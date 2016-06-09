//#ifndef __CUDACC__
//#define __CUDACC__
//#endif

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
//#include <cuda.h>
#include <stdint.h>
//#include <device_functions.h>
//#include <cuda_runtime_api.h>
//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"

#define PARAM 3 //input (bloom filter), input (reads), output
#define BASES 4
#define BLOCK_DIM 16
#define DATA_PER_THREAD 10
#define READS_LENGTH 35
#define L 10

/************ ERROR HANDLING *****************************/
void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

/**************** GLOBAL VARIABLE ************************/

/* Allocate vector of bases on device memory,
 * used to substitute character in the attempt of correcting reads
 * In order to reduce acces latency it will be stored on shared memory
 */
__shared__ char bases[BASES];

/****************  DEVICE FUNCTIONS  *********************/
__device__ ushort2 matrix_maximum(unsigned short int ** v) {
   unsigned short int i, j, maximum=0;
   ushort2 couple;
   for(i=0; i<READS_LENGTH; i++) {
      for(j=0;j<BASES;j++) {
         if(v[i][j] > maximum) {
            maximum = v[i][j];
            couple.x=i;
            couple.y=j;
         }
      }
   }
   
   return couple;
}

__device__ void gpu_strcpy(char * a, char * b) {
   unsigned short int i=0;
   do {
      a[i] = b[i];
      i++;
   } while (b[i-1] != '\0');
}


/**************** KERNEL DEFINITION ************************/

/***************** VOTING KERNEL ***************************/

//TODO

/***************** FIXING KERNEL ***************************/

__global__ void fixing(char ** reads, unsigned short int *** voting_matrix_array, unsigned int inputDim, uint64_t * gpu_hashed_spectrum) {
   ushort2 trim_indexes;
   trim_indexes.x = 0; //Starting index of longest substring
   trim_indexes.y = 0; //End index of longest substring
   short int i, h;
   
   int idx = blockIdx.x * blockDim.x + threadIdx.x;
   
   /* If first thread in the block should copy in shared memory
    * the variables from the global memory
    * 
    */
   if(threadIdx.x == 0) {
      bases[0] = 'A';
      bases[1] = 'C';
      bases[2] = 'G';
      bases[3] = 'T';
   }
   __syncthreads();
   
   char * read;
   unsigned short int ** voting_matrix;
   char rc[READS_LENGTH+1];
   unsigned short int j;
   char tuple[L + 1];
   bool corrected_flag, trimmed_flag;
   ushort2 couple;
   
   for(h=0; h<DATA_PER_THREAD; h++) {
      if(idx + inputDim/DATA_PER_THREAD * h >= inputDim)
         return;
      read = (char *)reads[idx + inputDim/DATA_PER_THREAD * h];
      voting_matrix = (unsigned short int **)voting_matrix_array[idx];
      
      couple = matrix_maximum(voting_matrix);
      
      if(voting_matrix[couple.x][couple.y] == 0) {
         continue; //read is already correct
      }
      
      for(i=0; i<couple.x; i++) {
         rc[i] = read[i];
      }
      rc[i] = bases[couple.y];
      for(i=couple.x+1; i<READS_LENGTH; i++) {
         rc[i] = read[i];
      }
      
      corrected_flag=1;
      trimmed_flag=0;
      
      for(j=0;j < (READS_LENGTH-(L+1)); j++) {
         //Create tuple
         for(i=j; i<j + L; i++) {
            tuple[i-j] = read[j];
         }
         tuple[L] = '\0';
         if( !(1/* TODO - query bloom filter for tuple */) ) {
            corrected_flag = 0;
           /* Check for trimming
            * If current subsequence is longer than previous one then update
            * Else the longest subsequent is already stored
            */
            if( (j+8 - trim_indexes.y+2) > (trim_indexes.y - trim_indexes.x) ) {
               trim_indexes.x = trim_indexes.y+2;
               trim_indexes.y = j + L - 1;
            }
         }
         else {
            trimmed_flag = 1;
         }
      }
      
      if(corrected_flag) {
         gpu_strcpy(reads[idx + inputDim/DATA_PER_THREAD * h], rc); //Return corrected read
         continue;
      }
      
      if(trimmed_flag) {
         //Trim read
         for(i=trim_indexes.x, j=0; i<trim_indexes.y; i++, j++) {
            read[j] = rc[i];
         }
         read[j] = '\0';
         gpu_strcpy(reads[idx + inputDim/DATA_PER_THREAD * h], read);
         continue;
      }
      
      //Uncorrect read, return empty string
      reads[idx + inputDim/DATA_PER_THREAD * h][0] = '\0';
   }
}



/********************** MAIN **************************/

int main (int argc, char * argv[]) {
   int i, j, h;
   
   if(argc < PARAM+1) {
      fprintf(stdout, "Error: parameters\n");
      exit(1);
   }
   
   /************** INITIALIZATION **************/

   /* Allocate and store the spectrum on the host,
    * reading it from the input file
    * 
    */
   uint64_t * hashed_spectrum;
   unsigned int spectrum_size;
   
   FILE * spectrumFP = fopen(argv[1], "r");
   fscanf(spectrumFP, "%u\n", &spectrum_size);
   
   if(!(hashed_spectrum = (uint64_t *)malloc(sizeof(uint64_t) * spectrum_size))) {
      fprintf(stdout, "Error: allocation\n");
      exit(1);
   }
   
   for(i=0; i<spectrum_size; i++) {
      fscanf(spectrumFP, "%lu", &hashed_spectrum[i]);
      //printf("%lu\n", hashed_spectrum[i]);
   }
   
   fclose(spectrumFP);
   
   unsigned int inputDim;
   
   FILE * readsFP = fopen(argv[2], "r");
   fscanf(readsFP, "%u\n", &inputDim);
   
   /* Allocate and read from input file the sequence reads,
    * on host memory
    * 
    */
   char ** reads;
   if(!(reads = (char **)malloc(sizeof(char *) * inputDim))) {
      fprintf(stdout, "Error: allocation\n");
      exit(1);
   }
   
   for(i=0; i<inputDim; i++) {
      if(!(reads[i] = (char *)malloc(sizeof(char) * (READS_LENGTH+1)))) {
         fprintf(stdout, "Error: allocation\n");
         exit(1);
      }
   }
   
   for(i=0; i<inputDim; i++) {
      //fgets(reads[i], READS_LENGTH+1, readsFP);
      fscanf(readsFP, "%s\n", reads[i]);
      //printf("%s\n", reads[i]);
   }
   
   fclose(readsFP);
   
   /************* CUDA ALLOCATION ***************/
   
   /* Allocate spectrum on texture memory
    * Inlcude memcopy of the spectrum data
    * 
    * 
    */
   uint64_t * gpu_hashed_spectrum;
   HANDLE_ERROR(cudaMalloc((void **)&gpu_hashed_spectrum, sizeof(uint64_t *) * spectrum_size));
   for(i=0; i<spectrum_size; i++) {
      HANDLE_ERROR(cudaMalloc((void **)&gpu_hashed_spectrum[i], sizeof(uint64_t)));
      HANDLE_ERROR(cudaMemcpy((void *)&gpu_hashed_spectrum[i], (const void *)&hashed_spectrum[i], sizeof(char), cudaMemcpyHostToDevice));
   }
    
   /* Allocate reads on device memory as gpu_reads
    * Include memcopy of already filled data
    * 
    */
   char ** gpu_reads;
   HANDLE_ERROR(cudaMalloc((void **)&gpu_reads, inputDim * sizeof(char *)));
   for(i=0; i<inputDim; i++) {
      HANDLE_ERROR(cudaMalloc((void **)&(gpu_reads[i]), (READS_LENGTH + 1) * sizeof(char)));
   }
   for(i=0; i<inputDim; i++) {
      HANDLE_ERROR(cudaMemcpy(gpu_reads[i], reads[i], sizeof(char) * (READS_LENGTH + 1), cudaMemcpyHostToDevice));
   }
   
   /* Initialize voting_matrix
    * Then proceed with the allocation on device memory of gpu_voting_matrix
    * Include memcopy of the initialized matrix
    * 
    */
   unsigned short int *** voting_matrix;
   if(!(voting_matrix = (unsigned short int ***)malloc(sizeof(unsigned short int **) * inputDim))) {
      fprintf(stdout, "Error: allocation\n");
      exit(1);
   }
   for(i=0; i<inputDim; i++) {
      if(!(voting_matrix[i] = (unsigned short int **)malloc(sizeof(unsigned short int *) * READS_LENGTH))) {
         fprintf(stdout, "Error: allocation\n");
         exit(1);
      }
      for(j=0; j<READS_LENGTH; j++) {
         if(!(voting_matrix[i][j] = (unsigned short int *)malloc(sizeof(unsigned short int) * BASES))) {
            fprintf(stdout, "Error: allocation\n");
            exit(1);
         }
      }
   }
   
   //Allocate voting_matrix on device as gpu_voting_matrix
   unsigned short int *** gpu_voting_matrix;
   HANDLE_ERROR(cudaMalloc((void **)&gpu_voting_matrix, sizeof(unsigned short int ***) * inputDim));
   for(i=0; i<inputDim; i++) {
      HANDLE_ERROR(cudaMalloc((void **)&gpu_voting_matrix[i], sizeof(unsigned short int **) * READS_LENGTH));
      for(j=0; j<inputDim; j++) {
         HANDLE_ERROR(cudaMalloc((void **)&gpu_voting_matrix[i][j], sizeof(unsigned short int *) * BASES));
         for(h=0; h<BASES; h++) {
            HANDLE_ERROR(cudaMalloc((void **)&gpu_voting_matrix[i][j][h], sizeof(unsigned short int)));
            //Copy initialized matrix
            HANDLE_ERROR(cudaMemcpy((void *)&gpu_voting_matrix[i][j][h], (const void *)&voting_matrix[i][j][h], sizeof(unsigned short int), cudaMemcpyHostToDevice));
         }
      }
   }
   
   /************* VOTING ***************/
   
   /* Already allocated and initializated voting_matrix and gpu_voting_matrix;
    * reads and gpu_reads allocated and filled with data;
    * inputDim is defined as gpu_inputDim;
    */
   
   
   /************* FIXING ***************/
   
   /* Assume gpu_voting_matrix already computed and stored in device memory,
    * as gpu_reads. Do not remove from the device memory
    * 
    */
   
   //Execute kernel
   fixing <<< inputDim/BLOCK_DIM/DATA_PER_THREAD, BLOCK_DIM >>> (gpu_reads, gpu_voting_matrix, inputDim, gpu_hashed_spectrum);
   
   /************ RETRIEVE RESULT ********/
   
   /* Need to retrieve only the gpu_reads, gpu_voting_matrix is not needed anymore
    * 
    * 
    */
    for(i=0; i<inputDim; i++) {
      HANDLE_ERROR(cudaMemcpy(reads[i], gpu_reads[i], sizeof(char) * (READS_LENGTH + 1), cudaMemcpyDeviceToHost));
    } 
   
   
   /************* MEM FREE *****************/
   /* Only not freed yet */
   
   for(i=0; i<inputDim; i++) {
      cudaFree(gpu_reads[i]);
      free(reads[i]);
   }
   cudaFree((void **)&gpu_reads);
   cudaFree((void **)&gpu_voting_matrix);
   free(voting_matrix);
   
   /************* WRITE OUT RESULTS ********/
   
   FILE * outFP;
   outFP = fopen(argv[3], "w+");
   
   for(i=0; i<inputDim; i++) {
      //Ignore discarded reads
      if(strcmp(reads[i], "") != 0) {
         fprintf(outFP, "%s\n", reads[i]);
      }
   }
   
   fclose(outFP);
   free(reads);
   
   return 0;
}
