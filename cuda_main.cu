#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cuda.h>

#define PARAM 6 //l, spectrum (bloom filter), dim, input, dim, output
#define BASES 4
#define BLOCK_DIM 16
#define READS_LENGHT 35

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
__shared__ unsigned short int gpu_l;
__shared__ unsigned int gpu_inputDim;

/* Allocate vector of bases on device memory,
 * used to substitute character in the attempt of correcting reads
 * In order to reduce acces latency it will be stored on shared memory
 */
__shared__ char bases[BASES];

/****************  DEVICE FUNCTIONS  *********************/
__device__ ushort2 matrix_maximum(unsigned short int ** v) {
   unsigned short int i, j, maximum=0;
   ushort2 couple;
   for(i=0; i<READS_LENGHT; i++) {
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

__global__ void fixing(char ** reads, unsigned short int *** voting_matrix_array) {
   ushort2 trim_indexes;
   trim_indexes.x = 0; //Starting index of longest substring
   trim_indexes.y = 0; //End index of longest substring
   
   int idx = blockIdx.x * blockDim.x + threadIdx.x;
   
   if(idx >= gpu_inputDim)
      return;
   
   char * read;
   unsigned short int ** voting_matrix;
   
   read = (char *)reads[idx];
   voting_matrix = (unsigned short int **)voting_matrix_array[idx];
   
   ushort2 couple = matrix_maximum(voting_matrix);
   
   if(voting_matrix[couple.x][couple.y] == 0) {
      return; //read is already correct
   }
   
   char rc[READS_LENGHT+1];
   unsigned short int i;
   for(i=0; i<couple.x; i++) {
      rc[i] = read[i];
   }
   rc[i] = bases[couple.y];
   for(i=couple.x+1; i<READS_LENGHT; i++) {
      rc[i] = read[i];
   }
   
   bool corrected_flag=1, trimmed_flag=0;
   unsigned short int j;
   char tuple[READS_LENGHT + 1]; //Memory waste - TODO
   /*
   char * tuple;
   if(cudaMalloc((void **)&tuple, (gpu_l + 1) * sizeof(char)) == cudaErrorMemoryAllocation) {
      fprintf(stdout, "Error: CUDA allocation\n");
      return; //How to exit???
   }
   */
   for(j=0;j < (READS_LENGHT-(gpu_l+1)); j++) {
      //Create tuple
      for(i=j; i<j+gpu_l; i++) {
         tuple[i-j] = read[j];
      }
      tuple[gpu_l] = '\0';
      if( !(1/* TODO - query bloom filter for tuple */) ) {
         corrected_flag = 0;
        /* Check for trimming
         * If current subsequence is longer than previous one then update
         * Else the longest subsequent is already stored
         */
         if( (j+8 - trim_indexes.y+2) > (trim_indexes.y - trim_indexes.x) ) {
            trim_indexes.x = trim_indexes.y+2;
            trim_indexes.y = j+gpu_l-1;
         }
      }
      else {
         trimmed_flag = 1;
      }
   }
   
   if(corrected_flag) {
      gpu_strcpy(read, rc); //Return corrected read
      return;
   }
   
   if(trimmed_flag) {
      //Trim read
      for(i=trim_indexes.x, j=0; i<trim_indexes.y; i++, j++) {
         read[j] = rc[i];
      }
      read[j] = '\0';
      return;
   }
   
   //Uncorrect read, return empty string
   read[0] = '\0';
}



/********************** MAIN **************************/

int main (int argc, char * argv[]) {
   int i, j, h;
   
   if(argc < PARAM+1) {
      fprintf(stdout, "Error: parameters\n");
      exit(1);
   }
   
   const unsigned short int l = atoi(argv[1]);
   
   /************** INITIALIZATION **************/
   
   /* Allocate basis vector
    * 
    */
   /* WARNING:
    * cuda_main.cu(157): warning: address of a __shared__ variable "bases"
    * cannot be directly taken in a host function
    * 
    * Also on lines 231, 237, 325 and 331
    */
   char host_bases[BASES] = {'A', 'C', 'G','T'};
   cudaMemcpyToSymbol((const char *)bases, (const void *)host_bases, sizeof(char) * BASES);

   /* Allocate and store the spectrum on the host,
    * reading it from the input file
    * 
    */
   //TODO Spectrum variable declaration
   
   FILE * spectrumFP = fopen(argv[2], "r");
   
   //TODO Read spectrum bloom filter
   
   fclose(spectrumFP);
   
   const unsigned int inputDim = atoi(argv[5]);
   
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
      if(!(reads[i] = (char *)malloc(sizeof(char) * (READS_LENGHT+1)))) {
         fprintf(stdout, "Error: allocation\n");
         exit(1);
      }
   }
   
   FILE * readsFP = fopen(argv[4], "r");
   
   for(i=0; i<inputDim; i++) {
      fgets(reads[i], READS_LENGHT+1, readsFP);
      printf("%s\n", reads[i]);
   }
   
   fclose(readsFP);
   
   /************* CUDA ALLOCATION ***************/
   
   /* Allocate spectrum (on texture memory?)
    * Inlcude memcopy of the spectrum data
    * 
    * 
    */
   //TODO -  use __constant__ qualifier
    
   /* Allocate reads on device memory as gpu_reads
    * Include memcopy of already filled data
    * 
    */
   char ** gpu_reads;
   HANDLE_ERROR(cudaMalloc((void **)&gpu_reads, inputDim * sizeof(char *)));
   for(i=0; i<inputDim; i++) {
      HANDLE_ERROR(cudaMalloc((void **)&(gpu_reads[i]), (READS_LENGHT + 1) * sizeof(char)));
   }
   for(i=0; i<inputDim; i++) {
      HANDLE_ERROR(cudaMemcpy(gpu_reads[i], reads[i], sizeof(char) * (READS_LENGHT + 1), cudaMemcpyHostToDevice));
   }
   
   /* Allocate inputDim on gpu memory as gpu_inputDim
    * Include memcopy
    * 
    */
   cudaMemcpyToSymbol((const char *)&gpu_inputDim, (const void *)&inputDim, sizeof(unsigned int));
   
   /* Allocate l on device memory
    * Include memcopy
    * 
    */
   cudaMemcpyToSymbol((const char *)&gpu_l, (const void *)&l, sizeof(unsigned short int));
   
   /* Initialize voting_matrix of zeros
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
      if(!(voting_matrix[i] = (unsigned short int **)malloc(sizeof(unsigned short int *) * READS_LENGHT))) {
         fprintf(stdout, "Error: allocation\n");
         exit(1);
      }
      for(j=0; j<READS_LENGHT; j++) {
         if(!(voting_matrix[i][j] = (unsigned short int *)malloc(sizeof(unsigned short int) * BASES))) {
            fprintf(stdout, "Error: allocation\n");
            exit(1);
         }
         for(h=0; h<BASES; h++) {
            voting_matrix[i][j][h] = 0;
         }
      }
   }
   
   //Allocate voting_matrix on device as gpu_voting_matrix
   unsigned short int *** gpu_voting_matrix;
   if(cudaMalloc((void **)&gpu_voting_matrix, sizeof(unsigned short int ***) * inputDim) == cudaErrorMemoryAllocation) {
      fprintf(stdout, "Error: CUDA allocation\n");
      exit(1);
   }
   for(i=0; i<inputDim; i++) {
      HANDLE_ERROR(cudaMalloc((void **)&gpu_voting_matrix[i], sizeof(unsigned short int **) * READS_LENGHT));
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
    * TODO spectrum allocation and filling (in proper section)
    */
   
   
   /************* FIXING ***************/
   
   /* Assume gpu_voting_matrix already computed and stored in device memory,
    * as gpu_reads. Do not remove from the device memory
    * 
    */
   
   //Execute kernel
   fixing <<< inputDim/BLOCK_DIM, BLOCK_DIM >>> (gpu_reads, gpu_voting_matrix);
   
   /************ RETRIEVE RESULT ********/
   
   /* Need to retrieve only the gpu_reads, gpu_voting_matrix is not needed anymore
    * 
    * 
    */
    for(i=0; i<inputDim; i++) {
      HANDLE_ERROR(cudaMemcpy(reads[i], gpu_reads[i], sizeof(char) * (READS_LENGHT + 1), cudaMemcpyDeviceToHost));
    } 
   
   
   /************* MEM FREE *****************/
   /* Only not freed yet */
   
   cudaFree((void **)&gpu_inputDim);
   for(i=0; i<inputDim; i++) {
      cudaFree(gpu_reads[i]);
      free(reads[i]);
   }
   cudaFree((void **)&gpu_reads);
   cudaFree((void **)&gpu_l);
   cudaFree((void **)&gpu_voting_matrix);
   free(voting_matrix);
   
   /************* WRITE OUT RESULTS ********/
   
   FILE * outFP;
   outFP = fopen(argv[6], "w");
   
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
