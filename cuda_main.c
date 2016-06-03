#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#define PARAM 6 //l, spectrum (bloom filter), dim, input, dim, output
#define BASES 4
#define BLOCK_DIM 16
#define READS_LENGHT 35

/****************  DEVICE FUNCTIONS  *********************/
__device__ ushort2 matrix_maximum(int ** v) {
   unsigned short int i, j, max=0;
   ushort2 couple;
   for(i=0; i<READS_LENGHT; i++) {
      for(j=0;j<BASES;j++) {
         if(v[i][j] > max) {
            max = v[i][j];
            couple.x=i;
            copule.y=j;
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


/***************** FIXING KERNEL ***************************/

__global__ void fixing(char ** reads, unsigned short int l, unsigned short int ** voting_matrix, unsigned int inputDim) {
   ushort2 couple = matrix_maximum(voting_matrix);
   ushort2 trim_indexes;
   trim_indexes.x = 0; //Starting index of longest substring
   trim_indexes.y = 0; //End index of longest substring
   
   int idx = blockIdx.x * blockDim.x + threadIdx.x;
   
   if(idx >= inputDim)
      return;
   
   __shared__ char read[L+1];
   
   read = reads[idx];
   
   if(voting_matrix[couple.x][couple.y] == 0) {
      return; //read is already correct
   }
   
   __shared__ char rc[READS_LENGHT+1];
   unsigned short int i;
   for(i=0; i<couple.x; i++) {
      rc[i] = read[i];
   }
   rc[i] = bases(copule.y);
   for(i=couple.x+1; i<READS_LENGHT; i++) {
      rc[i] = read[i];
   }
   
   bool corrected_flag=TRUE, trimmed_flag=FALSE;
   unsigned short int j;
   __shared__ char tuple[l+1];
   for(j=0;j < (READS_LENGHT-(l+1)); j++) {
      //Create tuple
      for(i=j; i<j+l; i++) {
         tuple[i-j] = read[j];
      }
      tuple[l] = '\0';
      if( !(/* query bloom filter for tuple */) ) {
         corrected_flag = FALSE;
        /* Check for trimming
         * If current subsequence is longer than previous one then update
         * Else the longest subsequent is already stored
         */
         if( (j+8 - trim_indexes.y+2) > (trim_indexes.y - trim_indexes.x) ) {
            trim_indexes.x = trim_indexes.y+2;
            trim_indexes.y = j+l-1;
         }
      }
      else {
         trimmed_flag = TRUE;
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
   int l;
   int i;
   
   if(argc < PARAM+1) {
      fprintf(stdout, "Error: parameters\n");
      exit(1);
   }
   
   l = atoi(argv[1]);

   FILE * spectrumFP = fopen(argv[2], "r");
   
   //TODO Read spectrum bloom filter
   
   fclose(spectrumFP);
   
   unsigned int inputDim = atoi(argv[5]);
   
   char ** reads;
   if(!(reads = (char **)malloc(sizeof(char *) * inputDim))) {
      frpintf(stdout, "Error: allocation\n");
      exit(1);
   }
   
   for(i=0; i<inputDim; i++) {
      if(!(reads[i] = (char *)malloc(sizeof(char) * (READS_LENGHT+1)))) {
         frpintf(stdout, "Error: allocation\n");
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
   
   //TODO Allocate spectrum on texture memory
   
   //Allocate reads on device memory
   char *** gpu_reads;
   if(cudaMalloc(&gpu_reads, inputDim * sizeof(char *)) == cudaErrorMemoryAllocation) {
      frpintf(stdout, "Error: CUDA allocation\n");
      exit(1);
   }
   for(i=0; i<inputDim; i++) {
      if(cudaMalloc(&(gpu_reads[i]), (READS_LENGHT + 1) * sizeof(char)) == cudaErrorMemoryAllocation) {
         fprintf(stdout, "Error: CUDA allocation\n");
         exit(1);
      }
   }
   for(i=0; i<inputDim; i++) {
      cudaMemcpy(gpuReads[i], reads[i], sizeof(char) * (READS_LENGHT + 1), cudaMemcpyHostToDevice);
   }
   
   //Allocate inputDim on gpu memory
   if(cudaMalloc(&gpu_inputDim, sizeof(unsigned int)) == cudaErrorMemoryAllocation) {
      fprintf(stdout, "Error: CUDA allocation\n");
      exit(1);
   }
   cudaMemcpy(gpu_inputDim, inputDim, sizeof(unsigned int), cudaMemcpyHostToDevice);
   
   //Allocate l on device memory
   unsigned short int * gpu_l;
   if(cudaMalloc(&gpu_l, sizeof(unsigned short int)) == cudaErrorMemoryAllocation) {
      fprintf(stdout, "Error: CUDA allocation\n");
      exit(1);
   }
   cudaMemcpy(gpu_l, l, sizeof(unsigned short int), cudaMemcpyHostToDevice);
   
   /************* VOTING ***************/
   
   
   
   /************* FIXING ***************/
   
   __device__ const thrust::device_vector<char> bases(BASES) = {'A', 'C', 'G'.,'T'};
   
   //Execute kernel
   fixing <<< inputDim/BLOCK_DIM, BLOCK_DIM >>> (gpu_reads, gpu_l, gpu_voting_matrix, gpu_inputDim);
   
   
   /************* MEM FREE *****************/
   /* Only not freed yet */
   
   cudaFree(gpu_inputDim);
   for(i=0; i<inputDim; i++) {
      cudaFree(gpu_reads[i]);
      free(reads[i]);
   }
   cudaFree(gpu_reads);
   cudaFree(gpu_l);
   free(reads);
   
   return 0;
}
