#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
//#include <cuda.h>
#include <stdint.h>
//#include <device_functions.h>
//#include <cuda_runtime_api.h>
//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"

#ifdef SILENT
   #define printf(...)
#endif

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

/******************** MAIN *******************************/

int main (int argc, char * argv[]) {
   int i, j;
   
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
   char * reads;
   if(!(reads = (char *)malloc(sizeof(char) * (READS_LENGTH+1) * inputDim))) {
      fprintf(stdout, "Error: allocation\n");
      exit(1);
   }
   char cons;
   
   for(i=0; i<inputDim; i++) {
      for(j=0; j<READS_LENGTH; j++) {
         fscanf(readsFP, "%c", reads + (READS_LENGTH+1) * i + j);
      }
      *(reads + (READS_LENGTH+1) * i + j) = '\0';
      //printf("%s\n", reads + (READS_LENGTH+1)*i);
      fscanf(readsFP, "%c", &cons);
   }
   
   fclose(readsFP);
   
   /************* CUDA ALLOCATION ***************/
   
   /* Allocate spectrum on texture memory
    * Inlcude memcopy of the spectrum data
    * 
    * 
    */
   uint64_t * gpu_hashed_spectrum;
   
   /* Allocate reads on device memory as gpu_reads
    * Include memcopy of already filled data
    * 
    */
   char * gpu_reads;
   HANDLE_ERROR(cudaMalloc(&gpu_reads, sizeof(char) * inputDim * (READS_LENGTH+1)));
   HANDLE_ERROR(cudaMemcpy(gpu_reads, reads, sizeof(char) * inputDim * (READS_LENGTH+1), cudaMemcpyHostToDevice));
   
   /************ KERNEL **************************/
   
   
   
   
   /*********** READ BACK ************************/
   
   HANDLE_ERROR(cudaMemcpy(reads, gpu_reads, sizeof(char) * inputDim * (READS_LENGTH+1), cudaMemcpyDeviceToHost));
   
   /*********** WRITE OUT ************************/
   
   FILE * outFP = fopen(argv[3], "w+");
   
   for(i=0; i<inputDim; i++) {
      for(j=0; j<READS_LENGTH; j++) {
         fprintf(outFP, "%c", *(reads + (READS_LENGTH+1) * i + j));
      }
      fprintf(outFP, "\n");
   }
   
   fclose(outFP);
   
   free(reads);
   free(hashed_spectrum);
   
   return 0;
}
