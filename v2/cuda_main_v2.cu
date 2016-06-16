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

/**************** PROTOTYPES *****************************/
__device__ int CheckBit(uint64_t *, unsigned long, unsigned int);
__device__ int CheckHash(uint64_t *, char *, unsigned int);
__device__ uint64_t MurmurHash64A (const void *, int, unsigned int);

/**************** GLOBAL VARIABLE ************************/

/* Allocate vector of bases on device memory,
 * used to substitute character in the attempt of correcting reads
 * In order to reduce acces latency it will be stored on shared memory
 */
__shared__ char bases[BASES];
__shared__ unsigned int hashseeds[8];

/****************  DEVICE FUNCTIONS  *********************/

/* @requires v poining to actual matrix
 * 
 */
__device__ ushort2 matrix_maximum(unsigned short int * v) {
   unsigned short int i, j, maximum=0;
   ushort2 couple;
   couple.x=0;
   couple.y=0;
   for(i=0; i<READS_LENGTH; i++) {
      for(j=0;j<BASES;j++) {
         if(*(v + i*BASES + j) >= maximum) {
            maximum = *(v + i*BASES + j);
            couple.x=i;
            couple.y=j;
         }
      }
   }
   
   return couple;
}

/* @requires v poining to actual string
 * 
 */
__device__ void gpu_strcpy(char * a, char * b) {
   unsigned short int i=0;
   do {
      *(a+i) = *(b+i);
      i++;
   } while (*(b+i-1) != '\0' && i<READS_LENGTH);
}

__device__ void gpu_strncpy(char * a, char* b, unsigned short int l){
    unsigned short int i=0;
    do{
        *(a+i) = *(b+i);
        i++;
    }
    while(i<l && *(b+i-1) != '\0');
}

/************************* KERNEL ****************************************************/

/***************** VOTING KERNEL ***************************/

__global__ void voting(char * reads, unsigned short int * voting_matrix_array, unsigned int inputDim, uint64_t * gpu_hashed_spectrum, unsigned int spectrum_size){ 
 
	char * read;
	char subseq[L+1], subcpy[L+1]; //substring of dimension L, copy with substituted character
	unsigned short int * voting_matrix;
	int j, c;
	unsigned short h;
	
	int idx = blockDim.x * blockIdx.x + threadIdx.x;
	
	if( threadIdx.x == 0){
	    bases[0] = 'A';
	    bases[1] = 'C';
	    bases[2] = 'G';
	    bases[3] = 'T';
	    
	    hashseeds[0] = 19;
	    hashseeds[1] = 98638;
	    hashseeds[2] = 8229873;
	    hashseeds[3] = 42;
	    hashseeds[4] = 3275875;
	    hashseeds[5] = 6757572;
	    hashseeds[6] = 678;
	    hashseeds[7] = 1234567890;
	}
	__syncthreads();
	
	
	for(h=0; h<DATA_PER_THREAD; h++){
	
	    if(inputDim/DATA_PER_THREAD*h + idx >= inputDim)
	        return;
	        
	    /* Get data
	     *
	     */
        voting_matrix = voting_matrix_array + (inputDim/DATA_PER_THREAD*h + idx)*BASES*(READS_LENGTH);
        read = reads + (inputDim/DATA_PER_THREAD*h + idx)*(READS_LENGTH+1);
	    
	    /* Initialize voting matrix to zeroes
	     *
	     */ 
	    for(j=0; j<READS_LENGTH; j++)
		    for(c=0; c<BASES; c++)
			    *(voting_matrix + j*BASES + c) = 0;
	    
	    /* Voting procedure
	     * 
	     */
	    for(j=0; j<(READS_LENGTH-L+1); j++){
		    gpu_strncpy(subseq, read+j, L);  //take substring of length L
		    subseq[L] = '\0';		    
		    if(!(CheckHash(gpu_hashed_spectrum, subseq, spectrum_size))){  //begin voting procedure for non-solid tuples
			    for(int p=0; p<L; p++){
				    for(c=0; c<BASES; c++){
					    if(subseq[p] != bases[c]){
						    gpu_strcpy(subcpy, subseq); 
						    subcpy[p]=bases[c];	 //substitute current base with any other one		   
							if(CheckHash(gpu_hashed_spectrum, subcpy, spectrum_size))
							    *(voting_matrix + (j+p)*BASES + c) = *(voting_matrix + (j+p)*BASES + c) + 1;  //increment number of solid tuples with that correction
					    }
				    }
			    }
		    }
	    }
	    
	    
	    
	}
}

/***************** FIXING KERNEL ***************************/

__global__ void fixing(char * reads, unsigned short int * voting_matrix_array, unsigned int inputDim, uint64_t * gpu_hashed_spectrum, unsigned int spectrum_size) {
   unsigned short int i,h,k;
   
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
      
      hashseeds[0] = 19;
	  hashseeds[1] = 98638;
	  hashseeds[2] = 8229873;
	  hashseeds[3] = 42;
	  hashseeds[4] = 3275875;
	  hashseeds[5] = 6757572;
	  hashseeds[6] = 678;
	  hashseeds[7] = 1234567890;
   }
   __syncthreads();
   
   char * read_ptr;
   char read[READS_LENGTH+1];
   unsigned short int * voting_matrix;
   char rc[READS_LENGTH+1];
   unsigned short int j;
   char tuple[L+1];
   bool corrected_flag, trimmed_flag;
   ushort2 couple;
   
   for(h=0; h<DATA_PER_THREAD; h++) {
      if(idx + inputDim/DATA_PER_THREAD * h >= inputDim)
         return;
      read_ptr = reads + (idx + inputDim/DATA_PER_THREAD * h) * (READS_LENGTH+1);
      gpu_strcpy(read, read_ptr);
      voting_matrix = voting_matrix_array + (idx + inputDim/DATA_PER_THREAD * h) * READS_LENGTH * BASES;
      
      couple = matrix_maximum(voting_matrix);
      
      if(*(voting_matrix + couple.x * BASES + couple.y) == 0) {
         continue; //read is already correct
      }
      
      gpu_strcpy(rc, read);
      rc[READS_LENGTH] = '\0';
      rc[couple.x] = bases[couple.y];
      
      corrected_flag=1;
      trimmed_flag=0;
      
      for(j=0;j < (READS_LENGTH-L+1); j++) {
         //Create tuple
         gpu_strncpy(tuple, rc+j, L);
         tuple[L] = '\0';
         if( !(CheckHash(gpu_hashed_spectrum, tuple, spectrum_size)) ) { /* Query bloom filter for tuple */
            corrected_flag = 0;
         }
         else {
            trimmed_flag = 1;
         }
      }
      
      if(corrected_flag) {
         gpu_strcpy(read_ptr, rc); //Return corrected read
         continue;
      }
      
      bool flag=0;
      int max_length=0;
      if(trimmed_flag) {
         //Trim read
         for(k=0; k < READS_LENGTH-L+1; k++) {
            for(j=READS_LENGTH-k; j>L+1 && !flag; j--) {
               gpu_strncpy(rc, read+k, j);
               rc[j] = '\0';
               //Check each tuple of shorter string
               flag=1;
               for(i=0; i<j-L+1; i++) {
                  gpu_strncpy(tuple, read+i, L);
                  tuple[L] = '\0';
                  flag = flag && (CheckHash(gpu_hashed_spectrum, tuple, spectrum_size));
               }
               //Copy trimmered read
               if(flag) {
                  if(j-k > max_length) {
                     gpu_strncpy(read_ptr, rc, j);
                     *(read_ptr+j-k+1) = '\0';
                     max_length = j-k;
                  }
               }
            }
            flag=0;
         }
         continue;
      }
      
      //Uncorrect read, return empty string
      *(read_ptr) = '\0';
   }
}

/******************** MAIN ***********************************************************/

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
   HANDLE_ERROR(cudaMalloc(&gpu_hashed_spectrum, sizeof(uint64_t) * spectrum_size));
   HANDLE_ERROR(cudaMemcpy(gpu_hashed_spectrum, hashed_spectrum, sizeof(uint64_t) * spectrum_size, cudaMemcpyHostToDevice));
   
   /* Allocate reads on device memory as gpu_reads
    * Include memcopy of already filled data
    * 
    */
   char * gpu_reads;
   HANDLE_ERROR(cudaMalloc(&gpu_reads, sizeof(char) * inputDim * (READS_LENGTH+1)));
   HANDLE_ERROR(cudaMemcpy(gpu_reads, reads, sizeof(char) * inputDim * (READS_LENGTH+1), cudaMemcpyHostToDevice));
   
   /* Allocate space gpu_voting_matrix
    * 
    */
   unsigned short int * gpu_voting_matrix;
   HANDLE_ERROR(cudaMalloc(&gpu_voting_matrix, sizeof(unsigned short int) * inputDim * READS_LENGTH * BASES));
   
   
   /************ KERNEL **************************/
   
   /************* VOTING ***************/
   
   /* Already allocated and initializated voting_matrix and gpu_voting_matrix;
    * reads and gpu_reads allocated and filled with data;
    * inputDim is defined as gpu_inputDim;
    */
   
   //Execute kernel
   voting <<< inputDim/BLOCK_DIM/DATA_PER_THREAD, BLOCK_DIM >>> (gpu_reads, gpu_voting_matrix, inputDim, gpu_hashed_spectrum, spectrum_size);   
   
   /*
   //DEBUG
   
   unsigned short int * v;
   v = (unsigned short int *)malloc(sizeof(unsigned short int) * inputDim * READS_LENGTH * BASES);
   HANDLE_ERROR(cudaMemcpy(v, gpu_voting_matrix, sizeof(unsigned short int) * inputDim * READS_LENGTH * BASES, cudaMemcpyDeviceToHost));
   
   for(int k=0; k<inputDim; k++) {
      for(i=0; i<BASES; i++) {
         for(j=0; j<READS_LENGTH; j++) {
            printf("%3d", *(v + j*BASES + i));
            //if (*(v + j*BASES + i)) printf("%3d", *(v + j*BASES + i));
         }
         printf("\n");
      }
      printf("\n\n\n");
   }
   
   //END DEBUG
   */
   
   
   
   /************* FIXING ***************/
   
   /* Assume gpu_voting_matrix already computed and stored in device memory,
    * as gpu_reads. Do not remove from the device memory
    * 
    */
   
   //Execute kernel
   fixing <<< inputDim/BLOCK_DIM/DATA_PER_THREAD, BLOCK_DIM >>> (gpu_reads, gpu_voting_matrix, inputDim, gpu_hashed_spectrum, spectrum_size);
   
   
   /*********** READ BACK ************************/
   
   HANDLE_ERROR(cudaMemcpy(reads, gpu_reads, sizeof(char) * inputDim * (READS_LENGTH+1), cudaMemcpyDeviceToHost));
      
   /*
   //DEBUG hashed_spectrum
   HANDLE_ERROR(cudaMemcpy(hashed_spectrum, gpu_hashed_spectrum, sizeof(uint64_t) * spectrum_size, cudaMemcpyDeviceToHost));
   for(i=0; i<spectrum_size;i++)
      printf("%lu\n", *(hashed_spectrum + i));
   */
   
   /*********** WRITE OUT ************************/
   
   FILE * outFP = fopen(argv[3], "w+");
   
   for(i=0; i<inputDim; i++) {
      if( 1/* *(reads + (READS_LENGTH+1) * i) != '\0' */) { //If not discarded
         for(j=0; j<READS_LENGTH && *(reads + (READS_LENGTH+1) * i + j) != '\0'; j++) {
            fprintf(outFP, "%c", *(reads + (READS_LENGTH+1) * i + j));
         }
         fprintf(outFP, "\n");
      }
   }
   
   fclose(outFP);
   
   free(reads);
   free(hashed_spectrum);
   cudaFree(gpu_reads);
   cudaFree(gpu_hashed_spectrum);
   cudaFree(gpu_voting_matrix);
   
   return 0;
}


/***************************** UTILITY FUNCTION **************************************************/
//Check that the bit is in the filter

__device__ int CheckBit(uint64_t* filter, unsigned long i, unsigned int n) {
	unsigned long k = i % n;
	unsigned long pos = i % 64;
	uint64_t bit = 1, res;
	bit = bit << pos;
	res = filter[k] & bit;
	if(res !=0)
		return 1;
	else
		return 0;
}

//Hashes the read and checks
__device__ int CheckHash(uint64_t* filter, char* read, unsigned int n) {
	int flag = 1;
	unsigned long i;
	int j;

	for(j=0; j<8 && flag == 1; j++){
		i = MurmurHash64A(read, L, hashseeds[j]);
		flag = CheckBit(filter, i, n);
	}

	return flag;
}


/****************** HASHES *******************/

//MurmurHash
__device__ uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end)
	{
		uint64_t k = *data++;

		k *= m; 
		k ^= k >> r; 
		k *= m; 
		
		h ^= k;
		h *= m; 
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7)
	{
	case 7: h ^= (uint64_t)data2[6] << 48;
	case 6: h ^= (uint64_t)data2[5] << 40;
	case 5: h ^= (uint64_t)data2[4] << 32;
	case 4: h ^= (uint64_t)data2[3] << 24;
	case 3: h ^= (uint64_t)data2[2] << 16;
	case 2: h ^= (uint64_t)data2[1] << 8;
	case 1: h ^= (uint64_t)data2[0];
	        h *= m;
	};
 
	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
} 
