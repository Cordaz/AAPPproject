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

#define MSEED 7127		//static seed for murmur function

//Spooky utils
#define SSEED 5449		//static seed for spooky function
#define SC_CONST 0xdeadbeefdeadbeefLL   
#define SC_NUMVARS 12              
#define SC_BLOCKSIZE (8 * SC_NUMVARS)
#define SC_BUFSIZE (2 * SC_BLOCKSIZE)

struct spooky_state
{
	uint64_t m_data[2 * SC_NUMVARS];
	uint64_t m_state[SC_NUMVARS];
	size_t m_length;
	unsigned char m_remainder;
};

//City utils
typedef uint8_t uint8;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef struct _uint128 uint128;
struct _uint128 {
  uint64 first;
  uint64 second;
};

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
__device__ unsigned long djb2_hash(unsigned char *);
__device__ uint64_t MurmurHash64A (const void *, int, unsigned int);
__device__ uint64_t APHash(char *, unsigned int);
__device__ uint64_t fnvhash(char *);
__device__ uint64_t SDBMHash(char *, unsigned int);
__device__ uint64_t RSHash(char *, unsigned int);

//City prototypes
#ifndef CITY_HASH_H_
#define CITY_HASH_H_
__device__ uint64 CityHash64(const char *buf, size_t len);
__device__ static uint64 UNALIGNED_LOAD64(const char *);
__device__ static uint32 UNALIGNED_LOAD32(const char *);
__device__ static uint64 Fetch64(const char *);
__device__ static uint32 Fetch32(const char *);
__device__ uint64 CityHash64(const char *, size_t);
__device__ static uint64 HashLen0to16(const char *, size_t);
__device__ static uint64 HashLen17to32(const char *, size_t);
__device__ static uint64 HashLen33to64(const char *, size_t);
__device__ static uint64 HashLen16(uint64, uint64);
__device__ static inline uint64 Hash128to64(const uint128);
__device__ static uint64 Rotate(uint64, int);
__device__ static uint64 RotateByAtLeast1(uint64, int);
__device__ static uint64 ShiftMix(uint64);
__device__ uint128 WeakHashLen32WithSeeds6(uint64, uint64, uint64, uint64, uint64, uint64);
__device__ uint128 WeakHashLen32WithSeeds(const char *, uint64, uint64);
#endif  // CITY_HASH_H_

//Spooky prototypes
__device__ static inline uint64_t rot64(uint64_t x, int k);
__device__ static inline void mix(const uint64_t *,	uint64_t *, uint64_t *, uint64_t *, uint64_t *,	uint64_t *, uint64_t *, uint64_t *,  uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *);
__device__ static inline void endPartial(uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *);
__device__ static inline void end(uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *, uint64_t *,	uint64_t *,	uint64_t *,	uint64_t *,	uint64_t *,	uint64_t *,	uint64_t *);
__device__ static inline void short_mix(uint64_t *,	uint64_t *,	uint64_t *,	uint64_t *);
__device__ static inline void short_end(uint64_t *,	uint64_t *,	uint64_t *,	uint64_t *);
__device__ void spooky_shorthash(const void *, size_t ,	uint64_t *,	uint64_t *);
__device__ void spooky_init(struct spooky_state *, uint64_t, uint64_t);
__device__ void spooky_update(struct spooky_state *, const void *, size_t);
__device__ void spooky_final(struct spooky_state *, uint64_t *, uint64_t *);
__device__ uint64_t spooky_hash64(const void *,	size_t,	uint64_t);
__device__ void spooky_hash128(const void *, size_t, uint64_t *, uint64_t *);

/**************** GLOBAL VARIABLE ************************/

/* Allocate vector of bases on device memory,
 * used to substitute character in the attempt of correcting reads
 * In order to reduce acces latency it will be stored on shared memory
 */
__shared__ char bases[BASES];

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

	for(int k=1; k<=8 && flag; k++){
		
		switch (k){
		
		case 1:
			i = djb2_hash((unsigned char*)read);
			break;
		case 2:
			i = MurmurHash64A(read, L, MSEED);
			break;
		case 3:
			i = APHash(read, L);
			break;
		case 4:
			i = CityHash64(read, L);
			break;
		case 5:
			i = spooky_hash64(read, L, SSEED);
			break;
		case 6:
			i = fnvhash(read);
			break;
		case 7:
			i = SDBMHash(read, L);
			break;
		case 8:
			i = RSHash(read, L);
			break;
		
		default: break;
		}
		
		flag = CheckBit(filter, i, n);
	}

	return flag;
}


/****************** HASHES *******************/

const uint64_t FNV_PRIME    = 1099511628211;				
const uint64_t OFFSET_BASIS = 14695981039346656037UL;

//djb2 hash
__device__ unsigned long djb2_hash(unsigned char *str){
    unsigned long hash = 5381;
    int c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

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

//AP hash
__device__ uint64_t APHash(char* str, unsigned int length) {
	uint64_t hash = 0xAAAAAAAA;
	unsigned int i = 0;

	for (i = 0; i < length; str++, i++)
	{
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ (*str) * (hash >> 3)) :
			(~((hash << 11) + ((*str) ^ (hash >> 5))));
	}

	return hash;
}

//FNV hash
__device__ uint64_t fnvhash(char * string){
   uint64_t hash = OFFSET_BASIS;
   for(uint8_t * c = (uint8_t*)string; *c != 0; ++c){
      hash ^= *c;
      hash *= FNV_PRIME;
   }
   return hash;
}

//SDBM hash
__device__ uint64_t SDBMHash(char* str, unsigned int length) {
	uint64_t hash = 0;
	unsigned int i = 0;

	for (i = 0; i < length; str++, i++)
	{
		hash = (*str) + (hash << 6) + (hash << 16) - hash;
	}

	return hash;
}

//RS hash
__device__ uint64_t RSHash(char* str, unsigned int length) {
	unsigned int b = 378551;
	unsigned int a = 63689;
	uint64_t hash = 0;
	unsigned int i = 0;

	for (i = 0; i < length; str++, i++)
	{
		hash = hash * a + (*str);
		a = a * b;
	}

	return hash;
}

/*************** CITY HASH **********************/

#define Uint128Low64(x) 	(x).first
#define Uint128High64(x)	(x).second

#if !defined(WORDS_BIGENDIAN)

#define uint32_in_expected_order(x) (x)
#define uint64_in_expected_order(x) (x)

#else

#ifdef _MSC_VER
#include <stdlib.h>
#define bswap_32(x) _byteswap_ulong(x)
#define bswap_64(x) _byteswap_uint64(x)

#elif defined(__APPLE__)
// Mac OS X / Darwin features
#include <libkern/OSByteOrder.h>
#define bswap_32(x) OSSwapInt32(x)
#define bswap_64(x) OSSwapInt64(x)

#else
#include <byteswap.h>
#endif

#define uint32_in_expected_order(x) (bswap_32(x))
#define uint64_in_expected_order(x) (bswap_64(x))

#endif  // WORDS_BIGENDIAN

#if !defined(LIKELY)
#if HAVE_BUILTIN_EXPECT
#define LIKELY(x) (__builtin_expect(!!(x), 1))
#else
#define LIKELY(x) (x)
#endif
#endif

static const uint64 k0 = 0xc3a5c85c97cb3127ULL;
static const uint64 k1 = 0xb492b66fbe98f273ULL;
static const uint64 k2 = 0x9ae16a3b2f90404fULL;
static const uint64 k3 = 0xc949d7c7509e6557ULL;

__device__ static uint64 UNALIGNED_LOAD64(const char *p) {
  uint64 result;
  memcpy(&result, p, sizeof(result));
  return result;
}

__device__ static uint32 UNALIGNED_LOAD32(const char *p) {
  uint32 result;
  memcpy(&result, p, sizeof(result));
  return result;
}

__device__ static uint64 Fetch64(const char *p) {
  return uint64_in_expected_order(UNALIGNED_LOAD64(p));
}

__device__ static uint32 Fetch32(const char *p) {
  return uint32_in_expected_order(UNALIGNED_LOAD32(p));
}

__device__ uint64 CityHash64(const char *s, size_t len) {
  if (len <= 32) {
    if (len <= 16) {
      return HashLen0to16(s, len);
    } else {
      return HashLen17to32(s, len);
    }
  } else if (len <= 64) {
    return HashLen33to64(s, len);
  }

  // For strings over 64 bytes we hash the end first, and then as we
  // loop we keep 56 bytes of state: v, w, x, y, and z.
  uint64 x = Fetch64(s + len - 40);
  uint64 y = Fetch64(s + len - 16) + Fetch64(s + len - 56);
  uint64 z = HashLen16(Fetch64(s + len - 48) + len, Fetch64(s + len - 24));
  uint64 temp;
  uint128 v = WeakHashLen32WithSeeds(s + len - 64, len, z);
  uint128 w = WeakHashLen32WithSeeds(s + len - 32, y + k1, x);
  x = x * k1 + Fetch64(s);

  // Decrease len to the nearest multiple of 64, and operate on 64-byte chunks.
  len = (len - 1) & ~(size_t)(63);
  do {
    x = Rotate(x + y + v.first + Fetch64(s + 8), 37) * k1;
    y = Rotate(y + v.second + Fetch64(s + 48), 42) * k1;
    x ^= w.second;
    y += v.first + Fetch64(s + 40);
    z = Rotate(z + w.first, 33) * k1;
    v = WeakHashLen32WithSeeds(s, v.second * k1, x + w.first);
    w = WeakHashLen32WithSeeds(s + 32, z + w.second, y + Fetch64(s + 16));
    temp = z;
    z = x;
    x = temp;
    s += 64;
    len -= 64;
  } while (len != 0);
  return HashLen16(HashLen16(v.first, w.first) + ShiftMix(y) * k1 + z,
                   HashLen16(v.second, w.second) + x);
}


__device__ static uint64 HashLen0to16(const char *s, size_t len) {
  if (len > 8) {
    uint64 a = Fetch64(s);
    uint64 b = Fetch64(s + len - 8);
    return HashLen16(a, RotateByAtLeast1(b + len, len)) ^ b;
  }
  if (len >= 4) {
    uint64 a = Fetch32(s);
    return HashLen16(len + (a << 3), Fetch32(s + len - 4));
  }
  if (len > 0) {
    uint8 a = s[0];
    uint8 b = s[len >> 1];
    uint8 c = s[len - 1];
    uint32 y = (uint32)(a) + ((uint32)(b) << 8);
    uint32 z = len + ((uint32)(c) << 2);
    return ShiftMix(y * k2 ^ z * k3) * k2;
  }
  return k2;
}


__device__ static uint64 HashLen17to32(const char *s, size_t len) {
  uint64 a = Fetch64(s) * k1;
  uint64 b = Fetch64(s + 8);
  uint64 c = Fetch64(s + len - 8) * k2;
  uint64 d = Fetch64(s + len - 16) * k0;
  return HashLen16(Rotate(a - b, 43) + Rotate(c, 30) + d,
                   a + Rotate(b ^ k3, 20) - c + len);
}

__device__ static uint64 HashLen33to64(const char *s, size_t len) {
  uint64 z = Fetch64(s + 24);
  uint64 a = Fetch64(s) + (len + Fetch64(s + len - 16)) * k0;
  uint64 b = Rotate(a + z, 52);
  uint64 c = Rotate(a, 37);
  a += Fetch64(s + 8);
  c += Rotate(a, 7);
  a += Fetch64(s + 16);
  uint64 vf = a + z;
  uint64 vs = b + Rotate(a, 31) + c;
  a = Fetch64(s + 16) + Fetch64(s + len - 32);
  z = Fetch64(s + len - 8);
  b = Rotate(a + z, 52);
  c = Rotate(a, 37);
  a += Fetch64(s + len - 24);
  c += Rotate(a, 7);
  a += Fetch64(s + len - 16);
  uint64 wf = a + z;
  uint64 ws = b + Rotate(a, 31) + c;
  uint64 r = ShiftMix((vf + ws) * k2 + (wf + vs) * k0);
  return ShiftMix(r * k0 + vs) * k2;
}

__device__ static uint64 HashLen16(uint64 u, uint64 v) {
  uint128 result;
  result.first = u;
  result.second = v;
  return Hash128to64(result);
}

__device__ static inline uint64 Hash128to64(const uint128 x) {
  // Murmur-inspired hashing.
  const uint64 kMul = 0x9ddfea08eb382d69ULL;
  uint64 a = (Uint128Low64(x) ^ Uint128High64(x)) * kMul;
  a ^= (a >> 47);
  uint64 b = (Uint128High64(x) ^ a) * kMul;
  b ^= (b >> 47);
  b *= kMul;
  return b;
}

__device__ static uint64 Rotate(uint64 val, int shift) {
  // Avoid shifting by 64: doing so yields an undefined result.
  return shift == 0 ? val : ((val >> shift) | (val << (64 - shift)));
}

__device__ static uint64 RotateByAtLeast1(uint64 val, int shift) {
  return (val >> shift) | (val << (64 - shift));
}

__device__ static uint64 ShiftMix(uint64 val) {
  return val ^ (val >> 47);
}

__device__ uint128 WeakHashLen32WithSeeds6(
    uint64 w, uint64 x, uint64 y, uint64 z, uint64 a, uint64 b) {
  a += w;
  b = Rotate(b + a + z, 21);
  uint64 c = a;
  a += x;
  a += y;
  b += Rotate(a, 44);

  uint128 result;
  result.first = (uint64) (a + z);
  result.second = (uint64) (b + c);
  return result;
}

__device__ uint128 WeakHashLen32WithSeeds(
    const char* s, uint64 a, uint64 b) {
  return WeakHashLen32WithSeeds6(Fetch64(s),
                                Fetch64(s + 8),
                                Fetch64(s + 16),
                                Fetch64(s + 24),
                                a,
                                b);
}


/**************** SPOOKY HASH *******************/
#include <stddef.h>
#ifdef HAVE_CONFIG_H
#  include <config.h>
#  ifndef HAVE_ALIGNED_ACCESS_REQUIRED
#     define ALLOW_UNALIGNED_READS 1
#  else
#     define ALLOW_UNALIGNED_READS 0
#  endif
#else
#  if defined(__i386__) || defined(__x86_64__) // add more architectures here
#     define ALLOW_UNALIGNED_READS 1
#  else
#     define ALLOW_UNALIGNED_READS 0
#  endif
#endif /* HAVE_CONFIG_H */

#include <memory.h>

__device__ static inline uint64_t rot64(uint64_t x, int k)
{
	return (x << k) | (x >> (64 - k));
}

__device__ static inline void mix
(
	const uint64_t *data,
	uint64_t *s0, uint64_t *s1, uint64_t *s2,  uint64_t *s3,
	uint64_t *s4, uint64_t *s5, uint64_t *s6,  uint64_t *s7,
	uint64_t *s8, uint64_t *s9, uint64_t *s10, uint64_t *s11
)
{
	*s0 += data[0];		*s2 ^= *s10;	*s11 ^= *s0;	*s0 = rot64(*s0, 11);	*s11 += *s1;
	*s1 += data[1];		*s3 ^= *s11;	*s0 ^= *s1;		*s1 = rot64(*s1, 32);	*s0 += *s2;
	*s2 += data[2];		*s4 ^= *s0;		*s1 ^= *s2;		*s2 = rot64(*s2, 43);	*s1 += *s3;
	*s3 += data[3];		*s5 ^= *s1;		*s2 ^= *s3;		*s3 = rot64(*s3, 31);	*s2 += *s4;
	*s4 += data[4];		*s6 ^= *s2;		*s3 ^= *s4;		*s4 = rot64(*s4, 17);	*s3 += *s5;
	*s5 += data[5];		*s7 ^= *s3;		*s4 ^= *s5;		*s5 = rot64(*s5, 28);	*s4 += *s6;
	*s6 += data[6];		*s8 ^= *s4;		*s5 ^= *s6;		*s6 = rot64(*s6, 39);	*s5 += *s7;
	*s7 += data[7];		*s9 ^= *s5;		*s6 ^= *s7;		*s7 = rot64(*s7, 57);	*s6 += *s8;
	*s8 += data[8];		*s10 ^= *s6;	*s7 ^= *s8;		*s8 = rot64(*s8, 55);	*s7 += *s9;
	*s9 += data[9];		*s11 ^= *s7;	*s8 ^= *s9;		*s9 = rot64(*s9, 54);	*s8 += *s10;
	*s10 += data[10];	*s0 ^= *s8;		*s9 ^= *s10;	*s10 = rot64(*s10, 22);	*s9 += *s11;
	*s11 += data[11];	*s1 ^= *s9;		*s10 ^= *s11;	*s11 = rot64(*s11, 46);	*s10 += *s0;
}

__device__ static inline void endPartial
(
	uint64_t *h0, uint64_t *h1, uint64_t *h2,  uint64_t *h3,
	uint64_t *h4, uint64_t *h5, uint64_t *h6,  uint64_t *h7,
	uint64_t *h8, uint64_t *h9, uint64_t *h10, uint64_t *h11
)
{
	*h11+= *h1;		*h2 ^= *h11;	*h1 = rot64(*h1, 44);
	*h0 += *h2;		*h3 ^= *h0;		*h2 = rot64(*h2, 15);
	*h1 += *h3;		*h4 ^= *h1;		*h3 = rot64(*h3, 34);
	*h2 += *h4;		*h5 ^= *h2;		*h4 = rot64(*h4, 21);
	*h3 += *h5;		*h6 ^= *h3;		*h5 = rot64(*h5, 38);
	*h4 += *h6;		*h7 ^= *h4;		*h6 = rot64(*h6, 33);
	*h5 += *h7;		*h8 ^= *h5;		*h7 = rot64(*h7, 10);
	*h6 += *h8;		*h9 ^= *h6;		*h8 = rot64(*h8, 13);
	*h7 += *h9;		*h10^= *h7;		*h9 = rot64(*h9, 38);
	*h8 += *h10;	*h11^= *h8;		*h10= rot64(*h10, 53);
	*h9 += *h11;	*h0 ^= *h9;		*h11= rot64(*h11, 42);
	*h10+= *h0;		*h1 ^= *h10;	*h0 = rot64(*h0, 54);
}

__device__ static inline void end
(
	uint64_t *h0,	uint64_t *h1,	uint64_t *h2,	uint64_t *h3,
	uint64_t *h4,	uint64_t *h5,	uint64_t *h6,	uint64_t *h7,
	uint64_t *h8,	uint64_t *h9,	uint64_t *h10,	uint64_t *h11
)
{
	endPartial(h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11);
	endPartial(h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11);
	endPartial(h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11);
}

__device__ static inline void short_mix
(
	uint64_t *h0,
	uint64_t *h1,
	uint64_t *h2,
	uint64_t *h3
)
{
	*h2 = rot64(*h2, 50);	*h2 += *h3;  *h0 ^= *h2;
	*h3 = rot64(*h3, 52);	*h3 += *h0;  *h1 ^= *h3;
	*h0 = rot64(*h0, 30);	*h0 += *h1;  *h2 ^= *h0;
	*h1 = rot64(*h1, 41);	*h1 += *h2;  *h3 ^= *h1;
	*h2 = rot64(*h2, 54);	*h2 += *h3;  *h0 ^= *h2;
	*h3 = rot64(*h3, 48);	*h3 += *h0;  *h1 ^= *h3;
	*h0 = rot64(*h0, 38);	*h0 += *h1;  *h2 ^= *h0;
	*h1 = rot64(*h1, 37);	*h1 += *h2;  *h3 ^= *h1;
	*h2 = rot64(*h2, 62);	*h2 += *h3;  *h0 ^= *h2;
	*h3 = rot64(*h3, 34);	*h3 += *h0;  *h1 ^= *h3;
	*h0 = rot64(*h0, 5);	*h0 += *h1;  *h2 ^= *h0;
	*h1 = rot64(*h1, 36);	*h1 += *h2;  *h3 ^= *h1;
}

__device__ static inline void short_end
(
	uint64_t *h0,
	uint64_t *h1,
	uint64_t *h2,
	uint64_t *h3
)
{
	*h3 ^= *h2;  *h2 = rot64(*h2, 15);  *h3 += *h2;
	*h0 ^= *h3;  *h3 = rot64(*h3, 52);  *h0 += *h3;
	*h1 ^= *h0;  *h0 = rot64(*h0, 26);  *h1 += *h0;
	*h2 ^= *h1;  *h1 = rot64(*h1, 51);  *h2 += *h1;
	*h3 ^= *h2;  *h2 = rot64(*h2, 28);  *h3 += *h2;
	*h0 ^= *h3;  *h3 = rot64(*h3, 9);   *h0 += *h3;
	*h1 ^= *h0;  *h0 = rot64(*h0, 47);  *h1 += *h0;
	*h2 ^= *h1;  *h1 = rot64(*h1, 54);  *h2 += *h1;
	*h3 ^= *h2;  *h2 = rot64(*h2, 32);  *h3 += *h2;
	*h0 ^= *h3;  *h3 = rot64(*h3, 25);  *h0 += *h3;
	*h1 ^= *h0;  *h0 = rot64(*h0, 63);  *h1 += *h0;
}

__device__ void spooky_shorthash
(
	const void *message,
	size_t length,
	uint64_t *hash1,
	uint64_t *hash2
)
{
	uint64_t buf[2 * SC_NUMVARS];
	union
	{
		const uint8_t *p8;
		uint32_t *p32;
		uint64_t *p64;
		size_t i;
	} u;
	size_t remainder;
	uint64_t a, b, c, d;
	u.p8 = (const uint8_t *)message;

	if (!ALLOW_UNALIGNED_READS && (u.i & 0x7))
	{
		memcpy(buf, message, length);
		u.p64 = buf;
	}

	remainder = length % 32;
	a = *hash1;
	b = *hash2;
	c = SC_CONST;
	d = SC_CONST;

	if (length > 15)
	{
		const uint64_t *endp = u.p64 + (length/32)*4;

		// handle all complete sets of 32 bytes
		for (; u.p64 < endp; u.p64 += 4)
		{
			c += u.p64[0];
			d += u.p64[1];
			short_mix(&a, &b, &c, &d);
			a += u.p64[2];
			b += u.p64[3];
		}

		// Handle the case of 16+ remaining bytes.
		if (remainder >= 16)
		{
			c += u.p64[0];
			d += u.p64[1];
			short_mix(&a, &b, &c, &d);
			u.p64 += 2;
			remainder -= 16;
		}
	}

	// Handle the last 0..15 bytes, and its length
	d = ((uint64_t)length) << 56;
	switch (remainder)
	{
		case 15:
			d += ((uint64_t)u.p8[14]) << 48;
		case 14:
			d += ((uint64_t)u.p8[13]) << 40;
		case 13:
			d += ((uint64_t)u.p8[12]) << 32;
		case 12:
			d += u.p32[2];
			c += u.p64[0];
			break;
		case 11:
			d += ((uint64_t)u.p8[10]) << 16;
		case 10:
			d += ((uint64_t)u.p8[9]) << 8;
		case 9:
			d += (uint64_t)u.p8[8];
		case 8:
			c += u.p64[0];
			break;
		case 7:
			c += ((uint64_t)u.p8[6]) << 48;
		case 6:
			c += ((uint64_t)u.p8[5]) << 40;
		case 5:
			c += ((uint64_t)u.p8[4]) << 32;
		case 4:
			c += u.p32[0];
			break;
		case 3:
			c += ((uint64_t)u.p8[2]) << 16;
		case 2:
			c += ((uint64_t)u.p8[1]) << 8;
		case 1:
			c += (uint64_t)u.p8[0];
			break;
		case 0:
			c += SC_CONST;
			d += SC_CONST;
	}
	short_end(&a, &b, &c, &d);
	*hash1 = a;
	*hash2 = b;
}

__device__ void spooky_init
(
	struct spooky_state *state,
	uint64_t seed1,
	uint64_t seed2
)
{
	state->m_length = 0;
	state->m_remainder = 0;
	state->m_state[0] = seed1;
	state->m_state[1] = seed2;
}

__device__ void spooky_update
(
	struct spooky_state *state,
	const void *message,
	size_t length
)
{
	uint64_t h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11;
	size_t newLength = length + state->m_remainder;
	uint8_t remainder;
	union
	{
		const uint8_t *p8;
		uint64_t *p64;
		size_t i;
	} u;
	const uint64_t *endp;

	// Is this message fragment too short?  If it is, stuff it away.
	if (newLength < SC_BUFSIZE)
	{
		memcpy(&((uint8_t *)state->m_data)[state->m_remainder], message, length);
		state->m_length = length + state->m_length;
		state->m_remainder = (uint8_t)newLength;
		return;
	}

	// init the variables
	if (state->m_length < SC_BUFSIZE)
	{
		h0 = h3 = h6 = h9  = state->m_state[0];
		h1 = h4 = h7 = h10 = state->m_state[1];
		h2 = h5 = h8 = h11 = SC_CONST;
	}
	else
	{
		h0 = state->m_state[0];
		h1 = state->m_state[1];
		h2 = state->m_state[2];
		h3 = state->m_state[3];
		h4 = state->m_state[4];
		h5 = state->m_state[5];
		h6 = state->m_state[6];
		h7 = state->m_state[7];
		h8 = state->m_state[8];
		h9 = state->m_state[9];
		h10 = state->m_state[10];
		h11 = state->m_state[11];
	}
	state->m_length = length + state->m_length;

	// if we've got anything stuffed away, use it now
	if (state->m_remainder)
	{
		uint8_t prefix = SC_BUFSIZE-state->m_remainder;
		memcpy(&(((uint8_t *)state->m_data)[state->m_remainder]), message, prefix);
		u.p64 = state->m_data;
		mix(u.p64, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);
		mix(&u.p64[SC_NUMVARS], &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);
		u.p8 = ((const uint8_t *)message) + prefix;
		length -= prefix;
	}
	else
	{
		u.p8 = (const uint8_t *)message;
	}

	// handle all whole blocks of SC_BLOCKSIZE bytes
	endp = u.p64 + (length/SC_BLOCKSIZE)*SC_NUMVARS;
	remainder = (uint8_t)(length-((const uint8_t *)endp - u.p8));
	if (ALLOW_UNALIGNED_READS || (u.i & 0x7) == 0)
	{
		while (u.p64 < endp)
		{
			mix(u.p64, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);
			u.p64 += SC_NUMVARS;
		}
	}
	else
	{
		while (u.p64 < endp)
		{
			memcpy(state->m_data, u.p8, SC_BLOCKSIZE);
			mix(state->m_data, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);
			u.p64 += SC_NUMVARS;
		}
	}

	// stuff away the last few bytes
	state->m_remainder = remainder;
	memcpy(state->m_data, endp, remainder);

	// stuff away the variables
	state->m_state[0] = h0;
	state->m_state[1] = h1;
	state->m_state[2] = h2;
	state->m_state[3] = h3;
	state->m_state[4] = h4;
	state->m_state[5] = h5;
	state->m_state[6] = h6;
	state->m_state[7] = h7;
	state->m_state[8] = h8;
	state->m_state[9] = h9;
	state->m_state[10] = h10;
	state->m_state[11] = h11;
}

__device__ void spooky_final
(
	struct spooky_state *state,
	uint64_t *hash1,
	uint64_t *hash2
)
{
	uint64_t h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11;
	const uint64_t *data = (const uint64_t *)state->m_data;
	uint8_t remainder = state->m_remainder;

	// init the variables
	if (state->m_length < SC_BUFSIZE)
	{
		spooky_shorthash(state->m_data, state->m_length, hash1, hash2);
		return;
	}

	h0 = state->m_state[0];
	h1 = state->m_state[1];
	h2 = state->m_state[2];
	h3 = state->m_state[3];
	h4 = state->m_state[4];
	h5 = state->m_state[5];
	h6 = state->m_state[6];
	h7 = state->m_state[7];
	h8 = state->m_state[8];
	h9 = state->m_state[9];
	h10 = state->m_state[10];
	h11 = state->m_state[11];

	if (remainder >= SC_BLOCKSIZE)
	{
		// m_data can contain two blocks; handle any whole first block
		mix(data, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);
		data += SC_NUMVARS;
		remainder -= SC_BLOCKSIZE;
	}

	// mix in the last partial block, and the length mod SC_BLOCKSIZE
	memset(&((uint8_t *)data)[remainder], 0, (SC_BLOCKSIZE-remainder));

	((uint8_t *)data)[SC_BLOCKSIZE-1] = remainder;
	mix(data, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);

	// do some final mixing
	end(&h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);

	*hash1 = h0;
	*hash2 = h1;
}

__device__ uint64_t spooky_hash64
(
	const void *message,
	size_t length,
	uint64_t seed
)
{
	uint64_t hash1 = seed;
	spooky_hash128(message, length, &hash1, &seed);
	return hash1;
}

__device__ void spooky_hash128
(
	const void *message,
	size_t length,
	uint64_t *hash1,
	uint64_t *hash2
)
{
	uint64_t h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11;
	uint64_t buf[SC_NUMVARS];
	uint64_t *endp;
	union
	{
		const uint8_t *p8;
		uint64_t *p64;
		uintptr_t i;
	} u;
	size_t remainder;

	if (length < SC_BUFSIZE)
	{
		spooky_shorthash(message, length, hash1, hash2);
		return;
	}

	h0 = h3 = h6 = h9  = *hash1;
	h1 = h4 = h7 = h10 = *hash2;
	h2 = h5 = h8 = h11 = SC_CONST;

	u.p8 = (const uint8_t *)message;
	endp = u.p64 + (length/SC_BLOCKSIZE)*SC_NUMVARS;

	// handle all whole blocks of SC_BLOCKSIZE bytes
	if (ALLOW_UNALIGNED_READS || (u.i & 0x7) == 0)
	{
		while (u.p64 < endp)
		{
			mix(u.p64, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);
			u.p64 += SC_NUMVARS;
		}
	}
	else
	{
		while (u.p64 < endp)
		{
			memcpy(buf, u.p64, SC_BLOCKSIZE);
			mix(buf, &h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);
			u.p64 += SC_NUMVARS;
		}
	}

	// handle the last partial block of SC_BLOCKSIZE bytes
	remainder = (length - ((const uint8_t *)endp-(const uint8_t *)message));
	memcpy(buf, endp, remainder);
	memset(((uint8_t *)buf)+remainder, 0, SC_BLOCKSIZE-remainder);
	((uint8_t *)buf)[SC_BLOCKSIZE-1] = remainder;
	mix(buf, &h0 , &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);

	// do some final mixing
	end(&h0, &h1, &h2, &h3, &h4, &h5, &h6, &h7, &h8, &h9, &h10, &h11);
	*hash1 = h0;
	*hash2 = h1;
}
