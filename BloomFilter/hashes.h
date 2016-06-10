#ifdef ONLYC
#define __device__
#endif

__device__ unsigned long djb2_hash(unsigned char *str);

__device__ uint64_t MurmurHash64A (const void * key, int len, unsigned int seed);

__device__ uint64_t APHash(char* str, unsigned int length);

__device__ uint64_t fnvhash(char * string);

__device__ uint64_t SDBMHash(char* str, unsigned int length);

__device__ uint64_t RSHash(char* str, unsigned int length);
