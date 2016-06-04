#include <stdint.h>
#include "hashes.h"

const uint64_t FNV_PRIME    = 1099511628211;				
const uint64_t OFFSET_BASIS = 14695981039346656037UL;

//djb2 hash
unsigned long djb2_hash(unsigned char *str){
    unsigned long hash = 5381;
    int c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

//MurmurHash
uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
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
uint64_t APHash(char* str, unsigned int length) {
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
uint64_t fnvhash(char * string){
   uint64_t hash = OFFSET_BASIS;
   for(uint8_t * c = (uint8_t*)string; *c != 0; ++c){
      hash ^= *c;
      hash *= FNV_PRIME;
   }
   return hash;
}

//SDBM hash
uint64_t SDBMHash(char* str, unsigned int length) {
	uint64_t hash = 0;
	unsigned int i = 0;

	for (i = 0; i < length; str++, i++)
	{
		hash = (*str) + (hash << 6) + (hash << 16) - hash;
	}

	return hash;
}

//RS hash
uint64_t RSHash(char* str, unsigned int length) {
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
