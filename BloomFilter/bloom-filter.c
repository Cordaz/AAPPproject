#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "city.h"
#include "spooky.h"

#ifdef SILENT
#define printf(...)
#endif

#define SEED 75489		//define another one for spooky
const uint64_t FNV_PRIME    = 1099511628211;				
const uint64_t OFFSET_BASIS = 14695981039346656037UL;

void SetBit(uint64_t* filter, unsigned long i, unsigned long n);
void HashRead(uint64_t* filter, char* read, unsigned long n);

unsigned long djb2_hash(unsigned char *str);
uint64_t MurmurHash64A (const void * key, int len, unsigned int seed);
uint64_t APHash(char* str, unsigned int length)
uint64_t fnvhash(char * string);
uint64_t SDBMHash(char* str, unsigned int length);
uint64_t RSHash(char* str, unsigned int length);


int main (int argc, char* argv[]){		//File name (read), spectrum dimension, length of reads, file name (write)

	FILE* fp, * bf;
	uint64_t* bloom;				//pointer to filter
	unsigned long n;				//spectrum dimension 
	int l; 						//length of read
	char* seq; 					//Read sequence
	
	n = atoi(argv[2]);
	l = atoi(argv[3]);

	//Allocate memory for Bloom Filter
	if(!(bloom = (uint64_t*)malloc(n*sizeof(uint64_t)))){		//assigns 64*n bits
		fprintf(stdout, "Error: Not enough memory\n");
		exit(1);
	}
	for(int i=0; i<n; i++)
		bloom[i]=0;
	
	//Try file
	if(!(fp = fopen(argv[1], "r"))){
		fprintf(stdout, "Error: File not found\n");
		exit(1);
	}
	fgets(seq, l, fp);
	if(feof(fp)){
		fprintf(stdout, "Warning: Empty file\n");
		exit(1);
	}
	
	//Fill up the filter
	while(!feof(fp)){
		HashRead(bloom, seq, n);
		
		fgets(seq, l, fp);
	}
	fclose(fp);
	
	bf = fopen(argv[4], "w+");
	for(int k=0; k<n; k++)
		fprintf(bf, "%lu\n", bloom[k]);
	fclose(bf);
	
	free(bloom);
	
	return 0;
}

//Sets bit into bloom filter
void SetBit(uint64_t* filter, unsigned long i, unsigned long n){
	unsigned long k = i % n;
	unsigned long pos = i % 64;
	
	printf(", cell: %lu, bit: %lu\n", k, pos);
	
	filter[k] |= 1 << pos;
}

//Takes the read and feeds it to hash functions
void HashRead(uint64_t* filter, char* read, unsigned long n){
	unsigned long i;
	
	i = djb2_hash((unsigned char*)read);
	printf("djb2: %lu ", i);
	SetBit(filter, i, n);

	i = MurmurHash64A(read, strlen(read), SEED);
	printf("Murmurhash: %lu", i);
	SetBit(filter, i, n);
	
	i = hash64shift(read);
	printf("hash64: %lu", i);
	SetBit(filter, i, n);
	
	i = CityHash64(read, strlen(read));
	printf("CityHash: %lu", i);
	SetBit(filter, i, n);
	
	i = spooky_hash64(read, strlen(read), SEED);
	printf("SpookyHash: %lu", i);
	SetBit(filter, i, n);
	
	i = fnvhash(read);
	printf("FNVhash: %lu", i);
	SetBit(filter, i, n);
	
	i = SDBMHash(read, strlen(read));
	printf("SDBMhash: %lu", i);
	SetBit(filter, i, n);
	
	i = RSHash(read, strlen(read));
	printf("RShash: %lu", i);
	SetBit(filter, i, n);
	
}

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