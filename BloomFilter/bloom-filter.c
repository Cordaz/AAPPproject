#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "hashes.h"
#include "city.h"
#include "spooky.h"

#ifdef SILENT
#define printf(...)
#endif

#define L 10			//length of reads
#define MSEED 7127		//static seed for murmur function
#define SSEED 5449		//static seed for spooky function

void SetBit(uint64_t* filter, unsigned long i, unsigned long n);
void HashRead(uint64_t* filter, char* read, unsigned long n);


int main (int argc, char* argv[]){		//File name (read), file name (write), spectrum dimension

	FILE* fp, * bf;
	uint64_t* bloom;				//pointer to filter
	unsigned long n;				//spectrum dimension 
	char* seq; 					//Read sequence
	char cons;					//consume \n	

	n = atoi(argv[3]);

	//Allocate memory for Bloom Filter
	if(!(bloom = (uint64_t*)malloc(n*sizeof(uint64_t)))){		//assigns 64*n bits
		fprintf(stdout, "Error: Not enough memory\n");
		exit(1);
	}
	for(int i=0; i<n; i++)
		bloom[i]=0;
	seq = (char*)malloc(L*sizeof(char));	

	//Try file
	if(!(fp = fopen(argv[1], "r"))){
		fprintf(stdout, "Error: File not found\n");
		exit(1);
	}
	fgets(seq, L+1, fp);
	//printf("%s\n", seq);
	if(feof(fp)){
		fprintf(stdout, "Warning: Empty file\n");
		exit(1);
	}
	
	//Fill up the filter
	while(!feof(fp)){
		HashRead(bloom, seq, n);
		cons = fgetc(fp);
		fgets(seq, L+1, fp);
		//printf("%s\n", seq);
	}
	fclose(fp);
	
	bf = fopen(argv[2], "w+");
	for(int k=0; k<n; k++){
		fprintf(bf, "%lu\n", bloom[k]);
	}	
	fclose(bf);
	
	free(bloom);
	
	return 0;
}

//Sets bit into bloom filter
void SetBit(uint64_t* filter, unsigned long i, unsigned long n){
	unsigned long k = i % n;
	unsigned long pos = i % 64;
	uint64_t bit = 1;
	
	bit = bit << pos;
	//printf(", cell: %lu, bit: %lu\n", k, pos);
	
	filter[k] = filter[k] | bit;
}

//Takes the read and feeds it to hash functions
void HashRead(uint64_t* filter, char* read, unsigned long n){
	unsigned long i;
	
	i = djb2_hash((unsigned char*)read);
	//printf("djb2: %lu ", i);
	SetBit(filter, i, n);

	i = MurmurHash64A(read, strlen(read), MSEED);
	//printf("Murmurhash: %lu", i);
	SetBit(filter, i, n);
	
	i = APHash(read, strlen(read));
	//printf("APHash: %lu", i);
	SetBit(filter, i, n);
	
	i = CityHash64(read, strlen(read));
	//printf("CityHash: %lu", i);
	SetBit(filter, i, n);
	
	i = spooky_hash64(read, strlen(read), SSEED);
	//printf("SpookyHash: %lu", i);
	SetBit(filter, i, n);
	
	i = fnvhash(read);
	//printf("FNVhash: %lu", i);
	SetBit(filter, i, n);
	
	i = SDBMHash(read, strlen(read));
	//printf("SDBMhash: %lu", i);
	SetBit(filter, i, n);
	
	i = RSHash(read, strlen(read));
	//printf("RShash: %lu", i);
	SetBit(filter, i, n);
	
}
