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

#define MSEED 7127		//static seed for murmur function
#define SSEED 5449		//static seed for spooky function

void SetBit(uint64_t* filter, unsigned long i, unsigned long n);
void HashRead(uint64_t* filter, char* read, unsigned long n);
int CheckBit(uint64_t* filter, unsigned long i, unsigned long n);
int CheckHash(uint64_t* filter, char* read, unsigned long n);


int main (int argc, char* argv[]){		//File name (read), spectrum dimension, length of reads, file name (write)

	FILE* fp, * bf;
	uint64_t* bloom;				//pointer to filter
	unsigned long n;				//spectrum dimension 
	int l; 						//length of read
	char* seq; 					//Read sequence
	char cons;					//consume \n
	
	n = atoi(argv[2]);
	l = atoi(argv[3])+1;

	//Allocate memory for Bloom Filter
	if(!(bloom = (uint64_t*)malloc(n*sizeof(uint64_t)))){		//assigns 64*n bits
		fprintf(stdout, "Error: Not enough memory\n");
		exit(1);
	}
	seq = (char*)malloc(10*sizeof(char));
	for(int i=0; i<n; i++)
		bloom[i]=0;
	
	//Try file
	if(!(fp = fopen(argv[1], "r"))){
		fprintf(stdout, "Error: File not found\n");
		exit(1);
	}
	fgets(seq, l, fp);
	//printf("%s\n", seq);
	if(feof(fp)){
		fprintf(stdout, "Warning: Empty file\n");
		exit(1);
	}
	
	//Fill up the filter
	while(!feof(fp)){
		HashRead(bloom, seq, n);
		cons = fgetc(fp);
		fgets(seq, l, fp);
		//printf("%s\n", seq);
	}
	fclose(fp);
	

	//Check hash
	if(!(fp = fopen(argv[4], "r"))){
		fprintf(stdout, "Error: File not found\n");
		exit(1);
	}
	fgets(seq, l, fp);
	cons = fgetc(fp);
	if(feof(fp)){
		fprintf(stdout, "Warning: Empty file\n");
		exit(1);
	}
	while(!feof(fp)){

		if(CheckHash(bloom, seq, n))
			printf("%s belongs to the spectrum\n", seq);
		else
			printf("%s not in spectrum\n", seq);
		
		fgets(seq, l, fp);
		cons = fgetc(fp);
	}
	fclose(fp);
	
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

//Check that the bit is in the filter
int CheckBit(uint64_t* filter, unsigned long i, unsigned long n){
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
int CheckHash(uint64_t* filter, char* read, unsigned long n){
	int flag = 1;
	unsigned long i;

	for(int k=1; k<=8 && flag; k++){
		
		switch (k){
		
		case 1:
			i = djb2_hash((unsigned char*)read);
			break;
		case 2:
			i = MurmurHash64A(read, strlen(read), MSEED);
			break;
		case 3:
			i = APHash(read, strlen(read));
			break;
		case 4:
			i = CityHash64(read, strlen(read));
			break;
		case 5:
			i = spooky_hash64(read, strlen(read), SSEED);
			break;
		case 6:
			i = fnvhash(read);
			break;
		case 7:
			i = SDBMHash(read, strlen(read));
			break;
		case 8:
			i = RSHash(read, strlen(read));
			break;
		
		default: break;
		}
		
		flag = CheckBit(filter, i, n);
	}

	return flag;
}
