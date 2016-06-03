#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef SILENT
#define printf(...)
#endif

#define SEED 75489

void SetBit(uint64_t* filter, unsigned long i, unsigned long n);
void HashRead(uint64_t* filter, char* read, unsigned long n);
unsigned long djb2_hash(unsigned char *str);
uint64_t MurmurHash64A (const void * key, int len, unsigned int seed);


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
		fprintf(bf, "%lu", bloom[k]);
	fclose(bf);
	
	free(bloom);
	
	return 0;
}

//Sets bit into bloom filter
void SetBit(uint64_t* filter, unsigned long i, unsigned long n){
	unsigned long k = i % n;
	unsigned long pos = i % 64;
	
	filter[k] |= 1 << pos;
}

//Takes the read and feeds it to hash functions
void HashRead(uint64_t* filter, char* read, unsigned long n){
	unsigned long i;
	
	i = djb2_hash((unsigned char*)read);
	SetBit(filter, i, n);
	
	printf("djb2: %lu ", i);

	i = MurmurHash64A(read, strlen(read), SEED);
	SetBit(filter, i, n);

	printf("Murmurhash: %lu\n", i);
}

//Simple hash
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
