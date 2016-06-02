#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef SILENT
#define printf(...)
#endif

void SetBit(uint64_t* filter, unsigned long i, unsigned long n);
void HashRead(uint64_t* filter, char* read, unsigned long n);
unsigned long djb2_hash(unsigned char *str);


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
	bloom = 0;
	
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
	fprintf(bf, "%s", (char*)bloom);
	fclose(bf);
	
	free(bloom);
	
	return 0;
}

//Sets bit into bloom filter
void SetBit(uint64_t* filter, unsigned long i, unsigned long n){
	unsigned long k = i/(64*n);
	unsigned long pos = i%(64*n);
	
	filter[k] |= 1 << pos;
}

//Takes the read and feeds it to hash functions
void HashRead(uint64_t* filter, char* read, unsigned long n){
	unsigned long i;
	
	i = djb2_hash((unsigned char*)read);
	SetBit(filter, i, n);
}

//Simple hash
unsigned long djb2_hash(unsigned char *str){
    unsigned long hash = 5381;
    int c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}
