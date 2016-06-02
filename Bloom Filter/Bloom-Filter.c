#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef SILENT
#define printf(...)
#endif

void SetBit(uint64_t* filter, int i);
void HashRead(uint64_t* filter, char* read);
unsigned long djb2_hash(unsigned char *str);


int main (int argc, char* argv[]){		//File name (read), spectrum dimension, length of reads, file name (write)

	FILE* fp, bf;
	uint64_t* bloom;
	int  i;								//Filter index
	char* seq; 						//Read sequence
	
	
	//Allocate memory for Bloom Filter
	if(!(bloom = (uint64_t*)malloc(argv[2]*sizeof(uint64_t))){		//assigns 64*n bits
		fprintf(stdout, "Error: Not enough memory\n");
		exit(1);
	}
	bloom = 0;
	
	//Try file
	if(!(fp = fopen(argv[1], r))){
		fprintf(stdout, "Error: File not found\n");
		exit(1);
	}
	fgets(seq, argv[3], fp);
	if(seq == EOF){
		fprintf(stdoud, "Warning: Empty file\n");
		exit(1);
	}
	
	//Fill up the filter
	while(seq != EOF){
		HashRead(bloom, seq, argv[2]);
		
		fgets(seq, argv[3], fp);
	}
	fclose(fp);
	
	bf = fopen(argv[2], w+);
	fputs(bloom, bf);
	fclose(bf);
	
	return 0;
}

void SetBit(uint64_t* filter, unsigned long i, unsigned long n){
	unsigned long k = i/(64*n);
	unsigned long pos = i%(64*n);
	
	filter[k] |= 1 << pos;
}

void HashRead(uint64_t* filter, char* read, unsigned long n){
	unsigned long i;
	
	i = djb2_hash((unsigned char)read);
	SetBit(filter, i, argv[2]);
}

//Simple hash
unsigned long djb2_hash(unsigned char *str){
    unsigned long hash = 5381;
    int c;

    while (c = *str++)
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}