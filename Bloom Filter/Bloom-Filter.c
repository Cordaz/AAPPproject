#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef SILENT
#define printf(...)
#endif

void SetBit(char* filter, int i);


int main (int argc, char* argv[]){		//File name (read), spectrum dimension, length of reads, file name (write)

	FILE* fp, bf;
	uint64_t* bloom;
	int  i;								//Filter index
	STRING seq; 						//Read sequence
	
	
	//Allocate memory for Bloom Filter
	bloom = (uint64_t*)malloc(argv[2]*sizeof(uint64_t));		//assigns 64*n bits, 8*n bytes
	bloom = 0;
		
	fp = fopen(argv[1], r);
	
	if(fp == NULL){
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
		//feed to hash
		
		//bit manipulation on filter
		
		fgets(seq, argv[3], fp);
	}
	fclose(fp);
	
	bf = fopen(argv[2], w+);
	fputs(bloom, bf);
	fclose(bf);
	
	return 0;
}

void SetBit(uint64_t* filter, int i, int n){
	int k = i/(64*n);
	int pos = i%(64*n);
	
	filter[k] |= 1 << pos;
}