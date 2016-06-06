#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "hashes.h"
#include "spooky.h"
#include "city.h"

#define BUF 255			//length of the filter cells on file
#define L 11			//length of the reads
#define MSEED 7127		//static seed for murmur function
#define SSEED 5449		//static seed for spooky function

int CheckBit(uint64_t* filter, unsigned long i, unsigned long n);
int CheckHash(uint64_t* filter, char* read, unsigned long n);

int main(int argc, char* argv[]){

	uint64_t* bloom;
	unsigned long cell;
	FILE* bf, * fp;
	unsigned long n;
	char* seq;
	char cons;
	
	n = atoi(argv[2]);

	if(!(bloom = (uint64_t*)malloc(n*sizeof(uint64_t)))){		
		fprintf(stdout, "Error: Not enough memory\n");
		exit(1);
	}

	if(!(bf = fopen(argv[1], "r"))){
		fprintf(stdout, "Error: File not found\n");
		exit(1);
	}
	fscanf(bf, "%lu", &cell);

	if(feof(bf)){
		fprintf(stdout, "Warning: Empty file\n");
		exit(1);
	}

	//Load bloom filter
	for(int i=0; i<n && !feof(bf); i++){
		bloom[i] = cell;
		fscanf(bf, "%lu", &cell);
	}
	fclose(bf);
	
	/*DEBUG
	for(int k=0; k<n; k++)
		fprintf(stdout, "%lu\n", bloom[k]);
	*/
	
	//Checking part

	//with file
	if(!(fp = fopen(argv[3], "r"))){
		fprintf(stdout, "Error: File not found\n");
		exit(1);
	}
	fgets(seq, L, fp);
	if(feof(fp)){
		fprintf(stdout, "Warning: Empty file\n");
		exit(1);
	}
	while(!feof(fp)){
		cons = fgetc(fp);

		if(CheckHash(bloom, seq, n))
			printf("%s belongs to the spectrum\n", seq);
		else
			printf("%s not in spectrum\n", seq);
		
		fgets(seq, L, fp);
	}
	fclose(fp);
	
	free(bloom);
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
