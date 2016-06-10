#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "BloomFilter/hashes.h"
#include "BloomFilter/city.h"
#include "BloomFilter/spooky.h"

#define PARAM 3 //input file (*.seq), output file (bloom filter), output file (read list)
#define READS_LENGTH 35
#define L 10 //l-tuple lenght
#define M 6 //multiplicity
#define SPECTRUM_MAX_SIZE 1048576 //BASES^L
#define BUFFER 50

#define A 0 //00
#define C 1 //01
#define G 2 //10
#define T 3 //11

#define INTSIZE 16
#define INT64SIZE 64
#define BASES 4

#define MSEED 7127		//static seed for murmur function
#define SSEED 5449		//static seed for spooky function

#ifdef SILENT
#define printf(...)
#endif

typedef struct list_s {
   char str[READS_LENGTH + 1];
   struct list_s * next;
   } list_t;

/***************** PROTOTYPES ********************/
int check_read(char *);
void get_tuple(char *, char *, int);

void increment_counter(uint8_t *, char *);
void extract_tuple(unsigned int, char *);

uint64_t * bloom_filter (char **, unsigned int);
void SetBit(uint64_t *, int, unsigned int);
void HashRead(uint64_t *, char *, unsigned int);

void rm_newline(char *);
void list_to_array(list_t *, char **);
list_t * add_elem(list_t *, char *);
void free_list(list_t *);

/******************** MAIN ***********************/
int main (int argc, char * argv[]) {
   if(argc < PARAM + 1) {
      fprintf(stdout, "Error: parameters\n");
      exit(1);
   }
   
   unsigned int i;
   
   list_t * reads_list = NULL;
   
   FILE * in;
   char buf[BUFFER+1];
   char tuple[L+1];
   unsigned int dim_reads=0;
   uint8_t counter[SPECTRUM_MAX_SIZE];
   for(i=0; i<SPECTRUM_MAX_SIZE; i++)
      counter[i] = 0;
   
   in = fopen(argv[1], "r");
   
   fgets(buf, BUFFER, in);
   while(!feof(in)) {
      if(buf[0] != '>') {
         rm_newline(buf);
         if(check_read(buf)) {
            //printf("Read: %s, tuple:\n", buf);
            //Store read
            reads_list = add_elem(reads_list, buf);
            dim_reads++;
            
            //Check tuple
            for(i=0; i< (READS_LENGTH - L + 1); i++) {
               get_tuple(buf, tuple, i);
               //printf("\t%s\n", tuple);
               //Increment counter of tuple occurences
               increment_counter(counter, tuple);
            }
         }
      }
      fgets(buf, BUFFER, in);
   }
   
   fclose(in);
   
   printf("Input file read\n");
   
   //Store reads in an array
   char ** reads;
   if(!(reads = (char **)malloc(sizeof(char *) * dim_reads))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   for(i=0; i<dim_reads; i++) {
      if(!(reads[i] = (char *)malloc(sizeof(char) * (READS_LENGTH+1)))) {
         fprintf(stdout, "Error allocation\n");
         exit(1);
      }
   }
   
   list_to_array(reads_list, reads);
   printf("Freeing reads_list\n");
   free_list(reads_list);
   printf("Freed reads_list\n");
   
   
   //Compute spectrum (Extract freqeunt tuples from counter)
   list_t * tuples_list = NULL;
   unsigned int spectrum_size=0;
   
   for(i=0; i<SPECTRUM_MAX_SIZE; i++) {
      if(counter[i] >= M) {
         //It's frequent, extract
         extract_tuple(i, tuple);
         tuples_list = add_elem(tuples_list, tuple);
         spectrum_size++;
      }
   }
   printf("Spectrum computed\n");
   
   char ** tuples;
   if(!(tuples = (char **)malloc(sizeof(char *) * spectrum_size))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   for(i=0; i<spectrum_size; i++) {
      if(!(tuples[i] = (char *)malloc(sizeof(char) * (L+1)))) {
         fprintf(stdout, "Error allocation\n");
         exit(1);
      }
   }
   list_to_array(tuples_list, tuples);
   printf("Freeing tuples_list\n");
   free_list(tuples_list);
   printf("Freed tuples_list\n");
   
   /*
   printf("DEBUG: printing array of tuples (dim: %u)\n", spectrum_size);
   /*
   for(i=0;i<spectrum_size;i++) {
      printf("%s\n", tuples[i]);
   }
   */
   
   //Apply bloom filter
   printf("Applying bloom filter\n");
   uint64_t * bloom = bloom_filter(tuples, spectrum_size);
   
   
   
   //Write out
   FILE * out;
   //Bloom filter
   out = fopen(argv[2], "w+");
   
   //Print size
   fprintf(out, "%u\n", spectrum_size);
   
   //Print data
   for(i=0; i<spectrum_size; i++) {
      fprintf(out, "%lu\n", bloom[i]);
   }
   
   fclose(out);
   //Reads
   out = fopen(argv[3], "w+");
   
   //Print size
   fprintf(out, "%u\n", dim_reads);
   
   //Print data
   for(i=0; i<dim_reads; i++) {
      fprintf(out, "%s\n", reads[i]);
   }
   
   fclose(out);
   
   for(i=0;i<dim_reads;i++)
      free(reads[i]);
   free(reads);
   for(i=0;i<spectrum_size;i++)
      free(tuples[i]);
   free(tuples);
}

/**************** PREPROCESSING ******************/

int check_read(char * read) {
   int i;
   for(i=0; i<READS_LENGTH+1; i++) {
      if(read[i] == 'N')
         return 0;
   }
   return 1;
}

void get_tuple(char * read, char * tuple, int start) {
   int i;
   for(i=0; i<L; i++) {
      tuple[i] = read[start+i];
   }
   tuple[L] = '\0';
}

/***************** SPECTRUM **********************/

void increment_counter(uint8_t * counter, char * tuple) {
   int j;
   unsigned int arrayIndex=0;
   
   for(j=0; j<L; j++) {
      arrayIndex = arrayIndex << 2; //Move bit to left
      switch(tuple[j]) {
         case 'A':
            //Should add A (0)
            break;
         case 'C':
            arrayIndex += C;
            break;
         case 'G':
            arrayIndex += G;
            break;
         case 'T':
            arrayIndex += T;
            break;
      }
   }
   
   //Counter update procedure
   //Check if necessary (< M)
   if(counter[arrayIndex] < M) {
      counter[arrayIndex]++;
   }
}

void extract_tuple(unsigned int coded, char * tuple) {
   unsigned int mask;
   unsigned int tmp;
   int i;
   
   for(i=0;i<L;i++) {
      //Create mask
      mask = 3; //11 binary
      mask = mask << (i*2);
      tmp = coded & mask;
      tmp = tmp >> (i*2);
      switch(tmp) {
         case A:
            tuple[L-1-i] = 'A';
            break;
         case C:
            tuple[L-1-i] = 'C';
            break;
         case G:
            tuple[L-1-i] = 'G';
            break;
         case T:
            tuple[L-1-i] = 'T';
            break;
      }
   }
   tuple[L] = '\0'; //Terminate string
   //printf("Decoded %s from %u\n", tuple, coded);
}

/************** BLOOM FILTER *********************/

uint64_t * bloom_filter (char ** spectrum, unsigned int n) {
   unsigned int i;
	uint64_t * bloom;				//pointer to filter

	//Allocate memory for Bloom Filter
	if(!(bloom = (uint64_t*)malloc(n*sizeof(uint64_t)))){		//assigns 64*n bits
		fprintf(stdout, "Error: Not enough memory\n");
		exit(1);
	}
	for(i=0; i<n; i++)
		bloom[i]=0;
	
	//Fill up the filter
	for(i=0; i<n; i++) {
		HashRead(bloom, spectrum[i], n);
	}
	
	return bloom;
}

//Sets bit into bloom filter
void SetBit(uint64_t* filter, int i, unsigned int n){
	int k = i % n;
	int pos = i % 64;
	uint64_t bit = 1;
	
	bit = bit << pos;
	//printf(", cell: %lu, bit: %lu\n", k, pos);
	
	filter[k] = filter[k] | bit;
}

//Takes the read and feeds it to hash functions
void HashRead(uint64_t* filter, char* read, unsigned int n){
	int i;
	
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

/****************** UTILITY **********************/

void list_to_array(list_t * head, char ** str_array) {
   list_t * tmp = head;
   int i=0;
   while(tmp) {
      strcpy(str_array[i], tmp->str);
      //printf("%s\n", str_array[i]);
      tmp = tmp->next;
      i++;
   }
}

list_t * add_elem(list_t * head, char * str) {
   if(!head) {
      if(!(head = (list_t *)malloc(sizeof(list_t)))) {
         fprintf(stdout, "Error: allocation\n");
         exit(1);
      }
      strcpy(head->str, str);
      head->next = NULL;
      return head;
   }
   
   list_t * new;
   if(!(new = (list_t *)malloc(sizeof(list_t)))) {
      fprintf(stdout, "Error: allocation\n");
      exit(1);
   }
   strcpy(new->str, str);
   new->next=head;
   return new;
}

void rm_newline(char * str) {
   int i;
   for(i=0; i<strlen(str) && str[i] != '\n'; i++)
      ;
   str[i] = '\0';
}

void free_list(list_t * head) {
   list_t * curr;
   while ( (curr = head) ) {
      head = head->next;
      free(curr);
   }
}
