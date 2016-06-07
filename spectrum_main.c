#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "BloomFilter/hashes.h"
#include "BloomFilter/city.h"
#include "BloomFilter/spooky.h"

#define PARAM 2 //input file (*.seq), output file (bloom filter)
#define READS_LENGHT 35
#define L 10 //l-tuple lenght
#define M 6 //multiplicity
#define SPECTRUM_MAX_SIZE 1048576 //BASES^L
#define BUFFER 100

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

typedef struct read_list_s {
   char str[READS_LENGHT];
   struct read_list_s * next;
} read_list;

void preprocess_input(char *, char **, int *, char **, int *);
void get_tuple(char *, char *, int);
int check_read(char *);

int compute_spectrum (char **, int, char **);

uint64_t * bloom_filter (char **, int);
void SetBit(uint64_t *, int, int);
void HashRead(uint64_t *, char *, int);

void free_list(read_list *);
read_list * extract_list(read_list *, char *);
read_list * insert_list(read_list *, char *);
read_list * create_list();
void rm_newline(char *);

/********************************* MAIN **************************************************************************************************************************/
int main (int argc, char * argv[]) {
   int i;
   FILE * out;
   
   if(argc < PARAM + 1) {
      fprintf(stdout, "Error: parameters\n");
      exit(1);
   }
      
   char ** reads;
   char ** tuples;
   int dim_reads;
   int dim_tuples;
   
   printf("Starting\n");
   
   //TODO precompute dimensions, can't allocate in function
   // OR work with list and restore array after
   
   //Acquire reads and L-tuples
   preprocess_input(argv[1], reads, &dim_reads, tuples, &dim_tuples);
   printf("Preprocessed");
   
   char ** spectrum;
   int spectrum_size;
   
   //Compute the spectrum
   spectrum_size = compute_spectrum(tuples, dim_tuples, spectrum);
   printf("Spectrum computed\n");
   
   //Apply bloom filter
   uint64_t * hashed_spectrum = bloom_filter(spectrum, spectrum_size);
   printf("Hashed\n");
   
	out = fopen(argv[2], "w+");
	for(i=0; i < spectrum_size; i++){
		fprintf(out, "%lu\n", hashed_spectrum[i]);
	}	
	fclose(out);
	free(hashed_spectrum);
	free(reads);
   free(tuples);
   free(spectrum);
}

/************************************ SPECTRUM PREPROCESSING ******************************************************************************************************/
void preprocess_input (char * inputFile, char ** reads, int * dim_reads_ptr, char ** tuples, int * dim_tuples_ptr) {
   int i;
   
   char * tuple;
   if(!(tuple = (char *)malloc(sizeof(char) * (L+1)))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   char * read;
   if(!(read = (char *)malloc(sizeof(char) * BUFFER))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   
   read_list * list_reads = NULL, * list_tuples = NULL;
   read_list * last_read = NULL, * last_tuple = NULL;
   
   list_reads = create_list();
   last_read = list_reads;
   list_tuples = create_list();
   last_tuple = list_tuples;
   
   int dim_reads = 0;
   int dim_tuples = 0;
   
   FILE * in = fopen(inputFile, "r");
   
   printf("Reading input\n");
   fgets(read, BUFFER, in);
   while(!feof(in)) {
      rm_newline(read);
      if(read[0] != '>') {
         //printf("Read: %s\n", read);
         if(check_read(read)) {
            last_read = insert_list(last_read, read);
            dim_reads++;
            //printf("Read: %s\n", read);
            for(i=0; i< (READS_LENGHT - L + 1); i++) {
               get_tuple(read, tuple, i);
               last_tuple = insert_list(last_tuple, tuple);
               dim_tuples++;
               //printf("\tTuple %d: %s\n", i, tuple);
            }
         }
      }
      fgets(read, BUFFER, in);
   }
   
   fclose(in);
   printf("Reads and tuples imported and cleaned, store in array\n");
   
   read_list * list_r = list_reads->next;
   read_list * list_t = list_tuples->next;
   
   /* Store all the reads in the proper dimension array
    * 
    * 
    */
   if(!(reads = (char **)malloc(sizeof(char *) * dim_reads))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   for(i=0; i<dim_reads; i++) {
      if(!(reads[i] = (char *)malloc(sizeof(char) * (READS_LENGHT+1)))) {
         fprintf(stdout, "Error allocation\n");
         exit(1);
      }
      //Allocated, store value
      list_reads = extract_list(list_r, reads[i]);
      //printf("%s\n", reads[i]);
   }
   
   /* Store all the tuples in the proper dimension array
    * 
    * 
    */
   if(!(tuples = (char **)malloc(sizeof(char *) * dim_tuples))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   for(i=0; i<dim_tuples; i++) {
      if(!(tuples[i] = (char *)malloc(sizeof(char) * (L+1)))) {
         fprintf(stdout, "Error allocation\n");
         exit(1);
      }
      //Allocated, store value
      list_tuples = extract_list(list_t, tuples[i]);
      //printf("%s\n", tuples[i]);
   }
   
   printf("Reads and tuples stored in the array\n");
   
   *dim_reads_ptr = dim_reads;
   *dim_tuples_ptr = dim_tuples;
   printf("Reads: %d, Tuples: %d\n", *dim_reads_ptr, *dim_tuples_ptr);
   
   /* TODO Not working!!!!
   printf("Freeing %x\n", list_reads);
   free_list(list_reads);
   printf("Freed\nFreeing %x\n", list_tuples);
   free_list(list_tuples);
   printf("Freed\n");
   */
   free(read);
   free(tuple);
}

int check_read(char * read) {
   int i;
   for(i=0; i<READS_LENGHT; i++) {
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

/****************************** SPECTRUM COMPUTATON **************************************************************************************/
int compute_spectrum (char ** tuples, int dim_tuples, char ** spectrum) {
   int i, j;
   
  /* Counter of occurences
   * The sequence is represented encoded in 2 bits for each character
   * So, from the address of the array cell we can find the original sequence.
   */
   uint8_t * counter;
   uint64_t arrayIndex;
   if((!(counter = (uint8_t *)malloc(sizeof(uint8_t) * SPECTRUM_MAX_SIZE)))) {
      printf("Error: allocation\n");
      exit(1);
   }
   for(arrayIndex=0; arrayIndex<SPECTRUM_MAX_SIZE; arrayIndex++) {
      counter[arrayIndex] = 0;
   }
   printf("Spectrum computation: initializated\n");
   
   for(i=0; i<dim_tuples; i++) {
      //Tokenize tuple
      printf("Tokenizing %s\n", tuples[i]);
      arrayIndex = 0;
      for(j=0; j<L; j++) {
         arrayIndex = arrayIndex << 2; //Move bit to left
         switch(tuples[i][j]) {
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
      printf("\t%x\n", arrayIndex);
      
      //Counter update procedure
      //Check if necessary (< M)
      if(counter[arrayIndex] < M) {
         counter[arrayIndex]++;
      }
   }
   //End of counting section
   printf("Counted\n");
   
   uint64_t tmp;
   uint64_t mask;
   char tuple[L+1];
   int dim=0;
   //Start of writing output section
   //Use list as support
   read_list * list = NULL, * last = NULL;
   list = create_list();
   last = list;
   for(arrayIndex=0; arrayIndex < SPECTRUM_MAX_SIZE; arrayIndex++) {
      if(counter[arrayIndex] >= M) {
         for(i=0;i<L;i++) {
            //Create mask
            mask = 3; //11 binary
            mask = mask << (i*2);
            tmp = arrayIndex & mask;
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
         last = insert_list(last, tuple);
         dim++;
      }
   }
   
   list = list->next;
   
   /* Store all the tuple of the spectrum in the proper dimension array
    * 
    * 
    */
   if(!(spectrum = (char **)malloc(sizeof(char *) * dim))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   for(i=0; i<dim; i++) {
      if(!(spectrum[i] = (char *)malloc(sizeof(char) * (L+1)))) {
         fprintf(stdout, "Error allocation\n");
         exit(1);
      }
      //Allocated, store value
      list = extract_list(list, spectrum[i]);
   }
   
   return dim;
}

/****************************** BLOOM FILTER ************************************************************************************************/
uint64_t * bloom_filter (char ** spectrum, int n) {
   int i;
	uint64_t * bloom;				//pointer to filter

	//Allocate memory for Bloom Filter
	if(!(bloom = (uint64_t*)malloc(n*sizeof(uint64_t)))){		//assigns 64*n bits
		fprintf(stdout, "Error: Not enough memory\n");
		exit(1);
	}
	for(int i=0; i<n; i++)
		bloom[i]=0;
	
	//Fill up the filter
	for(i=0; i<n; i++) {
		HashRead(bloom, spectrum[i], n);
	}
	
	return bloom;
}

//Sets bit into bloom filter
void SetBit(uint64_t* filter, int i, int n){
	int k = i % n;
	int pos = i % 64;
	uint64_t bit = 1;
	
	bit = bit << pos;
	//printf(", cell: %lu, bit: %lu\n", k, pos);
	
	filter[k] = filter[k] | bit;
}

//Takes the read and feeds it to hash functions
void HashRead(uint64_t* filter, char* read, int n){
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

/****************************** UTILITY *****************************************************************************************************/

void rm_newline(char * str) {
   int i;
   for(i=0; i<strlen(str) && str[i] != '\n'; i++)
      ;
   str[i] = '\0';
}

read_list * extract_list(read_list * head, char * str) {
   read_list * tmp;
   if(head) {
      tmp = head;
      head = head->next;
      strcpy(str,tmp->str);
      return head;
   }
   strcpy(str,"ERROR");
   return NULL;
}

read_list * insert_list(read_list * last, char * read) {
   read_list * elem;
   
   while(last->next) {
      last = last->next;
   }
   
   if(!(elem = (read_list *)malloc(sizeof(read_list)))) {
      fprintf(stdout, "Error: allocation\n");
      exit(1);
   }
   elem->next = NULL;
   strcpy(elem->str, read);
   last->next = elem;
   
   return last;
}

read_list * create_list() {
   read_list * head;
   if(!(head = (read_list *)malloc(sizeof(read_list)))) {
      fprintf(stdout, "Error: allocation\n");
      exit(1);
   }
   head->next = NULL;
   strcpy(head->str,"");
}

void free_list(read_list * t) {
   if(t->next) {
      free_list(t->next);
   }
   if(t) {
      free(t);
   }
}
