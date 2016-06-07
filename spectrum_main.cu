#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#define PARAM 2 //input file (*.seq), output file
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

#ifdef SILENT
#define printf(...)
#endif

typedef struct read_list_s {
   char str[READS_LENGHT];
   struct read_list_s * next;
} read_list;

void free_list(read_list *);
read_list * extract_list(read_list *, char *);
void insert_list(read_list *, char *);
void get_tuple(char *, char *, int);
void rm_newline(char *);
int check_read(char *);
void preprocess_input(FILE *, char **, int *, char **, int *);
int compute_spectrum (char **, int, char **);

/********************************* MAIN **************************************************************************************************************************/
int main (int argc, char * argv[]) {
   FILE * in, * out;
   
   if(argc < PARAM + 1) {
      fprintf(stdout, "Error: parameters\n");
      exit(1);
   }
   
   in = fopen(argv[1], "r");
   
   char ** reads;
   char ** tuples;
   int * dim_reads;
   int * dim_tuples;
   
   preprocess_input(in, reads, dim_reads, tuples, dim_tuples);
   
   char ** spectrum;
   int spectrum_size;
   
   spectrum_size = compute_spectrum(tuples, *dim_tuples, spectrum);
   
   //Apply bloom filter TODO
   uint64_t * hashed_spectrum;
}

/************************************ SPECTRUM PREPROCESSING ******************************************************************************************************/
void preprocess_input (FILE * in, char ** reads, int * dim_reads, char ** tuples, int * dim_tuples) {
   int i;
   
   char * tuple;
   if(!(tuple = (char *)malloc(sizeof(char) * L))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   char * read;
   if(!(read = (char *)malloc(sizeof(char) * BUFFER))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   
   read_list * list_reads = 0, * list_tuples=0;
   
   *dim_reads=0;
   *dim_tuples=0;
   
   fgets(read, BUFFER, in);
   while(!feof(in)) {
      rm_newline(read);
      printf("%s\n", read);
      if(read[0] != '>') {
         if(check_read(read)) {
            insert_list(list_reads, read);
            (*dim_reads)++;
            for(i=0; i< (READS_LENGHT - L + 1); i++) {
               get_tuple(read, tuple, i);
               insert_list(list_tuples, tuple);
               (*dim_tuples)++;
            }
         }
      }
      fgets(read, BUFFER, in);
   }
   
   /* Store all the reads in the proper dimension array
    * 
    * 
    */
   if(!(reads = (char **)malloc(sizeof(char *) * *dim_reads))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   for(i=0; i<*dim_reads; i++) {
      if(!(reads[i] = (char *)malloc(sizeof(char) * (READS_LENGHT+1)))) {
         fprintf(stdout, "Error allocation\n");
         exit(1);
      }
      //Allocated, store value
      list_reads = extract_list(list_reads, reads[i]);
   }
   
   /* Store all the tuples in the proper dimension array
    * 
    * 
    */
   if(!(tuples = (char **)malloc(sizeof(char *) * *dim_tuples))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   for(i=0; i<*dim_tuples; i++) {
      if(!(reads[i] = (char *)malloc(sizeof(char) * (L+1)))) {
         fprintf(stdout, "Error allocation\n");
         exit(1);
      }
      //Allocated, store value
      list_tuples = extract_list(list_tuples, tuples[i]);
   }
   
   free_list(list_tuples);
   free_list(list_reads);
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
   
   for(i=0; i<dim_tuples; i++) {
      //Tokenize tuple
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
      
      //Counter update procedure
      //Check if necessary (< M)
      if(counter[arrayIndex] < M) {
         counter[arrayIndex]++;
      }
   }
   //End of counting section
   
   uint64_t tmp;
   uint64_t mask;
   char tuple[L+1];
   int dim=0;
   //Start of writing output section
   //Use list as support
   read_list * list;
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
         insert_list(list, tuple);
         dim++;
      }
   }
   
   
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

/****************************** UTILITY *****************************************************************************************************/

void rm_newline(char * str) {
   int i;
   for(i=0; i<strlen(str) && str[i] != '\n'; i++)
      ;
   str[i] = '\0';
}

read_list * extract_list(read_list * head, char * str) {
   read_list * tmp;
   if(head != 0) {
      tmp = head;
      head = head->next;
      strcpy(str,tmp->str);
   }
   return head;
}

void insert_list(read_list * list, char * read) {
   read_list * head = list, * elem;
   
   if(head == 0) {
      if(!(head = (read_list *)malloc(sizeof(read_list)))) {
         fprintf(stdout, "Error: allocation\n");
         exit(1);
      }
      head->next = 0;
      strcpy(head->str, read);
      return;
   }
   
   while(head->next != 0) {
      head = head->next;
   }
   
   if(!(elem = (read_list *)malloc(sizeof(read_list)))) {
      fprintf(stdout, "Error: allocation\n");
      exit(1);
   }
   elem->next = 0;
   strcpy(elem->str, read);
   head->next = elem;
}

void free_list(read_list * head) {
   if(head->next == 0) {
      free(head);
      return;
   }
   free_list(head->next);
   free(head);
}
