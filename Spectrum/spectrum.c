#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define A 0 //00
#define C 1 //01
#define G 2 //10
#define T 3 //11

#define INTSIZE 16
#define INT64SIZE 64
#define PARAM 4 //l, m, input file, output file
#define BASES 4.0 //Need double

#ifdef SILENT
#define printf(...)
#endif

void rm_newline(char *);

/* Compute the spectrum T_m,l
 * 
 */
int main (int argc, char * argv[]) {
   int l; //lenght of l-tuple
   int m; //multiplicity
   int i, j;
   const double spectrum_size = 1048576; //pow(BASES, (double)l);
   
   //Read parameters
   if(argc < PARAM+1) {
      fprintf(stdout, "Error: parameters\n");
      exit(1);
   }
   
   l = atoi(argv[1]);
   m = atoi(argv[2]);
   char * inputFile;
   if(!(inputFile = (char *)malloc(sizeof(char) * strlen(argv[3])))) {
      printf("Error: allocation\n");
      exit(1);
   }
   strcpy(inputFile, argv[3]);
   
   char * outputFile;
   if(!(outputFile = (char *)malloc(sizeof(char) * strlen(argv[4])))) {
      printf("Error: allocation\n");
      exit(1);
   }
   strcpy(outputFile, argv[4]);
   
  /* Counter of occurences
   * The sequence is represented encoded in 2 bits for each character
   * So, from the address of the array cell we can find the original sequence.
   */
   uint8_t * counter;
   uint64_t arrayIndex;
   if((!(counter = (uint8_t *)malloc(sizeof(uint8_t) * spectrum_size)))) {
      printf("Error: allocation\n");
      exit(1);
   }
   for(arrayIndex=0; arrayIndex<spectrum_size; arrayIndex++) {
      counter[arrayIndex] = 0;
   }
   
   //Addressing the counter
   char * tuple;
   if(!(tuple = (char *)malloc(sizeof(char) * (l+1)))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   
   //Iterate over input file
   FILE * in;
   in = fopen(inputFile, "r");
   
   fgets(tuple, l+1, in);
   while (!feof(in)) {
      rm_newline(tuple);
      //Tokenize tuple
      arrayIndex = 0;
      for(i=0; i<l; i++) {
         arrayIndex = arrayIndex << 2; //Move bit to left
         switch(tuple[i]) {
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
      
      
      /*
      //DEBUG
      if(arrayIndex == 0)
         fprintf(stdout, "Tuple: %s, index: %x, count: %d\n", tuple, arrayIndex, counter[arrayIndex]);
      */
      
      //Counter update procedure
      //Check if necessary (< m)
      if(counter[arrayIndex] < m) {
         counter[arrayIndex]++;
      }
      printf("Tuple: %s, index: %x, count: %d\n", tuple, arrayIndex, counter[arrayIndex]);
      //tuple read, read next
      fgets(tuple, l+1, in);
      fgets(tuple, l+1, in);
   }
   
   fclose(in);
   free(inputFile);
   //End of counting section
   
   /*
   //DEBUG
   for(arrayIndex=0; arrayIndex < spectrum_size; arrayIndex++) {
      fprintf(stdout, "Index: %x, counter: %d\n", arrayIndex, counter[arrayIndex]);
   }
   //END DEBUG
   */
   
   FILE * out;
   out = fopen(outputFile, "w");
   
   /*
   //DEBUG
   FILE * debugFP;
   debugFP = fopen("debug.txt", "w");
   //END DEBUG
   */
   
   uint64_t tmp;
   uint64_t mask;
   //Start of writing output section
   for(arrayIndex=0; arrayIndex < spectrum_size; arrayIndex++) {
      if(counter[arrayIndex] >= m) {
         printf("%x, count: %d", arrayIndex, counter[arrayIndex]);
         for(i=0;i<l;i++) {
            //Create mask
            mask = 3; //11 binary
            mask = mask << (i*2);
            tmp = arrayIndex & mask;
            tmp = tmp >> (i*2);
            switch(tmp) {
               case A:
                  tuple[l-1-i] = 'A';
                  break;
               case C:
                  tuple[l-1-i] = 'C';
                  break;
               case G:
                  tuple[l-1-i] = 'G';
                  break;
               case T:
                  tuple[l-1-i] = 'T';
                  break;
            }
         }
         tuple[l] = '\0'; //Terminate string
         fprintf(out, "%s\n", tuple);
      }
      
      /*
      //DEBUG
      else {
         for(i=0; i<l; i++) {
            //Create mask
            mask = 3; //11 binary
            mask = mask << (i*2);
            tmp = arrayIndex & mask;
            tmp = tmp >> (i*2);
            switch(tmp) {
               case A:
                  tuple[l-1-i] = 'A';
                  break;
               case C:
                  tuple[l-1-i] = 'C';
                  break;
               case G:
                  tuple[l-1-i] = 'G';
                  break;
               case T:
                  tuple[l-1-i] = 'T';
                  break;
            }
         }
         tuple[l] = '\0'; //Terminate string
         fprintf(debugFP, "%s\n", tuple);
      }
      //END DEBUG
      */
      
   }
   
   /*
   //DEBUG
   fclose(debugFP);
   //END DEBUG
   */
   
   fclose(out);
   free(tuple);
   free(counter);
   free(outputFile);
   
   
   return 0;
}

void rm_newline(char * str) {
   int i;
   for(i=0; i<strlen(str) && str[i] != '\n'; i++)
      ;
   str[i] = '\0';
}
