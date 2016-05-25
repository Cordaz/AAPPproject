#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#define A 0 //00
#define C 1 //01
#define G 2 //10
#define T 3 //11

#define INTSIZE 16
#define INT64SIZE 64
#define PARAM 4 //l, m, input file, output file

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
   unsigned int spectrum_size = 262144; // = (4^10) / 4 for now assuming l=10 fixed
   
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
   * So, from the address of the array cell and the frame (1,2,3 or 4)
   * we can find the original sequence.
   */
   uint16_t * counter;
   uint64_t arrayIndex;
   if((!(counter = (uint16_t *)malloc(sizeof(uint16_t) * spectrum_size)))) {
      printf("Error: allocation\n");
      exit(1);
   }
   for(arrayIndex=0; arrayIndex<spectrum_size; arrayIndex++) {
      counter[arrayIndex] = 0;
   }
   
   //Addressing the counter
   unsigned short int frame;
   uint64_t adder;
   uint16_t count;
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
      //Tokenize touple
      arrayIndex = 0;
      frame = 0;
      for(i=0; i<l-1; i++) {
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
         if(i != l-2)
            arrayIndex = arrayIndex << 2; //Move bit to left (if not last)
      }
      //The last one indexes the frame
      switch(tuple[l-1]) {
         case 'A':
            frame=A * 4;
            break;
         case 'C':
            frame=C * 4;
            break;
         case 'G':
            frame=G * 4;
            break;
         case 'T':
            frame=T * 4;
            break;
      }
      printf("Index: %lu, Frame: %u", arrayIndex, frame/4);
      //Counter update procedure
      //Check if necessary (< m)
      //Can't exceed 4 bits
      count = counter[arrayIndex];
      count = count << (INTSIZE - frame -4);
      count = count >> (INTSIZE - 4);
      //Now count contains the wanted count
      if(count < m) {
         //Preparing adder
         adder=1;
         adder = adder << frame; //Move the 1 of multiple of 4 position in the correct frame
         counter[arrayIndex] += adder; //Add one;
      }
      printf(", Count: %d\n", count);
      //tuple read, read next
      fgets(tuple, l+1, in);
   }
   
   fclose(in);
   free(inputFile);
   //End of counting section
   
   /*
   //DEBUG: print count
   for(arrayIndex=0; arrayIndex<spectrum_size; arrayIndex++) {
      for(j=0; j<4; j++) {
         count = counter[arrayIndex];
         count = count << (INTSIZE - (j*4) - 4);
         count = count >> (INTSIZE - 4);
         fprintf(stdout, "%3d", count);
      }
      fprintf(stdout, "\n");
   }
   */
   
   
   FILE * out;
   out = fopen(outputFile, "w");
   uint64_t tmp;
   //Start of writing output section
   for(arrayIndex=0; arrayIndex < spectrum_size; arrayIndex++) {
      //For each index 4 sequence must be extracted
      for(i=0;i<l-1;i++) {
         tmp = arrayIndex;
         //tmp = tmp << (INT64SIZE - i*2 - 2); //TODO - NOT WORKING...
         //tmp = tmp >> (INT64SIZE - 2);
         tmp = tmp >> (i*2);
         tmp = tmp % 4;
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
      for(j=0; j<4; j++) {
         count = counter[arrayIndex];
         count = count << (INTSIZE - (j*4) - 4);
         count = count >> (INTSIZE - 4);
         if(count >= m) {
            switch(j) {
               case A:
                  tuple[l-1] = 'A';
                  break;
               case C:
                  tuple[l-1] = 'C';
                  break;
               case G:
                  tuple[l-1] = 'G';
                  break;
               case T:
                  tuple[l-1] = 'T';
                  break;
            }
            //tuple now is the reconstructed sequence
            
            //Ready to write out
            fprintf(out, "%s\n", tuple);
         }
      }
      //All 4 sequence indexed by "arrayIndex" writed
   }
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
