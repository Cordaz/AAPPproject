#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define SPECTRUM_SIZE 262144 // = (4^10) / 4 for now assuming l=10 fixed

#define A 0 //00
#define C 1 //01
#define G 2 //10
#define T 3 //11

#define INTSIZE 16 //Compiler dependet... Fixed size? For now assume 16 bits

#ifdef SILENT
#define printf(...)
#endif

int main (int argc, char * argv[]) {
   int l=10; //lenght of l-tuple
   int m=6; //multiplicity
   int i;
   
  /* Counter of occurences
   * The sequence is represented encoded in 2 bits for each character
   * So, from the address of the array cell and the frame (1,2,3 or 4)
   * we can find the original sequence.
   */
   unsigned int counter[SPECTRUM_SIZE];
   
   //Addressing the counter
   unsigned long int arrayIndex;
   unsigned short int frame;
   unsigned long int adder;
   int count;
   char * tuple;
   if(!(tuple = (char *)malloc(sizeof(char) * l))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   
   //Should read tuple from the file and iterate
   
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
      arrayIndex = arrayIndex << 2; //Move bit to left
   }
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
   
   //Counter update procedure
   //Check if necessary (< m)
   count = counter[arrayIndex];
   count = count << INTSIZE - frame;
   count = count >> INTSIZE;
   //Now count contains the wanted count
   if(count < m) {
      //Preparing adder
      adder=1;
      adder = adder << frame; //Move the 1 of multiple of 4 position in the correct frame
      counter[arrayIndex] += adder; //Add one;
   }
   
   //End of counting section
   
   //Start of writing output section
   //Need to reconstruct the counted sequence from the encoding
   
   return 0;
}
