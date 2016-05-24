#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define PARAM 3 //l, input file, output file
#define READS_LENGHT 35
#define BUFFER 100

#ifdef SILENT
#define printf(...)
#endif

void get_tuple(char *, char *, int, int);
void rm_newline(char *);

/* Spectrum preprocessing, remove useless line and compute the sequences of lenght l
 * 
 */
int main (int argc, char * argv[]) {
   int i;
   
   if(argc < PARAM+1) {
      fprintf(stdout, "Error: paramters\n");
      exit(1);
   }
   
   int l;
   l = atoi(argv[1]);
   
   char * inputFile;
   if(!(inputFile = (char *)malloc(sizeof(char) * strlen(argv[2])))) {
      printf("Error: allocation\n");
      exit(1);
   }
   strcpy(inputFile, argv[2]);
   
   char * outputFile;
   if(!(outputFile = (char *)malloc(sizeof(char) * strlen(argv[3])))) {
      printf("Error: allocation\n");
      exit(1);
   }
   strcpy(outputFile, argv[3]);
   
   FILE * in, * out;
   in = fopen(inputFile, "r");
   out = fopen(outputFile, "w");
   
   char * tuple;
   if(!(tuple = (char *)malloc(sizeof(char) * l))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   char * read;
   if(!(read = (char *)malloc(sizeof(char) * BUFFER))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   
   fgets(read, BUFFER, in);
   while(!feof(in)) {
      rm_newline(read);
      printf("%s\n", read);
      if(read[0] != '>') {
         for(i=0; i< (READS_LENGHT - l + 1); i++) {
            get_tuple(read, tuple, i, l);
            //Write out
            fprintf(out, "%s\n", tuple);
         }
      }
      fgets(read, BUFFER, in);
   }
   
   fclose(in);
   fclose(out);
   free(inputFile);
   free(outputFile);
   free(tuple);
}

void get_tuple(char * read, char * tuple, int start, int l) {
   int i;
   for(i=0; i<l; i++) {
      tuple[i] = read[start+i];
   }
   tuple[l] = '\0';
}

void rm_newline(char * str) {
   int i;
   for(i=0; i<strlen(str) && str[i] != '\n'; i++)
      ;
   str[i] = '\0';
}
