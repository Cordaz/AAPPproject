#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define PARAM 2 //input file, output file
#define READS_LENGHT 35
#define BUFFER 100

#ifdef SILENT
#define printf(...)
#endif

int check_read(char *);
void rm_newline(char *);

/* 
 * 
 */
int main (int argc, char * argv[]) {
   int i;
   
   if(argc < PARAM+1) {
      fprintf(stdout, "Error: paramters\n");
      exit(1);
   }
   
   char * inputFile;
   if(!(inputFile = (char *)malloc(sizeof(char) * strlen(argv[1])))) {
      printf("Error: allocation\n");
      exit(1);
   }
   strcpy(inputFile, argv[1]);
   
   char * outputFile;
   if(!(outputFile = (char *)malloc(sizeof(char) * strlen(argv[2])))) {
      printf("Error: allocation\n");
      exit(1);
   }
   strcpy(outputFile, argv[2]);
   
   FILE * in, * out;
   in = fopen(inputFile, "r");
   out = fopen(outputFile, "w");
   
   char * read;
   if(!(read = (char *)malloc(sizeof(char) * BUFFER))) {
      fprintf(stdout, "Error allocation\n");
      exit(1);
   }
   
   int flag;
   
   fgets(read, BUFFER, in);
   while(!feof(in)) {
      rm_newline(read);
      printf("%s\n", read);
      if(read[0] != '>') {
         flag = check_read(read);
         //Write out
         if(flag) //If is complete, otherwise discard
            fprintf(out, "%s\n", read);
      }
      fgets(read, BUFFER, in);
   }
   
   fclose(in);
   fclose(out);
   free(inputFile);
   free(outputFile);
}

int check_read(char * read) {
   int i;
   for(i=0; i<READS_LENGHT; i++) {
      if(read[i] == 'N')
         return 0;
   }
   return 1;
}

void rm_newline(char * str) {
   int i;
   for(i=0; i<strlen(str) && str[i] != '\n'; i++)
      ;
   str[i] = '\0';
}
