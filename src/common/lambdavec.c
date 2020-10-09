#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "lambdavec.h"


/******************************************************************************/
/*                             Lambda Vector handling                         */
/******************************************************************************/

LambdaVec *LambdaVectorCreate(int size)
{
   LambdaVec *lvec;

   lvec = malloc(sizeof(LambdaVec));
   if (lvec==NULL)  return(NULL);

   lvec->NumLambdas = size;
   lvec->Lambda = malloc(size*sizeof(float));
   if (lvec->Lambda==NULL)
   {
     free(lvec);
      return(NULL);
   }
   return(lvec);
}

void LambdaVectorDelete(LambdaVec *lvec)
{
  free(lvec->Lambda);
  free(lvec);
} 

LambdaVec *LambdaVectorRead(char *fname, double imScale)
{
   FILE *infile;
   LambdaVec *l;
   int size;
   int c;

   infile = fopen(fname, "r");
   if (infile==NULL)  return(NULL);
   fscanf(infile, "LambdaVector\n");

   fscanf(infile, "%d\n", &size);

   l = LambdaVectorCreate(size);

   if(l)
   {

      for(c=0;c<size;c++)
	{   
	  fscanf(infile,"%f ",&l->Lambda[c]); 
	  l->Lambda[c] /= imScale;
	}  
   } 
   fclose(infile);
   return(l);
} 
