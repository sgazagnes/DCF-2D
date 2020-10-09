#ifndef LAMBDAVEC_H
#define LAMBDAVEC_H

typedef struct LambdaVec LambdaVec;

struct LambdaVec
{
   int NumLambdas;
   float *Lambda;    
}; 
/******************************************************************************/
/*                             Lambda Vector handling                         */
/******************************************************************************/

LambdaVec *LambdaVectorCreate(int size);

void LambdaVectorDelete(LambdaVec *lvec);

LambdaVec *LambdaVectorRead(char *fname, double imScale);

#endif
