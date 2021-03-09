#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>


/* screening
Author: Seal Souvik
Mar 8, 2021
CSPH
*/

/* const int N, const int K, const size_t Num */

void screening_loop(double* output, double* r1, double* s1, double* b1, double* r2, double* s2, double* b2, const int p, const int n1, const int n2)
{
    double T1ij, T2ij;
    size_t idij, idii, idjj, idjmini;
    for(size_t i = 0; i<(p-1); i++)
     {
      for(size_t j = (i+1); j<p; j++)
       {
        double temp = i*p+j;
        idij = (size_t) temp;
        temp = i*(p+1);
        idii = (size_t) temp;
        temp = j*(p+1);
        idjj = (size_t) temp;
        temp = (j-1)*p+i;
        idjmini = (size_t) temp;
        T1ij = (r1[idij]+s1[i]*b1[idij]+s1[j]*b1[idjmini])/(r1[idii]*r1[idjj]);
        T2ij = (r2[idij]+s2[i]*b2[idij]+s2[j]*b2[idjmini])/(r2[idii]*r2[idjj]);
        output[idij] = (T1ij-T2ij)/sqrt((1+(b1[idij]*b1[idij])*r1[idii]/r1[idjj])/(r1[idii]*r1[idjj]*n1)+(1+(b2[idij]*b2[idij])*r2[idii]/r2[idjj])/(r2[idii]*r2[idjj]*n2));
       }
     }
}

