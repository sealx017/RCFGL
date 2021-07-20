#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>


/* screening
Author: Seal Souvik
Mar 21, 2021
CSPH
*/

bool screening(double* v, const int K, const double lam, const double rho, const double EPS)
{
    double sumT = 0;
    double threshold = rho;
    for(int k = 0; k<K-1;k++)
    {
        sumT += v[k];
        threshold += lam;
        if (fabs(sumT) + EPS> threshold )
            return false;
    }
    sumT += v[K-1];
    threshold += lam - rho;
    if (fabs(sumT)  + EPS> threshold)
        return false;

    for(int k = 1; k < K-1; k++)
        for(int j = k; j<K-1;j++)
        {
            sumT = 0;
            threshold = (j-k+1)*lam + 2*rho;
            for(int i = k; i<j+1;i++)
                sumT += v[i];
            if (fabs(sumT) + EPS> threshold )
                return false;
        }

    sumT = 0;
    threshold = rho;
    for(int k = K-1; k>0;k--)
    {
        sumT += v[k];
        threshold += lam;
        if (fabs(sumT) + EPS> threshold )
            return false;
    }
    return true;

}


size_t globalScreening(double* G, double* adj, const int N, const int K, const double lam, const double rho)
{
    size_t indxNum = 0;
    size_t indx,N2 = N*N;
    double *v = new double[K];
    for(int i = 0; i<N; i++)
        for(int j = i+1; j<N;j++)
        {
            indx = i*K+j*N*K;
            for(int k=0; k<K;k++)
                v[k] = G[indx + k];

            if(screening(v,K,lam,rho, 0.0))
                continue;

            adj[i+j*N] = adj[j+i*N] = 1;
            indxNum++;
        }
    delete [] v;
    return indxNum;
}


