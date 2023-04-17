#ifndef CS_H
#define CS_H
// #include<stdio.h>
// #include<stdlib.h>
#include<iostream>

using namespace std;
#define CS_CSC(A) (A&&(A->nz==-1))
#define CS_TRIPLET(A) (A&&(A->nz>=0))
#define CS_MAX(x, y) (x > y ? x : y)
#define CS_MIN(x, y) (x < y ? x : y)
#define CS_VER 0
#define CS_SUBVER 0
#define CS_SUBSUB 0
#define CS_DATE "20230414"
#define CS_COPYRIGHT "liqun.su"

#define HEAD(k, j) (ata ? head[k] : j)
#define NEXT(J) (ata ? next[J] : -1)

typedef struct cs_sparse
{
    int nzmax;   /*maximum number of entries*/
    int m;       /*number of rows*/
    int n;       /*number of columns*/
    int *p;      /*column ptr(size n+1) or col indices(size nzmax)*/
    int *i;      /*row indices, size nzmax*/
    double *x;   /*numerical values, size nzmax*/
    int nz;      /*#of entries in triplet maxtrix, -1 for compress-col*/
} cs;
extern "C" {
void show_cs_details(cs *A);


void *cs_malloc(int n, size_t size);


void *cs_calloc(int n, size_t size);


void *cs_free(void *p);


void *cs_realloc(void *p, int n, size_t size, int *ok);

cs *cs_spfree(cs *A);

cs *cs_spalloc(int m , int n, int nzmax, int values, int triplet);

int cs_sprealloc(cs *A, int nzmax);



int cs_entry(cs *T, int i, int j, double x);



cs *cs_done(cs *C, void *w, void *x, int ok);
int *cs_idone(int *p, cs *C, void *w, int ok);
double cs_cumsum(int *p, int *c, int n);

cs *cs_compress(const cs *T);



int cs_dupl(cs *A);


int cs_fkeep(cs *A, int (*fkeep)(int, int, double, void *), void *other);

static int cs_nonzero(int i, int j, double aij, void *other);

int cs_dropzeros(cs *A);


cs *cs_load(FILE *f);
int cs_print(const cs *A, int brief);
}
#endif