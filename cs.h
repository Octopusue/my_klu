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
#define CS_VER 2
#define CS_SUBVER 9
#define CS_SUBSUB 1
#define CS_DATE "20230414"
#define CS_COPYRIGHT "liqun.su"
#define csi int

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


typedef struct cs_symbolic
{
    int *pinv;
    int *q;
    int *parent;
    int *cp;
    int *leftmost;
    int m2;
    double lnz;
    double unz;
} css;

typedef struct cs_numeric
{
    cs *L;
    cs *U;
    int *pinv;
    double *B;

}csn;
typedef struct cs_dmperm_results    /* cs_dmperm or cs_scc output */
{
    csi *p ;        /* size m, row permutation */
    csi *q ;        /* size n, column permutation */
    csi *r ;        /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    csi *s ;        /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    csi nb ;        /* # of blocks in fine dmperm decomposition */
    csi rr [5] ;    /* coarse row decomposition */
    csi cc [5] ;    /* coarse column decomposition */
} csd ;

css *cs_schol(int order, const cs *A);

extern "C" {
void show_cs_details(cs *A);


void *cs_malloc(int n, size_t size);


void *cs_calloc(int n, size_t size);


void *cs_free(void *p);
csd *cs_dfree (csd *D);
csd *cs_dalloc (csi m, csi n);
csd *cs_ddone (csd *D, cs *C, void *w, csi ok);
void *cs_realloc(void *p, int n, size_t size, int *ok);

cs *cs_spfree(cs *A);
css *cs_sfree(css *S);
csn *cs_nfree(csn *N);
csn *cs_ndone(csn *N, cs *C, void *w, void *x, int ok);


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