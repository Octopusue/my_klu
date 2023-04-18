#ifndef CS_OP_H
#define CS_OP_H
#include "cs.h"
#include <math.h>
extern "C"{

/**************************
* Basic operation of Matrix
* *************************/
double cs_norm(const cs *A);

cs *cs_transpose(const cs *A, int values);

int cs_gaxpy(const cs *A, const double *x, double *y);

int cs_scatter(const cs* A, int j, double beta, int *w, double *x, int mark, cs *C, int nz);

cs *cs_multiply(const cs* A, const cs* B);

cs *cs_add(const cs* A, const cs* B, double alpha, double beta);

/**************************
* Calculating permutations of matrix
* *************************/
int cs_pvec(const int *p, const double *b, double *x, int n);

int cs_ipvec(const int *p, const double *b, double *x, int n);

int *cs_pinv(int const *p, int n);

cs *cs_permute(const cs* A, const int *pinv, const int *q, int values);

cs *cs_symperm(const cs* A, const int* pinv, int values);

int *cs_amd(int order, const cs *A);
/**************************
* Solving triangular systems
* Mx=b;
* first assume x=b;
* *************************/

/*dense right hand side*/

int cs_lsolve(const cs *L, double *x);

int cs_ltsolve(const cs *L, double *x);

int cs_usolve(const cs *U, double *x);

int cs_utsolve(const cs *U, double *x);

/*sparse right hand side*/
#define CS_FLIP(i) (-(i) - 2) /*convert index to negtive index*/
#define CS_UNFLIP(i) ((i < 0) ? CS_FLIP(i) : i) /*convert negtive index into original index*/
#define CS_MARKED(w, j) (w[j] < 0)
#define CS_MARK(w, j) {w[j] = CS_FLIP(w[j]);}
int reachr(const cs *L, const cs *B, int *xi, int *w);

void dfsr(int j, const cs *L, int *top, int *xi, int *w);

int cs_reach(cs *G, const cs *B, int k, int *xi, const int *pinv);
int cs_dfs(int j, cs *G, int top, int *xi, int *pstack, const int *pinv);
int cs_spsolve(cs *G, const cs *B, int k, int *xi, double *x, const int *pinv, int lo);

/*Cholesky factorization*/
int *cs_etree(const cs *A, int ata);
int cs_ereach(const cs *A, int k, const int *parent, int *s, int *w);

int *cs_post(const int * parent, int n);
int cs_tdfs(int j, int k, int *head, const int *next, int *post, int *stack);

void firstdesc(int n, int *parent, int *post, int *first, int *level);

int *rowcnt(cs *A, int *parent, int *post);
int cs_leaf(int i, int j, const int *first, int *maxfirst, int *prevleaf, int *ancestor, int *jleaf);


int *cs_counts(const cs *A, const int *parent, const int *post, int ata);

/*Up-Looking Cholesky*/
csn *cs_chol(const cs *A, const css *S);
/*Left-Looking Cholesky*/
}



#endif