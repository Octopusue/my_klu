#include "cs_op.h"
#include "cs.h"
#include <math.h>
cs *cs_transpose(const cs *A, int values)
{
    int p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w;
    double *Cx, *Ax;
    cs *C;

    if (!CS_CSC(A)) return NULL;


    m=A->m;n=A->n; 
    Ap=A->p; Ai=A->i;Ax=A->x;
    C = cs_spalloc(n, m, Ap[n], values&&Ax, 0);
    w = (int *)cs_calloc(m, sizeof(int));

    if (!C || !w) return (cs_done(C, w, NULL, 0));
    Cp=C->p; Ci=C->i; Cx=C->x;
    for (p=0;p<Ap[n];p++) w[Ai[p]]++;
    cs_cumsum(Cp, w, m);

    for (j=0;j<n;j++)
    {
        for (p=Ap[j]; p<Ap[j+1]; p++)
        {
            q = w[Ai[p]]++;
            Ci[q]=j;
            if (Cx) Cx[q]=Ax[p];
        }
    }
    return (cs_done(C, w, NULL, 1));

}


int cs_gaxpy(const cs *A, const double *x, double *y)
{
    int p, j, n, *Ap, *Ai;
    
    double *Ax;

    if (!CS_CSC(A) || !x || !y) return (0);

    n=A->n; Ap=A->p; Ai=A->i; Ax=A->x;
    /* p point first index of columns in i */
    for (j=0; j<n; j++)
    {
        for (p = Ap[j]; p<Ap[j+1]; p++)
        {
            y[Ai[p]] += Ax[p] * x[j];
        }
    }
    return 1;

}

int cs_scatter(const cs* A, int j, double beta, int *w, double *x, int mark, cs *C, int nz)
{
    int i, p, *Ap, *Ai, *Ci;
    double *Ax;

    if (!CS_CSC(A) || !w || !CS_CSC(C)) return -1;

    Ap = A->p; Ai = A ->i; Ax=A->x; Ci=C->i;
    for (p=Ap[j]; p<Ap[j+1];p++)
    {
        /* A row index is C row index */
        i = Ai[p];
        /*************************************** 
        * w[i] is zero initially, 
        * to make sure when we doing j columns of B,
        * w[i]<mark only work once
        * **************************************/
        if (w[i] < mark)
        {
            /* if i row of C is not calculated, calculate*/
            w[i] = mark;
            Ci[nz] = i;
            nz++;
            if (x) x[i] = beta * Ax[p];
        }
        else {
            /* i row of C is calculated, sum new solution*/
            if (x) x[i] += beta * Ax[p];
        }
    }
    return nz;

}
cs *cs_multiply(const cs* A, const cs* B)
{
    int p, j, nz=0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi;
    double *x, *Bx, *Cx;
    cs *C;

    if (!CS_CSC(A) || !CS_CSC(B)) return NULL;

    m=A->m; anz=A->p[A->n];
    n=B->n;Bp=B->p; Bi=B->i; Bx=B->x; bnz=Bp[n];
    w = (int *)cs_calloc(m, sizeof(int));
    values=(A->x !=NULL)&&(Bx!=NULL);
    x = values ? (double *)cs_malloc(m, sizeof(double)) : NULL;//tmp var for every columns after multiply
    C = cs_spalloc(m, n, anz+bnz, values, 0);
    
    if (!C || !w || (values && !x)) return (cs_done(C, w, x, 0));

    Cp = C->p;
    for (j = 0; j<n; j++)
    {
        if (nz+m>C->nzmax && !cs_sprealloc(C, 2*(C->nzmax) + m ))
        {
            return (cs_done(C, w, x, 0));
        }

        Ci=C->i; Cx=C->x;
        
        Cp[j] = nz; //start position of j columns

        for (p = Bp[j]; p<Bp[j+1];p++)
        {   
            /******************************************** 
            * Bi[p] is i row of B, and it need multiply i columns of A 
            * nz of C is accumulated
            * *******************************************/
            nz = cs_scatter(A, Bi[p], Bx? Bx[p] : 1, w, x, j+1, C, nz);
        }
        /*we have the j columns of C*/
        if (values) {
            for (p=Cp[j]; p<nz; p++)
            {
                Cx[p] = x[Ci[p]];
            }
        }


    }

    Cp[n] = nz;
    cs_sprealloc(C, 0); 
    return (cs_done(C, w, x, 1));


}


cs *cs_add(const cs* A, const cs* B, double alpha, double beta)
{
    int p, j, nz=0,anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values;
    double *x, *Bx, *Cx;
    cs *C;
    if (!CS_CSC(A) || !CS_CSC(B)) return NULL;

    m = A->m; anz = A->p[A->n];
    n = B->n; Bp = B->p; Bx = B->x; bnz = Bp[n];

    w = (int *)cs_calloc(m, sizeof(int));
    values = (A->x != NULL) && (Bx != NULL);
    x = values? (double *)cs_malloc(m, sizeof(double)) : NULL;

    C = cs_spalloc(m, n, anz + bnz, values, 0);

    if (!C || !w || (values && !x)) return (cs_done(C, w, x, 0));

    Cp = C->p; Ci = C->i; Cx = C->x;

    for (j=0;j<n;j++)
    {
        Cp[j] = nz;
        nz = cs_scatter(A, j, alpha, w, x, j+1, C, nz);
        nz = cs_scatter(B, j, beta, w, x, j+1, C, nz);
        if (values){
            for (p = Cp[j]; p < nz; p++)
            {
                Cx[p] = x[Ci[p]];
            }
        }
    }
    Cp[n]=nz;
    cs_sprealloc(C, 0);
    return (cs_done(C, w, x, 1));
}

int cs_pvec(const int *p, const double *b, double *x, int n)
{
    /* permutation vector : permutate position of k and p[k]*/
    
    if (!x || !b) return 0;
    for (int k = 0; k<n; k++)
    {
        x[k] = b[p ? p[k]: k];
    }

    return 1;
}

int cs_ipvec(const int *p, const double *b, double *x, int n)
{
    /* permutation vector : inverse permutate position of k and p[k]*/
    if (!x || !b) return 0;

    for (int k = 0; k < 0; k++)
    {
        x[p ? p[k] : k] = b[k];
    }
    return 1;
}

int *cs_pinv(int const *p, int n)
{
    int k, *pinv;
    if (!p) return NULL;
    pinv = (int *)cs_malloc(n, sizeof(int));
    if (!pinv) return NULL;
    for (k = 0; k < n; k++)
    {
        pinv[p[k]] = k;  /*inverse the permutation*/
    }
    return pinv;
}

cs *cs_permute(const cs* A, const int *pinv, const int *q, int values)
{
    /*C = PAQ, P permute rows, Q permute columns
    * here input pinv is inverse permutation vector of p */
    int t, j, k, nz=0, m, n, *Ap, *Ai, *Cp, *Ci;
    double *Cx, *Ax;
    cs *C;
    if (!CS_CSC(A)) return NULL;

    m = A->m; n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
    C = cs_spalloc(m, n, Ap[n], values && Ax != NULL, 0 );

    if (!C) return (cs_done(C, NULL, NULL, 0));
    /*generate new csc matrix: Cp->Cx->Ci*/
    for (k = 0; k < n; k++)
    {
        Cp[k] = nz;                             /* column k of C is column q[k] of A */
        j = q ? q[k] : k;                       /* column k permute to column j*/
        for (t = Ap[j]; t < Ap[j+1]; t++)
        {
            if (Cx) Cx[nz] = Ax[t];             /* row i of A if row pinv[i] of C */
            Ci[nz] = pinv ? pinv[Ai[t]] : Ai[t]; /* row Ai[t] permute to row pinv[Ai[t]]*/
            nz++;
        }

    }
    Cp[n] = nz;
    return (cs_done(C, NULL, NULL, 1));
}

cs *cs_symperm(const cs* A, const int* pinv, int values)
{
    //THIS IS NOT FULLY UNDERSTAND!!!!!!!!!!!!!!
    /****************************************************
    * symmetric matrix permutation:
    * this function computes C=PAP for symmetric matrix A 
    * whose uppper triangular part is stored
    * returning C in the same format
    * ***************************************************/
    int i, j,p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w;
    double *Cx, *Ax;
    cs *C;
    if (!CS_CSC(A)) return NULL;

    n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
    C = cs_spalloc(n, n, Ap[n], values && Ax != NULL, 0 );
    w = (int *)cs_calloc(n, sizeof(int)); /*record count of every column*/
    if (!C) return (cs_done(C, w, NULL, 0));
    Cp = C->p; Ci = C->i; Cx = C->x;
    /*generate new csc matrix: Cp->Cx->Ci*/
    for (j=0;j<n;j++)
    {
        j2 = pinv ? pinv[j] : j;
        for (p = Ap[j2]; p < Ap[j+1]; p++)
        {
            i = Ai[p];
            if (i>j) continue; /*only consider uppper triangular */
            i2 = pinv ? pinv[i] : i;
            w[CS_MAX(i2, j2)]++;  /* if permutation conducted to lower triangular place it to uppper part*/

        }
    }

    cs_cumsum(Cp, w, n); /*compute column pointers of C*/

    for (j = 0; j<n;j++)
    {
        j2 = pinv? pinv[j]: j;
        for (p=Ap[j]; p<Ap[j+1];p++)
        {
            i = Ai[p];
            if (i>j) continue;
            i2 = pinv? pinv[i] : i;

            q = w[CS_MAX(i2, j2)];

            w[CS_MAX(i2, j2)]++;

            Ci[q] = CS_MIN(i2, j2);
            if (Cx) Cx[q]=Ax[p];

        }
    }

    return (cs_done(C, w, NULL, 1));


}

double cs_norm(const cs *A)
{
    int p, j, n, *Ap;
    double *Ax, norm=0, s;
    if (!CS_CSC(A) || !A->x) return -1;

    n = A->n; Ap = A->p; Ax = A->x;

    for (j = 0; j < n; j++)
    {
        for (s = 0, p=Ap[j]; p<Ap[j+1]; p++)
        {
            s += fabs(Ax[p]);
            norm = CS_MAX(norm, s);
        }
    }
    return norm;

}
int cs_lsolve(const cs *L, double *x)
{
    int p, j, n, *Lp, *Li;
    double *Lx;
    if (!CS_CSC(L) || !x) return 0;
    n = L->n; Lp = L->p; Li = L->i; Lx = L->x;
    for (j =0; j<n; j++)
    {
        x[j]/=Lx[Lp[j]];
        for (p = Lp[j]+1; p < Lp[j+1]; p++)
        {
            x[Li[p]] -= Lx[p] * x[j];
        }
    } 
    return 1;
}
int cs_ltsolve(const cs *L, double *x)
{
    int p, j, n, *Lp, *Li;
    double *Lx;
    if (!CS_CSC(L) || !x) return 0;
    n = L->n; Lp = L->p; Li = L->i; Lx = L->x;
    for (j = n-1; j>=0; j--)
    {
        
        for (p = Lp[j]+1; p < Lp[j+1]; p++)
        {
            x[j] -= Lx[p] * x[Li[p]];
        }
        x[j] /= Lx[Lp[j]];
    } 
    return 1;
}

int cs_usolve(const cs *U, double *x)
{
    int p, j, n, *Up, *Ui;
    double *Ux;
    if (!CS_CSC(U) || !x) return 0;
    n = U->n; Up = U->p; Ui = U->i; Ux = U->x;
    for (j = n-1; j>=0; j--)
    {
        x[j]/=Ux[Up[j]];
        for (p = Up[j]; p < Up[j+1]-1; p++)
        {
            x[Ui[p]] -= Ux[p] * x[j];
        }
    } 
    return 1;
}

int cs_utsolve(const cs *U, double *x)
{
    int p, j, n, *Up, *Ui;
    double *Ux;
    if (!CS_CSC(U) || !x) return 0;
    n = U->n; Up = U->p; Ui = U->i; Ux = U->x;
    for (j = 0; j<n; j++)
    {
        
        for (p = Up[j]; p < Up[j+1]-1; p++)
        {
            x[j] -= Ux[p] * x[Ui[p]];
        }
        x[j] /= Ux[Up[j+1] -1];
    } 
    return 1;
}



int reachr(const cs *L, const cs *B, int *xi, int *w){
    int p, n = L->n;
    int top = n;
    for (p = B->p[0]; p<B->p[1]; p++) /*for each i in pattern of b*/
    {
        if (w[B->i[p]]!=1)
        {
            dfsr(B->i[p], L, &top, xi, w);
        }

    }
    return top;
}

void dfsr(int j, const cs *L, int *top, int *xi, int *w){
    int p;
    w[j] = 1;                               /*mark node j*/
    for (p = L->p[j]; p<L->p[j+1]; p++)     /*for each i in L (:, j)*/
    {
        if (w[L->i[p]]!=1)                  /*if i is unmarked*/
        {
            dfsr(L->i[p], L, top, xi, w);   /*start a dfs at i*/
        }

    }
    xi[--(*top)] = j;                       /*push j onto the stack*/
}

int cs_reach(cs *G, const cs *B, int k, int *xi, const int *pinv)
{
    int p, n, top, *Bp, *Bi, *Gp;
    if (!CS_CSC(G) || !CS_CSC(B) || !xi) return -1;

    n = G->n; Bp = B->p; Bi = B->i; Gp = G->p;
    top = n;

    for (p = Bp[k]; p<Bp[k+1]; p++)   /* k is 0 in this situation */
    {
        if (!CS_MARKED(Gp, Bi[p]))
        {
            top = cs_dfs(Bi[p], G, top, xi, xi+n, pinv);

        }
    }
    for (p = top; p<n; p++)
    {
        CS_MARK(Gp, xi[p]);/*restore G*/
    }
    return top;
}


int cs_dfs(int j, cs *G, int top, int *xi, int *pstack, const int *pinv)
{
    int i, p, p2, done, jnew, head =0, *Gp, *Gi;
    if (!CS_CSC(G) || !xi || !pstack) return (-1);

    Gp = G->p; Gi = G->i; 

    xi[0] = j;                                              /*initialize the recursion stack*/

    while (head >=0)
    {
        j = xi[head];                                   /*dfs push new node into stack*/
        jnew = pinv ? pinv[j] : j;                      /*permutation vector give new node*/

        if (!CS_MARKED(Gp, j))                          /*if the node is not marked*/
        {       
            CS_MARK(Gp, j);                             /*mark the node*/
            pstack[head] = (jnew < 0) ? 0 : CS_UNFLIP(Gp[jnew]);/*record recursion order*/
        }
        done = 1;                                       /*node j done if no unvisited neighbors*/

        p2 = (jnew < 0) ? 0 : CS_UNFLIP(Gp[jnew+1]);    

        for (p = pstack[head]; p < p2; p++)     
        {
            i = Gi[p];                                  /*examine all neighbors of j*/
            if (CS_MARKED(Gp, i)) continue;
            
            pstack[head] = p;                           /*push neighbors of j into recurison stack, record next column*/
            xi[++head] = i;                             /*push next rows*/
            done = 0;
        } 
        if (done)                                       /* if no more deep node exists, start to load nonzero node*/
        {
            head--;
            xi[--top] = j;                              /* head definitly will be less than top*/
        }

    }
    return top;
}


int cs_spsolve(cs *G, const cs *B, int k, int *xi, double *x, const int *pinv, int lo){
    /*lo: non-zero G=L, lo: zero G=U*/
    int j, J, p, q, px, top, n, *Gp, *Gi, *Bp, *Bi;
    double *Gx, *Bx;
    if (!CS_CSC(G) || !CS_CSC(B) || !xi || !x) return -1;
    Gp = G->p; Gi = G->i; Gx = G->x; n = G->n;
    Bp = B->p; Bi = B->i; Bx = B->x;
    top = cs_reach(G, B, k, xi, pinv);
    
    for (p = top; p < n; p++) x[xi[p]] = 0;             /*clear x*/
    for (p = Bp[k]; p< Bp[k+1]; p++) x[Bi[p]] = Bx[p];  /* make x=b */ 

    for (px = top; px < n; px++)
    {
        j = xi[px];
        J = pinv ? pinv[j] : j;
        if (J<0) continue;

        x[j] /= Gx[lo? Gp[J] : (Gp[J+1] - 1)]; /* x(j) /= G(i, j) */
        p = lo ? (Gp[J] +1) : (Gp[J]);
        q = lo ? (Gp[J+1]) : (Gp[J+1] -1);
        for (; p<q; p++)
        {
            x[Gi[p]] -= Gx[p] * x[j];
        }
    }
    return top;

}

int *cs_etree(const cs* A, int ata)
{
    int i, k, p, m, n, inext, *Ap, *Ai, *w, *parent, *ancestor, *prev;

    if (!CS_CSC(A)) return NULL; 

    m = A->m; n = A->n; Ap = A->p; Ai = A->i;
    parent = (int *) cs_malloc(n, sizeof(int));
    w  = (int *)cs_malloc(n +(ata ? m : 0), sizeof(int));
    if (!w || !parent) return (cs_idone(parent, NULL, w, 0));

    ancestor = w; prev = w + m;

    if (ata) for (i = 0; i < m; i++) prev[i] = -i;

    for (k = 0; k < n; k++)
    {
        parent[k] = -1;
        ancestor[k] = -1;
        for (p = Ap[k]; p < Ap[k+1]; p++) /*searching path on the uppper triangle matrix*/
        {
            i = ata ? prev[Ai[p]] : Ai[p];
            for (; i != -1 && i < k; i=inext) 
            {
                inext = ancestor[i];
                ancestor[i] = k;
                if (inext == -1) parent[i] = k; /*former root of i link to current node k*/
            }
            if (ata) prev[Ai[p]] = k;
        }

    }

    return (cs_idone(parent, NULL, w, i));

}

int cs_ereach(const cs *A, int k, const int *parent, int *s, int *w)
{
    int i, p, n, len, top, *Ap, *Ai;
    if (!CS_CSC(A) || ! parent || !s || !w) return -1;

    n = A->n; top = n; Ap = A->p; Ai = A->i;
    CS_MARK(w, k);

    for (p = Ap[k]; p < Ap[k+1]; p++)
    {
        i = Ai[p];
        if (i>k) continue;
        for (len=0; !CS_MARKED(w, i); i = parent[i]) /*traverse up etree*/
        {
            s[len] = i;
            len++;
            CS_MARK(w, i);
        }
        while(len>0)
        {
            s[--top] = s[--len];
        }
    }
    for (p = top; p < n; p++) CS_MARK(w, s[p]);
    CS_MARK(w, k);
    return top;
}

int *cs_post(const int * parent, int n){
    int j, k = 0, *post, *w, *head, *next, *stack;

    if (!parent) return NULL;

    post = (int *) cs_malloc(n, sizeof(int));
    w = (int *) cs_malloc(3*n, sizeof(int));

    if (!w || !post) return (cs_idone(post, NULL, w, 0));

    head = w; next = w+n; stack = w + 2 * n;

    for (j = 0; j < n; j++) head[j] = -1;

    for (j = n-1; j>=0; j--)        /*traverse node in reverse order*/
    {
        if (parent[j]==-1) continue;
        next[j] = head[parent[j]];  /**/
        head[parent[j]] = j;        /*head carry child of j*/
    }

    for (j = 0; j < n; j++)
    {
        if (parent[j] != -1) continue;
        k = cs_tdfs(j, k, head, next, post, stack); /*k is the top of post list */
    }
    return (cs_idone(post, NULL, w, 1));


}
int cs_tdfs(int j, int k, int *head, const int *next, int *post, int *stack){
    int i, p, top = 0;
    if (!head || !next || !post || !stack) return -1;

    while (top >= 0)
    {
        p = stack[top];
        i = head[p];            /*get child of p*/
        if (i == -1)
        {
            top--;
            post[k++] = p;
        }
        else
        {
            head[p] = next[i]; /*remove i form child of p, get another child of p*/
            stack[++top] = i; 
        }
    }
    return k;
}