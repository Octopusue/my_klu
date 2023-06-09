#include "cs_op.h"
#include "cs.h"
#include <math.h>

cs *cs_sort(const cs *A)
{
    cs *AT = cs_transpose(A, 1);
    return cs_transpose(AT, 1);
}
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
        
            Ci[q=w[Ai[p]]++]=j;
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
    C = cs_spalloc(m, n, anz+bnz, values, 0);// initiate enough mem;
    
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
            * B are iterating by rows and columns, downward for every column;
            * Bi[p] is k which need to contract with A[i][k];
            * once jump out of p iteration, C[][j] becomes settled
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
    Cp = C->p; Ci = C->i; Cx = C->x;
    for (k = 0; k < n; k++)
    {
        Cp[k] = nz;                             /* column k of C is column q[k] of A */
        j = q ? q[k] : k;                       /* column k permute to column j*/
        for (t = Ap[j]; t < Ap[j+1]; t++)
        {
            if (Cx) Cx[nz] = Ax[t];             /* row i of A if row pinv[i] of C */
            Ci[nz++] = pinv ? pinv[Ai[t]] : Ai[t]; /* row Ai[t] permute to row pinv[Ai[t]]*/
            
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
    w  = (int *)cs_malloc(n + (ata ? m : 0), sizeof(int));
    if (!w || !parent) 
        return (cs_idone(parent, NULL, w, 0));

    ancestor = w; prev = w + n;

    if (ata) for (i = 0; i < m; i++) prev[i] = -1;

    for (k = 0; k < n; k++)
    {
        parent[k] = -1;
        ancestor[k] = -1;/*ancestor keeps it's ancestor beyond it's father*/
        for (p = Ap[k]; p < Ap[k+1]; p++) /*searching path on the uppper triangle matrix*/
        {
            i = ata ? prev[Ai[p]] : Ai[p];
            for (; i != -1 && i < k; i=inext) /*traverse upper trangle matrix*/
            {
                inext = ancestor[i];
                ancestor[i] = k;                /*ancestor become large node, path compression*/
                if (inext == -1) parent[i] = k; /*former root of i link to current node k*/
            }
            /************ SEARCH FOR ELIMINATION TREE WHEN A IS NOT SYMMETIC *******************
            * if A is not sym pos denfinate matrix; we calculate tree of A'A
            * A(i, 0) may have values do not exist in A(0, i);
            * we have give up triangle search more node
            * It is column eliminate tree
            * for example ancestor of 6 maybe 1;
            * the next time node 6 need to append to tree
            * if 1 node does not attach to the tree, it will append 
            * 
            * what if node 6 does not show up?
            * even ata==0, node 6 become it own tree
            * ata==1, we always meet node 6 in iteration!!!!
            ************************************************************************************/
            if (ata) prev[Ai[p]] = k;           
        }

    }
    
   
    return (cs_idone(parent, NULL, w, i));

}

int cs_ereach(const cs *A, int k, const int *parent, int *s, int *w)
{
    /*traverse kth subtree from leaf to its root*/

    int i, p, n, len, top, *Ap, *Ai;
    if (!CS_CSC(A) || ! parent || !s || !w) return -1;

    n = A->n; top = n; Ap = A->p; Ai = A->i;
    CS_MARK(w, k);

    for (p = Ap[k]; p < Ap[k+1]; p++)
    {
        i = Ai[p];
        if (i>k) continue;                              /*still use upper triangular part of A*/
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
    /************************************************************
    * It is algrithm post order traverse the tree 
    * eliminate tree ordered by parent ptr，child node point parent node。
    * postordering tree ， it gives node orders; 
    * *********************************************************/
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
        next[j] = head[parent[j]];  /* find j‘s father‘s other child*/
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
    stack[0] = j;
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

void firstdesc(int n, int *parent, int *post, int *first, int *level)
{
    /*this function disjoint the elinimation tree into pieces
    * the smallest node of each path, must be leaves for some kth row subree
    * first descendant means the leaf node */
    int len, i, k, r, s;
    for (i = 0; i < n; i++) first[i] = -1;

    for (k = 0; k < n; k++)
    {
        i = post[k];
        len = 0;
        /******************************************** 
        * r == -1 represents root node, 
        * first[r] != -1 represents node r has visited
        * ********************************************/
        for (r = i; r!=-1 && first[r] == -1; r = parent[r], len++)
        {
            first[r] = k;
        }
        len += (r==-1) ? -1 : level[r];         /* root node : visited (end of path) */

        for (s = i; s != r; s = parent[s]) level[s] = len--;/*level represent length from root node*/
    }
}



int *rowcnt(cs *A, int *parent, int *post){
    /* row count is the first time traverse A matrix in lower triangular matrix
    * the kth row count traverse leaf to its k node or least common node of this subtree;
    * the number of node of the subtree is count for kth row
    *  
    */
    int i, j, k, len, s, p, jprev, q, n, sparent, jleaf;
    int  *Ap, *Ai, *maxfirst, *ancestor, *preleaf, *w, *first, *level, *rowcount;

    n = A->n; Ap = A->p; Ai = A->i;

    w = (int *) cs_malloc(5*n, sizeof(int));
    ancestor = w; maxfirst = w + n; preleaf = w + 2 * n; 
    first = w + 3 * n; level = w + 4 * n;

    rowcount = (int *) cs_malloc(n, sizeof(int));

    firstdesc(n, parent, post, first, level);

    for (i = 0; i < n; i++)
    {
        rowcount[i] = 1;
        preleaf[i] = -1;
        maxfirst[i] = -1;
        ancestor[i] = i;
    }

    for (k = 0; k < n; k++)
    {
        j = post[k];
        for (p = Ap[j]; p < Ap[j+1]; p++)
        {
            i = Ai[p];
            /*find leaf node add its level to kth subtree*/
            q = cs_leaf(i, j, first, maxfirst, preleaf, ancestor, &jleaf);/*low triangular tree*/
            if (jleaf){
                /* from small node, make sure leaf of i is start from j, node smaller than j is cut*/
                rowcount[i] += level[j] - level[q]; /*every j leaf node provides i row count*/
            }

        }
        if (parent[j] != -1){
            ancestor[j] = parent[j];
        }
    }
    cs_free(w);
    return rowcount;
}
int cs_leaf(int i, int j, const int *first, int *maxfirst, int *prevleaf, int *ancestor, int *jleaf){
    int q, s, sparent, jprev;

    if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return -1;
    *jleaf = 0;
    if (i <= j || first[j] <= maxfirst[i]) return -1;   /* j is not leaf of i */
    /*j is a leaf of i subtree*/
    printf("%d is leaf of %d \n", j, i);
    maxfirst[i] = j;
    jprev = prevleaf[i];                                /*previous leaf of i subtree*/
    prevleaf[i] = j;

    *jleaf = (jprev == -1) ? 1 : 2;                     /*j is first leaf of i or j is subsequent leaf of i*/

    if (*jleaf ==1) return i;
    /* **********************************************
    * find q root, it won't traverse to the last root, 
    * because it has not been defined outside this function
    * **************************************************/
    for (q = jprev; q != ancestor[q]; q = ancestor[q]); 

    for (s = jprev; s != q; s = sparent)
    {
        sparent = ancestor[s];
        ancestor[s] = q;                                /*path compression of path to its root*/
    }
    return q;                                           /* least common ancestor (jprev, j)*/
}

static void init_ata(cs *AT, const int *post, int *w, int **head, int **next){
    /*the elimination tree of A'A is constructed, it can be recover from post
    * if AT is low triangular matrix, post is equal to head */
    int i, k, p, m = AT->n, n = AT->m, *ATp=AT->p, *ATi = AT->i;
    *head = w + 4 * n; *next = w + 5 * n + 1;
    for (k =0; k<n; k++) w[post[k]] = k;
    for (i =0; i<m;i++)
    {
        for (k = n, p = ATp[i]; p < ATp[i+1]; p++)
            k = CS_MIN(k, w[ATi[p]]);
        (*next)[i] = (*head)[k];    //
        (*head)[k] = i;

        for (int it=0;it<n;it++)
            cout<<(*next)[it]<<'\t';
        cout<<endl;
        for (int it=0;it<n;it++)
            cout<<(*head)[it]<<'\t';
        cout<<endl<<endl;
    }
    /*head [i] == j : max column j having columns i*/
}
int *cs_counts(const cs *A, const int *parent, const int *post, int ata){
    /***************************************************************
     * 
     * col count L track of how many times j appears in any row subtree
     * ata is 1 : (A^T)A  is computed; ata is 0: A is computed
     * 
     * *************************************************************/
    int i, j, k, n, m, J, s, p, q, jleaf;
    int *ATp, *ATi, *maxfirst, *prevleaf, *ancestor, *head=NULL, *next=NULL, *colcount, *w, *first, *delta;
    cs *AT;
    if (!CS_CSC(A) || !parent || !post) return NULL;

    m = A->m; n = A->n;
    s = 4 * n + (ata ? (n+m+1) : 0);
    delta = colcount = (int *) cs_malloc(n, sizeof(int));
    w = (int *)cs_malloc(s, sizeof(int));
    /*when we transpose A matrix, we can use row count alogrithm to calculate column count*/
    AT = cs_transpose(A, 0);// IT SHOULD BE value = 1? no need
    // cs_print(AT, 1);
    if (!AT || !colcount || !w) return (cs_idone(colcount, AT, w, 0));

    ancestor = w; maxfirst = w + n; prevleaf = w + 2 * n; first = w + 3 * n;
    for (k =0; k < s; k++) w[k] = -1;

    for (k =0; k<n; k++)
    {
        j = post[k];
        delta[j] = (first[j] == -1) ? 1 : 0;        /* delta[j] == 1 if j is a leaf */
        for (; j != -1 && first[j] == -1; j = parent[j]) first[j] = k;
    }

    ATp = AT->p; ATi = AT ->i;
    cout<<endl;
    if (ata) init_ata(AT, post, w, &head, &next);
    for (int it=0;it<n;it++)
        cout<<post[it]<<'\t';
    cout<<endl;
    for (int it=0;it<n;it++)
        cout<<parent[it]<<'\t';
    cout<<endl;
    for (int it=0;it<n;it++)
        cout<<head[it]<<'\t';
    cout<<endl;
    for (int it=0;it<n;it++)
        cout<<next[it]<<'\t';
    cout<<endl;
    cs_print(AT, 0);
    for (i = 0; i < n; i++) ancestor[i] = i;
    for (k = 0; k < n; k++) 
    {
        j = post[k];
        if (parent[j] != -1) delta[parent[j]]--; /*if j is leaf delta[]*/
        for (J = HEAD(k, j); J!=-1; J=NEXT(J))/*Head jump from large node to small node !!!*/
        {   
            /*Head Next for A'A situation, */
            /**************************************************************
             * HEAD(k, j) has same function of post[k] -> tell us column order to tranverse; 
             * Head jump from large node to small node !!!
             * The outcome of Head + Next  == The outcome of post !!! 
             * when you print head and next you will find out !!!
             * At the beginning of Head point to large node which has nonzero in uppper triangular!!!
             * the Next point to the small node which has nonzero in low triangular !!!
             * only in this way we consider every possiable node in matrix!!!!!!
             * as done in etree where ata == 1, it traverse every nonzero node
             * **********************************************************/
            for (p = ATp[J]; p < ATp[J+1]; p++)
            {
                
                i = ATi[p];
                printf("k = %d, j = %d, J = %d, p = %d, i = %d, ATp[J]= %d\n", k, j, J, p, i, ATp[J+1]);
                /*******************************************************************
                 * although it calculate low triangular matrix of AT,
                 * it actually calculate upper triangular matrix of original matrix A;
                 *******************************************************************/
                q = cs_leaf(i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
                if (jleaf >= 1)
                    delta[j]++;         /*j become some leaf node*/
                if (jleaf == 2) 
                    delta[q]--;         /*q become a least common node*/
            }
        }
        if (parent[j] != -1) ancestor[j] = parent[j];
    }
    for (j = 0; j < n; j++){
        if (parent[j] != -1) colcount[parent[j]] += colcount[j];
    }
    for (int it=0;it<n;it++)
        cout<<colcount[it]<<'\t';
    cout<<endl;
    return (cs_idone(colcount, AT, w, 1));/*return L R^T*/

}

csn *cs_chol(const cs *A, const css *S)
{
    double d, lki, *Lx, *x, *Cx;
    int top, i, p, k, n, *Li, *Lp, *cp, *pinv, *s, *c, *parent, *Cp, *Ci;
    cs *L, *C, *E;
    csn *N;
    if (!CS_CSC(A) || !S || !S->cp || !S->parent) return (NULL);

    n = A->n;
    N = (csn *)cs_calloc(1, sizeof(csn));
    c = (int *)cs_malloc(2*n, sizeof(int));
    x = (double *)cs_malloc(n, sizeof(double));
    cp = S->cp; pinv = S->pinv; parent = S->parent;
    C = pinv? cs_symperm(A, pinv, 1) : ((cs *) A);
    E = pinv ? C : NULL; /*E is alias of A, or a copy E = A(p, p)*/

    if (!N || !c || !x || !C) return (cs_ndone(N, E, c, x, 0));

    s = c + n; 
    Cp = C->p; Ci = C->i; Cx = C->x;
    N->L = L = cs_spalloc(n, n, cp[n], 1, 0);
    if (!L) return (cs_ndone(N, E ,c, x, 0));
    Lp = L->p; Li = L->i; Lx = L->x;
    for (k = 0; k < n; k++) Lp[k] = c[k] = cp[k];
    for (k = 0; k < n; k++) /*compute L*L' = C   */
    {
        /*--------- Nonzero pattern of L(k,:) ---------------------------------*/
        top = cs_ereach(C, k, parent, s, c);        /*find nonzero pattern of L(K, :)*/
        x[k] = 0;
        for (p = Cp[k]; p < Cp[k+1]; p++)
        {
            if (Ci[p] <= k)
                x[Ci[p]]= Cx[p];                    /* as Lx=b, x=b, in cholesky L(:k, k) = A(:k, k)*/
        }
        d = x[k];                                   /* save L(k, k)*/
        x[k] = 0;
        /*--------- Triangular solve ---------------------------------*/
        for (; top < n; top++)
        {
            i = s[top];                             /*s[top .. n-1] is pattern of L(k, :)*/
            lki = x[i]/Lx[Lp[i]];                   /* as Lx = b , solve x[i] = x[i] / L[i][i] */
            x[i] = 0;
            for (p = Lp[i]+1; p < c[i]; p++){
                x[Li[p]] -= Lx[p] *lki;
            }
            d -= lki*lki;
            p = c[i]++;
            Li[p] = k;
            Lx[p] = lki;                            /*store L(k, i in set(0~k))*/
        }
        /*--------- Compute L(k, k) ---------------------------------*/
        if (d <= 0) return (cs_ndone(N, E, c, x ,0));
        p = c[k]++;
        Li[p] = k;
        Lx[p] = sqrt(d);
    } 
    Lp[n] = cp[n];
    return (cs_ndone(N, E, c, x, 1));

}


double cs_house(double *x, double *beta, int n){
    /*it ensures s = norm(x, 2). do not change norm(x)
    It computes beta and s and overwrites x with v*/
    double s, sigma = 0;
    int i;
    if (!x || !beta) return -1;

    for (i = 1; i < n; i++)/*only calculate entries except first one*/
        sigma += x[i]*x[i];

    if (sigma==0)
    {
        s = fabs(x[0]);
        (*beta) = (x[0] <= 0) ? 2 : 0;
        x[0] = 1;
    }
    else{
        s = sqrt(x[0] * x[0] + sigma);
        x[0] = (x[0] <= 0) ? (x[0] - s) : (-sigma/(x[0] + s));/*if x[0]<0 H should change it sign*/
        (*beta) = -1./(s*x[0]);
    }

    return s;
}
int cs_happly(const cs *V, int i, double beta, double *x){
    /******************************************************* 
     * happly means H apply
     * it apploes a HouseHolder reflection to a dense vector x, 
     * where v is sparse
     * it overwrites x with Hx = x - v*(beta*v'*x) 
     * 
     * HouseHolder Reflection convert column vector into new coorrdinate
     * 
     * ******************************************************/
    int p, *Vp, *Vi;
    double *Vx, tau = 0;

    if (!CS_CSC(V) || !x) return 0;
    Vp = V->p; Vi = V->i; Vx = V->x;
    for (p = Vp[i]; p < Vp[i+1]; p++)
    {
        tau += Vx[p] * x[Vi[p]];
    }
    tau *= beta;
    for (p = Vp[i]; p < Vp[i+1]; p++)
    {
        x[Vi[p]] -= Vx[p] * tau;
    }
    return 1;

}

static int cs_vcount(const cs *A, css *S)
{
    /*find the column counts of V matrix that holds the Householder vector*/
    int i, k, p, pa, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i, 
        *next, *head, *tail, *nque, *pinv, *leftmost, *w, *parent = S->parent;

    S->pinv = pinv = (int *)cs_malloc(n+m, sizeof(int));
    S->leftmost = leftmost = (int *)cs_malloc(m, sizeof(int));
    w = (int *)cs_malloc(m + 3 * n, sizeof(int));

    if (!pinv || !w || !leftmost) 
    {
        cs_free(w);     /* pinv and left most free later*/
        return (0);
    }

    next = w; head = w + m; tail = w + m + n; nque = w + m + 2 * n;

    for (k = 0; k < n; k++)
    {
        head[k] = -1;
        tail[k] = -1;
        nque[k] = 0;
    }
    for (i = 0; i < m; i++)
        leftmost[i] = -1;
    for (k = n-1; k >= 0; k--)
    {
        for (p = Ap[k]; p < Ap[k+1]; p++)
        {
            leftmost[Ai[p]] = k;            /*smallest kth non-zero node of each row*/
        }
    }
    /*split min col into clique*/
    for (i = m-1; i >= 0; i--)              /*scan rows in reverse order*/
    {
        pinv[i] = -1;
        k = leftmost[i];

        if (k == -1) continue;
        if (nque[k]++ == 0) 
            tail[k] = i;    /* first row in queue k, also clique of each leftmost[k]*/
        next[i] = head[k];                  /*just like col count link rows at certain order*/
        head[k] = i;                        /*head is row order to traverse*/
    }
    
    S->lnz = 0;
    S->m2 = m;
    cout<<"it"<<endl;
    for (int it=0;it<n;it++)
        cout<<it<<'\t';
    cout<<endl;
    cout<<"leftmost"<<endl;
    for (int it=0;it<n;it++)
        cout<<leftmost[it]<<'\t';
    cout<<endl;
    cout<<"head"<<endl;
    for (int it=0;it<n;it++)
        cout<<head[it]<<'\t';
    cout<<endl;
    cout<<"next"<<endl;
    for (int it=0;it<n;it++)
        cout<<next[it]<<'\t';
    cout<<endl;
    cout<<"nque"<<endl;
    for (int it=0;it<n;it++)
        cout<<nque[it]<<'\t';
    cout<<endl;
    cout<<"tail"<<endl;
    for (int it=0;it<n;it++)
        cout<<tail[it]<<'\t';
    cout<<endl;
    for (k = 0; k < n; k++)                 /*find row permutation and nnz(V)*/
    {   
        /* A～V 父结点是子节点的集合，如果把矩阵变成梯形矩阵，
        * 那么就不再需要把子结点的集合加到父节点上，
        * 也就是说，只要行交换变成梯形矩阵以后，non-zero A 就是non-zero V*/
        cout<<endl<<"iteration: "<<k<<endl;
        i = head[k];                        /*remove row i from queue k*/
        S->lnz++;                           /* count V(k, k) as non-zero*/
        if (i < 0) i = S->m2++;             /* add a fictitious row*/
        pinv[i] = k;                        /* associate row i with V(:, k) ???*/
        if (--nque[k] <= 0) continue;       /*skip if V(k+1:m, k)*/
        S->lnz += nque[k];                  /**/
        
        if ((pa = parent[k]) != -1)         /*eq~(5.3): k node contributes father node of k */
        {
            cout<<pa<<"=parent["<<k<<"]"<<endl; 
            if (nque[pa] == 0) tail[pa] = tail[k];    
            cout<<"tail["<<pa<<"] = tail["<<k<<"]"<<endl;         
            next[tail[k]] = head[pa];
            cout<<"next["<<tail[k]<<"] = head["<<pa<<"]"<<endl; 
            head[pa] = next[i];
            cout<<"head["<<pa<<"] = next["<<i<<"]"<<endl; 
            nque[pa] += nque[k];
            
        }
        cout<<"it"<<endl;
        for (int it=0;it<n;it++)
            cout<<it<<'\t';
        cout<<endl;
        cout<<"leftmost"<<endl;
        for (int it=0;it<n;it++)
            cout<<leftmost[it]<<'\t';
        cout<<endl;
        cout<<"pinv"<<endl;
        for (int it=0;it<n;it++)
            cout<<pinv[it]<<'\t';
        cout<<endl;
        cout<<"head"<<endl;
        for (int it=0;it<n;it++)
            cout<<head[it]<<'\t';
        cout<<endl;
        cout<<"next"<<endl;
        for (int it=0;it<n;it++)
            cout<<next[it]<<'\t';
        cout<<endl;
        cout<<"nque"<<endl;
        for (int it=0;it<n;it++)
            cout<<nque[it]<<'\t';
        cout<<endl;
        cout<<"tail"<<endl;
        for (int it=0;it<n;it++)
            cout<<tail[it]<<'\t';
        cout<<endl;
        
    }
    /*************************************************************
    * pinv[i] == k represents i row of A makes k row of V nonzero!!!!!!!!!!
    * pinv is core! nonzero node i of A after k iteration become nonzero node pinv[i], 
    * therefore V pattern is determined by pinv[i], where i is non-zero pattern of A;
    * the pinv makes PA matrix has non-zero at dignal entries
    * The row permutation pinv are ment to make matrix strong hall
    * 还是通过pinv 把矩阵拍成了阶梯矩阵
    * *************************************************************/
    for (i = 0; i < m; i++)
        if (pinv[i] < 0)
            pinv[i] = k++;
    
    // cs *permA = cs_permute(A, pinv, NULL, 1);
    // cs *tanspermA = cs_transpose(permA, 0);
    cs_free(w);
    return 1;
}
css *cs_sqr(int order, const cs* A, int qr){
    int n, k, ok = 1, *post;
    css *S;
    if (!CS_CSC(A)) return NULL;

    n = A->n;
    S = (css *) cs_calloc(1, sizeof(css));

    if (!S) return NULL;
    S->q = cs_amd(order, A);
    
    if (order && !S->q) return (cs_sfree(S));

    if (qr)
    {
        cs *C = order ? cs_permute(A, NULL, S->q, 0): (cs *)A;
        S->parent = cs_etree(C, 1);
        post = cs_post(S->parent, n);
        S->cp = cs_counts(C, S->parent, post, 1);
        cs_free(post);
        ok = C && S->parent && S->cp && cs_vcount(C, S);
        if (ok)
            for (S->unz = 0, k = 0; k < n; k++)
                S->unz += S->cp[k];
        ok = ok && (S->lnz >=0) && (S->unz >= 0);
        if (order)
            cs_spfree(C);
    }
    else{
        S->unz = 4 *(A->p[n]) + n; /*for LU factorization only*/
        S->lnz = S->unz;
    }
    return (ok? S:cs_sfree(S));

}
csn *cs_qr(const cs*A, const css *S){
    double *Rx, *Vx, *Ax, *Beta, *x;
    int i, k, p, m, n, vnz, p1, top, m2, len, col, rnz, 
        *s, *leftmost, *Ap, *Ai, *parent, *Rp, *Ri, *Vp, *Vi, *w, *pinv, *q;

    cs *R, *V;
    csn *N;

    if (!CS_CSC(A) || !S) return (NULL);

    m = A->m; n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
    q = S->q; parent = S->parent; pinv = S->pinv; m2 = S->m2;
    vnz = S->lnz; rnz = S->unz; leftmost = S->leftmost;
    
    w = (int *)cs_malloc(m2 + n, sizeof(int));
    x = (double *)cs_malloc(m2, sizeof(double));
    N = (csn *)cs_calloc(1, sizeof(csn));

    if (!w || !x || !N) return (cs_ndone(N, NULL, w, x, 0));

    s = w + m2;
    for (k = 0; k < m2; k++)
        x[k] = 0;
    N->L = V = cs_spalloc(m2, n, vnz, 1, 0);        /*V is triangular matrix, */
    N->U = R = cs_spalloc(m2, n, rnz, 1, 0);
    N->B = Beta = (double *)cs_malloc(n, sizeof(double));
    if (!R || !V || !Beta) return (cs_ndone(N, NULL, w, x, 0));

    Rp = R->p; Ri = R->i; Rx = R->x;
    Vp = V->p; Vi = V->i; Vx = V->x;

    for (i = 0; i <m2; i++) w[i] = -1;

    rnz = 0; vnz = 0;
    cout<<"pinv"<<endl;
    for (int it=0;it<n;it++)
        cout<<pinv[it]<<'\t';
    cout<<endl;
    for (k = 0; k < n; k++)                     /*tranverse column*/
    {
        /* calculate every column*/
        Rp[k] = rnz;
        Vp[k] = p1 = vnz;
        w[k] = k;
        Vi[vnz++] = k;                          /*digal entries*/
        top = n;
        col = q ? q[k] : k;                     /*column permutation*/
        for (p = Ap[col]; p < Ap[col+1]; p++)   /*traverse row*/
        {
            /******************************/
            i = leftmost[Ai[p]];                /*find min col*/
            //cout<<"Ai[p] = "<<Ai[p]<<"  i = "<<i<<endl;
            for (len = 0; w[i] != k; i = parent[i])/*traverse the whole clique for asymmetric matrix*/
            {
                s[len++] = i;
                w[i] = k;                       /*save max subtree, and mark i node*/
            }
            
            while (len > 0)
                s[--top] = s[--len];            /*push path into stack*/
            /**************
             * x[pinv[Ai[p]]] = Ax[p] 
             * in other word non-zero pattern of  A :
             * A_{ik}^{k-1 iteration} = pinv[A_{ik}]
             * ****************/
            i = pinv[Ai[p]];                    /*Is is how we determine V pattern : change row into trapezium*/
            x[i] = Ax[p];                       /*x[i] = A[i][col]*/
            cout<<"COL: "<<col<<"Ai[p] = "<<Ai[p]<<" i = "<<i<<"  Ax[p] = "<<Ax[p]<<endl;
            if (i > k && w[i] <k)               /*pattern of V(:, k) = x(k+1:m), upper triangular of A does not contribute to V*/
            {
                Vi[vnz++] = i;                  /* V pattern == leftmost A pattern */
                w[i] = k;
            }
        }
        printf("vnz = %d\n", vnz);
        for (p = top; p < n; p++)               /*k subtree traverse, 0th subtree skip*/
        {
            i = s[p];                           /*R(i, k) is nonzero*/
            cs_happly(V, i, Beta[i], x);        /*this make x multiply over and over eq~(5.1)!!!!!!*/
            Ri[rnz] = i;                        /*R nonzero pattern determined by Reach(min col) Tree */
            Rx[rnz++] = x[i];
            
            //printf("%d, %g", i, Rx[rnz-1]);

            x[i] = 0;
            if (parent[i] == k)
                vnz = cs_scatter(V, i, 0, w, NULL, k, V, vnz);//*k columns of V are set to zero 
        }
        //cs_print(R, 0);
        for (p = p1; p < vnz; p++)
        {
            Vx[p] = x[Vi[p]];                   /*V now represent A^{k-1}_{*,k}*/
            cout<<"Vi[p] = "<<Vi[p]<<"  Vx[p] "<<Vx[p]<<endl;
            x[Vi[p]] = 0;
        }
        Ri[rnz] = k;
        Rx[rnz++] = cs_house(Vx+p1, Beta+k, vnz-p1);/* V now become actual V */
        cout<<vnz-p1<<endl;
        //cs_print(V, 0);

    }

    Rp[n] = rnz;
    Vp[n] = vnz;
    return cs_ndone(N, NULL, w, x, 1);
}


/* clear w */
static csi cs_wclear (csi mark, csi lemax, csi *w, csi n)
{
    csi k ;
    if (mark < 2 || (mark + lemax < 0))
    {
        for (k = 0 ; k < n ; k++) if (w [k] != 0) w [k] = 1 ;
        mark = 2 ;
    }
    return (mark) ;     /* at this point, w [0..n-1] < mark holds */
}


/* keep off-diagonal entries; drop diagonal entries */
static csi cs_diag (csi i, csi j, double aij, void *other) { return (i != j) ; }

/* p = amd(A+A') if symmetric is true, or amd(A'A) otherwise */
csi *cs_amd (csi order, const cs *A)  /* order 0:natural, 1:Chol, 2:LU, 3:QR */
{
    cs *C, *A2, *AT ;
    csi *Cp, *Ci, *last, *W, *len, *nv, *next, *P, *head, *elen, *degree, *w,
        *hhead, *ATp, *ATi, d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1,
        k2, k3, jlast, ln, dense, nzmax, mindeg = 0, nvi, nvj, nvk, mark, wnvi,
        ok, cnz, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q, n, m, t ;
    csi h ;
    /* --- Construct matrix C ----------------------------------------------- */
    if (!CS_CSC (A) || order <= 0 || order > 3) return (NULL) ; /* check */
    AT = cs_transpose (A, 0) ;              /* compute A' */
    if (!AT) return (NULL) ;
    m = A->m ; n = A->n ;
    dense = CS_MAX (16, 10 * sqrt ((double) n)) ;   /* find dense threshold */
    dense = CS_MIN (n-2, dense) ;
    if (order == 1 && n == m)
    {
        C = cs_add (A, AT, 0, 0) ;          /* C = A+A' */
    }
    else if (order == 2)
    {
        ATp = AT->p ;                       /* drop dense columns from AT */
        ATi = AT->i ;
        for (p2 = 0, j = 0 ; j < m ; j++)
        {
            p = ATp [j] ;                   /* column j of AT starts here */
            ATp [j] = p2 ;                  /* new column j starts here */
            if (ATp [j+1] - p > dense) continue ;   /* skip dense col j */
            for ( ; p < ATp [j+1] ; p++) ATi [p2++] = ATi [p] ;
        }
        ATp [m] = p2 ;                      /* finalize AT */
        A2 = cs_transpose (AT, 0) ;         /* A2 = AT' */
        C = A2 ? cs_multiply (AT, A2) : NULL ;  /* C=A'*A with no dense rows */
        cs_spfree (A2) ;
    }
    else
    {
        C = cs_multiply (AT, A) ;           /* C=A'*A */
    }
    cs_spfree (AT) ;
    if (!C) return (NULL) ;
    cs_fkeep (C, &cs_diag, NULL) ;          /* drop diagonal entries */
    Cp = C->p ;
    cnz = Cp [n] ;
    P = (csi *)cs_malloc (n+1, sizeof (csi)) ;     /* allocate result */
    W = (csi *)cs_malloc (8*(n+1), sizeof (csi)) ; /* get workspace */
    t = cnz + cnz/5 + 2*n ;                 /* add elbow room to C */
    if (!P || !W || !cs_sprealloc (C, t)) return (cs_idone (P, C, W, 0)) ;
    len  = W           ; nv     = W +   (n+1) ; next   = W + 2*(n+1) ;
    head = W + 3*(n+1) ; elen   = W + 4*(n+1) ; degree = W + 5*(n+1) ;
    w    = W + 6*(n+1) ; hhead  = W + 7*(n+1) ;
    last = P ;                              /* use P as workspace for last */
    /* --- Initialize quotient graph ---------------------------------------- */
    for (k = 0 ; k < n ; k++) len [k] = Cp [k+1] - Cp [k] ;
    len [n] = 0 ;
    nzmax = C->nzmax ;
    Ci = C->i ;
    for (i = 0 ; i <= n ; i++)
    {
        head [i] = -1 ;                     /* degree list i is empty */
        last [i] = -1 ;
        next [i] = -1 ;
        hhead [i] = -1 ;                    /* hash list i is empty */
        nv [i] = 1 ;                        /* node i is just one node */
        w [i] = 1 ;                         /* node i is alive */
        elen [i] = 0 ;                      /* Ek of node i is empty */
        degree [i] = len [i] ;              /* degree of node i */
    }
    mark = cs_wclear (0, 0, w, n) ;         /* clear w */
    elen [n] = -2 ;                         /* n is a dead element */
    Cp [n] = -1 ;                           /* n is a root of assembly tree */
    w [n] = 0 ;                             /* n is a dead element */
    /* --- Initialize degree lists ------------------------------------------ */
    for (i = 0 ; i < n ; i++)
    {
        d = degree [i] ;
        if (d == 0)                         /* node i is empty */
        {
            elen [i] = -2 ;                 /* element i is dead */
            nel++ ;
            Cp [i] = -1 ;                   /* i is a root of assembly tree */
            w [i] = 0 ;
        }
        else if (d > dense)                 /* node i is dense */
        {
            nv [i] = 0 ;                    /* absorb i into element n */
            elen [i] = -1 ;                 /* node i is dead */
            nel++ ;
            Cp [i] = CS_FLIP (n) ;
            nv [n]++ ;
        }
        else
        {
            if (head [d] != -1) last [head [d]] = i ;
            next [i] = head [d] ;           /* put node i in degree list d */
            head [d] = i ;
        }
    }
    while (nel < n)                         /* while (selecting pivots) do */
    {
        /* --- Select node of minimum approximate degree -------------------- */
        for (k = -1 ; mindeg < n && (k = head [mindeg]) == -1 ; mindeg++) ;
        if (next [k] != -1) last [next [k]] = -1 ;
        head [mindeg] = next [k] ;          /* remove k from degree list */
        elenk = elen [k] ;                  /* elenk = |Ek| */
        nvk = nv [k] ;                      /* # of nodes k represents */
        nel += nvk ;                        /* nv[k] nodes of A eliminated */
        /* --- Garbage collection ------------------------------------------- */
        if (elenk > 0 && cnz + mindeg >= nzmax)
        {
            for (j = 0 ; j < n ; j++)
            {
                if ((p = Cp [j]) >= 0)      /* j is a live node or element */
                {
                    Cp [j] = Ci [p] ;       /* save first entry of object */
                    Ci [p] = CS_FLIP (j) ;  /* first entry is now CS_FLIP(j) */
                }
            }
            for (q = 0, p = 0 ; p < cnz ; ) /* scan all of memory */
            {
                if ((j = CS_FLIP (Ci [p++])) >= 0)  /* found object j */
                {
                    Ci [q] = Cp [j] ;       /* restore first entry of object */
                    Cp [j] = q++ ;          /* new pointer to object j */
                    for (k3 = 0 ; k3 < len [j]-1 ; k3++) Ci [q++] = Ci [p++] ;
                }
            }
            cnz = q ;                       /* Ci [cnz...nzmax-1] now free */
        }
        /* --- Construct new element ---------------------------------------- */
        dk = 0 ;
        nv [k] = -nvk ;                     /* flag k as in Lk */
        p = Cp [k] ;
        pk1 = (elenk == 0) ? p : cnz ;      /* do in place if elen[k] == 0 */
        pk2 = pk1 ;
        for (k1 = 1 ; k1 <= elenk + 1 ; k1++)
        {
            if (k1 > elenk)
            {
                e = k ;                     /* search the nodes in k */
                pj = p ;                    /* list of nodes starts at Ci[pj]*/
                ln = len [k] - elenk ;      /* length of list of nodes in k */
            }
            else
            {
                e = Ci [p++] ;              /* search the nodes in e */
                pj = Cp [e] ;
                ln = len [e] ;              /* length of list of nodes in e */
            }
            for (k2 = 1 ; k2 <= ln ; k2++)
            {
                i = Ci [pj++] ;
                if ((nvi = nv [i]) <= 0) continue ; /* node i dead, or seen */
                dk += nvi ;                 /* degree[Lk] += size of node i */
                nv [i] = -nvi ;             /* negate nv[i] to denote i in Lk*/
                Ci [pk2++] = i ;            /* place i in Lk */
                if (next [i] != -1) last [next [i]] = last [i] ;
                if (last [i] != -1)         /* remove i from degree list */
                {
                    next [last [i]] = next [i] ;
                }
                else
                {
                    head [degree [i]] = next [i] ;
                }
            }
            if (e != k)
            {
                Cp [e] = CS_FLIP (k) ;      /* absorb e into k */
                w [e] = 0 ;                 /* e is now a dead element */
            }
        }
        if (elenk != 0) cnz = pk2 ;         /* Ci [cnz...nzmax] is free */
        degree [k] = dk ;                   /* external degree of k - |Lk\i| */
        Cp [k] = pk1 ;                      /* element k is in Ci[pk1..pk2-1] */
        len [k] = pk2 - pk1 ;
        elen [k] = -2 ;                     /* k is now an element */
        /* --- Find set differences ----------------------------------------- */
        mark = cs_wclear (mark, lemax, w, n) ;  /* clear w if necessary */
        for (pk = pk1 ; pk < pk2 ; pk++)    /* scan 1: find |Le\Lk| */
        {
            i = Ci [pk] ;
            if ((eln = elen [i]) <= 0) continue ;/* skip if elen[i] empty */
            nvi = -nv [i] ;                      /* nv [i] was negated */
            wnvi = mark - nvi ;
            for (p = Cp [i] ; p <= Cp [i] + eln - 1 ; p++)  /* scan Ei */
            {
                e = Ci [p] ;
                if (w [e] >= mark)
                {
                    w [e] -= nvi ;          /* decrement |Le\Lk| */
                }
                else if (w [e] != 0)        /* ensure e is a live element */
                {
                    w [e] = degree [e] + wnvi ; /* 1st time e seen in scan 1 */
                }
            }
        }
        /* --- Degree update ------------------------------------------------ */
        for (pk = pk1 ; pk < pk2 ; pk++)    /* scan2: degree update */
        {
            i = Ci [pk] ;                   /* consider node i in Lk */
            p1 = Cp [i] ;
            p2 = p1 + elen [i] - 1 ;
            pn = p1 ;
            for (h = 0, d = 0, p = p1 ; p <= p2 ; p++)    /* scan Ei */
            {
                e = Ci [p] ;
                if (w [e] != 0)             /* e is an unabsorbed element */
                {
                    dext = w [e] - mark ;   /* dext = |Le\Lk| */
                    if (dext > 0)
                    {
                        d += dext ;         /* sum up the set differences */
                        Ci [pn++] = e ;     /* keep e in Ei */
                        h += e ;            /* compute the hash of node i */
                    }
                    else
                    {
                        Cp [e] = CS_FLIP (k) ;  /* aggressive absorb. e->k */
                        w [e] = 0 ;             /* e is a dead element */
                    }
                }
            }
            elen [i] = pn - p1 + 1 ;        /* elen[i] = |Ei| */
            p3 = pn ;
            p4 = p1 + len [i] ;
            for (p = p2 + 1 ; p < p4 ; p++) /* prune edges in Ai */
            {
                j = Ci [p] ;
                if ((nvj = nv [j]) <= 0) continue ; /* node j dead or in Lk */
                d += nvj ;                  /* degree(i) += |j| */
                Ci [pn++] = j ;             /* place j in node list of i */
                h += j ;                    /* compute hash for node i */
            }
            if (d == 0)                     /* check for mass elimination */
            {
                Cp [i] = CS_FLIP (k) ;      /* absorb i into k */
                nvi = -nv [i] ;
                dk -= nvi ;                 /* |Lk| -= |i| */
                nvk += nvi ;                /* |k| += nv[i] */
                nel += nvi ;
                nv [i] = 0 ;
                elen [i] = -1 ;             /* node i is dead */
            }
            else
            {
                degree [i] = CS_MIN (degree [i], d) ;   /* update degree(i) */
                Ci [pn] = Ci [p3] ;         /* move first node to end */
                Ci [p3] = Ci [p1] ;         /* move 1st el. to end of Ei */
                Ci [p1] = k ;               /* add k as 1st element in of Ei */
                len [i] = pn - p1 + 1 ;     /* new len of adj. list of node i */
                h = ((h<0) ? (-h):h) % n ;  /* finalize hash of i */
                next [i] = hhead [h] ;      /* place i in hash bucket */
                hhead [h] = i ;
                last [i] = h ;              /* save hash of i in last[i] */
            }
        }                                   /* scan2 is done */
        degree [k] = dk ;                   /* finalize |Lk| */
        lemax = CS_MAX (lemax, dk) ;
        mark = cs_wclear (mark+lemax, lemax, w, n) ;    /* clear w */
        /* --- Supernode detection ------------------------------------------ */
        for (pk = pk1 ; pk < pk2 ; pk++)
        {
            i = Ci [pk] ;
            if (nv [i] >= 0) continue ;         /* skip if i is dead */
            h = last [i] ;                      /* scan hash bucket of node i */
            i = hhead [h] ;
            hhead [h] = -1 ;                    /* hash bucket will be empty */
            for ( ; i != -1 && next [i] != -1 ; i = next [i], mark++)
            {
                ln = len [i] ;
                eln = elen [i] ;
                for (p = Cp [i]+1 ; p <= Cp [i] + ln-1 ; p++) w [Ci [p]] = mark;
                jlast = i ;
                for (j = next [i] ; j != -1 ; ) /* compare i with all j */
                {
                    ok = (len [j] == ln) && (elen [j] == eln) ;
                    for (p = Cp [j] + 1 ; ok && p <= Cp [j] + ln - 1 ; p++)
                    {
                        if (w [Ci [p]] != mark) ok = 0 ;    /* compare i and j*/
                    }
                    if (ok)                     /* i and j are identical */
                    {
                        Cp [j] = CS_FLIP (i) ;  /* absorb j into i */
                        nv [i] += nv [j] ;
                        nv [j] = 0 ;
                        elen [j] = -1 ;         /* node j is dead */
                        j = next [j] ;          /* delete j from hash bucket */
                        next [jlast] = j ;
                    }
                    else
                    {
                        jlast = j ;             /* j and i are different */
                        j = next [j] ;
                    }
                }
            }
        }
        /* --- Finalize new element------------------------------------------ */
        for (p = pk1, pk = pk1 ; pk < pk2 ; pk++)   /* finalize Lk */
        {
            i = Ci [pk] ;
            if ((nvi = -nv [i]) <= 0) continue ;/* skip if i is dead */
            nv [i] = nvi ;                      /* restore nv[i] */
            d = degree [i] + dk - nvi ;         /* compute external degree(i) */
            d = CS_MIN (d, n - nel - nvi) ;
            if (head [d] != -1) last [head [d]] = i ;
            next [i] = head [d] ;               /* put i back in degree list */
            last [i] = -1 ;
            head [d] = i ;
            mindeg = CS_MIN (mindeg, d) ;       /* find new minimum degree */
            degree [i] = d ;
            Ci [p++] = i ;                      /* place i in Lk */
        }
        nv [k] = nvk ;                      /* # nodes absorbed into k */
        if ((len [k] = p-pk1) == 0)         /* length of adj list of element k*/
        {
            Cp [k] = -1 ;                   /* k is a root of the tree */
            w [k] = 0 ;                     /* k is now a dead element */
        }
        if (elenk != 0) cnz = p ;           /* free unused space in Lk */
    }
    /* --- Postordering ----------------------------------------------------- */
    for (i = 0 ; i < n ; i++) Cp [i] = CS_FLIP (Cp [i]) ;/* fix assembly tree */
    for (j = 0 ; j <= n ; j++) head [j] = -1 ;
    for (j = n ; j >= 0 ; j--)              /* place unordered nodes in lists */
    {
        if (nv [j] > 0) continue ;          /* skip if j is an element */
        next [j] = head [Cp [j]] ;          /* place j in list of its parent */
        head [Cp [j]] = j ;
    }
    for (e = n ; e >= 0 ; e--)              /* place elements in lists */
    {
        if (nv [e] <= 0) continue ;         /* skip unless e is an element */
        if (Cp [e] != -1)
        {
            next [e] = head [Cp [e]] ;      /* place e in list of its parent */
            head [Cp [e]] = e ;
        }
    }
    for (k = 0, i = 0 ; i <= n ; i++)       /* postorder the assembly tree */
    {
        if (Cp [i] == -1) k = cs_tdfs (i, k, head, next, P, w) ;
    }
    return (cs_idone (P, C, W, 1)) ;
}


csi *cs_randperm (csi n, csi seed)
{
    csi *p, k, j, t ;
    if (seed == 0) return (NULL) ;      /* return p = NULL (identity) */
    p = (csi *)cs_malloc (n, sizeof (csi)) ;   /* allocate result */
    if (!p) return (NULL) ;             /* out of memory */
    for (k = 0 ; k < n ; k++) p [k] = n-k-1 ;
    if (seed == -1) return (p) ;        /* return reverse permutation */
    srand (seed) ;                      /* get new random number seed */
    for (k = 0 ; k < n ; k++)
    {
        j = k + (rand ( ) % (n-k)) ;    /* j = rand integer in range k to n-1 */
        t = p [j] ;                     /* swap p[k] and p[j] */
        p [j] = p [k] ;
        p [k] = t ;
    }
    return (p) ;
}

static csi cs_bfs (const cs *A, csi n, csi *wi, csi *wj, csi *queue,
    const csi *imatch, const csi *jmatch, csi mark)
{
    
    csi *Ap, *Ai, head = 0, tail = 0, j, i, p, j2 ;
    cs *C ;

    for (i =0; i < A->n; i++) 
        cout<<"("<<i<<","<<jmatch[i]<<")"<<endl;


    for (j = 0 ; j < n ; j++)           /* place all unmatched nodes in queue */
    {
        if (imatch [j] >= 0) continue ; /* skip j if matched */
        wj [j] = 0 ;                    /* j in set C0 (R0 if transpose) */
        queue [tail++] = j ;            /* place unmatched col j in queue */
    }
    if (tail == 0) return (1) ;         /* quick return if no unmatched nodes */
    C = (mark == 1) ? ((cs *) A) : cs_transpose (A, 0) ;
    if (!C) return (0) ;                /* bfs of C=A' to find R3,C3 from R0 */
    Ap = C->p ; Ai = C->i ;


    for (int it = head; it < tail; it++)
    {
        j = queue[it];
        cout<<"("<<imatch[j]<<","<<j<<")"<<endl;
    } 
        
    while (head < tail)                 /* while queue is not empty */
    {
        j = queue [head++] ;            /* get the head of the queue */
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;
            cout<<"("<<i<<","<<jmatch[i]<<")"<<endl;
            if (wi [i] >= 0) continue ; /* skip if i is marked */
            wi [i] = mark ;             /* i in set R1 (C3 if transpose) */
            j2 = jmatch [i] ;           /* traverse alternating path to j2 */
            if (wj [j2] >= 0) continue ;/* skip j2 if it is marked */
            wj [j2] = mark ;            /* j2 in set C1 (R3 if transpose) */
            queue [tail++] = j2 ;       /* add j2 to queue */
        }
    }
    if (mark != 1) cs_spfree (C) ;      /* free A' if it was created */
    return (1) ;
}

/* collect matched rows and columns into p and q */
static void cs_matched (csi n, const csi *wj, const csi *imatch, csi *p, csi *q,
    csi *cc, csi *rr, csi set, csi mark)
{
    csi kc = cc [set], j ;
    csi kr = rr [set-1] ;
    for (j = 0 ; j < n ; j++)
    {
        if (wj [j] != mark) continue ;      /* skip if j is not in C set */
        p [kr++] = imatch [j] ;
        q [kc++] = j ;
    }
    cc [set+1] = kc ;
    rr [set] = kr ;
}

/* collect unmatched rows into the permutation vector p */
static void cs_unmatched (csi m, const csi *wi, csi *p, csi *rr, csi set)
{
    csi i, kr = rr [set] ;
    for (i = 0 ; i < m ; i++) if (wi [i] == 0) p [kr++] = i ;
    rr [set+1] = kr ;
}

/* return 1 if row i is in R2 */
static csi cs_rprune (csi i, csi j, double aij, void *other)
{
    csi *rr = (csi *) other ;
    return (i >= rr [1] && i < rr [2]) ;
}
/* find the strongly connected components of a square matrix */
csd *cs_scc (cs *A)     /* matrix A temporarily modified, then restored */
{
    csi n, i, k, b, nb = 0, top, *xi, *pstack, *p, *r, *Ap, *ATp, *rcopy, *Blk ;
    cs *AT ;
    csd *D ;
    if (!CS_CSC (A)) return (NULL) ;                /* check inputs */
    n = A->n ; Ap = A->p ;
    D = cs_dalloc (n, 0) ;                          /* allocate result */
    AT = cs_transpose (A, 0) ;                      /* AT = A' */
    xi = (csi *)cs_malloc (2*n+1, sizeof (csi)) ;          /* get workspace */
    if (!D || !AT || !xi) return (cs_ddone (D, AT, xi, 0)) ;
    Blk = xi ; rcopy = pstack = xi + n ;
    p = D->p ; r = D->r ; ATp = AT->p ;
    top = n ;
    for (i = 0 ; i < n ; i++)   /* first dfs(A) to find finish times (xi) */
    {
        if (!CS_MARKED (Ap, i)) top = cs_dfs (i, A, top, xi, pstack, NULL) ;
    }
    for (i = 0 ; i < n ; i++) CS_MARK (Ap, i) ; /* restore A; unmark all nodes*/
    top = n ;
    nb = n ;
    for (k = 0 ; k < n ; k++)   /* dfs(A') to find strongly connnected comp */
    {
        i = xi [k] ;            /* get i in reverse order of finish times */
        if (CS_MARKED (ATp, i)) continue ;  /* skip node i if already ordered */
        r [nb--] = top ;        /* node i is the start of a component in p */
        top = cs_dfs (i, AT, top, p, pstack, NULL) ;
    }
    r [nb] = 0 ;                /* first block starts at zero; shift r up */
    for (k = nb ; k <= n ; k++) r [k-nb] = r [k] ;
    D->nb = nb = n-nb ;         /* nb = # of strongly connected components */
    for (b = 0 ; b < nb ; b++)  /* sort each block in natural order */
    {
        for (k = r [b] ; k < r [b+1] ; k++) Blk [p [k]] = b ;
    }
    for (b = 0 ; b <= nb ; b++) rcopy [b] = r [b] ;
    for (i = 0 ; i < n ; i++) p [rcopy [Blk [i]]++] = i ;
    return (cs_ddone (D, AT, xi, 1)) ;
}

/* Given A, compute coarse and then fine dmperm */
csd *cs_dmperm (const cs *A, csi seed)
{
    csi m, n, i, j, k, cnz, nc, *jmatch, *imatch, *wi, *wj, *pinv, *Cp, *Ci,
        *ps, *rs, nb1, nb2, *p, *q, *cc, *rr, *r, *s, ok ;
    cs *C ;
    csd *D, *scc ;
    /* --- Maximum matching ------------------------------------------------- */
    if (!CS_CSC (A)) return (NULL) ;            /* check inputs */
    m = A->m ; n = A->n ;
    D = cs_dalloc (m, n) ;                      /* allocate result */
    if (!D) return (NULL) ;
    p = D->p ; q = D->q ; r = D->r ; s = D->s ; cc = D->cc ; rr = D->rr ;
    jmatch = cs_maxtrans (A, seed) ;            /* max transversal */
    imatch = jmatch + m ;                       /* imatch = inverse of jmatch */
    if (!jmatch) return (cs_ddone (D, NULL, jmatch, 0)) ;
    /* --- Coarse decomposition --------------------------------------------- */
    wi = r ; wj = s ;                           /* use r and s as workspace */
    for (j = 0 ; j < n ; j++) wj [j] = -1 ;     /* unmark all cols for bfs */
    for (i = 0 ; i < m ; i++) wi [i] = -1 ;     /* unmark all rows for bfs */
    cs_bfs (A, n, wi, wj, q, imatch, jmatch, 1) ;       /* find C1, R1 from C0*/
    for (int it = 0; it < m; it++)
        cout<<wi[it]<<"\t";
    cout<<endl;
    for (int it = 0; it < n; it++)
        cout<<wj[it]<<"\t";
    cout<<endl;
    ok = cs_bfs (A, m, wj, wi, p, jmatch, imatch, 3) ;  /* find R3, C3 from R0*/
    for (int it = 0; it < m; it++)
        cout<<wi[it]<<"\t";
    cout<<endl;
    for (int it = 0; it < n; it++)
        cout<<wj[it]<<"\t";
    cout<<endl;
    if (!ok) return (cs_ddone (D, NULL, jmatch, 0)) ;
    cs_unmatched (n, wj, q, cc, 0) ;                    /* unmatched set C0 */
    cs_matched (n, wj, imatch, p, q, cc, rr, 1, 1) ;    /* set R1 and C1 */
    cs_matched (n, wj, imatch, p, q, cc, rr, 2, -1) ;   /* set R2 and C2 */
    cs_matched (n, wj, imatch, p, q, cc, rr, 3, 3) ;    /* set R3 and C3 */
    cs_unmatched (m, wi, p, rr, 3) ;                    /* unmatched set R0 */
    cs_free (jmatch) ;
    /* --- Fine decomposition ----------------------------------------------- */
    pinv = cs_pinv (p, m) ;         /* pinv=p' */
    if (!pinv) return (cs_ddone (D, NULL, NULL, 0)) ;
    C = cs_permute (A, pinv, q, 1) ;/* C=A(p,q) (it will hold A(R2,C2)) todo*/
    show_raw_matrix(C, "matchedM");
    cs_free (pinv) ;
    if (!C) return (cs_ddone (D, NULL, NULL, 0)) ;
    Cp = C->p ;
    nc = cc [3] - cc [2] ;          /* delete cols C0, C1, and C3 from C */
    if (cc [2] > 0) for (j = cc [2] ; j <= cc [3] ; j++) Cp [j-cc[2]] = Cp [j] ;
    C->n = nc ;
    if (rr [2] - rr [1] < m)        /* delete rows R0, R1, and R3 from C */
    {
        cs_fkeep (C, cs_rprune, rr) ;
        cnz = Cp [nc] ;
        Ci = C->i ;
        if (rr [1] > 0) for (k = 0 ; k < cnz ; k++) Ci [k] -= rr [1] ;
    }
    C->m = nc ;
    scc = cs_scc (C) ;              /* find strongly connected components of C*/
    if (!scc) return (cs_ddone (D, C, NULL, 0)) ;
    /* --- Combine coarse and fine decompositions --------------------------- */
    ps = scc->p ;                   /* C(ps,ps) is the permuted matrix */
    rs = scc->r ;                   /* kth block is rs[k]..rs[k+1]-1 */
    nb1 = scc->nb  ;                /* # of blocks of A(R2,C2) */
    for (k = 0 ; k < nc ; k++) wj [k] = q [ps [k] + cc [2]] ;
    for (k = 0 ; k < nc ; k++) q [k + cc [2]] = wj [k] ;
    for (k = 0 ; k < nc ; k++) wi [k] = p [ps [k] + rr [1]] ;
    for (k = 0 ; k < nc ; k++) p [k + rr [1]] = wi [k] ;
    nb2 = 0 ;                       /* create the fine block partitions */
    r [0] = s [0] = 0 ;
    if (cc [2] > 0) nb2++ ;         /* leading coarse block A (R1, [C0 C1]) */
    for (k = 0 ; k < nb1 ; k++)     /* coarse block A (R2,C2) */
    {
        r [nb2] = rs [k] + rr [1] ; /* A (R2,C2) splits into nb1 fine blocks */
        s [nb2] = rs [k] + cc [2] ;
        nb2++ ;
    }
    if (rr [2] < m)
    {
        r [nb2] = rr [2] ;          /* trailing coarse block A ([R3 R0], C3) */
        s [nb2] = cc [3] ;
        nb2++ ;
    }
    r [nb2] = m ;
    s [nb2] = n ;
    D->nb = nb2 ;
    cs_dfree (scc) ;
    return (cs_ddone (D, C, NULL, 1)) ;
}

// #include "cs.h"
/* find an augmenting path starting at column k and extend the match if found */
static void cs_augment (csi k, const cs *A, csi *jmatch, csi *cheap, csi *w,
        csi *js, csi *is, csi *ps)
{
    csi found = 0, p, i = -1, *Ap = A->p, *Ai = A->i, head = 0, j ;
    
    js [0] = k ;                        /* start with just node k in jstack */
    while (head >= 0)
    {
        /* --- Start (or continue) depth-first-search at node j ------------- */
        j = js [head] ;                 /* get j from top of jstack */
        if (w [j] != k)                 /* 1st time j visited for kth path */
        {
            w [j] = k ;                 /* mark j as visited for kth path */
            
            for (p = cheap [j] ; p < Ap [j+1] && !found ; p++)
            {
                i = Ai [p] ;            /* try a cheap assignment (i,j) */
                found = (jmatch [i] == -1) ;
            }
            cheap [j] = p ;             /* start here next time j is traversed*/
            if (found)
            {
                is [head] = i ;         /* column j matched with row i */
                break ;                 /* end of augmenting path */
            }
            ps [head] = Ap [j] ;        /* no cheap match: start dfs for j */
        }
        /* --- Depth-first-search of neighbors of j ------------------------- */
        for (p = ps [head] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;                /* consider row i */
            cout<<k<<"   "<<w [jmatch [i]]<<"   "<<jmatch [i]<<"   "<<i<<endl;
            if (w [jmatch [i]] == k) continue ; /* skip jmatch [i] if marked */
            ps [head] = p + 1 ;         /* pause dfs of node j */
            is [head] = i ;             /* i will be matched with j if found */
            js [++head] = jmatch [i] ;  /* start dfs at column jmatch [i] */
            break ;
        }
        if (p == Ap [j+1]) head-- ;     /* node j is done; pop from stack */
    }                                   /* augment the match if path found: */
    if (found) for (p = head ; p >= 0 ; p--) jmatch [is [p]] = js [p] ;
    cout<<k<<endl;
    for (i =0; i < A->n; i++) 
        cout<<"("<<i<<","<<jmatch[i]<<")"<<'\t';
    cout<<endl;
    
}

/* find a maximum transveral */
csi *cs_maxtrans (const cs *A, csi seed)  /*[jmatch [0..m-1]; imatch [0..n-1]]*/
{
    
    csi i, j, k, n, m, p, n2 = 0, m2 = 0, *Ap, *jimatch, *w, *cheap, *js, *is,
        *ps, *Ai, *Cp, *jmatch, *imatch, *q ;
    cs *C ;
    
    if (!CS_CSC (A)) return (NULL) ;                /* check inputs */
    n = A->n ; m = A->m ; Ap = A->p ; Ai = A->i ;
    w = jimatch = (csi *)cs_calloc (m+n, sizeof (csi)) ;   /* allocate result */
    if (!jimatch) return (NULL) ;
    for (k = 0, j = 0 ; j < n ; j++)    /* count nonempty rows and columns */
    {
        n2 += (Ap [j] < Ap [j+1]) ;
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            w [Ai [p]] = 1 ;
            k += (j == Ai [p]) ;        /* count entries already on diagonal */
        }
    }
    if (k == CS_MIN (m,n))              /* quick return if diagonal zero-free */
    {
        jmatch = jimatch ; imatch = jimatch + m ;
        for (i = 0 ; i < k ; i++) jmatch [i] = i ;
        for (      ; i < m ; i++) jmatch [i] = -1 ;
        for (j = 0 ; j < k ; j++) imatch [j] = j ;
        for (      ; j < n ; j++) imatch [j] = -1 ;
        return (cs_idone (jimatch, NULL, NULL, 1)) ;
    }
    for (i = 0 ; i < m ; i++) m2 += w [i] ;
    C = (m2 < n2) ? cs_transpose (A,0) : ((cs *) A) ; /* transpose if needed */
    if (!C) return (cs_idone (jimatch, (m2 < n2) ? C : NULL, NULL, 0)) ;
    n = C->n ; m = C->m ; Cp = C->p ;
    jmatch = (m2 < n2) ? jimatch + n : jimatch ;
    imatch = (m2 < n2) ? jimatch : jimatch + m ;
    w = (csi *)cs_malloc (5*n, sizeof (csi)) ;             /* get workspace */
    if (!w) return (cs_idone (jimatch, (m2 < n2) ? C : NULL, w, 0)) ;
    cheap = w + n ; js = w + 2*n ; is = w + 3*n ; ps = w + 4*n ;
    for (j = 0 ; j < n ; j++) cheap [j] = Cp [j] ;  /* for cheap assignment */
    for (j = 0 ; j < n ; j++) w [j] = -1 ;          /* all columns unflagged */
    for (i = 0 ; i < m ; i++) jmatch [i] = -1 ;     /* nothing matched yet */
    q = cs_randperm (n, seed) ;                     /* q = random permutation */
    for (k = 0 ; k < n ; k++)   /* augment, starting at column q[k] */
    {
        cs_augment (q ? q [k]: k, C, jmatch, cheap, w, js, is, ps) ;
    }
    cs_free (q) ;
    for (j = 0 ; j < n ; j++) imatch [j] = -1 ;     /* find row match */
    for (i = 0 ; i < m ; i++) if (jmatch [i] >= 0) imatch [jmatch [i]] = i ;


    for (i =0; i < n; i++) 
        cout<<"("<<i<<","<<jmatch[i]<<")"<<endl;

    return (cs_idone (jimatch, (m2 < n2) ? C : NULL, w, 1)) ;
}
