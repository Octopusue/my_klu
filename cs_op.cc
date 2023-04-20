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
    /*it ensures s = norm(x, 2). It computes beta and s and overwrites x with v*/
    double s, sigma = 0;
    int i;
    if (! x || !beta) return -1;

    for (i = 1; i < n; i++)
        sigma += x[i]*x[i];

    if (sigma==0)
    {
        s = fabs(x[0]);
        (*beta) = (x[0] <= 0) ? 2 : 0;
        x[0] = 1;
    }
    else{
        s = sqrt(x[0] * x[0] + sigma);
        x[0] = (x[0] <= 0) ? (x[0] - s) : (-sigma/(x[0] + s));
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

    next = w; head = w + n; tail = w + m + n; nque = w + +2 * n;

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
            leftmost[Ai[p]] = k;
        }
    }
    for (i = m-1; i >= 0; i--)
    {
        pinv[i] = -1;
        k = leftmost[i];
        if (k == -1) continue;
        if (nque[k]++ == 0) tail[k] = i;
        next[i] = head[k];
        head[k] = i;
    }
    S->lnz = 0;
    S->m2 = m;
    for (k = 0; k < n; k++)
    {
        i = head[k];
        S->lnz++;
        if (i < 0) i = S->m2++;
        pinv[i] = k;
        if (--nque[k] <= 0) tail[pa] = tail[k];
        next[tail[k]] = head[pa];
        head[pa] = next[i];
        nque[pa] += nque[k];
    }
    for (i = 0; i < m; i++)
        if (pinv[i] < 0)
            pinv[i] = k++;
    cs_free(w);
    return 1;


}
// css *cs_sqr(int order, const cs* A, int qr){
//     int n, k, ok = 1, *post;
//     css *S;
//     if (!CS_CSC(A)) return NULL;

//     n = A->n;
//     S = (css *) cs_calloc(1, sizeof(css));

//     if (!S) return NULL;
//     S->q = cs_amd(order, A);
//     if (order && !S->q) return (cs_sfree(S));

//     if (qr)
//     {
//         cs *C = order ? cs_permute(A, NULL, S->q, 0): (cs *)A;
//         S->parent = cs_etree(C, 1);
//         post = cs_post(S->parent, n);
//         S->cp = cs_counts(C, S->parent, post, 1);
//         cs_free(post);
//         ok = C && S->parent && S->cp && cs_vcount(C, S);
//         if (ok)
//             for (S->unz = 0, k = 0; k < n; k++)
//                 S->unz += S->cp[k];
//         ok = ok && (S->lnz >=0) && (S->unz >= 0);
//         if (order)
//             cs_spfree(C);
//     }
//     else{
//         S->unz = 4 *(A->p[n]) + n; /*for LU factorization only*/
//         S->lnz = S->unz;
//     }
//     return (ok? S:cs_sfree(S));

// }
// csn *cs_qr(const cs*A, const css *S){
//     double *Rx, *Vx, *Ax, *Beta, *x;
//     int i, k, p, m, n, vnz, p1, top, m2, len, col, rnz, 
//         *s, *leftmost, *Ap, *Ai, *parent, *Rp, *Ri, *Vp, *Vi, *w, *pinv, *q;

//     cs *R, *V;
//     csn *N;

//     if (!CS_CSC(A) || !S) return (NULL);

//     m = A->m; n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
//     q = S->q; parent = S->parent; pinv = S->pinv; m2 = S->m2;
//     vnz = S->lnz; rnz = S->unz; leftmost = S->leftmost;
    
//     w = (int *)cs_malloc(m2 + n, sizeof(int));
//     x = (double *)cs_malloc(m2, sizeof(double));
//     N = (csn *)cs_calloc(1, sizeof(csn));

//     if (!w || !x || !N) return (cs_ndone(N, NULL, w, x, 0));

//     s = w + m2;
//     for (k = 0; k < m2; k++)
//         x[k] = 0;
//     N->L = V = cs_spalloc(m2, n, vnz, 1, 0);
//     N->U = R = cs_spalloc(m2, n, rnz, 1, 0);
//     N->B = Beta = (double *)cs_malloc(n, sizeof(double));
//     if (!R || !V || !Beta) return (cs_ndone(N, NULL, w, x, 0));

//     Rp = R->p; Ri = R->i; Rx = R->x;
//     Vp = V->p; Vi = V->i; Vx = V->x;

//     for (i = 0; i <m2; i++) w[i] = -1;
//     rnz = 0; vnz = 0;
//     for (k = 0; k < n; k++)
//     {
//         Rp[k] = rnz;
//         Vp[k] = p1 = vnz;
//         w[k] = k;
//         Vi[vnz++] = k;
//         top = n;
//         col = q ? q[k] : k;
//         for (p = Ap[col]; p < Ap[col+1]; p++)
//         {
//             i = leftmost[Ai[p]];
//             for (len = 0; w[i] != k; i = parent[i])
//             {
//                 s[len++] = i;
//                 w[i] = k;
//             }
            
//             while (len > 0)
//                 s[--top] = s[--len];
//             i = pinv[Ai[p]];
//             x[i] = Ax[p];
//             if (i > k && w[i] <k)
//             {
//                 Vi[vnz++] = i;
//                 w[i] = k;
//             }
//         }
//         for (p = top; p < n; p++)
//         {
//             i = s[p];
//             cs_happly(V, i, Beta[i], x);
//             Ri[rnz] = i;
//             Rx[rnz++] = x[i];
//             x[i] = 0;
//             if (parent[i] == k)
//                 vnz = cs_scatter(V, i, 0, w, NULL, k, V, vnz);
//         }

//         for (p = p1; p < vnz; p++)
//         {
//             Vx[p] = x[Vi[p]];
//             x[Vi[p]] = 0;
//         }
//         Ri[rnz] = k;
//         Rx[rnz++] = cs_house(Vx+p1, Beta+k, vnz-p1);
    

//     }

//     Rp[n] = rnz;
//     Vp[n] = vnz;
//     return cs_ndone(N, NULL, w, x, 1);
// }
