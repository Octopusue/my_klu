
#include<iostream>
#include"cs.h"
#include"cs_op.h"
using namespace std;

void show_cs_details(cs *A)
{
    cout<<"nzmax: "<<A->nzmax<<endl;
    cout<<"m: "<<A->m<<" n: "<<A->n<<endl;
    int pSize=A->nz>=0? A->nzmax:A->n+1;
    for (int i=0;i<pSize;i++)
    {
        cout<<"p: "<<A->p[i]<<'\t';
    }
    cout<<endl;
    
    for (int i=0;i<A->nzmax;i++)
    {
        cout<<"i: "<<A->i[i]<<'\t';
    }
    cout<<endl;
    if (A->x)
    {
        for (int i=0;i<A->nzmax;i++)
        {
            cout<<"x: "<<A->x[i]<<'\t';
        }
        cout<<endl;
    }
    
    

    cout<<"nz: "<<A->nz<<endl;
    cout<<"*****************"<<endl;
}

cs *ccs2tri(const cs *A){
    int n = A->n, m =A->m, *Ap = A->p, *Ai = A->i, p;
    double *Ax = A->x;

    cs *triM = cs_spalloc(m, n, A->nzmax, 0, 1);
    for (int j = 0; j < n; j++)
    {
        p = Ap[j];
        for (int i = Ai[p]; i < Ai[p+1]; i++)
        {
            cs_entry(triM, i, j, Ax[p]);
        }
    }
    return triM;
}
int show_raw_matrix(const cs *A, const char *filename){
    
    FILE *fout = fopen(filename, "w+"); 
    int n = A->n, m =A->m, *Ap = A->p, *Ai = A->i, p;
    double *Ax = A->x;

    double **rawM, *colM;
    rawM = (double **)malloc(m*sizeof(double *));
    if (!rawM) return 0;

    for (int i = 0; i < m;i++)
    {
        rawM[i] = (double *)calloc(n, sizeof(double));
        if (!rawM[i]) 
        {
            for (int j = 0; j < i; j++)
            {
                free(rawM[j]);
            }
            free(rawM);
            return 0;
        }
    }

    for (int j = 0; j < n; j++)
    {
        p = Ap[j];
        for (int i = Ai[p]; i < Ai[p+1]; i++)
        {
            rawM[i][j] = Ax[p];
        }
    }
    for (int i = 0; i < m ; i++)
    {
        cout<<"|\t";
        fprintf(fout, "|\t");
        for (int j = 0; j < n ; j++)
        {
            cout<<rawM[i][j]<<'\t';
            fprintf(fout, "%g\t", rawM[i][j]);
        }
        cout<<"|\t"<<endl;
        fprintf(fout,"| \t\n");
    }
    fclose(fout);
    return 1;
}


void *cs_malloc(int n, size_t size)
{
    
    return malloc(n * size);
}

void *cs_calloc(int n, size_t size)
{
    /* auto 0 values*/
    return calloc(n, size);
}

void *cs_free(void *p)
{
    if (p) 
    {
        free (p);
    }

    return NULL;
}

void *cs_realloc(void *p, int n, size_t size, int *ok)
{
    void *pnew;
    pnew=realloc(p, n*size);
    *ok=(pnew!=NULL);
    return ((*ok)? pnew:p);
}
cs *cs_spfree(cs *A)
{
    if (!A) return NULL;

    cs_free(A->p);
    cs_free(A->i);
    cs_free(A->x);
    return (cs*)cs_free(A);
}
css *cs_sfree(css *S)
{
    if (!S) return NULL;
    cs_free(S->pinv);
    cs_free(S->q);
    cs_free(S->parent);
    cs_free(S->cp);
    cs_free(S->leftmost);
    return (css *)(cs_free(S));
}
/* allocate a cs_dmperm or cs_scc result */
csd *cs_dalloc (csi m, csi n)
{
    csd *D ;
    D = (csd *)cs_calloc (1, sizeof (csd)) ;
    if (!D) return (NULL) ;
    D->p = (csi *)cs_malloc (m, sizeof (csi)) ;
    D->r = (csi *)cs_malloc (m+6, sizeof (csi)) ;
    D->q = (csi *)cs_malloc (n, sizeof (csi)) ;
    D->s = (csi *)cs_malloc (n+6, sizeof (csi)) ;
    return ((!D->p || !D->r || !D->q || !D->s) ? cs_dfree (D) : D) ;
}
csd *cs_ddone (csd *D, cs *C, void *w, csi ok)
{
    cs_spfree (C) ;                     /* free temporary matrix */
    cs_free (w) ;                       /* free workspace */
    return (ok ? D : cs_dfree (D)) ;    /* return result if OK, else free it */
}
/* free a cs_dmperm or cs_scc result */
csd *cs_dfree (csd *D)
{
    if (!D) return (NULL) ;     /* do nothing if D already NULL */
    cs_free (D->p) ;
    cs_free (D->q) ;
    cs_free (D->r) ;
    cs_free (D->s) ;
    return ((csd *) cs_free (D)) ;  /* free the csd struct and return NULL */
}

cs *cs_spalloc(int m , int n, int nzmax, int values, int triplet)
{
    cs *A = (cs*)cs_calloc(1, sizeof(cs));
    if (!A) return NULL;

    A->m=m;
    A->n=n;
    A->nzmax = nzmax;
    A->nz=triplet? 0:-1;
    A->p = (int *)cs_calloc(triplet? nzmax : n+1, sizeof(int));
    A->i = (int *)cs_calloc(nzmax, sizeof(int));
    A->x = values? (double *)cs_malloc(nzmax, sizeof(double)):NULL;
    return (!A->p ||!A->i || (values && !A->x)) ? cs_spfree(A) : A;

}
int cs_sprealloc(cs *A, int nzmax)
{
    int ok, oki, okj=1, okx=1;
    if (!A) return 0;
    /*******************
    * for csc matrix, 
    * the number of entries is last entries' order;
    * ******************/
    if (nzmax<=0) nzmax = CS_CSC(A)? (A->p[A->n]) : A->nz;

    A->i = (int *)cs_realloc(A->i, nzmax, sizeof(int), &oki);
    if (CS_TRIPLET(A)) A->p = (int *)cs_realloc(A->p, nzmax, sizeof(int), &okj);
    if (A->x) A->x = (double *)cs_realloc(A->x, nzmax, sizeof(double), &okx);

    ok = (oki && okj && okx);
    if (ok) A->nzmax=nzmax;
    return ok;

}


int cs_entry(cs *T, int i, int j, double x)
{
    /* this is triplet*/
    if (!T || i<0 || j<0 ) return 0;

    if (T->nz>=T->nzmax && !cs_sprealloc(T, 2*(T->nzmax))) return 0;

    if (T->x) 
    {
        T->x[T->nz]=x;
    }
    else{
        T->x = (double *)cs_calloc(T->nzmax, sizeof(double));
        T->x[T->nz]=x;
    }
    T->i[T->nz]=i;
    T->p[T->nz++]=j;
    T->m = T->m>i+1? T->m : i+1;
    T->n = T->n>i+1? T->n : i+1;
    return 1;


}


cs *cs_done(cs *C, void *w, void *x, int ok)
{
    cs_free(w);
    cs_free(x);
    return (ok? C : cs_spfree(C));
}
int *cs_idone(int *p, cs *C, void *w, int ok)
{
    cs_spfree(C);
    cs_free(w);
    return (ok ? p : (int *)cs_free(p));
}

csn *cs_nfree(csn *N){
    if (!N) return NULL;
    cs_spfree(N->L);
    cs_spfree(N->U);
    cs_free(N->pinv);
    cs_free(N->B);
    return (csn *)(cs_free(N));
}
csn *cs_ndone(csn *N, cs *C, void *w, void *x, int ok){
    cs_spfree(C);
    cs_free(w);
    cs_free(x);
    return (ok? N : cs_nfree(N));
}
double cs_cumsum(int *p, int *c, int n)
{
    int i, nz=0;
    double nz2=0;
    if (!p || !c) return -1;
    for (i=0;i<n;i++)
    {
        p[i]=nz;
        nz+=c[i];
        nz2+=c[i];
        c[i]=p[i];
    }
    p[n]=nz;
    return nz2;
}
cs *cs_compress(const cs *T)
{   
    /*fix order todo*/
    int m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj;
    double *Cx, *Tx;
    cs *C;
    if (!CS_TRIPLET(T)) return NULL;

    m=T->m; n=T->n; 
    Ti= T->i; Tj=T->p;
    Tx=T->x;nz=T->nz;

    C = cs_spalloc(m, n, nz, Tx != NULL, 0);
    w = (int *)cs_calloc(n, sizeof(int));

    if (!C || !w) return (cs_done(C, w, NULL, 0));

    Cp=C->p;Ci=C->i;Cx=C->x;

#if 0
    for (int it=0;it<n;it++)
    {
        cout<<w[it]<<"\t";
    }
    cout<<endl;
    cout<<"***********"<<endl;
#endif    


    for (k=0;k<nz;k++)
    {
        w[Tj[k]]++; /*count number of entries on each columns*/
    }
    cs_cumsum(Cp, w, n); /*set column pointer*/

    /*the CSC matrix may have same entries and they are required to be summed*/
    for (k=0; k<nz; k++)
    {
        //cout<<k<<"  "<<Ti[k]<<"  "<<Tj[k]<<"  "<<w[Tj[k]]<<endl;
        Ci[p = w[Tj[k]]++] = Ti[k]; /* triplet entries are not in order*/
        if (Cx) Cx[p] = Tx[k];
    }
    return (cs_done(C, w, NULL, 1));
}



int cs_dupl(cs *A)
{
    /**************************************
    * triplet matrix is not convert from original matrix.
    * it allows duplicates entries.
    * therefore, we must sum values with same entries in csc matrix.
    ***************************************/
    int i, j, p, q, nz=0, n, m, *Ap, *Ai, *w;
    double *Ax;
    if (!CS_CSC(A)) return 0;
    m = A->m; 
    n=A->n; Ap=A->p; Ai=A->i; 
    Ax=A->x;
    w=(int *)cs_malloc(m, sizeof(int));
    if (!w) return 0;
    for (i=0;i<m;i++) w[i]=-1;
    for (j=0;j<n;j++)
    {
        q = nz;
        for (p=Ap[j]; p<Ap[j+1];p++)
        {
            /******************************************* 
            * q record new index of entry on each columns
            * p record old index of entry on each columns
            * here is iteration of row on each columns
            * *******************************************/

            i = Ai[p];
            /* if new index of entry */
            if (w[i]>=q)
            {
                Ax[w[i]]+=Ax[p];/*add duplicates entries to first one*/
            }
            else
            {
                w[i]=nz;
                Ai[nz]=i;
                Ax[nz]=Ax[p];
                nz++;
            }
        }
        Ap[j]=q; /*start row index ptr point to the new begin of row*/
    }
    Ap[n]=nz;
    cs_free(w);
    return (cs_sprealloc(A, 0));

}

int cs_fkeep(cs *A, int (*fkeep)(int, int, double, void *), void *other)
{
    /* ********************************
     * remove some entries in cs matrix, 
     * delete extra mem from the back of matrix;
     * *********************************/
    int j, p, nz=0, n, *Ap, *Ai;
    double *Ax;
    if (!CS_CSC(A) || !fkeep) return -1;
    n = A->n; Ap=A->p; Ai=A->i;Ax=A->x;
    for (j=0; j<n;j++)
    {
        p = Ap[j];   /* record original row start prt*/
        Ap[j] = nz;  /* give column ptr new row start index*/

        for ( ; p<Ap[j+1]; p++)
        {
            if (fkeep(Ai[p], j, Ax? Ax[p]:1, other))
            {
                /*old ptr to new ptr*/
                if (Ax) 
                {
                    Ax[nz] = Ax[p];
                }
                Ai[nz] = Ai[p];
                nz++; /*re index new entries, move entries in back to front*/
            }   
        }
    }

    Ap[n]=nz;
    cs_sprealloc(A, 0);
    return nz;
}

static int cs_nonzero(int i, int j, double aij, void *other)
{
    return (aij!=0);
}

int cs_dropzeros(cs *A)
{
    return (cs_fkeep(A, &cs_nonzero, NULL));
}


cs *cs_load(FILE *f)
{
    int i, j;
    double x;
    cs *T;
    if (!f) NULL;
    T = cs_spalloc(0, 0, 1, 1, 1);
    while (fscanf(f, "%d%d%lg\n", &i, &j, &x)==3)
    {
        if (!cs_entry(T, i, j, x)) return cs_spfree(T);
    }
    return T;
}
int cs_print(const cs *A, int brief)
{
    int p, j, m , n, nzmax, nz, *Ap, *Ai;
    double *Ax;

    if (!A) {
        printf(("(null)\n"));
        return 0;
    }

    m = A->m; n = A->n; Ap=A->p; Ai = A->i; Ax = A->x;
    nzmax = A->nzmax; nz=A->nz;

    printf("CSparse Version %d %d %d %s, %s\n",
        CS_VER, CS_SUBVER, CS_SUBSUB, CS_DATE, CS_COPYRIGHT);

    if (nz<0)
    {
        printf("%d-by-%d, nzmax: %d nnz: %d, 1-norm: %g\n",
            m, n, nzmax, Ap[n], cs_norm(A));
        for (j = 0; j<n;j++)
        {
            printf("    col %d : locations %d to %d\n",
                j, Ap[j], Ap[j+1]-1);
            for (p = Ap[j]; p<Ap[j+1];p++)
            {
                printf("    %d : %g\n", Ai[p], Ax ? Ax[p] : 1);
                if (brief && p > 20)
                {
                    printf("    ...\n");
                    return 1;
                }
            }
        }
        
    }
    else{
            printf("triplet: %d-by-%d, nzmax: %d nnz: %d\n", 
                m, n, nzmax, nz);
            for (p=0;p<nz;p++)
            {
                printf("    %d %d : %g\n", Ai[p], Ap[p], Ax? Ax[p] : 1);
                if (brief &&p>20)
                {
                    printf("    ...\n");
                    return 1;
                }
            }   
    }
    return 1;

    

}

css *cs_schol(int order, const cs *A){
    int n, *c, *post, *P;
    cs *C;
    css *S;
    if (!CS_CSC(A)) return NULL;
    n = A->n;
    S = (css *)cs_calloc(1, sizeof(css));

    if (!S) return NULL;
    //P = cs_amd(order, A);   /*find a permutation P so that PAP^T has fewer nonzero in its factorization than A*/
    P = NULL;//todo
    S->pinv = cs_pinv(P, n);
    cs_free(P);
    if (order && !S->pinv) return (cs_sfree(S));
    C = cs_symperm(A, S->pinv, 0);
    S->parent = cs_etree(C, 0);
    post = cs_post(S->parent, n);
    c = cs_counts(C, S->parent, post, 0);
    cs_free(post);
    cs_spfree(C);
    S->cp = (int *)cs_malloc(n+1, sizeof(int));
    S->unz = S->lnz = cs_cumsum(S->cp, c, n);
    cs_free(c);
    return (S->lnz>=0? S: cs_sfree(S));
}