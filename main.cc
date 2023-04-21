#include<iostream>
#include"cs.h"
#include"cs_op.h"
#include<vector>
#include<stdlib.h>
#include<limits.h>
#include<stdio.h>
using namespace std;
void show_matrix(vector<vector<double>> m)
{
    for (int i=0;i<m.size();i++)
    {
        for (int j=0;j<m[i].size();j++)
        {
            cout<<m[i][j]<<'\t';
        }
        cout<<endl;
    }
}
int get_nzmax(vector<vector<double>> m)
{
    int nz=0;
    for (int i=0;i<m.size();i++)
    {
        for (int j=0;j<m[i].size();j++)
        {
            if(m[i][j]!=0)
                nz++;
        }
        
    }
    return nz;
}

void convert2triplet(cs * A, vector<vector<double>> m)
{
    for (int i=0;i<m.size();i++)
    {
        for (int j=0;j<m[i].size();j++)
        {
            if(m[i][j]!=0)
                cs_entry(A, i, j, m[i][j]);
        }
        
    }

}

cs *get_csMatrix(vector<vector<double>> rawMatrix, size_t dim)
{
    /*generate triplet matrix*/
    int nzmax = get_nzmax(rawMatrix);
    cs *triMatrix = cs_spalloc(dim, dim, nzmax, 0, 1);
    // show_cs_details(triMatrix);
    /* convert triplet matrix to csc matrix*/
    convert2triplet(triMatrix, rawMatrix);
    show_cs_details(triMatrix);
    cs *csMatrix = cs_compress(triMatrix);
    
    return csMatrix;
}
int main()
{
   
   
    
    FILE *fp;
    fp = fopen("qr_tri", "r");
    cs *triMatrix = cs_load(fp);
    cs *A = cs_compress(triMatrix);
    cs *AT = cs_transpose(A, 1);
    cs *ATA = cs_multiply(AT, A);
    A = cs_sort(A);
    css *S = cs_sqr(0, A, 1);
    csn *SN = cs_qr(A, S);

    // cs_spfree(A);
    // A = AT;

    cout<<A->nzmax<<"  "<<ATA->nzmax<<endl;
    int *post, *parent;
    parent = cs_etree(A, 1);
    post = cs_post(parent, A->n);
    // int n = A->n;
    // for (int it=0;it<n;it++)
    //     cout<<it<<'\t';
    // cout<<endl;
    // for (int it=0;it<n;it++)
    //     cout<<parent[it]<<'\t';
    // cout<<endl;

    // cs_free(parent);
    // cs_free(post);
    // parent = cs_etree(ATA, 0);
    // post = cs_post(parent, ATA->n);
    
    // for (int it=0;it<n;it++)
    //     cout<<it<<'\t';
    // cout<<endl;
    // for (int it=0;it<n;it++)
    //     cout<<parent[it]<<'\t';
    // cout<<endl;


    // for (int it=0;it<n;it++)
    //     cout<<post[it]<<'\t';
    // cout<<endl;
    
    cs_counts(A, parent, post, 1);
    

    return 0;
}