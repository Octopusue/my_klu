#include<iostream>
#include"cs.h"
#include"cs_op.h"
#include<vector>
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
    size_t dim=4;
    vector<vector<double>> rawMatrix;
    vector<double> m1 = {4.5, 0, 3.2, 0}, m2={3.1, 2.9, 0, 0.9}, m3={0, 1.7, 3.0, 0}, m4{3.5, 0.4, 0, 1.0};
    rawMatrix.push_back(m1);rawMatrix.push_back(m2);rawMatrix.push_back(m3);rawMatrix.push_back(m4);

    show_matrix(rawMatrix);

    cs *csMatrix1 = get_csMatrix(rawMatrix, dim);
    cs *csMatrix2 = get_csMatrix(rawMatrix, dim);

    csMatrix2 = cs_transpose(csMatrix2, 1);

    show_cs_details(csMatrix1);
    show_cs_details(csMatrix2);

    cs_print(csMatrix1, 1);

    // /*add duplicates entries*/
    // cs_entry(triMatrix, 2, 2, 0.05);
    // show_cs_details(triMatrix);
    // cs_spfree(csMatrix);
    // csMatrix = cs_compress(triMatrix);
    // show_cs_details(csMatrix);

    return 0;
}