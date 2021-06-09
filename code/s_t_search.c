
#include <GraphBLAS.h>
#include <stdio.h>
#include "readmtx.h"
#include "s_t_bfs.h"

int main(int argc, char** argv)
{
    GrB_init ( GrB_NONBLOCKING );
    GrB_Index n;            //No. of vertices
    GrB_Index edges;        //No. of edges
    GrB_Index* row_indeces;
    GrB_Index* col_indeces;
    double* values;
    GrB_Matrix A = NULL;    //Input graph
    GrB_Index s = 0; //Default source
    GrB_Index t = -1;


    readMtx(argv[1], &n, &edges, &row_indeces, &col_indeces, &values);
    CHECK( GrB_Matrix_new (&A, GrB_FP64, n, n) );
    CHECK( GrB_Matrix_build (A, row_indeces, col_indeces, values, edges, GrB_PLUS_FP64) );

    s_t_bfs(A,&s,&t);

    GrB_finalize ( ) ;
}