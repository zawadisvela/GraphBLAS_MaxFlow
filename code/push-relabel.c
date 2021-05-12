
#include <GraphBLAS.h>
#include <stdio.h>
#include "readmtx.h"


#define CHECK(x) \
{   \
    int err_code = x; \
    if (err_code != 0) { \
        printf("ERROR! %s:%d:%s - errcode: %d\n",__FILE__, __LINE__, __func__, err_code); \
    } \
}

#define NO_MASK GxB_DEFAULT
#define NO_ACCUM GxB_DEFAULT
#define DEFAULT_DESC GxB_DEFAULT

#define DOUBLE_PREC 1

int main (int argc, char **argv)
{
    if(argc < 4) {
        printf("Usage: ./edmund-karp filename.mtx <s> <t>\n");
        return 1;
    }
    //GrB_init ( GrB_NONBLOCKING )
    //Blocking better for debuggbing
    GrB_init ( GrB_BLOCKING ) ;



    GrB_Index n;
    GrB_Index edges;
    GrB_Index* row_indeces;
    GrB_Index* col_indeces;
    double* values;
    GrB_Matrix A = NULL;

    readMtx(argv[1], &n, &edges, &row_indeces, &col_indeces, &values);

    //Old hard-coded matrix
    /*
    n = 6;
    edges = 6;
    GrB_Index row_indeces[6] = {0,0,1,1,3,2} ;
    GrB_Index col_indeces[6] = {1,2,3,2,2,5} ;
    double values[6] = {4,2,3,1,5,7};
    */

    CHECK( GrB_Matrix_new (&A, GrB_FP64, n, n) );
    GxB_print(A, GxB_SHORT);

    CHECK( GrB_Matrix_build (A, row_indeces, col_indeces, values, edges, GrB_PLUS_FP64) );

    GxB_print(A, GxB_SHORT);

    GrB_Matrix M = NULL;
    GrB_Matrix_new(&M, GrB_BOOL, n, n);

    GrB_Index source = atoi(argv[2]), sink = atoi(argv[3]); //TODO set sink to something sensible
    printf("SOURCE:%ld, SINK:%ld\n", source, sink);
    if (sink > n) {
        printf("Out of bounds sink y'all ERROR to the MAX!\n");
        return 1;
    }

    GrB_Matrix R = NULL;
    CHECK( GrB_Matrix_dup(&R, A) );

    GrB_Vector index_ramp;
    CHECK( GrB_Vector_new(&index_ramp, GrB_INT32, n) );
    for (int i = 0; i < n; i++)
        CHECK( GrB_Vector_setElement(index_ramp, i, i) );
    printf("Index ramp:\n");
    //GxB_print(index_ramp, GxB_SHORT);

    //printf("Success? %s\n", get_augmenting_path(A, 0, 5, M)?"yes! :D":"no :(");

// Unsure if needed. Maybe some edge cases?
    int count = 0;
    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, A);
    while (get_augmenting_path(R, source, sink, index_ramp, M) && count++ < nvals)
    {
        printf("----------- Iteration: %d -----------\n", count);
        //sprintf("\nPath M:\n");
        //GxB_print(M, GxB_SHORT);

        GrB_Matrix P;
#if DOUBLE_PREC
        GrB_Matrix_new(&P, GrB_FP64, n, n);
        GrB_eWiseMult(P, NO_MASK, NO_ACCUM, GrB_TIMES_FP64, M, R, DEFAULT_DESC);
#endif
        printf("\nPath:\n");
        GxB_print(P, GxB_SHORT);
        GrB_Vector test;
        GrB_Vector_new(&test, GrB_FP64, n);
        //GrB_reduce(scalar, accum op, reduction monoid, vector/matrix, descriptor)
        //matrix to vector reduction also has a mask argument, not optional, but can be NULL
        CHECK(GrB_reduce(&delta_f_global, NO_ACCUM, GrB_MIN_MONOID_FP64, P, DEFAULT_DESC));
        //printf("Gamma value; %lf\n", delta_f_global);

        GrB_UnaryOp apply_delta_op;
#if DOUBLE_PREC
        CHECK( GrB_UnaryOp_new(&apply_delta_op, apply_delta, GrB_FP64, GrB_BOOL) );
#endif

        //GxB_print(M, GxB_SHORT);

        GrB_Descriptor transpose_a;
        GrB_Descriptor_new(&transpose_a);
        GrB_Descriptor_set(transpose_a, GrB_INP0, GrB_TRAN);

        CHECK( GrB_Matrix_apply(P, NO_MASK, NO_ACCUM, apply_delta_op, M, transpose_a) );
        //GxB_print(P, GxB_SHORT);
        CHECK( GrB_Matrix_apply(P, NO_MASK, GrB_PLUS_FP64, GrB_AINV_FP64, P, transpose_a) );
        //GxB_print(P, GxB_SHORT);

        //GxB_print(R, GxB_SHORT);
        GrB_eWiseAdd(R, NO_MASK, NO_ACCUM, GxB_PLUS_FP64_MONOID, R, P, DEFAULT_DESC);
        //GxB_print(R, GxB_SHORT);

        GrB_Descriptor replace;
        GrB_Descriptor_new(&replace);
        GrB_Descriptor_set(replace, GrB_OUTP, GrB_REPLACE);

        GrB_apply(R, R, NO_ACCUM, GrB_IDENTITY_FP64, R, replace); //Remove zero-edges
        //GxB_print(R, GxB_SHORT);
    }

    double total_flow;
    for (GrB_Index i = 0; i < n; i++) {
        double capacity;
        double residual;

        if(GrB_Matrix_extractElement(&capacity, A, source, i) == GrB_NO_VALUE)
            capacity = 0;
        if(GrB_Matrix_extractElement(&residual, R, source, i) == GrB_NO_VALUE)
            residual = 0;
        //printf("Found edge (%ld,%ld) w cap: %lf and residual: %lf. Adding %lf to total flow\n", source, i, capacity, residual, capacity-residual);
        total_flow += capacity-residual;
    }
    //Could use elementwise multiplication between graph and R to figure out the min-cut?
    //Intersection between graph and residuals complement?
    GrB_Descriptor mask_complement;
    GrB_Descriptor_new(&mask_complement);
    GrB_Descriptor_set(mask_complement, GrB_MASK, GrB_COMP);

    GrB_Matrix min_cut;
    GrB_Matrix_new(&min_cut, GrB_FP64, n, n);
    GrB_apply(min_cut, R, NO_ACCUM, GrB_IDENTITY_FP64, A, mask_complement);

    GrB_Vector reachable;
    GrB_Vector_new(&reachable, GrB_BOOL, n);
    mincut_bfs(&reachable, R, source);

    printf("Min-cut edges:\n");
    double val;
    double min_cut_value = 0;
    for (GrB_Index i = 0; i < n; i++) {
        if(GrB_Vector_extractElement(&val, reachable, i) != GrB_NO_VALUE) {
            for (GrB_Index j = 0; j < n; j++) {
                if(GrB_Matrix_extractElement(&val, min_cut, i, j) != GrB_NO_VALUE) {
                    printf("(%ld,%ld)\t%lf\n", i, j, val);
                    min_cut_value += val;
                }
            }
        }
    }
    //Are there edge-cases where flow will go back into source? No, becasue BFS always begins in the source
    printf("Max flow: %lf\n", total_flow);
    printf("Min cut: %lf\n", min_cut_value);
    GrB_finalize ( ) ;
}
