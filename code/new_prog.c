

#define CHECK(x) \
{   \
    int err_code = x; \
    if (err_code != 0) \
        printf("Error! %s:%d:%s - errcode: %d\n",__FILE__, __LINE__, __func__, err_code); \
}

#include <GraphBLAS.h>
#include <stdio.h>

int main (int argc, char **argv)
{
    GrB_init ( GrB_NONBLOCKING ) ;
    //GrB_init ( GrB_BLOCKING ) ;

    printf("Hello world\n");
    int64_t nrows = 4 ;
    int64_t ncols = 4 ;
    GrB_Matrix C = NULL ;
    //GrB_Index *I = (GrB_Index *) malloc (len * sizeof (GrB_Index));
    //GrB_Index *J = (GrB_Index *) malloc (len * sizeof (GrB_Index));
    //double *X = malloc (len * sizeof(double)) ;
    size_t ntuples = 6;
    GrB_Index I[6] = {0,0,2,2,2,3} ;
    GrB_Index J[6] = {1,2,0,1,3,1} ;
    int8_t X[6] = {4,1,3,2,5,6} ;
    GrB_Type xtype = GrB_INT8 ;
    GrB_BinaryOp xop = GrB_PLUS_INT8 ;
    CHECK(GrB_Matrix_new (&C, xtype, nrows, ncols));
    CHECK(GrB_Matrix_build (C, I, J, X, ntuples, xop)) ;
    GxB_print(C, GxB_COMPLETE);

    GrB_Vector Y = NULL;
    GrB_Vector_new(&Y, xtype, nrows);
    //GrB_Vector_setElement(VECTOR, VALUE, INDEX);
    GrB_Vector_setElement(Y, 1, 0);
    GxB_print(Y, GxB_COMPLETE);

    GrB_Vector f = NULL;
    GrB_Vector_new(&f, xtype, nrows);
    //TODO: construct semring
    //GrB_vxm(Y, NULL, NULL, C, )

    GrB_Matrix A;     // input graph, treated as if boolean in semiring
    GrB_Index s = 0;           // starting node of the BFS

    //--------------------------------------------------------------------------
    // set up the semiring and initialize the vector v
    //--------------------------------------------------------------------------

    GrB_Index n = 4;                          // # of nodes in the graph
    GrB_Vector q = NULL ;                  // nodes visited at each level
    GrB_Vector v = NULL ;                  // result vector
    GrB_Descriptor desc = NULL ;           // Descriptor for vxm

    CHECK(GrB_Matrix_new (&A, GrB_INT8, n, n));
    CHECK(GrB_Matrix_build (A, I, J, X, ntuples, xop)) ;

    CHECK (GrB_Vector_new (&v, GrB_INT32, n)) ;    // q<int32_t> v(n) = 0

    //These two lines just force a fully instantiated dense vector v. why?
    CHECK (GrB_Vector_assign_INT32 (v, NULL, NULL, 0, GrB_ALL, n, NULL) );// make dense. Zawadi: Why??
    GrB_Vector_nvals (&n, v) ;             // finish pending work on v

    GrB_Vector_new (&q, GrB_BOOL, n) ;     // Vector<bool> q(n) = false
    GrB_Vector_setElement_BOOL (q, true, s) ;   // q[s] = true, false elsewhere

    GrB_Descriptor_new (&desc) ;
    GrB_Descriptor_set (desc, GrB_MASK, GrB_COMP) ;     // invert the mask
    GrB_Descriptor_set (desc, GrB_OUTP, GrB_REPLACE) ;  // clear q first

    printf("after setu-up print\n");

    GxB_print(v, GxB_COMPLETE);

    GrB_Semiring any_and_semiring = NULL;
    GrB_Semiring_new(&any_and_semiring, GxB_ANY_BOOL_MONOID, GxB_LAND_BOOL);
    //--------------------------------------------------------------------------
    // BFS traversal and label the nodes
    //--------------------------------------------------------------------------

    bool successor = true ; // true when some successor found
    for (int32_t level = 1 ; successor && level <= n ; level++)
    {
        // v<q> = level, using vector assign with q as the mask
        GrB_Vector_assign_INT32 (v, q, NULL, level, GrB_ALL, n, NULL) ;
        printf("round %d:\n", level);

        GxB_print(v, GxB_COMPLETE);
        // q<!v> = q ||.&& A ; finds all the unvisited
        // successors from current q, using !v as the mask
        GrB_vxm (q, v, NULL, any_and_semiring, q, A, desc) ;

        // successor = ||(q)
        GrB_Vector_reduce_BOOL (&successor, NULL, GrB_LOR_MONOID_BOOL, q, NULL) ;
    }

    // make v sparse
    GrB_Descriptor_set (desc, GrB_MASK, GxB_DEFAULT) ;  // mask not inverted
    GrB_Vector_assign (v, v, NULL, v, GrB_ALL, n, desc) ;

    GxB_print(v, GxB_COMPLETE);


    GrB_finalize ( ) ;
}
