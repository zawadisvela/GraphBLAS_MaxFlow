
#include <GraphBLAS.h>
#include <stdio.h>
#include "readmtx.h"

#ifdef DEBUG
#define CHECK(x) \
{   \
    int err_code = x; \
    if (err_code != 0) { \
        printf("ERROR! %s:%d:%s - errcode: %d\n",__FILE__, __LINE__, __func__, err_code); \
    } \
}
#else
#define CHECK(x) x
#endif

#define NO_MASK GxB_DEFAULT
#define NO_ACCUM GxB_DEFAULT
#define DEFAULT_DESC GxB_DEFAULT

#define DOUBLE_PREC 1


double delta_f_global = 0 ;
void apply_delta (void *result, const void *element)
{
    (* ((double *) result)) = delta_f_global ;
}

//Perform multiple simultaneous BFS searches to determine s-t
void set_s_t_bfs(
    const GrB_Matrix R,
    GrB_Index *s,
    GrB_Index *t
)
{
    /*
    GrB_Index n;
    GrB_Matrix level = NULL;
    GrB_Matrix frontiers;
    GrB_Descriptor desc = NULL;

    CHECK( GrB_Matrix_nrows(&n, R) );

    CHECK( GrB_Matrixr_new(&level, GrB_INT32, n, n) );

    GrB_Vector_setElement(level, 0, s);

    GrB_Vector_new(&frontier, GrB_BOOL, n);
    GrB_Vector_setElement(frontier, true, s);

    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_COMP);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

    bool successor = true;
    while (successor) {
        GrB_vxm(frontier, reachable, NO_ACCUM, GxB_LOR_LAND_BOOL, frontier, R, desc);
        GrB_reduce(&successor, NO_ACCUM, GrB_LOR_MONOID_BOOL, frontier, DEFAULT_DESC);
        GrB_assign(reachable, frontier, NULL, true, GrB_ALL, n, NULL);
    }
    */

}

void mincut_bfs(
    GrB_Vector* result,
    const GrB_Matrix R,
    GrB_Index s
)
{
    GrB_Index n;
    CHECK( GrB_Matrix_nrows(&n, R) );

    GrB_Vector reachable = *result;
    CHECK( GrB_Vector_setElement(reachable, true, s) );

    GrB_Vector frontier;
    CHECK( GrB_Vector_new(&frontier, GrB_BOOL, n) );
    CHECK( GrB_Vector_setElement(frontier, true, s) );

    GrB_Descriptor desc = NULL;
    CHECK( GrB_Descriptor_new(&desc) );
    CHECK( GrB_Descriptor_set(desc, GrB_MASK, GrB_COMP) );
    CHECK( GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE) );

    bool successor = true;
    while (successor) {
        CHECK( GrB_vxm(frontier, reachable, NO_ACCUM, GxB_LOR_LAND_BOOL, frontier, R, desc) );
        CHECK( GrB_reduce(&successor, NO_ACCUM, GrB_LOR_MONOID_BOOL, frontier, DEFAULT_DESC) );
        CHECK( GrB_assign(reachable, frontier, NULL, true, GrB_ALL, n, NULL) );
    }
}

bool get_augmenting_path(
//        GrB_Matrix  const graph,
        GrB_Matrix  const R,
        GrB_Index   source,
        GrB_Index   sink,
        GrB_Vector  const index_ramp,
        GrB_Matrix  M
)
{
    GrB_Index n;
    CHECK( GrB_Matrix_nrows(&n, R) );

/*
    GrB_Matrix R = NULL;
    CHECK( GrB_Matrix_dup(&R, graph) );
*/
    GrB_Vector parent_list = NULL;
    CHECK( GrB_Vector_new(&parent_list, GrB_INT32, n) );
    CHECK( GrB_Vector_setElement(parent_list, source, source) );

    GrB_Vector wavefront = NULL;
    CHECK( GrB_Vector_new(&wavefront, GrB_INT32, n) );
    CHECK( GrB_Vector_setElement(wavefront, 0xDEADBEEF, source) ); //Value doesn't matter, only existence

    int sink_parent;
    GrB_Index wavefront_nvals;
    CHECK( GrB_Vector_nvals(&wavefront_nvals, wavefront) );
    //while ((!parent_list.hasElement(sink)) && (wavefront.nvals() > 0))

    GrB_Descriptor desc = NULL ;           // Descriptor for vxm
    CHECK( GrB_Descriptor_new (&desc) );
    CHECK( GrB_Descriptor_set (desc, GrB_MASK, GrB_COMP) );     // invert the mask
    CHECK( GrB_Descriptor_set (desc, GrB_MASK, GrB_STRUCTURE) );     // Don't care about values, only structure of mask
    CHECK( GrB_Descriptor_set (desc, GrB_OUTP, GrB_REPLACE) );  // clear q first

    while (((GrB_Vector_extractElement(&sink_parent, parent_list, sink) == GrB_NO_VALUE)) && (wavefront_nvals > 0))
    {

        //printf("\nCurrent frontier:\n");
        //GxB_print(wavefront, GxB_SHORT);

        // convert all stored values to their column index
        CHECK( GrB_eWiseMult(wavefront,
            NO_MASK, NO_ACCUM, //mask, accumulate op
            GrB_FIRST_INT32,
            index_ramp, wavefront,
            DEFAULT_DESC) ); //descriptor

        ////GxB_print(wavefront, GxB_SHORT);

        // First because we are left multiplying wavefront rows
        // Masking out the parent list ensures wavefront values do not
        // overlap values already stored in the parent list
/*
        printf("VxM\n");
        printf("w\n:");
        //GxB_print(wavefront, GxB_SHORT);
        printf("mask\n:");
        //GxB_print(parent_list, GxB_SHORT);
        printf("u\n:");
        //GxB_print(wavefront, GxB_SHORT);
        printf("A\n:");
        //GxB_print(R, GxB_SHORT);
        printf("desc\n:");
        //GxB_print(desc, GxB_SHORT);
*/
        CHECK( GrB_vxm(wavefront, //c
            //NULL,
            parent_list, //mask <- automatically removes self-loops
            NO_ACCUM, //no accumulate/default
#if DOUBLE_PREC
            //GxB_MIN_FIRST_FP64, //semiring. First to extract parent. Why min? Could be ANY?
            GxB_ANY_FIRST_FP64, //Should yield same result, and leaves more implementation freedom
#endif
            wavefront, // v
            R, // M
            desc) ); //descriptor

        //printf("RESULT w\n:");
        ////GxB_print(wavefront, GxB_SHORT);
        printf("\npre-apply parent list\n");
        GxB_print(parent_list, GxB_SHORT);
        CHECK( GrB_Vector_apply(parent_list, //w
            NO_MASK, //mask
            GrB_PLUS_INT32, //accum
            GrB_IDENTITY_INT32, //unary op. w<M> = accum (w, op (u)). Nor sure this identity is the same as for indexes
            wavefront, //u
            DEFAULT_DESC) ); // descriptor
        printf("post-apply parent list\n\n");
        GxB_print(parent_list, GxB_SHORT);
/*

        CHECK( GrB_eWiseAdd(parent_list, //w
            NO_MASK, //mask
            GrB_PLUS_INT32, //accum
            GrB_IDENTITY_INT32, //unary op. w<M> = accum (w, op (u)). Nor sure this identity is the same as for indexes
            wavefront, //u
            DEFAULT_DESC) ); // descriptor
*/
        CHECK( GrB_Vector_nvals(&wavefront_nvals, wavefront) );

        //printf("\nparent_list:\n");
        ////GxB_print(parent_list, GxB_SHORT);
    }

    if ((GrB_Vector_extractElement(&sink_parent, parent_list, sink) == GrB_NO_VALUE))
    {
        return false;
    }

    // Extract path from source to sink from parent list (reverse traverse)
    // build a mask
//    M.clear();

    CHECK( GrB_Matrix_clear(M) );
    GrB_Index curr_vertex = sink;
    while (curr_vertex != source)
    {
    //    grb::IndexType parent(parent_list.extractElement(curr_vertex));
        GrB_Index parent;
        CHECK( GrB_Vector_extractElement(&parent, parent_list, curr_vertex) );

        //M.setElement(parent, curr_vertex, true);
        CHECK( GrB_Matrix_setElement(M, true, parent, curr_vertex) );

        curr_vertex = parent;
    }

    return true;
}



int main (int argc, char **argv)
{
    if(argc < 4) {
        printf("Usage: ./edmund-karp filename.mtx <s> <t>\n");
        return 1;
    }
    GrB_init ( GrB_NONBLOCKING );
    //Blocking better for debuggbing
    //GrB_init ( GrB_BLOCKING ) ;



    GrB_Index n;            //No. of vertices
    GrB_Index edges;        //No. of edges
    GrB_Index* row_indeces;
    GrB_Index* col_indeces;
    double* values;
    GrB_Matrix A = NULL;    //Input graph
    GrB_Matrix M = NULL;    //Mask for path
    GrB_Matrix P = NULL;    //Weighted path
    GrB_Matrix R = NULL;    //Residual network
    GrB_Vector index_ramp;  //v[i] = i. For extracting parent in BFS
    GrB_Index source = -1;
    GrB_Index sink = -1;

    readMtx(argv[1], &n, &edges, &row_indeces, &col_indeces, &values);
    CHECK( GrB_Matrix_new (&A, GrB_FP64, n, n) );
    CHECK( GrB_Matrix_build (A, row_indeces, col_indeces, values, edges, GrB_PLUS_FP64) );
#ifdef DEBUG
    GxB_print(A, GxB_SHORT);
#endif
    GrB_Matrix_new(&M, GrB_BOOL, n, n);
    GrB_Matrix_new(&P, GrB_FP64, n, n);

    CHECK( GrB_Matrix_dup(&R, A) );

    CHECK( GrB_Vector_new(&index_ramp, GrB_INT32, n) );
    for (int i = 0; i < n; i++)
        CHECK( GrB_Vector_setElement(index_ramp, i, i) );
#ifdef DEBUG
    printf("Index ramp:\n");
    GxB_print(index_ramp, GxB_SHORT);
#endif
    //printf("Success? %s\n", get_augmenting_path(A, 0, 5, M)?"yes! :D":"no :(");

    source = atoi(argv[2]);
    sink = atoi(argv[3]); //TODO set sink to something sensible
    printf("SOURCE:%ld, SINK:%ld\n", source, sink);
    if (sink > n) {
        printf("Out of bounds sink y'all ERROR to the MAX!\n");
        return 1;
    }

    int count = 0;
    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, A);
    while (get_augmenting_path(R, source, sink, index_ramp, M) && count++ < nvals)
    {

        printf("----------- Iteration: %d -----------\n", count);
        //sprintf("\nPath M:\n");
        //GxB_print(M, GxB_SHORT);

        GrB_eWiseMult(P, NO_MASK, NO_ACCUM, GrB_TIMES_FP64, M, R, DEFAULT_DESC);
#ifdef DEBUG
        printf("\nPath:\n");
        GxB_print(P, GxB_SHORT);
#endif
        //GrB_reduce(scalar, accum op, reduction monoid, vector/matrix, descriptor)
        //matrix to vector reduction also has a mask argument, not optional, but can be NULL
        CHECK(GrB_reduce(&delta_f_global, NO_ACCUM, GrB_MIN_MONOID_FP64, P, DEFAULT_DESC));
        //printf("Gamma value; %lf\n", delta_f_global);

        GrB_UnaryOp apply_delta_op;
        CHECK( GrB_UnaryOp_new(&apply_delta_op, apply_delta, GrB_FP64, GrB_BOOL) );

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

    double total_flow = 0;
    for (GrB_Index i = 0; i < n; i++) {
        double capacity = 0;
        double residual = 0;

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
