
#include <GraphBLAS.h>
#include <stdio.h>
#include <omp.h>
//#include <time.h>
#include "readmtx.h"
#include "s_t_bfs.h"
#include <sys/time.h>
#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)


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


static int runs = -1;
static int run = 0;

void get_mincut(
    GrB_Matrix *C,
    const GrB_Matrix A,
    const GrB_Matrix R,
    GrB_Index s
)
{
#ifdef DEBUG
    printf(">>>Beginning get_mincut()\n");
    GxB_print(A, GxB_SHORT);
    GxB_print(R, GxB_SHORT);
#endif
    GrB_Vector reachable;
    GrB_Vector frontier;
    GrB_Descriptor desc = NULL;
    GrB_Descriptor transpose_b = NULL;
    GrB_Vector s_cut;
    GrB_Vector t_cut;
    GrB_Matrix cut_extract = NULL;

    GrB_Index n;
    CHECK( GrB_Matrix_nrows(&n, A) );

    GrB_Vector_new(&reachable, GrB_INT32, n);
    CHECK( GrB_Vector_setElement(reachable, 1, s) );

    CHECK( GrB_Vector_new(&frontier, GrB_BOOL, n) );
    CHECK( GrB_Vector_setElement(frontier, true, s) );

    CHECK( GrB_Descriptor_new(&desc) );
    CHECK( GrB_Descriptor_set(desc, GrB_MASK, GrB_COMP) );
    CHECK( GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE) );

    bool successor = true;
    while (successor) {
    //for (int i = 0; i < 3; i++) {
        //GxB_print(frontier, GxB_SHORT);
        CHECK( GrB_vxm(frontier, reachable, NO_ACCUM, GxB_LOR_LAND_BOOL, frontier, R, desc) );
        CHECK( GrB_reduce(&successor, NO_ACCUM, GrB_LOR_MONOID_BOOL, frontier, DEFAULT_DESC) );
        CHECK( GrB_assign(reachable, frontier, NULL, true, GrB_ALL, n, NULL) );
    }

    #ifdef DEBUG
        GxB_print(reachable, GxB_SHORT);
    #endif

    GrB_Vector_new(&t_cut, GrB_INT32, n);
    CHECK( GrB_vxm(t_cut, reachable, NO_ACCUM, GxB_ANY_FIRSTJ_INT32, reachable, A, desc) );

    CHECK( GrB_Descriptor_new(&transpose_b) );
    CHECK( GrB_Descriptor_set(transpose_b, GrB_INP1, GrB_TRAN) );

    GrB_Vector_new(&s_cut, GrB_INT32, n);
    CHECK( GrB_vxm(s_cut, reachable, NO_ACCUM, GxB_ANY_FIRSTJ_INT32, t_cut, A, transpose_b) );

#ifdef DEBUG
    GxB_print(s_cut, GxB_SHORT);
    GxB_print(t_cut, GxB_SHORT);
#endif

    GrB_Index s_len  = -1;
    CHECK( GrB_Vector_nvals(&s_len, s_cut) );
    GrB_Index *s_indices = malloc(s_len * sizeof(GrB_Index));
    GrB_Index *s_vals = malloc(s_len * sizeof(GrB_Index));

    GrB_Index t_len  = -1;
    CHECK( GrB_Vector_nvals(&t_len, t_cut) );
    GrB_Index *t_indices = malloc(t_len * sizeof(GrB_Index));
    GrB_Index *t_vals = malloc(t_len * sizeof(GrB_Index));

    CHECK( GrB_Vector_extractTuples(s_indices, s_vals, &s_len, s_cut) );
    CHECK( GrB_Vector_extractTuples(t_indices, t_vals, &t_len, t_cut) );

    GrB_Matrix_new(&cut_extract, GrB_FP64, s_len, t_len);

    CHECK( GrB_extract(cut_extract, NO_MASK, NO_ACCUM, A, s_indices, s_len, t_indices, t_len, DEFAULT_DESC) );

    GrB_Index cut_len = -1;
    GrB_Matrix_nvals(&cut_len, cut_extract);

    GrB_Index *cut_s = malloc(cut_len * sizeof(GrB_Index));
    GrB_Index *cut_t = malloc(cut_len * sizeof(GrB_Index));
    double *cut_val = malloc(cut_len * sizeof(GrB_Index));

    GrB_Matrix_extractTuples(cut_s, cut_t, cut_val, &cut_len, cut_extract);

    for(int i = 0; i < cut_len; i++) {
        cut_s[i] = s_indices[cut_s[i]];
        cut_t[i] = t_indices[cut_t[i]];
    }

    GrB_Matrix_new(C, GrB_FP64, n, n);
    CHECK( GrB_Matrix_build (*C, cut_s, cut_t, cut_val, cut_len, GrB_PLUS_FP64) );

    GrB_free(&reachable);
    GrB_free(&frontier);
    GrB_free(&desc);
    GrB_free(&transpose_b);
    GrB_free(&s_cut);
    GrB_free(&t_cut);
    GrB_free(&cut_extract);
    free(s_indices);
    free(s_vals);
    free(t_indices);
    free(t_vals);
    free(cut_s);
    free(cut_t);
    free(cut_val);
    #ifdef DEBUG
        printf(">>>End get_mincut()\n");
    #endif
}

bool get_augmenting_path(
//        GrB_Matrix  const graph,
        GrB_Matrix  const R,
        GrB_Index   source,
        GrB_Index   sink,
        GrB_Matrix  M
)
{
    GrB_Index n;
    GrB_Vector parent_list = NULL;
    GrB_Vector frontier = NULL;
    GrB_Descriptor vxm_desc = NULL ;     // Descriptor for vxm
    int sink_parent = -1;
    GrB_Index frontier_nvals;
    bool found_path = true;

    CHECK( GrB_Matrix_nrows(&n, R) );

    CHECK( GrB_Vector_new(&parent_list, GrB_INT32, n) );
    CHECK( GrB_Vector_new(&frontier, GrB_INT32, n) );

    CHECK( GrB_Descriptor_new (&vxm_desc) );
    CHECK( GrB_Descriptor_set (vxm_desc, GrB_MASK, GrB_COMP) );     // invert the mask
    CHECK( GrB_Descriptor_set (vxm_desc, GrB_MASK, GrB_STRUCTURE) );     // Don't care about values, only structure of mask
    CHECK( GrB_Descriptor_set (vxm_desc, GrB_OUTP, GrB_REPLACE) );  // clear q first

    int depth = 0;
    if (run == runs-1) {
        printf("Starting BFS\n");
    }

    CHECK( GrB_Vector_setElement(frontier, 1, source) ); //Value irrelevant
    CHECK( GrB_Vector_setElement(parent_list, source, source) ); //Source is it's own parent
    CHECK( GrB_Vector_nvals(&frontier_nvals, frontier) );

    while ((frontier_nvals > 0) && (GrB_Vector_extractElement(&sink_parent, parent_list, sink) == GrB_NO_VALUE))
    {
        if (run == runs-1) {
            printf("Depth: %d, frontier size: %ld \n", depth, frontier_nvals);
        }
        depth++;

        #if defined(DEBUG) || defined(PROFILE) || defined(DETERMINISTIC)
        CHECK( GrB_vxm(frontier, parent_list, NO_ACCUM, GxB_MIN_FIRSTJ_INT32, frontier, R, vxm_desc) );
        #else
        CHECK( GrB_vxm(frontier, parent_list, NO_ACCUM, GxB_ANY_FIRSTJ_INT32, frontier, R, vxm_desc) );
        #endif

        //accum needed because otherwise the operation will act as a replace.
        //As with all graphblas operations, if only one operand is present it will be passed through
        //This could have been any binary operator. see top of p. 171 in the API
        CHECK( GrB_Vector_apply(parent_list, NO_MASK, GrB_PLUS_INT32, GrB_IDENTITY_INT32, frontier, DEFAULT_DESC) );

        CHECK( GrB_Vector_nvals(&frontier_nvals, frontier) );
    }

    if ((GrB_Vector_extractElement(&sink_parent, parent_list, sink) == GrB_NO_VALUE))
    {
        found_path = false;
    } else {
        found_path = true;
        // Extract path from source to sink from parent list (reverse traverse)
        // build a mask
        CHECK( GrB_Matrix_clear(M) );
        GrB_Index curr_vertex = sink;
        GrB_Index parent;
        while (curr_vertex != source)
        {
            CHECK( GrB_Vector_extractElement(&parent, parent_list, curr_vertex) );
            CHECK( GrB_Matrix_setElement(M, true, parent, curr_vertex) );
            curr_vertex = parent;
        }
    }

    GrB_free(&parent_list);
    GrB_free(&frontier);
    GrB_free(&vxm_desc);

    return found_path;
}


int main (int argc, char **argv)
{
    printf("Usage: ./edmund-karp filename.mtx (<runs>) (<s>) (<t>)\n");
    if(argc < 2) {
        return 1;
    }
    //Blocking better for debuggbing
    #if defined(DEBUG) || defined(PROFILE) || defined(DETERMINISTIC)
    GrB_init ( GrB_BLOCKING ) ;
    #else
    GrB_init ( GrB_NONBLOCKING );
    #endif

    size_t threads = omp_get_max_threads();
    printf("max threads:%ld\n", threads);

    GrB_Index n = 0;            //No. of vertices
    GrB_Index edges = 0;        //No. of edges
    GrB_Index* row_indeces;
    GrB_Index* col_indeces;
    double* values;
    GrB_Matrix A = NULL;    //Input graph
    GrB_Matrix M = NULL;    //Mask for path
    GrB_Matrix P = NULL;    //Weighted path
    GrB_Matrix R = NULL;    //Residual network
    GrB_Descriptor transpose_a;
    GrB_Descriptor replace;
    GrB_Vector index_ramp;  //v[i] = i. For extracting parent in BFS
    GrB_Index source = -1;
    GrB_Index sink = -1;
    struct timeval start, stop;
    double build_matrix_time = 0;
    double total_run_time = 0;
    double bfs_time = 0;
    double augment_time = 0;
    double flow_calc_time = 0;
    double min_cut_time = 0;

    gettimeofday ( &start, NULL );

    readMtx(argv[1], &n, &edges, &row_indeces, &col_indeces, &values);
    CHECK( GrB_Matrix_new (&A, GrB_FP64, n, n) );
    CHECK( GrB_Matrix_build (A, row_indeces, col_indeces, values, edges, GrB_PLUS_FP64) );

    gettimeofday(&stop, NULL);
    build_matrix_time = (WALLTIME(stop)-WALLTIME(start));
    printf("Read&build matrix time: %lf\n", build_matrix_time);

    GrB_Descriptor_new(&transpose_a);
    GrB_Descriptor_set(transpose_a, GrB_INP0, GrB_TRAN);

    GrB_Descriptor_new(&replace);
    GrB_Descriptor_set(replace, GrB_OUTP, GrB_REPLACE);

    free(row_indeces);
    free(col_indeces);
    free(values);

#ifdef DEBUG
    GxB_print(A, GxB_SHORT);
#endif
    GrB_Matrix_new(&M, GrB_BOOL, n, n);
    GrB_Matrix_new(&P, GrB_FP64, n, n);

    if(argc < 5) {
        s_t_bfs(A, &source, &sink);
        //return 0;
    } else {
        source = atoi(argv[3]);
        sink = atoi(argv[4]); //TODO set sink to something sensible
        printf("SOURCE:%ld, SINK:%ld\n", source, sink);
        if (sink > n) {
            printf("Out of bounds sink y'all ERROR to the MAX!\n");
            return 1;
        }
    }

    if(argc < 3) {
        runs = 3;
    } else {
        runs = atoi(argv[2]);
    }
    printf("Runs:%d\n", runs);

    for (; run < runs; run++)
    {


        total_run_time = 0;
        bfs_time = 0;
        augment_time = 0;
        flow_calc_time = 0;
        min_cut_time = 0;
        gettimeofday(&start, NULL);

        printf("\n-------Run %d-------\n", run+1);

        CHECK( GrB_Matrix_dup(&R, A) );
        size_t count = 1;
        GrB_Index nvals;
        GrB_Matrix_nvals(&nvals, A);

        gettimeofday(&stop, NULL);
        total_run_time += (WALLTIME(stop)-WALLTIME(start));
        gettimeofday(&start, NULL);

        printf("----------- Iteration: %ld -----------\n", count);
        while (get_augmenting_path(R, source, sink, M) && count++ < nvals)
        {
            gettimeofday(&stop, NULL);
            total_run_time += (WALLTIME(stop)-WALLTIME(start));
            bfs_time += (WALLTIME(stop)-WALLTIME(start));
            gettimeofday(&start, NULL);

            //Extract the capacity values coinciding with the path
            CHECK( GrB_eWiseMult(P, NO_MASK, NO_ACCUM, GrB_SECOND_FP64, M, R, DEFAULT_DESC) );
#ifdef DEBUG
            printf("\nPath:\n");
            GxB_print(P, GxB_SHORT);
#endif
            if(run == runs-1){
                GrB_Index path_len = 0;
                GrB_Matrix_nvals(&path_len, P);
                printf("\nPath length: %ld\n", path_len);
            }
            //Reduction to get the delta, the minimum capacity along the path
            //matrix to vector reduction also has a mask argument, not optional, but can be NULL
            double delta_f = 0;
            CHECK(GrB_reduce(&delta_f, NO_ACCUM, GrB_MIN_MONOID_FP64, P, DEFAULT_DESC));

            //Assign the negated delta first using the path mask to get the negative capacity changes...
            CHECK( GrB_assign(P, M, NO_ACCUM, -delta_f, GrB_ALL, n, GrB_ALL, n, DEFAULT_DESC) );
            //...then apply the transposed negated result matrix to itself to get the positive capacity changes

            CHECK( GrB_Matrix_apply(P, NO_MASK, GrB_PLUS_FP64, GrB_AINV_FP64, P, transpose_a) );

            //Add values to residual, reducing/increasing capacities as defined
            GrB_eWiseAdd(R, NO_MASK, NO_ACCUM, GxB_PLUS_FP64_MONOID, R, P, DEFAULT_DESC);

             //Remove zero-edges
            GrB_apply(R, R, NO_ACCUM, GrB_IDENTITY_FP64, R, replace);
#ifdef DEBUG
            GxB_print(R, GxB_SHORT);
#endif

            gettimeofday(&stop, NULL);
            total_run_time += (WALLTIME(stop)-WALLTIME(start));
            augment_time += (WALLTIME(stop)-WALLTIME(start));
            gettimeofday(&start, NULL);

            if(count < 20 || run == runs-1)
                printf("----------- Iteration: %ld -----------\n", count);
            else if(count == 20)
                printf("...\n\n");
        }

        gettimeofday(&stop, NULL);
        total_run_time += (WALLTIME(stop)-WALLTIME(start));
        gettimeofday(&start, NULL);

        double total_flow = 0;
        for (GrB_Index i = 0; i < n; i++) {
            double capacity = 0;
            double residual = 0;

            if(GrB_Matrix_extractElement(&capacity, A, source, i) == GrB_NO_VALUE)
                capacity = 0;
            if(GrB_Matrix_extractElement(&residual, R, source, i) == GrB_NO_VALUE)
                residual = 0;

            if(capacity != 0){
                total_flow += capacity-residual;
                //if(run == runs-1)
                    //printf("Found edge (%ld,%ld) w cap: %lf and residual: %lf. Adding %lf to total flow\n", source, i, capacity, residual, capacity-residual);
            }
            //printf("Flow now at %d\n", total_flow);
        }

        //End after determining flow, but before correctnes
        gettimeofday(&stop, NULL);
        total_run_time += (WALLTIME(stop)-WALLTIME(start));
        flow_calc_time += (WALLTIME(stop)-WALLTIME(start));
        gettimeofday(&start, NULL);


        printf("s-t: %ld-%ld. Max flow: %lf\n", source, sink, total_flow);
        printf("Total run time: %lf\n", total_run_time);
        printf("Time spent finding augmenting paths (BFS): %lf\n", bfs_time);
        printf("Time computing and augmenting flow: %lf\n", augment_time);
        printf("Time spent summing total flow: %lf\n", flow_calc_time);

        if (run == runs-1) {
            printf("\nCorrectness and min-cut:\n");

            double total_source_flow = total_flow;
            printf("Total flow out of source: %lf\n", total_source_flow);

            double total_sink_flow = 0;
            for (GrB_Index i = 0; i < n; i++) {
                double capacity = 0;
                double residual = 0;

                if(GrB_Matrix_extractElement(&capacity, A, i, sink) == GrB_NO_VALUE)
                    capacity = 0;
                if(GrB_Matrix_extractElement(&residual, R, i, sink) == GrB_NO_VALUE)
                    residual = 0;
                total_sink_flow += capacity-residual;
            }
            printf("Total flow into sink: %lf\n", total_sink_flow);

            GrB_Matrix min_cut;
            get_mincut(&min_cut, A, R, source);

            double min_cut_value = -1;
            CHECK( GrB_reduce(&min_cut_value, NO_ACCUM, GrB_PLUS_MONOID_FP64, min_cut, DEFAULT_DESC) );

            printf("Min-cut sum: %lf\n", min_cut_value);
#ifdef DEBUG
            GxB_print(min_cut, GxB_SHORT);
#endif
            min_cut_time += (WALLTIME(start)-WALLTIME(stop));

            printf("Min-cut time: %lf\n", min_cut_time);
            printf("Total number of iterations: %ld\n", count);

            GrB_free(&min_cut);
        }
    }
    GrB_finalize ( ) ;
}
