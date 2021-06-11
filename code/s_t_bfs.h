#include <GraphBLAS.h>

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

void s_t_bfs(
    GrB_Matrix A,
    GrB_Index *s,
    GrB_Index *t
)
{
    GrB_Index n;
    CHECK( GrB_Matrix_nrows(&n, A) );

    GrB_Index num_reachable = 0;

    GrB_Descriptor desc = NULL;
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_COMP);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_STRUCTURE);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

    GrB_Vector level = NULL;
    GrB_Vector frontier = NULL;

    CHECK( GrB_Vector_new(&level, GrB_INT32, n) );
    CHECK( GrB_Vector_new(&frontier, GrB_BOOL, n) );

    int tries = 0;

    while(num_reachable < n/2 && tries++ < n/2) {
        *s = rand() % n;

        printf("Try:%d\tSource:%ld\n", tries, *s);

        GrB_Vector_setElement(frontier, true, *s);

        bool successor = true;
        int depth = 0;
        while (successor) {
    #ifdef DEBUG
            printf("current depth: %d \n", depth);
    #endif
            GrB_Vector_assign_INT32     // C<Mask>(I,J) = accum (C(I,J),x)
            (
                level,                   // input/output matrix for results
                frontier,          // optional mask for C, unused if NULL
                NO_ACCUM,       // optional accum for Z=accum(C(I,J),x)
                depth,                   // scalar to assign to C(I,J)
                GrB_ALL,             // row indices
                n,                  // number of indices
                DEFAULT_DESC       // descriptor for C and Mask
            );
            GrB_vxm(frontier, level, NO_ACCUM, GxB_LOR_LAND_BOOL, frontier, A, desc);
            GrB_reduce(&successor, NO_ACCUM, GrB_LOR_MONOID_BOOL, frontier, DEFAULT_DESC);
            depth++;
        }
    #ifdef DEBUG
        GxB_print(level, GxB_SHORT);
    #endif
        GrB_Index deepest = -1;
        CHECK( GrB_reduce(&deepest, NO_ACCUM, GrB_MAX_MONOID_INT32, level, DEFAULT_DESC) );
        printf("Deepest:%ld\n", deepest);

        int val = -1;
        for(int j = 0; j < n; j++){
            if(GrB_Vector_extractElement(&val, level, j) != GrB_NO_VALUE && val == depth-1){
                *t = j;
                break;
            }
        }

        GrB_Vector_nvals(&num_reachable, level);
        printf("Number reachable vertices:%ld\n\n", num_reachable);

        GrB_Vector_clear(frontier);
        GrB_Vector_clear(level);
    }
    printf("----------------\n");
    if(num_reachable >= n/2){
        printf("\ns-t: %ld-%ld\n", *s, *t);
    }else{
        printf("\nCould not find source reaching at least half of graph\n");
    }
    GrB_free(&frontier);
    GrB_free(&level);
}
