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


static float dampening_factor = 1.2;
void dampen (void *result, const void *element)
{
    (* ((float *) result)) = (* ((float *) element)) / dampening_factor;
}


void s_t_bfs(
    const GrB_Matrix A,
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

    GrB_Vector weights = NULL;
    GrB_Vector depths = NULL;
    GrB_Vector frontier = NULL;

    CHECK( GrB_Vector_new(&weights, GrB_FP32, n) );
    CHECK( GrB_Vector_new(&depths, GrB_INT32, n) );
    CHECK( GrB_Vector_new(&frontier, GrB_FP32, n) );

    int tries = 0;
    bool single_try = false;

    if((long)(*s) < 0){
        *s = 0; //rand() % n;
    }else{
        single_try = true;
    }

    while (true) {
        printf("Try number:%d\tSource:%ld\n", tries, *s);

        GrB_Vector_setElement(frontier, 1, *s);
        GrB_Vector_setElement(weights, 1, *s);

        GrB_UnaryOp dampen_op;
        CHECK( GrB_UnaryOp_new(&dampen_op, dampen, GrB_FP32, GrB_FP32) );


        bool successor = true;
        int depth = 0;
        GrB_Index frontier_nvals = -1;

        GrB_Vector_nvals(&frontier_nvals, frontier);
        printf("Depth: %d, frontier size: %ld \n", depth++, frontier_nvals);
        while (successor) {
    #ifdef DEBUG
            printf("current depth: %d \n", depth);
    #endif
            GrB_vxm(frontier, weights, NO_ACCUM, GxB_PLUS_FIRST_FP32, frontier, A, desc);
            GrB_reduce(&successor, NO_ACCUM, GrB_LOR_MONOID_BOOL, frontier, DEFAULT_DESC);

            //printf("\n>>>>>pre-apply:\n");
            //GxB_print(frontier, GxB_SHORT);
            CHECK( GrB_Vector_apply(frontier, NO_MASK, NO_ACCUM, dampen_op, frontier, DEFAULT_DESC) );
            //printf(">>>>>post-apply:\n");
            //GxB_print(frontier, GxB_SHORT);
            //printf(">>>>>>>>>>>>>>>\n");



            GrB_eWiseAdd(weights, NO_MASK, NO_ACCUM, GxB_PLUS_FP32_MONOID, weights, frontier, DEFAULT_DESC);

            depth++;
            GrB_Vector_assign_INT32(depths, frontier, NO_ACCUM, depth, GrB_ALL, n, DEFAULT_DESC);

            GrB_Vector_nvals(&frontier_nvals, frontier);
            printf("Depth: %d, frontier size: %ld \n", depth, frontier_nvals);
        }
    #ifdef DEBUG
        GxB_print(weights, GxB_SHORT);
    #endif
        float max_weight = -1;
        CHECK( GrB_reduce(&max_weight, NO_ACCUM, GrB_MAX_MONOID_FP32, weights, DEFAULT_DESC) );
        printf("Max weight:%f", max_weight);

        int t_depth = -1;
        float val = -1;
        for(int j = 0; j < n; j++){
            if(GrB_Vector_extractElement(&val, weights, j) != GrB_NO_VALUE && val == max_weight){
                *t = j;
                GrB_Vector_extractElement(&t_depth, depths, j);
                break;
            }
        }
        printf(" - depth: %d", t_depth);
        printf(" - dampening_factor: %f\n", dampening_factor);

        printf("\nCandidate s-t: %ld-%ld\n", *s, *t);
        GrB_Vector_nvals(&num_reachable, weights);
        printf("Number reachable vertices:%ld/%ld\n\n", num_reachable, n);

        char answer;
        printf("New search? (y/n)\n");
        scanf("\n%c", &answer);
        printf("The answer you gave: %c\n", answer);
        if(answer == 'y' || answer == 'Y') {
            printf("Same source? (y/n)\n");
            scanf(" %c", &answer);
            printf("The answer you gave: %c\n", answer);
            if(answer == 'n' || answer == 'N') {
                *s = rand() % n;
                printf("New source: %ld\n", *s);
            }
            printf("Enter new dampening factor. Current dampening factor = %f\n", dampening_factor);
            scanf("%f", &dampening_factor);
            GrB_Vector_clear(frontier);
            GrB_Vector_clear(depths);
            GrB_Vector_clear(weights);
        } else {
            break;
        }
/*
        if(num_reachable < n/3 && tries++ < n/2 && !single_try) {
            *s = rand() % n;
            GrB_Vector_clear(frontier);
            GrB_Vector_clear(weights);
        } else {
            break;
        }
*/
    }
    printf("----------------\n");
    if(num_reachable >= n/3){
        printf("\ns-t: %ld-%ld\n", *s, *t);
    }else{
        printf("\ns-t: %ld-%ld\n", *s, *t);
        printf("\nCould not find source reaching at least half of graph\n");
    }
    GrB_free(&frontier);
    GrB_free(&depths);
    GrB_free(&weights);
}
