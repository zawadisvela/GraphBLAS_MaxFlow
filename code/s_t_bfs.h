#include <GraphBLAS.h>
#include <string.h>

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


static double dampening_factor = 1.2;
void dampen (void *result, const void *element)
{
    (* ((double *) result)) = (* ((double *) element)) / dampening_factor;
}


void s_t_bfs(
    const GrB_Matrix A,
    GrB_Index *s,
    GrB_Index *t
)
{
    GrB_Index n;
    GrB_Descriptor vxm_desc = NULL;
    GrB_Vector weights = NULL;
    GrB_Vector depths = NULL;
    GrB_Vector frontier = NULL;
    GrB_Matrix A_t = NULL;
    GrB_UnaryOp dampen_op;

    GrB_Index num_reachable = 0;
    bool single_try = false;

    CHECK( GrB_Matrix_nrows(&n, A) );

    GrB_Descriptor_new(&vxm_desc);
    GrB_Descriptor_set(vxm_desc, GrB_MASK, GrB_COMP);
    GrB_Descriptor_set(vxm_desc, GrB_MASK, GrB_STRUCTURE);
    GrB_Descriptor_set(vxm_desc, GrB_OUTP, GrB_REPLACE);

    CHECK( GrB_Vector_new(&weights, GrB_FP64, n) );
    CHECK( GrB_Vector_new(&depths, GrB_INT32, n) );
    CHECK( GrB_Vector_new(&frontier, GrB_FP64, n) );

    GrB_Matrix_new(&A_t, GrB_BOOL, n, n);
    GrB_transpose(A_t, NO_MASK, NO_ACCUM, A, DEFAULT_DESC);

    CHECK( GrB_UnaryOp_new(&dampen_op, dampen, GrB_FP64, GrB_FP64) );

    if((long)(*s) < 0){
        *s = 0; //rand() % n; Start w 0, then move to random
    }else{
        single_try = true;
    }

    GrB_Index current_source = *s;
    GrB_Matrix current_graph = A;

    GrB_Index s_prev = 0;
    GrB_Index t_prev = 0;

    int iter = -1;
    int max_iter = 5;

    while (true) {
        //Force minimum 1 iteration
        s_prev = *s+1;
        t_prev = *t+1;
        for (iter = 0; iter < max_iter && (s_prev != *s || t_prev != *t); iter++) {
            printf("------------------------\n" );
            printf("Iteration:%d\t", iter);

            if (iter % 2 == 0) {
                printf("Forward search from %ld\n", *s);
                current_graph = A;
                current_source = *s;
            }else{
                printf("Backward search from %ld\n", *t);
                current_graph = A_t;
                current_source = *t;
            }
            s_prev = *s;
            t_prev = *t;
            printf("\n");

            GrB_Vector_clear(frontier);
            GrB_Vector_clear(depths);
            GrB_Vector_clear(weights);

            GrB_Vector_setElement(frontier, 1, current_source);
            GrB_Vector_setElement(weights, 1, current_source);
            GrB_Vector_setElement(depths, 0, current_source);

            bool successor = true;
            int depth = 0;
            GrB_Index frontier_nvals = -1;

            GrB_Vector_nvals(&frontier_nvals, frontier);
            //while (successor) {
            while (frontier_nvals > 0) {
                if(depth < 500)
                    printf("Depth: %d, frontier size: %ld \n", depth, frontier_nvals);
                else if (depth == 500)
                    printf("...\n");
                depth++;
                //GxB_print(weights, GxB_SHORT);
                //GrB_eWiseAdd(weights, NO_MASK, NO_ACCUM, GxB_PLUS_FP64_MONOID, weights, frontier, DEFAULT_DESC);

                //GxB_print(weights, GxB_SHORT);
                //GxB_print(depths, GxB_SHORT);
                //GrB_Vector frontier_2;
                //GrB_Vector_new(&frontier_2, GrB_FP64, n);

                GrB_vxm(frontier, weights, NO_ACCUM, GxB_PLUS_FIRST_FP64, frontier, current_graph, vxm_desc);

                CHECK( GrB_Vector_apply(frontier, NO_MASK, NO_ACCUM, dampen_op, frontier, DEFAULT_DESC) );

                CHECK( GrB_Vector_apply(weights, NO_MASK, GrB_PLUS_FP64, GrB_IDENTITY_FP64, frontier, DEFAULT_DESC) );

                GrB_Vector_assign_INT32(depths, frontier, NO_ACCUM, depth, GrB_ALL, n, DEFAULT_DESC);

                //GrB_reduce(&successor, NO_ACCUM, GrB_LOR_MONOID_BOOL, frontier, DEFAULT_DESC);

                GrB_Vector_nvals(&frontier_nvals, frontier);
            }
        #ifdef DEBUG
            GxB_print(weights, GxB_SHORT);
        #endif
            double max_weight = -1;
            CHECK( GrB_reduce(&max_weight, NO_ACCUM, GrB_MAX_MONOID_FP64, weights, DEFAULT_DESC) );
            printf("Max weight:%lf", max_weight);

            int max_depth = -1;
            GrB_Index max_index;
            double val = -1;
            for(int j = 0; j < n; j++){
                if(GrB_Vector_extractElement(&val, weights, j) != GrB_NO_VALUE && val == max_weight){
                    max_index = j;
                    GrB_Vector_extractElement(&max_depth, depths, j);
                    break;
                }
            }
            printf(" - depth: %d", max_depth);
            printf(" - dampening_factor: %lf\n", dampening_factor);

            if(current_source == *s) {
                *t = max_index;
            }else{
                *s = max_index;
            }

            printf("\nCandidate s-t: %ld-%ld\n", *s, *t);
            GrB_Vector_nvals(&num_reachable, weights);
            printf("Number reachable vertices:%ld/%ld\n\n", num_reachable, n);
        }

        printf("\nFinished in %d/%d iterations\n\n", iter, max_iter);

        printf("Finished w these states: iter = %d, max_iter = %d\n s_prev = %ld, *s = %ld\n t_prev = %ld, *t = %ld\n",
            iter, max_iter, s_prev, *s, t_prev, *t);

        char answer[255];
        int successful_scans = -1;
        printf("Enter new search: <dampening> <source>/\"random\"\nEnter anything else to abort\n" );
        successful_scans = scanf("%lf %s", &dampening_factor, answer);
        printf("Answers: %lf %s\n", dampening_factor, answer);
        if(successful_scans < 2){
            break;
        }else if(strcmp(answer, "random") == 0){
            *s = rand() % n;
            printf("Random source: %ld!!!!\n", *s);
        }else{
            GrB_Index new_source = -1;
            successful_scans = sscanf(answer, "%ld\n", &new_source);
            if(successful_scans == 1){
                *s = new_source;
            }else{
                printf("Failed to read new source, using previous\n");
            }
        }
        *t = *s;
        /*
        printf("New search? (y/n)\n");
        scanf("\n%c", &answer);
        printf("The answer you gave: %c\n", answer);
        if(answer == 'y' || answer == 'Y') {
            printf("Same source? y/n/R(andom)\n");
            scanf(" %c", &answer);
            printf("The answer you gave: %c\n", answer);
            if(answer == 'r' || answer == 'R') {
                *s = rand() % n;
                printf("New source: %ld\n", *s);
            }else if(answer == 'n' || answer == 'N'){
                printf("Enter source y/n/R(andom)\n");
                scanf(" %d", &new_source);
                *s = new_source;
            }
            printf("Enter new dampening factor. Current dampening factor = %lf\n", dampening_factor);
            scanf("%lf", &dampening_factor);
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
