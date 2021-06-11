///Index ramp based bfs stuff



CHECK( GrB_Vector_new(&index_ramp, GrB_INT32, n) );
for (int i = 0; i < n; i++)
        CHECK( GrB_Vector_setElement(index_ramp, i, i) );
#ifdef DEBUG
    printf("Index ramp:\n");
    GxB_print(index_ramp, GxB_SHORT);
#endif


// convert all stored values to their column index
CHECK( GrB_eWiseMult(wavefront,
    NO_MASK, NO_ACCUM, //mask, accumulate op
    GrB_FIRST_INT32,
    index_ramp, wavefront,
    DEFAULT_DESC) ); //descriptor


CHECK( GrB_vxm(wavefront, //c
    parent_list, //mask, ignores preciously found vertices and automatically removes self-loops
    NO_ACCUM, //no accumulate/default
    #if defined(DEBUG) || defined(PROFILE)
    GxB_MIN_FIRST_FP64, //semiring. First to extract parent. Why min? Could be ANY?
    #else
    GxB_ANY_FIRST_FP64, //Should yield same result, and leaves more implementation freedom
    #endif
    wavefront,
    R, // M
    desc) ); //descriptor



//Old garbage min-cut extraction

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

#if defined(PROFILE) || defined(DEBUG)
        GrB_Vector reachable = NULL;
        GrB_Vector_new(&reachable, GrB_BOOL, n);
        mincut_bfs(&reachable, R, source);

        printf("Min-cut edges:\n");
        double val;
        double min_cut_value = 0;
        for (GrB_Index i = 0; i < n; i++) {
            if(GrB_Vector_extractElement(&val, reachable, i) != GrB_NO_VALUE) {
                for (GrB_Index j = 0; j < n; j++) {
                    if(GrB_Vector_extractElement(&val, reachable, j) == GrB_NO_VALUE
                        && GrB_Matrix_extractElement(&val, A, i, j) != GrB_NO_VALUE)
                    {
                        printf("(%ld,%ld)\t%lf\n", i, j, val);
                        min_cut_value += val;
                    }
                }
            }
        }
        printf("Min cut: %lf\n", min_cut_value);

#endif



//Todo: make readmtx generic

void readMtx(char* filename, GrB_Index* n, GrB_Index* edges, GrB_Index** I, GrB_Index** J, void** V)

    if (is_pattern) {
        {};
    } else if (is_integer) {
	*V = malloc(sizeof(int)*e);
    } else if (is_real) {
        *V = malloc(sizeof(double)*e);
    } else {
        printf ("BIG MF ERROR HERE\n");
    }

    //Continued TODO: Make generic
        	if (is_pattern) {
                {};
            } else if (is_integer) {
                int v;
                sscanf(line, "%d %d %d", &i, &j, &v);
                (*(int**)V)[count] = v;
            } else if (is_real) {
    	    double v;
                sscanf(line, "%d %d %lf", &i, &j, &v);
                printf("-->extracted %d %d %lf\n", i, j, v);
                (*(double**)V)[count] = v;
            }


//Tested whether a transposed mask could be used to replace apply with assign, but it did not work
GrB_Matrix _P = NULL;
GrB_Matrix_dup(&_P, P);

GxB_print(P, GxB_SHORT);

GxB_print(_P, GxB_SHORT);

GrB_Descriptor transpose_mask;
CHECK( GrB_Descriptor_new(&transpose_mask) );
printf(">>>Setting transpose mask\n");
CHECK( GrB_Descriptor_set(transpose_mask, GrB_MASK, GrB_TRAN));
printf(">>>Transpose mask has been set\n");

CHECK( GrB_assign(_P, M, NO_ACCUM, delta_f, GrB_ALL, n, GrB_ALL, n, transpose_mask) );
GxB_print(_P, GxB_SHORT);
return false;

//Code for apply-based setting of delta-values in the path
#ifdef APPLY
            //Construct path and inverse, negative path
            printf(">>> We doing apply now");
            GrB_UnaryOp apply_delta_op;
            CHECK( GrB_UnaryOp_new(&apply_delta_op, apply_delta, GrB_FP64, GrB_BOOL) );
            CHECK( GrB_Matrix_apply(P, NO_MASK, NO_ACCUM, apply_delta_op, M, transpose_a) );
            CHECK( GrB_Matrix_apply(P, NO_MASK, GrB_PLUS_FP64, GrB_AINV_FP64, P, transpose_a) );
#else

double delta_f_global = 0 ;
void apply_delta (void *result, const void *element)
{
    (* ((double *) result)) = delta_f_global ;
}


//Matrix-based all-to-all bfs search for s & t. too time-consuming
//Perform multiple simultaneous BFS searches to determine s-t
void set_s_t_bfs(
    const GrB_Matrix A,
    GrB_Index *s,
    GrB_Index *t
)
{
    GrB_Index n;
    GrB_Matrix levels = NULL;
    GrB_Matrix frontiers = NULL;
    GrB_Descriptor desc = NULL;

    CHECK( GrB_Matrix_nrows(&n, A) );

    CHECK( GrB_Matrix_new(&levels, GrB_INT32, n, n) );
    CHECK( GrB_Matrix_new(&frontiers, GrB_INT32, n, n) );

if (false) { //Can enable all-to-all search by using GrB_extract afterwards
    for(int i = 0; i < n; i++){
        GrB_Matrix_setElement(frontiers, true, i, i);
    }
} else {
    for(int i = 0; i < 1; i++){
        GrB_Matrix_setElement(frontiers, true, i, i);
    }
}
    GrB_Descriptor_new(&desc);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_COMP);
    GrB_Descriptor_set(desc, GrB_MASK, GrB_STRUCTURE);
    GrB_Descriptor_set(desc, GrB_OUTP, GrB_REPLACE);

    bool successor = true;
    int level = 0;
    while (successor) {
        GrB_Matrix_assign_INT32     // C<Mask>(I,J) = accum (C(I,J),x)
        (
            levels,                   // input/output matrix for results
            frontiers,          // optional mask for C, unused if NULL
            NO_ACCUM,       // optional accum for Z=accum(C(I,J),x)
            level,                   // scalar to assign to C(I,J)
            GrB_ALL,             // row indices
            n,                   // number of row indices
            GrB_ALL,             // column indices
            n,                   // number of column indices
            DEFAULT_DESC       // descriptor for C and Mask
        );
        GrB_mxm(frontiers, levels, NO_ACCUM, GxB_LOR_LAND_BOOL, frontiers, A, desc);
        GrB_reduce(&successor, NO_ACCUM, GrB_LOR_MONOID_BOOL, frontiers, DEFAULT_DESC);
        level++;
        printf("current depth: %d \n", level);
    }
#ifdef DEBUG
    GxB_print(levels, GxB_SHORT);
#endif
    GrB_Vector deepest;
    CHECK( GrB_Vector_new(&deepest, GrB_INT32, n) );
    CHECK( GrB_reduce(deepest, NO_MASK, NO_ACCUM, GrB_MAX_MONOID_INT32, levels, DEFAULT_DESC) );
    printf("deepest:\n", deepest);
    GxB_print(deepest, GxB_SHORT);

    int val = -1;
    for(int i = 0; i < n; i++){
        if(GrB_Vector_extractElement(&val, deepest, i) != GrB_NO_VALUE && val == level-1){
            *s = i;
            for(int j = 0; j < n; j++){
                if(GrB_Matrix_extractElement(&val, levels, i, j) != GrB_NO_VALUE && val == level-1){
                    *t = j;
                    break;
                }
            }
            break;
        }
    }
    printf("s:%ld, t:%ld\n", *s, *t);

    GrB_Index num_reachable = -1;
    GrB_Vector_nvals(&num_reachable, levels);
    printf("Number reachable vertices:%ld\n", num_reachable);
}

//ewise addition for applying frontier to parent list in bfs, less efficient than apply
CHECK( GrB_eWiseAdd(parent_list, //w
    NO_MASK, //mask
    //GrB_PLUS_INT32, //accum
    NO_ACCUM,
    GrB_PLUS_INT32, //unary op. w<M> = accum (w, op (u)). Nor sure this identity is the same as for indexes
    parent_list,
    wavefront, //u
    DEFAULT_DESC) ); // descriptor

//Used for locating all filled edges, not necessary for finding min-cut

GrB_Descriptor mask_complement;
GrB_Descriptor_new(&mask_complement);
GrB_Descriptor_set(mask_complement, GrB_MASK, GrB_COMP);

GrB_Matrix min_cut;
GrB_Matrix_new(&min_cut, GrB_FP64, n, n);
GrB_apply(min_cut, R, NO_ACCUM, GrB_IDENTITY_FP64, A, mask_complement);



//Used to check the incoming flow into source, unecessary

if(GrB_Matrix_extractElement(&capacity, A, i, source) == GrB_NO_VALUE)
    capacity = 0;
if(GrB_Matrix_extractElement(&residual, R, i, source) == GrB_NO_VALUE)
    residual = 0;
else
    CHECK( GrB_Matrix_extractElement(&residual, R, i, source) );

if (i==1 && source == 0) {
    CHECK( GrB_Matrix_extractElement(&residual, R, i, source) );
    CHECK( GrB_Matrix_extractElement(&val, R, 1, 0) );
    printf("Double checking edge 1,0. Residual:%lf, val:%lf\n", residual, val);

printf("capacity along (%ld, %ld): %lf\n", i, source, capacity);
printf("residual along (%ld, %ld): %lf\n", i, source, residual);
printf("\n");

//crude remove self-loops in Matrix, unneeded
printf("Gonna build myself a matrix out of row_indeces:0x%lx, col_indeces:0x%lx, values:0x%lx, n:%ld\n", (uint64_t)row_indeces, (uint64_t)col_indeces, (uint64_t)values, n);
for (int i = 0; i < edges; i++) {
    //printf("%ld - %ld = %ld\n",row_indeces[i],col_indeces[i],row_indeces[i]-col_indeces[i]);
    /*if (row_indeces[i] == col_indeces[i]) {
        printf("boom?\n");
        //col_indeces[i] += 1;
        //col_indeces[i] = col_indeces[i] % n;
    }*/
    //values[i] = 1;
    printf("%ld,%ld - %lf\n",row_indeces[i],col_indeces[i],values[i]);
}

//Old hard-coded matrix
    GrB_Index n = 6;
    GrB_Index I[6] = {0,0,1,1,3,2} ;
    GrB_Index J[6] = {1,2,3,2,2,5} ;
    GrB_Matrix A = NULL;
    GrB_Index s = 0;

#if DOUBLE_PREC
    double X[6] = {4,2,3,1,5,7};               // # of nodes in the graph
    CHECK( GrB_Matrix_new (&A, GrB_FP64, n, n) );
    CHECK( GrB_Matrix_build (A, I, J, X, n, GrB_PLUS_FP64) );
#endif



//Old BFS-stuff

CHECK (GrB_Vector_new (&v, GrB_INT32, n)) ;

GrB_Vector_new (&q, GrB_BOOL, n) ;     // Vector<bool> q(n) = false
GrB_Vector_setElement_BOOL (q, true, s) ;   // q[s] = true, false elsewhere

GrB_Descriptor_new (&desc) ;
GrB_Descriptor_set (desc, GrB_MASK, GrB_COMP) ;     // invert the mask
GrB_Descriptor_set (desc, GrB_OUTP, GrB_REPLACE) ;  // clear q first

GrB_Semiring_new(&any_and_semiring, GxB_ANY_BOOL_MONOID, GxB_LAND_BOOL);

//These two lines just force a fully instantiated dense vector v. why?
//CHECK (GrB_Vector_assign_INT32 (v, NULL, NULL, 0, GrB_ALL, n, NULL) );// make dense. Zawadi: Why??
//GrB_Vector_nvals (&n, v) ;             // finish pending work on v


bool successor = true ; // true when some successor found
for (int32_t level = 1 ; successor && level <= n ; level++)
{
    // v<q> = level, using vector assign with q as the mask
    GrB_Vector_assign_INT32 (v, q, NULL, level, GrB_ALL, n, NULL) ;
    printf("\nROUND %d:\n\n", level);

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

GrB_Vector x;
GrB_Vector_new(&x, GrB_INT32, n);

CHECK( GrB_eWiseMult(x, GxB_DEFAULT, GxB_DEFAULT, GrB_TIMES_INT32, v, v, GxB_DEFAULT) );

GxB_print(x, GxB_COMPLETE);

CHECK( GrB_eWiseMult(x, GxB_DEFAULT, GxB_DEFAULT, GrB_FIRST_INT32, x, v, GxB_DEFAULT) );

GxB_print(x, GxB_COMPLETE);


//Failed attempt at using ewise add to add the negative edge changes
GrB_Descriptor transpose_b;
GrB_Descriptor_new(&transpose_b);
GrB_Descriptor_set(transpose_b, GrB_INP1, GrB_TRAN);
GrB_Matrix G;
GrB_Matrix_new(&G, GrB_INT8, n, n);
GxB_print(G, GxB_SHORT);

CHECK (GrB_eWiseAdd(G, NO_MASK, NO_ACCUM, GrB_PLUS_INT8, M, M, transpose_b));
GxB_print(G, GxB_SHORT);
CHECK (GrB_eWiseAdd(G, NO_MASK, NO_ACCUM, GrB_MINUS_INT8, M, M, transpose_b));
GxB_print(G, GxB_SHORT);

CHECK (GrB_eWiseAdd(G, NO_MASK, NO_ACCUM, GrB_PLUS_INT8, M, M, GxB_DEFAULT));
GxB_print(G, GxB_SHORT);
CHECK (GrB_eWiseAdd(G, NO_MASK, NO_ACCUM, GrB_MINUS_INT8, M, M, GxB_DEFAULT));
GxB_print(G, GxB_SHORT);

CHECK (GrB_eWiseAdd(G, NO_MASK, NO_ACCUM, GrB_MINUS_INT8, M, M, transpose_b));
CHECK (GrB_eWiseAdd(G, NO_MASK, NO_ACCUM, GrB_MINUS_INT8, G, M, transpose_b));
CHECK (GrB_eWiseAdd(G, NO_MASK, NO_ACCUM, GrB_MINUS_INT8, G, M, transpose_b));
GxB_print(G, GxB_SHORT);
