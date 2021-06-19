
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GraphBLAS.h>

void readMtx(char* filename, GrB_Index* n, GrB_Index* edges, GrB_Index** I, GrB_Index** J, double** V)
{
    FILE* fp = fp = fopen(filename, "r");
    size_t line_size = 0;
    char* line = NULL;

    getline(&line, &line_size, fp);
    char delim[] = " \t\r\n\v\f";
    char* word = strtok(line, delim);

    bool is_symmetric=false, is_integer=false, is_real=false, is_pattern=false;
    while(word != NULL){
        if(strcmp(word, "symmetric") == 0) {
            is_symmetric = true;
        } else if(strcmp(word, "integer") == 0) {
            is_integer = true;
	    } else if(strcmp(word, "real") == 0) {
            is_real = true;
        } else if(strcmp(word, "pattern") == 0) {
            is_pattern = true;
        }
        word = strtok(NULL, delim);
    }

    printf("Graph properties:%s%s%s%s\n",
        is_symmetric?" symmetric,":"",
        is_pattern?" pattern,":"",
        is_real?" real,":"",
        is_integer?" integer,":"");

    while(line[0] == '%')
        getline(&line, &line_size, fp);

    size_t m=0, e=0;
    sscanf(line, "%ld %ld %ld", n, &m, &e);
    if (is_symmetric) {
        e *= 2;
    }
    printf("Graphsize: %ld %ld %ld\n", *n, m, e);
    *I = malloc(sizeof(GrB_Index)*e);
    *J = malloc(sizeof(GrB_Index)*e);
    *V = malloc(sizeof(double)*e);

    GrB_Index i, j;
    size_t count = 0;
    size_t total_updates = 1000;
    size_t update_count = (e/(total_updates*2))*2;//Need even number
    if(update_count == 0) {
        update_count = 2;
    }
    printf("Determined update frequency: %ld\n", update_count);
    while(getline(&line, &line_size, fp) != -1){
        #ifdef DEBUG
        if(count % update_count == 0){
            //printf("Reading:  %.4f%%\r", 100*count/(double)e);
        }
        #endif
        double v = 0;

        sscanf(line, "%ld %ld %lf", &i, &j, &v);

        if(v == 0) {
            v = rand() % (int)1000 + 40; //new edge weights
            //v = rand() % (int)1.7976931348623157E+308 + 1; //new edge weights
        }else if(v < 0){
            v = -v; //Need positive values for maxflow
        }

        (*V)[count] = v;
        (*I)[count] = i-1;
        (*J)[count] = j-1;

        if (is_symmetric) {
            count++;
            (*V)[count] = v;
            (*I)[count] = j-1;
            (*J)[count] = i-1;
        }

        count++;
    }
    if(e==count) {
        printf("Everything is fine!\n");
    } else {
        printf("\ne:%ld, count:%ld. \neverything is not fine.....\n", e, count);
    }

    *edges = e;
}
