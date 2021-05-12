
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <GraphBLAS.h>

void readMtx(char* filename, GrB_Index* n, GrB_Index* edges, GrB_Index** I, GrB_Index** J, void** V)
//void readMtx(char* filename, GrB_Index* n, GrB_Index* edges, GrB_Index** I, GrB_Index** J, double** V)
{
    FILE* fp = fp = fopen(filename, "r");
    size_t line_size = 0;
    char* line = NULL;

    getline(&line, &line_size, fp);
    char delim[] = " \t\r\n\v\f";
    char* word = strtok(line, delim);

    bool is_symmetric, is_integer, is_real, is_pattern;
    while(word != NULL){
        if(strcmp(word, "symmetric") == 0)
            is_symmetric = true;
        else if(strcmp(word, "integer") == 0)
            is_integer = true;
	else if(strcmp(word, "real") == 0)
            is_real = true;
        else if(strcmp(word, "pattern") == 0)
            is_pattern = true;
        word = strtok(NULL, delim);
    }

    while(line[0] == '%')
        getline(&line, &line_size, fp);

    //printf("Line len: %ld\nContents: %s", line_size, line);
    int m=0, e=0;
    sscanf(line, "%ld %d %d", n, &m, &e);
    //printf("Exracted V. n:%ld m:%d e:%d\n", *n, m, e);

    *I = malloc(sizeof(GrB_Index)*e);
    *J = malloc(sizeof(GrB_Index)*e);

    *V = malloc(sizeof(double)*e); //
/* Todo: make generic
    if (is_pattern) {
        {};
    } else if (is_integer) {
	*V = malloc(sizeof(int)*e);
    } else if (is_real) {
        *V = malloc(sizeof(double)*e);
    } else {
        printf ("BIG MF ERROR HERE\n");
    }
*/

    GrB_Index i, j;
    int count = 0;
    while(getline(&line, &line_size, fp) != -1){
    //    printf("READING: %s", line);

        double v = 0;
        sscanf(line, "%ld %ld %lf", &i, &j, &v);
        if(v == 0)
            v = rand() % (int)1.7976931348623157E+308 + 1;
        else if(v < 0)
            v = -v;
        (*(double**)V)[count] = v;

        (*I)[count] = i-1;
        (*J)[count] = j-1;
/* TODO: Make generic
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
*/
        count++;
    }

    if(e==count) {
        printf("Everything is fine!\n");
        //for (int i = 0; i < e; i++)
        //    printf("V[%d]:%lf\n", i, (*(double**)V)[i]);
    }
    else
        printf("e:%d, count:%d. everything is not fine.....\n", e, count);

    *edges = e;
}
