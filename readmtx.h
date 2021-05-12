
#include <stdio.h>

void readMtx(char* filename, void* values, void* I, void* J)
{
    FILE *fp = fp = fopen(filename, "r");
    char buff[255];
    fscanf(fp, "%s", buff);
    printf("%s", buff);
}
