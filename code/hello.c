
#include "readmtx.h"


int main (int argc, char **argv)
{
    unsigned long* i;
    unsigned long* j;
    void* val;
    readMtx(argv[1], &i, &j, &val);
}
