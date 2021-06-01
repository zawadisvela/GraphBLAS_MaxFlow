

common_flags = -lgraphblas -fopenmp

all:

edmund-karp: edmund-karp.c readmtx.h
	gcc edmund-karp.c $(common_flags) -o edmund-karp

debug: edmund-karp.c readmtx.h
	gcc edmund-karp.c -D DEBUG $(common_flags) -Og -o edmund-karp-debug


profile: edmund-karp.c readmtx.h
	gcc edmund-karp.c -D PROFILE $(common_flags) -Og -o edmund-karp-prof

.PHONY: run
run: edmund-karp
	./edmund-karp
