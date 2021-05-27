


all:

edmund-karp: edmund-karp.c readmtx.h
	gcc edmund-karp.c -lgraphblas -o edmund-karp

debug: edmund-karp.c readmtx.h
	gcc edmund-karp.c -D DEBUG -lgraphblas -Og -o edmund-karp


profile: edmund-karp.c readmtx.h
	gcc edmund-karp.c -lgraphblas -Og -o edmund-karp-prof

.PHONY: run
run: edmund-karp
	./edmund-karp
