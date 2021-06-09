

common_flags = -lgraphblas -fopenmp

deps = readmtx.h s_t_bfs.h

all:

edmund-karp: edmund-karp.c $(deps)
	gcc edmund-karp.c $(common_flags) -o edmund-karp

debug: edmund-karp.c $(deps)
	gcc edmund-karp.c -D DEBUG $(common_flags) -Og -o edmund-karp-debug


profile: edmund-karp.c $(deps)
	gcc edmund-karp.c -D PROFILE $(common_flags) -Og -o edmund-karp-prof

apply: edmund-karp.c $(deps)
	gcc edmund-karp.c -D DEBUG -D APPLY $(common_flags) -Og -o edmund-karp-apply


s_t_search: s_t_search.c $(deps)
	gcc s_t_search.c $(common_flags) -o s_t_search
