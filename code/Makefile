all: mst

mst: mst-skeleton.c mst-solution.c
	mpicc -O3 mst-skeleton.c -lm -o mst

smpi: mst-skeleton.c mst-solution.c
	smpicc -O3 mst-skeleton.c -lm -o mst

clean:
	rm mst
