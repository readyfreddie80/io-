#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>


#define MASTER 0

void io(int rank, int size, int l, int a, int b, int N) {

    int seed;
    int *result  = calloc(l * l * a * b, sizeof(int));

    if (rank == MASTER){
            srand(time(NULL));
            int *seeds = (int*)malloc(size * sizeof(int));
            for (int i = 0; i < size; i++)
                seeds[i] = (int)rand();
            MPI_Scatter(seeds, 1, MPI_INT, &seed, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
            free(seeds);
    }
    else {
            MPI_Scatter(NULL, 1, MPI_INT, &seed, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    }

    srand(seed);
    int ab = a * b;
    for (int i = 0; i < N; i++) {
		int x = rand() % l;
		int y = rand() % l;
		int r = rand() % ab;
		result[(y * l + x) * ab + r] += 1;
	}

	MPI_File data;
	MPI_Datatype contiguous, view;

	MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &data);
	MPI_File_set_size(data, 0);

	MPI_Type_contiguous(l * ab, MPI_INT, &contiguous);
	MPI_Type_create_resized(contiguous, 0, a * l * ab * sizeof(int), &view);
	MPI_Type_commit(&view);
	MPI_File_set_view(data, ((rank / a) * a * l + rank % a) * l * ab * sizeof(int), MPI_INT, view, "native", MPI_INFO_NULL);
	MPI_File_write_all(data, result, l * l * ab, MPI_INT, MPI_STATUS_IGNORE);

	MPI_File_close(&data);

	free(result);

}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int l = atoi(argv[1]);
    int a = atoi(argv[2]);
    int b = atoi(argv[3]);
    int N = atoi(argv[4]);


    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start, end;

    if(rank == MASTER) {
        start =  MPI_Wtime();
    }
    io(rank, size, l, a, b, N);

    if(rank == MASTER) {
        end = MPI_Wtime();
        double delta = end - start;
        FILE *f = fopen("stats.txt", "w");
        fprintf(f, "%d %d %d %d %fs", l, a, b, N, delta);

        fclose(f);
    }

    MPI_Finalize();
    return 0;
}
