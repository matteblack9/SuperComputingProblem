#include "Boyer-Moore.h"
#include <unistd.h>
#include <string.h>
#include <mpi.h>

#define MAX_STRING_LENGTH 256

int64_t get_filesize(char *target_path)
{
	int64_t buffsize;
	int32_t fd = open(target_path, O_RDONLY);
	if(fd == -1)
		return -1;
	struct stat fd_stat;
	fstat(fd, &fd_stat);
	buffsize = fd_stat.st_size;
	return buffsize;
}

char* read_targetfile(char *target, int64_t target_length, char *target_path)
{
	FILE *fp = fopen((char*)target_path, "r");
	*target = '\0';

	if (fp != NULL) {
		size_t newLen = fread(target, sizeof(char), target_length, fp);
		if (newLen == 0) {
			printf("\nError: Cannot read target file [ %s ]\n", target_path);
			MPI_Finalize();
			exit(1);
		} else {
			target[newLen] = '\0'; /* Just to be safe. */
		}
	}
	fclose(fp);
	return target;
}

int32_t main(int32_t argc, char *argv[])
{
	int32_t rankID = 0, nRanks = 1;
	char rankName[MAX_STRING_LENGTH];
	gethostname(rankName, MAX_STRING_LENGTH);
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rankID);
	MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

	char *target_path = argv[1], *target = NULL;
	int64_t target_length = 0;
	target_length = get_filesize(target_path);
	if(target_length < 0)
	{
		printf("\nError: Cannot read target file [ %s ]\n", target_path);
		MPI_Finalize();
		exit(-1);
	}
	if(rankID == 0)
	{
		printf("--------------------------------------------------\n");
		printf("- Read target: [ %s ]\n", target_path);
		target = (char*)malloc(sizeof(char)*target_length);
		double read_time = 0;
		read_time -= MPI_Wtime();
		read_targetfile(target, target_length, target_path);
		read_time += MPI_Wtime();
		printf("- Target length: %ld (read time: %lf secs)\n", target_length, read_time);
		printf("--------------------------------------------------\n");
	}

	char *pattern = argv[2];
	int64_t pattern_length = 0;
	if(pattern == NULL)
	{
		printf("\nError: Cannot read pattern [ %s ]\n", pattern);
		free(target);
		MPI_Finalize();
		exit(-1);
	}
	pattern_length = strlen(pattern);
	if(rankID == 0)
	{
		printf("- Pattern: [ %s ]\n", pattern);
		printf("- Pattern length: %ld\n", pattern_length);
		printf("--------------------------------------------------\n");
	}
	int32_t* BCS = (int32_t*)malloc(ALPHABET_LEN * sizeof(int32_t));
	int32_t* GSS = (int32_t*)malloc(pattern_length * sizeof(int32_t));;
	make_BCS(BCS, pattern, pattern_length);
	make_GSS(GSS, pattern, pattern_length);

	int64_t found_count = 0;
	double search_time = 0;
	if(rankID == 0)
	{
		search_time -= MPI_Wtime();
	}
// DO NOT EDIT UPPER CODE //
//==============================================================================================================//




	found_count = do_search(target, target_length, 0, target_length, pattern, pattern_length, BCS, GSS);
	if(found_count < 0)
	{
		free(target);
		free(BCS);
		free(GSS);
		MPI_Finalize();
		exit(-1);
	}

//==============================================================================================================//
// DO NOT EDIT LOWER CODE //
	if(rankID == 0)
	{
		search_time += MPI_Wtime();
		printf("- Found count: %ld\n", found_count);
		printf("--------------------------------------------------\n");
		printf("- Time: %lf secs\n", search_time);
		printf("--------------------------------------------------\n");
	}

	free(target);
	free(BCS);
	free(GSS);
	MPI_Finalize();

	return 0;
}