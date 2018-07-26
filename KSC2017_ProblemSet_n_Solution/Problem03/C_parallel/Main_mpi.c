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

	int64_t mpi_found_count = 0;
	char* chunk = NULL;
	if(argv[3] == NULL)
	{
		printf("\nError: Check chunk size [ %s ]\n", argv[3]);
		free(target);
		free(BCS);
		free(GSS);
		MPI_Finalize();
		exit(-1);
	}

	if(rankID == 0) printf("\ttarget_length = %ld\n", target_length);
	int64_t nChunksPerRank = atoi(argv[3]);
	//각 rank에 몇개의 문자열 덩어리를 줄것인가 결정
	int64_t nTotalChunks = (nRanks-1) * nChunksPerRank; 
	//rank 0은 검사하지않으므로 nRanks - 1
	//문자열은 총 nTotalChunks개로 쪼개진다.
	if(rankID == 0) printf("\tnTotalChunks = %ld\n", nTotalChunks);
	
	int64_t overlap_length = (pattern_length - 1) * (nTotalChunks - 1);
	//쪼개진 덩어리 중 마지막 1개는 겹치는 부분이 없으므로 nTotalChunks - 1
	//문자열은 최악의 경우 pattern의 첫글자가 하나의 코어 
	//나머지 글자가 하나의 코어에 분배되는 경우이므로 pattern_length - 1

	if(rankID == 0) printf("\toverlap_length = %ld\n", overlap_length);
	int64_t quotient = (target_length + overlap_length) / nTotalChunks; 
	//각 코어당 최악의 경우를 방지하기 위해 덩어리마다 pattern_length - 1을 추가
	//즉 target_length + overlap_length가 되고 이를 정해진 nChunksPerRank씩
	//각 코어에 분배하기 위하여 nTotalChunks로 나누어 준다.

	if(rankID == 0) printf("\tquotient = %ld\n", quotient);
	int64_t remainder = (target_length + overlap_length) - (quotient * nTotalChunks);
	//나누는 경우에 나누어 떨어지지 않는 경우가 있으므로 나머지를 따로 처리해준다.
	
	if(rankID == 0) printf("\tremainder = %ld\n\n", remainder);

	int64_t chunkID = 0;
	int64_t* chunk_length = (int64_t*)malloc((nTotalChunks+1)*sizeof(int64_t)); 
	int64_t* chunk_start_idx = (int64_t*)malloc((nTotalChunks+1)*sizeof(int64_t)); 
	//remainder의 경우를 위해 nTotalChunks에 + 1 을 한다.
	
	int64_t i;
	for(i=0; i<nTotalChunks; i++)
		chunk_length[i] = quotient;
	for(i=0; i<remainder; i++)
		chunk_length[i] += 1;
	chunk_start_idx[0] = 0;
	for(i=1; i<nTotalChunks; i++)
		chunk_start_idx[i] = chunk_start_idx[i-1] + chunk_length[i-1] - (pattern_length-1); 
	//마지막에 - (pattern_length - 1) 을 해줌으로서 첫 번째 chunk를 제외하고
	//모든 chunk는 이전 chunk의 마지막 문자열의 -4번째 포인터를 chunk_start_idx로 가진다.
	//따라서 첫번째 chunk를 제외한 모든 chunk 이전 chunk의 마지막 4글자를 무조건 포함한다.

	chunk_start_idx[nTotalChunks] = 0;
	chunk_length[nTotalChunks] = 0;

	//chunk가 끝났다는 것을 표시하기 위해 nTotalChunk + 1번째 chunk의
	//start idx 와 length는 모두 0으로 지정한다.

	MPI_Request MPI_req[2];
	MPI_Status MPI_stat[2];
	int32_t MPI_tag = 0;
	int32_t request_rankID = -1;
	if(rankID == 0)
	{
		int64_t nFinishRanks = 0;
		while(nFinishRanks < nRanks-1)
		{
			MPI_Recv(&request_rankID, 1, MPI_INT32_T, MPI_ANY_SOURCE, MPI_tag, MPI_COMM_WORLD, &MPI_stat[0]);
			MPI_Isend(&target[chunk_start_idx[chunkID]], chunk_length[chunkID], MPI_CHAR, request_rankID, chunkID, MPI_COMM_WORLD, &MPI_req[1]);
			printf("\trequest_rankID = %d\n", request_rankID);
			printf("\tchunkID = %ld\n", chunkID);
			printf("\tchunk_start_idx[chunkID] = %ld\n", chunk_start_idx[chunkID]);
			printf("\ttarget[chunk_start_idx[chunkID] = %c\n\n", target[chunk_start_idx[chunkID]]);
			if(chunkID < nTotalChunks)
				chunkID++;
			else
				nFinishRanks++;
		}
	}
	else
	{
		chunk = (char *)malloc(chunk_length[0] * sizeof(char));
		int64_t chunk_found_count = 0;
		int64_t call_count = 0;
		while(chunkID < nTotalChunks)
		{
			MPI_Isend(&rankID, 1, MPI_INT32_T, 0, MPI_tag, MPI_COMM_WORLD, &MPI_req[0]);
			MPI_Recv(chunk, chunk_length[0], MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &MPI_stat[1]);
			printf("\trank = %d chunk = %s\n", rankID, chunk);
			chunkID = MPI_stat[1].MPI_TAG;
			if(chunkID < nTotalChunks)
			{
				chunk_found_count = do_search(chunk, target_length, 0, chunk_length[chunkID], pattern, pattern_length, BCS, GSS);
				if(found_count < 0)
				{
					free(chunk);
					free(BCS);
					free(GSS);
					free(chunk_length);
					free(chunk_start_idx);
					MPI_Finalize();
					exit(-1);
				}
				mpi_found_count += chunk_found_count;
				call_count++;
			}
		}
		printf("- [%02d: %s] call_count: %ld\n", rankID, rankName, call_count);
	}
	MPI_Reduce(&mpi_found_count, &found_count, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
	free(chunk);
	free(chunk_length);
	free(chunk_start_idx);

//==============================================================================================================//
// DO NOT EDIT LOWER CODE //
	if(rankID == 0)
	{
		search_time += MPI_Wtime();
		printf("- Found_count: %ld\n", found_count);
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
