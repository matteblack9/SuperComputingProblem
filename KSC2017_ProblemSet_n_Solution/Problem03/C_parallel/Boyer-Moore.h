#ifndef __BOYER_MOORE_H__
#define __BOYER_MOORE_H__

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#define ALPHABET_LEN 256
#define NOT_FOUND pattern_length
#define MAX(a, b) ((a < b) ? b : a)

void	make_BCS(int32_t *BCS, char *pattern, int64_t pattern_length);
void	make_GSS(int32_t *GSS, char *pattern, int64_t pattern_length);
int32_t	is_prefix(char *word, int64_t wordlen, int64_t pos);
int64_t	get_suffix_length(char *word, int64_t wordlen, int64_t pos);
int64_t	do_search(char *target, int64_t target_length, int64_t search_idx, int64_t search_length, char *pattern, int64_t pattern_length, int32_t *BCS, int32_t *GSS);

#endif