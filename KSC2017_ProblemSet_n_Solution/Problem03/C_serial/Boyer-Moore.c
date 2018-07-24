#include "Boyer-Moore.h"

void make_BCS(int32_t *BCS, char *pattern, int64_t pattern_length)
{
	int64_t i;
	for(i=0; i < ALPHABET_LEN; i++)
	{
		BCS[i] = NOT_FOUND;
	}

	for(i=0; i < pattern_length-1; i++)
	{
		BCS[(int32_t)pattern[i]] = pattern_length-1 - i;
	}
}

void make_GSS(int32_t *GSS, char *pattern, int64_t pattern_length)
{
	int64_t last_prefix_index = pattern_length-1;
	int64_t i;
	for(i=pattern_length-1; i>=0; i--)
	{
		if(is_prefix(pattern, pattern_length, i+1))
		{
			last_prefix_index = i + 1;
		}
		GSS[i] = last_prefix_index + (pattern_length-1 - i);
	}
	for(i=0; i<pattern_length-1; i++)
	{
		int64_t suffix_length = get_suffix_length(pattern, pattern_length, i);
		if(pattern[i - suffix_length] != pattern[pattern_length-1 - suffix_length])
		{
			GSS[pattern_length-1 - suffix_length] = pattern_length-1 - i + suffix_length;
		}
	}
}

int32_t is_prefix(char *pattern, int64_t pattern_length, int64_t pattern_idx)
{
	int64_t suffix_length = pattern_length - pattern_idx;
	int64_t i;
	for(i=0; i<suffix_length; i++)
	{
		if(pattern[i] != pattern[pattern_idx+i])
		{
			return 0;
		}
	}
	return 1;
}

int64_t get_suffix_length(char *pattern, int64_t pattern_length, int64_t pattern_idx)
{
	int64_t i;
	for(i=0; (pattern[pattern_idx-i] == pattern[pattern_length-1 - i]) && (i < pattern_idx); i++);
	return i;
}

int64_t do_search(char *target, int64_t target_length, int64_t search_idx, int64_t search_length, char *pattern, int64_t pattern_length, int32_t *BCS, int32_t *GSS)
{
	int64_t target_cmp_idx, pattern_cmp_idx;
	int64_t found_count = 0;

	if(pattern_length == 0)
	{
		printf("\nError: pattern_length = 0\n");
		return -1;
	}

	target_cmp_idx = pattern_length-1;
	while(target_cmp_idx < search_length)
	{
		pattern_cmp_idx = pattern_length-1;
		while((pattern_cmp_idx >= 0) && (target[search_idx + target_cmp_idx] == pattern[pattern_cmp_idx]))
		{
			--target_cmp_idx;
			--pattern_cmp_idx;
		}
		if(pattern_cmp_idx < 0)
		{
			target_cmp_idx += GSS[0];
			found_count++;
		}
		else
		{
			target_cmp_idx += MAX(BCS[(int32_t)target[search_idx + target_cmp_idx]], GSS[pattern_cmp_idx]);
		}

		if(target_cmp_idx >= search_length)
		{
			break;
		}

		if(target[search_idx + target_cmp_idx] < 0)
		{
			printf("\nError: target index error - (%ld vs %ld + %ld)\n", target_length, search_idx, target_cmp_idx);
			return -1;
		}
	}
	return found_count;
}