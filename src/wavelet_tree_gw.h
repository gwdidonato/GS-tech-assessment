#define WAVELET_TREE_H
#include "bitvector_bin_gw.c"
#include <stdio.h>
#include <stdbool.h>

/*
Each node has a bitvector struct and two children.
*/
typedef struct _data {
	bitvector_bin_gw *vector;
	//unsigned char *compressed_data;
	//int compressed_data_dimension;
	struct _data *childZero;
	struct _data *childOne;
	char *alphabetZero;
	char *alphabetOne;
	int alphabet_size_zero;
	int alphabet_size_one;
} wavelet_node_gw;


wavelet_node_gw * create_tree(unsigned char [],int ,char [], int , int, int);
wavelet_node_gw * create_tree2(unsigned char [],int ,char [], int , int, int, int*);
bool isvalueinarray(char , char [], int );
int posvalueinarray(char , char [], int);
int w_rank(char, int , wavelet_node_gw *);
char w_access(wavelet_node_gw *, int );
void free_wavelet_tree(wavelet_node_gw *);
