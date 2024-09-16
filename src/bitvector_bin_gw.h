#define BITVECTOR_BIN_GW_H
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <zlib.h>
#include <math.h>
#define NUM_BITS_PER_BYTE 8
#define member_size(type, member) sizeof(((type *)0)->member)

typedef struct _bitvector_bin{
	unsigned int *first;
	unsigned int *second;
	//unsigned int **third; //third table not used in this case.
	int l; //first block dim (should be multiple of k)
	int k; //second block dim
	int dim;
} bitvector_bin; //RRR data structure for bitvector (no compression, just fast rank)

typedef struct _bitvector_bin_gw{
	unsigned int *rank_sum;
	unsigned int *pos_sum;
	unsigned char *_class; // can be compressed
	unsigned char *offset; //bitvector -> fields of different length
	unsigned int *bit_table;
	unsigned int *class_offset;
	int l; //superblock dim (should be multiple of k)
	int k; //block dim
	int dim;
} bitvector_bin_gw; //RRR data structure for bitvector


void insert(unsigned char * , int );
int read_(unsigned char *, int);
int rank(unsigned char *, int, int);
int rank_int32(unsigned int , int , int );
int rank_struct(bitvector_bin_gw *, int, int);
unsigned int rank_interval(unsigned char *, int, int, int );
int select_bin(unsigned char *, int , int );
int binary_search(unsigned int *, int , int , int );
int binary_search_approx(unsigned int *, int , int , int );
int binary_search_approx_iter(unsigned int *, int , int , int );
int select_rank_based(unsigned char *, bitvector_bin_gw *, int );
unsigned int popcount(unsigned int );
int binomialCoeff(int , int );
void swap(unsigned int *, unsigned int *);
void bubbleSort_class(unsigned int *, unsigned int *, int );
void bubbleSort_range(unsigned int *, int , int );
unsigned int* populate_bits(unsigned int* , int );
bitvector_bin_gw * populate_bitvector_bin_gw_rank(bitvector_bin_gw *, unsigned char *, int , int , int );
void free_bitvector_gw(bitvector_bin_gw *);
