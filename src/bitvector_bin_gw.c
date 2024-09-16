/*
Bitvector & RRR with bits

Third table is not used as it seems to be inefficient when using unsigned char
*/

#include "bitvector_bin_gw.h"
#include <stdio.h>
#include <time.h>

/*
Assumption: at the beginning each unsigned char is set to 0
*/
void insert(unsigned char * vector, int index){
	int shift_index = index % NUM_BITS_PER_BYTE;
	int vector_index = index / NUM_BITS_PER_BYTE;

	//vector[vector_index] |= (1 << 7-shift_index);
	vector[vector_index] |= (1 << shift_index);

}

int read_(unsigned char *vector, int index){
	int shift_index =  index % NUM_BITS_PER_BYTE;
	int vector_index = index / NUM_BITS_PER_BYTE;

	//return (vector[vector_index] >> (7-shift_index)) & 1;
	return (vector[vector_index] >> shift_index) & 1;
}

int rank(unsigned char *vector, int index, int dim){ //rank with popcount could me more efficient?
	int count = 0;
	for (int i = 0; i < index && i < dim; i++){
		if(read_(vector, i) == 1)
			count ++;
	}

	return count;
}

int rank_int32(unsigned int value, int index, int block_dim){ //rank with popcount could me more efficient?
	int count = 0;
	for (int i = block_dim-1; i > block_dim-1-index && i >= 0; i--){
		if((value>>i)&1)
			count ++;
	}
	return count;
}

unsigned int rank_interval(unsigned char *vector, int index_s, int index_e, int dim){
	unsigned int count = 0;

	for(int i = index_s; i < index_e && i < dim; i++){
		if(read_(vector,i) == 1)
			count++;
	}

	return count;
}


int rank_struct(bitvector_bin_gw * vector, int index, int dim){ //index è considerato a partire da 1...
	int count = 0;
	int i;
	if (index>dim)
		return vector -> rank_sum[dim/vector->l+1];

	if(index % vector -> l == 0){
		return vector -> rank_sum[index/vector->l];
	}else if(index % vector -> k == 0){
		for (i = index/vector->l * vector->l/vector->k; i < index/vector->k; i++){
			count += (int)vector->_class[i];
		}
		return vector -> rank_sum[index/vector->l] + count;
	}else{
		unsigned int posix_off = vector->pos_sum[index/vector->l];
		for (i = index/vector->l * vector->l/vector->k; i < index/vector->k; i++){
			count += (int)vector->_class[i];
			posix_off += ceil(log(binomialCoeff(vector->k,(int)vector -> _class[i]))/log(2));
		}

		unsigned int offset = 0;
		for(i = 0; i < ceil(log(binomialCoeff(vector->k,(int)vector -> _class[index/vector->k]))/log(2)); i++){
			offset |= (read_(vector->offset,posix_off+i) << (int)(ceil(log(binomialCoeff(vector->k,(int)vector -> _class[index/vector->k]))/log(2)) - 1 -i));
		}

		int c_off = vector->class_offset[vector->_class[index/vector->k]];
		count += rank_int32(vector->bit_table[c_off+offset], index%vector->k, vector->k);

		return vector -> rank_sum[index/vector->l] + count;
	}
}


/*
returns the position of the ith_one one, or -1 in case there are not so many 1.
*/
int select_bin(unsigned char *vector, int ith_one, int dim){
	int count = 0;
	for(int i = 0; i < dim; i++){
		if(read_(vector, i) == 1){
			count++;
		}
		if(count == ith_one){
			return i;
		}
	}

	return -1; //there is not
}



/*
return -1 in case not found, the position otherwise
*/

int binary_search(unsigned int *index, int val, int min, int max){
	int min_val = index[min];
	int max_val = index[max];

	if(val > max_val || val < min_val){
		return -1;
	}else if (max - min < 2){
		if(val == index[min]){
			return min;
		}else{
			return -1;
		}
	}else{
		int half = (min + max) / 2;
		if(val >=index[half]){
			return binary_search(index, val, half, max);
		}else{
			return binary_search(index, val, 0, half);
		}
	}
}

int binary_search_approx(unsigned int *index, int val, int min, int max){
	int min_val = index[min]; // lower bound
	int max_val = index[max];

	if(val > max_val || val < min_val){
		return -1;
	}else if(max - min < 2){
		if( val > index[min]){
			return min;
		}
		if(val <= index[min]){
			return min - 1;
		}
		/*if(val > index[min]){
			return min;
		}*/
	}

	int half = (min + max) / 2;
	if(val >= index[half]){
		return binary_search_approx(index, val, half, max);
	}else{
		return binary_search_approx(index, val, 0, half);
	}
}

int binary_search_approx_iter(unsigned int *index, int val, int min, int max){
	int it_min = min;
	int it_max = max;

	do {
	int min_val = index[it_min];
	int max_val = index[it_max];
	int half = (it_min + it_max) / 2;

		if(val > max_val || val < min_val){
			return -1;
		}else if(it_max - it_min < 2){
			if( val > index[it_min]){
				return it_min;
			}
			if(val <= index[it_min]){
				return it_min - 1;
			}
			/*if(val > index[min]){
				return min;
			}*/
		}

		if(val >= index[half]){
			it_min = half;
		}else{
			it_max=half;
		}
	}while(it_max - it_min >= 2);
	return -1;
}

unsigned int popcount(unsigned int value){
	unsigned int count = 0;
	for (count = 0; value != 0; value >>=1)
		if (value & 1)
			count ++;
	return count;
}


int binomialCoeff(int n, int k){
    int res = 1;

    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;

    for (int i = 0; i < k; ++i){
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
}


void swap(unsigned int *xp, unsigned int *yp){
    unsigned int temp = *xp;
    *xp = *yp;
    *yp = temp;
}


void bubbleSort_class(unsigned int _class[], unsigned int arr[], int end){
   	int i, j;
   	bool swapped;
   	for (i = 0; i < end; i++){
    	swapped = false;
     	for (j = 0; j < end-i; j++){
        	if (_class[j] > _class[j+1]){
           		swap(&arr[j], &arr[j+1]);
           		swap(&_class[j], &_class[j+1]);
           		swapped = true;
        	}
     	}

    	if (swapped == false)
        	break;
  	}
}


void bubbleSort_range(unsigned int arr[], int start, int end){
   int i, j;
   bool swapped;
   for (i = start; i < end; i++){
     swapped = false;
     for (j = start; j < end-i; j++){
        if (arr[j] > arr[j+1]){
           swap(&arr[j], &arr[j+1]);
           swapped = true;
        }
     }

     if (swapped == false)
        break;
   }
}


unsigned int* populate_bits(unsigned int* bits, int block_dim){

	int dim2 = pow(2,block_dim);
	unsigned int pcount[dim2];
	unsigned int i;
	for(i = 0; i < dim2; i++){
		bits[i]=i;
		pcount[i]=popcount(i);
	}

	//sort per class
	bubbleSort_class(pcount, bits, dim2-1);

	int start = 0;
	int end = 0;

	//sort inside class
	for(i=1; i<=block_dim; i++){
		start = 1;
		for (int j = 0; j < i; j++)
			start += binomialCoeff(block_dim,j);
		end = start - 1 + binomialCoeff(block_dim,i);
		bubbleSort_range(bits, start, end);
	}

/* //testing correctness
	for (int i = 0; i < dim2; ++i){
		printf("%d\t", pcount[i]);
		for (int j = block_dim - 1; j >= 0; j--){
			printf("%d ", (int)((bits[i]>>j)&1));
		}
		printf("\n");
	}
*/
	return bits;

}


/*
This function populates the bitvector_bin structure.
Given an unsigned char composed of 0s and 1s, the dimension of the vector and the dimensions of the first and second block,
it creates the binary vector and the first and second indexes.
*/
bitvector_bin_gw * populate_bitvector_bin_gw_rank(bitvector_bin_gw * vector, unsigned char *structure, int dimension, int superblock_dim, int block_dim){
	//allocate memory for the structure
	int i,j;
	unsigned long long int memory = 0;
	unsigned long long int shared = 0;
	if(vector = (bitvector_bin_gw *)malloc(sizeof(bitvector_bin_gw))){
		memory+=(int)sizeof(bitvector_bin_gw);
		vector -> l = superblock_dim;
		vector -> k = block_dim;
		vector -> dim = dimension;

		//allocate rank_sum and class
		if(vector -> rank_sum = (unsigned int *)calloc((dimension/superblock_dim+2), sizeof(unsigned int))){
			memory+=(int)((dimension/superblock_dim+2)*sizeof(unsigned int));
			if(vector -> _class = (unsigned char *)calloc((dimension/block_dim+1), sizeof(unsigned char))){
				memory+=(int)((dimension/block_dim+1)*sizeof(unsigned char));

				//populate rank_sum
				vector -> rank_sum[0] = 0;
				for(i = 1; i <= dimension/superblock_dim +1; i++)
					vector -> rank_sum[i] = vector -> rank_sum[i-1] + rank_interval(structure, (i-1)*superblock_dim, (i)*superblock_dim, dimension);

				//populate class
				if(superblock_dim % block_dim != 0){
					printf("Error, second block is not multiple of the first \n");
					exit(EXIT_FAILURE);
				}

				for(i = 0; i < dimension / block_dim +1; i++){
					vector -> _class[i] = (unsigned char)rank_interval(structure, i*block_dim, (i+1)*block_dim, dimension);
				}

				//calculate offsets bitvector length in bits
				int len_offset = 0;
				for (i = 0; i < dimension/block_dim +1; ++i)
					len_offset += ceil(log(binomialCoeff(block_dim,(int)vector -> _class[i]))/log(2));
				printf("\nlen_offset = %d bits \n\n", len_offset);
				//allocate and initialize offsets bitvector
				if(vector -> offset = (unsigned char *)calloc(len_offset/8 + 1, sizeof(unsigned char))){
					memory+=(int)((len_offset/8 + 1) * sizeof(unsigned char));
					for (i = 0; i < len_offset/8+1; ++i)
						vector->offset[i] = 0;
					//allocate pos_sum
					if(vector -> pos_sum = (unsigned int *)calloc((dimension/superblock_dim + 1), sizeof(unsigned int))){
						memory+=(int)((dimension/superblock_dim + 1)* sizeof(unsigned int));
						//allocate class_offset
						if(vector -> class_offset = (unsigned int*)calloc(block_dim+1, sizeof(unsigned int))){
							shared+=(int)((block_dim+1)* sizeof(short int));

							//compute class offsets
							clock_t cl_off0 = clock();
							vector -> class_offset[0] = 0;
							for (i = 1; i < block_dim+1; ++i)
								vector->class_offset[i] = vector->class_offset[i-1] + binomialCoeff(block_dim,i-1);
							clock_t cl_off1 = clock();
							printf("class_offset time: %f s\n\n", (double)(cl_off1 - cl_off0)/CLOCKS_PER_SEC);


							//allocate bitvectors table
							if(vector -> bit_table = (unsigned int*)calloc(pow(2,block_dim), sizeof(unsigned int))){
								shared+=(int)(pow(2,block_dim)* sizeof(unsigned int));

								//populate bitvectors table
								clock_t bit0 = clock();
								vector -> bit_table = populate_bits(vector -> bit_table, block_dim);
								clock_t bit1 = clock();
								printf("bit_table time: %f s\n\n", (double)(bit1 - bit0)/CLOCKS_PER_SEC);

								unsigned int temp_bv = 0;
								unsigned int posix = 0;
								unsigned int temp_off = 0;
								//populate offsets bitvector and pos_sum
								for (i = 0; i < dimension/block_dim+1; ++i){
									temp_bv = 0;
									temp_off = 0;
									for (j = 0; j < block_dim && j < dimension - i * block_dim; j++){
										temp_bv |= read_(structure,i*block_dim+j) << (block_dim - 1 - j);
									}
									//printf("temp_bv(%d) %d\n", i, temp_bv); //checked: it's correct!

									for (j = vector->class_offset[vector->_class[i]]; j < vector->class_offset[vector->_class[i]+1]; j++){
										if (vector->bit_table[j] == temp_bv){
											temp_off = j - vector->class_offset[vector->_class[i]];
											//printf("temp_of(%d) %d\n", i, temp_off); //checked: it's correct!
											for (int k = 0; k < ceil(log(binomialCoeff(block_dim,(int)vector -> _class[i]))/log(2)); k++){
												if((temp_off >> (int)(ceil(log(binomialCoeff(block_dim,(int)vector -> _class[i]))/log(2))-k-1)) & 1)
													insert(vector->offset,k+posix);
											}
											break;
										}
									}
									//populate pos_sum
									if(!(i%(superblock_dim/block_dim)))
										vector->pos_sum[i/(superblock_dim/block_dim)]=posix;

									posix += ceil(log(binomialCoeff(block_dim,(int)vector -> _class[i]))/log(2));

								}

							}else{
								printf("Error while allocating bit_table \n");
								exit(EXIT_FAILURE);
							}
						}else{
							printf("Error while allocating class_offset \n");
							exit(EXIT_FAILURE);
						}
					}else{
						printf("Error while allocating pos_sum \n");
						exit(EXIT_FAILURE);
					}
				}else{
					printf("Error while allocating offset \n");
					exit(EXIT_FAILURE);
				}
			}else{
				printf("Error while allocating class index \n");
				exit(EXIT_FAILURE);
			}
		}else{
			printf("Error while allocating rank_sum index \n");
			exit(EXIT_FAILURE);
		}
	}else{
		printf("Error while allocating bitvector \n");
		exit(EXIT_FAILURE);
	}
	printf("Memory Occupancy = %llu\n", memory);
	printf("Shared Memory = %llu\n", shared);
	return vector;
}


void free_bitvector_gw(bitvector_bin_gw * vector){
	free(vector -> rank_sum);
	free(vector -> pos_sum);
	free(vector -> _class);
	free(vector -> offset);
	free(vector);
}
