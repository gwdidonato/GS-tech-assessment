
#include "wavelet_tree_gw.h"

wavelet_node_gw * create_tree2(unsigned char *original_letters,int original_letters_dim, char alphabet[], int first_block_dim, int second_block_dim, int alphabet_size, int* dollar_posix){

	if(alphabet_size <2){
		return NULL;
	}
	int i = 0;
	int j;
	wavelet_node_gw *temp = NULL;

	if(!(temp = (wavelet_node_gw *) malloc(sizeof(wavelet_node_gw)))){
		printf("ERROR: Unable to allocate memory for the wavelet_node \n");
		exit(EXIT_FAILURE);
	}

	//memory allocated here, the children are init to null
	temp -> childZero = NULL;
	temp -> childOne = NULL;

	//now populate the vector and create its indexes..
	int counter = 0;
	for(i = 0; i < original_letters_dim; i++){
		if(isvalueinarray(original_letters[i], alphabet, alphabet_size) == true){
			counter++;
		}
		else if(original_letters[i]==4){
			*dollar_posix=i;
		}
	}

	int dimension = counter;

	//set who's zero and who's one
	/*
		this works only with alphabet's sizes 2^x
	*/
	if(!(temp -> alphabetZero = (char*)malloc(sizeof(char) * alphabet_size/2))){
		printf("Error allocating\n");
		exit(EXIT_FAILURE);
	}
	if(!(temp -> alphabetOne = (char*)malloc(sizeof(char) * alphabet_size/2))){
		printf("Error allocating\n");
		exit(EXIT_FAILURE);
	}

	temp -> alphabet_size_zero = 0;
	for(i = 0; i < alphabet_size/2; i++){
		temp -> alphabetZero[i] = alphabet[i];
		temp -> alphabet_size_zero++;
	}

	temp -> alphabet_size_one = 0;
	for(i = 0, j = alphabet_size/2; j < alphabet_size; i++, j++){
		temp -> alphabetOne[i] = alphabet[j];
		temp -> alphabet_size_one++;
	}

	unsigned char * compressed_data;
	if(!(compressed_data = (unsigned char *) calloc(dimension/ NUM_BITS_PER_BYTE + 1, sizeof(unsigned char)))){
		printf("ERROR while allocating compressed \n");
		exit(EXIT_FAILURE);
	}
	for(i = 0, j = 0; i < original_letters_dim; i++){
		if(isvalueinarray(original_letters[i], temp -> alphabetZero, temp -> alphabet_size_zero) == true){
			j++;
			//printf("%c ", original_letters[i]);
		}
		else if(isvalueinarray(original_letters[i], temp -> alphabetOne, temp -> alphabet_size_one) == true){
			insert(compressed_data, j);
			j++;
			//printf("%c ", original_letters[i]);
		}
	}
	int compressed_data_dimension = dimension;

	/*printf("\n dimension %d comnpressed string: \n",dimension);
	for(int i = 0; i < dimension; i++)
		printf("%d ", read_(temp -> compressed_data,i));

	printf("\n");*/

	temp -> vector = populate_bitvector_bin_gw_rank(temp -> vector, compressed_data, dimension, first_block_dim, second_block_dim);

	temp -> childZero = create_tree2(original_letters, original_letters_dim,temp->alphabetZero, first_block_dim, second_block_dim, temp->alphabet_size_zero, dollar_posix);
	temp -> childOne = create_tree2(original_letters, original_letters_dim,temp->alphabetOne, first_block_dim, second_block_dim, temp->alphabet_size_one, dollar_posix);

	return temp;
}

bool isvalueinarray(char val, char *arr, int size){
    int i;
    for (i=0; i < size; i++) {
        if (arr[i] == val)
            return true;
    }
    return false;
}

/*
check if a value is present in the array. In positive case, it returns the position
-1, otherwise.
*/
int posvalueinarray(char val, char *arr, int size){
	int i=0;
	for(i = 0; i < size; i++){
		if(arr[i] == val)
			return i;
	}
	return -1;
}

int w_rank(char letter, int pos, wavelet_node_gw *head){
	//1. what character is represented by letter in the first level?
	wavelet_node_gw *temp;
	int rank = pos;

	for(temp = head; temp != NULL; /*update below*/){
		short flag_one = 1;
		int bit = posvalueinarray(letter,temp->alphabetOne, temp->alphabet_size_one);
		if(bit == -1){
			bit = posvalueinarray(letter, temp->alphabetZero, temp->alphabet_size_zero);
			flag_one = 0;
		}

		if(bit != 0 && bit != 1){
			printf("Can't find letter in alphabet \n");
			return -1;
		}


		//rank the value in the first level.
		int rank_temp =  rank_struct(temp-> vector, rank, temp->vector->dim);

		if(flag_one == 0){
			rank_temp = rank - rank_temp;
		}

		rank = rank_temp;
		//printf("RANK: %d \n", rank);

		if(flag_one == 1){
			temp = temp -> childOne;
		}
		else{
			temp = temp -> childZero;
		}


	}

	return rank;

}

void free_wavelet_tree(wavelet_node_gw *head){
	if(head == NULL)
		return;

	free_wavelet_tree(head -> childOne);
	free_wavelet_tree(head -> childZero);

	free_bitvector_gw(head->vector);
	//free(head->compressed_data);
	free(head->alphabetOne);
	free(head->alphabetZero);

	return;

}
