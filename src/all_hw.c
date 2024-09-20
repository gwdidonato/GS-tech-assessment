#include "kseq.h"
#include "wavelet_tree_gw.c"

#define ALPHABET_SIZE 4 // A,C,G,T
#define N 5				// A,C,G,T,$

#define A 0
#define C 1
#define G 2
#define T 3


KSEQ_INIT(gzFile, gzread)

typedef struct{
	char * BWT; 		//BWT of the reference genome
	int *suffix_array; 	//suffix array (positions of the suffixes in the origal string)
	int n; 				//regerence length
	int c[N]; 			//FM_index1: for each character in the BWT alphabet ($), a counter of lexicographically smaller characters
} reference;

typedef struct{
	char *qseq; //query sequence
	int start; 	//index for backward search
	int end; 	//index for backward search
	int m; 		//query length
} query;


void insert_BWT (reference*,char*);
void compute_c(reference*);
void insert_reads(query*, kseq_t*);
void bw_search(reference*, query*, int*, wavelet_node_gw*, char*);
void locate(reference*, query*, FILE*);
void rev_comp(query*, int);
void locate_inv(reference*, query*, FILE*);


int main(int argc, char *argv[]){
	reference ref;
	query r;
	FILE *fresults;
	gzFile freads;
	kseq_t *seq;
	int i, l;
	int dollar_posix = 0;
	char alphabet[ALPHABET_SIZE] = {A, C, G, T};
	// default values for block and superblock
	int BLOCK_DIM = 15;
	int SUPERBLOCK_DIM = 150;

	// read block and superblock
	sscanf (argv[4],"%d",&BLOCK_DIM);
	sscanf (argv[4],"%d",&SUPERBLOCK_DIM);

	SUPERBLOCK_DIM *= BLOCK_DIM;
	printf("%d - %d\n", BLOCK_DIM, SUPERBLOCK_DIM);

	r.qseq=NULL;
	//populate reference struct
	insert_BWT(&ref,  argv[2]);
	//allocate wavelet tree
	wavelet_node_gw *head_tree = NULL;
	//populate wavelet tree of the BWT
	head_tree = create_tree2(ref.BWT, ref.n, alphabet, SUPERBLOCK_DIM, BLOCK_DIM, ALPHABET_SIZE, &dollar_posix);
	//compute c index
	compute_c(&ref);

	//open output file
	fresults = fopen(argv[3], "w");
	if(fresults==NULL){
		perror("Error while opening results file");
		exit(1);
	}

	//open queries file
	freads = gzopen(argv[1], "r");
	if(freads==NULL){
		perror("Error while opening reads file");
		exit(3);
	}

	seq = kseq_init(freads);
	//iteratively read queries and map them to the reference and its reverse complement
	while ((l = kseq_read(seq)) >= 0) {
		//load query
		insert_reads(&r, seq);

		//backward search on reference stored as wavelet tree
		bw_search(&ref, &r, &dollar_posix, head_tree, seq->name.s); // bw_search function takes as input a pointer to ref (&ref), not ref itself.
		//write the number of matches
		fprintf(fresults, "> %s | [%d bp] | + %d:\n", seq->name.s, (int)(seq->seq.l), (r.end - r.start +1));
		//locate the matches' positions
		locate(&ref, &r, fresults);
		//write the matches' positions
		fprintf(fresults, "\n");

		//compute reverse complement of the query, overwriting it
		rev_comp(&r,(int)(seq->seq.l));

		//backward search on reference stored as wavelet tree
		bw_search(&ref, &r, &dollar_posix, head_tree, seq->name.s); //Removed &head_tree. head_tree is already a pointer. &head_tree is the memory address of the pointer.
		//write the number of reverse matches
		fprintf(fresults, "> %s | [%d bp] | - %d:\n", seq->name.s, (int)(seq->seq.l), (r.end - r.start +1));
		//locate the reverse matches' positions
		locate_inv (&ref, &r, fresults);
		//write the reverse matches' positions
		fprintf(fresults, "\n\n");
	}

	fclose(fresults);

	//free reference
	free(ref.BWT);
	free(ref.suffix_array);

	//free query
	free(r.qseq);

	//free wavelet tree
	free_wavelet_tree(head_tree);

  	return 0;
}


// populate BWT from file
// (si può considerare corretta)
void insert_BWT (reference* ref, char* filename){
	FILE *fBWT;
	int i,j=0;
	int dim=0;
	fBWT = fopen(filename, "r"); // da file prendo l'output della BWT.

	if( fBWT==NULL){
		perror("Error while opening BWT file");
		exit(2);
	}

	while (getc(fBWT) != '\n')
		dim++;

	ref->n = dim;
	ref->BWT = (char*)(calloc(dim + 1, sizeof(char)));
	ref->suffix_array = (int*)(calloc(dim, sizeof(int)));

	fseek(fBWT, 0L, SEEK_SET);
	fscanf(fBWT, "%s", ref->BWT);

	while(!feof(fBWT) && j < dim){
		fscanf(fBWT, "%d ", &(ref->suffix_array[j]));
		j++;
	}


	for(i=0; i<ref->n; i++){
		switch(ref->BWT[i]){
			case('A'):
				ref->BWT[i]=A;
				break;
			case('C'):
				ref->BWT[i]=C;
				break;
			case('G'):
				ref->BWT[i]=G;
				break;
			case('T'):
				ref->BWT[i]=T;
				break;
			case('a'):
				ref->BWT[i]=A;
				break;
			case('c'):
				ref->BWT[i]=C;
				break;
			case('g'):
				ref->BWT[i]=G;
				break;
			case('t'):
				ref->BWT[i]=T;
				break;
			default:
				ref->BWT[i]=4;
				break;
		}
	}

	fclose(fBWT);

}

// for each char of the alphabet, count the number of lexicographically smaller elements
void compute_c (reference* ref){
	int i=0;

	//inizialization
	for(i=0; i<N; i++)
		ref->c[i]=0;
	//count occurrences for each symbol in the BWT
	for(i=0; i<ref->n; i++)
		++(ref->c[ref->BWT[i]]); // Add BWT to loop through each character of the transformation and update the count.

	//for each character in the BWT alphabet ($), count lexicographically smaller characters
	ref->c[T]=(ref->c[G])+(ref->c[C])+(ref->c[A])+1;
	ref->c[G]=(ref->c[C])+(ref->c[A])+1;
	ref->c[C]=(ref->c[A])+1;
	ref->c[A]=1;
	ref->c[4]=0;
}

// insert reads using a number for each character, according to the defined alphabet
// (si può considerare corretta)
void insert_reads (query* r, kseq_t* seq){
	int i=0;
	char * old_b=r->qseq;
    // add 1 to compensate for possible not null termination
	r->qseq = (char*)(malloc((seq->seq.l * sizeof(char))+1));
    // if the old ptr is set, free the old one
	if (old_b!=NULL)
		free(old_b);
    // copy the string with at most seq->seq.l and add the final '\0'
	strncpy(r->qseq, seq->seq.s, seq->seq.l);
    r->qseq[seq->seq.l] = '\0';
	r->start=0;
	r->end=0;
	r->m = seq->seq.l;

	if(r->qseq[i] < 'a'){
		for(i=0; i<r->m; i++){
		switch(r->qseq[i]){
			case('A'):
				r->qseq[i]=A;
				break;

			case('C'):
				r->qseq[i]=C;
				break;

			case('G'):
				r->qseq[i]=G;
				break;

			case('T'):
				r->qseq[i]=T;
				break;
			}
		}
	}
	else{
		for(i=0; i<r->m; i++){
		switch(r->qseq[i]){

			case('a'):
				r->qseq[i]=A;
				break;

			case('c'):
				r->qseq[i]=C;
				break;

			case('g'):
				r->qseq[i]=G;
				break;

			case('t'):
				r->qseq[i]=T;
				break;
			}
        }
    }
}

// execute the backward search using the wavelet tree
// (si può considerare corretta)
void bw_search (reference* ref, query* r, int* dollar_posix, wavelet_node_gw* tree, char* name){
	int val;
	int i=1;

	// starting from last char
	val = r->qseq[r->m-1];

	// initializing start & end
	r->start = ref->c[val];
	if(val == T)
		r->end = ref->n-1;
	else
		r->end = ref->c[(val+1)]-1;

	// iterative backward search
	while(i < r->m && r->start <= r->end){
		i++;

		val = r->qseq[r->m - i];

		if (r->end<*dollar_posix){
			r->start = ref->c[val] + w_rank(val,r->start,tree);
			r->end = ref->c[val] + w_rank(val,r->end+1,tree)-1;
		}else{
			if (r->start<*dollar_posix)
				r->start = ref->c[val] + w_rank(val,r->start,tree);
			else
				r->start = ref->c[val] + w_rank(val,r->start-1,tree);
			r->end = ref->c[val] + w_rank(val,r->end,tree)-1;
		}
		if(r->end > (ref->n))
			r->end = (ref->n)-1;
	}
}

// print the positions of the matches (suffix array elements)
void locate (reference* ref, query* r, FILE* fresults){
	int i;

	if(fresults==NULL){
		perror("Errore nel passaggio del file");
		exit(6);
	}

	for(i = r->start; i <= r->end; i++){ // Add parenthesis
		fprintf(fresults, "\t%d", ref->suffix_array[i]); // Add ref to point to suffix array
	}
		
}

// compute the reverse complement of a string
void rev_comp(query* r, int dim){
	char tmp;
	int i;

	r->start=0;
	r->end=0;
	r->m = dim;

	// swap the elements
	for (i = 0; i < dim / 2; ++i){ // Do not loop through all dimension of query array, but only half.
		tmp=r->qseq[i];
		r->qseq[i]=r->qseq[r->m-1-i];
		r->qseq[r->m-1-i]=tmp;
	}
	// compute the reverse
	for (i = 0; i < dim; ++i){
		r->qseq[i]=3-(r->qseq[i]);
	}
}

// print the positions of the reverse matches (suffix array elements)
void locate_inv (reference* ref, query* r, FILE* fresults){
	int i;

	if(fresults==NULL){
		perror("Errore nel passaggio del file");
		exit(7);
	}

	for(i = r->start; i <= r->end; i++){
		fprintf(fresults, "\t%d", (ref->n)-(ref->suffix_array[i])-1); // Add file pointer and parenthesis.
	}
			
}