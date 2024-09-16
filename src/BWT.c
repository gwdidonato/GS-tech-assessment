#include "kseq.h"
#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <time.h>
KSEQ_INIT(gzFile, gzread)

typedef struct{
	char *a;	//stringa in cui andrò a memorizzare il genoma, su cui riscriverò l'output della BWT//
	int n; 		//dimensione effettiva del genoma che l'utente andrà ad inserire//
	int k;		//la variabile in cui salvo la lunghezza dei k-meri.
	int *suffix;//l'array in cui salverò il suffix array 
} genoma;

void inserisci_genoma(genoma*, char*);
void de_brujin (genoma*, FILE*);
void make_k (genoma*, char*);
void quickSort( genoma*, char*, char*, char*, int, int);
int partition(genoma*, char*, char*, char*, int, int);
void stampa_matrice(genoma*, char*);
void output(genoma*, char*, FILE*);

int main(int argc, char *argv[]){
	genoma d;
	// popolo il genoma
	inserisci_genoma(&d, argv[1]);

	//leggo parametro k
	sscanf (argv[2],"%d",&(d.k));

	// apro file di output
	FILE *fd;
	fd=fopen(argv[3], "w");
	if(fd==NULL){
		perror("Errore nell'apertura del file di output");
		exit(1);
	}

	// calcolo e salvo il suffix array
	de_brujin (&d, fd);

	// chiudo il file di output
	fclose(fd);

	return 0;
}

// funzione per la lettura del genoma
void inserisci_genoma (genoma* d, char* filename){
	int i,l;
	gzFile stringa;
	kseq_t *seq;

	// apro file di input
	stringa = gzopen(filename, "r");

	if(stringa==NULL){ //Controllo
		perror("Errore nell'apertura del file\n");
		exit(1);
	}

	// leggo genoma in seq tramite kseq
	seq = kseq_init(stringa);
	l = kseq_read(seq);

	// aggiungo 1 pe considerare carattere speciale '$'
	d->n=(int)(seq->seq.l) + 1;

	// alloco spazio per la stringa
	d->a = (char*)(malloc((d->n) * sizeof(char)));

	// alloco spazio per suffix array
	d->suffix = (int*)(malloc((d->) * sizeof(int)));

	// copio in a la sequenza in seq.s
	strncpy(d->a, seq->seq.s, (d->n)-1);
	// aggiungo il carattere speciale in ultima posizione
	d->a[(d->)-1]='$';

	// libero memoria in seq e chiudo il file.
	kseq_destroy(seq);
	gzclose(stringa);

	// inizializzo il suffix array (numerazione 1-based invece che 0-based)
	for(i=1; i<(d->n)+1; i++)
		d->suffix[i-1] = i;
}

// funzione che popola una matrice mat (n righe, k+1 colonne) linearizzata
// ogni riga di mat è una finestra di k elementi presi scorrendo sulla sequenza genomica
// per ogni riga, l'ultimo elemento è il carattere antecedente la finestra campionata
void make_k (genoma* d, char* mat){
	int i, j;

	// popolo separatamente la prima riga perché l'ultimo carattere della prima riga è l'ultimo carattere dell'intera sequenza
	for(j=0; j<(d->k); j++)
		*(mat+j)=d->a[j];
	*(mat+(d->k))=d->a[(d->n)-1]; //in coda, l'elemento precedente

	// popolo tutte le altre righe della matrice
	for(i=1; i<(d->n); i++){
		for(j=0; j<(d->); j++)
			*(mat+i*(d->k+1)+j)=d->a[(i+j)%(d->n)]; //(i+j)%(d->n) perchè quando finisco la stringa le lettere devo prenderle dall'inizio
		*(mat+i*(d->k+1)+(d->k))=d->a[(i-1)];
	}
}

// funzione che crea il suffix array ordinando lessicograficamente la matrice
void de_brujin (genoma* d, FILE* fd){
	// alloco la matrice n*(k+1)
	char *mat = (char*)(malloc(sizeof(char) * ((d->n))*((d->k)+1)));
	// alloco stringa temporanea necessaria per algoritmo di sorting
	char *tmp = (char*)(malloc(sizeof(char) * (d->k)+1));
	// alloco stringa pivot necessaria per algoritmo di sorting
    char *pivot = (char*)(malloc(sizeof(char)* (d->k)+1));

	// creo la matrice mat
	make_k(d, mat);

	// stampo la matrice per controllo
	//stampa_matrice (d, mat);

	clock_t ordina0 = clock();
	// ordino la matrice mat tramite algoritmo quickSort
	quickSort (d, mat, tmp, pivot, 0, d->n-1);
	clock_t ordina1 = clock();

	// controllo esistenza file di output
	if(fd==NULL){
		perror("\n\nErrore nel passaggio del file\n\n");
		exit(2);
	}

	// stampo a video tempo di sorting
	printf("Sorting time: %f s\n\n", (double)(ordina1-ordina0)/CLOCKS_PER_SEC);

	// stampo la matrice ordinata per controllo
	//stampa_matrice (d, mat);

	clock_t output0 = clock();
	// salvo l'output nel rispettivo file
	output (d, mat, fd);
	clock_t output1 = clock();

	printf("Output time: %f s\n\n", (double)(output1-output0)/CLOCKS_PER_SEC);
	printf("\n\n\t\t\t\t********************************\n\n\n");

	// libero memoria allocata
	free(met);
	free(pivot);
   	free(tmp);
}

// funzione per stampare a video la matrice
void stampa_matrice (genoma* d, char* mat){
	int i, j;

	for(i=0; i<(d->n); i++){
		printf("Riga %d ", d->suffix[i]);
		for(j=0; j<(d->k+1); j++)
			printf("%c ", *(mat+i*(d->k+1)+j));
		printf("\n");
	}
	printf("\n\n");
}

// funzione che ordina alfabeticamente la matrice precedente, sovrascrivendola attraverso l'uso di una stringa di appoggio
// (si può dare per assodato sia corretta)
void quickSort( genoma *d, char* mat , char* tmp, char* pivot, int l, int r){
	int j;

	if( l < r ){
		// dividi et impera - continuo a chiamare quicksort fintanto che gli elementi da ordinare non sono lunghi un solo carattere
    	j = partition(d, mat, tmp, pivot, l, r);
    	quickSort(d, mat, tmp, pivot, l, j-1);
    	quickSort(d, mat, tmp, pivot, j+1, r);
	}

}

// funzione che implementa la vera logica di quickSort 
// (si può dare per assodato sia corretta)
int partition(genoma* d, char* mat, char* tmp, char* pivot, int l, int r) {
	int i, j, sx, h, z, ok;

	strncpy(pivot, mat+l*(d->k+1), (d->k+1));

	i = l;
	j = r+1;

	while(1)
	{ //Il do while viene sempre eseguito almeno una volta. Prima esegue il corpo delle istruzioni,
	//e poi testa la condizione nel while. Questo significa che il primo elemento, che è anche il pivot
	//non viene mai testato e che il primo che viene testato dall'indice j è l'ultima riga.

	ok=1;
	do{
		++i;
		if(strncmp(mat+i*(d->k+1), pivot, d->k)==0 && i <= r ){
			h=(d->suffix[i])%(d->n);
			z=(d->suffix[l])%(d->n);
			while (d->a[h] == d->a[z]){ //Finchè non trovi due caratteri lesicograficamente diversi, incrementa
				h++;
				z++;
			}
			if(d->a[h] > d->a[z]){ //Se il carattere per cui differiscono è tale da invertire le righe, invertile
				ok=0;
			}
		}
		else if(strncmp(mat+i*(d->k+1), pivot, d->k)>0)
			ok=0;
	}while( ok && i < r );

	ok=1;
	do{
		--j;
		if(strncmp(mat+j*(d->k+1), pivot, d->k)==0 && j > l){
			h=(d->suffix[j])%(d->n);
			z=(d->suffix[l])%(d->n);
			while (d->a[h] == d->a[z]){ //Finchè non trovi due caratteri lesicograficamente diversi, incrementa
				h++;
				z++;
			}
			if(d->a[h] < d->a[z]){ //Se il carattere per cui differiscono è tale da invertire le righe, invertile
				ok=0;
			}
		}
		if(strncmp(mat+j*(d->k+1), pivot, d->k)==0 && j==l)
			ok=0;
		else if(strncmp(mat+j*(d->k+1), pivot, d->k)<0)
			ok=0;
	}while( ok && j > l );


	if( i >= j ) break;

		else{
			strncpy(tmp, mat+i*(d->k+1), (d->k+1));
			strncpy(mat+i*(d->k+1), mat+j*(d->k+1), (d->k+1));
			strncpy(mat+j*(d->k+1), tmp, (d->k+1));

				sx = d->suffix[i];
				d->suffix[i] = d->suffix[j];
				d->suffix[j] = sx;
		}
	}

	if ( l != j ){
		strncpy(tmp, mat+l*(d->k+1), (d->k+1));
		strncpy(mat+l*(d->k+1), mat+j*(d->k+1), (d->k+1));
		strncpy(mat+j*(d->k+1), tmp, (d->k+1));
	}

	sx = d->suffix[l];
	d->suffix[l] = d->suffix[j];
	d->suffix[j] = sx;

	return j; //è la nuova posizione del pivot. Il pivot è l'unico elemento che di volta in volta viene
	//ordinato.
}

// funzione per la stampa a video dell'output
// il file di output contiene nella prima riga la BWT, e nella seconda il suffix array (le posizioni di ciascuna lettera della BWT nel genoma iniziale)
void output (genoma* d, char* mat, FILE* fd){
	int i, k, x;

	// per ciascuno degli n elementi del genoma (a), sovrascrivo inserendo il relativo carattere della BWT (ultimo carattere della rispettica riga)
	for(i=0; i<(d->); i++)
		d->a[i]=*(mat+i*(d->k+1)+d->k);

	if(fd==NULL){
		perror("Errore nel passaggio del file");
		exit(3);
	}

	// scrivo BWT
	fprintf(fd, "%s", d->a);
	fprintf(fd, "\n");

	//scrivo Suffix Array
	for(i=0; i< ; i++)
		fprintf(fd, "%d ", d->suffix[i]);

}