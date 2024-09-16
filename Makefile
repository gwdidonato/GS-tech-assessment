BWT:
	gcc -O3 src/BWT.c -o BWT -lz -lm -w

MAP:
	gcc -O3 all_hw.c -o MAP -lm -w

clean:
	rm BWT MAP