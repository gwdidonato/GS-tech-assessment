BWT:
	gcc -O3 src/BWT.c -o BWT -lz -lm -w

#ERROR: not using zlib library
#SOLUTION: add -lz to compile
MAP:
	gcc -O3 src/all_hw.c -o MAP -lz -lm -w

clean:
	rm BWT MAP
