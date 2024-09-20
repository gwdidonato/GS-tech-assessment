BWT:
	gcc -O3 src/BWT.c -o BWT -lz -lm -w

MAP:
# Change path to correct and add -lz to link executable to zlib library. 
	gcc -O3 src/all_hw.c -o MAP -lz -lm -w 

clean:
	rm BWT MAP