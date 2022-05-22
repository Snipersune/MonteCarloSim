INCLUDES = -I../g2/include -L../g2
LIBS = -lm -lg2 -lgd -lX11
CFLAGS = -std=c99 -DPLOT -O3 

SRC = ../src/ran.c ../src/config.c ../src/ising.c ../src/my_queue.c

ising:
	gcc $(SRC) -o ising $(CFLAGS) $(INCLUDES) $(LIBS)