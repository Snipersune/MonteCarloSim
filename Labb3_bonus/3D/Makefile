.DEFAULT_GOAL := perc
INCLUDES = -I../g2/include -L../g2
LIBS = -lm -lg2 -lgd -lX11
CFLAGS = -std=c99 -DTHREE_D -O3

SRC = ../src/ran.c ../src/perc.c ../src/my_queue.c ../src/utils.c

perc:
	gcc $(SRC) -o perc $(CFLAGS) $(INCLUDES) $(LIBS)

summary:
	gcc ../src/summary.c -o summary $(CFLAGS) $(INCLUDES) $(LIBS)