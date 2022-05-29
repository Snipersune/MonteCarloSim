.DEFAULT_GOAL := rwalk
INCLUDES = -I../g2/include -L../g2
LIBS = -lm -lg2 -lgd
CFLAGS = -std=c99 -O3

SRC = ../src/ran.c ../src/utils.c ../src/rwalk.c

rwalk:
	gcc $(SRC) -o rwalk $(CFLAGS) $(INCLUDES) $(LIBS)

summary:
	gcc ../src/summary.c -o summary $(CFLAGS) $(INCLUDES) $(LIBS)