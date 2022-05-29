#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "my_queue.h"


int MAX_size;
int *intArray;
int front = 0;
int rear = -1;
int itemCount = 0;

void setQueueMAXSize(int size){
    MAX_size = size;
}

int init_queue(){
    intArray = malloc(MAX_size*sizeof(int));
    if(intArray == NULL){
        fprintf(stderr, "Could not initiate queue! Program terminated.\n");
        return 0;
    }
    return 1;
}

void free_queue(){
    free(intArray);
}

int peek() {
   return intArray[front];
}

int isEmpty() {
   return itemCount == 0;
}

int isFull() {
   return itemCount == MAX_size;
}

int size() {
   return itemCount;
}  

void insert(int data) {

   if(!isFull()) {
	
      if(rear == MAX_size-1) {
         rear = -1;            
      }       

      intArray[++rear] = data;
      itemCount++;
   }
}

int removeData() {
   int data = intArray[front++];
	
   if(front == MAX_size) {
      front = 0;
   }
	
   itemCount--;
   return data;  
}

void clearQueue() {
   while (!isEmpty())
      removeData();
}