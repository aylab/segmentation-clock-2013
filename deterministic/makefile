CC = g++
CFLAGS = -Wall -g
TARGET = deterministic
OBJS = main.o input_functions.o functions.o output_functions.o

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm *.o *~

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)