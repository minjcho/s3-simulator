CC      = gcc
CFLAGS  = -Wall -Wextra -O2 -std=c11
LDFLAGS = -lm

TARGET  = mini_s3
SRCS    = main.c gf256.c erasure.c
OBJS    = $(SRCS:.c=.o)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

main.o: main.c gf256.h erasure.h
erasure.o: erasure.c erasure.h gf256.h
gf256.o: gf256.c gf256.h

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean run
