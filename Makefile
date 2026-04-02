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

TEST    = mini_s3_test

main.o: main.c gf256.h erasure.h
erasure.o: erasure.c erasure.h gf256.h
gf256.o: gf256.c gf256.h
test.o: test.c gf256.h erasure.h

$(TEST): test.o gf256.o erasure.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

test: $(TEST)
	./$(TEST)

clean:
	rm -f $(OBJS) test.o $(TARGET) $(TEST)

.PHONY: clean run test
