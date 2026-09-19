CC = gcc-16
CFLAGS = -O2 -Wall
LDFLAGS = -lm

jacobi_solve: jacobi_solve.c system_gen.c system_gen.h
	$(CC) $(CFLAGS) -o $@ jacobi_solve.c system_gen.c $(LDFLAGS)

clean:
	rm -f jacobi_solve
