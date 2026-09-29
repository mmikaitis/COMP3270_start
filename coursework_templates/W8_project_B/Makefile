CC = gcc
MPICC = mpicc
CFLAGS = -O2 -std=c99 -Wall -Wextra
LDFLAGS = -lm
COMMON = fft.c
HEADERS = fft.h

.PHONY: all serial clean
all: fft_serial fft_mpi
serial: fft_serial

fft_serial: fft_serial.c $(COMMON) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ fft_serial.c $(COMMON) $(LDFLAGS)

fft_mpi: fft_mpi.c $(COMMON) $(HEADERS)
	$(MPICC) $(CFLAGS) -o $@ fft_mpi.c $(COMMON) $(LDFLAGS)

clean:
	rm -f fft_serial fft_mpi fft_serial.exe fft_mpi.exe
