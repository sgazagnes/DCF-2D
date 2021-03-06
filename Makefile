# --------------------------------------------# Compiler options

CC=mpicc
CFLAGS  =  -pedantic -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wstrict-prototypes -Wwrite-strings -Wcast-qual -Wswitch-default -Wunreachable-code -Wformat=2 -Winit-self -march=native -g -std=gnu99
LDFLAGS	= 

LIBS    = -lm -lfreeimage -lhdf5 -lcfitsio -pg -L /home/simon/lib/lib  
COMPILE = $(CC) -c $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) $(LIBS)
LINK    = $(CC) $(LDFLAGS) -o $@



# -----------------------------------------   Folders
SRC=src
BUILD=build
OBJ=$(BUILD)/obj
COMMON=$(SRC)/common
AREA=$(SRC)/area
DAP=$(SRC)/dap
PAT=$(SRC)/pattern

# --------------------------------------------
all: area dap pattern

area: $(COMMON)/*.c $(AREA)/*.c 
	$(CC) $(CFLAGS) -I /home/simon/lib/include -I $(COMMON) -I $(AREA)  -o $@ $^  $(LIBS)

dap:  $(COMMON)/*.c $(DAP)/*.c 
	$(CC) $(CFLAGS) -I /home/simon/lib/include -I $(COMMON) -I $(DAP)  -o $@ $^  $(LIBS)

pattern:  $(COMMON)/*.c $(PAT)/*.c 
	$(CC) $(CFLAGS) -I /home/simon/lib/include -I $(COMMON) -I $(PAT)  -o $@ $^  $(LIBS)

mpitest: mpitest.c
	$(CC) $(CFLAGS) -o mpitest mpitest.c  $(LIBS)


clean:
	rm area dap pattern
