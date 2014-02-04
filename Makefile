# Copyright 2010, 2011 Paul Chote
# This file is part of Puoko-nui, which is free software. It is made available
# to you under the terms of version 3 of the GNU General Public License, as
# published by the Free Software Foundation. For more information, see LICENSE.

CC = gcc
LINKER = gfortran
CFLAGS = -g -c -Wall -pedantic -Dlinux --std=c99 -D_XOPEN_SOURCE -D_BSD_SOURCE
LFLAGS = -lcfitsio -lxpa -lgsl

LFLAGS += -L/usr/local/lib -L/usr/X11R6/lib -lX11 -lcpgplot -lpgplot -lpng


SRC = main.c
OBJ = $(SRC:.c=.o)


harmonics: $(OBJ)
	$(LINKER) -o $@ $(OBJ) $(LFLAGS)

clean:
	-rm $(OBJ) harmonics

.SUFFIXES: .c
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
