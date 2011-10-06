# Copyright 2010, 2011 Paul Chote
# This file is part of Puoko-nui, which is free software. It is made available
# to you under the terms of version 3 of the GNU General Public License, as
# published by the Free Software Foundation. For more information, see LICENSE.

CC = gcc
LINKER = gfortran
CFLAGS = -g -c -Wall -pedantic -Dlinux --std=c99 -D_XOPEN_SOURCE -D_BSD_SOURCE -I/usr/local/Cellar/gsl/1.15/include/gsl
LFLAGS = -lcfitsio -lxpa -lgsl

CFLAGS += -I/sw/lib/pgplot
LFLAGS += -L/usr/X11R6/lib -lX11 -L/sw/lib -laquaterm -Wl,-framework -Wl,Foundation -L/sw/lib/pgplot -lcpgplot -lpgplot -lpng


SRC = main.c
OBJ = $(SRC:.c=.o)


harmonics: $(OBJ)
	$(LINKER) -o $@ $(OBJ) $(LFLAGS)

clean:
	-rm $(OBJ) harmonics

.SUFFIXES: .c
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
