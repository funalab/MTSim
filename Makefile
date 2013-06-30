PROG = mtsim
HEADERS = $(wildcard *.h)
SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(wildcard *.c))
#OBJS = color_setting.o display_setting.o draw_graphs.o free.o lubksb_double.o ludcmp_double.o math_func.o nrutil.o ran1_double.o save_logs.o solve_1D_double.o store_graphs.o main.o
CC = gcc
#CFLAGS = -Wall -Wextra -O0 -g -I/usr/X11R6/include
#CFLAGS = -DDEBUG_PRINT -Wall -O0 -g -I/usr/X11R6/include
CFLAGS = -Wall -O0 -g -I/usr/X11R6/include
LDFLAGS = -L/usr/X11R6/lib -lX11 -lm

.PHONY: all
all: TAGS $(PROG)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $<

$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean: 
	rm -f $(PROG) $(OBJS)

TAGS: $(SRCS) $(HEADERS)
	ctags -Re *.c *.h
