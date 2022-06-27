PROG = mtsim
HEADERS = $(wildcard *.h)
SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(wildcard *.c))
CC = gcc
#CFLAGS = -Wall -Wextra -O0 -g -I/usr/X11R6/include
CFLAGS = -Wall -O0 -g -I/usr/X11R6/include
LDFLAGS = -L/usr/X11R6/lib -lX11 -lm

debug := false
ifeq ($(debug),true)
	CFLAGS += -DDEBUG_PRINT
else ifeq ($(debug),false)
	CFLAGS += 
else
 $(error debug=... must be "true" or "false")
endif

.PHONY: all
all: TAGS $(PROG)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -DRadS=$(CRADS) -c $<

$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean: 
	rm -f $(PROG) $(OBJS)

TAGS: $(SRCS) $(HEADERS)
	ctags -Re *.c *.h
