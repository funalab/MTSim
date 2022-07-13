PROG = mtsim
HEADERS = $(wildcard *.h)
SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(wildcard *.c))

CC = gcc

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
		CFLAGS += -Wall -O0 -g -I/opt/local/include
		LDFLAGS = -L/opt/local/lib -lX11 -lm
else
		CFLAGS += -Wall -O0 -g -I/usr/X11R6/include
		LDFLAGS = -L/usr/X11R6/lib -lX11 -lm
endif

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
	$(CC) $(CFLAGS) -c $<

$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(PROG) $(OBJS)

TAGS: $(SRCS) $(HEADERS)
	ctags -Re *.c *.h
