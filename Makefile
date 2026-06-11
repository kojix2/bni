CC ?= cc
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

CPPFLAGS ?=
CFLAGS ?= -O2 -g -Wall -Wextra -std=c11
LDFLAGS ?=
LDLIBS ?=

CFLAGS += -D_POSIX_C_SOURCE=200809L -D_FILE_OFFSET_BITS=64 -Iinclude

ifdef HTSLIB_DIR
CPPFLAGS += -I$(HTSLIB_DIR)
LDFLAGS += -L$(HTSLIB_DIR)
LDLIBS += -lhts -lz -lbz2 -llzma -lm -lpthread
else
HTSLIB_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
HTSLIB_LIBS := $(shell pkg-config --libs htslib 2>/dev/null || echo -lhts -lz -lbz2 -llzma -lm -lpthread)
CPPFLAGS += $(HTSLIB_CFLAGS)
LDLIBS += $(HTSLIB_LIBS)
endif

SRC = \
	src/main.c \
	src/bni_util.c \
	src/bni_format.c \
	src/bni_index.c \
	src/bni_get.c \
	src/bni_stats.c \
	src/bni_check.c

OBJ = $(SRC:.c=.o)
BIN = bni

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

install: $(BIN)
	mkdir -p $(DESTDIR)$(BINDIR)
	cp $(BIN) $(DESTDIR)$(BINDIR)/$(BIN)

clean:
	rm -f $(OBJ) $(BIN)

.PHONY: all install clean
