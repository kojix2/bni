CC ?= cc
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
LIBDIR ?= $(PREFIX)/lib
INCLUDEDIR ?= $(PREFIX)/include

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

LIB_SRC = \
	src/bni_util.c \
	src/bni_format.c \
	src/bni_index.c \
	src/bni_get.c \
	src/bni_stats.c \
	src/bni_check.c

CLI_SRC = src/main.c
LIB_OBJ = $(LIB_SRC:.c=.o)
CLI_OBJ = $(CLI_SRC:.c=.o)
PIC_OBJ = $(LIB_SRC:.c=.pic.o)
BIN = bni
STATIC_LIB = libbni.a

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SHARED_LIB = libbni.dylib
SHARED_LDFLAGS = -dynamiclib -install_name @rpath/$(SHARED_LIB)
else
SHARED_LIB = libbni.so
SHARED_LDFLAGS = -shared -Wl,-soname,$(SHARED_LIB)
endif

ifneq ($(filter -static,$(LDFLAGS)),)
DEFAULT_LIBS = $(STATIC_LIB)
else
DEFAULT_LIBS = $(STATIC_LIB) $(SHARED_LIB)
endif

all: $(BIN) $(DEFAULT_LIBS)

$(BIN): $(CLI_OBJ) $(LIB_OBJ)
	$(CC) $(LDFLAGS) -o $@ $(CLI_OBJ) $(LIB_OBJ) $(LDLIBS)

$(STATIC_LIB): $(LIB_OBJ)
	ar rcs $@ $(LIB_OBJ)

$(SHARED_LIB): $(PIC_OBJ)
	$(CC) $(LDFLAGS) $(SHARED_LDFLAGS) -o $@ $(PIC_OBJ) $(LDLIBS)

$(LIB_OBJ) $(CLI_OBJ) $(PIC_OBJ): include/bni.h src/bni_internal.h

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

%.pic.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -fvisibility=hidden -DBNI_BUILD_SHARED -c -o $@ $<

install: all
	mkdir -p $(DESTDIR)$(BINDIR)
	mkdir -p $(DESTDIR)$(LIBDIR)
	mkdir -p $(DESTDIR)$(INCLUDEDIR)
	cp $(BIN) $(DESTDIR)$(BINDIR)/$(BIN)
	cp $(DEFAULT_LIBS) $(DESTDIR)$(LIBDIR)/
	cp include/bni.h $(DESTDIR)$(INCLUDEDIR)/bni.h

clean:
	rm -f $(LIB_OBJ) $(CLI_OBJ) $(PIC_OBJ) $(BIN) $(STATIC_LIB) $(SHARED_LIB)

.PHONY: all install clean
