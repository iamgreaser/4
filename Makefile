
BINNAME = 4

OBJS = box.o level.o vec.o main.o
INCLUDES = common.h

CFLAGS = -fopenmp -g -Wall -O2 -msse2 `sdl-config --cflags` $(CFLAGS_EXTRA)
LDFLAGS = -fopenmp -g -Wall -lm `sdl-config --libs` $(LDFLAGS_EXTRA)

RM = rm
RM_F = $(RM) -f

all: $(BINNAME)

clean:
	$(RM_F) $(OBJS)

clean-bin:
	$(RM_F) $(BINNAME)

$(BINNAME): $(OBJS)
	$(CC) -o $(BINNAME) $(OBJS) $(LDFLAGS)

%.o: %.c $(INCLUDES)
	$(CC) -c -o $@ $(CFLAGS) $<

