
BINNAME = 4

OBJS = box.o kd.o level.o render.o sphere.o vec.o main.o
INCLUDES = common.h

SDL_CF = `sdl-config --cflags` 
LIBS = -lm `sdl-config --libs` 
CFLAGS = -fopenmp -g -fno-strict-aliasing -Wall -O2 -msse2 -mfpmath=sse $(SDL_CF) $(CFLAGS_EXTRA)
LDFLAGS = -fopenmp -g -Wall $(LIBS) $(LDFLAGS_EXTRA)

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

