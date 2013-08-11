/*
4 (pronounced "sardinelauncher")
a first person shooter in 4 dimensions

Copyright (c) 2013 Ben "GreaseMonkey" Russell

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
claim that you wrote the original software. If you use this software
in a product, an acknowledgment in the product documentation would be
appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be
misrepresented as being the original software.

3. This notice may not be removed or altered from any source
distribution.

*/

// this header (or its stupidly-named equivalents) will be vital.
// if you don't have it... i honestly can't remember the really stupid names they use.
// i only remember that mmintrin is MMX and xmmintrin is the first SSE.
// pretty sure ymmintrin is NOT AVX (SSE2 perhaps?).
#include <immintrin.h>

// we'll need this for an extra speed boost on multicore systems :)
#include <omp.h>

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <errno.h>

#include <math.h>

#include <SDL.h>

#ifndef __SSE2__
#error "SSE2 is required to compile + run this. Use the -msse2 flag in GCC."
#endif


enum
{
	SHP_PAIR,

	SHP_ADD,
	SHP_SUB,
	SHP_DBL,
};

enum
{
	F_XN, F_YN, F_ZN, F_WN,
	F_XP, F_YP, F_ZP, F_WP,
};

__attribute__((aligned(16)))
typedef union v4f
{
	struct { float x, y, z, w; } v;
	struct { float b, g, r, a; } c;
	struct { float s, t, r, q; } t;
	float a[4];
	__m128 m;
} v4f_t;

// and by "box" we mean "tesseract"
typedef struct box box_t;
struct box
{
	v4f_t v0, v1; // min/max points in 4D space
	v4f_t color;
	box_t *c[2]; // children
	box_t *p; // parent
	int op, depth;
	int ptroffs;
};

__attribute__((aligned(16)))
typedef union m4f
{
	struct { v4f_t x, y, z, w; } v;
	struct { v4f_t b, g, r, a; } c;
	struct { v4f_t s, t, r, q; } t;
	v4f_t a[4];
} m4f_t;

__attribute__((aligned(16)))
typedef struct camera
{
	m4f_t m;
	v4f_t o;
} camera_t;

extern int fps_counter;
extern int fps_next_tick;
extern box_t *root;

extern SDL_Surface *screen;
extern uint32_t *rtbuf;
extern int rtbuf_width;
extern int rtbuf_height;
extern int rtbuf_scale;
extern int quitflag;

extern camera_t cam;

// box.c
float box_volume(box_t *b);
void box_depth_change_up(box_t *b);
void box_free(box_t *b);
void box_free_tree(box_t *b);
box_t *box_new(v4f_t *v0, v4f_t *v1, v4f_t *color, int op);
void box_print(box_t *b, int tabs);
box_t *box_inject(box_t *r, box_t *b);
int box_in(box_t *box, v4f_t *p);
int box_pair_touch(box_t *b0, box_t *b1);
box_t *box_in_tree(box_t *box, v4f_t *p, box_t **ignore, int ignore_count);
void box_normal(box_t *box, v4f_t *p, v4f_t *n, int inside);
float box_crosses(box_t *box, v4f_t *p, v4f_t *vi, int *inside);

// level.c
box_t *level_load_fp(FILE *fp);
box_t *level_load_fname(const char *fname);

// vec.c
void vec_norm(v4f_t *v);
void mat_ident(m4f_t *m);
void mat_mul(m4f_t *a, m4f_t *b);
uint32_t color_vec_sse(__m128 v);
void cam_rotate_by(float axz, float ayz, float axw, float ayw);

