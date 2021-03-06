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

#include <enet/enet.h>

#include <SDL.h>

#ifdef WIN32
#ifdef stderr
#undef stderr
#endif
#define stderr stdout
#endif

#ifndef __SSE2__
#error "SSE2 is required to compile + run this. Use the -msse2 flag in GCC."
#endif


enum
{
	SHP_SCENE_LIST,

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
};

typedef struct sphere
{
	v4f_t v;
	v4f_t color;
	float r;
	int plidx;
} sphere_t;

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

typedef struct kd kd_t;
struct kd
{
	// outer box (for SHP_SUB hypervolumes)
	// only useful for the fringes of a space
	box_t *obox[2];

	// associated box (i.e. what fills this space)
	box_t *box;

	// children
	//
	// CONSTRAINT: If c[0] == NULL, then c[1] == NULL.
	//
	// CONSTRAINT: If c[0] != NULL, then box == NULL.
	//
	// When c[1] != NULL, this is a pair.
	// When c[0] != NULL, but c[1] == NULL, this is an axis change.
	// When c[0] == NULL, this is a leaf.
	//
	// The following two constraints apply where box == NULL and c[1] != NULL.
	//
	// CONSTRAINT: c[0]->p == c[1]->p == this
	//
	// CONSTRAINT:
	// (c[0]->axis == c[1]->axis == this->axis
	//  && c[0]->v[1] == c[1]->v[0]
	//    && c[0]->v[0] == this->v[0]
	//    && c[1]->v[1] == this->v[1])
	//
	// The following two constraints apply where box == NULL, but c[1] == NULL.
	//
	// CONSTRAINT: c[0]->p == this
	//
	// CONSTRAINT: c[0]->axis == this->axis + 1
	//
	kd_t *c[2];

	// parent
	kd_t *p;

	// start/end marks
	//
	// CONSTRAINT: v[0] < v[1]
	float v[2];

	// current axis where 0,1,2,3 = X,Y,Z,W
	//
	// CONSTRAINT: parent->axis <= this->axis
	int axis;

	//
	// ACCELERATION VARIABLES
	// note, any of these can be NULL during tree construction.
	//

	// axis parent
	//
	// CONSTRAINT: if axis == 0, ap == NULL
	//
	// CONSTRAINT: if axis > 0,
	//   ap->axis == this->axis - 1
	//   && ap is an axis change node
	kd_t *ap;

	// adjacency list
	kd_t *adj[8];
};

typedef struct player
{
	camera_t cam;
	v4f_t lv;
	char name[33];
	uint8_t grounded;
	uint8_t vflip;
	char pad1[1];
	float grav_v;
	char pad2[5];
	uint8_t pid;
	int8_t team;
	uint8_t magic; // 0xC4, or 0xC9 if this is the current player, or 0x4C if not in game, or 0x66 if ingame and dead, or 0x69 if also current player
} __attribute__((__packed__)) player_t;

// don't reach this unless you're at a lan with a gigabit connection to a switch with a 100 gigabit interconnect
#define PLAYERS_MAX 128

extern int fps_counter;
extern int fps_next_tick;
extern box_t *root;
extern sphere_t *sroot;
extern int sroot_len;

extern int net_mode;
extern ENetHost *h_server;
extern ENetHost *h_client;
extern ENetPeer *h_toserver;
extern ENetPeer *h_player[PLAYERS_MAX];

extern ENetAddress a_server;
extern ENetAddress a_client;

extern SDL_Surface *screen;
extern uint32_t *rtbuf;
extern int rtbuf_width;
extern int rtbuf_height;
extern int rtbuf_scale;
extern int quitflag;

extern int bseed;
extern int mbutts;
extern int mrelease;

extern player_t players[];
extern int cplr;

// box.c
float box_volume(box_t *b);
void box_depth_change_up(box_t *b);
void box_free(box_t *b);
void box_free_tree(box_t *b);
box_t *box_new(v4f_t *v0, v4f_t *v1, v4f_t *color, int op);
void box_print(box_t *b, int tabs);
box_t *box_inject(box_t *r, box_t *b);
int box_in(const box_t *box, const v4f_t *p);
int box_pair_touch(box_t *b0, box_t *b1);
box_t *box_in_tree(box_t *box, const v4f_t *p, box_t **ignore, int ignore_count);
void box_normal(box_t *box, v4f_t *p, v4f_t *n, int inside);
float box_crosses(const box_t *box, const v4f_t *p, const v4f_t *vi, int *inside, int *side);

// game.c
void game_player_tick(player_t *pl, float dt);
int game_input(player_t *pl, float dt);

// kd.c
void kd_free(kd_t *kd);
void kd_free_down(kd_t *kd);
void kd_print(const kd_t *kd, int tabs);
const kd_t *kd_in_tree(const kd_t *kd, const v4f_t *p, v4f_t *b0, v4f_t *b1);
kd_t *kd_add_box(kd_t *kd, box_t *box);
void kd_accelerate(kd_t *kd);

// level.c
box_t *level_load_fp(FILE *fp);
box_t *level_load_fname(const char *fname);

// net.c
char *net_new_client(const char *shost, const char *sport, int toourselves);
void net_new_server(const char *sport, const char *fname);
void net_update_client(void);
void net_update_server(void);

// render.c
void render_main(void);

// sphere.c
sphere_t *sphere_list_add(sphere_t *l, int *llen, const v4f_t *v, float r, const v4f_t *color, int plidx);
void sphere_list_rm(sphere_t *l, int *llen, int idx);
void sphere_normal(const sphere_t *s, const v4f_t *p, v4f_t *normal);
float sphere_trace(sphere_t *l, const int *llen, const v4f_t *p, const v4f_t *v, float md, sphere_t **rets);

// vec.c
void vec_norm(v4f_t *v);
void mat_ident(m4f_t *m);
void mat_mul(m4f_t *a, m4f_t *b);
uint32_t color_vec_sse(__m128 v);
void cam_rotate_by(camera_t *cam, float axz, float ayz, float axw, float ayw, float axy);

// main.c
void *malloc16(size_t len);
void *realloc16(void *p, size_t len);
void free16(void *p);
void refresh_fps(void);
float trace_box(box_t *box, const v4f_t *p, const v4f_t *v, v4f_t *color, box_t **retbox, int *inside, float md, int *side);
uint32_t trace_pixel(camera_t *cam, box_t *bstart, float sx, float sy, const v4f_t *dirx, const v4f_t *diry, const v4f_t *dirz, const v4f_t *dirw, int *seed);
void player_sphere_clear(void);
void player_sphere_add(player_t *pl);

