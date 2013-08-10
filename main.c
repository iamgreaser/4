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

typedef union m4f
{
	struct { v4f_t x, y, z, w; } v;
	struct { v4f_t b, g, r, a; } c;
	struct { v4f_t s, t, r, q; } t;
	v4f_t a[4];
} m4f_t;

typedef struct camera
{
	m4f_t m;
	v4f_t o;
} camera_t;

#define MAX_BOX 10000

box_t *root = NULL;

SDL_Surface *screen = NULL;
uint32_t *rtbuf = NULL;
int rtbuf_width = 320;
int rtbuf_height = 240;
int rtbuf_scale = 3;
int quitflag = 0;

camera_t cam;

void vec_norm(v4f_t *v)
{
	// TODO: newton-raphson to improve accuracy?
	v4f_t rsp;
	rsp.m = _mm_mul_ps(v->m, v->m);
	float rsf = rsp.v.x + rsp.v.y + rsp.v.z + rsp.v.w;
	__m128 rs = _mm_rsqrt_ss(_mm_load_ss(&rsf));
	__m128 rsc = _mm_shuffle_ps(rs, rs, 0x00);
	v->m = _mm_mul_ps(v->m, rsc);
}

void mat_ident(m4f_t *m)
{
	m->v.x.m = _mm_set_ps(0, 0, 0, 1);
	m->v.y.m = _mm_set_ps(0, 0, 1, 0);
	m->v.z.m = _mm_set_ps(0, 1, 0, 0);
	m->v.w.m = _mm_set_ps(1, 0, 0, 0);
}

void mat_mul(m4f_t *a, m4f_t *b)
{
	// TODO: replace this with an SSE version
	int i,j,k;
	float v;

	m4f_t c;

	for(i=0; i<4; i++)
	for(j=0; j<4; j++)
	{
		v = 0.0f;

		for(k=0; k<4; k++)
			v += b->a[i].a[k] * a->a[k].a[j];

		c.a[i].a[j] = v;
	}

	for(i=0; i<4; i++)
	for(j=0; j<4; j++)
		a->a[i].a[j] = c.a[i].a[j];
}

uint32_t color_vec_sse(__m128 v)
{
	const __m128 mc_255 = _mm_set1_ps(255.0f);
	const __m128 mc_0 = _mm_set1_ps(0.0f);
	const __m128i mc_0i = _mm_set1_epi32(0);
	const __m128 mc_1 = _mm_set1_ps(1.0f);
	const __m128i mc_255i = _mm_set1_epi32(0x000000FF);

	// clamp
	v = _mm_min_ps(v, mc_1);
	v = _mm_max_ps(v, mc_0);

	// mul
	v = _mm_mul_ps(v, mc_255);

	// convert
	__m128i vi = _mm_cvtps_epi32(v);

	// pack
	vi = (__m128i)_mm_and_ps((__m128)vi, (__m128)mc_255i); // XXX is this step necessary?
	vi = _mm_packs_epi32(vi, mc_0i);
	vi = _mm_packus_epi16(vi, mc_0i);

	// return lowest dword
	return _mm_cvtsi128_si32(vi);
}

void cam_rotate_by(float axy, float axw, float ayw)
{
	m4f_t rotmat;
	float s, c;

	// aiming for left-hand coordinates - y up, x right, z forward, ana-kata wherever the fuck it goes

	// rotate around xw
	mat_ident(&rotmat);
	s = sinf(axw);
	c = cosf(axw);
	rotmat.v.y.m = _mm_set_ps(0, -s, c, 0);
	rotmat.v.z.m = _mm_set_ps(0, c, s, 0);
	mat_mul(&cam.m, &rotmat);

	// rotate around yw
	mat_ident(&rotmat);
	s = sinf(ayw);
	c = cosf(ayw);
	rotmat.v.x.m = _mm_set_ps(0, -s, 0, c);
	rotmat.v.z.m = _mm_set_ps(0, c, 0, s);
	mat_mul(&cam.m, &rotmat);

	// rotate around xy
	mat_ident(&rotmat);
	s = sinf(axy);
	c = cosf(axy);
	rotmat.v.z.m = _mm_set_ps(-s, c, 0, 0);
	rotmat.v.w.m = _mm_set_ps(c, s, 0, 0);
	mat_mul(&cam.m, &rotmat);
}

float box_volume(box_t *b)
{
	v4f_t v;

	v.m = _mm_sub_ps(b->v1.m, b->v0.m);

	return v.v.x * v.v.y * v.v.z * v.v.w;
}

void box_depth_change_up(box_t *b, int d)
{
	if(b == NULL)
		return;
	
	b->depth += d;
	return box_depth_change_up(b->p, d);
}

void box_depth_change_down(box_t *b, int d)
{
	if(b == NULL)
		return;
	
	b->depth += d;

	box_depth_change_down(b->c[0], d);
	return box_depth_change_down(b->c[1], d);
}

box_t *box_new(v4f_t *v0, v4f_t *v1, v4f_t *color, int op)
{
	box_t *b = malloc(sizeof(box_t));

	b->v0.m = v0->m;
	b->v1.m = v1->m;
	if(color != NULL)
		b->color.m = color->m;
	b->op = op;

	b->c[0] = b->c[1] = NULL;
	b->p = NULL;
	b->depth = 0;

	return b;
}

void box_print(box_t *b, int tabs)
{
	if(b == NULL)
		return;
	
	int i;

	for(i = 0; i < tabs; i++)
		printf(" ");

	printf("%i [%f %f %f %f] -> [%f %f %f %f]\n"
		, b->depth
		, b->v0.v.x
		, b->v0.v.y
		, b->v0.v.z
		, b->v0.v.w
		, b->v1.v.x
		, b->v1.v.y
		, b->v1.v.z
		, b->v1.v.w);
	
	tabs += 2;
	box_print(b->c[0], tabs);
	box_print(b->c[1], tabs);
}

box_t *box_inject(box_t *r, box_t *b)
{
	// if this is all there is then return our box
	if(r == NULL)
		return b;
	
	// is r a pair, as opposed to a leaf?
	if(r->op == SHP_PAIR)
	{
		// is b fully in r?
		v4f_t contv;
		contv.m = _mm_and_ps(
			_mm_cmpge_ps(b->v0.m, r->v0.m),
			_mm_cmpge_ps(r->v1.m, b->v1.m));
		
		if(_mm_movemask_ps(contv.m) == 15)
		{
			// contained. check best candidate for child.
			int ci = 0;

			if(r->c[0]->depth >= r->c[1]->depth)
				ci = 1;

			// apply.
			r->c[ci] = box_inject(r->c[ci], b);

			// return.
			return r;
		}
	}

	// create a parent node to encapsulate our two nodes
	v4f_t nv0, nv1;
	nv0.m = _mm_min_ps(b->v0.m, r->v0.m);
	nv1.m = _mm_max_ps(b->v1.m, r->v1.m);
	box_t *p = box_new(&nv0, &nv1, NULL, SHP_PAIR);

	// encapsulate!
	p->c[0] = r;
	p->c[1] = b;
	p->p = r->p;
	r->p = b->p = p;
	p->depth = p->c[p->c[0]->depth > p->c[1]->depth ? 0 : 1]->depth + 1;
	box_depth_change_up(p->p, 1);
	return p;
}

int box_in(box_t *box, v4f_t *p)
{
	v4f_t cmp;

	cmp.m = _mm_and_ps(
		_mm_cmpge_ps(p->m, box->v0.m),
		_mm_cmpge_ps(box->v1.m, p->m));
	
	return (_mm_movemask_ps(cmp.m) == 15);
}

void box_normal(box_t *box, v4f_t *p, v4f_t *n, int inside)
{
	// get centre
	__m128 cent = _mm_mul_ps(_mm_add_ps(box->v0.m, box->v1.m), _mm_set1_ps(0.5f));

	// get abs offset
	__m128 offs = _mm_sub_ps(p->m, cent);
	__m128 offs_abs = _mm_and_ps(offs, (__m128)_mm_set1_epi32(0x7FFFFFFF));

	// build mask
	__m128 mask = _mm_and_ps(
		_mm_and_ps(
			_mm_cmpge_ps(offs_abs, _mm_shuffle_ps(offs_abs, offs_abs, 0x4E)),
			_mm_cmpge_ps(offs_abs, _mm_shuffle_ps(offs_abs, offs_abs, 0x93))
		), 
		_mm_cmpge_ps(offs_abs, _mm_shuffle_ps(offs_abs, offs_abs, 0x39))
	);

	// get masked value and normalise
	if(inside)
		n->m = _mm_and_ps(mask, offs);
	else
		n->m = _mm_and_ps(mask, _mm_sub_ps(_mm_setzero_ps(), offs));

	vec_norm(n);
}

float box_crosses(box_t *box, v4f_t *p, v4f_t *vi, int *inside)
{
	// copy vi so we don't mutilate it
	v4f_t vz;
	v4f_t *v;
	vz.m = vi->m;
	v = &vz;

	// determine what faces are "in" or "out"
	v4f_t cmp0, cmp1;
	cmp0.m = _mm_cmpge_ps(p->m, box->v0.m);
	cmp1.m = _mm_cmpge_ps(box->v1.m, p->m);

	int outmask = (_mm_movemask_ps(cmp0.m) | (_mm_movemask_ps(cmp1.m) << 4));

	// check if we are inside this box (?)
	int start_inside = (outmask == 0xFF);
	
	// get deltas
	__m128 d0 = _mm_sub_ps(box->v0.m, p->m);
	__m128 d1 = _mm_sub_ps(box->v1.m, p->m);

	// mask them
	v4f_t m0, m1;
	m0.m = _mm_cmpgt_ps(_mm_mul_ps(v->m, d0), _mm_setzero_ps());
	m1.m = _mm_cmpgt_ps(_mm_mul_ps(v->m, d1), _mm_setzero_ps());

	int colmask = (_mm_movemask_ps(m0.m) | (_mm_movemask_ps(m1.m) << 4));

	// check if we can hit an in face
	if((!start_inside) && (colmask & ~outmask) == 0)
		return -1.0f;
	
	// get abs velocity
	v->m = _mm_max_ps(_mm_and_ps(v->m, (__m128)_mm_set1_epi32(0x7FFFFFFF)), _mm_set1_ps(0.00001f));

	// get abs deltas
	d0 = _mm_and_ps(d0, (__m128)_mm_set1_epi32(0x7FFFFFFF));
	d1 = _mm_and_ps(d1, (__m128)_mm_set1_epi32(0x7FFFFFFF));

	// get base time to plane
	__m128 tb0 = _mm_div_ps(d0, v->m);
	__m128 tb1 = _mm_div_ps(d1, v->m);

	// loop
	for(;;)
	{
		// get masked out time to plane
		// add a large value for stuff we can't touch
		v4f_t t0, t1;
		t0.m = _mm_add_ps(tb0, _mm_andnot_ps(m0.m, _mm_set1_ps(1000000000.0f)));
		t1.m = _mm_add_ps(tb1, _mm_andnot_ps(m1.m, _mm_set1_ps(1000000000.0f)));

		// find the first thing we hit
		int idx0, idx1;

		if(t0.v.x < t0.v.y && t0.v.x < t0.v.z && t0.v.x < t0.v.w)
			idx0 = F_XN;
		else if(t0.v.y < t0.v.z && t0.v.y < t0.v.w)
			idx0 = F_YN;
		else if(t0.v.z < t0.v.w)
			idx0 = F_ZN;
		else
			idx0 = F_WN;

		if(t1.v.x < t1.v.y && t1.v.x < t1.v.z && t1.v.x < t1.v.w)
			idx1 = F_XP;
		else if(t1.v.y < t1.v.z && t1.v.y < t1.v.w)
			idx1 = F_YP;
		else if(t1.v.z < t1.v.w)
			idx1 = F_ZP;
		else
			idx1 = F_WP;
		
		// get index + distance
		int idx;
		float minv;
		if(t0.a[idx0] <= t1.a[idx1&3])
		{
			idx = idx0;
			minv = t0.a[idx0];
		} else {
			idx = idx1;
			minv = t1.a[idx1&3];
		}

		// check if we started off inside this object - if so, give distance to plane NOW.
		if(start_inside)
		{
			if(inside != NULL) *inside = 1;
			return minv;
		}

		// check if this plane is an out plane - if so, FAIL.
		if(((1<<idx)&(~outmask)&colmask) == 0)
			return -1.0f;

		// remove this plane from colmask and add it to outmask
		colmask &= ~(1<<idx);
		outmask |= (1<<idx);

		// if outmask is 0xFF, we've FINALLY reached the inside of the volume - RETURN DISTANCE
		if((outmask&0xFF) == 0xFF)
		{
			if(inside != NULL) *inside = 0;
			return minv;
		}

		// mask this out
		if(idx < 4)
			m0.a[idx] = 0.0f;
		else
			m1.a[idx&3] = 0.0f;

		// loop around!
	}
}

float trace_box(box_t *box, v4f_t *p, v4f_t *v, v4f_t *color, box_t **retbox, int *inside, float md)
{
	// trace against this box
	v4f_t rb;
	rb.m = p->m;
	int ins;
	float d = (box->op == SHP_PAIR && box_in(box, p)
		? 0.0f
		: box_crosses(box, &rb, v, &ins));
	int doret = 0;

	if(d >= 0.0f && d <= md)
	{
		switch(box->op)
		{
			case SHP_PAIR: {
				// find nearest of the two
				v4f_t r0, r1;
				int ins0, ins1;
				r0.m = p->m;
				r1.m = p->m;

				box_t *retbox0 = NULL;
				float d0 = trace_box(box->c[0], &r0, v, color, &retbox0, &ins0, md);
				if(d0 > 0.0f && d0 <= md)
				{
					doret = 1;
					*retbox = retbox0;
					p->m = r0.m;
					if(inside != NULL) *inside = ins0;
					md = d0;
				}

				box_t *retbox1 = NULL;
				float d1 = trace_box(box->c[1], &r1, v, color, &retbox1, &ins1, md);
				if(d1 > 0.0f && d1 <= md)
				{
					doret = 1;
					*retbox = retbox1;
					p->m = r1.m;
					if(inside != NULL) *inside = ins1;
					md = d1;
				}
			} break;
			default:
				if(d > 0.0f && d <= md)
				{
					color->m = box->color.m;
					md = d;
					doret = 1;
					if(inside != NULL) *inside = ins;
					p->m = _mm_add_ps(p->m, _mm_mul_ps(
							v->m, _mm_set1_ps(d)));
					*retbox = box;
				}
				break;
		}
	}

	return doret ? md : -1.0f;
}

uint32_t trace_pixel(float sx, float sy)
{
	v4f_t f;

	f.m = _mm_add_ps(
		cam.m.v.z.m,
		_mm_add_ps(
			_mm_mul_ps(cam.m.v.x.m, _mm_set1_ps(sx)),
			_mm_mul_ps(cam.m.v.y.m, _mm_set1_ps(sy))
		));
	
	// normalise
	vec_norm(&f);

	// trace recursively
	v4f_t color;
	color.m = _mm_add_ps(_mm_mul_ps(f.m, _mm_set1_ps(0.5f)), _mm_set1_ps(0.5f));

	v4f_t p;
	box_t *box = NULL;
	p.m = cam.o.m;
	int inside;
	trace_box(root, &p, &f, &color, &box, &inside, 100.0f);

	if(box != NULL)
	{
		// calculate normal
		v4f_t n;
		box_normal(box, &p, &n, inside);

		// calculate diffuse (against one point for now)
		v4f_t l, dc;
		//printf("%f %f %f %f\n", p.v.x, p.v.y, p.v.z, p.v.w);
		l.m = cam.o.m;
		l.m = _mm_sub_ps(p.m, l.m);
		vec_norm(&l);
		//printf("%f %f %f %f\n", l.v.x, l.v.y, l.v.z, l.v.w);
		//printf("%f %f %f %f\n", n.v.x, n.v.y, n.v.z, n.v.w);
		dc.m = _mm_mul_ps(l.m, n.m);
		float diff = dc.v.x + dc.v.y + dc.v.z + dc.v.w;
		if(diff < 0.0f) diff = -diff;
		//printf("%f\n", diff);

		// multiply diffuse
		color.m = _mm_mul_ps(color.m, _mm_set1_ps(diff));
	}

	return color_vec_sse(color.m);
}

void render_screen(void)
{
	int x,y;
	uint32_t *p = rtbuf;
	SDL_LockSurface(screen);
	float sxd = (2.0f)/(rtbuf_width/2.0f);
	float syd = (2.0f)/(rtbuf_width/2.0f);
	float sxi = -sxd*rtbuf_width/2.0f;
	float syi = -syd*rtbuf_height/2.0f;

	float sy = syi;
	/*
	printf("%f %f %f %f\n", cam.m.v.x.v.x, cam.m.v.x.v.y, cam.m.v.x.v.z, cam.m.v.x.v.w);
	printf("%f %f %f %f\n", cam.m.v.y.v.x, cam.m.v.y.v.y, cam.m.v.y.v.z, cam.m.v.y.v.w);
	printf("%f %f %f %f\n", cam.m.v.z.v.x, cam.m.v.z.v.y, cam.m.v.z.v.z, cam.m.v.z.v.w);
	*/
	//printf("%f %f %f %f\n", cam.m.v.z.v.x, cam.m.v.z.v.y, cam.m.v.z.v.z, cam.m.v.z.v.w);
	for(y = 0; y < rtbuf_height; y++)
	{
		float sx = sxi;

		for(x = 0; x < rtbuf_width; x++)
		{
			*(p++) = trace_pixel(sx, sy);
			sx += sxd;
		}

		sy += syd;
	}

	// scale up onto screen
	uint32_t *d0, *d1, *d2;
	d0 = screen->pixels;
	d1 = (uint32_t *)((screen->pixels) + screen->pitch);
	d2 = (uint32_t *)((screen->pixels) + screen->pitch*2);
	int pitch = screen->pitch*3/4 - screen->w;

	SDL_LockSurface(screen);

	p = rtbuf;
	for(y = 0; y < rtbuf_height; y++)
	{
		for(x = 0; x < rtbuf_width; x++)
		{
			uint32_t v = *(p++);
			*(d0++) = v;
			*(d0++) = v;
			*(d0++) = v;
			*(d1++) = v;
			*(d1++) = v;
			*(d1++) = v;
			*(d2++) = v;
			*(d2++) = v;
			*(d2++) = v;
		}

		d0 += pitch;
		d1 += pitch;
		d2 += pitch;
	}

	SDL_UnlockSurface(screen);

	SDL_Flip(screen);
}

void cam_init(void)
{
	cam.o.m = _mm_set_ps(0, 0, 0, 0);
	mat_ident(&cam.m);
}

void level_init(void)
{
	v4f_t v0, v1, color;
	box_t **blist = malloc(sizeof(box_t *) * 100000);

	int i;
	int bcount = 0;

	for(i = 0; i < 8; i++)
	{
		v0.m = _mm_set_ps(-1.0f, -1.0f, -1.0f, -1.0f);
		v1.m = _mm_set_ps(1.0f, 1.0f, 1.0f, 1.0f);

		float offs_f = i < 4 ? 3.0f : -3.0f;
		__m128 offs = _mm_load_ss(&offs_f);
		if(i&1) offs = _mm_shuffle_ps(offs, offs, 0x93);
		if(i&2) offs = _mm_shuffle_ps(offs, offs, 0x4E);

		v0.m = _mm_add_ps(v0.m, offs);
		v1.m = _mm_add_ps(v1.m, offs);

		color.m = _mm_set_ps(1.0f, 1.0f, 1.0f, 1.0f);

		blist[i] = box_new(&v0, &v1, &color, SHP_ADD);
		bcount++;
	}

	for(i = 0; i < bcount-1; i += 2, bcount++)
		blist[bcount] = box_inject(blist[i], blist[i+1]);
	
	root = blist[bcount-1];
	free(blist);
	box_print(root, 0);
}

void render_main(void)
{
	SDL_Event ev;

	quitflag = 0;
	float vaxy = 0.0f;
	float vaxw = 0.0f;
	float vayw = 0.0f;
	float vx = 0.0f;
	float vy = 0.0f;
	float vz = 0.0f;
	float vw = 0.0f;
	while(!quitflag)
	{
		render_screen();

		const float vas = 0.03f;
		cam_rotate_by(vaxy*vas, vaxw*vas, vayw*vas);

		const float vs = 0.1f;
		cam.o.m = _mm_add_ps(cam.o.m, _mm_mul_ps(cam.m.v.x.m, _mm_set1_ps(vx*vs)));
		cam.o.m = _mm_add_ps(cam.o.m, _mm_mul_ps(cam.m.v.y.m, _mm_set1_ps(vy*vs)));
		cam.o.m = _mm_add_ps(cam.o.m, _mm_mul_ps(cam.m.v.z.m, _mm_set1_ps(vz*vs)));
		cam.o.m = _mm_add_ps(cam.o.m, _mm_mul_ps(cam.m.v.w.m, _mm_set1_ps(vw*vs)));

		while(SDL_PollEvent(&ev))
		switch(ev.type)
		{
			case SDL_QUIT:
				quitflag = 1;
				break;

			case SDL_KEYUP:
			switch(ev.key.keysym.sym)
			{
				case SDLK_w:
				case SDLK_s:
					vz = 0.0f;
					break;
				case SDLK_a:
				case SDLK_d:
					vx = 0.0f;
					break;
				case SDLK_q:
				case SDLK_e:
					vw = 0.0f;
					break;
				case SDLK_LCTRL:
				case SDLK_SPACE:
					vy = 0.0f;
					break;
					
				case SDLK_u:
				case SDLK_o:
					vaxy = 0.0f;
					break;
				case SDLK_i:
				case SDLK_k:
					vaxw = 0.0f;
					break;
				case SDLK_j:
				case SDLK_l:
					vayw = 0.0f;
					break;
				default:
					break;
			} break;

			case SDL_KEYDOWN:
			switch(ev.key.keysym.sym)
			{
				case SDLK_s:
					vz = -1.0f;
					break;
				case SDLK_a:
					vx = -1.0f;
					break;
				case SDLK_q:
					vw = -1.0f;
					break;
				case SDLK_SPACE:
					vy = -1.0f;
					break;
				case SDLK_w:
					vz = +1.0f;
					break;
				case SDLK_d:
					vx = +1.0f;
					break;
				case SDLK_e:
					vw = +1.0f;
					break;
				case SDLK_LCTRL:
					vy = +1.0f;
					break;

				case SDLK_u:
					vaxy = -1.0f;
					break;
				case SDLK_i:
					vaxw = -1.0f;
					break;
				case SDLK_j:
					vayw = -1.0f;
					break;
				case SDLK_o:
					vaxy = +1.0f;
					break;
				case SDLK_k:
					vaxw = +1.0f;
					break;
				case SDLK_l:
					vayw = +1.0f;
					break;
				default:
					break;
			} break;
		}
	}
}

int main(int argc, char *argv[])
{
	SDL_Init(SDL_INIT_TIMER | SDL_INIT_AUDIO | SDL_INIT_VIDEO);
	SDL_WM_SetCaption("4 (pronounced \"sardinelauncher\")", NULL);
	screen = SDL_SetVideoMode(rtbuf_width * rtbuf_scale, rtbuf_height * rtbuf_scale, 32, 0);
	rtbuf = malloc(rtbuf_width * rtbuf_height * 4);

	cam_init();
	level_init();
	render_main();

	free(rtbuf);

	return 0;
}

