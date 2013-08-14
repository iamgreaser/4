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

#include "common.h"

#define MAX_BOX 10000

int bseed = 1;

// util functions
void *align16(void *p, size_t len)
{
	int offs = sizeof(size_t) + sizeof(int);
	char *q = p + offs;
	while((((long)q) & 15) != 0)
	{
		q++;
		offs++;
	}

	int *r = (int *)q;
	r[-1] = offs-sizeof(int);

	size_t *s = (size_t *)(&r[-1]);
	s--;
	*s = len;

	return (void *)r;
}

void *unalign16(void *p)
{
	int *q = p;
	q--;
	char *r = (char *)q;
	r -= *q;
	return (void *)r;
}

size_t getsize16(void *p)
{
	char *q = p;
	q -= sizeof(int) + sizeof(size_t);
	return *(size_t *)q;
}

void *malloc16(size_t len)
{
	return align16(malloc(len + 16 + 2*sizeof(int)), len);
}

void free16(void *p)
{
	free(unalign16(p));
}

void *realloc16(void *p, size_t len)
{
	// naive method 1: align16(realloc(unalign16(p), len), len) - DOES NOT WORK.
	// naive method 2 should.

	if(len <= 0)
	{
		if(p != NULL)
			free16(p);

		return NULL;
	}

	// allocate a little bit more.
	// just so, y'know, we can be a bit lazy in other places.
	size_t mlen = len*3/2+1;
	if(mlen < len) mlen = len;

	if(p == NULL)
		return malloc16(mlen);
	
	// compare lengths; if <= current, don't expand.
	size_t clen = getsize16(p);
	if(clen >= len)
		return p;
	
	// malloc, memcpy, free.
	void *np = malloc16(mlen);
	memcpy(np, p, clen);
	free16(p);

	return np;
}

int myrand(int *seed)
{
	int v = (*seed * 1103515245) + 12345;
	v &= 0x7FFFFFFF;
	*seed = v;
	return v;
}

__m128 myrand_sse_mul(int *seed, float mul)
{
	__m128 v = _mm_set_ps(
		(float)(myrand(seed)%65539),
		(float)(myrand(seed)%65539),
		(float)(myrand(seed)%65539),
		(float)(myrand(seed)%65539));
	
	v = _mm_sub_ps(v, _mm_set1_ps(65539.0f / 2.0f));
	v = _mm_mul_ps(v, _mm_set1_ps(mul * 2.0f / 65539.0f));
	return v;
}

// TODO: move these two to the camera (well, a player entity, that is)
float grav_v = 0.0f;
int grounded = 0;

int mbutts = 0;
int mrelease = 0;

int fps_counter = 0;
int fps_next_tick = 0;
box_t *root = NULL;
kd_t *kdroot = NULL;

sphere_t *sroot = NULL;
int sroot_len = 0;

SDL_Surface *screen = NULL;
int rtbuf_width = 320;
int rtbuf_height = 240;
int rtbuf_scale = 3;
int quitflag = 0;

player_t players[PLAYERS_MAX];
int cplr = 0;

void refresh_fps(void)
{
	char buf[128];
	sprintf(buf, "4 (pronounced \"sardinelauncher\") - FPS: %i", fps_counter);
	SDL_WM_SetCaption(buf, NULL);
	fps_counter = 0;
	fps_next_tick += 1000;
}

float trace_into_box(box_t **retbox, const v4f_t *p, const v4f_t *v, float md, int *side)
{
	box_t *box = *retbox;
	float d = 0.0f;

	// check if we're in this box.
	// (skipping SHP_PAIR as box_crosses isn't the fastest op ever)
	if(box->op != SHP_PAIR && !box_in(box, p))
	{
		// nope. can we trace to it?
		//d = box_crosses_outside(box, p, v);
		d = box_crosses(box, p, v, NULL, side);

		if(d < -1.0f || d >= md)
		{
			// nope.
			return -1.0f;
		}
	}

	// check if pair
	if(box->op == SHP_PAIR)
	{
		// yes. trace children.
		box_t *b0 = box->c[0];
		box_t *b1 = box->c[1];
		float d0 = trace_into_box(&b0, p, v, md, side);
		if(d0 >= 0.0f)
			md = d0;
		float d1 = trace_into_box(&b1, p, v, md, side);

		if(d1 >= 0.0f && d1 < d0)
		{
			*retbox = b1;
			return d1;
		} else {
			if(d0 >= 0.0f) *retbox = b0;
			return d0;
		}
	} else {
		// no. return distance.
		*retbox = box;
		return d;
	}
}

#define MAX_IGN 100
float trace_box(box_t *box, const v4f_t *p, const v4f_t *v, v4f_t *color, box_t **retbox, int *inside, float md, int *side)
{
	int i; (void)i;
	box_t *ign_l[MAX_IGN];
	int ign_c = 0;
	int ins;
	v4f_t tp;
	float ad = 0.0f;
	float td;

	// just in case we have an empty level!
	if(box == NULL)
		return -1.0f;
	
	// find the box we're tracing from within
	box_t *obox = box;
	box = box_in_tree(box, p, ign_l, ign_c);

	// store trace position
	tp.m = p->m;

	// are we in anything?
	if(box == NULL)
	{
		// trace inwards.
		box = obox;
		td = trace_into_box(&box, p, v, md, side);

		// did we go anywhere?
		if(td < 0.0f || td >= md)
		{
			// nope.
			return -1.0f;
		} else {
			// yep.
			// return the point.
			if(retbox != NULL) *retbox = box;
			if(inside != NULL) *inside = 0;
			if(color != NULL) color->m = box->color.m;
			return td;
		}
	}

	// trace onwards.
	for(;;)
	{
		// trace to the box end.
		td = box_crosses(box, &tp, v, &ins, side);

		// ASSERTION.
		if(td < 0.0f)
			// FIXME: this sometimes happens.
			return -1.0f;
			//abort();

		// check if we have exceeded our distance.
		ad += td;
		if(ad >= md)
			// we have. fail.
			return -1.0f;

		// moving right along...
		tp.m = _mm_add_ps(p->m, _mm_mul_ps(_mm_set1_ps(ad), v->m));
		//tp.m = _mm_add_ps(tp.m, _mm_mul_ps(_mm_set1_ps(td), v->m));
		
#if 0
		// remove old entries from ignore list
		// XXX: this actually slows it down. let's NOT do this.
		for(i = 0; i < ign_c; i++)
			if(!box_in(ign_l[i], &tp))
				ign_l[i--] = ign_l[ign_c--];
#endif

		// add our current box to the ignore list
		ign_l[ign_c++] = box;
		if(ign_c > 100)
			// FIXME: this sometimes happens.
			return -1.0f;
			//abort();

		// check to see if we're in another box.
		box_t *nbox = box_in_tree(box, &tp, ign_l, ign_c);
		if(nbox == NULL || nbox->op != SHP_SUB)
		{
			// nope. let's end it here.
			if(nbox == NULL)
				nbox = box;

			if(retbox != NULL) *retbox = nbox;
			if(inside != NULL) *inside = ins;
			if(color != NULL) color->m = nbox->color.m;
			return ad;
		}

		box = nbox;
	}
}

void trace_main(int bouncerem, box_t *bstart, box_t *bign, const v4f_t *p, const v4f_t *f, const v4f_t *cf, v4f_t *color, float md, int *seed)
{
	int side;
	int inside;
	box_t *box = bign;
	sphere_t *s = NULL;

	if(bouncerem-- <= 0)
		return;

	v4f_t rcolor;
	rcolor.m = color->m;
	float d1 = trace_box(bstart, p, f, &rcolor, &box, &inside, md, &side);
	float d2 = sphere_trace(sroot, &sroot_len, p, f, (d1 > 0.0f ? d1 : md), &s);

	if(d1 >= 0.0f || d2 > 0.00001f)
	{
		int use_sphere = (d2 > 0.00001f && (d1 < 0.0f || d2 < d1));
		float d = (use_sphere ? d2 : d1);

		// move position
		v4f_t np;
		np.m = _mm_add_ps(p->m, _mm_mul_ps(
			f->m, _mm_set1_ps(d-0.01f)));

		// calculate normal
		v4f_t n;
		if(use_sphere)
		{
			rcolor.m = s->color.m;
			n.m = _mm_sub_ps(s->v.m, np.m);
			vec_norm(&n);
		} else {
			n.m = _mm_setzero_ps();
			n.a[side&3] = ((side & 4) ? 1.0f : -1.0f) * (inside ? 1.0f : -1.0f);
		}

		// calculate diffuse (against one point for now)
		v4f_t dc;
		dc.m = _mm_mul_ps(f->m, n.m);
		float diff = dc.v.x + dc.v.y + dc.v.z + dc.v.w;

		// bounce
		v4f_t ncolor, fb, fb2;
		fb.m = f->m;
		fb.a[side&3] *= -1.0f;
		ncolor.m = _mm_setzero_ps();
		const float blurfac = 0.02f;
		fb2.m = _mm_add_ps(fb.m, myrand_sse_mul(seed, blurfac));
		trace_main(bouncerem, box, NULL, &np, &fb2, cf, &ncolor, md - d, seed);

		// multiply diffuse
		const float refl = 0.3f;
		color->m = _mm_add_ps(
			_mm_mul_ps(ncolor.m, _mm_set1_ps(refl)),
			_mm_mul_ps(rcolor.m, _mm_set1_ps(diff*(1.0f-refl)))
		);

	}
}

__m128 trace_pixel_step(camera_t *cam, box_t *bstart, float sx, float sy, const v4f_t *dirx, const v4f_t *diry, const v4f_t *dirz, const v4f_t *dirw, int *seed)
{
	v4f_t p, f;

	const float bluraper = 0.00f;
	const float blurdist = 1.5f;

	if(bluraper == 0.0f)
	{
		p.m = cam->o.m;
		f.m = _mm_add_ps(
			dirz->m,
			_mm_add_ps(
				_mm_mul_ps(dirx->m, _mm_set1_ps(sx)),
				_mm_mul_ps(diry->m, _mm_set1_ps(sy))
			));
		
		// normalise
		vec_norm(&f);
	} else {
		// give a slight offset to cause a focal blur
		__m128 bluroffs = myrand_sse_mul(seed, bluraper);
		__m128 focus = _mm_add_ps(p.m, _mm_mul_ps(_mm_set1_ps(blurdist), dirz->m));
		p.m = _mm_add_ps(cam->o.m, bluroffs);
		f.m = _mm_sub_ps(focus, p.m);
		vec_norm(&f);

		f.m = _mm_add_ps(
			f.m,
			_mm_add_ps(
				_mm_mul_ps(dirx->m, _mm_set1_ps(sx)),
				_mm_mul_ps(diry->m, _mm_set1_ps(sy))
			));
		vec_norm(&f);
	}
	
	// disabled because it spoils the suprise
	/*
	// slightly distort for w
	float wdist = 0.1f;
	f.m = _mm_add_ps(
		f.m,
		_mm_mul_ps(dirw->m, _mm_set1_ps(sqrtf(sx*sx + sy*sy)*wdist))
	);
	*/

	// trace recursively
	v4f_t color;
	color.m = _mm_add_ps(_mm_mul_ps(f.m, _mm_set1_ps(0.5f)), _mm_set1_ps(0.5f));

	trace_main(3, bstart, NULL, &p, &f, &f, &color, 100.0f, seed);

	return color.m;
}

uint32_t trace_pixel(camera_t *cam, box_t *bstart, float sx, float sy, const v4f_t *dirx, const v4f_t *diry, const v4f_t *dirz, const v4f_t *dirw, int *seed)
{
	__m128 color;
	__m128 cmin, cmax;
	v4f_t cdiff;
	int i;

	const int bcmin = 1;
	const int bcmax = 1;
	const float bcap = 1.0f;
	/*
	const int bcmin = 3;
	const int bcmax = 10;
	const float bcap = 0.5f;
	*/

	const int bcvary = bcmax - bcmin;

	int bc = 1;

	color = trace_pixel_step(cam, bstart, sx, sy, dirx, diry, dirz, dirw, seed);
	cmin = cmax = color;

	for(i = 0; i < bcmin-1; i++, bc++)
	{
		__m128 ncolor = trace_pixel_step(cam, bstart, sx, sy, dirx, diry, dirz, dirw, seed);
		cmin = _mm_min_ps(cmin, ncolor);
		cmax = _mm_max_ps(cmax, ncolor);
		color = _mm_add_ps(color, ncolor);
	}

	cdiff.m = _mm_sub_ps(cmax, cmin);

	float bdiff = (cdiff.a[0]
		+ cdiff.a[1]
		+ cdiff.a[2]);
	
	bdiff /= bcap;
	if(bdiff > 1.0f)
		bdiff = 1.0f;
	
	bdiff = 1.0f - bdiff;
	bdiff = bdiff * bdiff;
	bdiff = 1.0f - bdiff;
	int bcvc = (int)(bcvary*bdiff + 0.9f);
	for(i = 0; i < bcvc; i++, bc++)
		color = _mm_add_ps(color, trace_pixel_step(cam, bstart, sx, sy, dirx, diry, dirz, dirw, seed));
	
	color = _mm_mul_ps(color, _mm_set1_ps(1.0f/bc));
	return color_vec_sse(color);
}

void cam_init(void)
{
	int i;
	for(i = 0; i < PLAYERS_MAX; i++)
	{
		camera_t *cam = &(players[i].cam);
		cam->o.m = _mm_set_ps(0, 0, 0, 0);
		mat_ident(&cam->m);
	}
}

void level_init(const char *fname)
{
	root = level_load_fname(fname);
	//kdroot = kd_add_box(kdroot, root);
	//kd_accelerate(kdroot);
	//kd_print(kdroot, 0);

	v4f_t v, c;
	v.m = _mm_set_ps(0,1,0,0); c.m = _mm_set_ps(1,1,0,1); sroot = sphere_list_add(sroot, &sroot_len, &v, 0.5f, &c);
	v.m = _mm_set_ps(0,1.7f,0,0); c.m = _mm_set_ps(1,1,0,1); sroot = sphere_list_add(sroot, &sroot_len, &v, 0.3f, &c);
	v.m = _mm_set_ps(0.5f,1.5f,0,0); c.m = _mm_set_ps(1,1,0,1); sroot = sphere_list_add(sroot, &sroot_len, &v, 0.2f, &c);
	v.m = _mm_set_ps(0,1.5f,0,0.5f); c.m = _mm_set_ps(1,1,0,1); sroot = sphere_list_add(sroot, &sroot_len, &v, 0.2f, &c);
	v.m = _mm_set_ps(-0.5f,1.5f,0,0); c.m = _mm_set_ps(1,1,0,1); sroot = sphere_list_add(sroot, &sroot_len, &v, 0.2f, &c);
	v.m = _mm_set_ps(0,1.5f,0,-0.5f); c.m = _mm_set_ps(1,1,0,1); sroot = sphere_list_add(sroot, &sroot_len, &v, 0.2f, &c);
}

int main(int argc, char *argv[])
{
	const char *fname;
	if(argc <= 1)
	{
		printf("%s dat/shortcut.4lv (or something like that).\n", argv[0]);
		fname = "dat/shortcut.4lv";
	} else {
		fname = argv[1];
	}

	SDL_Init(SDL_INIT_TIMER | SDL_INIT_AUDIO | SDL_INIT_VIDEO);
	SDL_WM_SetCaption("4 (pronounced \"sardinelauncher\")", NULL);
	screen = SDL_SetVideoMode(rtbuf_width * rtbuf_scale, rtbuf_height * rtbuf_scale, 32, 0);

	cam_init();
	level_init(fname);

	SDL_WM_GrabInput(1);
	SDL_ShowCursor(0);
	render_main();
	SDL_ShowCursor(1);
	SDL_WM_GrabInput(0);

	return 0;
}

