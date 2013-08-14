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
	char *q = p + sizeof(size_t) + sizeof(int);
	int offs = 0;
	while((((long)q) & 15) != 0)
	{
		q++;
		offs++;
	}

	int *r = (int *)q;
	r--;
	*r = offs-4;

	size_t *s = (size_t *)r;
	s--;
	*s = len;

	r++;
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

SDL_Surface *screen = NULL;
int rtbuf_width = 320;
int rtbuf_height = 240;
int rtbuf_scale = 3;
int quitflag = 0;

camera_t cam;

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

	if(bouncerem-- <= 0)
		return;

	v4f_t rcolor;
	rcolor.m = color->m;
	float d = trace_box(bstart, p, f, &rcolor, &box, &inside, md, &side);

	if(d >= 0.0f)
	{
		// move position
		v4f_t np;
		np.m = _mm_add_ps(p->m, _mm_mul_ps(
			f->m, _mm_set1_ps(d-0.01f)));

		// calculate normal
		v4f_t n;
		n.m = _mm_setzero_ps();
		n.a[side&3] = ((side & 4) ? 1.0f : -1.0f) * (inside ? 1.0f : -1.0f);

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

__m128 trace_pixel_step(box_t *bstart, float sx, float sy, const v4f_t *dirx, const v4f_t *diry, const v4f_t *dirz, const v4f_t *dirw, int *seed)
{
	v4f_t p, f;

	const float bluraper = 0.00f;
	const float blurdist = 1.5f;

	if(bluraper == 0.0f)
	{
		p.m = cam.o.m;
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
		p.m = _mm_add_ps(cam.o.m, bluroffs);
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

uint32_t trace_pixel(box_t *bstart, float sx, float sy, const v4f_t *dirx, const v4f_t *diry, const v4f_t *dirz, const v4f_t *dirw, int *seed)
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

	color = trace_pixel_step(bstart, sx, sy, dirx, diry, dirz, dirw, seed);
	cmin = cmax = color;

	for(i = 0; i < bcmin-1; i++, bc++)
	{
		__m128 ncolor = trace_pixel_step(bstart, sx, sy, dirx, diry, dirz, dirw, seed);
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
		color = _mm_add_ps(color, trace_pixel_step(bstart, sx, sy, dirx, diry, dirz, dirw, seed));
	
	color = _mm_mul_ps(color, _mm_set1_ps(1.0f/bc));
	return color_vec_sse(color);
}

void render_viewport(int vx, int vy, int w, int h, v4f_t *dx, v4f_t *dy, v4f_t *dz, v4f_t *dw)
{
	int y;

	float sxd = (2.0f)/(w/2.0f);
	float syd = (2.0f)/(w/2.0f);
	float sxi = -sxd*w/2.0f;
	float syi = -syd*h/2.0f;

	// find starting box
	box_t *bstart = box_in_tree(root, &cam.o, NULL, 0);
	if(bstart == NULL) bstart = root;

	// scale up onto screen
	#pragma omp parallel for
	for(y = 0; y < h; y++)
	{
		int x;
		uint32_t *d0, *d1, *d2;
		__attribute__((aligned(16))) uint32_t vs[4];

		int seed = bseed ^ (y*6325);

		d0 = (uint32_t *)((screen->pixels) + screen->pitch*(0+3*(y+vy)) + 4*3*vx);
		d1 = (uint32_t *)((screen->pixels) + screen->pitch*(1+3*(y+vy)) + 4*3*vx);
		d2 = (uint32_t *)((screen->pixels) + screen->pitch*(2+3*(y+vy)) + 4*3*vx);

		float sy = syi + syd*y;

		float sx = sxi;

		for(x = 0; x < w-3 && (((long)d0)&15) != 0; x++)
		{
			uint32_t v = trace_pixel(bstart, sx, sy, dx, dy, dz, dw, &seed);

			*(d0++) = v;
			*(d0++) = v;
			*(d0++) = v;
			*(d1++) = v;
			*(d1++) = v;
			*(d1++) = v;
			*(d2++) = v;
			*(d2++) = v;
			*(d2++) = v;

			sx += sxd;
		}

		for(; x < w-3; x += 4)
		{
			vs[0] = trace_pixel(bstart, sx, sy, dx, dy, dz, dw, &seed); sx += sxd;
			vs[1] = trace_pixel(bstart, sx, sy, dx, dy, dz, dw, &seed); sx += sxd;
			vs[2] = trace_pixel(bstart, sx, sy, dx, dy, dz, dw, &seed); sx += sxd;
			vs[3] = trace_pixel(bstart, sx, sy, dx, dy, dz, dw, &seed); sx += sxd;

			__m128i pb = _mm_load_si128((__m128i *)vs);

			// stream!
			__m128i p0 = _mm_shuffle_epi32(pb, 0x40);
			__m128i p1 = _mm_shuffle_epi32(pb, 0xA5);
			_mm_store_si128((__m128i *)d0, p0);
			_mm_store_si128((__m128i *)d1, p0);
			_mm_store_si128((__m128i *)d2, p0);
			d0 += 4; d1 += 4; d2 += 4;

			__m128i p2 = _mm_shuffle_epi32(pb, 0xFE);
			_mm_store_si128((__m128i *)d0, p1);
			_mm_store_si128((__m128i *)d1, p1);
			_mm_store_si128((__m128i *)d2, p1);
			d0 += 4; d1 += 4; d2 += 4;

			_mm_store_si128((__m128i *)d0, p2);
			_mm_store_si128((__m128i *)d1, p2);
			_mm_store_si128((__m128i *)d2, p2);
			d0 += 4; d1 += 4; d2 += 4;
		}

		for(; x < w; x++)
		{
			uint32_t v = trace_pixel(bstart, sx, sy, dx, dy, dz, dw, &seed);

			*(d0++) = v;
			*(d0++) = v;
			*(d0++) = v;
			*(d1++) = v;
			*(d1++) = v;
			*(d1++) = v;
			*(d2++) = v;
			*(d2++) = v;
			*(d2++) = v;

			sx += sxd;
		}

		sy += syd;
	}
}

void cam_init(void)
{
	cam.o.m = _mm_set_ps(0, 0, 0, 0);
	mat_ident(&cam.m);
}

void level_init(const char *fname)
{
	root = level_load_fname(fname);
	//kdroot = kd_add_box(kdroot, root);
	//kd_accelerate(kdroot);
	//kd_print(kdroot, 0);
}

void render_screen(void)
{
	int y0 = 0;
	int y1 = rtbuf_height*3/4;
	int y2 = rtbuf_height;

	int x0 = 0;
	int x1 = rtbuf_width/3;
	int x2 = rtbuf_width*2/3;
	int x3 = rtbuf_width;

	v4f_t nx,  nz, nw;
	nx.m = _mm_mul_ps(cam.m.v.x.m, _mm_set1_ps(-1.0f));
	//ny.m = _mm_mul_ps(cam.m.v.y.m, _mm_set1_ps(-1.0f));
	nz.m = _mm_mul_ps(cam.m.v.z.m, _mm_set1_ps(-1.0f));
	nw.m = _mm_mul_ps(cam.m.v.w.m, _mm_set1_ps(-1.0f));

	SDL_LockSurface(screen);

	if(mbutts & 4)
	{
		render_viewport(x0, y0, x3-x0, y1-y0, &nx, &cam.m.v.y, &cam.m.v.w, &cam.m.v.z);
		render_viewport(x0, y1, x1-x0, y2-y1, &nx, &cam.m.v.y, &nz, &cam.m.v.w);
		render_viewport(x1, y1, x2-x1, y2-y1, &cam.m.v.x, &cam.m.v.y, &nw, &cam.m.v.z);
		render_viewport(x2, y1, x3-x2, y2-y1, &cam.m.v.x, &cam.m.v.y, &cam.m.v.z, &cam.m.v.w);
	} else {
		render_viewport(x0, y0, x3-x0, y1-y0, &cam.m.v.x, &cam.m.v.y, &cam.m.v.z, &cam.m.v.w);
		render_viewport(x0, y1, x1-x0, y2-y1, &cam.m.v.x, &cam.m.v.y, &nw, &cam.m.v.z);
		render_viewport(x1, y1, x2-x1, y2-y1, &nx, &cam.m.v.y, &nz, &cam.m.v.w);
		render_viewport(x2, y1, x3-x2, y2-y1, &nx, &cam.m.v.y, &cam.m.v.w, &cam.m.v.z);
	}

	SDL_UnlockSurface(screen);

	SDL_Flip(screen);

	fps_counter++;
	if(SDL_GetTicks() >= fps_next_tick)
		refresh_fps();
	
	v4f_t b0, b1;
	const kd_t *kd = NULL;//kd_in_tree(kdroot, &cam.o, &b0, &b1);
	if(kd != NULL)
	{
		printf("%p: %f %f %f %f -> %f %f %f %f\n",
			kd,
			b0.v.x,
			b0.v.y,
			b0.v.z,
			b0.v.w,
			b1.v.x,
			b1.v.y,
			b1.v.z,
			b1.v.w);
	}
}

void render_main(void)
{
	fps_counter = 0;
	fps_next_tick = SDL_GetTicks() + 1000;

	SDL_Event ev;

	quitflag = 0;
	float vaxz = 0.0f;
	float vayz = 0.0f;
	float vaxw = 0.0f;
	float vayw = 0.0f;
	float vaxy = 0.0f;
	float vx = 0.0f;
	float vy = 0.0f;
	float vz = 0.0f;
	float vw = 0.0f;
	while(!quitflag)
	{
		render_screen();

		const float vas = 0.002f;
		cam_rotate_by(vaxz*vas, vayz*vas, vaxw*vas, vayw*vas, vaxy*vas);
		vaxz = vayz = vaxw = vayw = vaxy = 0.0f;

		const float vs = 0.2f;

		// trace motion
		v4f_t no, tno, tv;
		no.m = _mm_setzero_ps();
		if(mbutts & 4)
		{
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam.m.v.x.m, _mm_set1_ps(-vx*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam.m.v.y.m, _mm_set1_ps(vy*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam.m.v.z.m, _mm_set1_ps(vw*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam.m.v.w.m, _mm_set1_ps(vz*vs)));
		} else {
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam.m.v.x.m, _mm_set1_ps(vx*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam.m.v.y.m, _mm_set1_ps(vy*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam.m.v.z.m, _mm_set1_ps(vz*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam.m.v.w.m, _mm_set1_ps(vw*vs)));
		}

		// add gravity
		grav_v += 0.30f*vs*vs;
		if(grav_v > 0.5f)
			grav_v = 0.5f;
		no.v.y = grav_v;

		// check distance
		v4f_t tv2;
		tv2.m = _mm_mul_ps(no.m, no.m);
		//float md = sqrtf(vx*vx + vy*vy + vz*vz + vw*vw)*vs;
		float md = tv2.v.x + tv2.v.y + tv2.v.z + tv2.v.w;

		if(md > 0.0000001f)
		{
			// normalise for direction
			tv.m = no.m;
			vec_norm(&tv);

			// cast a ray
			float r = 0.2f;
			tno.m = cam.o.m;
			int side;

			// trace away
			for(;;)
			{
				float d = trace_box(root, &tno, &tv, NULL, NULL, NULL, md + r, &side);
				
				//
				// apply collision
				//

				if(d < 0.0f)
				{
					// no collision. jump to point.
					cam.o.m = _mm_add_ps(cam.o.m,
						_mm_mul_ps(_mm_set1_ps(md), tv.m));

					break;
				} else {
					// we've hit a plane. slide back.
					if(side == F_YP)
					{
						if(grav_v >= 0.0f)
						{
							grav_v = 0.0f;
							grounded = 1;
						}
					} else if(side == F_YN) {
						if(grav_v <= 0.0f)
							grav_v = 0.0f;
					}

					float dd = md - (d - r);

					cam.o.m = _mm_add_ps(cam.o.m,
						_mm_mul_ps(_mm_set1_ps(md - dd), tv.m));

					// mask out velocity.
					tv.a[side&3] = 0.0f;

					// reduce distance.
					md = dd;
				}
			}
		}

		while(SDL_PollEvent(&ev))
		switch(ev.type)
		{
			case SDL_QUIT:
				quitflag = 1;
				break;

			case SDL_MOUSEBUTTONDOWN:
				mbutts |= (1<<(ev.button.button-1));
				if(ev.button.button == 4)
				{
					vaxy += 40.0f;
				} else if(ev.button.button == 5) {
					vaxy -= 40.0f;
				}
				break;

			case SDL_MOUSEBUTTONUP:
				mbutts &= ~(1<<(ev.button.button-1));
				break;

			case SDL_MOUSEMOTION:
				if(mrelease)
					break;
				if(mbutts & 4)
				{
					vayz += ev.motion.xrel;
					vaxz += ev.motion.yrel;
				} else {
					vayw += ev.motion.xrel;
					vaxw += ev.motion.yrel;
				}
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
				case SDLK_f:
				case SDLK_r:
					vy = 0.0f;
					break;

				default:
					break;
			} break;

			case SDL_KEYDOWN:
			switch(ev.key.keysym.sym)
			{
				case SDLK_ESCAPE:
					quitflag = 1;
					break;

				case SDLK_F5:
					if(mrelease)
					{
						SDL_WM_GrabInput(1);
						SDL_ShowCursor(0);
					} else {
						SDL_ShowCursor(1);
						SDL_WM_GrabInput(0);
					}
					mrelease = !mrelease;
					break;

				case SDLK_SPACE:
					if(grounded)
					{
						grav_v = -0.3f;
						grounded = 0;
					}
					break;
				case SDLK_s:
					vz = -1.0f;
					break;
				case SDLK_a:
					vx = -1.0f;
					break;
				case SDLK_q:
					vw = -1.0f;
					break;
				case SDLK_r:
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
				case SDLK_f:
					vy = +1.0f;
					break;
				
				default:
					break;
			} break;
		}
	}
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

