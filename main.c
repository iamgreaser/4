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

int mbutts = 0;
int mrelease = 0;

int fps_counter = 0;
int fps_next_tick = 0;
box_t *root = NULL;

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

float trace_into_box(box_t **retbox, v4f_t *p, v4f_t *v, float md)
{
	box_t *box = *retbox;
	float d = 0.0f;

	// check if we're in this box.
	// (skipping SHP_PAIR as box_crosses isn't the fastest op ever)
	if(box->op != SHP_PAIR && !box_in(box, p))
	{
		// nope. can we trace to it?
		//d = box_crosses_outside(box, p, v);
		d = box_crosses(box, p, v, NULL);

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
		float d0 = trace_into_box(&b0, p, v, md);
		if(d0 >= 0.0f)
			md = d0;
		float d1 = trace_into_box(&b1, p, v, md);

		if(d0 < 0.0f || (d1 >= 0.0f && d1 < d0))
		{
			if(d1 >= 0.0f) *retbox = b1;
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
float trace_box(box_t *box, v4f_t *p, v4f_t *v, v4f_t *color, box_t **retbox, int *inside, float md)
{
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
		td = trace_into_box(&box, p, v, md);

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
		td = box_crosses(box, &tp, v, &ins);

		// ASSERTION.
		if(td < 0.0f)
			abort();

		// check if we have exceeded our distance.
		ad += td;
		if(ad >= md)
			// we have. fail.
			return -1.0f;
		
		// moving right along...
		// TODO: clear the ignore buffer more efficiently
		ign_l[ign_c++] = box;
		if(ign_c > 100)
			abort();
		tp.m = _mm_add_ps(p->m, _mm_mul_ps(_mm_set1_ps(ad), v->m));
		//tp.m = _mm_add_ps(tp.m, _mm_mul_ps(_mm_set1_ps(td), v->m));

		// check to see if we're in another box.
		box_t *nbox = box_in_tree(root, &tp, ign_l, ign_c);
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

uint32_t trace_pixel(float sx, float sy, v4f_t *dirx, v4f_t *diry, v4f_t *dirz)
{
	v4f_t f;

	f.m = _mm_add_ps(
		dirz->m,
		_mm_add_ps(
			_mm_mul_ps(dirx->m, _mm_set1_ps(sx)),
			_mm_mul_ps(diry->m, _mm_set1_ps(sy))
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
	float d = trace_box(root, &p, &f, &color, &box, &inside, 100.0f);

	if(d >= 0.0f)
	{
		// move position
		p.m = _mm_add_ps(p.m, _mm_mul_ps(
			f.m, _mm_set1_ps(d)));

		// calculate normal
		v4f_t n;
		box_normal(box, &p, &n, inside);

		// calculate diffuse (against one point for now)
		v4f_t dc;
		dc.m = _mm_mul_ps(f.m, n.m);
		float diff = dc.v.x + dc.v.y + dc.v.z + dc.v.w;
		//if(diff < 0.0f) diff = -diff;

		// multiply diffuse
		color.m = _mm_mul_ps(color.m, _mm_set1_ps(diff));
	}

	return color_vec_sse(color.m);
}

void render_viewport(int vx, int vy, int w, int h, v4f_t *dx, v4f_t *dy, v4f_t *dz)
{
	int y;

	float sxd = (2.0f)/(w/2.0f);
	float syd = (2.0f)/(w/2.0f);
	float sxi = -sxd*w/2.0f;
	float syi = -syd*h/2.0f;

	// scale up onto screen
	#pragma omp parallel for
	for(y = 0; y < h; y++)
	{
		int x;
		uint32_t *d0, *d1, *d2;

		d0 = (uint32_t *)((screen->pixels) + screen->pitch*(0+3*(y+vy)) + 4*3*vx);
		d1 = (uint32_t *)((screen->pixels) + screen->pitch*(1+3*(y+vy)) + 4*3*vx);
		d2 = (uint32_t *)((screen->pixels) + screen->pitch*(2+3*(y+vy)) + 4*3*vx);

		float sy = syi + syd*y;

		float sx = sxi;

		for(x = 0; x < w; x++)
		{
			uint32_t v = trace_pixel(sx, sy, dx, dy, dz);

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
	box_print(root, 0);
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
		render_viewport(x0, y0, x3-x0, y1-y0, &nx, &cam.m.v.y, &cam.m.v.w);
		render_viewport(x0, y1, x1-x0, y2-y1, &nx, &cam.m.v.y, &nz);
		render_viewport(x1, y1, x2-x1, y2-y1, &cam.m.v.x, &cam.m.v.y, &nw);
		render_viewport(x2, y1, x3-x2, y2-y1, &cam.m.v.x, &cam.m.v.y, &cam.m.v.z);
	} else {
		render_viewport(x0, y0, x3-x0, y1-y0, &cam.m.v.x, &cam.m.v.y, &cam.m.v.z);
		render_viewport(x0, y1, x1-x0, y2-y1, &cam.m.v.x, &cam.m.v.y, &nw);
		render_viewport(x1, y1, x2-x1, y2-y1, &nx, &cam.m.v.y, &nz);
		render_viewport(x2, y1, x3-x2, y2-y1, &nx, &cam.m.v.y, &cam.m.v.w);
	}

	SDL_UnlockSurface(screen);

	SDL_Flip(screen);

	fps_counter++;
	if(SDL_GetTicks() >= fps_next_tick)
		refresh_fps();
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
	float vx = 0.0f;
	float vy = 0.0f;
	float vz = 0.0f;
	float vw = 0.0f;
	while(!quitflag)
	{
		render_screen();

		const float vas = 0.002f;
		cam_rotate_by(vaxz*vas, vayz*vas, vaxw*vas, vayw*vas);
		vaxz = vayz = vaxw = vayw = 0.0f;

		const float vs = 0.02f;

		// check velocity
		if(vx*vx + vy*vy + vz*vz + vw*vw > 0.000001f)
		{
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

			// normalise for direction
			tv.m = no.m;
			vec_norm(&tv);

			// cast a ray
			float r = 0.2f;
			float md = sqrtf(vx*vx + vy*vy + vz*vz + vw*vw)*vs;
			tno.m = cam.o.m;
			float d = trace_box(root, &tno, &tv, NULL, NULL, NULL, md + r);
			
			(void)d;
			if(d >= 0.0f) md = d - r;

			cam.o.m = _mm_add_ps(cam.o.m,
				_mm_mul_ps(_mm_set1_ps(md), tv.m));
		}

		while(SDL_PollEvent(&ev))
		switch(ev.type)
		{
			case SDL_QUIT:
				quitflag = 1;
				break;

			case SDL_MOUSEBUTTONDOWN:
				mbutts |= (1<<(ev.button.button-1));
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
				case SDLK_LCTRL:
				case SDLK_SPACE:
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

