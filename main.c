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

int fps_counter = 0;
int fps_next_tick = 0;
box_t *root = NULL;

SDL_Surface *screen = NULL;
uint32_t *rtbuf = NULL;
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
	box = box_in_tree(box, p, ign_l, ign_c);

	// if this fails, we need to trace to the first surface we hit...
	// TODO? not sure if this is necessary.
	// at the moment, we're just going to fail, because you shouldn't be inside a SHP_ADD,
	// and should always be inside at least one SHP_SUB.
	if(box == NULL)
		return -1.0f;

	// otherwise, trace onwards.
	tp.m = p->m;
	for(;;)
	{
		// trace to the box end.
		td = box_crosses(box, &tp, v, &ins);
		if(td < 0.0f) // wat, pretty sure this should just crash and burn as an assertion
			return -1.0f;

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
		//tp.m = _mm_add_ps(p->m, _mm_mul_ps(_mm_set1_ps(ad), v->m));
		tp.m = _mm_add_ps(tp.m, _mm_mul_ps(_mm_set1_ps(td), v->m));

		// check to see if we're in another box.
		box_t *nbox = box_in_tree(root, &tp, ign_l, ign_c);
		if(nbox == NULL || nbox->op != SHP_SUB)
		{
			// nope. let's end it here.
			if(nbox == NULL)
				nbox = box;

			*retbox = nbox;
			*inside = ins;
			color->m = nbox->color.m;
			return ad;
		}

		box = nbox;
	}
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

	fps_counter++;
	if(SDL_GetTicks() >= fps_next_tick)
		refresh_fps();
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

void render_main(void)
{
	fps_counter = 0;
	fps_next_tick = SDL_GetTicks() + 1000;

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
	if(argc <= 1)
	{
		printf("%s dat/midbox.4lv (or something like that).\n", argv[0]);
		return 99;
	}

	SDL_Init(SDL_INIT_TIMER | SDL_INIT_AUDIO | SDL_INIT_VIDEO);
	SDL_WM_SetCaption("4 (pronounced \"sardinelauncher\")", NULL);
	screen = SDL_SetVideoMode(rtbuf_width * rtbuf_scale, rtbuf_height * rtbuf_scale, 32, 0);
	rtbuf = malloc(rtbuf_width * rtbuf_height * 4);

	cam_init();
	level_init(argv[1]);
	render_main();

	free(rtbuf);

	return 0;
}

