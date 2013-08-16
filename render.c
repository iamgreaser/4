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

int last_tick = 0;

void render_viewport(camera_t *cam, int vx, int vy, int w, int h, v4f_t *dx, v4f_t *dy, v4f_t *dz, v4f_t *dw)
{
	int y;

	float sxd = (2.0f)/(w/2.0f);
	float syd = (2.0f)/(w/2.0f);
	float sxi = -sxd*w/2.0f;
	float syi = -syd*h/2.0f;

	// find starting box
	box_t *bstart = box_in_tree(root, &cam->o, NULL, 0);
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
			uint32_t v = trace_pixel(cam, bstart, sx, sy, dx, dy, dz, dw, &seed);

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
			vs[0] = trace_pixel(cam, bstart, sx, sy, dx, dy, dz, dw, &seed); sx += sxd;
			vs[1] = trace_pixel(cam, bstart, sx, sy, dx, dy, dz, dw, &seed); sx += sxd;
			vs[2] = trace_pixel(cam, bstart, sx, sy, dx, dy, dz, dw, &seed); sx += sxd;
			vs[3] = trace_pixel(cam, bstart, sx, sy, dx, dy, dz, dw, &seed); sx += sxd;

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
			uint32_t v = trace_pixel(cam, bstart, sx, sy, dx, dy, dz, dw, &seed);

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

void render_screen(camera_t *cam)
{
	int y0 = 0;
	int y1 = rtbuf_height*3/4;
	int y2 = rtbuf_height;

	int x0 = 0;
	int x1 = rtbuf_width/3;
	int x2 = rtbuf_width*2/3;
	int x3 = rtbuf_width;

	v4f_t nx,  nz, nw;
	nx.m = _mm_mul_ps(cam->m.v.x.m, _mm_set1_ps(-1.0f));
	//ny.m = _mm_mul_ps(cam->m.v.y.m, _mm_set1_ps(-1.0f));
	nz.m = _mm_mul_ps(cam->m.v.z.m, _mm_set1_ps(-1.0f));
	nw.m = _mm_mul_ps(cam->m.v.w.m, _mm_set1_ps(-1.0f));

	SDL_LockSurface(screen);

	if(mbutts & 4)
	{
		render_viewport(cam, x0, y0, x3-x0, y1-y0, &nx, &cam->m.v.y, &cam->m.v.w, &cam->m.v.z);
		render_viewport(cam, x0, y1, x1-x0, y2-y1, &nx, &cam->m.v.y, &nz, &cam->m.v.w);
		render_viewport(cam, x1, y1, x2-x1, y2-y1, &cam->m.v.x, &cam->m.v.y, &nw, &cam->m.v.z);
		render_viewport(cam, x2, y1, x3-x2, y2-y1, &cam->m.v.x, &cam->m.v.y, &cam->m.v.z, &cam->m.v.w);
	} else {
		render_viewport(cam, x0, y0, x3-x0, y1-y0, &cam->m.v.x, &cam->m.v.y, &cam->m.v.z, &cam->m.v.w);
		render_viewport(cam, x0, y1, x1-x0, y2-y1, &cam->m.v.x, &cam->m.v.y, &nw, &cam->m.v.z);
		render_viewport(cam, x1, y1, x2-x1, y2-y1, &nx, &cam->m.v.y, &nz, &cam->m.v.w);
		render_viewport(cam, x2, y1, x3-x2, y2-y1, &nx, &cam->m.v.y, &cam->m.v.w, &cam->m.v.z);
	}

	SDL_UnlockSurface(screen);

	SDL_Flip(screen);

	fps_counter++;
	if(SDL_GetTicks() >= fps_next_tick)
		refresh_fps();
	
	v4f_t b0, b1;
	const kd_t *kd = NULL;//kd_in_tree(kdroot, &cam->o, &b0, &b1);
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
	int i;

	fps_counter = 0;
	fps_next_tick = SDL_GetTicks() + 1000;

	quitflag = 0;
	float dt = 0.00f;
	last_tick = SDL_GetTicks();
	while(!quitflag)
	{
		player_t *pl = &players[cplr >= 0 ? cplr : 127];
		camera_t *cam = &(pl->cam);

		render_screen(cam);

		if(net_mode == 1)
			net_update_client();
		else if(net_mode == 3)
			net_update_server();

		quitflag = quitflag || game_input(pl, dt);
		
		player_sphere_clear();
		for(i = 0; i < PLAYERS_MAX; i++)
		{
			player_t *pl = &players[i];
			if(pl->magic == 0xC4)
				player_sphere_add(pl);
		}

		for(i = 0; i < PLAYERS_MAX; i++)
			game_player_tick(&players[i], dt);

		int next_tick = SDL_GetTicks();
		dt = ((float)(next_tick - last_tick))/1000.0f;
		last_tick = next_tick;
	}
}

