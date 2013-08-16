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

void game_player_tick(player_t *pl, float dt)
{
	if(pl->magic != 0xC4 && pl->magic != 0xC9 && pl->magic != 0x66 && pl->magic != 0x69)
		return;
	
	camera_t *cam = &(pl->cam);

	float vs = 0.15f*100.0f*dt;

	// trace motion
	v4f_t no, tno, tv;
	no.m = _mm_setzero_ps();
	if(pl->magic != 0x66 && pl->magic != 0x69)
	{
		if(pl->vflip)
		{
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam->m.v.x.m, _mm_set1_ps(-pl->lv.v.x*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam->m.v.y.m, _mm_set1_ps(pl->lv.v.y*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam->m.v.z.m, _mm_set1_ps(pl->lv.v.w*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam->m.v.w.m, _mm_set1_ps(pl->lv.v.z*vs)));
		} else {
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam->m.v.x.m, _mm_set1_ps(pl->lv.v.x*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam->m.v.y.m, _mm_set1_ps(pl->lv.v.y*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam->m.v.z.m, _mm_set1_ps(pl->lv.v.z*vs)));
			no.m = _mm_add_ps(no.m, _mm_mul_ps(cam->m.v.w.m, _mm_set1_ps(pl->lv.v.w*vs)));
		}
	}

	// add gravity
	pl->grav_v += 0.30f*vs*vs;
	if(pl->grav_v > 0.5f)
		pl->grav_v = 0.5f;
	no.v.y = pl->grav_v;

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
		float r = 0.5f;
		tno.m = cam->o.m;
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
				cam->o.m = _mm_add_ps(cam->o.m,
					_mm_mul_ps(_mm_set1_ps(md), tv.m));

				break;
			} else {
				// we've hit a plane. slide back.
				if(side == F_YP)
				{
					if(pl->grav_v >= 0.0f)
					{
						pl->grav_v = 0.0f;
						pl->grounded = 1;
					}
				} else if(side == F_YN) {
					if(pl->grav_v <= 0.0f)
						pl->grav_v = 0.0f;
				}

				float dd = md - (d - r);

				cam->o.m = _mm_add_ps(cam->o.m,
					_mm_mul_ps(_mm_set1_ps(md - dd), tv.m));

				// mask out velocity.
				tv.a[side&3] = 0.0f;

				// reduce distance.
				md = dd;
			}
		}
	}
}

int game_input(player_t *pl, float dt)
{
	int quitflag = 0;

	camera_t *cam = &(pl->cam);

	SDL_Event ev;

	float vaxz = 0.0f;
	float vayz = 0.0f;
	float vaxw = 0.0f;
	float vayw = 0.0f;
	float vaxy = 0.0f;

	vaxz = vayz = vaxw = vayw = vaxy = 0.0f;

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
				pl->lv.v.z = 0.0f;
				break;
			case SDLK_a:
			case SDLK_d:
				pl->lv.v.x = 0.0f;
				break;
			case SDLK_q:
			case SDLK_e:
				pl->lv.v.w = 0.0f;
				break;
			case SDLK_f:
			case SDLK_r:
				pl->lv.v.y = 0.0f;
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
				if(pl->grounded)
				{
					pl->grav_v = -0.3f;
					pl->grounded = 0;
				}
				break;
			case SDLK_s:
				pl->lv.v.z = -1.0f;
				break;
			case SDLK_a:
				pl->lv.v.x = -1.0f;
				break;
			case SDLK_q:
				pl->lv.v.w = -1.0f;
				break;
			case SDLK_r:
				pl->lv.v.y = -1.0f;
				break;
			case SDLK_w:
				pl->lv.v.z = +1.0f;
				break;
			case SDLK_d:
				pl->lv.v.x = +1.0f;
				break;
			case SDLK_e:
				pl->lv.v.w = +1.0f;
				break;
			case SDLK_f:
				pl->lv.v.y = +1.0f;
				break;
			
			default:
				break;
		} break;
	}

	const float vas = 0.002f;
	cam_rotate_by(cam, vaxz*vas, vayz*vas, vaxw*vas, vayw*vas, vaxy*vas);
	pl->vflip |= mbutts & 4;

	return quitflag;
}

