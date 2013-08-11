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

void cam_rotate_by(float axz, float ayz, float axw, float ayw)
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

	// rotate around xz
	mat_ident(&rotmat);
	s = sinf(axz);
	c = cosf(axz);
	rotmat.v.y.m = _mm_set_ps(-s, 0, c, 0);
	rotmat.v.w.m = _mm_set_ps(c, 0, s, 0);
	mat_mul(&cam.m, &rotmat);

	// rotate around yz
	mat_ident(&rotmat);
	s = sinf(ayz);
	c = cosf(ayz);
	rotmat.v.x.m = _mm_set_ps(s, 0, 0, c);
	rotmat.v.w.m = _mm_set_ps(c, 0, 0, -s);
	mat_mul(&cam.m, &rotmat);

	// correct up vector by rotating around zw
	mat_ident(&rotmat);
	s = cam.m.v.x.v.y;
	c = sqrtf(1.0f - s*s);
	rotmat.v.x.m = _mm_set_ps(0, 0, -s, c);
	rotmat.v.y.m = _mm_set_ps(0, 0, c, s);
	mat_mul(&cam.m, &rotmat);
}


