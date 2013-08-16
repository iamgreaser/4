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

sphere_t *sphere_list_add(sphere_t *l, int *llen, const v4f_t *v, float r, const v4f_t *color, int plidx)
{
	// do reallocation
	(*llen)++;
	l = realloc16(l, sizeof(sphere_t) * *llen);

	// get current sphere
	sphere_t *csph = &l[(*llen)-1];

	// fill it in
	csph->v.m = v->m;
	csph->r = r;
	csph->plidx = plidx;
	csph->color.m = color->m;

	// return our new list
	return l;
}

void sphere_list_rm(sphere_t *l, int *llen, int idx)
{
	sphere_t *csph = l;

	// assertion!
	if(idx < 0 || idx >= *llen)
		abort();

	// scrub this sphere out
	int mlen = ((*llen)-1) - idx;
	if(mlen > 0)
		memmove(csph, csph+1, mlen*sizeof(sphere_t));

	// decrement length
	(*llen)--;
}

float sphere_trace_one(sphere_t *s, const v4f_t *p, const v4f_t *v, float md)
{
	(void)md;

	// based off some code i wrote for iceball
	// which was based off some math i wrote

	// get delta and delta squared
	__m128 delta = _mm_sub_ps(s->v.m, p->m);
	v4f_t delta2;
	delta2.m = _mm_mul_ps(delta, delta);

	// get distance
	float dd = delta2.v.x + delta2.v.y + delta2.v.z + delta2.v.w;
	float ddl = sqrtf(dd);

	// cut it short
	if(ddl - s->r > md)
		return -1.0f;

	// get dot product for some particular side
	v4f_t dotb;
	dotb.m = _mm_mul_ps(delta, v->m);
	float dotk = dotb.v.x + dotb.v.y + dotb.v.z + dotb.v.w;

	// make sure we're not hitting on the back side of the ray
	if(dotk <= 0.0f)
		return -1.0f;
	
	// get the radius of intersection
	float sr2 = s->r * s->r;
	float cr2 = dd - dotk*dotk;

	// now we have to come up with our own stuff.
	// check if the intersection radius is in range.
	if(cr2 > sr2)
		return -1.0f; // nope!
	
	// this is where i use a piece of paper and do weird semi-Euclidean shit.
	float sh2 = sr2 - cr2;
	float sh = sqrtf(sh2);
	
	float bo = ddl-sh;

	// if this is < 0.0f... well, it doesn't matter.
	// we don't need to trace from within the inside.
	
	float d = sqrtf(cr2 + bo*bo);

	return d;
}

void sphere_normal(const sphere_t *s, const v4f_t *p, v4f_t *normal)
{
	normal->m = _mm_sub_ps(p->m, s->v.m);
	vec_norm(normal);
}

float sphere_trace(sphere_t *l, const int *llen, const v4f_t *p, const v4f_t *v, float md, sphere_t **rets)
{
	sphere_t *fs = NULL;
	float fd = -1.0f;
	int i;

	for(i = 0; i < *llen; i++)
	{
		sphere_t *cs = &l[i];
		float cd = sphere_trace_one(cs, p, v, md);

		if(cd > 0.0f && cd < md && (fd < 0.0f || cd < fd))
		{
			fs = cs;
			fd = cd;
		}
	}

	if(fs == NULL)
		return -1.0f;
	
	if(rets != NULL) *rets = fs;
	return fd;
}

