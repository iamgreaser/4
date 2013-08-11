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

float box_volume(box_t *b)
{
	v4f_t v;

	v.m = _mm_sub_ps(b->v1.m, b->v0.m);

	return v.v.x * v.v.y * v.v.z * v.v.w;
}

void box_depth_change_up(box_t *b)
{
	if(b == NULL)
		return;
	
	b->depth = b->c[(b->c[0]->depth > b->c[1]->depth ? 0 : 1)]->depth + 1;
	return box_depth_change_up(b->p);
}

void box_free(box_t *b)
{
	uint8_t *p = (uint8_t *)b;
	p -= b->ptroffs;
	free(b);
}

void box_free_tree(box_t *b)
{
	if(b == NULL)
		return;
	
	box_free_tree(b->c[0]);
	box_free_tree(b->c[1]);
	b->c[0] = b->c[1] = NULL;
	box_free(b);
}

box_t *box_new(v4f_t *v0, v4f_t *v1, v4f_t *color, int op)
{
	box_t *b = malloc(sizeof(box_t) + 16);
	
	int ptroffs = 0;

	while((((long)b)&15) != 0)
	{
		b = (box_t *)(((uint8_t *)b) + 1);
		ptroffs++;
	}

	b->ptroffs = ptroffs;

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
	box_depth_change_up(p->p);
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

int box_pair_touch(box_t *b0, box_t *b1)
{
	v4f_t cmp;

	cmp.m = _mm_or_ps(
		_mm_cmpge_ps(b0->v1.m, b1->v0.m),
		_mm_cmpge_ps(b1->v1.m, b0->v0.m));
	
	return (_mm_movemask_ps(cmp.m) != 0);

}

box_t *box_in_tree(box_t *box, v4f_t *p, box_t **ignore, int ignore_count)
{
	int i;

	// check if null
	if(box == NULL)
		return NULL;
	
	// check if in this box,
	if(!box_in(box, p))
		return NULL;
	
	// check if we're a pair
	if(box->op != SHP_PAIR)
	{
		// check if we're ignoring this box
		for(i = 0; i < ignore_count; i++)
			if(box == ignore[i])
				return NULL;

		// return.
		return box;
	}
	
	// traverse children
	box_t *b = box_in_tree(box->c[0], p, ignore, ignore_count);
	if(b == NULL)
		b = box_in_tree(box->c[1], p, ignore, ignore_count);
	
	// return
	return b;
}

void box_normal(box_t *box, v4f_t *p, v4f_t *n, int inside)
{
	// get dims + centre
	__m128 dims = _mm_sub_ps(box->v1.m, box->v0.m);
	__m128 cent = _mm_mul_ps(_mm_add_ps(box->v0.m, box->v1.m), _mm_set1_ps(0.5f));

	// get abs offset
	__m128 offs = _mm_sub_ps(p->m, cent);
	__m128 offs_abs = _mm_and_ps(offs, (__m128)_mm_set1_epi32(0x7FFFFFFF));
	offs_abs = _mm_div_ps(offs_abs, dims);

	// scale offset

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

	// TODO: suggested optimisation to apply:
	// if max(in_dist) < min(out_dist), then we are inside the box.

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


