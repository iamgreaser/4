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

int box_in(const box_t *box, const v4f_t *p)
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

const box_t *box_in_tree_down(const box_t *box, const v4f_t *p, box_t **ignore, int ignore_count)
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
	const box_t *b = box_in_tree_down(box->c[0], p, ignore, ignore_count);
	if(b == NULL)
		b = box_in_tree_down(box->c[1], p, ignore, ignore_count);
	
	// return
	return b;
}

const box_t *box_in_tree_up(const box_t *box, const v4f_t *p, box_t **ignore, int ignore_count, const box_t *prev)
{
	int i;

	// check if null
	if(box == NULL || box->op == SHP_SCENE_LIST)
		return NULL;
	
	// check if in this box
	if(!box_in(box, p))
		return box_in_tree_up(box->p, p, ignore, ignore_count, box);
	
	// check if we're a pair
	if(box->op != SHP_PAIR)
	{
		// check if we're ignoring this box
		int fail = 0;
		for(i = 0; i < ignore_count; i++)
			if(box == ignore[i])
			{
				fail = 1;
				break;
			}

		if(!fail)
			return box; // just return us.
	}
	
	// trace down
	const box_t *b = NULL;
	if(box->c[0] != prev)
		b = box_in_tree_down(box->c[0], p, ignore, ignore_count);
	if(b == NULL && box->c[1] != prev)
		b = box_in_tree_down(box->c[1], p, ignore, ignore_count);
	if(b == NULL)
		b = box_in_tree_up(box->p, p, ignore, ignore_count, box);
	
	return b;
}

box_t *box_in_tree(box_t *box, const v4f_t *p, box_t **ignore, int ignore_count)
{
	return (box_t *)box_in_tree_up(box, p, ignore, ignore_count, NULL);
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

float box_crosses(const box_t *box, const v4f_t *p, const v4f_t *vi, int *inside, int *side)
{
	// the general idea:
	// define the "out" faces to be in the direction of the ray velocity.
	// define the "in" faces to be the other faces.
	// if the largest "in" face distance is smaller than the smallest "out" face distance,
	// return the largest "in" face distance - we've hit the box.
	// otherwise, return -1.0f.

	// copy vi so we don't mutilate it
	v4f_t v;
	v.m = vi->m;

	// determine direction
	__m128 dirneg = _mm_cmpgt_ps(_mm_setzero_ps(), v.m);

	// get in/out values
	__m128 fi, fo, fd, fdm;
	fd = _mm_sub_ps(box->v1.m, box->v0.m);
	fdm = _mm_and_ps(dirneg, fd);
	fi = _mm_add_ps(box->v0.m, fdm);
	fo = _mm_sub_ps(box->v1.m, fdm);

	// clamp v to not be in the +/- epsilon range
	__m128 vflip = _mm_and_ps(dirneg, (__m128)_mm_set1_epi32(-0x80000000));
	v.m = _mm_xor_ps(vflip,
		_mm_max_ps(
			_mm_set1_ps(0.00000001f),
			_mm_xor_ps(vflip, v.m)
		));
	
	if(_mm_movemask_ps(_mm_xor_ps(vflip, v.m)) != 0)
		abort();

	// offset the values
	// doing a max so we're actually using -out
	fi = _mm_sub_ps(fi, p->m);
	fo = _mm_sub_ps(p->m, fo);

	// get times
	// encode the indices on the lowest 3 bits of the fractional part
	__m128 ti, to;
	ti = _mm_xor_ps(
		_mm_and_ps(dirneg, (__m128)_mm_set1_epi32(4)),
		_mm_or_ps((__m128)_mm_set_epi32(3, 2, 1, 0),
			_mm_and_ps((__m128)_mm_set1_epi32(~0x7),
				_mm_div_ps(fi, v.m)
			)
		)
	);
	to = _mm_xor_ps(
		_mm_and_ps(dirneg, (__m128)_mm_set1_epi32(4)),
		_mm_or_ps((__m128)_mm_set_epi32(7, 6, 5, 4),
			_mm_and_ps((__m128)_mm_set1_epi32(~0x7),
				_mm_div_ps(fo, v.m)
			)
		)
	);

	/*
	printf("%08X %08X %08X %08X\n",
		(*(0+(int *)&to)),
		(*(1+(int *)&to)),
		(*(2+(int *)&to)),
		(*(3+(int *)&to)));
	*/

	// get relevant face indices + distances
	v4f_t c0;
	c0.m = _mm_max_ps(
	        _mm_movehl_ps(to, ti),
	        _mm_movelh_ps(ti, to)
	);

	float dci = (c0.a[0] > c0.a[1] ? c0.a[0] : c0.a[1]);
	float dco = -(c0.a[2] > c0.a[3] ? c0.a[2] : c0.a[3]);
	int idxi = (*(int *)(float *)&dci) & 7;
	int idxo = (*(int *)(float *)&dco) & 7;

	// return appropriately
	if(dco < dci)
		return -1.0f;
	
	if(dci <= 0.0f)
	{
		if(inside != NULL) *inside = 1;
		if(side != NULL) *side = idxo+(v.a[idxo] < 0.0f ? 0 : 4);
		return dco;
	} else {
		if(inside != NULL) *inside = 0;
		if(side != NULL) *side = idxi+(v.a[idxi] < 0.0f ? 4 : 0);
		return dci;
	}
}


