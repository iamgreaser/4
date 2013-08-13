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

kd_t *kd_new(int axis, float v0, float v1, box_t *box, box_t *obox0, box_t *obox1)
{
	kd_t *kd = malloc(sizeof(kd_t));

	kd->axis = axis;
	kd->v[0] = v0;
	kd->v[1] = v1;
	kd->box = box;
	kd->obox[0] = obox0;
	kd->obox[1] = obox1;

	kd->c[0] = kd->c[1] = kd->p = NULL;

	return kd;
}

kd_t *kd_dup(kd_t *kd, kd_t *p)
{
	if(kd == NULL)
		return NULL;
	
	kd_t *ret = malloc(sizeof(kd_t));
	memcpy(ret, kd, sizeof(kd_t));

	ret->c[0] = kd_dup(ret->c[0], ret);
	ret->c[1] = kd_dup(ret->c[1], ret);
	ret->p = p;

	return ret;
}

void kd_free(kd_t *kd)
{
	if(kd == NULL)
		return;
	
	free(kd);
}

void kd_free_down(kd_t *kd)
{
	if(kd == NULL)
		return;
	
	kd_free_down(kd->c[0]);
	kd_free_down(kd->c[1]);
	kd_free(kd);
}

int kd_equal(const kd_t *k0, const kd_t *k1, int cmpv)
{
	// null check
	if(k0 == NULL)
		return k1 == NULL;
	if(k1 == NULL)
		return 0;
	
	// basic constants
	if(k0->axis != k1->axis) return 0;
	if(k0->box != k1->box) return 0;

	if(cmpv)
	{
		if(k0->v[0] != k1->v[0]) return 0;
		if(k0->v[1] != k1->v[1]) return 0;
	}

	// children
	if(k0->c[0] != NULL && k1->c[0] == NULL) return 0;
	if(k0->c[0] == NULL && k1->c[0] != NULL) return 0;
	
	// trace down
	return (k0->c[0] == NULL || kd_equal(k0->c[0], k1->c[0], cmpv))
		&& (k1->c[1] == NULL || kd_equal(k0->c[1], k1->c[1], cmpv));
}

kd_t *kd_contract_down(kd_t *kd, kd_t *cmpkd)
{
	// some simple checks
	if(kd == NULL || cmpkd == NULL || kd->axis != cmpkd->axis)
		return kd;
	
	// check if out of range
	if(kd->v[0] > cmpkd->v[1] || kd->v[1] < cmpkd->v[0])
	{
		// this node no longer exists.
		kd_free_down(kd);
		return NULL;
	}
	
	// contract
	if(kd->v[0] < cmpkd->v[0])
		kd->v[0] = cmpkd->v[0];
	if(kd->v[1] > cmpkd->v[1])
		kd->v[1] = cmpkd->v[1];

	// check node type
	if(kd->c[1] != NULL)
	{
		// pair. check each child.
		kd_t *c0 = kd_contract_down(kd->c[0], cmpkd);
		kd_t *c1 = kd_contract_down(kd->c[1], cmpkd);

		// are these both non-null?
		if(c0 != NULL && c1 != NULL)
		{
			// combine and return ourselves.
			kd->c[0] = c0;
			kd->c[1] = c1;
			c0->p = kd;
			c1->p = kd;

			return kd;
		} else {
			// one is NULL, the other has a value.
			kd = c0;
			if(kd == NULL)
				kd = c1;

			// assertion...
			if(kd == NULL)
				abort();

			// return the child node that we have.
			return kd;
		}
	}

	// otherwise, axis change OR leaf.
	// return ourselves.
	return kd;
}

void kd_expand_down(kd_t *kd, int child)
{
	if(kd == NULL || kd->c[child] == NULL || kd->axis != kd->c[child]->axis)
		return;
	
	kd->c[child]->v[child] = kd->v[child];

	return kd_expand_down(kd->c[child], child);
}

kd_t *kd_split(kd_t *kd, float split)
{
	// is the split in range?
	if(split <= kd->v[0] || split >= kd->v[1])
		// nope. return.
		return kd;
	
	// make split.
	kd_t *p = kd_new(kd->axis, kd->v[0], kd->v[1], kd->box, kd->obox[0], kd->obox[1]);
	p->p = kd->p;

	// duplicate.
	kd_t *dup = kd_dup(kd, p);

	// set axes.
	dup->v[0] = kd->v[1] = split;

	// contract down.
	dup = kd_contract_down(dup, dup);
	kd = kd_contract_down(kd, kd);

	// set children.
	p->c[0] = kd;
	p->c[1] = dup;
	p->c[0]->p = p;
	p->c[1]->p = p;

	// return "parent".
	return p;
}

kd_t *kd_add_box_step(int axis, kd_t *kd, box_t *box)
{
	//printf("%i %016llX %016llX\n", axis, kd, box);

	// check: do we even have a box?
	if(box == NULL)
		return kd;

	// check: do we need to create a new node?
	if(kd == NULL)
	{
		// yes. create a new kd node.
		kd = kd_new(axis, box->v0.a[axis], box->v1.a[axis],
			(axis != 3 || box->op == SHP_SUB ? NULL : box),
			(box->op == SHP_ADD ? NULL : box),
			(box->op == SHP_ADD ? NULL : box));

		// check if we're at axis level 3.
		if(axis == 3)
			return kd;

		// nope. create an axis change child.
		kd->c[0] = kd_add_box_step(axis + 1, NULL, box);
		kd->c[0]->p = kd;

		return kd;
	}

	// check: do we not reach the -ve limit?
	if(box->v0.a[axis] > kd->v[0])
	{
		// split this node and recurse.
		kd = kd_split(kd, box->v0.a[axis]);
		return kd_add_box_step(axis, kd->c[1], box);
	}

	// check: do we not reach the +ve limit?
	if(box->v1.a[axis] < kd->v[1])
	{
		// split this node and recurse.
		kd = kd_split(kd, box->v1.a[axis]);
		return kd_add_box_step(axis, kd->c[0], box);
	}

	// check: is this a leaf?
	if(kd->c[0] == NULL)
	{
		// replace this node.
		kd_free_down(kd);
		return kd_add_box_step(axis, NULL, box);
	}

	// then this is an axis switch.
	// step down.
	kd->c[0] = kd_add_box_step(axis + 1, kd->c[0], box);
	kd->c[0]->p = kd;
	return kd;
}

kd_t *kd_add_box(kd_t *kd, box_t *box)
{
	if(box == NULL)
		return NULL;
	
	if(box->op == SHP_PAIR)
	{
		kd = kd_add_box(kd, box->c[0]);
		kd = kd_add_box(kd, box->c[1]);
		return kd;
	} else {
		return kd_add_box_step(0, kd, box);
	}
}

