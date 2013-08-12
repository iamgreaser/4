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

kd_t *kd_new(int axis, float v0, float v1, box_t *box, box_t *obox)
{
	kd_t *kd = malloc(sizeof(kd_t));

	kd->axis = axis;
	kd->v[0] = v0;
	kd->v[1] = v1;
	kd->box = box;
	kd->obox = obox;

	kd->c[0] = kd->c[1] = kd->p;

	return kd;
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

void kd_expand_down(kd_t *kd, int child)
{
	if(kd == NULL || kd->c[child] == NULL || kd->axis != kd->c[child]->axis)
		return;
	
	kd->c[child]->v[child] = kd->v[child];

	return kd_expand_down(kd->c[child], child);
}

kd_t *kd_add_box_step(int axis, kd_t *kd, box_t *box)
{
	// check: do we need to create a new node?
	if(kd == NULL)
	{
		// yes. create a new kd node.
		kd = kd_new(axis, box->v0.a[axis], box->v1.a[axis],
			(axis != 3 || box->op == SHP_SUB ? NULL : box),
			(box->op == SHP_ADD ? NULL : box));

		// check if we're at axis level 3.
		if(axis == 3)
			return kd;

		// nope. create an axis change child.
		kd->c[0] = kd_add_box_step(axis + 1, kd->c[0], box);
		kd->c[0]->p = kd;

		return kd;
	}

	// check: are we on a leaf node (axis 3 only) or axis change node?
	// (we can merge the nodes later)
	if(kd->c[1] == NULL)
	{
		// check: do we hit the -ve limit? if we pass it, it's OK.
		if(box->v0.a[axis] <= kd->v[0])
		{
			// no. split this node.
			kd_t *old_kd = kd;
			kd = kd_new(axis, kd->v[0], box->v0.a[axis], NULL, kd->obox);

			// link it up correctly
			kd->p = old_kd->p;
			kd->c[0] = old_kd;
			old_kd->p = kd;

			// trace further.
			// TODO: recursion
		}
	}

	// check: are we on a 
}

kd_t *kd_add_box(kd_t *kd, box_t *box)
{
	if(box == NULL)
		return;
	
	if(box->op == SHP_PAIR)
	{
		kd = kd_add_box(kd, box->c[0]);
		kd = kd_add_box(kd, box->c[1]);
		return kd;
	} else {
		return kd_add_box_step(0, kd, box);
	}
}

