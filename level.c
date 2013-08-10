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

int level_parse_vector(FILE *fp, v4f_t *dest)
{
	char det;

	if(fscanf(fp, " %c", &det) != 1)
	{
		fprintf(stderr, "level_parse_vector: expected '<', EOF found instead\n");
		return 0;
	}

	if(det != '<')
	{
		fprintf(stderr, "level_parse_vector: expected '<', '%c' found instead\n", det);
		return 0;
	}

	v4f_t v;
	int r = fscanf(fp, " %f , %f , %f , %f >", &v.v.x, &v.v.y, &v.v.z, &v.v.w);
	if(r != 4)
	{
		perror("level_parse_vector(actual parsing)");
		return 0;
	}

	dest->m = v.m;
	return 1;
}

box_t *level_parse_group(FILE *fp)
{
	box_t *r = NULL;

	char det;

	do
	{
		if(fscanf(fp, " %c", &det) != 1)
		{
			fprintf(stderr, "level_parse_group: expected '}', EOF found instead\n");
			return NULL;
		}

		switch(det)
		{
			case '}':
				break; // falls through to do/while block

			case '{': {
				box_t *b = level_parse_group(fp);

				if(b == NULL)
				{
					box_free_tree(r);
					return NULL;
				}

				r = box_inject(r, b);
			} break;
			
			case '#':
				fscanf(fp, "%*[^\r\n]");
				break;

			case '-':
			case '+': {
				v4f_t v0, v1, color;

				if(!(level_parse_vector(fp, &v0)
					&& level_parse_vector(fp, &v1)
					&& level_parse_vector(fp, &color)))
				{
					box_free_tree(r);
					return NULL;
				}

				if(fscanf(fp, " ;") == EOF)
				{
					perror("level_parse_group(expected ';')");
					box_free_tree(r);
					return NULL;
				}
				
				int typ = SHP_ADD;

				if(det == '-')
					typ = SHP_SUB;

				box_t *b = box_new(&v0, &v1, &color, typ);
				r = box_inject(r, b);
			} break;

			default:
				fprintf(stderr, "level_parse_group: expected '}', '%c' found instead\n", det);
				box_free_tree(r);
				return NULL;
		}
	} while(det != '}');

	if(r == NULL)
		fprintf(stderr, "level_parse_group: cannot have an empty group!\n");
		// follows through

	return r;
}

box_t *level_load_fp(FILE *fp)
{
	// check header
	char hbuf[16];
	if(fscanf(fp, "%15[^ \r\n\t{]", hbuf) != 1 || strcmp(hbuf, "4level"))
	{
		fprintf(stderr, "level_load_fp: not a 4level file\n");
		return NULL;
	}

	// load a group
	if(fscanf(fp, " %c", &hbuf[0]) != 1 || hbuf[0] != '{')
	{
		fprintf(stderr, "level_load_fp: expected '{' after \"4level\"\n");
		return NULL;
	}

	return level_parse_group(fp);
}

box_t *level_load_fname(const char *fname)
{
	FILE *fp = fopen(fname, "rb");

	if(fp == NULL)
	{
		perror("level_load_fname");
		return NULL;
	}

	box_t *b = level_load_fp(fp);
	
	fclose(fp);

	return b;
}


