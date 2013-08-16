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

int net_send_state = 0;
int net_send_spawn = -1;

ENetHost *h_server;
ENetHost *h_client;
ENetPeer *h_toserver;
ENetPeer *h_player[PLAYERS_MAX];

ENetAddress a_server;
ENetAddress a_client;

const char *net_fname_map = NULL;

void net_player_broadcast(player_t *pl, int feedback)
{
	int i;

	for(i = 0; i < PLAYERS_MAX; i++)
	{
		if(h_player[i] == NULL)
			continue;

		if(i == pl->pid && !feedback)
			continue;

		if(pl->magic != 0x4C)
		{
			int lmagic = pl->magic;
			if(i == pl->pid) pl->magic = (pl->magic == 0x66 || pl->magic == 0x69 ? 0x69 : 0xC9);
			else pl->magic = (pl->magic == 0x69 || pl->magic == 0x66 ? 0x66 : 0xC4);
			ENetPacket *p = enet_packet_create(pl, sizeof(player_t), ENET_PACKET_FLAG_RELIABLE);
			pl->magic = lmagic;
			enet_peer_send(h_player[i], 0, p);
		}
	}
}

void net_player_find_spawn(player_t *pl)
{
	// TODO: find a location!
	pl->magic = (pl->magic == 0x69 ? 0xC9 : 0xC4);
}

void net_player_spawn(player_t *pl)
{
	net_player_find_spawn(pl);
	net_player_broadcast(pl, 1);
}

void net_player_hurt(int pid0, int pid1)
{
	if(pid0 >= PLAYERS_MAX)
		return;
	if(pid1 >= PLAYERS_MAX)
		return;
	
	player_t *pl0 = &players[pid0];
	player_t *pl1 = &players[pid1];

	printf("pid %i killed %i\n", pid0, pid1);
	pl1->magic = (pl1->pid == cplr ? 0x69 : 0x66);
	net_player_broadcast(pl1, 1);
}

int net_find_player(ENetPeer *p)
{
	int i;

	for(i = 0; i < PLAYERS_MAX; i++)	
		if(h_player[i] == p)
			return i;
	
	return -1;
}

int net_alloc_player(ENetPeer *h, ENetEvent *e)
{
	int i;

	for(i = 0; i < PLAYERS_MAX; i++)
	{
		if(players[i].magic == 0x4C)
		{
			printf("New player at index %i\n", i);
			players[i].magic = 0xC4;
			h_player[i] = h;
			enet_peer_send(h, 0, enet_packet_create(net_fname_map, strlen(net_fname_map)+1, ENET_PACKET_FLAG_RELIABLE));
			net_player_spawn(&(players[i]));
			return i;
		}
	}

	return -1;
}

char *net_new_client(const char *shost, const char *sport, int toourselves)
{
	ENetEvent ev;

	a_client.port = atoi(sport);
	enet_address_set_host(&a_client, shost);

	// you're not REALLY considering running this on the open internet at this stage, are you?
	h_client = enet_host_create(NULL, 1, 1, 0, 0);
	if(h_client == NULL)
		abort();
	
	// connect!
	h_toserver = enet_host_connect(h_client, &a_client, 1, 0x4DF95);

	if(toourselves) return NULL;

	if(h_toserver == NULL || enet_host_service(h_client, &ev, 5000) <= 0 || ev.type != ENET_EVENT_TYPE_CONNECT)
	{
		printf("Couldn't connect!\n");
		enet_peer_reset(h_toserver);
		fflush(stdout);
		exit(1);
	}

	// yeah ok we're in WHO CARES

	// await filename packet
	for(;;)
	{
		ENetEvent ev;
		while(enet_host_service(h_client, &ev, 5000) > 0)
		{
			if(ev.type != ENET_EVENT_TYPE_RECEIVE)
				break;

			// blatantly obvious exploit right there guys
			// didn't i say don't use this on the open internet?
			return strdup((char *)(ev.packet->data));
		}

		printf("Couldn't finish connecting!\n");
		enet_peer_reset(h_toserver);
		fflush(stdout);
		exit(1);
	}

	net_send_spawn = net_send_state = SDL_GetTicks();
}

void net_new_server(const char *sport, const char *fname)
{
	net_fname_map = fname;
	a_server.host = ENET_HOST_ANY;
	a_server.port = atoi(sport);

	h_server = enet_host_create(&a_server, 128, 1, 0, 0);
	if(h_server == NULL)
		abort();

	// connect to ourselves.
	net_new_client("localhost", sport, 1);
}

void net_update_client(void)
{
	int tupd = SDL_GetTicks();
	int dupd = tupd - net_send_state;

	if(dupd >= 0)
	{
		if(cplr != -1)
		{
			ENetPacket *p = enet_packet_create(&players[cplr], sizeof(player_t), ENET_PACKET_FLAG_RELIABLE);
			enet_peer_send(h_toserver, 0, p);
			net_send_state = tupd + 200;
		}
	}

	if(net_send_spawn == -1 && cplr != -1)
	{
		if(players[cplr].magic == 0x66 || players[cplr].magic == 0x69)
		{
			net_send_spawn = SDL_GetTicks() + 3000;
		}
	}

	if(net_send_spawn != -1)
	{
		dupd = tupd - net_send_spawn;
		if(dupd >= 0)
		{
			net_player_find_spawn(&players[cplr]);
			net_send_spawn = -1;
		}
	}

	ENetEvent ev;
	while(enet_host_service(h_client, &ev, 0) > 0)
	switch(ev.type)
	{
		case ENET_EVENT_TYPE_NONE:
			// what the fuck?
			break;
		case ENET_EVENT_TYPE_CONNECT:
			break;
		case ENET_EVENT_TYPE_RECEIVE:
			if(ev.packet->dataLength == sizeof(player_t))
			{
				player_t *pl = (player_t *)(ev.packet->data);
				if(pl->pid < PLAYERS_MAX)
				{
					memcpy(&players[pl->pid], pl, sizeof(player_t));
					if(pl->magic == 0xC9 || pl->magic == 0x69)
						cplr = pl->pid;
				}
			}
			break;
		case ENET_EVENT_TYPE_DISCONNECT:
			break;
	}
}

void net_update_server(void)
{
	ENetEvent ev;
	while(enet_host_service(h_server, &ev, 0) > 0)
	switch(ev.type)
	{
		case ENET_EVENT_TYPE_NONE:
			// what the fuck?
			break;
		case ENET_EVENT_TYPE_CONNECT:
			if(ev.data != 0x4DF95)
			{
				enet_peer_reset(ev.peer);
				break;
			}

			net_alloc_player(ev.peer, &ev);
			break;
		case ENET_EVENT_TYPE_RECEIVE:
			// oh look, another exploit! ;D
			// let's not be COMPLETELY useless.
			if(ev.packet->dataLength == sizeof(player_t))
			{
				player_t *pl = (player_t *)ev.packet->data;
				if(pl->pid < PLAYERS_MAX)
				{
					if(pl->pid != cplr)
						memcpy(&players[pl->pid], pl, sizeof(player_t));
					net_player_broadcast(pl, 0);
				}
			} else if(ev.packet->dataLength == 2) {
				int pid0 = (int)(ev.packet->data[0]);
				int pid1 = (int)(ev.packet->data[1]);
				net_player_hurt(pid0, pid1);
			}
			break;
		case ENET_EVENT_TYPE_DISCONNECT:
			// finally, we're going to leave players lying around.
			break;
	}

	net_update_client();
}

