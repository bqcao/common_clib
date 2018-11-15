//
//  DBSCN.h
//  Created by Baoqiang Cao on 11/13/18.
//  Copyright © 2018 Baoqiang Cao. All rights reserved.
//
// Reference: https://en.wikipedia.org/wiki/DBSCAN
//

#ifndef DBSCN_h
#define DBSCN_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct point {
    int x;
    int y;
    int assigned;
    int cluster;
} point;

typedef struct queue_node {
    point pt;
    struct queue_node *next;
} queue_node;

typedef struct link_list_node {
    int len;
    int value;
    struct link_list_node *next;
} link_list_node;

void link_list_push(link_list_node **head, int v);
int link_list_pop(link_list_node **head);
int link_list_len(link_list_node *head);
void link_list_free(link_list_node **head);
void set_cluster_w_link_list(point *p, link_list_node *lhead, int cluster, int assigned);

double distance(point A, point B);
void init_visit(point *p, int num_p);
int dbscan(point *p, int num_p, int min_c, double eps);

#endif /* DBSCN_h */
