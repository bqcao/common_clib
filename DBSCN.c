//
//  DBSCN.c
//  Created by Baoqiang Cao on 11/13/18.
//  Copyright © 2018 Baoqiang Cao. All rights reserved.
//
// Reference: https://en.wikipedia.org/wiki/DBSCAN
//
//

#include "DBSCN.h"

//design filo stack, mistakenly called que. -Baoqiang

int que_push(queue_node **que, point A) {
    queue_node *newp = malloc(sizeof(queue_node));
    int status = 1;
    if (newp == NULL){
        //failed to allocate memory
        status = -1;
    } else {
        newp->pt.x = A.x;
        newp->pt.y = A.y;
        newp->pt.assigned = A.assigned;
        newp->pt.cluster = A.cluster;
        newp->next = NULL;
        if (*que == NULL) {
            *que = newp;
        } else {
            newp->next = *que;
            *que = newp;
        }
    }
    return(status);
}

queue_node *que_pop(queue_node **que) {
    //Be sure to free the returned queue_node from que_pop() call.
    queue_node* ret = NULL;
    if (*que != NULL) {
        ret = *que;
        *que = (*que)->next;
        ret->next = NULL;
    }
    return ret;
}

void print_queue_free(queue_node *que) {
    queue_node* cque = que;
    while (cque != NULL) {
        printf("cque.x=%d, cque.y=%d\n", cque->pt.x, cque->pt.y);
        cque = cque->next;
    }
}

void queue_free(queue_node *que) {
    while (que != NULL) {
        queue_node *ret;
        ret = que;
        que = que->next;
        free(ret);
    }
}

int queue_len(queue_node *que) {
    int cnt = 0;
    while (que != NULL) {
        que = que->next;
        cnt++;
    }
    return(cnt);
}

void set_cluster(queue_node *que, int c) {
    while (que != NULL) {
        que->pt.cluster = c;
        que = que->next;
    }
}

void set_assigned(queue_node *que, int a) {
    while (que != NULL) {
        que->pt.assigned = a;
        que = que->next;
    }
}

void link_list_push(link_list_node **head, int v) {
    link_list_node *nd = malloc(sizeof(link_list_node));
    nd->value = v;
    nd->next = NULL;
    if (*head == NULL) {
        nd->len = 1;
        *head = nd;
    } else {
        nd->next = *head;
        nd->len = (*head)->len+1;
        *head = nd;
    }
}

int link_list_pop(link_list_node **head) {
    int i = -1; // default is a error return, -1.
    link_list_node *nd;
    if (*head != NULL) {
        i = (*head)->value;
        nd = *head;
        *head = (*head)->next;
        free(nd);
    }
    return(i);
}

int link_list_len(link_list_node *head) {
    int cnt = 0;
    while (head != NULL) {
        head = head->next;
        cnt++;
    }
    return cnt;
}

void link_list_free(link_list_node **head) {
    link_list_node *nd;
    if (*head != NULL) {
        nd = *head;
        *head = (*head)->next;
        free(nd);
    }
}

void set_cluster_w_link_list(point *p, link_list_node *lhead, int cluster, int assigned) {
    while (lhead != NULL) {
        (p + (lhead->value))->cluster = cluster;
        if ( (p + (lhead->value))->assigned == -1) {
            // those already assigned, don't bother to re-assign, this is particulary for those core points
            // which already are assigned to 2.
            (p + (lhead->value))->assigned = assigned;
        }
        lhead = lhead->next;
    }
}

void link_list_print_cluster(link_list_node *lhead, point *p) {
    int cnt = 0;
    while (lhead != NULL) {
        printf("%d: index=%d, x=%d, y=%d, assigned=%d, cluster=%d\n", cnt, lhead->value, (p+lhead->value)->x, (p+lhead->value)->y,(p+lhead->value)->assigned,(p+lhead->value)->cluster);
        lhead = lhead->next;
        cnt++;
    }
}


double distance(point A, point B) {
    return sqrt(pow(A.x - B.x, 2.0) + pow(A.y - B.y, 2.0));
}

void init_visit(point *p, int num_p) {
    for (int i = 0; i < num_p; i++) {
        (p+i)->assigned = -1;
        (p+i)->cluster = -1;
    }
}

/*
 *point *p storing all the input points
 num_p, is the count of the input points
 min_c, minimum number of points within a cluster
 eps, the radius of which all points within are defined as a cluster
 return:
        number of cluster
        the cluster assignment are done in *p,
            p->assigned = -1 means, not a cluster member
            p->assigned = 0 means it's a member of a cluster, but it is small island, its members are less than min_c
            p->assigned = 1 means it'a member of a legiment cluster, it's at the outer margin of the cluster
            p->assigned = 2 means it'a member of a legiment cluster, it's at the core center of the cluster
 */
int dbscan (point *p, int num_p, int min_c, double eps) {
    int clusters = 0;
    int i,j;
    init_visit(p, num_p);
    
    link_list_node *link_same_cluster = NULL;
    link_list_node *link_NN = NULL;
    for (i = 0; i < num_p; i++) {
        if ((p+i)->assigned == -1) { //not visited yet
            link_list_push(&link_same_cluster, i);
            (p+i)->assigned = 0;
            int cnt_neighbors = 0;
            for (j = 0; j < num_p; j++) {
                if ( (p+j)->assigned == -1 && distance(*(p+i), *(p+j)) < eps ) {
                    link_list_push(&link_same_cluster, j);
                    link_list_push(&link_NN, j);
                    (p+j)->assigned = 0;
                    cnt_neighbors++;
                }
            }
            if (cnt_neighbors > min_c) {
                //core points of a cluster
                (p+i)->assigned = 2;
                cnt_neighbors = 0;
            }
            
            while (link_NN != NULL) {
                int qi = link_list_pop(&link_NN);
                if (qi > 0) {
                    cnt_neighbors = 0;
                    for (j = 0; j < num_p; j++) {
                        if ( (p+j)->assigned == -1 && distance(*(p+qi), *(p+j)) < eps ) {
                            link_list_push(&link_same_cluster, j);
                            link_list_push(&link_NN, j);
                            (p+j)->assigned = 0;
                            cnt_neighbors++;
                        }
                    }
                    if (cnt_neighbors > min_c) {
                        //core points of a cluster
                        (p+qi)->assigned = 2;
                        cnt_neighbors = 0;
                    }
                }
            }
            
            if (link_list_len(link_same_cluster) >= min_c) {
                set_cluster_w_link_list(p, link_same_cluster, clusters, 1); //fully connected in clusters
                clusters++;
            } else {
                set_cluster_w_link_list(p, link_same_cluster, -1, 0); // small islands
            }
            //print
            //link_list_print_cluster(link_same_cluster, p);
            
            link_list_free(&link_same_cluster);
            link_same_cluster = NULL;
        }
    }
    return(clusters);
}
