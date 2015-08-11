/*
 *  treecompare2.c
 *  
 *
 *  Created by Chris Creevey on Sat Oct 04 2003.
 *  Copyright (c) 2003 - 2008 Chris Creevey. All rights reserved.
 *
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>


/****** Define  ********/

#define FUNDAMENTAL_NUM 1000
#define TREE_LENGTH 1000
#define TAXA_NUM 50
#define NAME_LENGTH 50



#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TEMP
#define TEMP 2
#endif



/******** Structure definitions *************/

struct node {
	

	char level[100]; /* this will be called "bacteria" or Prevotella" etc depending on the name of the level */
	int  **abundance; /* contains a summary of the totoals of all samples for this level */

	struct node *parent;
	struct node *daughter;
	struct node *prev_sibling;
	struct node *next_sibling;
	


	} node_type;
	
	



