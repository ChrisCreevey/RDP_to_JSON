/*
 *  RDP_to_JSON.c
 *  
 *
 *  Created by Chris Creevey on Tuesday Aug 11 2015.
 *  Copyright (c) 2015 Chris Creevey. All rights reserved.
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
	

FILE *infile = '\0';
struct node *root = '\0';
char c = '\0', OTU_name[1000];
int found=FALSE, cutoff=0, i=0, j=0;

	
int main(int argc, char *argv[])
{


if(argc < 3)
    {
    printf(" the usage of this program is:\n\n\nRDP_to_JSON RDP.file %cutoff\n\n");
    exit(1);
    }

/* 1) open the RDP file */
if((infile = fopen(argv[2], "r")) == '\0')		/* check to see if the file is there */
	{								/* Open the source tree file */
	printf("Error: Cannot open file %s\n", argv[2]);
	exit(1);
	}
cufoff=atoi(argv[3]);

/* 2) read in the taxonomy lines one by one and create nodes in the taxonomy tree as necessary */
while(!feof(infile) && !finished)
	{
	found=FALSE;
	/* read in file until blank line is found */
	while ( !feof(infile) && !found) 
		{
		if (c=(getc(infile)) == '\n' || c == '\r')
			{
			if(c=(getc(infile)) == '\n' || c == '\r')
				found=TRUE;
			}
		}
	if (!found)
		{	
		printf("error no taxonomy information found\n");
		exit(1);
		}

	/* now read in the lines */
	/* This iare in the following format: */
	/* OTU1;+;Root;100%;Bacteria;100%;"Bacteroidetes";99%;"Bacteroidia";90%;"Bacteroidales";90%;"Prevotellaceae";71%;Prevotella;67% */

	/* read in the OTU name */
	OTU_name[0] = '\0'; i=0;
	while(!feof(infile) && (c=(getc(infile)) != ';'))
		{
		OTU_name[i] = c;
		i++;
		}
	OTU_name[i] = '\0';

	/* itierate through each level on the assignment traversigin the tree as you go and adding nodes where necessary */


	}


}



