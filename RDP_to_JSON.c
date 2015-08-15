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


#define NAME_LENGTH 1000



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
	
/* function definitions */
void print_JSON(struct node * pos, int num_tabs);
struct node* find_name(struct node * pos, char * present_name);
struct node * add_node(struct node * pos, char * present_name);
struct node * make_node(void);
struct node * find_last_daughter(struct node *p);
void memory_error(int);
void free_memory(void);
void deconstruct_tree(struct node * this);


FILE *infile = '\0', *outfile = '\0', *total_file = '\0';
struct node *root = '\0', *pos = '\0', *new_node = '\0', *tmp_pos = '\0';
char c = '\0', *OTU_name = '\0', *present_name = '\0', *present_pid = '\0', **otu_names = '\0', number[100];
int found=FALSE, cutoff=0, i=0, j=0, pid=0, num_otus, *otu_totals = '\0', value=0;

	
int main(int argc, char *argv[])
{


if(argc < 3)
    {
    printf(" the usage of this program is:\n\n\nRDP_to_JSON RDP.file otu_count.file %cutoff\n\n");
    exit(1);
    }

OTU_name = malloc(NAME_LENGTH*sizeof(char));
if(!OTU_name) memory_error(1);

present_name = malloc(NAME_LENGTH*sizeof(char));
if(!present_name) memory_error(2);

present_pid = malloc(NAME_LENGTH*sizeof(char));
if(!present_pid) memory_error(3);


OTU_name[0] = present_name[0] = present_pid[0] = '\0';


/* 1) open the RDP file */
if((infile = fopen(argv[1], "r")) == '\0')		/* check to see if the file is there */
	{								
	printf("Error: Cannot open file %s\n", argv[1]);
	exit(1);
	}
if((outfile = fopen("taxonomy.json", "w")) == '\0')		/* check to see if the file is there */
	{								
	printf("Error: Cannot open file taxonomy.json\n");
	exit(1);
	}
if((total_file = fopen(argv[2], "r")) == '\0')		/* check to see if the file is there */
	{								
	printf("Error: Cannot open file %s\n", argv[2]);
	exit(1);
	}
cutoff=atoi(argv[3]);

num_otus =0;
while(!feof(total_file))
	{
	if(c=getc(total_file) == '\n' || c == '\r') num_otus++;
	}
rewind(total_file);
/*printf("number of OTUs = %d\n", num_otus);*/

otu_totals = malloc(num_otus*sizeof(int));
if(!otu_totals) memory_error(4);


otu_names = malloc(num_otus*sizeof(char*));
if(!otu_names) memory_error(5);

for(i=0; i<num_otus; i++)
	{	
	otu_totals[i] =0;
	otu_names[i] = malloc(1000*sizeof(char));
	if(!otu_names[i]) memory_error(6);

	otu_names[i][0] = '\0';
	}
i=0; 
while(!feof(total_file))
	{
	j=0;
	while(!feof(total_file) && (c = getc(total_file)) != '\t')
		{
		if(!feof(total_file))
			{	
			otu_names[i][j] = c;
			j++;
			}
		}
	if(!feof(total_file))
		{		
		otu_names[i][j]='\0';
		/*printf("%s\t", otu_names[i]);*/
		j=0;
		}
	while(!feof(total_file) && (c = getc(total_file)) != '\n' && c != '\r')
		{	
		if(!feof(total_file))
			{
			number[j]=c;
			j++;
			}
		}
	if(!feof(total_file))
		{
		number[j] = '\0';
		otu_totals[i] = atoi(number);
		/*printf("%d\n", otu_totals[i]);*/
		/*printf("i=%d %s %d\n", i, otu_names, otu_totals);*/
		i++;
		}
	}

fclose(total_file);
/*printf("finished reading total file\n");*/
/* read in the header information to get to the taxonomly lines */
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
		free_memory();
		exit(1);
		}


/* 2) read in the taxonomy lines one by one and create nodes in the taxonomy tree as necessary */
while(!feof(infile))
	{

	/* now read in the lines */
	/* This is in the following format: */
	/* OTU1;+;Root;100%;Bacteria;100%;"Bacteroidetes";99%;"Bacteroidia";90%;"Bacteroidales";90%;"Prevotellaceae";71%;Prevotella;67% */

	/* read in the OTU name */
	OTU_name[0] = '\0'; i=0;
	while(!feof(infile) && (c=getc(infile)) != ';')
		{
		OTU_name[i] = c;
		i++;
		}
	OTU_name[i] = '\0';

	/* read the "+" or "-" and to the "root" level */
	while(!feof(infile) && (c=getc(infile)) != '+' && c != '-');

	c=(getc(infile)); /* read the ";" */

	/* iterate through each level on the assignment traversigin the tree as you go and adding nodes where necessary */

	pos = root;

	while(!feof(infile) && c!='\n' && c!='\r')
		{	
		present_name[0] = '\0';
		pid=0;
		i=0;
		while(!feof(infile) && (c=getc(infile)) != ';' && c != '\n' && c != '\r')
			{
			present_name[i]=c;
			if(c != '"') i++;
			}
		
		present_name[i] = '\0';
		/*printf("present_name = %s\n", present_name);*/
		i=0;
		while(!feof(infile) && (c=getc(infile)) != '%' && c != '\n' && c != '\r')
			{
			present_pid[i]=c;
			
			i++;
			}

		present_pid[i] = '\0';
		pid = atoi(present_pid);
		/*printf("present pid = %d\n", pid);*/
		if(!feof(infile))
			{
			c=getc(infile); /* read the ';' after the '%' */
			}
		if(pid >= cutoff)
			{	
			if(root == '\0') /* if this is the start of the tree */
				{
				/*printf("starting tree\n");*/
				new_node = make_node();
				strcpy(new_node->level, present_name);
				pos = new_node;
				root = pos;
				}
			else 
				{
				/* check to see if this name exists as a child of the present node */
				if(strcmp(pos->level, present_name) != 0) /* if the present node is not the name found (i.e. root) */
					{
					if(pos->daughter != '\0' ) /* if it has daughters */
						{
						
						if((tmp_pos = find_name(pos->daughter, present_name)) == '\0') /* if the name is not present in the daughters */
							{
								/*printf("name doesn't exists as daughter\n");*/
							new_node = make_node();
							pos = find_last_daughter(pos->daughter);
							/*printf("last daughter is %s\n", pos->level);*/
							strcpy(new_node->level, present_name);
							pos->next_sibling = new_node;
							new_node->prev_sibling = pos;
							pos = new_node;
							}
						else
							{
								/*printf("name exists as daughter %s\n", tmp_pos->level);*/
							pos = tmp_pos; /* if the node already exists, then just assign the pos to be pointing to it in the tree */
							}
						}
					else /* this position has no daughters */
						{	
							/*printf("name doesn't have daughters yet... creating first\n");*/
						new_node = make_node();
						strcpy(new_node->level, present_name);
						pos->daughter = new_node;
						new_node->parent = pos;
						pos = new_node;
						}
					}
				else
					{
					/*printf("%s is present node...skipping\n", present_name );*/
					}
				}
			}
		else /* the %ID is below the cutoff */
			{
			while(!feof(infile) && (c=getc(infile)) != '\n' && c != '\r'); /* read to the end of the line skipping everything else */
			}
		}
	/* now add the OTU as a child of the current position */
	if(!feof(infile))
		{
		new_node = make_node();
		strcpy(new_node->level, OTU_name);
		if(pos->daughter != '\0') /* if this already has a daughter OTU */
			{
			pos = find_last_daughter(pos->daughter);
			pos->next_sibling = new_node;
			new_node->prev_sibling = pos;
			}
		else /* if this is the first daughter OTU */
			{	
			pos->daughter = new_node;
			new_node->parent = pos;
			}
		pos = new_node;
		}

	}	
	fclose(infile);
/* 3) print out the tree of taxonomy as a JSON format */
	print_JSON(root, 0);
	fclose(outfile);
/* 4) clen up the memory */



}


/* prints the JSON structure from the inputted taxonomy details, traveersing the tree in a depth-first manner */
void print_JSON(struct node * pos, int num_tabs)
	{
	int i;
	for(i=0; i<num_tabs; i++) fprintf(outfile, "\t");
	fprintf(outfile, "{\n");
	for(i=0; i<num_tabs; i++) fprintf(outfile, "\t");
	fprintf(outfile, "\"name\": \"%s\",", pos->level);
	
	if(pos->daughter != '\0')
		{
		fprintf(outfile, "\n");
		for(i=0; i<num_tabs; i++) fprintf(outfile, "\t");
		fprintf(outfile, "\"children\" : [\n");
		print_JSON(pos->daughter, num_tabs+1);
		for(i=0; i<num_tabs; i++) fprintf(outfile, "\t");
		fprintf(outfile, "]\n");

		}
	else
		{
		value=0;
		for(i=0; i<num_otus; i++)
			{
			if(strcmp(pos->level, otu_names[i]) == 0)
				{	
				value = otu_totals[i];
				i=num_otus;
				}
			}
		fprintf(outfile, " \"size\": %d\n", value);
		}

	if(pos->next_sibling !='\0')
		{
		for(i=0; i<num_tabs; i++) fprintf(outfile, "\t");
		fprintf(outfile, "},\n");
		print_JSON(pos->next_sibling, num_tabs);
		}
	else
		{
		for(i=0; i<num_tabs; i++) fprintf(outfile, "\t");
		fprintf(outfile, "}\n");
		}
	}


/* finds if a node exists in the chilren of the current node that match the name given, if not, returns NULL, otherwise returns the position of the matched node */
struct node* find_name(struct node * pos, char * name)
	{
	struct node *result = '\0';

	if (strcmp(pos->level, name) == 0)
		{
		result=pos;
		/*printf("found match\n");*/
		}
	else
		{	
		if(pos->next_sibling != '\0')
			{
			result = find_name(pos->next_sibling, name);
			}
		}
	return(result);
	}


/* Adds the node with the given name as a child to the provided position*/
struct node * add_node(struct node * pos, char * present_name)
	{
	struct node *new = '\0';

	new = malloc(sizeof(node_type));
	if(!new) memory_error(10);

	new->level[0] = '\0';
	strcpy(new->level, present_name);
	new->abundance = '\0';

	new->parent;

	}

/* This makes the taxon structure when we need it so I don't have to keep typing the assignments all the time */
struct node * make_node(void)
	{
	struct node *new = '\0';

	new = malloc(sizeof(node_type));
	if(!new) memory_error(11);

	new->level[0] = '\0';
	new->abundance = '\0';

	new->parent = '\0';
	new->daughter = '\0';
	new->prev_sibling = '\0';
	new->next_sibling = '\0';


	return(new);

	}



void memory_error(int place)
	{
	printf("out of memory at %d\n", place);
	free_memory();
	exit(1);
	}



struct node * find_last_daughter(struct node *p)
	{
	struct node * res ='\0';
	if(p->next_sibling != '\0')
		res = find_last_daughter(p->next_sibling);
	else
		res = p;

	return(res);
	}

void free_memory(void)
	{
	if(OTU_name != '\0')
		{
		free(OTU_name);
		OTU_name = '\0';
		}
	if(present_name != '\0')
		{
		free(present_name);
		present_name = '\0';
		}
	if(present_pid != '\0')
		{
		free(present_pid);
		present_pid = '\0';
		}
	if(otu_names != '\0')
		{
		for(i=0; i<num_otus; i++)
			{
			if(otu_names[i] != '\0')
				{	
				free(otu_names[i]);
				otu_names[i] = '\0';
				}
			}
		otu_names = '\0';
		}
	if(otu_totals != '\0')
		{
		free(otu_totals);
		otu_totals = '\0';
		}
	if(root!= '\0')
		{
		deconstruct_tree(root);
		root = '\0';
		}
	}

void deconstruct_tree(struct node * this)
	{
	if(this->daughter != '\0')
		deconstruct_tree(this->daughter);
	if(this->next_sibling != '\0')
		deconstruct_tree(this->next_sibling);

	free(this);
	this='\0';
	}




