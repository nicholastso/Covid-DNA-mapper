/* DNA sequence mapper:
 *
 * Skeleton code written by Jianzhong Qi, April 2022
 * Edited by: [Add your name and student ID here]
 *
 */

#include <stdio.h>

#define STAGE_NUM_ONE 1						  /* stage numbers */ 
#define STAGE_NUM_TWO 2
#define STAGE_NUM_THREE 3
#define STAGE_NUM_FOUR 4
#define STAGE_NUM_FIVE 5
#define STAGE_HEADER "Stage %d\n==========\n" /* stage header format string */

#define MAX_READ_ID_LENGTH 100				  /* maximum read ID length */
#define MAX_READ_LENGTH 100					  /* maximum read length */
#define MAX_NUM_READS 100					  /* maximum number of reads */
#define MAX_REF_LENGTH 1000					  /* maximum reference DNA length */

typedef char read_id_t[MAX_READ_ID_LENGTH+1]; /* a read ID */
typedef char read_t[MAX_READ_LENGTH+1];		  /* a read */
typedef char score_t[MAX_READ_LENGTH+1];	  /* quality scores of a read */
typedef char ref_t[MAX_REF_LENGTH+1];		  /* a reference DNA sequence */

/****************************************************************/

/* function prototypes */
void take_one_read(read_t one_read, score_t scores_of_one_read);
void print_stage_header(int stage_num);

void stage_one(read_t one_read, score_t scores_of_one_read);
void stage_two(read_t reads[], score_t scores[], int *num_reads);
void stage_three(read_t reads[], score_t scores[], int num_reads);
void stage_four(ref_t ref_sequence);
void stage_five(read_t reads[], score_t scores[], int num_reads, 
	ref_t ref_sequence);

/* add your own function prototypes here */
/****************************************************************/

/* main function controls all the action, do NOT modify this function */
int
main(int argc, char *argv[]) {
	/* to hold all input reads and quality scores */
	read_t reads[MAX_NUM_READS];	
	score_t scores[MAX_NUM_READS];
	/* to hold the number of input reads */
	int num_reads = 0;	
	/* to hold the input reference sequence */
	ref_t ref_sequence;
	
	/* stage 1: process one read */
	stage_one(reads[0], scores[0]); 
	num_reads++;
	
	/* stage 2: process all reads */
	stage_two(reads, scores, &num_reads);
	
	/* stage 3: mask bases with high error probability */ 
	stage_three(reads, scores, num_reads);
	
	/* stage 4: process reference sequence */
	stage_four(ref_sequence);
	
	/* stage 5: map reads to the reference sequence */
	stage_five(reads, scores, num_reads, ref_sequence);
	
	/* all done; take some rest */
	return 0;
}

/* print stage header given stage number */
void 
print_stage_header(int stage_num) {
	printf(STAGE_HEADER, stage_num);
}

/****************************************************************/
/* add your code below */

/* process a read record */
void 
take_one_read(read_t one_read, score_t scores_of_one_read) {
	/* modify this function for stage 2 */
	read_id_t id;
	scanf("%s", id);
	scanf("%s", one_read);
	getchar();
	getchar(); // skip '+' and '\n'
	scanf("%s", scores_of_one_read);
}

/* stage 1: process one read */
void 
stage_one(read_t one_read, score_t scores_of_one_read) {
	/* print stage header */
	print_stage_header(STAGE_NUM_ONE);
	
	/* add code for stage 1 */
	
	/* process first read */
}

/* stage 2: process all reads */
void 
stage_two(read_t reads[], score_t scores[], int *num_reads) {
	/* add code for stage 2 */
}

/* stage 3: mask bases with high error probability */ 
void 
stage_three(read_t reads[], score_t scores[], int num_reads) {
	/* add code for stage 3 */
}

/* stage 4: process reference sequence */
void stage_four(ref_t ref_sequence) {
	/* add code for stage 4 */
}

/* stage 5: map reads to the reference sequence */
void 
stage_five(read_t reads[], score_t scores[], int num_reads, 
	ref_t ref_sequence) {
	/* add code for stage 5  */
}
