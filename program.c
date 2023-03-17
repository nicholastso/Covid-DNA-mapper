#include <stdio.h>
#include <string.h>
#include <math.h>

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
#define MAX_READ_RECORD 100					  /* maximum number of records */
#define SCORE_MIN 33
#define SCORE_MAX 73						  /* valid base quality score range */
#define MAX_ERROR_PROB 0.15					  /* error limit before we mask that base */
#define MASKED_BASE '*'						  /* mask character representation */

typedef char read_id_t[MAX_READ_ID_LENGTH+1]; /* a read ID */
typedef char read_t[MAX_READ_LENGTH+1];		  /* a read */
typedef char score_t[MAX_READ_LENGTH+1];	  /* quality scores of a read */
typedef char ref_t[MAX_REF_LENGTH+1];		  /* a reference DNA sequence */

void process_one_read(read_t one_read, score_t scores_of_one_read);
void print_stage_header(int stage_num);
void stage_one(read_t one_read, score_t scores_of_one_read);
void stage_two(read_t reads[], score_t scores[], int *num_reads);
void stage_three(read_t reads[], score_t scores[], int num_reads);
void stage_four(ref_t ref_sequence);
void stage_five(read_t reads[], score_t scores[], int num_reads, 
	ref_t ref_sequence);
int index_of_base_with_smallest_quality_score(score_t scores_of_one_read);
int smallest_average(score_t scores[], int num_reads);
double average_quality_score(score_t scores_of_one_read);
void replace_char(read_t *one_read, score_t scores_of_one_read);
double error_probability(char score);

int main(int argc, char *argv[]) {

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
	
	return 0;

}

void print_stage_header(int stage_num) {

	printf(STAGE_HEADER, stage_num);

}

/* process a read record */
void process_one_read(read_t one_read, score_t scores_of_one_read) {

	read_id_t id;
    scanf("%s", id);

    if (id[0] == '@') {
        scanf("%s", one_read);
        getchar();
        getchar();   // skip '+' and '\n'
        scanf("%s", scores_of_one_read);
	}
	else if (id[0] == '#') {
        one_read[0] = '\0';
    }

}

/* stage 1: process one read */
void stage_one(read_t one_read, score_t scores_of_one_read) {

	print_stage_header(STAGE_NUM_ONE);
	process_one_read(one_read, scores_of_one_read);

	int index=index_of_base_with_smallest_quality_score(
		scores_of_one_read);

	printf("Base with the samllest quality score: %c\n", one_read[index]);
	printf("Index: %d\n", index);

}

int index_of_base_with_smallest_quality_score(score_t scores) {

	int smallest=0;

	for(int i=1; scores[i]!='\0'; i++) {
		if(scores[i]<scores[smallest]) {
			smallest=i;
		}
	}

	return smallest;

}

/* stage 2: process all reads */
void 
stage_two(read_t reads[], score_t scores[], int *num_reads) {

    printf("\n");
    print_stage_header(STAGE_NUM_TWO);

	for(int i=1; i<MAX_READ_RECORD; i++) {
		process_one_read(reads[i], scores[i]);
		if(reads[i][0]=='\0') {
			break;
		}
        (*num_reads)++;
	}

	int smallest_average_score=smallest_average(scores, *num_reads);

    printf("Total number of reads: %d\n", *num_reads);
    printf("Smallest average quality score: %.2lf\n", average_quality_score(scores[smallest_average_score]));
    printf("Read with the smallest average quality score:\n");
    printf("%s\n", reads[smallest_average_score]);

}

int smallest_average(score_t scores[], int num_reads) {

	int smallest_average=0;
	int smallest_score=average_quality_score(scores[0]);

	for (int i = 1; i < num_reads; i++) {
        double avg_score = average_quality_score(scores[i]);
        if (avg_score < smallest_score) {
            smallest_average = i;
            smallest_score = avg_score;
        }
    }

    return smallest_average;

}

double average_quality_score(score_t scores) {

	double sum=0;
    double count=0;

	for(int i=0; scores[i]!='\0'; i++) {
		sum+=scores[i];
        count++;
	}

	return sum/count;

}

/* stage 3: mask bases with high error probability */ 
void 
stage_three(read_t reads[], score_t scores[], int num_reads) {

    print_stage_header(STAGE_NUM_THREE);

	for(int i=0; i<num_reads; i++) {
		for(int j=0; scores[i][j]!='\0'; j++) {
			if(error_probability(scores[i][j])>MAX_ERROR_PROB) {
				reads[i][j]=MASKED_BASE;
			}
		}
    }

	for(int i=0; i<num_reads; i++) {
		printf("%s\n", reads[i]);
	}
    
	printf("\n");

}

double error_probability(char score) {

    return pow(10,-(score-SCORE_MIN)/10);

}

/* stage 4: process reference sequence */
void stage_four(ref_t ref_sequence) {

	print_stage_header(STAGE_NUM_FOUR);
	scanf("%s", ref_sequence);
	int A=0, C=0, G=0, T=0, length=0;

	for(int i=0; ref_sequence[i]!='\0'; i++) {
		if (ref_sequence[i] == 'A') {
            A++;
        } else if (ref_sequence[i] == 'C') {
            C++;
        } else if (ref_sequence[i] == 'G') {
            G++;
        } else if (ref_sequence[i] == 'T') {
            T++;
        }
		length++;
	}

	printf("Length of the reference sequence: %d\n", length);
    printf("Number of A bases: %d\n", A);
    printf("Number of C bases: %d\n", C);
    printf("Number of G bases: %d\n", G);
    printf("Number of T bases: %d\n", T);
    printf("\n");

}

/* stage 5: map reads to the reference sequence */
void 
stage_five(read_t reads[], score_t scores[], int num_reads, 
	ref_t ref_sequence) {
		
}