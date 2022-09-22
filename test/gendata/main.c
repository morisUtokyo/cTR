#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include "MT.h"     // Use the Mersenne Twister.

#define NUM_MOSAIC_TRS 1000
#define ERROR_RATE 0 // 0.05

#define NUM_READS_IN_ONE_HAPLOTYPE 3
#define MAX_UNIT_OCC 20
#define MIN_UNIT_OCC 10

#define ST_LEN 100000

//#define DEBUG_errors

void print_units(char *st, char *unit, int n){
    for(int i=0; i<n; i++)
        sprintf(st, "%s%s", st, unit); // append the give unit n times.
}

char another_char(char c){
    char c1;
    for(;;){
        switch(genrand_int32()%4){
            case 0: c1='A';
            case 1: c1='C';
            case 2: c1='G';
            case 3: c1='T';
        }
        if(c!=c1) return(c1);
    }
}

int main(int argc, char *argv[])
{
    int num_units, num_reads_in_one_haplotype=NUM_READS_IN_ONE_HAPLOTYPE, max_unit_occ=MAX_UNIT_OCC, min_unit_occ=MIN_UNIT_OCC;
    int num_TRs = NUM_MOSAIC_TRS;
    int opt;
    while ((opt = getopt(argc, argv, "k:l:m:n:h:")) != -1) {
        switch(opt){
            case 'k':
                sscanf(optarg, "%d", &min_unit_occ);
                //fprintf(stderr, "Minimum number of unit occrrences = %d\n", min_unit_occ);
                break;
            case 'l':
                sscanf(optarg, "%d", &max_unit_occ);
                //fprintf(stderr, "Maximum number of unit occrrences = %d\n", max_unit_occ);
                break;
            case 'm':
                sscanf(optarg, "%d", &num_units);
                //fprintf(stderr, "Number of units = %d\n", num_units);
                break;
            case 'n':
                sscanf(optarg, "%d", &num_TRs);
                //fprintf(stderr, "Number of TRs = %d\n", num_TRs);
                break;
            case 'h':
                sscanf(optarg, "%d", &num_reads_in_one_haplotype);
                break;
            default:
                fprintf(stderr, "Usage: gen -f (number of units) units -n (number of TRs) -e (error rate) \n");
                exit(EXIT_FAILURE);
        }
    }
    
    if(max_unit_occ - min_unit_occ < 1){
        fprintf(stderr, "max_unit_occ must be greater than min_unit_occ\n");
        exit(EXIT_FAILURE);
    }

    int *randNums  = (int *) malloc( sizeof(int) * num_units );
    char *st = (char *) malloc( sizeof(char) * ST_LEN );

    // Generate reads for num_TRs individuals
    for(int j=0; j<num_TRs; j++){
        // Generate reads for two haplotypes, separately
        for(int h=0; h<2; h++){
            for (int i=0; i<num_units; i++)
                randNums[i] = (genrand_int32()%(max_unit_occ - min_unit_occ)) + min_unit_occ;
            sprintf(st, "");
            int start_units = argc - num_units; // The first unit is argv[i + start_units]
            for (int i=0; i<num_units; i++)
                print_units(st, argv[i + start_units], randNums[i]);            
            // Generate num_reads_in_one_haplotype reads for each haplotype
            for(int k=0; k<num_reads_in_one_haplotype; k++){
                // j is the individual ID, and k is the read ID.
                printf("> %d,%d,", j, k+h*num_reads_in_one_haplotype);
                for (int i=0; i<num_units; i++) // add mosaic TR pattern.
                    printf("(%s)%d", argv[i + start_units], randNums[i]);
                // Print the string of the mosaic TR.
                printf("\n%s\n", st);
            }
        }
    }

    free(randNums);
    free(st);
    return EXIT_SUCCESS;
}
