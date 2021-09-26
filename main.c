#include <sys/time.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include "ksw2.h"
#include "clusterTR.h"

void print_error_message(){
    fprintf(stderr, "cTR [-acs] [-f fasta file] [-d output directory] <fasta file name> \n");
    fprintf(stderr, "-f: Input a fasta file name or a list of loci\n");
    fprintf(stderr, "-d: Output the results to the specified directory \n");
    fprintf(stderr, "-a: Output a detailed analysis of groups\n");
    fprintf(stderr, "-c: Output a fasta file of representative reads of groups\n");
}

int main(int argc, char *argv[])
{
    //---------------------------------------------------------
    // Parse the command with arguments
    //---------------------------------------------------------
    char inputFile[1000];
    char directory[1000] = "";
    int print_analysis = 0;
    int print_centroid_fasta = 0;
    int print_CIGAR = 0;
    
    int opt;
    while ((opt = getopt(argc, argv, "f:d:acsm")) != -1) {
        switch(opt){
            case 'f':
                strcpy(inputFile,optarg);
                break;
            case 'd':
                strcat(directory, optarg);
                strcat(directory, "/");
                break;
            case 'a':
                print_analysis = 1;
                break;
            case 'c':
                print_centroid_fasta = 1;
                break;
            default:
                print_error_message();
                exit(EXIT_FAILURE);
        }
    }
    
    // The input file must be a fasta file of reads from a single locus
    comp_repCentroids(inputFile, directory, print_centroid_fasta, print_analysis, print_CIGAR);


    return 0;
}
