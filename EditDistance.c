#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include "ksw2.h"
#include "clusterTR.h"

// We obtained this function from https://github.com/lh3/ksw2/blob/master/ksw2.h and modified.
int align(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape, int print_CIGAR)
{
    int i, a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
    ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
    /*
    // Print the statistics of the alignment
    printf("%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t",
           ez.max_q, ez.max_t,
           ez.mqe, ez.mqe_t,
           ez.mte, ez.mte_q,
           ez.score,
           ez.m_cigar, ez.n_cigar,
           ez.reach_end);
    */
    if(print_CIGAR == 1){
        printf("CIGAR = ");
        for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
            printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
        putchar('\n');
    }
    
    free(ez.cigar); free(ts); free(qs);
    
    return(ez.score);
}

int **compute_edit_distance(int num_reads, int *individualReads, char **arg_reads, int *arg_readLen, int print_CIGAR)
{
    int sc_mch = MATCH_SCORE;
    int sc_mis = (-1) * sc_mch;
    int gapo = 0;
    int gape = 1;
    
    for(int i=0; i<num_reads; i++){
        dMat[i][i] = 0;
        for(int j=i+1; j<num_reads; j++){
            int score, distance;
            score = align(arg_reads[individualReads[i]], arg_reads[individualReads[j]], sc_mch, sc_mis, gapo, gape, print_CIGAR);
            distance = MIN( arg_readLen[individualReads[i]], arg_readLen[individualReads[j]] ) - score/sc_mch ;
            dMat[i][j] = distance;
            dMat[j][i] = distance;
        }
    }
#ifdef USE_BENCHMARK
    // Example given in https://en.wikipedia.org/wiki/Neighbor_joining
    num_reads = 5;
    dMat[0][0] = 0; dMat[0][1] = 5; dMat[0][2] = 9; dMat[0][3] = 9; dMat[0][4] = 8;
    dMat[1][0] = 5; dMat[1][1] = 0; dMat[1][2] = 10;dMat[1][3] = 10;dMat[1][4] = 9;
    dMat[2][0] = 9; dMat[2][1] = 10;dMat[2][2] = 0; dMat[2][3] = 8; dMat[2][4] = 7;
    dMat[3][0] = 9; dMat[3][1] = 10;dMat[3][2] = 8; dMat[3][3] = 0; dMat[3][4] = 3;
    dMat[4][0] = 8; dMat[4][1] = 9; dMat[4][2] = 7; dMat[4][3] = 3; dMat[4][4] = 0;
#endif
    return(dMat);
}


void simple_dump_dMat(int n, int *listReadIDs, int *arg_readLen, int **arg_dMat){
    fprintf(stderr, "Dump of distance matrix\n");
    fprintf(stderr, "readLen");
    for(int i=0; i<n; i++)
        fprintf(stderr, "\t%d", arg_readLen[listReadIDs[i]]);
    fprintf(stderr, "\nreadID");
    for(int i=0; i<n; i++)
        fprintf(stderr, "\t%d", listReadIDs[i]);
    fprintf(stderr, "\n");
    for(int i=0; i<n; i++){
        fprintf(stderr, "%d", listReadIDs[i]);
        for(int j=0; j<n; j++){
            fprintf(stderr, "\t%d", arg_dMat[i][j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}


void dump_dMat(int *a, int n, int *arg_readLen, int **arg_dMat){
    fprintf(stderr, "Dump of distance matrix\nlength");
    for(int i=0; i<n; i++)
        fprintf(stderr,"\t%d", arg_readLen[a[i]]);
    fprintf(stderr, "\n\n");
    fprintf(stderr, "readID");
    for(int i=0; i<n; i++)
        fprintf(stderr, "\t%d", a[i]);
    fprintf(stderr, "\n");
    for(int i=0; i<n; i++){
        fprintf(stderr, "%d", a[i]);
        for(int j=0; j<n; j++){
            fprintf(stderr, "\t%d", arg_dMat[a[i]][a[j]]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}
