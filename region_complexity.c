#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include "clusterTR.h"

// See https://en.wikipedia.org/wiki/Lempel-Ziv_complexity
//#define Debug_Lempel_Ziv
int Lempel_Ziv(char *s){
    int n;
    for(n=0; s[n]!='\0'; n++){}
    
    // Use 1-origin indexing
    char *tmp_s = malloc(sizeof(char)*(n+1));
    for(int j=0; j<n; j++){ tmp_s[j+1] = s[j]; }
    
    int i = 0;
    int LZC = 1;    // Lempel Ziv complexity
    int u = 1;
    int v = 1;
    int vmax = 1;
    while(u+v <= n){
        if(tmp_s[i+v] == tmp_s[u+v])
            v++;
        else{
            if(vmax < v){ vmax = v; }
            i++;
            if(i == u){     // The componet is exhaustive.
                #ifdef Debug_Lempel_Ziv
                fprintf( stderr, "u-1 = %d, u = %d, v = %d, i = %d, LZC = %d\n", u-1, u, v, i, LZC);
                #endif
                LZC++;
                u += vmax;
                v = 1;
                i = 0;
                vmax = v;
            }else{
                v = 1;
            }
        }
    }
    if(1 < v){
        #ifdef Debug_Lempel_Ziv
        fprintf( stderr, "u-1 = %d, u = %d, v = %d, i = %d, LZC = %d\n", u-1, u, v, i, LZC);
        #endif
        LZC++;
    }
    #ifdef Debug_Lempel_Ziv
    fprintf( stderr, "u+v-1 = %d, u = %d, v = %d, i = %d, LZC = %d\n", u+v-1, u, v, i, LZC);
    #endif
    free(tmp_s);
    return( LZC );
    //return( (double)LZC/n );
}




#define DEBUG_mosaic_repeats

void int2stringKmer(int intKmer, int k, char *Kmer){
    int msd; // The most significant digit
    int remainder = intKmer;
    for(int i=0; i<k; i++)
    {
        msd = remainder / ((int)pow(4,k-1-i));
        switch(msd){
            case 0: Kmer[i] = 'A'; break;
            case 1: Kmer[i] = 'C'; break;
            case 2: Kmer[i] = 'G'; break;
            case 3: Kmer[i] = 'T'; break;
            default: fprintf(stderr, "Illigal msd\n"); exit(1);
        }
        remainder = remainder % ((int)pow(4,k-1-i));
    }
    Kmer[k] = '\0';
}

void string2int(char *s, int *intKmers)
{
    int charCode;
    for(int i=0; s[i]!='\0'; i++){
        switch(s[i]){
            case 'A':
            case 'a':
                charCode = 0; break;
            case 'C':
            case 'c':
                charCode = 1; break;
            case 'G':
            case 'g':
                charCode = 2; break;
            case 'T':
            case 't':
                charCode = 3; break;
            default:
                fprintf(stderr, "Invalid character: %c \n", s[i]); exit(1);
        }
        intKmers[i] = charCode;
    }
}

int encodeKmer(int *intKmers, int i, int k){
    int answer = 0;
    for(int j=0; j<k; j++){
        answer = answer * 4 + intKmers[i+j];
    }
    return(answer);
}


void mosaic_repeats_sub(int k, int len, int *intKmers)
{
    int numKmers = (int)pow(4,k);
    int *freqKmers = malloc(sizeof(int)*numKmers);
    for(int i=0; i<numKmers; i++){ freqKmers[i] = 0; }
    for(int i = 0; i < (len-k+1); i++){
        freqKmers[ encodeKmer(intKmers, i, k) ]++;
    }
    
    double p = 1/pow(4,k);
    double avg = (len-k+1)*p;
    double sd = avg*(1-p);
    double threshold = avg + sd * 3;
    if(threshold < 2) threshold = 2;
    
    #ifdef DEBUG_mosaic_repeats
    fprintf(stderr, "\ninput length = %d, k = %d, avg = %f, sd = %f\n", len, k, avg, sd );
    //for(int i=0; i<numKmers; i++){ fprintf(stderr, "%d ", freqKmers[i]); }
    //fprintf(stderr, "\n");
    #endif
    
    char *Kmer= malloc(sizeof(char) * (k+1));
    int sum = 0;
    int uniqueKmers = 0;
    for(int i=0; i<numKmers; i++){
        if(freqKmers[i] >= threshold){
            sum += freqKmers[i];
            uniqueKmers++;
            int2stringKmer(i, k, Kmer);
            #ifdef DEBUG_mosaic_repeats
            fprintf(stderr, "%s\t%d\n", Kmer, freqKmers[i]);
            #endif
        }
    }
    #ifdef DEBUG_mosaic_repeats
    fprintf(stderr, "significant k-mers / k-mers = %d / %d, unique kmers = %d\n", sum, len-k+1, uniqueKmers);
    #endif
    free(freqKmers);
    free(Kmer);
}

void mosaic_repeats(char *s, int max_k){
    //int  intKmers[1000];
    int len;
    for(len=0; s[len]!='\0'; len++){}
    int *intKmers = malloc(sizeof(int)*len);
    string2int(s, intKmers);
    
    int k=0;
    for(double tmp_len = len; 4 < tmp_len; k++){ tmp_len = tmp_len/4; }
    // (k, len) = (0, < 4), (1, < 16), (2, < 64), (3, < 256), (4, < 1024)
    
    for(int tmp_k = k; tmp_k <= max_k; tmp_k++){
        mosaic_repeats_sub(tmp_k, len, intKmers);
    }
    
    free(intKmers);
}

