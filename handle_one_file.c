/*
 Copyright (c) 2019, Shinichi Morishita
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 The views and conclusions contained in the software and documentation are those
 of the authors and should not be interpreted as representing official policies,
 either expressed or implied, of the FreeBSD Project.
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "clusterTR.h"

int handle_one_file(char *inputFile){
    //---------------------------------------------------------------------------
    // Feed a string from a file, convert the string into a series of integers
    //---------------------------------------------------------------------------
    
    FILE *fp = fopen(inputFile, "r");
    if(fp == NULL){
        fprintf(stderr, "fatal error: cannot open %s\n", inputFile);
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
    
    char s[BLK];
    int i;
    char charCode;
    int cnt=0;
    int read_cnt = 0;
    int firstRead = 1;  // 1 means the first read.
    
    while (fgets(s, BLK, fp) != NULL) { // Feed a string of size BLK from fp into string s
        if( MAX_READ_LENGTH < cnt){
            fprintf(stderr, "fatal error: the length must be at most %i.\nread ID = %s\n", MAX_READ_LENGTH, readIDs[read_cnt]);
            free_global_variables_and_exit();
            exit(EXIT_FAILURE);
        }
        if( MAX_NUMBER_READS < read_cnt){
            fprintf(stderr, "fatal error: the number of reads must be at most %i. read_cnt=%d\n", MAX_NUMBER_READS, read_cnt);
            free_global_variables_and_exit();
            exit(EXIT_FAILURE);
        }
        if(s[0] == '>'){
            if(firstRead == 1){
                firstRead = 0;
            }else{  // Process the previous read if the current read is not the first one.
                reads[read_cnt][cnt] = '\0';
                readLen[read_cnt] = cnt;
                read_cnt++;
            }
            // Feed the header of the read.
            for(i=1; s[i]!='\0' && s[i]!='\n' && i<BLK; i++){
                readIDs[read_cnt][i-1] = s[i];
            }
            readIDs[read_cnt][i-1] = '\0';
            cnt = 0;
            
        }else{
            for(i=0; s[i]!='\0' && s[i]!='\n'; i++){
                switch(s[i]){
                    case 'A':
                    case 'a':
                        charCode = 'A'; break;
                    case 'C':
                    case 'c':
                        charCode = 'C'; break;
                    case 'G':
                    case 'g':
                        charCode = 'G'; break;
                    case 'T':
                    case 't':
                        charCode = 'T'; break;
                    default:
                        fprintf(stderr, "Invalid character: %c \n", s[i]); exit(EXIT_FAILURE);
                }
                reads[read_cnt][cnt] = charCode;
                cnt++;  // Count up here.
            }
        }
    }
    if(firstRead == 1)  // No annotation starting with ">" is found.
        return(0);      // The read count is zero.
    
    // Process the last read.
    reads[read_cnt][cnt] = '\0';
    readLen[read_cnt] = cnt;
    read_cnt++;
#ifdef DEBUG_feed
    fprintf(stderr, "read count = %d\n", read_cnt);
    for(int j=0; j < read_cnt; j++){
        fprintf(stderr, "%s\n%s\n", readIDs[j], reads[j]);
    }
#endif
    fclose(fp);
    return(read_cnt);
}
