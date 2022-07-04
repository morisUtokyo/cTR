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

// Key default parameters
#define MAX_DIFF_RATIO  0.01 // 0.03

#define MAX_NUMBER_READS    10000
#define MAX_READ_LENGTH     30000
#define MAX_NUMBER_READS_FROM_AN_INDIVIDUAL    200
#define MAX_NUMBER_INDIVIDUALS 1000

// Internal variables and data structures
#define MATCH_SCORE     1
#define BLK 4096

char **reads;   // all reads
char **readIDs; // all reads
int  *readLen;  // read length for all reads

int  **dMat;    // matric of distances between reads
int  **auxQ;   // auxiliary table Q for neighbor joining

FILE *fp_repCentroids, *fp_analysis, *fp_table;

// Features of the parent in NJnode
#define NO_PARENT   -2
// Features of left and right
#define TERMINAL_NODE    -1
// Features of the leaf in NJnode
#define LEAF    1
#define INTERNAL_NODE 2

typedef struct {
    int parent;
        // NO_PARENT if the node is the root.
    int left, right;
        // TERMINAL_NODE if the left (right) node is not in the parent's cluster.
    int d_parent;   // distance to the parent
    int leaf;       // LEAF or INTERNAL_NODE
    int diameter;   // The diameter, the largest distance between a pair of leaves
    int diameter_sum_len;  // The sum of the lengths of the two end reads that define the diameter
} NJnode;

NJnode *NJtree;  // array of nodes in the neighbor joining tree

// Representative centroids for each individual
typedef struct {
    int repReadID, size, diameter, radius, readLen;
    // repReadID has a local ID (0,1,...) in each individual, but we need to renumber repReadID when we merged reads from multiple individuals.
    char readName[100];
    int  members[MAX_NUMBER_INDIVIDUALS];
} oneCentroid;

oneCentroid *repCentroids;  // Representatives of centroids
int numRepCentroids;


typedef struct {
    int numCentroids; // 1 if homozygous, 2 if heterozygous
    oneCentroid group1, group2;
} centroidPair;


void malloc_global_variables();
void free_global_variables();
void free_global_variables_and_exit();
int handle_one_file(char *inputFile);
int **compute_edit_distance(int read_cnt, int *indReads, char **arg_reads, int *arg_readLen, int print_CIGAR);
void dump_dMat(int *a, int n, int *arg_readLen, int **arg_dMat);
void simple_dump_dMat(int n, int *listReadIDs, int *arg_readLen, int **arg_dMat);

int generate_NJtree_from_non_outliers(int read_cnt, int *individualReads, int *arg_readLen, int **arg_dMat);
int generate_NJtree_for_centroids(int numReads, int **arg_dMat);
void printNJtree(int root, int *individualReads, int *arg_readLen, NJnode *arg_NJtree, char **arg_readIDs);
void clustering_from_NJtree(int NJroot, int *centroidList, int *arg_readLen, NJnode *arg_NJtree, char **arg_readIDs, char **arg_reads, int **dMat_centroid);
centroidPair cluster_into_two_haplotypes(int NJroot, int NumNodes, int rootNode, int **arg_dMat, int *individualReads, int *arg_readLen, int num_reads);

void comp_repCentroids(char *fastaFileName, char *inputDirectory, char *outputDirectory, int print_centroid_fasta, int print_analysis, int print_CIGAR);

int Lempel_Ziv(char *s);


// External functions
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define DIFF(x, y) ((x) > (y) ? ((x) - (y)) : ((y) - (x)))
