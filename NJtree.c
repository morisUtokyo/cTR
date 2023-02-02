#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include "ksw2.h"
#include "clusterTR.h"

// -----------------------------------------------------------------------
// We implemented the method of bulding a neighbor joininng tree for a given distance matrix according to the description in:
// https://en.wikipedia.org/wiki/Neighbor_joining
// Recall that "given an additive distance matrix as input, neighbor joining is guaranteed to find the tree whose distances between taxa agree with it," though this condition is rare in practice.
// -----------------------------------------------------------------------
int neighbor_joining(int *a, int n, int new_node, int **arg_dMat){
    // a is the array of active leaves in the NJ tree.
    // n is the number of active leaves in l.
    if(n == 1){
        int root = new_node - 1;
        NJtree[root].d_parent = 0;
        return(root);
    }
    int *sum_distances = malloc(sizeof(int) * n );
    for(int i=0; i<n; i++){
        sum_distances[i] = 0;
        for(int j=0; j<n; j++)
            sum_distances[i] += arg_dMat[a[i]][a[j]];  // To access arg_dMat, use a[i] !
    }
    // Compute auxiliary table Q
    int min_i = 0;
    int min_j = 1;
    int minQ  = (n-2) * arg_dMat[a[min_i]][a[min_j]] - sum_distances[min_i] - sum_distances[min_j];
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            int tmpQ = (n-2) * arg_dMat[a[i]][a[j]] - sum_distances[i] - sum_distances[j];
            auxQ[a[i]][a[j]] = tmpQ;
            auxQ[a[j]][a[i]] = tmpQ;
            if(tmpQ < minQ){
                minQ = tmpQ; min_i = i; min_j = j;
    }}}
#ifdef DUMP_auxQ
    dump_auxQ(a, n);
#endif
    
    // Set the ID of the new node to "new_node"
    for(int i=0; i<n; i++){
        arg_dMat[new_node][a[i]] =
            ( arg_dMat[a[min_i]][a[i]] + arg_dMat[a[min_j]][a[i]] - arg_dMat[a[min_i]][a[min_j]] ) / 2;
        arg_dMat[a[i]][new_node] = arg_dMat[new_node][a[i]];
    }
    // Revise the distance between the new node and each of min_i and min_j
    int delta_min_i, delta_min_j;
    if(n>2){
        delta_min_i = arg_dMat[a[min_i]][a[min_j]] / 2 + (sum_distances[min_i] - sum_distances[min_j]) / (2*(n-2));
        delta_min_j = arg_dMat[a[min_i]][a[min_j]] - delta_min_i;
    }else{
        delta_min_i = arg_dMat[a[min_i]][a[min_j]];
        delta_min_j = arg_dMat[a[min_i]][a[min_j]];
    }
#ifdef DEBUG_NJ1
    printf("minQ = %d\ti = %d\tj = %d are replaced with %d.\td(%d,%d)=%d,\td(%d,%d)=%d.\n",
           minQ, a[min_i], a[min_j], new_node,
           a[min_i], new_node, delta_min_i,
           a[min_j], new_node, delta_min_j);
#endif
    arg_dMat[new_node][a[min_i]] = delta_min_i;
    arg_dMat[a[min_i]][new_node] = delta_min_i;
    arg_dMat[new_node][a[min_j]] = delta_min_j;
    arg_dMat[a[min_j]][new_node] = delta_min_j;
    arg_dMat[new_node][new_node] = 0;
    
    // Put "new node" to the NJ tree.
    NJtree[new_node].left   = a[min_i];
    NJtree[new_node].right  = a[min_j];
    NJtree[new_node].parent = NO_PARENT;  // default setting
    NJtree[new_node].leaf   = INTERNAL_NODE;    // an internal node
    
    // Put the distance to the parent "new_node"
    NJtree[a[min_i]].d_parent = delta_min_i;
    NJtree[a[min_i]].parent   = new_node;
    
    NJtree[a[min_j]].d_parent = delta_min_j;
    NJtree[a[min_j]].parent   = new_node;

    
    // Add the new node to the head of the active leaves.
    // This operation is quite important to implement the NJ method !!
    int *tmp_a = malloc(sizeof(int) * n );
    int j=0;
    tmp_a[j++] = new_node;  // The new node is added at the head.
    for(int i=0; i<n; i++)
        if(i != min_i && i != min_j)    // Remove  min_i and min_j
            tmp_a[j++] = a[i];          // "a" has n-1 active leaves.
    int root_node = neighbor_joining(tmp_a, n-1, new_node+1, arg_dMat);
    free(sum_distances);
    free(tmp_a);
    return(root_node);
}


int generate_NJtree(int *listReads, int numReads, int newNode, int **arg_dMat){
    
    for(int iReadID = 0; iReadID < numReads; iReadID++){
        NJtree[iReadID].left  = TERMINAL_NODE;
        NJtree[iReadID].right = TERMINAL_NODE;
        NJtree[iReadID].leaf  = LEAF;
        NJtree[iReadID].diameter = 0;
        NJtree[iReadID].diameter_sum_len = 0;
    }
    int NJroot = neighbor_joining(listReads, numReads, newNode, arg_dMat);
    return(NJroot);
}

int generate_NJtree_for_centroids(int numCentroids, int **arg_dMat){
    
    int *internalCentroids = malloc(sizeof(int)*numCentroids);
    for(int i=0; i<numCentroids; i++)
        internalCentroids[i]=i;
    int NJroot = generate_NJtree( internalCentroids, numCentroids, (numCentroids+1), dMat);
    free(internalCentroids);
    return(NJroot);
}

int outlier(int iReadID, int num_reads,
            int *individualReads, int *arg_readLen, int **arg_dMat){
    int outlier_flag = 1;
    int sum_len;
    for(int i=0; i<num_reads; i++)
        if(i != iReadID){
            //---------------------------------------------
            sum_len= arg_readLen[individualReads[iReadID]] +
                         arg_readLen[individualReads[i]];
            if( arg_dMat[iReadID][i]
               < ceil(sum_len * (double)MAX_DIFF_RATIO)  )
                outlier_flag = 0;
            //----------------------------------------------
        }
    return(outlier_flag);
}

int generate_NJtree_from_non_outliers(int num_reads, int *individualReads, int *arg_readLen, int **arg_dMat){
    // List non-outliers
    int *non_outliers = malloc(sizeof(int) * num_reads);
    int num_non_outliers = 0;
    for(int iReadID = 0; iReadID < num_reads; iReadID++){
        int readLen = arg_readLen[individualReads[iReadID]];
        if(outlier(iReadID, num_reads, individualReads, arg_readLen, arg_dMat) != 1){  // Remove outliers
            non_outliers[num_non_outliers++] = iReadID;
        }
    }
    if(num_non_outliers == 0)
        return(-1);
#ifdef DUMP_NJtree
    fprintf(stderr, "non outliers = ");
    for(int i=0; i<num_non_outliers; i++)
        fprintf(stderr, "%d ", non_outliers[i]);
    fprintf(stderr, "\nNum of non outliers = %d, read count = %d\n", num_non_outliers, num_reads);
#endif

    int NJroot = generate_NJtree(non_outliers, num_non_outliers, num_reads, arg_dMat);
    // The second last argument, num_reads, is the ID of the new node to be added.
    free(non_outliers);
    return(NJroot);
}

void dump_auxQ(int *a, int n){
    fprintf(stderr, "Dump of aux Q\ninternalReadID");
    for(int i=0; i<n; i++)
        fprintf(stderr, "\t%d", a[i]);
    fprintf(stderr, "\n");
    for(int i=0; i<n; i++){
        fprintf(stderr, "%d", a[i]);
        for(int j=0; j<n; j++){
            fprintf(stderr, "\t%d", auxQ[a[i]][a[j]]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

void printNJtree_sub(int root, int *individualReads, int *arg_readLen, int NumIndents, NJnode *arg_NJtree, char **arg_readIDs){
    if(root != TERMINAL_NODE){
        // Print the left tree
        printNJtree_sub( arg_NJtree[root].left,  individualReads, arg_readLen, NumIndents+1, arg_NJtree, arg_readIDs);
        // Print the internal node
        for(int i=0; i<NumIndents; i++)
            fprintf(fp_analysis, "  ");
        if(arg_NJtree[root].leaf == INTERNAL_NODE)
            fprintf(fp_analysis, "* %d\n", arg_NJtree[root].d_parent);
        else
            fprintf(fp_analysis, "%s (%d, %d)\n", arg_readIDs[individualReads[root]], arg_NJtree[root].d_parent,
                arg_readLen[individualReads[root]]);
        // Print the right tree
        printNJtree_sub( arg_NJtree[root].right, individualReads, arg_readLen, NumIndents+1, arg_NJtree, arg_readIDs);
    }
}
void printNJtree(int root, int *individualReads, int *arg_readLen, NJnode *arg_NJtree, char **arg_readIDs)
{
    printNJtree_sub(root, individualReads, arg_readLen, 0, arg_NJtree, arg_readIDs);
}
