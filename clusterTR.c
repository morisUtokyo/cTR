#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include "ksw2.h"
#include "clusterTR.h"

int inside(int node, int *insideNodes, int ith){
    insideNodes[ith] = node;
    int next_i = ith + 1;
    if( NJtree[node].leaf == INTERNAL_NODE ){
        next_i = inside( NJtree[node].left,  insideNodes, next_i );
        next_i = inside( NJtree[node].right, insideNodes, next_i );
    }
    return(next_i);
}

oneCentroid centroid(int *argReads, int numReads, int **arg_dMat){
    int min_radius = 10000000;
    int min_radius_i = 0;
    for(int i=0; i<numReads; i++){
        int radius = 0;
        for(int j=0; j<numReads; j++){
            if( radius < arg_dMat[argReads[i]][argReads[j]] )
            {
                radius = arg_dMat[argReads[i]][argReads[j]];
            }
        }
        if(radius < min_radius){
            min_radius   = radius;
            min_radius_i = i;
        }
    }
#ifdef DUMP_two_haplotypes
    fprintf(stderr, "centroid Read = %d, radius = %d, size of cluster = %d\n%", argReads[min_radius_i], min_radius, numReads);
#endif
    oneCentroid oneCent;
    oneCent.repReadID   = argReads[min_radius_i];
    oneCent.size        = numReads;
    oneCent.radius      = min_radius;
    return(oneCent);
}

int leaves_in_the_group(int subroot, NJnode *arg_NJtree, int *leaves, int numLeaves)
{
    int new_numLeaves = numLeaves;
    if(arg_NJtree[subroot].leaf == LEAF){
        leaves[new_numLeaves++] = subroot;
        return(new_numLeaves);
    }
    // The left subtree
    if(arg_NJtree[subroot].left != TERMINAL_NODE){
        new_numLeaves = leaves_in_the_group(arg_NJtree[subroot].left, arg_NJtree, leaves, new_numLeaves);
    }
    // The right subtree
    if(arg_NJtree[subroot].right != TERMINAL_NODE){
        new_numLeaves = leaves_in_the_group(arg_NJtree[subroot].right, arg_NJtree, leaves, new_numLeaves);
    }
    return(new_numLeaves);
}

void put_one_cluster_centroid(int subroot, NJnode *arg_NJtree, int *centroidList, char **arg_readIDs, char **arg_reads, int *arg_readLen, int **dMat_centroid)
{
    int leaves[MAX_NUMBER_READS];
    int numLeaves = 0;
    numLeaves = leaves_in_the_group(subroot, arg_NJtree, leaves, numLeaves);
    if(0 < numLeaves){
        oneCentroid Cent = centroid(leaves, numLeaves, dMat_centroid);
        Cent.repReadID   = centroidList[Cent.repReadID];
        Cent.diameter    = arg_NJtree[subroot].diameter;
        Cent.readLen     = readLen[Cent.repReadID];
        strcpy(Cent.readName, arg_readIDs[Cent.repReadID]);
        for(int i=0; i<numLeaves; i++)
            Cent.members[i] = centroidList[leaves[i]];
        
        repCentroids[numRepCentroids++] = Cent;
    }
}

void disconnect_subroot_of_NJtree_sub(int subroot, int left_node, int right_node, NJnode *arg_NJtree, int *centroidList, char **arg_readIDs, char **arg_reads, int *arg_readLen, int **dMat_centroid){
    
    int leaves_left[MAX_NUMBER_INDIVIDUALS];
    int numLeaves_left = 0;
    numLeaves_left = leaves_in_the_group(left_node, arg_NJtree, leaves_left, numLeaves_left);
    
    int leaves_right[MAX_NUMBER_INDIVIDUALS];
    int numLeaves_right = 0;
    numLeaves_right = leaves_in_the_group(right_node, arg_NJtree, leaves_right, numLeaves_right);
    
    // Compute the maxmimum distance between the left and right sets
    int tmp_diameter = -1;
    int tmp_diameter_sum_len;
    if(0 < numLeaves_left && 0 < numLeaves_right){
        for(int i=0; i<numLeaves_left; i++){
            for(int j=0; j<numLeaves_right; j++){
                if(tmp_diameter < dMat_centroid[leaves_left[i]][leaves_right[j]])
                {
                    tmp_diameter = dMat_centroid[leaves_left[i]][leaves_right[j]];
                    tmp_diameter_sum_len =
                        arg_readLen[leaves_left[i]] + arg_readLen[leaves_right[j]];
                }
            }
        }
    }
    if(tmp_diameter < arg_NJtree[left_node].diameter){
        tmp_diameter = arg_NJtree[left_node].diameter;
        tmp_diameter_sum_len = arg_NJtree[left_node].diameter_sum_len;
    }
    if(tmp_diameter < arg_NJtree[right_node].diameter){
        tmp_diameter = arg_NJtree[right_node].diameter;
        tmp_diameter_sum_len = arg_NJtree[right_node].diameter_sum_len;
    }
    // If the small diameter requirement is not met
    if(tmp_diameter > (int)ceil(tmp_diameter_sum_len * (double)MAX_DIFF_RATIO) )
    {
        int disconnected_node, other_node;
        // Disconnect the child with a greater distance
        // The distances can be 0 and are equal.
        // In this case, disconnect the subtee with a greater diameter.
        if(arg_NJtree[left_node].d_parent == arg_NJtree[right_node].d_parent){
            if(arg_NJtree[left_node].diameter > arg_NJtree[right_node].diameter){
                disconnected_node   = left_node;
                other_node          = right_node;
                arg_NJtree[subroot].left = TERMINAL_NODE;
            }else{
                disconnected_node   = right_node;
                other_node          = left_node;
                arg_NJtree[subroot].right = TERMINAL_NODE;
            }
        }else if(arg_NJtree[left_node].d_parent > arg_NJtree[right_node].d_parent){
            disconnected_node   = left_node;
            other_node          = right_node;
            arg_NJtree[subroot].left = TERMINAL_NODE;
        }else{
            disconnected_node   = right_node;
            other_node          = left_node;
            arg_NJtree[subroot].right = TERMINAL_NODE;
        }
        // Print the disconnected node
        put_one_cluster_centroid(disconnected_node, arg_NJtree, centroidList, arg_readIDs, arg_reads, arg_readLen, dMat_centroid);
        // Replace the diameter and other relevant values of the focal subroot with those of the non-disconnected node
        arg_NJtree[subroot].d_parent =
            arg_NJtree[subroot].d_parent + arg_NJtree[other_node].d_parent;
        arg_NJtree[subroot].diameter = arg_NJtree[other_node].diameter;
        arg_NJtree[subroot].diameter_sum_len = arg_NJtree[other_node].diameter_sum_len;
    }else{
        arg_NJtree[subroot].diameter = tmp_diameter;
        arg_NJtree[subroot].diameter_sum_len = tmp_diameter_sum_len;
    }
    // If subroot is the representative of a group, compute the centroid.
    if(arg_NJtree[subroot].parent == NO_PARENT)
        put_one_cluster_centroid(subroot, arg_NJtree, centroidList, arg_readIDs, arg_reads, arg_readLen, dMat_centroid);
}

void disconnect_subroot_of_NJtree(int subroot, NJnode *arg_NJtree, int *centroidList, char **arg_readIDs, char **arg_reads, int *arg_readLen, int **dMat_centroid)
{
    int left_node  = arg_NJtree[subroot].left;
    int right_node = arg_NJtree[subroot].right;
    
    if( left_node == TERMINAL_NODE && right_node == TERMINAL_NODE ){
        arg_NJtree[subroot].diameter = 0;
        arg_NJtree[subroot].diameter_sum_len = 2 * arg_readLen[subroot];
        return;
    }
    if( left_node != TERMINAL_NODE ){
        disconnect_subroot_of_NJtree(left_node, arg_NJtree, centroidList, arg_readIDs, arg_reads, arg_readLen, dMat_centroid);
    }
    if( right_node != TERMINAL_NODE ){
        disconnect_subroot_of_NJtree(right_node, arg_NJtree, centroidList, arg_readIDs, arg_reads, arg_readLen, dMat_centroid);
    }
    if( left_node != TERMINAL_NODE && right_node != TERMINAL_NODE ){
        disconnect_subroot_of_NJtree_sub(subroot, left_node, right_node, arg_NJtree, centroidList, arg_readIDs, arg_reads, arg_readLen, dMat_centroid);
    }
}

void clustering_from_NJtree(int NJroot, int *centroidList, int *arg_readLen, NJnode *arg_NJtree, char **arg_readIDs, char **arg_reads, int **dMat_centroid)
{
    #ifdef NJtree_centroids
    printNJtree(NJroot, centroidList, arg_readLen, arg_NJtree, arg_readIDs);
    #endif
    disconnect_subroot_of_NJtree(NJroot, arg_NJtree, centroidList, arg_readIDs, arg_reads, arg_readLen, dMat_centroid);
}

centroidPair one_haplotype(int rootNode, int NumNodes, int **arg_dMat, int *arg_readLen){
    int *Nodes   = malloc(sizeof(int)*NumNodes);
    NumNodes = inside(rootNode, Nodes, 0);
    int diameter = 0;
    int diameter_i = 0;
    int diameter_j = 0;;
    int first = 1;
    int *tmpReads   = malloc(sizeof(int)*NumNodes);
    int num_tmpReads = 0;
    for(int i=0; i<NumNodes; i++){
        if(NJtree[Nodes[i]].leaf == LEAF){
            tmpReads[num_tmpReads++] = Nodes[i];
            if(first == 1){
                diameter_i = i; first = 0;
            }
            for(int j=i+1; j<NumNodes ; j++){
                if(NJtree[Nodes[j]].leaf == LEAF){
                    int tmp_distance = arg_dMat[Nodes[i]][Nodes[j]];
                    if(diameter <= tmp_distance){
                        diameter = tmp_distance;
                        diameter_i = i;
                        diameter_j = j;
                    }
                }
            }
        }
    }
    centroidPair centPair;
    #ifdef DUMP_two_haplotypes
    fprintf(stderr, "\ndiameter = %d, read pair = (%d,%d)\n", diameter, Nodes[diameter_i], Nodes[diameter_j]);
    #endif
    if( diameter <= (int)ceil((arg_readLen[diameter_i] + arg_readLen[diameter_j]) * MAX_DIFF_RATIO) )
    {
        #ifdef DUMP_two_haplotypes
        fprintf(stderr, "No need to divide into two haplotypes\n");
        #endif
        centPair.numCentroids = 1;
        centPair.group1 = centroid(tmpReads, num_tmpReads, arg_dMat);
        free(Nodes);
        free(tmpReads);
        return(centPair);
    }else{
        centPair.numCentroids = 2;
        free(Nodes);
        free(tmpReads);
        return(centPair);
    }
}

int outside(int node, int insideRoot, int *outsideNodes, int ith){
    int next_i;
    if( node != insideRoot ){
        outsideNodes[ith] = node;
        next_i = ith + 1;
    }else
        return(ith);
    if( NJtree[node].leaf == INTERNAL_NODE){
        next_i = outside( NJtree[node].left,  insideRoot, outsideNodes, next_i );
        next_i = outside( NJtree[node].right, insideRoot, outsideNodes, next_i );
    }
    return(next_i);
}

int brother(int node){
    int parent = NJtree[node].parent;
    if( parent == NO_PARENT){
        return(-1); // This is the root and have no brothers.
    }else{
        if( NJtree[parent].left == node )
            return(NJtree[parent].right);
        else
            return(NJtree[parent].left);
    }
}


void inside_outside(int rootNode, int insideRootNode, int *insideNodes, int *outsideNodes){
    int NumInsideNodes  = inside(  insideRootNode, insideNodes,  0);
    insideNodes[NumInsideNodes] = -1; // The end of the nodes inside.
    int NumOutsideNodes;
    int parent = NJtree[insideRootNode].parent;
    if( parent == NO_PARENT){ // The parent is the root.
    //if( NJtree[parent].parent == NO_PARENT){ // The parent is the root.
        int bro = brother(insideRootNode);
        if( bro != -1){
            NumOutsideNodes = outside( bro, insideRootNode, outsideNodes, 0);
            outsideNodes[NumOutsideNodes] = -1; // The end of the nodes outside.
        }
    }else{
        NumOutsideNodes = outside( rootNode, insideRootNode, outsideNodes, 0);
        outsideNodes[NumOutsideNodes] = -1; // The end of the nodes outside.
    }
}

int min_diameter_node(int *Nodes, int NumNodes, int rootNode, int **arg_dMat){
    int *insideNodes = malloc(sizeof(int)*NumNodes);
    int *outsideNodes= malloc(sizeof(int)*NumNodes);
    int min_diameter = 10000000;
    int min_diameter_node = Nodes[0];
    
    for(int k=0; k < NumNodes; k++){
        int insideRootNode = Nodes[k];
        for(int i=0; i < NumNodes; i++){
            insideNodes[i]  = -1;
            outsideNodes[i] = -1;
        }
        inside_outside(rootNode, insideRootNode, insideNodes, outsideNodes);
        // Computer the diameter of insideNodes.
        int inside_diameter = 0;
        for(int i=0; ; i++){
            if(insideNodes[i] == -1)    break;  // Reach the end of nodes.
            if(NJtree[insideNodes[i]].leaf == LEAF)
                for(int j=i+1; ; j++){
                    if(insideNodes[j] == -1)    break;   // Reach the end of nodes.
                    if(NJtree[insideNodes[j]].leaf == LEAF)
                        inside_diameter = MAX(inside_diameter, arg_dMat[insideNodes[i]][insideNodes[j]]);
                }
        }
        // Computer the diameter of outsideNodes.
        int outside_diameter = 0;
        for(int i=0; ; i++){
            if(outsideNodes[i] == -1)    break;  // Reach the end of nodes.
            if(NJtree[outsideNodes[i]].leaf == LEAF)
                for(int j=i+1; ; j++){
                    if(outsideNodes[j] == -1)    break;  // Reach the end of nodes.
                    if(NJtree[outsideNodes[j]].leaf == LEAF)
                        outside_diameter = MAX(outside_diameter, arg_dMat[outsideNodes[i]][outsideNodes[j]]);
                }
        }
        int tmpD = MAX(inside_diameter, outside_diameter);
        if( min_diameter > tmpD ){
            min_diameter = tmpD;
            min_diameter_node = Nodes[k];
        }
    }
    free(insideNodes);
    free(outsideNodes);
#ifdef DUMP_two_haplotypes
    fprintf(stderr, "Max diamter of two haplotypes = %d, Min_diamster_node = %d\n", min_diameter, min_diameter_node);
#endif
    return(min_diameter_node);
}

centroidPair representatives_inside_outside(int insideRootNode, int *insideNodes, int *outsideNodes, int numNodes, int **arg_dMat, int *arg_readLen){
    int numInsideReads  = 0;
    int numOutsideReads = 0;
    int  *insideReads = malloc(sizeof(int) * numNodes);
    int *outsideReads = malloc(sizeof(int) * numNodes);
    
#ifdef DUMP_two_haplotypes
    fprintf(stderr, "inside = ");
#endif
    for(int i=0; ; i++){
        if(insideNodes[i] == -1){
            break;
        }else{
#ifdef DUMP_two_haplotypes
            if(NJtree[insideNodes[i]].leaf != LEAF)
                fprintf(stderr, "*");
            fprintf(stderr, "%d ", insideNodes[i]);
#endif
            if(NJtree[insideNodes[i]].leaf == LEAF)
                insideReads[numInsideReads++] = insideNodes[i];
        }
    }
#ifdef DUMP_two_haplotypes
    fprintf(stderr, "\n");
    fprintf(stderr, "outside = ");
#endif
    for(int i=0; ; i++){
        if(outsideNodes[i] == -1){
            break;
        }else{
#ifdef DUMP_two_haplotypes
            if(NJtree[outsideNodes[i]].leaf != LEAF)
                fprintf(stderr, "*");
            fprintf(stderr, "%d ", outsideNodes[i]);
#endif
            if(NJtree[outsideNodes[i]].leaf == LEAF)
                outsideReads[numOutsideReads++] = outsideNodes[i];
        }
    }
#ifdef DUMP_two_haplotypes
    fprintf(stderr, "\n");
    fprintf(stderr, "inside root = %d\tNum inside reads = %d\tNum outside reads = %d\n", insideRootNode, numInsideReads, numOutsideReads);
    
    int numAllReads = numInsideReads + numOutsideReads;
    int *allReads = malloc(sizeof(int) * numAllReads);
    for(int i=0; i < numInsideReads; i++)
        allReads[i] = insideReads[i];
    for(int i=0; i < numOutsideReads; i++)
        allReads[i + numInsideReads] = outsideReads[i];
    dump_dMat(allReads, numAllReads, arg_readLen, arg_dMat);
    free(allReads);
#endif
    
    centroidPair centPair;
    centPair.numCentroids = 2;
    centPair.group1 = centroid( insideReads, numInsideReads,  arg_dMat);
    centPair.group2 = centroid(outsideReads, numOutsideReads, arg_dMat);
    
    free(insideReads);
    free(outsideReads);
    return(centPair);
}

centroidPair cluster_into_two_haplotypes(int NJroot, int NumNodes, int rootNode, int **arg_dMat, int *individualReads, int *arg_readLen, int num_reads){
    
    // Initialize
    int *Nodes= malloc(sizeof(int)*(NJroot+1));
    for(int i=0; i<=NJroot; i++) Nodes[i] = i;
    
    // Generate a local list of read lengths
    int readLenInd[MAX_NUMBER_READS_FROM_AN_INDIVIDUAL];
    for(int i=0; i<num_reads; i++){
        readLenInd[i] = arg_readLen[individualReads[i]];
    }
    
    centroidPair centPair = one_haplotype(rootNode, NumNodes, arg_dMat, readLenInd);
    
    if( centPair.numCentroids != 1){    // Heterozygous
        int insideRootNode = min_diameter_node(Nodes, NumNodes, rootNode, arg_dMat);
        int *insideNodes = malloc(sizeof(int)*NumNodes);
        int *outsideNodes= malloc(sizeof(int)*NumNodes);
        inside_outside(rootNode, insideRootNode, insideNodes, outsideNodes);
        centPair = representatives_inside_outside(insideRootNode, insideNodes, outsideNodes, NumNodes, arg_dMat, readLenInd);
        
        free(insideNodes);
        free(outsideNodes);
    }
    // Free the nodes
    free(Nodes);
    return(centPair);
}
