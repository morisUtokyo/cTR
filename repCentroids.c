#include <sys/time.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include "ksw2.h"
#include "clusterTR.h"

void free_global_variables_and_exit(){
    // If any of the above global variables failed to be allocated in the heap, free other variables and exit.
    for(int i=0; i < MAX_NUMBER_READS; i++)
        if( reads[i] != NULL ){ free(reads[i]); }
    if(reads != NULL){ free(reads); }
    
    if(readLen != NULL){ free(readLen); }
    
    for(int i=0; i < MAX_NUMBER_READS; i++)
        if( readIDs[i] != NULL ){ free(readIDs[i]); }
    if(readIDs != NULL){ free(readIDs); }

    for(int i=0; i < MAX_NUMBER_READS; i++)
        if( dMat[i] != NULL ){ free(dMat[i]); }
    if(dMat != NULL){ free(dMat); }
    
    for(int i=0; i < MAX_NUMBER_READS; i++)
        if( auxQ[i] != NULL ){ free(auxQ[i]); }
    if(auxQ != NULL){ free(auxQ); }
    
    if(NJtree != NULL){ free(NJtree); }
    
    if(repCentroids != NULL){ free(repCentroids); }

    fprintf(stderr, "cannot allocate space for one of global variables in the heap.\n");
    exit(EXIT_FAILURE);
}

void malloc_global_variables(){
    // Allocate the main memory for global variables in the heap
    reads = malloc(sizeof(int *) * MAX_NUMBER_READS);
    if( reads == NULL ){ free_global_variables_and_exit(); }
    for(int i=0; i < MAX_NUMBER_READS; i++){
        reads[i] = malloc(sizeof(char) * (MAX_READ_LENGTH+1));
        if( reads[i] == NULL ){ free_global_variables_and_exit(); }
    }
    
    readLen = malloc(sizeof(int) * MAX_NUMBER_READS);
    if( readLen == NULL ){ free_global_variables_and_exit(); }
    
    readIDs = malloc(sizeof(int *) * MAX_NUMBER_READS);
    if( readIDs == NULL ){ free_global_variables_and_exit(); }
    for(int i=0; i < MAX_NUMBER_READS; i++){
        readIDs[i] = malloc(sizeof(char) * BLK);
        if( readIDs[i] == NULL ){ free_global_variables_and_exit(); }
    }

    dMat = malloc(sizeof(int *) * MAX_NUMBER_READS);
    if( reads == NULL ){ free_global_variables_and_exit(); }
    for(int i=0; i < MAX_NUMBER_READS; i++){
        dMat[i] = malloc(sizeof(int) * MAX_NUMBER_READS);
        if( dMat[i] == NULL ){ free_global_variables_and_exit(); }
    }
    
    auxQ = malloc(sizeof(int *) * MAX_NUMBER_READS);
    if( reads == NULL ){ free_global_variables_and_exit(); }
    for(int i=0; i < MAX_NUMBER_READS; i++){
        auxQ[i] = malloc(sizeof(int) * MAX_NUMBER_READS);
        if( auxQ[i] == NULL ){ free_global_variables_and_exit(); }
    }
    
    NJtree = malloc(sizeof(NJnode) * MAX_NUMBER_READS * 2);
    if( NJtree == NULL ){ free_global_variables_and_exit(); }
    
    repCentroids = malloc(sizeof(oneCentroid) * MAX_NUMBER_INDIVIDUALS);
    numRepCentroids = 0;
    if(repCentroids == NULL){ free_global_variables_and_exit(); }
    
#ifdef DEBUG_feed
    fprintf(stderr, "Succeeded in malooc (handle_one_file.c)\n");
#endif
    
}

void free_global_variables(){
    for(int i=0; i < MAX_NUMBER_READS; i++)
        if(reads[i] != NULL) free(reads[i]);
    if(reads != NULL)   free(reads);
    if(readLen != NULL) free(readLen);
    for(int i=0; i < MAX_NUMBER_READS; i++)
        if(readIDs[i] != NULL) free(readIDs[i]);
    if(readIDs != NULL) free(readIDs);

    for(int i=0; i < MAX_NUMBER_READS; i++)
        if(dMat[i] != NULL) free(dMat[i]);
    if(dMat != NULL)    free(dMat);
    for(int i=0; i < MAX_NUMBER_READS; i++)
        if(auxQ[i] != NULL) free(auxQ[i]);
    if(auxQ != NULL)    free(auxQ);
    if(NJtree != NULL)  free(NJtree);
    if(repCentroids != NULL) free(repCentroids);
#ifdef DEBUG_feed
    fprintf(stderr, "Succeeded in free (handle_one_file.c)\n");
#endif
}

int same_individual(char *readID1, char *readID2){
    // return 1 if ID1 and ID2 are from the same individual.
    int i=0;
    while(readID1[i] == readID2[i]){
        if(readID1[i] == ',')
            return(1);
        if(readID1[i] == '\0' || readID2[i] == '\0')
            return(0);  // When no comma is found, the individual ID is missing.
        i++;
    }
    return(0);
}

//#define DUMP_DISTANCE
centroidPair cluster_reads_of_an_individual(int num_reads, int *individualReads, char **arg_reads, int *arg_readLen, int print_CIGAR){
    
    centroidPair centPair;
    centPair.numCentroids = 0;
    if(num_reads > 1){
        //---------------------------------------------------------
        // Compute the edit distance between any pair of nodes
        // Use reads[][] and readLen[] to compute dMat[][]
        //---------------------------------------------------------
        int **dMat_individual = compute_edit_distance(num_reads, individualReads, arg_reads, arg_readLen, print_CIGAR);
        
        #ifdef DUMP_DISTANCE
        simple_dump_dMat(num_reads, individualReads, arg_readLen, dMat_individual);
        #endif

        //---------------------------------------------------------
        // Generate the neighbor joining (NJ) tree
        // Use readLen[] and dMat[][] to compute NJtree[]
        //---------------------------------------------------------
        int NJroot = generate_NJtree_from_non_outliers( num_reads, individualReads, arg_readLen, dMat_individual );
        #ifdef DUMP_NJtree
        printNJtree(NJroot, individualReads, arg_readLen, NJtree);
        #endif

        //---------------------------------------------------------
        // Divide into two haplotypes if the NJ tree is not empty.
        // Use dMat[][] and readLen[]
        //---------------------------------------------------------
        if(NJroot != -1){   // The NJ tree is not empty.
            // Group into two haplotypes
            centPair =  cluster_into_two_haplotypes(NJroot, NJroot+1, NJroot, dMat_individual, individualReads, readLen, num_reads);
        }
        
        //---------------------------------------------------------
        // Initialize dMat and NJroot
        //---------------------------------------------------------
        for(int i=0; i<num_reads; i++)
            for(int j=0; j<num_reads; j++)
                dMat[i][j] = 0;
        if(NJroot != -1){
            for(int i=0; i<=NJroot; i++){
                NJtree[i].parent = -1;
                NJtree[i].left   = -1;
                NJtree[i].right  = -1;
                NJtree[i].leaf   = LEAF;
            }
        }
    }
    return(centPair);
}

//#define Select_centroids
int putCentroidsSub(int *centroidList, int numCentroids, oneCentroid Cent, int *individualReads, char **arg_readIDs, int *arg_readLen){
 
    int next_i = numCentroids;
    int deleted = 1;
    if( 0 < Cent.size )
    //if( 1 < Cent.size )
    {
        // repReadID has a local ID (0,1,...) in each individual, but we need to renumber repReadID when we merged reads from multiple individuals.
        centroidList[next_i++] = individualReads[Cent.repReadID];
        deleted = 0;
    }
#ifdef Select_centroids
    if(deleted == 1)
        fprintf(stderr, "--");
    fprintf(stderr, "%s\t%d\t%d\t%d\t%d\n", arg_readIDs[individualReads[Cent.repReadID]], individualReads[Cent.repReadID], Cent.size, Cent.radius, arg_readLen[individualReads[Cent.repReadID]]);
#endif
    return(next_i);
}

int putCentroids(int *centroidList, int numCentroids, centroidPair centPair, int *individualReads, char **arg_readIDs, int *arg_readLen){
    
    int next_i = numCentroids;
    oneCentroid Cent;
    
    // Process the first group
    if(1 <= centPair.numCentroids){
        Cent = centPair.group1;
        next_i = putCentroidsSub(centroidList, next_i, Cent, individualReads, arg_readIDs, arg_readLen);
        //printf("First\t\t%s\t%d\t%d\t%d\n", arg_readIDs[individualReads[Cent.repReadID]], Cent.size, next_i, centPair.numCentroids);
    }
    // If the two haplotypes are homozygous, add the first group to the list as well
    if(1 == centPair.numCentroids){
        Cent = centPair.group1;
        next_i = putCentroidsSub(centroidList, next_i, Cent, individualReads, arg_readIDs, arg_readLen);
        //printf("Second homo\t%s\t%d\t%d\n", arg_readIDs[individualReads[Cent.repReadID]], Cent.size, next_i);
    }
    // Process the second group when we have two centroids
    if(2 == centPair.numCentroids){
        Cent = centPair.group2;
        next_i = putCentroidsSub(centroidList, next_i, Cent, individualReads, arg_readIDs, arg_readLen);
        //printf("Second hetero\t%s\t%d\t%d\n", arg_readIDs[individualReads[Cent.repReadID]], Cent.size, next_i);
    }
    return(next_i);
}

void print_centroid(FILE *fp, char *header, oneCentroid Cent){
    fprintf(fp, header);
    //fprintf(fp, "%d,%d,%s\n", Cent.size, Cent.readLen, Cent.readName);
    fprintf(fp, "GroupSize = %d, Diameter = %d, RadiusFromCentroid = %d, CentroidReadName = %s, CentroidReadLength = %d\n", Cent.size, Cent.diameter, Cent.radius,  Cent.readName, Cent.readLen );
}

int compCentroids(const void* a, const void* b) {
    if (((oneCentroid*)a)->size < ((oneCentroid*)b)->size) {
        return 1;
    } else {
        return -1;
    }
}

int top80p(oneCentroid *repCentroids, int numCentroids){
    int totalSize = 0;
    for(int i=0; i<numCentroids; i++){
        totalSize += repCentroids[i].size;
    }
    int tmpSize = 0;
    for(int i=0; i<numCentroids; i++){
        tmpSize += repCentroids[i].size;
        if((double)0.8 * totalSize < tmpSize)
            return(i+1);
    }
    return(numCentroids);
}

void comp_repCentroids(char *fastaFileName, char *inputDirectory, char *outputDirectory, int print_centroid_fasta, int print_analysis, int print_CIGAR){

    // Specify the output files
    char inputFile[1000]="";
    strcat(inputFile, inputDirectory);
    strcat(inputFile, fastaFileName);
    strcat(inputFile, ".fasta");
    
    if(print_centroid_fasta == 1){
        char file_repCentroids[100] = "";
        strcat(file_repCentroids, outputDirectory);
        strcat(file_repCentroids, fastaFileName);
        strcat(file_repCentroids, "_rep.fasta");
        fp_repCentroids = fopen(file_repCentroids, "w");
    }
    if(print_analysis == 1){
        char file_analysis[100] = "";
        strcat(file_analysis, outputDirectory);
        strcat(file_analysis, fastaFileName);
        strcat(file_analysis, "_analysis.txt");
        fp_analysis     = fopen(file_analysis, "w");
    }
    char file_table[100] = "";
    strcat(file_table, outputDirectory);
    strcat(file_table, fastaFileName);
    strcat(file_table, "_table.txt");
    fp_table     = fopen(file_table, "w");
    
    //---------------------------------------------------------
    // Feed reads from the input fasta file
    // Initialize reads[][], readIDs[][], and readLen[]
    //---------------------------------------------------------
    
    struct timeval s, e;
    float time_feed_reads, time_haplotyping, time_clustering;
    gettimeofday(&s, NULL);
    
    malloc_global_variables();
    
    int read_cnt = handle_one_file(inputFile);
    
    if(print_analysis == 1){
        fprintf(fp_analysis, "The input file name is %s.\nThe read count is %d.\n", inputFile, read_cnt);
        fflush(fp_analysis);
    }
    gettimeofday(&e, NULL);
    time_feed_reads = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
    //---------------------------------------------------------
    // Cluster reads from an individual into two haplotypes
    // Reads from an individual are put to readsInd and readLenInd
    //---------------------------------------------------------
    gettimeofday(&s, NULL);
    
    int focal = 0;
    int numReadsIndividual = 0;
    int individualReads[MAX_NUMBER_READS_FROM_AN_INDIVIDUAL];
    
    centroidPair centPair;
    int centroidList[MAX_NUMBER_INDIVIDUALS];
    int numCentroids = 0;
    
    for(int i=0; i<read_cnt; i++){
        if(same_individual(readIDs[focal], readIDs[i]) == 0){
            // We hit a read from another individual.
            centPair = cluster_reads_of_an_individual(numReadsIndividual, individualReads, reads, readLen, print_CIGAR);
            numCentroids = putCentroids(centroidList, numCentroids, centPair, individualReads, readIDs, readLen);
            focal = i;
            numReadsIndividual = 0;
        }
        // Even when focal == 0,  this works properly.
        individualReads[i-focal] = i;
        numReadsIndividual++;
    }
    centPair = cluster_reads_of_an_individual(numReadsIndividual, individualReads, reads, readLen, print_CIGAR);
    numCentroids = putCentroids(centroidList, numCentroids, centPair, individualReads, readIDs, readLen);
    
    gettimeofday(&e, NULL);
    time_haplotyping = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;

    //---------------------------------------------------------
    // centroidList[] has centroids from individuals.
    // numCentroids is the number of centroids.
    //---------------------------------------------------------
    gettimeofday(&s, NULL);
    
    // Assume that multiple centroids from individual haplotypes.
    if(1 < numCentroids){
        // Output representatives of centroids
        int **dMat_centroids = compute_edit_distance(numCentroids, centroidList, reads, readLen, print_CIGAR);
        #ifdef DUMP_DISTANCE
        simple_dump_dMat(numCentroids, centroidList, readLen, dMat_centroids);
        #endif
        int NJroot = generate_NJtree_for_centroids( numCentroids, dMat_centroids);
        clustering_from_NJtree(NJroot, centroidList, readLen, NJtree, readIDs, reads, dMat_centroids);
        //printNJtree(NJroot, centroidList, readLen, NJtree, readIDs);
        
        // Sort repCentroids
        qsort(repCentroids, numRepCentroids, sizeof(oneCentroid), compCentroids);
        
        // Print a fasta file of the reads for group representatives
        if(print_centroid_fasta == 1){
            for(int i=0; i<numRepCentroids; i++){
                oneCentroid Cent = repCentroids[i];
                print_centroid(fp_repCentroids, "> ", Cent);
                fprintf(fp_repCentroids, "%s\n", reads[Cent.repReadID]);
            }
        }
        // Print an analysis of groups
        if(print_analysis == 1){
            // Print the statistics and members of each cluster
            for(int i=0; i<numRepCentroids; i++){
                oneCentroid Cent = repCentroids[i];
                print_centroid(fp_analysis, "\n", Cent);
                for(int i=0; i<Cent.size; i++){
                    fprintf(fp_analysis, "(%s,%d) ", readIDs[Cent.members[i]], readLen[Cent.members[i]]);
                    fprintf(fp_table, "%s\t%s\n", readIDs[Cent.members[i]], Cent.readName);
                }
                fprintf(fp_analysis, "\n");
            }
            // Output the NJ tree of representatives of centroids
            int repCentroidsList[MAX_NUMBER_INDIVIDUALS];
            for(int i=0; i < numRepCentroids; i++)
                repCentroidsList[i] = repCentroids[i].repReadID;
            int **dMat_repCentroids = compute_edit_distance(numRepCentroids, repCentroidsList, reads, readLen, print_CIGAR);
            NJroot = generate_NJtree_for_centroids( numRepCentroids, dMat_repCentroids);
            fprintf(fp_analysis, "\nTh NJ tree of representative centroids\n");
            printNJtree(NJroot, repCentroidsList, readLen, NJtree, readIDs);
        }
    }
    //else{ fprintf(stderr, "No informative centroids are found in %s.\n", inputFile);}
    
    gettimeofday(&e, NULL);
    time_clustering = (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6;
    
    if(print_analysis == 1){
        fprintf(fp_analysis, "\n\nBreakdown of processing time\n");
        fprintf(fp_analysis, "Feed reads\t%f\n",  time_feed_reads);
        fprintf(fp_analysis, "Haplotyping\t%f\n", time_haplotyping);
        fprintf(fp_analysis, "Clustering\t%f\n",  time_clustering);
    }
    
    if(print_analysis == 1)         fclose(fp_analysis);
    if(print_centroid_fasta == 1)   fclose(fp_repCentroids);
    fclose(fp_table);
    
    /*
    // Print a brief statistics of representative centorids
    int num_top80Groups = top80p(repCentroids, numRepCentroids);
    printf( " %d %d", num_top80Groups, numRepCentroids);
    for(int i=0; i<numRepCentroids; i++){
        // Compute the Lempel-Ziv complexity of the rep. read
        int LZC = Lempel_Ziv(reads[repCentroids[i].repReadID]);
        printf( " %d %d %d", repCentroids[i].size, repCentroids[i].readLen, LZC);
    }
    printf("\n");
    */
    free_global_variables();
}
