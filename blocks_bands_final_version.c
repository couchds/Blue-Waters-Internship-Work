//
//  Parallel implementation of the Needleman-Wunsch sequence alignment algorithm.
//
//  Created by Daniel Couch
//

//  * Fixed indexing bugs (tricky to find index of sequence to look at... ended up having processes send their last index to the next.


#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define TAG 10
#define BIG 1000000
#define ROOT 0

int maxAux(int x, int y);

int max(int x, int y, int z);

/* Each process receives the number of rows it will store. */
int assign_num_rows(int my_rank, int num_procs, int m){
    int num_process_rows;
    /* If the number of process evenly divides the number of rows, give each process the same number. */
    if (m % num_procs == 0) {
        num_process_rows = m/num_procs;
    }
    else{
        /* Compute the number of rows left over after equal distribution. */
        int num_left = m % num_procs;
        /* How many processes will get the floor. */
        int smaller_procs = num_procs - num_left;
        /* Highest rank that will get more data (offset by 1). */
        int highest_rank = num_procs - smaller_procs;
        if (my_rank + 1 > highest_rank) {
            num_process_rows = (int)floor(m/num_procs);
        }
        else{
            num_process_rows = (int)floor(m/num_procs) + 1;
        }
    }
    return num_process_rows;
}

/*
 Process creates its part of the matrix.
 */
void create_matrix_rows(int m, int n, char *sequence_1, char *sequence_2, int *my_rows, int num_procs, int my_rank, int num_process_rows, int size_of_block, int match_value, int mismatch_value, int gap_penalty, MPI_Status *status){
    int i, j, k;
    /* Three values to compare. */
    int value_1 = 0;
    int value_2 = 0;
    int value_3 = 0;
    int number_of_blocks = (int) n/size_of_block;
    int counter = 0;
    int index;
    printf("%s\n", sequence_1);
    int sequence_1_index;
    int sequence_2_index;
    /* Size of last block. */
    int last_block_size;
    if (n % size_of_block != 0){
        /* block size doesn't divide number of columns -> last submatrices will have smaller block size. */
        last_block_size = n % size_of_block;
    }
    else{
        last_block_size = size_of_block;
    }
    
    if (my_rank == ROOT){
        for (i = 0; i <= number_of_blocks; i++) {
            if (i != number_of_blocks) {
                for (j = 0; j < size_of_block; j++) {
                    /* Set up the first elements of the submatrix (in the first row) */
                    my_rows[counter] = counter * gap_penalty;
                    counter++;
                }
                /* Then compute the rest of the submatrix (starting at 1 since first band is computed). */
                for (j = 1; j < num_process_rows; j++) {
                    for (k = 0; k < size_of_block; k++) {
                        /* index keeps track of where we are in the process's matrix. */
                        index = (n * j) + (i * size_of_block) + k;
                        sequence_1_index = j - 1;
                        /* Tell the next process where we are in the sequence. */
                        if (i == 0 && j == 1 && k == 0) {
                            int next_process_seq_start = num_process_rows - 1;
                            MPI_Send(&next_process_seq_start, 1, MPI_INT, my_rank + 1, 500, MPI_COMM_WORLD);
                        }
                        sequence_2_index = (i * size_of_block) + k - 1;
                        char sequence_1_char = sequence_1[sequence_1_index];
                        char sequence_2_char = sequence_2[sequence_2_index];
                        /* First value comes from the top. */
                        value_1 = my_rows[index - n] + gap_penalty;
                        /* Second value comes from the left. */
                        value_2 = my_rows[index - 1] + gap_penalty;
                        /* Third value comes from diagonal. */
                        if (sequence_1_char == sequence_2_char) {
                            value_3 = my_rows[index - n - 1] + match_value;
                        }
                        else{
                            value_3 = my_rows[index - n - 1] + mismatch_value;
                        }
                        int high_score = max(value_1, value_2, value_3);
                        my_rows[index] = high_score;
                        /* First column is filled with gaps. */
                        if (i == 0 && k == 0) {
                            my_rows[index] = value_1;
                        }
                    }
                }
            }
            /* Last sumbmatrix. */
            else if(i == number_of_blocks && last_block_size != size_of_block){
                for (j = 0; j < last_block_size; j++) {
                    /* Set up the first elements of the submatrix (in the first row) */
                    my_rows[counter] = counter * gap_penalty;
                    counter++;
                }
                for (j = 1; j < num_process_rows; j++) {
                    for (k = 0; k < last_block_size; k++) {
                        /* index keeps track of where we are in the process's matrix. */
                        index = (n * j) + (i * size_of_block) + k;
                        sequence_1_index = j - 1;
                        sequence_2_index = (i * size_of_block) + k - 1;
                        char sequence_1_char = sequence_1[sequence_1_index];
                        char sequence_2_char = sequence_2[sequence_2_index];
                        /* First value comes from the top. */
                        value_1 = my_rows[index - n] + gap_penalty;
                        /* Second value comes from the left. */
                        value_2 = my_rows[index - 1] + gap_penalty;
                        /* Third value comes from diagonal. */
                        if (sequence_1_char == sequence_2_char) {
                            value_3 = my_rows[index - n - 1] + match_value;
                        }
                        else{
                            value_3 = my_rows[index - n - 1] + mismatch_value;
                        }
                        int high_score = max(value_1, value_2, value_3);
                        my_rows[index] = high_score;
                    }
                }
            }
            if (i == number_of_blocks) {
                MPI_Send(&my_rows[((num_process_rows - 1) * n) + (i * size_of_block)], last_block_size, MPI_INT, my_rank + 1, TAG, MPI_COMM_WORLD);
            }
            else{
                MPI_Send(&my_rows[((num_process_rows - 1) * n) + (i * size_of_block)], size_of_block, MPI_INT, my_rank + 1, TAG, MPI_COMM_WORLD);
            }
        }
    }
    /* Computation by other processes. */
    else{
        int *diagonal_elements = malloc(number_of_blocks * sizeof(int));
        int my_sequence_start;
        /* Find out where to look in the columnar DNA sequence. */
        MPI_Recv(&my_sequence_start, 1, MPI_INT, my_rank - 1, 500, MPI_COMM_WORLD, status);
        /* All processes except last will send data to the next. */
        if (my_rank != num_procs - 1) {
            int next_process_seq_start;
            next_process_seq_start = my_sequence_start + num_process_rows;
            MPI_Send(&next_process_seq_start, 1, MPI_INT, my_rank + 1, 500, MPI_COMM_WORLD);
        }
        for (i = 0; i <= number_of_blocks; i++) {
            /* Each process receives last row of required submatrix. */
            int *received_block = malloc(size_of_block * sizeof(int));
            int *received_last_block = malloc(last_block_size * sizeof(int));
            /* Last block might receive a different size. */
            if (i == number_of_blocks) {
                MPI_Recv(received_last_block, last_block_size, MPI_INT, my_rank - 1, TAG, MPI_COMM_WORLD, status);
                for (j = 0; j < last_block_size; j++) {
                    printf("%d ", received_last_block[j]);
                }
                printf("\n");
            }
            else{
                MPI_Recv(received_block, size_of_block, MPI_INT, my_rank - 1, TAG, MPI_COMM_WORLD, status);
                for (j = 0; j < size_of_block; j++) {
                    printf("%d ", received_block[j]);
                }
                printf("\n");
            }
            /* Last element of the received block will be a diagonal element for the next submatrix. */
            diagonal_elements[i] = received_block[size_of_block - 1];
            int diagonal_element = diagonal_elements[i - 1];
            /* Computation of submatrix. */
            if (i != number_of_blocks) {
                for (j = 0; j < num_process_rows; j++) {
                    for (k = 0; k < size_of_block; k++) {
                        index = (n * j) + (i * size_of_block) + k;
                        /* Index of sequence1 that we're on is determined by adding current row to what we started on. */
                        sequence_1_index = my_sequence_start + j;
                        sequence_2_index = (i * size_of_block) + k - 1;
                        char sequence_1_char = sequence_1[sequence_1_index];
                        char sequence_2_char = sequence_2[sequence_2_index];
                        if (i == 1) {
                            printf("%c\n", sequence_2_char);
                        }
                        /* Compute gaps in first column. */
                        if (i == 0 && k == 0) {
                            if (j == 0){
                                my_rows[index] = received_block[k] + gap_penalty;
                            }
                            else{
                                my_rows[index] = my_rows[index - n] + gap_penalty;
                            }
                        }
                        /* On the first row of the submatrix, the received row is used. */
                        else if (j == 0){
                            /* On the top-left corner of the submatrix, the diagonal element isn't immediately received. Have to use stored diagonal element. */
                            if (k == 0) {
                                /* First value from top */
                                value_1 = received_block[k] + gap_penalty;
                                
                                /* Second value from left */
                                value_2 = my_rows[index - 1] + gap_penalty;
                                
                                if (sequence_1_char == sequence_2_char) {
                                    value_3 = diagonal_element + match_value;
                                }
                                else{
                                    value_3 = diagonal_element + mismatch_value;
                                }
                                
                                my_rows[index] = max(value_1, value_2, value_3);
                            }
                            /* On other entries of first row, */
                            else{
                                value_1 = received_block[k] + gap_penalty;
                                value_2 = my_rows[index - 1] + gap_penalty;
                                if (sequence_1_char == sequence_2_char) {
                                    value_3 = received_block[k - 1] + match_value;
                                }
                                else{
                                    value_3 = received_block[k - 1] + mismatch_value;
                                }
                                my_rows[index] = max(value_1, value_2, value_3);
                            }
                        }
                        else {
                            value_1 = my_rows[index - n] + gap_penalty;
                            value_2 = my_rows[index - 1] + gap_penalty;
                            if (sequence_1_char == sequence_2_char) {
                                value_3 = my_rows[index - n - 1] + match_value;
                            }
                            else{
                                value_3 = my_rows[index - n - 1] + mismatch_value;
                            }
                            my_rows[index] = max(value_1, value_2, value_3);
                        }
                    }
                }
            }
            /* Computation of last submatrix. */
            else if(i == number_of_blocks && last_block_size != size_of_block){
                for (j = 0; j < num_process_rows; j++) {
                    for (k = 0; k < last_block_size; k++) {
                        index = (n * j) + (i * size_of_block) + k;
                        sequence_1_index = my_sequence_start + j;
                        sequence_2_index = (i * size_of_block) + k - 1;
                        char sequence_1_char = sequence_1[sequence_1_index];
                        char sequence_2_char = sequence_2[sequence_2_index];
                        /* On the first row of the submatrix, the received row is used. */
                        if (j == 0){
                            /* On the top-left corner of the submatrix, the diagonal element isn't immediately received. Have to use stored diagonal element. */
                            if (k == 0) {
                                /* First value from top */
                                value_1 = received_last_block[k] + gap_penalty;
                                
                                /* Second value from left */
                                value_2 = my_rows[index - 1] + gap_penalty;
                                
                                if (sequence_1_char == sequence_2_char) {
                                    value_3 = diagonal_element + match_value;
                                }
                                else{
                                    value_3 = diagonal_element + mismatch_value;
                                }
                                
                                my_rows[index] = max(value_1, value_2, value_3);
                            }
                            /* On other entries of first row, */
                            else{
                                value_1 = received_last_block[k] + gap_penalty;
                                value_2 = my_rows[index - 1] + gap_penalty;
                                if (sequence_1_char == sequence_2_char) {
                                    value_3 = received_last_block[k - 1] + match_value;
                                }
                                else{
                                    value_3 = received_last_block[k - 1] + mismatch_value;
                                }
                                my_rows[index] = max(value_1, value_2, value_3);
                            }
                        }
                        else {
                            value_1 = my_rows[index - n] + gap_penalty;
                            value_2 = my_rows[index - 1] + gap_penalty;
                            if (sequence_1_char == sequence_2_char) {
                                value_3 = my_rows[index - n - 1] + match_value;
                            }
                            else{
                                value_3 = my_rows[index - n - 1] + mismatch_value;
                            }
                            my_rows[index] = max(value_1, value_2, value_3);
                        }
                    }
                }
            }
            /* All processes other than last one sends its last rows. */
            if (my_rank != num_procs - 1) {
                if (i == number_of_blocks) {
                    MPI_Send(&my_rows[((num_process_rows - 1) * n) + (i * size_of_block)], last_block_size, MPI_INT, my_rank + 1, TAG, MPI_COMM_WORLD);
                }
                else{
                    MPI_Send(&my_rows[((num_process_rows - 1) * n) + (i * size_of_block)], size_of_block, MPI_INT, my_rank + 1, TAG, MPI_COMM_WORLD);
                }
            }
        }
    }
}

/*  Used by the align function, calls itself recursively to help determine the optimal alignment. The parameter k keeps track of the size of the aligned sequences.  */
void recursive_align(char *seq1, char *seq2, int m, int n, int *matrix, int match_value, int mismatch_value, int space_penalty, char **seq1_alignment, char **seq2_alignment, int i, int j, int k){
    if (i == 0 && j == 0){
        /* Done with alignment! */
    }
    else{
        /* CASE: Movement from above (notice the decrement of i in the recursive call)
         Occurs if the element of the matrix being considered is equal to the element to the left of it plus the space penalty. */
        if (i > 0 && matrix[(i * n) + j] == matrix[((i - 1) * n) + j] + space_penalty){
            (*seq1_alignment)[k] = seq1[i - 1];
            (*seq2_alignment)[k] = '-';
            k++;
            recursive_align(seq1, seq2, m, n, matrix, match_value, mismatch_value, space_penalty, seq1_alignment, seq2_alignment, i - 1, j, k);
        }
        
        /* CASE: Movement from the left (decrement of j)
         Occurs if the element is equal to the element to the top of it plus the space penalty. */
        else if (j > 0 && matrix[(i * n) + j] == matrix[(i * n) + j - 1] + space_penalty) {
            (*seq1_alignment)[k] = '-';
            (*seq2_alignment)[k] = seq2[j - 1];
            k++;
            recursive_align(seq1, seq2, m, n, matrix, match_value, mismatch_value, space_penalty, seq1_alignment, seq2_alignment, i , j - 1, k);
        }
        
        /* CASE: Movement from diagonal (decrement of both i AND j)
         Occurs otherwise. */
        else{
            (*seq1_alignment)[k] = seq1[i - 1];
            (*seq2_alignment)[k] = seq2[j - 1];
            k++;
            recursive_align(seq1, seq2, m, n, matrix, match_value, mismatch_value, space_penalty, seq1_alignment, seq2_alignment, i - 1, j - 1, k);
        }
    }
}

/*  Just calls the recursive alignment function, sets up aligned sequences.  */
void align(char *seq1, char *seq2, int m, int n, int *matrix,int match_value, int mismatch_value, int space_penalty, char **seq1_alignment, char **seq2_alignment){
    
    recursive_align(seq1, seq2, m, n, matrix, match_value, mismatch_value, space_penalty, seq1_alignment, seq2_alignment, m - 1, n - 1, 0);
}

/*   This function is used to get rid of "junk" data in the aligned sequences, i.e. the 'z' placeholders.  */
char **clean_alignment(int sequence1_size, int sequence2_size, char *sequence1_aligned, char *sequence2_aligned){
    char *alignment1 = malloc(BIG);
    char *alignment2 = malloc(BIG);
    char **alignments = malloc(BIG);
    int i;
    /* k used to index the alignments. */
    int k = 0;
    
    /* The sequences are read in the reverse order (backtracking from the final value in the matrix yields sequences. */
    for (i = sequence1_size + sequence2_size; i >= 0; i--) {
        if (sequence1_aligned[i] == 'A' || sequence1_aligned[i] == 'C' || sequence1_aligned[i] == 'G' || sequence1_aligned[i] == 'T'|| sequence1_aligned[i] == '-') {
            alignment1[k] = sequence1_aligned[i];
            k++;
        }
    }
    k = 0;
    for (i = sequence1_size + sequence2_size; i >= 0; i--) {
        if (sequence2_aligned[i] == 'A' || sequence2_aligned[i] == 'C' || sequence2_aligned[i] == 'G' || sequence2_aligned[i] == 'T'|| sequence2_aligned[i] == '-') {
            alignment2[k] = sequence2_aligned[i];
            k++;
        }
    }
    alignments[0] = alignment1;
    alignments[1] = alignment2;
    return alignments;
}

int main(int argc, char * argv[]) {
    /* Initialize variables: */
    /* Parameters are hard coded. */
    int gap_penalty = -2;
    int match_value = 2;
    int mismatch_value = -1;
    int my_rank;
    int num_procs;
    int i, j;
    int sequence_1_length;
    int sequence_2_length;
    char sequence_1[BIG];
    char sequence_2[BIG];
    MPI_Status status;
    int block_size = atoi(argv[1]);
    
    /* Set up MPI. */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    /* Start the timing. */
    double start, finish;
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    
    /* Input (the root process is responsible for I/O) */
    FILE *fp1, *fp2;
    fp1 = fopen("sequence1.txt", "r");
    fp2 = fopen("sequence2.txt", "r");
    
    
    fscanf(fp1, "%s", sequence_1);
    fscanf(fp2, "%s", sequence_2);
    /* Get the length of the sequences. */
    sequence_1_length = 0;
    sequence_2_length = 0;
    
    /* Search for null terminator. */
    while (sequence_1[sequence_1_length] != '\0'){
        sequence_1_length++;
    }
    while (sequence_2[sequence_2_length] != '\0'){
        sequence_2_length++;
    }
    
    fclose(fp1);
    fclose(fp2);
    
    /* m and n are the dimensions (row-column order) of the scoring matrix. */
    int m = sequence_1_length + 1;
    int n = sequence_2_length + 1;
    
    
    /* Assign workload (how many rows each process gets): */
    int num_process_rows = assign_num_rows(my_rank, num_procs, m);
    /* Allocate memory for processes' rows. */
    int *my_rows = (int *)calloc(n * num_process_rows, sizeof(int));
    
    /* Processes calculate the values in their rows. */
    create_matrix_rows(m, n, sequence_1, sequence_2, my_rows, num_procs, my_rank, num_process_rows, block_size, match_value, mismatch_value, gap_penalty, &status);
    
    
    /* The root needs to know how many rows each process has. */
    int *all_num_rows = malloc(num_procs * sizeof(int));
    int *all_rows = malloc(m * n * sizeof(int));
    int *displacements = malloc(num_procs * sizeof(int));
    MPI_Gather(&num_process_rows, 1, MPI_INT, all_num_rows, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(all_num_rows, num_procs, MPI_INT, ROOT, MPI_COMM_WORLD);
    
    int *scoring_matrix = malloc(m * n * sizeof(int));
    int *receives = malloc(num_procs * sizeof(int));
    for (i = 0; i < num_procs; i++) {
        receives[i] = all_num_rows[i] * n;
    }
    int indexer = 0;
    for (i = 0; i < num_procs; i++) {
        displacements[i] = indexer;
        indexer += (n * all_num_rows[i]);
    }
    MPI_Gatherv(my_rows, n * num_process_rows, MPI_INT, scoring_matrix, receives, displacements, MPI_INT, ROOT, MPI_COMM_WORLD);
    
    if (my_rank == ROOT) {
        printf("\n");
        printf("Scoring matrix created.\n");
        for (i = 0; i < m; i++){
            for (j = 0; j < n; j++) {
                printf("%d  ", scoring_matrix[(i * n) + j]);
            }
            printf("\n");
        }
        
        //
        /*  The parallel part of the program is complete (other than timings). All that is left is the sequential part. */
        char *seq1_alignment = malloc(BIG);
        char *seq2_alignment = malloc(BIG);
        /*   Variables to store "raw" aligned sequences (containing placeholder characters):  */
        char *raw_seq1_aligned = malloc(sequence_1_length + sequence_2_length);
        char *raw_seq2_aligned = malloc(sequence_1_length + sequence_2_length);
        
        /*   Create placeholder characters in the string, to be replaced with appropriate nucleotides  */
        for (i = 0; i < sequence_1_length + sequence_2_length; i++) {
            raw_seq1_aligned[i] = 'z';
        }
        for (i = 0; i < sequence_1_length + sequence_2_length; i++) {
            raw_seq2_aligned[i] = 'z';
        }
        
        align(sequence_1, sequence_2, m, n, scoring_matrix, match_value, mismatch_value, gap_penalty, &raw_seq1_aligned, &raw_seq2_aligned);
        
        printf("\n");
        char **aligned_sequences = malloc(BIG);
        aligned_sequences = clean_alignment(sequence_1_length, sequence_2_length, raw_seq1_aligned, raw_seq2_aligned);
        
        printf("%s\n", aligned_sequences[0]);
        printf("%s\n", aligned_sequences[1]);
        
        free(aligned_sequences[0]);
        free(aligned_sequences[1]);
        free(aligned_sequences);
        free(scoring_matrix);
    }
    free(all_rows);
    free(my_rows);
    free(all_num_rows);
    
    /* Finish timing. */
    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    if (my_rank == ROOT) {
        printf("Time: %f\n", finish - start);
        printf("Input size: %d x %d\n", sequence_1_length, sequence_2_length);
        printf("Submatrix size: %d x %d\n", num_process_rows, block_size);
    }
    
    
    /* Finished using MPI. */
    MPI_Finalize();
    
    return 0;
}

/*  Determines and returns maximum b/t two integers.  */
int maxAux(int x, int y){
    return x > y ? x : y;
}

/*  Uses maxAux to determine and return maximum b/t three integers  */
int max(int x, int y, int z){
    return maxAux((maxAux(x, y)), z);
}
