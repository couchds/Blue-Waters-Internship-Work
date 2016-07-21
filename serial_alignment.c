//  Daniel Couch
//
//  Serial program that uses the strategy of dynamic programming in order to align two DNA sequences. A scoring matrix is created via the score function, and the optimal alignment is determined using the align function. Finally, the recAlign function determines the optimal alignment of two prefix sequences.
//
//  This is to be compared with the solution to the problem that uses parallel programming, shown later.
//
//
//  Adapted from: http://wofford-ecs.org/DataAndVisualization/GenomicSequence/index.htm
//
//

#include <stdio.h>
#include <stdlib.h>

#define STRING_SIZE 1000000

//////////////////////////////////////////////////////////////////////////////////////////////////////

// List of functions to be used:

void create_scoring_matrix(char *seq1, char *seq2, int match_value, int mismatch_value, int space_penalty, int m, int n, int matrix[m][n]);

void recursive_align(char *seq1, char *seq2, int m, int n, int matrix[m][n], int match_value, int mismatch_value, int space_penalty, char **seq1_alignment, char **seq2_alignment, int i, int j, int k);

void align(char *seq1, char *seq2, int m, int n, int matrix[m][n] ,int match_value, int mismatch_value, int space_penalty, char **seq1_alignment, char **seq2_alignment);

void display_scoring_matrix(int m, int n, int matrix[m][n]);

void display_alignment(char **sequence1_aligned, char **sequence2_aligned, int sequence1_size, int sequence2_size);

int max(int x, int y, int z);

/////////////////////////////////////////////////////////////////////////////////////////////////////


/*  Display the scoring matrix. Unfortunately, it doesn't work well for large matrices.  */
void display_scoring_matrix(int m, int n, int matrix[m][n]){
    int row, column;
    for (row = 0; row < m; row++) {
        for (column = 0; column < n; column++) {
            printf("%d  ", matrix[row][column]);
        }
        printf("\n");
    }
}

/*  Determine maximum between three integers.  */
int max(int a, int b, int c)
{
    int max_value = a;
    
    if (b > max_value) {
        max_value = b;
    }
    
    if (c > max_value) {
        max_value = c;
    }
    
    return  max_value;
}

/*  Based on the characters of strings seq1 and seq2, constructs a scoring matrix  */
void create_scoring_matrix(char *seq1, char *seq2, int match_value, int mismatch_value, int space_penalty, int m, int n, int matrix[m][n]){
    
    int i, j, diagonal_value;
    diagonal_value = 0;
    
    /* First, make the first column consist of the space penalties. */
    for (i = 0; i < m; i++) {
        matrix[i][0] = i * space_penalty;
    }
    /* The first row should also only consist of space penalties. */
    for (j = 0; j < n; j++) {
        matrix[0][j] = j * space_penalty;
    }
    
    /* Next, the matrix is completed. */
    for (i = 1; i < m; i++) {
        for (j = 1; j < n; j++) {
            
            /* Compare each nucleotide b/t the two sequences. */
            if ((seq1[i - 1]) == (seq2[j - 1])) {
                diagonal_value = match_value;
            }
            else{
                diagonal_value = mismatch_value;
            }
            /* Place the appropriate value in the scoring matrix. It should be the maximum of the possible values that can be entered. */
            matrix[i][j] = max(matrix[i - 1][j - 1] + diagonal_value,
                               matrix[i - 1][j] + space_penalty,
                               matrix[i][j - 1] + space_penalty);
        }
    }
}

/*  Used by the align function, calls itself recursively to help determine the optimal alignment. The parameter k keeps track of the size of the aligned sequences.  */
void recursive_align(char *seq1, char *seq2, int m, int n, int matrix[m][n], int match_value, int mismatch_value, int space_penalty, char **seq1_alignment, char **seq2_alignment, int i, int j, int k){
    if (i == 0 && j == 0){
        /* Done with alignment! */
    }
    else{
        
        /* CASE: Movement from above (notice the decrement of i in the recursive call)
         Occurs if the element of the matrix being considered is equal to the element to the left of it plus the space penalty. */
        if (i > 0 && matrix[i][j] == matrix[i - 1][j] + space_penalty){
            (*seq1_alignment)[k] = seq1[i - 1];
            (*seq2_alignment)[k] = '-';
            k++;
            recursive_align(seq1, seq2, m, n, matrix, match_value, mismatch_value, space_penalty, seq1_alignment, seq2_alignment, i - 1, j, k);
        }
        
        /* CASE: Movement from the left (decrement of j)
         Occurs if the element is equal to the element to the top of it plus the space penalty. */
        else if (j > 0 && matrix[i][j] == matrix[i][j - 1] + space_penalty) {
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
void align(char *seq1, char *seq2, int m, int n, int matrix[m][n],int match_value, int mismatch_value, int space_penalty, char **seq1_alignment, char **seq2_alignment){
    recursive_align(seq1, seq2, m, n, matrix, match_value, mismatch_value, space_penalty, seq1_alignment, seq2_alignment, m - 1, n - 1, 0);
}

/*   This function is used to get rid of "junk" data in the aligned sequences, i.e. the 'z' placeholders.  */
char **clean_alignment(int sequence1_size, int sequence2_size, char *sequence1_aligned, char *sequence2_aligned){
    char *alignment1 = malloc(STRING_SIZE);
    char *alignment2 = malloc(STRING_SIZE);
    char **alignments = malloc(STRING_SIZE);
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

int main() {
    int i;
    
    /* Hard coding the parameters. */
    int match_value = 2;
    int mismatch_penalty = -1;
    int gap_penalty = -2;
    
    /* Allocate memory for sequences. */
    char *seq1 = malloc(STRING_SIZE);
    char *seq2 = malloc(STRING_SIZE);
    
    /* File input: */
        FILE *fp1, *fp2;
        fp1 = fopen("sequence1.txt", "r");
        fp2 = fopen("sequence2.txt", "r");
        fscanf(fp1, "%s", seq1);
        fscanf(fp2, "%s", seq2);
    
    /*   Determine the length of the nucleotide sequences (using the string library).   */
    /* EDIT: At first, strlen was used for this... it seems that, with very large strings (~25,000) the function doesn't work so well. */
    int seq1_length = 0;
    int seq2_length = 0;
    
    /* Search for null terminator. */
    while (seq1[seq1_length] != '\0'){
        seq1_length++;
    }
    while (seq2[seq2_length] != '\0'){
        seq2_length++;
    }
    
    /* m and n are the dimensions of the scoring matrix, in row-column order (added 1 takes into account empty row/column with only gaps in each sequence. */
    int m = seq1_length + 1;
    int n = seq2_length + 1;
    
    /*   Variables to store aligned sequences (maximum possible size is seq1_length + seq2_length):  */
    char *seq1_aligned = malloc(seq1_length + seq2_length);
    char *seq2_aligned = malloc(seq1_length + seq2_length);
    
    /*   Variables to store "raw" aligned sequences (containing placeholder characters):  */
    char *raw_seq1_aligned = malloc(STRING_SIZE);
    char *raw_seq2_aligned = malloc(STRING_SIZE);
    
    /*   Create placeholder characters in the string, to be replaced with appropriate nucleotides  */
    for (i = 0; i < seq1_length + seq2_length; i++) {
        raw_seq1_aligned[i] = 'z';
    }
    for (i = 0; i < seq1_length + seq2_length; i++) {
        raw_seq2_aligned[i] = 'z';
    }
    
    /*  Allocate memory for the scoring matrix.  */
    int *scoring_matrix = malloc(m * n * sizeof(int));
    
    create_scoring_matrix(seq1, seq2, match_value, mismatch_penalty, gap_penalty, m, n, scoring_matrix);
    
    /*  Optional, can be removed.  */
    // printf("\nScoring matrix:\n");
    // display_scoring_matrix(m, n, scoring_matrix);
    
    /*   Get the aligned sequences.  */
    // align(seq1, seq2, m, n, scoring_matrix, match_value, mismatch_penalty, gap_penalty, &raw_seq1_aligned, &raw_seq2_aligned);
    
    char **aligned_sequences = malloc(STRING_SIZE);
    /*  Clean up sequences, store them in the array aligned_sequences  */
    // aligned_sequences = clean_alignment(seq1_length, seq2_length, raw_seq1_aligned, raw_seq2_aligned);
    
    /*  Display alignment. (only works well with smaller sequences)  */
    // printf("\nAlignments:\n");
    // printf("%s\n", aligned_sequences[0]);
    // printf("%s\n", aligned_sequences[1]);
    
    /* Free memory that was allocated to the heap.  */
    free(scoring_matrix);
    free(raw_seq1_aligned);
    free(raw_seq2_aligned);
    free(aligned_sequences[0]);
    free(aligned_sequences[1]);
    free(aligned_sequences);
    
    /*  Input size?  */
    printf("\nSequence 1 size: %d\n", seq1_length);
    printf("Sequence 2 size: %d\n", seq2_length);
    
    return 0;
}
