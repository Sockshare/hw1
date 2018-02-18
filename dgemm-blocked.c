/*
 Please include compiler name below (you may also include any other modules you would like to be loaded)
 
 COMPILER= gnu
 
 Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
 
 CC = cc
 OPT = -O3
 CFLAGS = -Wall -std=gnu99 $(OPT)
 MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
 LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
 
 */

const char* dgemm_desc = "Simple blocked dgemm.";
#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 64//Play around with the block size
#endif

#define min(a,b) (((a)<(b))?(a):(b))

/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, double* A, double* B, double* C)
{

    //Make a copy of the matrix here. Only a submatrix used for the multiplication. Both matrix A, B and C.
    //Same submatrix size
    static double subMatrix[BLOCK_SIZE*BLOCK_SIZE];//MAX size needed
    /* For each row i of A */
    for (int i = 0; i < M; ++i)
    /* For each column j of B */
        for (int j = 0; j < N; ++j)
        {
            /* Compute C(i,j) */
            double cij = C[i+j*lda];
            for (int k = 0; k < K; ++k)
                cij += A[i+k*lda] * B[k+j*lda];
            C[i+j*lda] = cij;
        }
}

//Back-up function

static void do_blockFast (int lda, int M, int N, int K, double* A, double* B, double* C)
{

    //Make a copy of the matrix here. Only a submatrix used for the multiplication. Both matrix A, B and C.
    //Same submatrix size
    static double subMatrix[BLOCK_SIZE*BLOCK_SIZE];//MAX size needed
    /* For each row i of A */


    for (int i = 0; i < M; i+=2)

        /* For each column j of B */
        for (int j = 0; j < N; j+=2)
        {
            /* Compute C(i,j) */
            double cij = C[i+j*lda];

            for (int k = 0; k < K; k+=2){

                cij += A[i+k*lda] * B[k+j*lda];
                cij += A[(i+1)+k*lda] * B[k+(j+1)*lda];
            }
            C[i+j*lda] = cij;
        }
}



/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm (int lda, double* A, double* B, double* C)
{
    /* For each block-row of A */
    for (int i = 0; i < lda; i += BLOCK_SIZE)//this matrix needs to be cached
    /* For each block-column of B */
        for (int j = 0; j < lda; j += BLOCK_SIZE)
        /* Accumulate block dgemms into block of C */
            for (int k = 0; k < lda; k += BLOCK_SIZE)
            {
                /* Correct block dimensions if block "goes off edge of" the matrix */
                //If any of these != blockSize, matrix != square
                int M = min (BLOCK_SIZE, lda-i);
                int N = min (BLOCK_SIZE, lda-j);
                int K = min (BLOCK_SIZE, lda-k);

                //if block is a square do:
                //if block is not a square do:
                
                /* Perform individual block dgemm */
                if(M !=BLOCK_SIZE || N !=BLOCK_SIZE  || K !=BLOCK_SIZE ){//test if the matrix is a square
                    //loop unroll inside here if the matrix is a square

                    do_blockFast(lda, M, N, K, A + i + k * lda, B + k + j * lda, C + i + j * lda);
                }
                else {
                    //Unroll matrix if the size != even
                    //eg: 3x3, 5x5, ect.
                    do_block(lda, M, N, K, A + i + k * lda, B + k + j * lda, C + i + j * lda);
                }
            }
}
