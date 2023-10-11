// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Code Originally written by Alex Lenkoski copyright 2013                                     |
//                                                                                                 |                                           |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

#ifndef gwish_H
#define gwish_H
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>


extern "C" {
	void util_es_to_A(int *es, int *A, int p);

    void transpose(int p_i, int p_j, double *A, double *At);

    void get_upper_triangle(int p,double *A_full,double *A_tri);

    void set_mat_identity(int p, double *A);

    void invert(int p, double *A, double *A_inv);

    double log_det(int p, double* A);

    void mult_mats(int p_i, int p_k, int p_j, double *A, double *B, double *C);

    void mult_square_mats(int p, double* A,double* B,double* C);

    void chol(int p,double *A);

    void get_cond_matrix(int p,int p_clique,int *clique_ID,
                 int *V_ID,double *Full, double *Result);

    void copy_mats(int p_i, int p_j, double *A, double *B);

    void make_sub_mat_dbl(int p, int p_sub, int *sub, double *A, double *B);

    void make_sub_mat_int(int p, int p_sub, int *sub, int **A, int *B);

    void get_complementary_set(int p, int p_clique, int *clique_ID, int *V_ID);

    int choose(int n, int m);

    double logchoose(int n, int m);

    int combinations_increment(int n, int r, int *c);

    void combinations_init(int n, int r, int *c);

    int test_add_var(int *A, int p, int *var_list, int list_length, int prop_var);

    int is_subclique(int *var_list, int list_size, int *clmat, int *cldims, int p);

    void add_clique(int *var_list_curr, int curr_length, int *clmat, int *cldims, int p);

    void list_can_augment(int *A, int p, int *var_list_curr, int curr_length, int *clmat, int *cldims);

    int get_cliques(int *A, int p, int *clmat, int *cldims);

    void IPF_MLE(int *CliqueMat, int *CliqueDims, int totcliques, 
             double *target, int p, double thresh, int *nonconverge_flag);

    double gwish_nc_complete(int delta, int p, double *D);

    double gwish_logC(int *A, int delta, double *T, int p);

    double gwish_norm_laplace(int p, int *A, int delta, double *D);

    // Murph: again, I don't want the normalizing constants.
    double gwish_exact_posterior(int p, int delta, int n,  double *D_post); // double *D_prior,

    double gwish_calculateLogPosterior(LPGraph graph, double *D_post, int delta,int n,int *nonconverge_flag);

    void gwish_get_psi_from_K(int p, int delta, double *D, double *K, double *Psi);

    double gwish_get_f_Tsq(int p, int *A, double *Psi);

    int FindDecomposableNeighbors(LPGraph graph,int* decNeighbors);
    
}

#endif