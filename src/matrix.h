// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C) 2012 - 2020  Reza Mohammadi                                                   |
//                                                                                                 |
//     This file is part of BDgraph package.                                                       |
//                                                                                                 |
//     BDgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

#ifndef matrix_H
#define matrix_H
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>

#include "util.h"

extern "C" {
	void sub_matrix( double A[], double sub_A[], int sub[], int *p_sub, int *p  );

	void sub_matrix_upper( double A[], double sub_A[], int sub[], int *p_sub, int *p  );

	void sub_row_mins( double A[], double sub_A[], int *sub, int *p );

	void sub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p );

	void sub_cols_mins( double A[], double sub_A[], int *row, int *col, int *p );

	void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p );

	void sub_matrices( double A[], double A11[], double A21[], double A22[], int *row, int *col, int *p );

	void sub_matrices_inv( double A[], double A11_inv[], double A21[], double A22[], int *row, int *col, int *p );
	
	void inverse( double A[], double A_inv[], int *p );

	void inverse_2x2( double B[], double B_inv[] );

	void cholesky( double A[], double U[], int *p );

	void log_determinant( double A[], double *det_A, int *p );
	
	void select_edge( double rates[], int *index_selected_edge, double *sum_rates, int *qp );

	void select_edge_ts( long double rates[], int *index_selected_edge, long double *sum_rates, int *qp );

	void select_multi_edges( double rates[], int index_selected_edges[], int *size_index, double *sum_rates, int *multi_update, int *qp );

				            
// - - - - - - - NEW for Lang codes - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
	// For Hermitian matrix
	void Hsub_row_mins( double A[], double sub_A[], int *sub, int *p );
		  
	void Hsub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p );

	void Hsub_matrices1( double A[], double A12[], double A22[], int *sub, int *p );

	void Hsub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p );
	
	void cinverse_2x2( double r_B[], double i_B[], double r_B_inv[], double i_B_inv[] );
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	
	void scale_free( int *G, int *p );
	
	void transfer_data( int r_data[], int data[], int *n, int *p, int *size_unique_data );
}

#endif
