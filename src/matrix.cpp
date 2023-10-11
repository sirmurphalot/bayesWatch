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
  
#include "matrix.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves square sub_matrix B (p_sub x p_sub), dictated by vector sub
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrix( double A[], double sub_A[], int sub[], int *p_sub, int *p )
{
	int i, j, ixp, subixp, psub = *p_sub, pdim = *p;
	
	for( i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[ i ] * pdim;
		
		for( j = 0; j < psub; j++ )
			sub_A[ ixp + j ] = A[ subixp + sub[ j ] ]; 
	}
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes symmetric matrix A (p x p) and 
// retrieves upper part of sub_matrix B (p_sub x p_sub), dictated by vector sub
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrix_upper( double A[], double sub_A[], int sub[], int *p_sub, int *p )
{
	int i, j, ixp, subixp, psub = *p_sub, pdim = *p;
			
	for( i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[ i ] * pdim;
			
		for( j = 0; j <= i; j++ )
			sub_A[ ixp + j ] = A[ subixp + sub[ j ] ]; 
	}
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves vector sub_A which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[j, -j] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int subj = *sub, pdim = *p, subxp = subj * pdim;

	memcpy( sub_A       , A + subxp           , sizeof( double ) * subj );		
	memcpy( sub_A + subj, A + subxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(2 x p-2) which is sub rows of matrix A, minus two elements
// Likes A[(i,j), -(i,j)] in R ONLY FOR SYMMETRIC MATRICES
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col, sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		sub_A[ l++ ] = A[ sub0p + i ]; 
		sub_A[ l++ ] = A[ sub1p + i ]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		sub_A[ l++ ] = A[ sub0p + i ]; 
		sub_A[ l++ ] = A[ sub1p + i ]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		sub_A[ l++ ] = A[ sub0p + i ]; 
		sub_A[ l++ ] = A[ sub1p + i ]; 
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(p-2 x 2) which is sub cols of matrix A, minus two elements
// Likes A[-(i,j), (i,j)] in R 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_cols_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int subi = *row, subj = *col, pdim = *p, p2 = pdim - 2, subixp = subi * pdim, subjxp = subj * pdim;

	memcpy( sub_A           , A + subixp           , sizeof( double ) * subi );		
	memcpy( sub_A + subi    , A + subixp + subi + 1, sizeof( double ) * ( subj - subi - 1 ) );	
	memcpy( sub_A + subj - 1, A + subixp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	

	memcpy( sub_A + p2           , A + subjxp           , sizeof( double ) * subi );		
	memcpy( sub_A + p2 + subi    , A + subjxp + subi + 1, sizeof( double ) * ( subj - subi - 1 ) );	
	memcpy( sub_A + p2 + subj - 1, A + subjxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes symmatric matrix A (p x p) and 
// retrieves A12(1x(p-1)) and A22((p-1)x(p-1))
// Like A12=A[j, -j], and A22=A[-j, -j] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, ixpdim, ixp1, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim, mpsub = pdim - psub - 1;
    int size_psub  = sizeof( double ) * psub;
    int size_mpsub = sizeof( double ) * mpsub;

	memcpy( A12,        A + subxp,            size_psub );	
	memcpy( A12 + psub, A + subxp + psub + 1, size_mpsub );	

	for( i = 0; i < psub; i++ )
	{	
		ixpdim = i * pdim;
		ixp1   = i * p1;
		
		memcpy( A22 + ixp1       , A + ixpdim           , size_psub );
		memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, size_mpsub );
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		ixpdim = i * pdim;
		ixp1   = ( i - 1 ) * p1;
		
		memcpy( A22 + ixp1       , A + ixpdim           , size_psub );
		memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, size_mpsub );
	}
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves A11(2x2), A12(2x(p-2)), and A22((p-2)x(p-2))
// Like A11=A[e, e], A12=A[e, -e], and A22=A[-e, -e] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	int i, j, ixp, ij, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;

	A11[ 0 ] = A[ sub0 * pdim + sub0 ];
	A11[ 1 ] = A[ sub0 * pdim + sub1 ];
	A11[ 2 ] = A11[ 1 ];                   // for symmetric matrices
	A11[ 3 ] = A[ sub1 * pdim + sub1 ];
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp = i * pdim;
		
		A12[ i + i     ] = A[ ixp + sub0 ];
		A12[ i + i + 1 ] = A[ ixp + sub1 ];
	
		for( j = 0; j < sub0; j++ )
			A22[ j * p2 + i ] = A[ ixp + j ];

		for( j = sub0 + 1; j < sub1; j++ )
		{
			ij = ixp + j;
			A22[ ( j - 1 ) * p2 + i ] = A[ ij ];
			A22[ i * p2 + j - 1     ] = A[ ij ];
		}
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[ ( j - 2 ) * p2 + i ] = A[ ij ];
			A22[ i * p2 + j - 2     ] = A[ ij ];
		}
	}
 
	for( i = sub0 + 1; i < sub1; i++ )
	{
		ixp = i * pdim;
		
		A12[ i + i - 2 ] = A[ ixp + sub0 ];
		A12[ i + i - 1 ] = A[ ixp + sub1 ];
	
		for( j = sub0 + 1; j < sub1; j++ )
			A22[ ( j - 1 ) * p2 + i - 1 ] = A[ ixp + j ];
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[ ( j - 2 ) * p2 + i - 1 ] = A[ ij ];
			A22[ ( i - 1 ) * p2 + j - 2 ] = A[ ij ];
		}
	}
	
	for( i = sub1 + 1; i < pdim; i++ )
	{
		ixp = i * pdim;
				
		A12[ i + i - 4 ] = A[ ixp + sub0 ];
		A12[ i + i - 3 ] = A[ ixp + sub1 ];
		
		for( j = sub1 + 1; j < pdim; j++ )
			A22[ ( j - 2 ) * p2 + i - 2 ] = A[ ixp + j ];
	}
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves A11_inv ( 2 x 2 ), A21 ( ( p - 2 ) x 2 ), and A22 ( ( p - 2 ) x ( p - 2 ) )
// Like A11_inv=inv ( A[ e, e ] ), A21 = A[ -e, e ], and A22 = A[ -e, -e ] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrices_inv( double A[], double A11_inv[], double A21[], double A22[], int *row, int *col, int *p )
{
	int i, ixp, ixp2, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;
	int sub0xp = sub0 * pdim, sub1xp = sub1 * pdim, sub0_plus = sub0 + 1, sub1_plus = sub1 + 1;
	
	double a11 = A[ sub0 * pdim + sub0 ];
	double a12 = A[ sub0 * pdim + sub1 ];
	double a22 = A[ sub1 * pdim + sub1 ];

	double det_A11 = a11 * a22 - a12 * a12;
	A11_inv[ 0 ]   = a22 / det_A11;
	A11_inv[ 1 ]   = - a12 / det_A11;
	A11_inv[ 2 ]   = A11_inv[ 1 ];
	A11_inv[ 3 ]   = a11 / det_A11;
	
	int size_sub0      = sizeof( double ) * sub0;
	int size_sub1_sub0 = sizeof( double ) * ( sub1 - sub0_plus );
	int size_pdim_sub0 = sizeof( double ) * ( pdim - sub1_plus );
	
	memcpy( A21           , A + sub0xp            , size_sub0 );		
	memcpy( A21 + sub0    , A + sub0xp + sub0_plus, size_sub1_sub0 );	
	memcpy( A21 + sub1 - 1, A + sub0xp + sub1_plus, size_pdim_sub0 );	

	memcpy( A21 + p2           , A + sub1xp            , size_sub0 );		
	memcpy( A21 + p2 + sub0    , A + sub1xp + sub0_plus, size_sub1_sub0 );	
	memcpy( A21 + p2 + sub1 - 1, A + sub1xp + sub1_plus, size_pdim_sub0 );	
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp  = i * pdim;
		ixp2 = i * p2;

		memcpy( A22 + ixp2           , A + ixp            , size_sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, size_sub1_sub0 );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, size_pdim_sub0 );	
	}
 
	for( i = sub0_plus; i < sub1; i++ )
	{
		ixp  = i * pdim;
		ixp2 = ( i - 1 ) * p2;

		memcpy( A22 + ixp2           , A + ixp            , size_sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, size_sub1_sub0 );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, size_pdim_sub0 );	
	}
	
	for( i = sub1_plus; i < pdim; i++ )
	{
		ixp  = i * pdim;
		ixp2 = ( i - 2 ) * p2;
				
		memcpy( A22 + ixp2           , A + ixp            , size_sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, size_sub1_sub0 );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, size_pdim_sub0 );		
	}
}
      
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// inverse function for symmetric positive-definite matrices (p x p)
// WARNING: Matrix you pass is overwritten with the result
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void inverse( double A[], double A_inv[], int *p )
{
//	int info, dim = *p;
//	char uplo = 'U';

	// creating an identity matrix
//	#pragma omp parallel for
//	for( int i = 0; i < dim; i++ )
//		for( int j = 0; j < dim; j++ )
//			A_inv[ j * dim + i ] = ( i == j );
//	
//	// LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
//	F77_NAME(dposv)( &uplo, &dim, &dim, A, &dim, A_inv, &dim, &info FCONE );
    
    Eigen::MatrixXd mat(1,1);
    mat.resize(*p,*p);
    int k = 0;
  for(int i = 0; i < *p; i++)
  {
    for(int j = 0; j < *p; j++)
    {
      mat(i,j) = A[k];
      k++;
    }
  }
//  printf("original matrix: \n");
//  std::cout << mat << std::endl;
//  Eigen::MatrixXd g = mat;
  Eigen::MatrixXd mat_inv = mat.inverse();
  k = 0;
  for(int i = 0; i < *p; i++)
    {
      for(int j = 0; j < *p; j++)
        {
          A_inv[k] = mat_inv(i,j);
          k++;
        }
    }
//  printf("inverted matrix: \n");
//  std::cout << mat_inv << std::endl;
//    
//  printf("their product is: \n");
//  std::cout << g*mat_inv << std::endl;
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// inverse function for symmetric (2 x 2)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void inverse_2x2( double B[], double B_inv[] )
{
	double detB = B[ 0 ] * B[ 3 ] - B[ 1 ] * B[ 1 ];
	B_inv[ 0 ]  = B[ 3 ] / detB;
	B_inv[ 1 ]  = - B[ 1 ] / detB;
	B_inv[ 2 ]  = B_inv[ 1 ];
	B_inv[ 3 ]  = B[ 0 ] / detB;
} 
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Cholesky decomposition of symmetric positive-definite matrix
// A = U' %*% U
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void cholesky( double A[], double U[], int *p )
{
//	char uplo = 'U';
	int dim = *p; //info, 
	
	memcpy( U, A, sizeof( double ) * dim * dim );
    
    Eigen::MatrixXd mat(1,1);
    mat.resize(*p,*p);
    int k = 0;
    for(int i = 0; i < *p; i++)
    {
        for(int j = 0; j < *p; j++)
        {
          mat(i,j) = A[k];
          k++;
        }
    }
    //  printf("original matrix: \n");
    //  std::cout << mat << std::endl;
    //  Eigen::MatrixXd g = mat;
    Eigen::MatrixXd L(mat.llt().matrixL());
    Eigen::MatrixXd UU =L.transpose();
    k = 0;
    for(int i = 0; i < *p; i++)
    {
      for(int j = 0; j < *p; j++)
        {
          U[k] = UU(i,j);
          k++;
        }
    }
    
	
//	F77_NAME(dpotrf)( &uplo, &dim, U, &dim, &info FCONE );	
	
//	#pragma omp parallel for
//	for( int i = 0; i < dim; i++ )
//		for( int j = 0; j < i; j++ )
//			U[ j * dim + i ] = 0.0;
    
//      printf("original matrix is: \n");
//      std::cout << mat << std::endl;
//    
//    printf("product matrix is: \n");
//      std::cout << UU.transpose()*UU << std::endl;
//    
//    printf("upper triangular matrix is: \n");
//      std::cout << UU << std::endl;
    
}
  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//  Determinant of a symmetric possitive-definite matrix ( A )
//       > > > > > > > > >  WARNING: Matrix you pass is overwritten < < < < < < < < < 
//  For any symmetric PD Matrix A, we have: |A| = |T| ^ 2, where T is cholesky decomposition of A. 
//  Thus, |T| = \prod_{i = 1}^p T_{ii}.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_determinant( double A[], double *det_A, int *p )
{
//	char uplo = 'U';
	int dim = *p; //info, , dim1 = dim + 1;

//	F77_NAME(dpotrf)( &uplo, &dim, A, &dim, &info FCONE );	
    
    Eigen::MatrixXd mat(1,1);
    mat.resize(*p,*p);
    int k = 0;
    for(int i = 0; i < *p; i++)
    {
        for(int j = 0; j < *p; j++)
        {
          mat(i,j) = A[k];
          k++;
        }
    }
    //  printf("original matrix: \n");
    //  std::cout << mat << std::endl;
    //  Eigen::MatrixXd g = mat;
    Eigen::MatrixXd L(mat.llt().matrixL());
    k = 0;
	
	double result = 0.0;
	for( int i = 0; i < dim; i++ ) result += log(L(i,i));

	*det_A = 2.0 * result;
}
        
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To select an edge for BDMCMC algorithm  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void select_edge( double rates[], int *index_selected_edge, double *sum_rates, int *qp )
{
	int qp_star = *qp;

	// rates = sum_sort_rates
	vector<double>cumulative_rates( qp_star, 0.0 );
	cumulative_rates[ 0 ] = rates[ 0 ];
	for( int i = 1; i < qp_star; i++ )
		cumulative_rates[ i ] = cumulative_rates[ i - 1 ] + rates[ i ];
	
	*sum_rates = cumulative_rates[ qp_star - 1 ];
	
	// GetRNGstate();
	// This is essentially choosing a random value in the vector, where values
	// at a higher rate are selected with a higher probability.  (super cool)
	double random_value = *sum_rates * unif_rand(); // Rf_runif( 0.0, *sum_rates );
	// PutRNGstate();

	//int counter = 0;
	//while( random_value > cumulative_rates[ counter ] )	++counter;
	//*index_selected_edge = counter;
	 
	// To start, find the subscript of the middle position.
	int lower_bound = 0;
	int upper_bound = qp_star - 1;
	int position    = upper_bound / 2;  // ( lower_bound + upper_bound ) / 2;

	// This is likely the fastest way to search this vector, since it does not require us to 
	// check every element (most of the time).
	while( upper_bound - lower_bound > 1 )
	{
		 //if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
		( cumulative_rates[ position ] > random_value ) ? upper_bound = position : lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	*index_selected_edge = ( cumulative_rates[ position ] < random_value ) ? ++position : position;
} 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// select_edge for bd_for_ts
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void select_edge_ts( long double rates[], int *index_selected_edge, long double *sum_rates, int *qp )
{
	int qp_star = *qp;

	// rates = sum_sort_rates
	vector<long double>cumulative_rates( qp_star, 0.0 );
	cumulative_rates[ 0 ] = rates[ 0 ];
	for( int i = 1; i < qp_star; i++ )
		cumulative_rates[ i ] = cumulative_rates[ i - 1 ] + rates[ i ];
	
	*sum_rates = cumulative_rates[ qp_star - 1 ];
	
	// GetRNGstate();
	long double random_value = *sum_rates * unif_rand();  // Rf_runif( 0.0, *sum_rates );
	// PutRNGstate();

	//int counter = 0;
	//while( random_value > cumulative_rates[ counter ] )	++counter;
	//*index_selected_edge = counter;
	 
	// To start, find the subscript of the middle position.
	int lower_bound = 0;
	int upper_bound = qp_star - 1;
	int position    = upper_bound / 2;  // ( lower_bound + upper_bound ) / 2;

	while( upper_bound - lower_bound > 1 )
	{
		 //if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
		( cumulative_rates[ position ] > random_value ) ? upper_bound = position : lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	*index_selected_edge = ( cumulative_rates[ position ] < random_value ) ? ++position : position;
} 
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To simultaneously select multiple edges for BDMCMC algorithm  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void select_multi_edges( double rates[], int index_selected_edges[], int *size_index, double *sum_rates, int *multi_update, int *qp )
{
	int i, qp_star = *qp, qp_star_1 = qp_star - 1;

	// rates = sum_sort_rates
	vector<double>cumulative_rates( qp_star, 0.0 );
	cumulative_rates[ 0 ] = rates[ 0 ];
	for ( int i = 1; i < qp_star; i++ )
		cumulative_rates[ i ] = cumulative_rates[ i - 1 ] + rates[ i ];
	
	double max_bound = cumulative_rates[ qp_star_1 ];
	
	// - - - - - - for first edge - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
	// To start, find the subscript of the middle position.
	int lower_bound = 0;
	int upper_bound = qp_star_1;
	int position    = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

	//GetRNGstate();
	double random_value = max_bound * unif_rand();
	//PutRNGstate();

	while( upper_bound - lower_bound > 1 )
	{
		//if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
		( cumulative_rates[ position ] > random_value ) ? upper_bound = position : lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	if ( cumulative_rates[position] < random_value ) ++position;
	index_selected_edges[0] = position;
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

	int counter = 1, same;
	//GetRNGstate();
	for( int it = 0; it < 200 * *multi_update; it++ )
	{
		if( counter == *multi_update ) break;
		
		random_value = max_bound * unif_rand();
	
		// To start, find the subscript of the middle position.
		lower_bound = 0;
		upper_bound = qp_star_1;
		position    = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

		while( upper_bound - lower_bound > 1 )
		{
			// if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
			( cumulative_rates[ position ] > random_value ) ? upper_bound = position : lower_bound = position;     
			
			position = ( lower_bound + upper_bound ) / 2;
		}
		
		if( cumulative_rates[position] < random_value ) ++position;
		
		same = 0;
		for( i = 0; i < counter; i++ )
			if( index_selected_edges[ i ] == position )
				++same;

		if( same == 0 ) index_selected_edges[ counter++ ] = position;
	}
	//PutRNGstate();

	*size_index = counter;
	*sum_rates  = max_bound;
} 
         

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// NEW for Lang codes for Hermitian matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void Hsub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int i, l = 0, subj = *sub, pdim = *p, subxp = subj * pdim;

	for( i = 0; i < subj; i++ )
		sub_A[ l++ ] = -A[ subxp + i ];
	
	for( i = subj + 1; i < pdim; i++ )
		sub_A[ l++ ] = -A[ subxp + i ];
}
      
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// For Hermitian matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void Hsub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col, sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		sub_A[ l++ ] = -A[ sub0p + i ]; 
		sub_A[ l++ ] = -A[ sub1p + i ]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		sub_A[ l++ ] = -A[ sub0p + i ]; 
		sub_A[ l++ ] = -A[ sub1p + i ]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		sub_A[ l++ ] = -A[ sub0p + i ]; 
		sub_A[ l++ ] = -A[ sub1p + i ]; 
	}
}
       
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// sub_matrices1 for Hermitian matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void Hsub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, ixpdim, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim, mpsub = pdim - psub - 1;

	for( i = 0; i < psub; i++ )
		A12[ i ] = -A[ subxp + i ];
	
	for( i = psub; i < pdim - 1; i++ )
		A12[ i ] = -A[ subxp + i + 1 ];

	for( i = 0; i < psub; i++ )
	{	
		ixpdim = i * pdim;
		memcpy( A22 + i * p1       , A + ixpdim           , sizeof( double ) * psub );
		memcpy( A22 + i * p1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		ixpdim = i * pdim;
		memcpy( A22 + ( i - 1 ) * p1       , A + ixpdim           , sizeof( double ) * psub);
		memcpy( A22 + ( i - 1 ) * p1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}
}
        
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// sub_matrices for Hermitian matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void Hsub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	int i, i1, i2, ixp, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;

	A11[ 0 ] = A[ sub0 * pdim + sub0 ];
	A11[ 1 ] = A[ sub0 * pdim + sub1 ];
	A11[ 2 ] = -A11[ 1 ];                   // for symmetric matrices
	A11[ 3 ] = A[ sub1 * pdim + sub1 ];
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp = i * pdim;
		
		A12[ i + i ]     = A[ ixp + sub0 ];
		A12[ i + i + 1 ] = A[ ixp + sub1 ];

		memcpy( A22 + i * p2,            A + ixp,            sizeof( double ) * sub0 );
		memcpy( A22 + i * p2 + sub0,     A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );	
	}
 
	for( i = sub0 + 1; i < sub1; i++ )
	{
		ixp = i * pdim;
		i1 = i - 1;

		A12[ i + i - 2 ] = A[ ixp + sub0 ];
		A12[ i + i - 1 ] = A[ ixp + sub1 ];

		memcpy( A22 + i1 * p2,            A + ixp,            sizeof( double ) * sub0 );
		memcpy( A22 + i1 * p2 + sub0,     A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i1 * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );	
	}
	
	for( i = sub1 + 1; i < pdim; i++ )
	{
		ixp = i * pdim;
		i2  = i - 2;
				
		A12[ i + i - 4 ] = A[ ixp + sub0 ];
		A12[ i + i - 3 ] = A[ ixp + sub1 ];

		memcpy( A22 + i2 * p2,            A + ixp,            sizeof( double ) * sub0 );
		memcpy( A22 + i2 * p2 + sub0,     A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i2 * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );		
	}
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// inverse function for Hermitian (2 x 2)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void cinverse_2x2( double r_B[], double i_B[], double r_B_inv[], double i_B_inv[] )
{
	double r_det = r_B[0] * r_B[3] - i_B[0] * i_B[3] - ( r_B[1] * r_B[1] + i_B[1] * i_B[1] );
	double i_det = r_B[0] * i_B[3] + i_B[0] * r_B[3];
	double mod   = r_det * r_det + i_det * i_det;
	
	r_B_inv[0] =  ( r_B[3] * r_det + i_B[3] * i_det ) / mod;
	i_B_inv[0] =  ( r_det * i_B[3] - r_B[3] * i_det ) / mod;
	r_B_inv[1] = -( r_B[1] * r_det + i_B[1] * i_det ) / mod;
	i_B_inv[1] = -( r_det * i_B[1] - r_B[1] * i_det ) / mod;
	r_B_inv[2] = -( r_B[1] * r_det - i_B[1] * i_det ) / mod;
	i_B_inv[2] =  ( r_det * i_B[1] + r_B[1] * i_det ) / mod;
	r_B_inv[3] =  ( r_B[0] * r_det + i_B[0] * i_det ) / mod;
	i_B_inv[3] =  ( r_det * i_B[0] - r_B[0] * i_det ) / mod;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// For generating scale-free graphs: matrix G (p x p) is an adjacency matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void scale_free( int *G, int *p )
{
	int i, j, tmp, dim = *p, p0 = 2;
	double random_value;
	std::vector<int> size_a( dim ); 

	for( i = 0; i < p0 - 1; i++ )
	{
		G[         i * dim + i + 1 ] = 1;
		G[ ( i + 1 ) * dim + i     ] = 1;
	}
		
	for( i = 0 ; i < p0 ; i++ ) size_a[ i ] = 2;
	for( i = p0; i < dim; i++ ) size_a[ i ] = 0;
	
	int total = 2 * p0;
	
	GetRNGstate();
	for( i = p0; i < dim; i++ )
	{
		random_value = (double) total * unif_rand();
	   
		tmp = 0;
		j   = 0;
		
		while( tmp < random_value && j < i ) 
			tmp += size_a[ j++ ];
		
		j--;
		
		G[ i * dim + j ] = 1;
		G[ j * dim + i ] = 1;
		
		total += 2;
		size_a[ j ]++;
		size_a[ i ]++;
	}
	PutRNGstate();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To transfer the raw discreate data 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void transfer_data( int r_data[], int data[], int *n, int *p, int *size_unique_data )
{
	int i, j, l, counter;
	
// - - tranfer each row of raw data to string - - - - - - - - - - - - - - - - - - - - - - - - - - -|
	vector<char> char_row( *p );             
	vector<string>all_patterns( *n );
	string *unique_patterns = new string[ *n ];
	
	for( i = 0; i < *n; i++ )
	{
		for( j = 0; j < *p; j++ )
			char_row[ j ] = r_data[ j * *n + i ] + '0';
		
		all_patterns[ i ] = string( char_row.begin(), char_row.end() );
	}

// - - find the unique string-rows - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	unique_patterns[0] = all_patterns[0];
	int length_unique_patterns = 1;
	for( i = 1; i < *n; i++ )
	{
		counter = 0;
		//for( j = 0; j < length_unique_patterns; j++ )
			//( all_patterns[i] == unique_patterns[j] ) ? j = length_unique_patterns : ++counter;					
		while( ( counter < length_unique_patterns ) and ( all_patterns[ i ] != unique_patterns[ counter ] ) )
			++counter;
		
		if( counter == length_unique_patterns )
			unique_patterns[ length_unique_patterns++ ] = all_patterns[ i ];
	}

// - - tranfer the data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
	int which_one = 0;
	for( l = 0; l < length_unique_patterns; l++ )  
	{
		counter = 0;
		for( i = 0; i < *n; i++ )
			if( all_patterns[ i ] == unique_patterns[ l ] ) 
			{
				counter++;
				which_one = i;
			}
			
		data[ *p * *n + l ] = counter;
		
		for( j = 0; j < *p; j++ )
			data[ j * *n + l ] = r_data[ j * *n + which_one ]; 
	}
	
	*size_unique_data = length_unique_patterns;
	
	delete[] unique_patterns;
}

