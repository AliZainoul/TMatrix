// classe vecteur
#include <iostream>
#include "TVector.hpp"
#include "full_mat_c.hpp"
#include "sparse_mat_c.hpp"

using namespace std;


int main(){

	cout.setf(ios::scientific, ios::floatfield);
	//cout.setf(ios::fixed, ios::floatfield);
	cout.precision(6);

	// first test for real number arithmetic

	cout << " ******** Testing for FullMtx ******** \n";

  int n = 4, m = 4;

  FullMtx<long double> mat1(n,m);
  for (int i = 0; i < n; i++)
  for (int j = 0; j < m; j++)
	mat1[i][j] =  n/(abs(i - j) + 1.0);
	cout << "mat1";
  cout << mat1;

	FullMtx<long double>matx(n,m,mat1.getmatrix());
	cout << "matx"<< endl;
	cout << matx<< endl;

  Vector<long double> vec1(n) ;
  for (int i = 0; i < n; i++) vec1[i] = i;
	cout << "vec1"<< endl;
	cout << vec1<< endl;

  Vector<long double> vecx = vec1;
  long double doot = dot<long double>(vec1, vecx);
	cout << "doot"<< endl;
	cout << doot<< endl;

	Vector<long double> vectrue(2);
	vectrue[0] = 9.09e-02;
	vectrue[1] = 6.36336e-01;
	cout << "vectrue" << endl;
	cout << vectrue<< endl;


	FullMtx<long double> matxx({{4,1},{1,3}});
	cout << "matxx"<< endl;
	cout << matxx<< endl;
  int prec = 1;
  int iter = 2;
  long double eps = 1.0e-14;
  Vector<long double> vec2(2,0.0);
	vec2[0] = 2;
	vec2[1] = 1;
	cout << "vec2"<< endl;
	cout << vec2 << endl;
	Vector<long double> vec2x(2,0.0);
	vec2x[0] = 1;
	vec2x[1] = 2;
	cout << "vec2x"<< endl;
	cout << vec2x << endl;

  int ret =  matxx.CG(vec2, vec2x, eps, iter, prec);
  if (ret == 0) cout << "CG returned successfully" << endl;
	cout << "The true solution to Ax = b is: " << endl << vectrue << endl;
	cout << "CGsolution =  " << endl << vec2 << endl;

  cout << iter << " iterations used. " ;
  cout << "Residual in CG= " << eps << " \n " ;
  cout << "Ture error in CG= " << (vectrue-vec2).maxnorm() << "\n" ;



	/*
  cout << "******** Sparse Matrix  ******** \n";
  SparseMtx<long double> sm1(n, n*n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      sm1[i*n + j] = mat1[i][j];
      sm1.getclm(i*n + j) = j;
    }
    sm1.getfnz(i) = i*n;
  }
  sm1.getfnz(n) = n*n;
	cout << "********************************" << endl;
  //cout << sm1;
	cout << "********************************" << endl;

  iter = n;
  eps = 1.0e-14;
  vec2.reset();

	cout << "********************************" << endl;
	cout << sm1*vec1;
	cout << "********************************" << endl;
	FullMtx<long double>	identity(n,n);
	for (unsigned int i=0; i<n ; i++)
	for (unsigned int j=0; j<n ; j++)
	if (i==j)	identity[i][j]= 1.0; else	identity[i][j] = 0.0;
	cout << "********************************" << endl;
	cout <<	identity;
	cout << "********************************" << endl;
	//for (unsigned int k=0;k<n-1;k++) cout << sm1[k] << endl;
	cout << "********************************" << endl;

  ret =  sm1.CG(vec2, sm1*vec1,eps,iter,prec);
  if (ret == 0) cout << "CG returned successfully\n";
	cout << "True solution is: " << vec1 << " ";
	cout << "Conjuguate Gradient solution =  " << vec2 << "\n";
	if (iter == 1)
	cout << iter << " iteration used. " ;
	else
	cout << iter << " iterations used. " ;
  cout << "Residual in CG = " << eps << "  " ;
  cout << "Ture error in CG= " << (vec2-vec1).maxnorm() << "\n\n" ;
	*/
}
