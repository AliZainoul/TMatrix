// classe vecteur
#include <iostream>
#include "TVector.hpp"
#include "full_mat_c.hpp"
#include "sparse_mat_c.hpp"

using namespace std;


int main(){

	cout.setf(ios::scientific, ios::floatfield);
	//cout.setf(ios::fixed, ios::floatfield);
	cout.precision(3);

	// first test for real number arithmetic

	cout << " ******** Testing for FullMtx ******** \n";

  int n = 10, m = 10;

  FullMtx<double> mat1(n,m);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++) {
      mat1[i][j] =  n/(abs(i - j) + 1.0);
    }
  }
  cout << mat1;

  Vector<double> vec1(n) ;
  for (int i = 0; i < n; i++) vec1[i] = i;
  cout << vec1;
  Vector<double> vecx = vec1;
  double alpha = dot<double>(vec1, vecx);
  cout << alpha;

    int prec = 1;
    int iter = n;
    double eps = 1.0e-14;
    Vector<double> vec2(n,0.0);
    int ret =  mat1.CG(vec2, mat1*vec1, eps, iter, prec);
    if (ret == 0) cout << "CG returned successfully\n";
	  cout << "true solutionn is: " << vec1 << "   ";
	  cout << "CGsolution =  " << vec2 << "\n";
  cout << iter << " iterations used. " ;
  cout << "Residual in CG= " << eps << " \n " ;
  cout << "Ture error in CG= " << (vec2-vec1).maxnorm() << "\n" ;




  cout << "******** Sparse Matrix  ******** \n";
  SparseMtx<double> sm1(n, n*n);
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
	FullMtx<double>	identity(n,n);
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

}
