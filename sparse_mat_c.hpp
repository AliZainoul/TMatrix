#ifndef SPARSE_MAT_C_H
#define SPARSE_MAT_C_H
#include "abstract_mat_c.hpp"
#include "TVector.hpp"
using namespace std;

// Compressed Sparse Row Format
template<class T>
class SparseMtx: public AbsMtx<T>
{
private:
  int lenth;    // Number of nonzero entries of the original matrix
  T* sra;       // Array for storing the non-zero entries
  int* clm;     // Column indexes in the original matrix of the entries in sra
  int* fnz;     // Position in sra of first non-zero entires of each row

public:
  // Constructor 0 : Constructor
  SparseMtx(int n, int m, T* t, int* c, int* f);
  // n: Number of rows (and columns) of the original matrix;
  // m: length of array sra for nonzero entries;
  // t: nonzero entries of the original matrix;
  // c: colunm indexes (in the original matrix) of entries in sra;
  // f: index in sra of first nonzero entry in each row;

  // Constructor 1 :  Copy Constructor
  SparseMtx(int n, int m);
  // Initialize all entries to zero
  // n: number of rows (and columns);
  // m: number of nonzero entries;

  // Constructor 2 :  Copy Constructor
  SparseMtx(const SparseMtx&);

  // Destructor
  ~SparseMtx(){ delete[] sra; delete[] fnz; delete[] clm; }

  T& operator[](int i) const { return sra[i]; }  // Subscripting
  int& getfnz(int i) const { return fnz[i]; }    // First !0 entry of each row
  int& getclm(int i) const { return clm[i]; }    // column index

  SparseMtx& operator=(const SparseMtx&);        // Overload of Operator '='
  Vector<T> operator*(const Vector<T>&) const;   // Matrix-Vector Product

  //template <class S>
  // friend ostream& operator<<(ostream&, const SparseMtx<S>&);
  // overload <<

};

// **************************
// IMPLEMENTATION
// **************************
template<class T>
SparseMtx<T>::SparseMtx(int n, int m, T* et, int* cn, int* da) {
  this->nrows = n;
  this->lenth = m;
  this->sra = new T [this->lenth];
  this->clm = new int [this->lenth];
  this->fnz = new int [this->nrows +1];

  for (int i =0; i< this->lenth; i++)
  {
    sra[i] = et[i];
    clm[i] = cn[i];
  }
  for (int i = 0; i <= this->nrows; i++) fnz[i] = da[i];
}

template<class T>
SparseMtx<T>::SparseMtx(int n, int m) {
  this->nrows = n;
  this->lenth = m;
  this->sra = new T [lenth];
  this->clm = new int [lenth];
  this->fnz = new int [this->nrows +1];

  for (int i =0; i< this->lenth; i++)
  {
    sra[i] = 0;
    clm[i] = 0;
  }
  for (int i =0; i <= this->nrows; i++)  fnz[i] = 0;
}

template<class T>
SparseMtx<T>::SparseMtx(const SparseMtx & mat) {           // copy constructor
  this->nrows = mat.nrows;
  this->lenth = mat.lenth;
  this->sra = new T [lenth];
  this->clm = new int [lenth];
  this->fnz = new int [this->nrows +1];

  for (int i = 0; i < this->lenth; i++) {
    sra[i] = mat[i];
    clm[i] = mat.clm[i];
  }
  for (int i =  0; i <= this->nrows; i++) fnz[i] = mat.fnz[i];
}

template<class T>
SparseMtx<T>& SparseMtx<T>::operator=(const SparseMtx & ssm) {
  if(this->nrows != ssm.nrows || this->lenth != ssm.lenth)
    error("bad matrix sizes in SparseMtx::operator=()");
  for (int i = 0; i < this->lenth; i++) {
    sra[i]  = ssm[i];
    clm[i]  = ssm.clm[i];
  }
  for (int i = 0; i <= this->nrows; i++) fnz[i] = ssm.fnz[i];
  return *this;
}

template<class T>
Vector<T> SparseMtx<T>::operator*(const Vector<T>& vec) const {
  if (this->nrows != vec.size())
    error("matrix and vector sizes do not match in SparseMtx::operator*()");
  Vector<T> tm(this->nrows);
  for (int i = 0; i < this->nrows; i++)
    for (int j = fnz[i]; j < fnz[i +1]; j++) tm[i] += sra[j]*vec[clm[j]];
  return tm;
}
/*
template<class T>
ostream& operator<<(ostream& s, const SparseMtx<T>& mat) {
  double alpha = 0.0;
  for (int i = 0; i < mat.nrows; i++)
  for (int j = 0; j < mat.nrows; j++)
  if
  else alpha = 0.0;
  s <<  alpha << " " << "\n";
  return s;

  int lenth;             // # of nonzero entries of the original matrix
  T* sra;                // array for storing the non-zero entries
  int* clm;              // column indexes in matrix of the entries in sra
  int* fnz;              // position in sra of first nonzero entires of each row
}*/
// **************************

#endif
