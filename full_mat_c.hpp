#ifndef FULL_MAT_C_H
#define FULL_MAT_C_H

#include "abstract_mat_c.hpp"
#include "TVector.hpp"
#include "error.hpp"

using namespace std;

template<class T>
class FullMtx: public AbsMtx<T>{

private:
  int ncols;   // Number of columns in the matrix
  T** mx;      // Entries of the matrix

public:
  // Constructor 0: n: Number of rows, m: Number of columns, T: Entries
  FullMtx(int n, int m, T**);
  // Constructor 1: All entries are set to 't'
  FullMtx(int n, int m, T t = 0);
  // constructeur 2: List Constructor - Matrix = {{1,2},{3,4}} = [1 2 ; 3 4]
  FullMtx(initializer_list<initializer_list<T>> lst);
  // Constructor 3: Copy Constructor
  FullMtx(const FullMtx&);

  // Destructor
  ~FullMtx(){
    for (int i = 0; i< this->nrows; i++) delete[]  mx[i];
    delete[] mx;
  }

  // Implement + as a member fcn and  implement - as a friend
  FullMtx& operator=(const FullMtx&);           // Overload of Operator '='
  Vector<T> operator*(const Vector<T>&) const;  // Matrix-Vector Product
  T* operator[](int i) const { return mx[i]; }  // Method that Returns i-th row
  T** getmatrix() const {return mx;}

  // Overload of Operator '<<'
  template <class S>
  friend ostream& operator<<(ostream&, const FullMtx<S>&);
};


// Implementation
// ********************
template<class T>
FullMtx<T>::FullMtx(int n, int m, T** dbp) {
  this->nrows = n;
  this->ncols = m;
  mx = new T* [this->nrows];
  for (int i =  0; i< this->nrows; i++) {
    mx[i] = new T [this->ncols];
    for (int j = 0; j < this->ncols; j++) mx[i][j] = dbp[i][j];
  }
}

template<class T>
FullMtx<T>::FullMtx(int n, int m, T a) {
  this->nrows = n;
  this->ncols = m;
  mx = new T* [this->nrows];
  for (int i =  0; i< this->nrows; i++) {
    mx[i] = new T [this->ncols];
    for (int j = 0; j < this->ncols; j++) mx[i][j] = a;
  }
}

template<class T>
FullMtx<T>::FullMtx(const FullMtx& mat) {
  this->nrows = mat.nrows;
  this->ncols = mat.ncols;
  mx = new T* [this->nrows];
  for (int i =  0; i< this->nrows; i++) {
    mx[i] = new T [this->ncols];
    for (int j = 0; j < this->ncols; j++) mx[i][j] = mat[i][j];
  }
}

template<class T>
FullMtx<T>::FullMtx(initializer_list<initializer_list<T>> lst){
  this->nrows = lst.size();
  this->ncols = (*lst.begin()).size();
  mx = new T* [this->nrows];
  int i=0;
  for (const auto& l : lst) {
      mx[i] = new T [this->ncols];
	  copy(l.begin(), l.end(),mx[i] );
	  i++;
  }
}

template<class T> // Usefull for testing small matrices
ostream& operator<<(ostream& s, const FullMtx<T>& mat) {
  for (int i =0; i< mat.nrows; i++) {
    s << "\n" << i << "-th row:    \n";
    for (int j =0; j< mat.ncols; j++)
    {
      s << mat.mx[i][j] << "  ";
      if (j%10 == 9) s << "\n";
    }
    s << "\n";
  }
  return s;
}

template<class T>
FullMtx<T> & FullMtx<T>::operator=(const FullMtx & mat) {
  if (this != &mat) {
    if (this->nrows !=mat.nrows || this->ncols !=mat.ncols)
    error("Bad matrix sizes in FullMtx::operator=()");
    for (int i = 0; i < this->nrows; i++)
      for (int j = 0; j < this->ncols; j++) mx[i][j]  = mat.mx[i][j];
  }
  return *this;
}

template<class T>
Vector<T> FullMtx<T>::operator*(const Vector<T> & vec) const {
  if (this->ncols != vec.size())
  error("Matrix and Vector sizes does not match in FullMtx::operator*()");
  Vector<T> tm(this->nrows);
  for (int i = 0; i < this->nrows; i++)
  for (int j = 0; j < this->ncols; j++) tm[i] += mx[i][j]*vec[j];
  return tm;
}
// ********************
#endif
