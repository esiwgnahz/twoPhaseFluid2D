#include "pointers.h"
#include <stdio.h>
//#include "string.h"
#include <math.h>
#include <algorithm>
//#include "point3d.h"
#ifndef M_PI
  /* Constants rounded for 21 decimals. */
  #define M_E         2.71828182845904523536
  #define M_LOG2E     1.44269504088896340736
  #define M_LOG10E    0.434294481903251827651
  #define M_LN2       0.693147180559945309417
  #define M_LN10      2.30258509299404568402
  #define M_PI        3.14159265358979323846
  #define M_PI_2      1.57079632679489661923
  #define M_PI_4      0.785398163397448309616
  #define M_1_PI      0.318309886183790671538
  #define M_2_PI      0.636619772367581343076
  #define M_1_SQRTPI  0.564189583547756286948
  #define M_2_SQRTPI  1.12837916709551257390
  #define M_SQRT2     1.41421356237309504880
  #define M_SQRT_2    0.707106781186547524401
#endif
  #define G_GRAVITY           9.80665 // м/с^2
  #define GAS_CONST        8.3144621 //[Дж/(моль*K)]

#if !defined M_SQRT_2 && defined M_SQRT1_2 // VC++2005
  #define M_SQRT_2 M_SQRT1_2
#endif

#if !defined M_1_SQRTPI // VC++2005
  #define M_1_SQRTPI  0.564189583547756286948
#endif

#ifndef RBF_FUNCTIONS_FLAG
enum RBF_FUNCTIONS
{
  MULTIQUADRIC_RBF         = 0,
  INVERSE_MULTIQUADRIC_RBF     = 1,
  INVERT_QUADRIC_RBF       = 2,
  GAUSS_RBF             = 3,
  THIN_PLATE_SPLINE_RBF         = 4,
  DISTANCE_RBF           = 5

  //MATERN_LINER_RBF       = 6,
  //MATERN_SQUARE_RBF       = 7,
  //BIGARMONIC_LAPLACE_RBF     = 8,
  //MACDONALD_RBF           = 9,
  //BESSEL_RBF               = 10
};
#endif

template <class T> inline T radToGrad(T angle){
  return angle/(M_PI/180);
}

template <class T> inline T gradToRad(T angle){
  return angle*(M_PI/180);
}

template <class T> inline int round(T value){
  return (value > 0) ? (int) (value + 0.5) : (int) (value - 0.5);
}

/* not used
inline double round(double value, int digits)
{
  double d=pow(10.0,digits);
  return round(value*d)/d;
}
*/

///////////////////////////
// Why following functions and classes are placed here?

struct APPR_STRUCT{
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;
};

//int porabol_approximation (double *a, double *b, double *c, APPR_STRUCT * ptStruct);

  //  classes declaration
template <class T> class TBasicArray;
template <class T> class TBasicArray2D;
template <class T> class TMatrix;
template <class T> class TVector;
class MCRadialBasis;

/****************************************************    TYPEDEF      *****************************************************/
typedef TBasicArray<int> MCBasicArray_int;
typedef TBasicArray<double> MCBasicArray;
typedef TBasicArray<char> MCBasicArray_char;
typedef TBasicArray<float> MCBasicArray_float;

typedef TBasicArray2D<int> MCBasicArray2D_int;
typedef TBasicArray2D<double> MCBasicArray2D;
typedef TBasicArray2D<char> MCBasicArray2D_char;
typedef TBasicArray2D<float> MCBasicArray2D_float;

typedef TVector<int> MCVector_int;
typedef TVector<double> MCVector;
typedef TVector<float> MCVector_float;

typedef TMatrix<int> MCMatrix_int;
typedef TMatrix<double> MCMatrix;
typedef TMatrix<float> MCMatrix_float;

  // Function to Inverse Matrix by it's LU Decomposition
template <class T>
TMatrix<T> MatrixInverse(const TMatrix<T>&);

//***************************************************************************************
template <class T> class TBasicArray
{
protected:
  int m_Err;
    MCAPtr<T> VectorItems;
    int nn;
public:
  int size() const
    { return nn; }
    T* pointer() const
    { Assert(!m_Err); return &(*VectorItems); }
    int  IsErr () const    // 0 - no Err
    { return m_Err; }

  TBasicArray();
  TBasicArray(const TBasicArray<T> &);
    explicit  TBasicArray(const int);
    template <class Y> TBasicArray(const int _nn, const Y& a)
    {
      Initialize(_nn);
      for (int i=0; i<nn; i++) {
          VectorItems[i] = (T) a;
      }
    }
    template <class Y> TBasicArray(const int _nn, const Y* a)
    {
      Initialize(_nn);
      for (int i=0; i<nn; i++) {
          VectorItems[i]= (T) a[i];
      }
    }
  virtual ~TBasicArray(){ };

  T& operator [] (const int i) const
    { return VectorItems[i]; }
   
  void Initialize (const int);
    void Initialize (const TBasicArray<T> &);
private:
  int operator == (const TBasicArray<T>&) const;
  int operator != (const TBasicArray<T>&) const;
  TBasicArray<T>& operator = (const TBasicArray<T>&);
  template <class Y> TBasicArray<T>& operator = (const Y& a){
    Assert(!m_Err);
    for (int i=0; i< nn; i++)
      VectorItems[i] = (T) a;
    return *this;
  }
};

template <class T> class TBasicArray2D
{
protected:
  int m_Err;
    int ni;  //row
    int nj;  //column
    MCAPtr< MCAPtr<T> > MatrixItems;
public:
    int row() const 
    { return ni; }
    int column() const
    { return nj; }
  MCAPtr<T>* pointer() const
    {/* Assert(!m_Err);*/ return &(*MatrixItems); }
    int IsErr() const
    { return m_Err; }    // 0 - no Error

  TBasicArray2D();
  TBasicArray2D(const TBasicArray2D<T> &);
    TBasicArray2D(const int,const int);
    template <class Y> TBasicArray2D(const int _ni,const int _nj, const Y& a)
    {
      Initialize(_ni,_nj);
      for (int i=0; i<ni; i++)
      for(int j=0; j<nj; j++)
      MatrixItems[i][j] = (T) a;
    }
    template <class Y> TBasicArray2D(const int _ni,const int _nj, const Y* a)
    {
      Initialize(_ni,_nj);
      for (int i=0; i<ni; i++)
      for(int j=0; j<nj; j++)
      MatrixItems[i][j] = (T) a[i][j];
    }
  virtual ~TBasicArray2D(){ };
  
    MCAPtr<T>& operator [] (const int i) const
    { return MatrixItems[i]; }

    void Initialize(const int,const int);
    void Initialize(const TBasicArray2D &);
private:
  int operator == (const TBasicArray2D<T>&) const;
  int operator != (const TBasicArray2D<T>&) const;
  TBasicArray2D<T>& operator = (const TBasicArray2D<T>&);
  template <class Y> TBasicArray2D<T>& operator = (const Y& a){
    Assert(!m_Err);
    for (int i=0; i< ni; i++)
        for(int j=0; j< nj; j++)
          MatrixItems[i][j] = (T) a;
    return *this;
  }
};

template <class T> class TVector : public TBasicArray<T>
{
public:
  TVector() : TBasicArray<T>() { }
  explicit  TVector(const int num) : TBasicArray<T>(num) { }
  TVector(const TVector &A) : TBasicArray<T>(A) { }
  template <class Y> TVector(const int _ni,const Y &a) : TBasicArray<T>(_ni, a) { }
  template <class Y> TVector(const int _ni,const Y* a) : TBasicArray<T>(_ni, a) { }
    ~TVector(){ };

  int operator == (const TVector<T>&) const;
  int operator != (const TVector<T>&) const;
  TVector<T>& operator = (const TVector<T>&);
    template <class Y> TVector<T>& operator = (const Y &a){
    Assert(!m_Err);
    for (int i=0; i< nn; i++)
      VectorItems[i] = (T) a;
    return *this;
    }
  TVector<T> operator + (const TVector<T>&) const;
  template <class Y> TVector<T> operator + (const Y &a) const {
    Assert(!m_Err);
        TVector<T> Temp(nn);
        for (int i=0; i< nn; i++)
      Temp[i] = (T) a + VectorItems[i];
        return Temp;
    }
    TVector<T> operator - (const TVector<T>&) const;
    template <class Y> TVector<T> operator - (const Y &a) const {
    Assert(!m_Err);
    TVector<T> Temp(nn);
    for (int i=0; i< nn; i++)
      Temp[i] = VectorItems[i] - (T) a;
    return Temp;
    }
    template <class Y> TVector<T> operator * (const Y &a) const {
    Assert(!m_Err);
      TVector<T> Temp(nn);
        for (int i=0; i< nn; i++)
            Temp[i] = (T) a * VectorItems[i];
        return Temp; 
    }
  TVector<T> operator * (const TMatrix<T>&) const ;
    template <class Y> TVector<T> operator / (const Y &a) const {
    Assert(!m_Err);
      TVector<T> Temp(nn);
        for (int i=0; i< nn; i++)
            Temp[i] = VectorItems[i]/((T) a);
        return Temp; 
    }
    template <class Y> TVector<T>& operator *= (const Y &a){
    Assert(!m_Err);
      for (int i=0; i<nn; i++)
          VectorItems[i] *= (T) a;
      return *this;
    }   
    template <class Y> TVector<T>& operator /= (const Y &a){
    Assert(!m_Err);
      for (int i=0; i<nn; i++)
          VectorItems[i] /= (T) a;
      return *this;
    }
    TVector<T>& operator += (const TVector<T>&);
    template <class Y> TVector<T>& operator += (const Y &a){
    Assert(!m_Err);
      for (int i=0; i<nn; i++)
          VectorItems[i] += (T) a;
      return *this;
    }
    TVector<T>& operator -= (const TVector<T>&);
    template <class Y> TVector<T>& operator -= (const Y &a){
    Assert(!m_Err);
      for (int i=0; i<nn; i++)
          VectorItems[i] -= (T) a;
      return *this;
    }
  
  T MultScalar (const TVector<T>&) const ;  // T = TVector * TVector
  TMatrix<T> MultMatrix(const TVector<T>&) const ;  // TMatrix = TVector * TVector
  template <class Y, class Z> friend TVector<Z> operator * (const Y &, const TVector<Z>& );
  template <class Y, class Z> friend TVector<Z> operator + (const Y &, const TVector<Z>& );
  template <class Y, class Z> friend TVector<Z> operator - (const Y &, const TVector<Z>& );

  T Norm() const ;   // Square root of the sum of squares
    void Normalize();
    T Max() const ;
    T Min() const ;
    void Print (const char * szFile) const ;
    void Probe () const ;    //create random T vector
};

template <class T> class TMatrix : public TBasicArray2D<T>
{
public:
    TMatrix() : TBasicArray2D<T>() { };
    TMatrix(const int _ni,const int _nj) : TBasicArray2D<T>(_ni,_nj) { }
    TMatrix(const TMatrix<T> &A) : TBasicArray2D<T>(A) { }
  template <class Y> TMatrix(const int _ni,const int _nj,const Y& a) : TBasicArray2D<T>(_ni,_nj,a) { }
    template <class Y> TMatrix(const int _ni,const int _nj,const Y* a) : TBasicArray2D<T>(_ni,_nj,a) { }
    ~TMatrix() { }

  int operator == (const TMatrix<T>&) const;
  int operator != (const TMatrix<T>&) const;
  TMatrix<T>& operator = (const TMatrix<T>&);
    template <class Y> TMatrix<T>& operator = (const Y &a){
    Assert(!m_Err);
      for (int i=0; i< ni; i++)
          for(int j=0; j< nj; j++)
              MatrixItems[i][j] = (T) a;
      return *this;
    }
    TMatrix<T> operator + (const TMatrix<T>&) const ;
    template <class Y> TMatrix<T> operator + (const Y &a) const {
    Assert(!m_Err);
      TMatrix<T> Temp(ni,nj);
      for (int i=0; i< ni; i++)
          for(int j=0; j< nj; j++)
            Temp[i][j] = (T) a + MatrixItems[i][j];
        return Temp;
    }
    TMatrix<T> operator - (const TMatrix<T>&) const;
  TMatrix<T> operator - () const;
    template <class Y> TMatrix<T> operator - (const Y &a) const {
    Assert(!m_Err);
        TMatrix<T> Temp(ni,nj);
      for (int i=0; i< ni; i++)
          for(int j=0; j< nj; j++)
            Temp[i][j] = MatrixItems[i][j] - (T) a;
      return Temp;
    }
    TMatrix<T> operator * (const TMatrix<T>&) const ;
  TVector<T> operator * (const TVector<T>&) const ;
  template <class Y> TMatrix<T> operator * (const Y &a) const {
    Assert(!m_Err);
        TMatrix<T> Temp(ni,nj);
      for (int i=0; i< ni; i++)
          for(int j=0; j< nj; j++)
            Temp[i][j] = ((T) a)*MatrixItems[i][j];
        return Temp; 
    }
    template <class Y> TMatrix<T> operator / (const Y &a) const {
    Assert(!m_Err);
        TMatrix<T> Temp(ni,nj);
      for (int i=0; i< ni; i++)
          for(int j=0; j< nj; j++)
            Temp[i][j] = MatrixItems[i][j]/((T) a);
        return Temp; 
    }
    template <class Y> TMatrix<T>& operator *= (const Y &a){
    Assert(!m_Err);
      for (int i=0; i< ni; i++)
          for(int j=0; j< nj; j++)
            MatrixItems[i][j] *= (T) a;
      return *this;
    }   
    template <class Y> TMatrix<T>& operator /= (const Y &a){
    Assert(!m_Err);
      for (int i=0; i< ni; i++)
          for(int j=0; j< nj; j++)
            MatrixItems[i][j] /= (T) a;
      return *this;
    }
    TMatrix<T>& operator += (const TMatrix<T>&);
    template <class Y> TMatrix<T>& operator += (const Y &a){
    Assert(!m_Err);
      for (int i=0; i< ni; i++)
          for(int j=0; j< nj; j++)
            MatrixItems[i][j] += (T) a;
      return *this;
    }
    TMatrix<T>& operator -= (const TMatrix<T>&);
    template <class Y> TMatrix<T>& operator -= (const Y &a){
    Assert(!m_Err);
      for (int i=0; i< ni; i++)
          for(int j=0; j< nj; j++)
            MatrixItems[i][j] -= (T) a;
      return *this;
    }
 
  template <class Y,class Z> friend TMatrix<Z> operator * (const Y &, const TMatrix<Z> &);
  template <class Y,class Z> friend TMatrix<Z> operator + (const Y &, const TMatrix<Z> &);
  template <class Y,class Z> friend TMatrix<Z> operator - (const Y &, const TMatrix<Z> &);
 
  int row() const 
    { return ni; }
  TVector<T> row(const int);  // return 'i' row
  int column() const
    { return nj; }
  TVector<T> column(const int);  // return 'j' column
    int IsSymmetric() const ;
    T Max() const ;
    T Min() const ;
    TMatrix Transpose() const ;
  T Trace() const ;
    void Print (const char * szFile) const ;
    void Probe () const ;    //create random T matrix

    // Ruzana
    void Inverse();  
    // C = A*B для выполнения меньшего количества операций  
    void Multiplication_and_Equating(const TMatrix<T>& A,const TMatrix<T>& B);
};

/****************************************************   TBasicArray    *****************************************************/
template <class T>
void TBasicArray<T> :: Initialize(const int _nn){
 // Assert(_nn!=0);
  m_Err=0;
  nn=_nn;

  VectorItems=new T[nn];
    memset (VectorItems, 0, sizeof(T)*nn);
}
template <class T>
void TBasicArray<T> :: Initialize(const TBasicArray<T> &A){
 // Assert(!A.IsErr());
  m_Err=0;
  nn = A.size();

  VectorItems = new T[nn];
  for (int i=0; i<nn; i++)
    VectorItems[i] = (T) A.VectorItems[i];
}
template <class T>
TBasicArray<T> :: TBasicArray()
  : m_Err(1),
  nn(0)
{ }
template <class T>
TBasicArray<T> :: TBasicArray(const int _nn)
{
  Initialize(_nn);
}
template <class T>
TBasicArray<T> :: TBasicArray(const TBasicArray<T> &A)
{
  Initialize(A);
}
template <class T>
int TBasicArray<T> :: operator == (const TBasicArray<T> &A) const
{
  if(nn!=A.size() || m_Err) return 0;
  for(int i=0; i < nn; i++)
    if(VectorItems[i]!=A[i])return 0;
  return 1;
}
template <class T>
int TBasicArray<T> :: operator != (const TBasicArray<T> &A) const
{
  if(nn!=A.size() || m_Err) return 1;
  for(int i=0; i < nn; i++)
    if(VectorItems[i]!=A[i])return 1;
  return 0;
}
template <class T>
TBasicArray<T>& TBasicArray<T> :: operator = (const TBasicArray<T> &A)
{
  Assert(!m_Err);
  Assert(nn==A.size());
  for (int i=0; i<nn; i++)
    VectorItems[i] = A.VectorItems[i];
  return *this;
}
/****************************************************   TBasicArray2D    *****************************************************/
template <class T>
void TBasicArray2D<T> :: Initialize(const int _ni, const int _nj){
  m_Err = 0;
  ni = _ni;
  nj = _nj;

  MatrixItems = new MCAPtr<T>[ni];
  for (int i=0; i<ni; i++) {
    MatrixItems[i] = new T[nj];
    memset (MatrixItems[i], 0, sizeof(T)*nj);
  }
}
template <class T>
void TBasicArray2D<T> :: Initialize(const TBasicArray2D<T> &A){
  //Assert(!A.IsErr());
  m_Err=0;
  ni = A.row();
  nj = A.column();

  MatrixItems = new MCAPtr<T>[ni];
  for (int i=0; i<ni; i++) {
    MatrixItems[i] = new T[nj];
  }
  for (int i=0; i< ni; i++)
    for (int j=0; j< nj; j++)
      MatrixItems[i][j] = A.MatrixItems[i][j];
}
template <class T>
TBasicArray2D<T> :: TBasicArray2D()
  : m_Err(1),
   ni(0),
   nj(0)
{ }
template <class T>
TBasicArray2D<T> :: TBasicArray2D(const int _ni, const int _nj)
{
  Initialize(_ni,_nj);
}
template <class T>
TBasicArray2D<T> :: TBasicArray2D(const TBasicArray2D<T> &A)
{
  Initialize(A);
}
template <class T>
int TBasicArray2D<T> :: operator == (const TBasicArray2D<T> &A) const
{
  if(ni!=A.row() || nj!=A.column() || m_Err) return 0;
  for(int i=0; i <ni; i++)
    for(int j=0; j<nj; j++)
      if(MatrixItems[i][j]!=A[i][j])return 0;
  return 1;
}
template <class T>
int TBasicArray2D<T> :: operator != (const TBasicArray2D<T> &A) const
{
  if(ni!=A.row() || nj!=A.column() || m_Err) return 1;
  for(int i=0; i <ni; i++)
    for(int j=0; j<nj; j++)
      if(MatrixItems[i][j]!=A[i][j])return 1;
  return 0;
}
template <class T>
TBasicArray2D<T>& TBasicArray2D<T> :: operator = (const TBasicArray2D<T> &A){
  Assert(!m_Err);
  Assert((ni==A.row()) && (nj==A.column()));
  for (int i=0; i< ni; i++)
    for (int j=0; j< nj; j++)
      MatrixItems[i][j] = A.MatrixItems[i][j];
  return *this;
}
/****************************************************   TVector    *****************************************************/
template <class Y, class Z> TVector<Z> operator * (const Y &a, const TVector<Z>& B){
  Assert(!B.IsErr());
  TVector<Z> Temp(B.size());
  for (int i =0; i< B.size(); i++)
    Temp[i] = ((Z) a)*B[i];
  return Temp;
    }
template <class Y, class Z> TVector<Z> operator + (const Y &a, const TVector<Z>& B){
    Assert(!B.IsErr());  
  TVector<Z> Temp(B.size());
      for (int i =0; i< B.size(); i++)
          Temp[i] = ((Z) a) + B[i];
      return Temp;
    }
template <class Y, class Z> TVector<Z> operator - (const Y &a, const TVector<Z>& B){
    Assert(!B.IsErr());  
  TVector<Z> Temp(B.size());
      for (int i =0; i< B.size(); i++)
      Temp[i] = ((Z) a) - B[i];
      return Temp;
    }

template <class T>
int TVector<T> :: operator == (const TVector<T> &A) const
{
  if(nn!=A.size() || m_Err) return 0;
  for(int i=0; i < nn; i++)
    if(VectorItems[i]!=A[i])return 0;
  return 1;
}
template <class T>
int TVector<T> :: operator != (const TVector<T> &A) const
{
  if(nn!=A.size() || m_Err) return 1;
  for(int i=0; i < nn; i++)
    if(VectorItems[i]!=A[i])return 1;
  return 0;
}
template <class T>
TVector<T>& TVector<T>::operator = (const TVector<T>& A){
 // Assert(!m_Err);
 // Assert(A.size()==nn);
  for (int i=0; i<nn; i++)
    VectorItems[i] = A[i];
  return *this;  
}
template <class T>
TVector<T> TVector<T>::operator + (const TVector<T>& A) const {
  Assert(!m_Err);
  Assert(nn==A.size());
  TVector<T> Temp(nn);
  for (int i=0; i< nn; i++)
    Temp[i] = A[i]+VectorItems[i];
  return Temp;
}
template <class T>
TVector<T> TVector<T>::operator - (const TVector<T>& A) const {
  Assert(!m_Err);
  Assert(nn==A.size());
  TVector<T> Temp(nn);
  for (int i=0; i< nn; i++)
    Temp[i] = VectorItems[i]-A[i];
  return Temp;
}
template <class T>
TVector<T> TVector<T>::operator * (const TMatrix<T>& B) const {
    //  TVector x  TMatrix
    Assert(!m_Err);
    Assert(!B.IsErr());
    Assert(nn==B.row());
  TVector<T> Temp(B.column());

  for (int i =0; i< Temp.size(); i++)
      Temp[i]=0.;
  for (int i =0; i< B.row(); i++)
      for(int j=0; j< B.column(); j++)
      Temp[j] += VectorItems[i]*B[i][j];
  return Temp;
}
template <class T>
TVector<T>& TVector<T>::operator += (const TVector<T>& A){
  Assert(!m_Err);
  Assert(nn==A.size());
  for (int i=0; i<nn; i++)
    VectorItems[i] += A[i];
  return *this;
}
template <class T>
TVector<T>& TVector<T>::operator -= (const TVector<T>& A){
  Assert(!m_Err);
  Assert(nn==A.size());
  for (int i=0; i<nn; i++)
    VectorItems[i] -= A[i];
  return *this;
}
template <class T>
T TVector<T>::MultScalar (const TVector<T>& B) const {
  Assert((!m_Err)&&(!B.IsErr()));
  Assert(nn==B.size());
  T res = 0.;
  for (int i =0; i< nn; i++)
    res += VectorItems[i]*B[i];
  return res;
}
template <class T>
TMatrix<T> TVector<T>::MultMatrix(const TVector<T>& A) const {
  Assert((!m_Err)&&(!A.IsErr()));
  TMatrix<T> Temp(nn,A.size());
  for(int i=0; i<nn; i++)
    for(int j=0; j<A.size(); j++)
      Temp[i][j]= VectorItems[i]*A[j];
  return Temp;
}
template <class T>
T TVector<T> :: Norm() const {
  Assert(!m_Err);
  T res = 0.;
  for (int i =0; i< nn; i++)
    res += VectorItems[i]*VectorItems[i];
  return sqrt(res);
}
template <class T>
void TVector<T> :: Normalize(){
  Assert(!m_Err);
  T norm = Norm2();
  for (int i =0; i< nn; i++)
    VectorItems[i] = VectorItems[i]/norm;
}
template <class T>
T TVector<T> :: Max() const {
  Assert(!m_Err);
  T max = VectorItems[0];
  for (int i =1; i< nn; i++)
    if(max < VectorItems[i]) max = VectorItems[i];
  return max;
}
template <class T>
T TVector<T> :: Min() const {
  Assert(!m_Err);
  T min = VectorItems[0];
  for (int i =1; i< nn; i++)
    if(min > VectorItems[i]) min = VectorItems[i];
  return min;
}
template <class T>
void TVector<T> :: Print (const char * szFile) const {
  //Assert(!m_Err);
  //MCFile file(szFile,_T("wt"));
  //for (int i=0; i < nn; i++)
  //  _ftprintf(file, _T("%20.10f \n"),VectorItems[i]);  //  !!! 20.10f
}
template <class T>
void TVector<T> :: Probe () const {
  Assert(!m_Err);
  long ltime = time(NULL);
  int stime = (unsigned) ltime/2;
  srand(stime); //srand(time(NULL));
  for(int i=0.; i< nn; i++)
    VectorItems[i]= (T) rand(); // % 100;
}
/****************************************************   TMatrix    *****************************************************/
template <class Y,class Z> TMatrix<Z> operator * (const Y &a, const TMatrix<Z> &B){
  Assert(!B.IsErr());
      TMatrix<Z> Temp(B.row(),B.column());
      for (int i=0; i< B.row(); i++)
          for(int j=0; j< B.column(); j++)
              Temp[i][j] = ((Z) a)*B[i][j];
      return Temp;
    }
template <class Y,class Z> TMatrix<Z> operator + (const Y &a, const TMatrix<Z> &B){
  Assert(!B.IsErr());
      TMatrix<Z> Temp(B.row(),B.column());
      for (int i=0; i< B.row(); i++)
          for(int j=0; j< B.column(); j++)
              Temp[i][j] = (Z) a + B[i][j];
      return Temp;
    }
template <class Y,class Z> TMatrix<Z> operator - (const Y &a, const TMatrix<Z> &B){
  Assert(!B.IsErr());
      TMatrix<Z> Temp(B.row(),B.column());
      for (int i=0; i< B.row(); i++)
          for(int j=0; j< B.column(); j++)
              Temp[i][j] = (Z) a - B[i][j];
      return Temp;
    }
template <class T>
TVector<T> TMatrix<T>::operator * (const TVector<T>& B) const {
  //  TMatrix x  TVector
  Assert(!m_Err);
  Assert(!B.IsErr());
  Assert(nj==B.size());
  TVector<T> Temp(ni);

  for (int i =0; i< Temp.size(); i++)
    Temp[i]=0.;
  for (int i =0; i< ni; i++)
    for(int j=0; j< nj; j++)
      Temp[i] += B[j]*MatrixItems[i][j];
  return Temp;
}
template <class T>
int TMatrix<T> :: operator == (const TMatrix<T> &A) const
{
  if(ni!=A.row() || nj!=A.column() || m_Err) return 0;
  for(int i=0; i <ni; i++)
    for(int j=0; j<nj; j++)
      if(MatrixItems[i][j]!=A[i][j])return 0;
  return 1;
}
template <class T>
int TMatrix<T> :: operator != (const TMatrix<T> &A) const
{
  if(ni!=A.row() || nj!=A.column() || m_Err) return 1;
  for(int i=0; i <ni; i++)
    for(int j=0; j<nj; j++)
      if(MatrixItems[i][j]!=A[i][j])return 1;
  return 0;
}
template <class T>
TMatrix<T>& TMatrix<T>::operator = (const TMatrix<T>& A){
 // Assert(!m_Err);
 // Assert((A.row()==ni)&&(A.column()==nj));
  for (int i=0; i< ni; i++)
    for (int j=0; j< nj; j++)
      MatrixItems[i][j] = A[i][j];
  return *this;
}
template <class T>
TMatrix<T> TMatrix<T>::operator * (const TMatrix<T> &B) const 
{
  //Assert(!m_Err);
 // Assert(nj==B.row());
 
  TMatrix<T> Temporary(ni,B.column());
  
  for (int i =0; i< ni; i++)
    for (int j =0; j< B.column(); j++)
      Temporary[i][j] =0.0;

  for (int i = 0; i< ni; i++)
    for (int j = 0; j< B.column(); j++)
      for (int l = 0; l< nj; l++)
        Temporary[i][j] += MatrixItems[i][l]*B[l][j];
  return Temporary;
}
template <class T>
TMatrix<T>& TMatrix<T>:: operator += (const TMatrix<T> &A)
{
  Assert(!m_Err);
  Assert((ni==A.row()) && (nj==A.column()));
  for (int i=0; i< ni; i++)
    for (int j=0; j< nj; j++)
      MatrixItems[i][j] += A[i][j];
    return *this;
}
template <class T> 
TMatrix<T>& TMatrix<T>:: operator -= (const TMatrix<T> &A)
{
  Assert(!m_Err);
  Assert((ni==A.row()) && (nj==A.column()));
  for (int i=0; i< ni; i++)
    for (int j=0; j< nj; j++)
      MatrixItems[i][j] -= A[i][j];
    return *this;
}
template <class T>
TMatrix<T> TMatrix<T>:: operator + (const TMatrix<T> &A) const 
{
  Assert(!m_Err);
  Assert((ni==A.row()) && (nj==A.column()));
  TMatrix<T> Temp(ni,nj);
  for (int i = 0; i< ni; i++)
    for (int j = 0; j< nj; j++)
    Temp[i][j] = A[i][j]+MatrixItems[i][j];
  return Temp;
}
template <class T>
TMatrix<T> TMatrix<T>:: operator - (const TMatrix<T> &A) const
{
  Assert(!m_Err);
  Assert((ni==A.row()) && (nj==A.column()));
  TMatrix<T> Temp(ni,nj);
  for (int i = 0; i< ni; i++)
    for (int j = 0; j< nj; j++)
    Temp[i][j] = MatrixItems[i][j]-A[i][j];
  return Temp;
}
template <class T>
TMatrix<T> TMatrix<T>:: operator - () const 
{
  Assert(!m_Err);
  TMatrix<T> Temp(ni,nj);
  for (int i = 0; i< ni; i++)
    for (int j = 0; j< nj; j++)
      Temp[i][j] = -MatrixItems[i][j];
  return Temp;
}
template <class T>
TVector<T> TMatrix<T>::row(const int i){
  Assert(!m_Err);
  TVector<T> Temp(nj);
  for(int j=0; j<nj; j++)
    Temp[j] = MatrixItems[i][j];
  return Temp;
}
template <class T>
TVector<T> TMatrix<T>::column(const int j){
  Assert(!m_Err);
  TVector<T> Temp(ni);
  for(int i=0; i<ni; i++)
    Temp[i] = MatrixItems[i][j];
  return Temp;
}
template <class T>
T TMatrix<T> :: Max() const 
{
  Assert(!m_Err);
  T max = MatrixItems[0][0];
  for (int i = 0; i< ni; i++)
    for (int j = 0; j< nj; j++)
      if(max<MatrixItems[i][j]) max = MatrixItems[i][j];
  return max;
}
template <class T>
T TMatrix<T> :: Min() const 
{
  Assert(!m_Err);
  T min = MatrixItems[0][0];
  for (int i = 0; i< ni; i++)
    for (int j = 0; j< nj; j++)
      if(min>MatrixItems[i][j]) min = MatrixItems[i][j];
  return min;
}
template <class T>
void TMatrix<T> :: Print (const char * szFile) const 
{
  //Assert(!m_Err);
  //MCFile file(szFile,_T("wt"));
  //for (int i=0; i < ni; i++)  {
  //  for (int j=0; j < nj; j++)  {
  //    _ftprintf(file, _T("%20.10f   "),MatrixItems[i][j]);  //  !!!  20.10f     
  //  }
  //_ftprintf(file,_T("\n"));
  //}
}

template <class T>
TMatrix<T> TMatrix<T> :: Transpose() const 
{
    Assert(!m_Err);
    TMatrix<T> Temp(nj,ni);
    for(int i=0; i< ni; i++)
        for(int j=0; j< nj; j++)
            Temp[i][j]=MatrixItems[j][i];
    return Temp;
}
template <class T>
T TMatrix<T> :: Trace() const {
  int size = ni;
  T sum = 0.;
  if(size > nj)size = nj;
  for(int i=0; i<size; i++)
    sum += MatrixItems[i][i];
  return sum;
}
template <class T>
int TMatrix<T> :: IsSymmetric() const {
  // Assert(!m_Err);
    for(int i=0; i< ni; i++)
        for(int j=0; j< nj; j++)
            if(MatrixItems[i][j]!=MatrixItems[j][i])return 0;
    return 1;
}

template <class T>
void TMatrix<T> :: Multiplication_and_Equating(const TMatrix<T>& A,const TMatrix<T>& B)
{
//  Assert(!m_Err);
//  Assert(ni == A.ni && nj==B.nj && A.nj == B.ni);
  int i,j,l;
  int A_i = A.row();
  int A_j = A.column();
  //int B_i = B.row();
  int B_j = B.column();


  for (i = 0;i<A_i;i++)
    for (j = 0;j<B_j;j++)    
      MatrixItems[i][j] =0.0;
       
  for (i = 0;i<A_i;i++)
    for (j = 0;j<B_j;j++)
      for (l = 0;l<A_j;l++)
        MatrixItems[i][j] = MatrixItems[i][j]+A[i][l]*B[l][j];
}

template <class T>
void TMatrix<T> :: Inverse()
{
    Assert(!m_Err);
 int i,j,k,l,N,N2;
 T det=1.0;
 T temp;
 N = this->ni;
 N2 = this->nj;
 if (N!= N2)
   return;
 T** a;
  a = (T**) calloc (N,sizeof(T *));
  for (i=0; i<N; i++) {
    a[i]= (T *)calloc (2*N,sizeof(T));}

 bool change;

 for(i=0;i<N;i++)
 {
   for(j=0;j<N;j++)
   {
     a[i][j]=this->MatrixItems[i][j];
     if(i==j) a[i][N+j]=1.0;
     else a[i][N+j]=0.0;
     
   }
 }

  for(i=0;i < N-1;i++)
  {
    if(fabs(a[i][i])<0.0000001)
    {
      change =false;
      for (k=i+1;k<=N-1;k++)
        {
          if(fabs(a[k][i])>0.0000001)
          {
            for (l=i;l<=2*N-1;l++)
            {
             temp= a[k][l];
               a[k][l]=a[i][l];
             a[i][l]=temp; 
             change=true;
            }
            if (change)
              goto break_label;

         }
        }

    }
    break_label:
     
    for(l=i+1;l<=2*N-1;l++)
    {

    a[i][l]=a[i][l]/a[i][i];
  
    } 

    a[i][i]=1;
   
    for(k=i+1;k < N;k++)
    {
      for(j=i+1;j <2*N;j++)
      {

          a[k][j]=a[k][j]-a[k][i]*a[i][j];

      }
      a[k][i]=0;
    }

  }
       
   for(l=N;l<=2*N-1;l++)
    a[N-1][l]/=a[N-1][N-1];
      // if(fabs(a[N-1][N-1])<0.0000001)
      //std::cout<<a[i][i]<<i<<std::endl;
     a[N-1][N-1]=1.0;
   
 for(i=N-2;i >=0; i--)
 {
  for(j=N-1;j!=i;j--)
  {
    for(l=j+1;l<=2*N-1;l++)
    {
    a[i][l]-=a[i][j]*a[j][l];
    }
    a[i][j]=0.0;
  }
 }

  for (i= 0;i<=N-1;i++)
 {
   for (j=0;j<=N-1;j++)
   {
    this->MatrixItems[i][j]=a[i][j+N];
   }
 }
 
  for ( i=0; i<N; i++)  {
    free (a[i]);
  }
  free (a);

 return;
}
template <class T>
void TMatrix<T> :: Probe () const 
{
  Assert(!m_Err);
  long ltime = time(NULL);
  int stime = (unsigned) ltime/2;
  srand(stime); //srand(time(NULL));
  for(int i=0.; i< ni; i++){
    for(int j=0; j< nj; j++){
      MatrixItems[i][j]= (T) rand(); // % 100;
    }
  }
}

  // Function to Inverse Matrix by it's LU Decomposition
template <class T>
TMatrix<T> MatrixInverse(const TMatrix<T> &ptMatrix){
  // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
  Assert(!ptMatrix.IsErr());
  int pivsign;
  int ni = ptMatrix.row();
  int nj = ptMatrix.column();
  TVector<T> Pivot_Vector(ni);
  TMatrix<T> LU_Matrix(ptMatrix);
  TMatrix<T> ResultMatrix(ni,nj);
  TMatrix<T> ImpactMatrix(ni,nj,0.);

  for (int i=0; i<ni; i++)
    ImpactMatrix[i][i] = 1;

  for (int i = 0; i < ni; i++) {
    Pivot_Vector[i] = i;
  }
  pivsign = 1;
  T *LUrowi = 0;
  TVector<T> LUcolj(ni);

  // Outer loop.
  for (int j = 0; j < nj; j++) {
    // Make a copy of the j-th column to localize references.
    for (int i = 0; i < ni; i++) {
      LUcolj[i] = LU_Matrix[i][j];
    }

    // Apply previous transformations.
    for (int i = 0; i < ni; i++) {
      LUrowi = LU_Matrix[i];
      // Most of the time is spent in the following dot product.
      int kmax = min(i,j);
      double s = 0.0;
      for (int k = 0; k < kmax; k++) {
        s += LUrowi[k]*LUcolj[k];
      }
      LUrowi[j] = LUcolj[i] -= s;
    }

    // Find pivot and exchange if necessary.
    int p = j;
    for (int i = j+1; i < ni; i++) {
      if (abs(LUcolj[i]) > abs(LUcolj[p])) {
        p = i;
      }
    }

    if (p != j) {
      int k=0;
      for (k = 0; k < nj; k++) {
        double t = LU_Matrix[p][k]; 
        LU_Matrix[p][k] = LU_Matrix[j][k]; 
        LU_Matrix[j][k] = t;
      }
      k = Pivot_Vector[p]; 
      Pivot_Vector[p] = Pivot_Vector[j]; 
      Pivot_Vector[j] = k;
      pivsign = -pivsign;
    }

    // Compute multipliers.
    if ((j < ni) && (LU_Matrix[j][j] != 0.0))
      for (int i = j+1; i < ni; i++)
        LU_Matrix[i][j] /= LU_Matrix[j][j];
  }

  // Copy right hand side with pivoting
  for (int i = 0; i < ni; i++) 
    for (int j = 0; j < ni; j++) 
      ResultMatrix[i][j] = ImpactMatrix[Pivot_Vector[i]][j];

  // Solve L*Y = B(piv,:)
  for (int k = 0; k < nj; k++) {
    for (int i = k+1; i < nj; i++) {
      for (int j = 0; j < ni; j++) {
        ResultMatrix[i][j] -= ResultMatrix[k][j]*LU_Matrix[i][k];
      }
    }
  }

  // Solve U*X = Y;
  for (int k = nj-1; k >= 0; k--) {
    for (int j = 0; j < ni; j++) {
      ResultMatrix[k][j] /= LU_Matrix[k][k];
    }
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < ni; j++) {
        ResultMatrix[i][j] -= ResultMatrix[k][j]*LU_Matrix[i][k];
      }
    }
  }
  return ResultMatrix;
}

class MCGaussLinearSystem
{
public:
  MCGaussLinearSystem  (int Num);
  ~MCGaussLinearSystem ();
  int m_nDim;                    // Order of matrix

  MCMatrix m_matrix;           // Line systm's matrix
  MCAPtr<double> ptResult;       // Vector of Results

  int  Setup (MCMatrix * ptMatrix, double * ptVector);
  int  Solve();
  void GetSolution(double * ptResultArray);
  void Print(const char * szFile);
 
private:
  MCGaussLinearSystem (const MCGaussLinearSystem &);
  //void operator=(const MCGaussLinearSystem &);
};

class MCLinearSystem
{
protected:
  // Order of LinearSystem (Matrix)
  int ni;    // 'num' of the equations
  int nj;    // 'num' of the variables
  //int Err;    // exist ptMatrix and ptImpact
public:
  int row()
  { return ni; }
  int column()
  { return nj; }
 // int IsErr() { return Err; }

  MCLinearSystem() : ni(0),nj(0)  { }
  explicit MCLinearSystem(const int,const int);    //  ptMatrix (ni x nj),  ptImpact (nj)
  MCLinearSystem(const MCMatrix&,const MCVector&);  //  (ptMatrix, ptImpact)
  MCLinearSystem(const MCMatrix&);  //  (ptMatrix)
  MCLinearSystem(const MCMatrix&,const MCMatrix&);  // (ptMatrix, ImpactMatrix)
  virtual ~MCLinearSystem () { }

  MCMatrix ptMatrix;       // Line systm's matrix
  MCVector ptImpact;       // Vector of Right-hand part of the linear system
  MCVector ptResult;       // Vector of Results
  
  // Matrix equations  A*X=B
  // A is mxn, X is nxk, B is mxk 
  MCMatrix ImpactMatrix;
  MCMatrix ResultMatrix;

  void Initialize (const MCMatrix&);  // ptMatrix
  void Initialize (const MCVector&);  // ptImpact
  void Initialize (const MCMatrix&,const MCVector&);  // (ptMatrix,ptImpact)
  void Initialize (const MCMatrix&,const MCMatrix&);  // (ptMatrix,ImpactMatrix)
  void Print_Matrix(const char * szFile);
  void Print_Impact(const char * szFile);
  void Print_Result(const char * szFile);
};

class MC_LU : public MCLinearSystem
{
/*  For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n unit lower triangular matrix L, an n-by-n upper
  triangular matrix U, and a permutation vector piv of length m so that A(piv,:) = L*U. If m < n, then L is m-by-m and U is
  m-by-n.
  The LU decompostion with pivoting always exists, even if the matrix is singular, so the constructor will never fail. The
  primary use of the LU decomposition is in the solution of square systems of simultaneous linear equations. This will fail if
  isNonsingular() returns false.
    JAMA/C++ Linear Algebra Library    */

  // for MCInteraction
  bool flag;  
  int index;  
  int infoID;

  int isDecompose;  // has LU_Decomposition been done?
  int pivsign;
  MCVector Pivot_Vector;
  MCMatrix permute_copy(const MCMatrix&,const MCVector&,const int,const int);
  MCVector permute_copy(const MCVector&,const MCVector&);
public:
  MC_LU() : flag(0), MCLinearSystem() { }
  explicit MC_LU (const int,const int);
  MC_LU(const MCMatrix&,const MCVector&);
  MC_LU(const MCMatrix&);
  MC_LU(const MCMatrix&,const MCMatrix&);
  ~MC_LU () { }

  MCMatrix LU_Matrix;    // internal storage of decomposition

  void LU_Decomposition (const int _infoID = NULL);
  void LU_Solve();
  void LU_SolveMatrix(const int _infoID = NULL);
  int isSingular(); // 1 - Singular
  double LU_Determinate();
  MCMatrix Get_U();
  MCMatrix Get_L();
  MCVector Get_Pivot()
  { return Pivot_Vector; }
};

class MC_QR : public MCLinearSystem
{
/*  For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
    orthogonal matrix Q and an n-by-n upper triangular matrix R so that
    A = Q*R.
    The QR decompostion always exists, even if the matrix does not have
    full rank, so the constructor will never fail.  The primary use of the
    QR decomposition is in the least squares solution of nonsquare systems
    of simultaneous linear equations.  This will fail if isFullRank()
    returns false.
    JAMA/C++ Linear Algebra Library    */
   
  // for MCInteraction
  bool flag;
  int index;
  int infoID;

  MCVector R_Diag_Vector;  // Diagonal Vector of R matrix
public:
  MC_QR() : flag(0), MCLinearSystem() { }
  explicit MC_QR (const int,const int);
  MC_QR(const MCMatrix&,const MCVector&);
  MC_QR(const MCMatrix&);
  MC_QR(const MCMatrix&,const MCMatrix&);
  ~MC_QR () { }

  MCMatrix QR_Matrix;  //  m x n matrix
 
  void QR_Decomposition (const int _infoID = NULL);
  void QR_Solve(const int _infoID = NULL);
  void QR_SolveMatrix(const int _infoID = NULL);
  int isFullRank(); // 1 - Full Rank
  MCMatrix Get_Q();
  MCMatrix Get_R();
  MCMatrix Get_Householder();
};

class MC_Cholesky : public MCLinearSystem  // Direct Method
{
/*    For a symmetric, positive definite matrix A, this function
   computes the Cholesky factorization, i.e. it computes a lower 
   triangular matrix L such that A = L*L'.
   If the matrix is not symmetric or positive definite, the function
   computes only a partial decomposition.  This can be tested with
   the is_spd() flag.

     Typical usage looks like:
    Array2D<double> A(n,n);
    Array2D<double> L;
     ... 
    Cholesky<double> chol(A);
    if (chol.is_spd())
      L = chol.getL();
      else
      cout << "factorization was not complete.\n";
  JAMA/C++ Linear Algebra Library    */

  // for MCInteraction
  bool flag;  
  int index;
  int infoID;
  int isspd;
public:
  MC_Cholesky() : isspd(0), MCLinearSystem() { }
  explicit MC_Cholesky (const int,const int,const int _spd=0);
  MC_Cholesky(const MCMatrix&,const MCVector&,const int _spd=0);
  MC_Cholesky(const MCMatrix&,const int _spd=0);
  MC_Cholesky(const MCMatrix&,const MCMatrix&,const int _spd=0);
  ~MC_Cholesky () { }

  MCMatrix Cholesky_Matrix;

  void Cholesky_Decomposition(const int _infoID = NULL);
  void Cholesky_Solve(const int _infoID = NULL);
  void Cholesky_SolveMatrix(const int _infoID = NULL);
  int isSPD()
  { return isspd; }
};

class MC_SVD
{
/* For an m-by-n matrix A with m >= n, the singular value decomposition is
   an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
   an n-by-n orthogonal matrix V so that A = U*S*V'.
   The singular values, sigma[k] = S[k][k], are ordered so that
   sigma[0] >= sigma[1] >= ... >= sigma[n-1].
   The singular value decompostion always exists, so the constructor will
   never fail.  The matrix condition number and the effective numerical
   rank can be computed from this decomposition.
    JAMA/C++ Linear Algebra Library    */

  int ni;
  int nj;
public:
  explicit MC_SVD  (const int,const int);
  MC_SVD(const MCMatrix&);
  virtual ~MC_SVD () { }

  MCMatrix ptMatrix;  //  Initial Matrix
  MCMatrix U_Matrix, V_Matrix; // Orthogonal matrixis
  MCVector S_Vector;           // Vector whis singular values

  void Initialize(const MCMatrix&);
  void SVD_Decomposition();
  MCMatrix Get_S();    // Get Matrix whis singular values on diagonal
  double Cond();
  int Rank();
};

class MC_Eigenvalue
{
/*  Computes eigenvalues and eigenvectors of a double (non-complex)
    matrix. 
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal. That is,
  the diagonal values of D are the eigenvalues, and
    V*V' = I, where I is the identity matrix.  The columns of V 
    represent the eigenvectors in the sense that A*V = V*D.
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the double eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    a + i*b, in 2-by-2 blocks, [a, b; -b, a].  That is, if the complex
    eigenvalues look like
          u + iv     .        .          .      .    .
            .      u - iv     .          .      .    .
            .        .      a + ib       .      .    .
            .        .        .        a - ib   .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
        then D looks like
            u        v        .          .      .    .
           -v        u        .          .      .    . 
            .        .        a          b      .    .
            .        .       -b          a      .    .
            .        .        .          .      x    .
            .        .        .          .      .    y
    This keeps V a double matrix in both symmetric and non-symmetric
    cases, and A*V = V*D.

  The matrix V may be badly conditioned, or even singular, so the 
  validity of the equation A = V*D*inverse(V) depends upon the 
  condition number of V.
    JAMA/C++ Linear Algebra Library    */
  
  int size;
  int Err;
  int issymmetric;
  double cdivr, cdivi;  // working variables
  MCVector ort;  // working storage for nonsymmetric algorithm
  void tred2();  // symmetric Householder reduction to tridiagonal form
  void tql2();  // symmetric tridiagonal QL algorithm
  void orthes();  // nonsymmetric reduction to Hessenberg form
  void cdiv(double, double, double, double);  // complex scalar division
  void hqr2();  // nonsymmetric reduction from Hessenberg to double Schur form
public:
  MC_Eigenvalue(const int);
  MC_Eigenvalue(MCMatrix&);
  virtual ~MC_Eigenvalue() { }

  MCMatrix ptMatrix;  // initial matrix

  MCVector Eigendouble;  // eigenvalues double part
  MCVector EigenImg;  // eigenvalues img part
  MCMatrix EigenMatrix;  // eigenvectors
  MCMatrix Hess;  // nonsymmetric Hessenberg form

  int IsSymmetric()
    { return ptMatrix.IsSymmetric(); }
    void Initialize(MCMatrix&);
  void EigenDecomposition();
  MCMatrix Get_D();    // Get block diagonal eigenvalue matrix
};
class Point3d {
public:
	float x;
	float y;
	float z;
	Point3d(){x = 0; y = 0; z = 0;}
	Point3d(float _x,float _y,float _z){ x= _x; y= _y; z = _z;}
	Point3d& operator += ( const Point3d & adder ) { x += adder.x; y += adder.y; z += adder.z; return *this; }
	Point3d& operator -= ( const Point3d & subber ) { x -= subber.x; y -= subber.y; z -= subber.z; return *this; }
	Point3d& operator *= ( float coeff ) { x *= coeff; y *= coeff; z *= coeff; return *this; }
	Point3d& operator /= ( float coeff );
	Point3d  operator -  () const { return Point3d(-x, -y, -z); }
	Point3d  operator *  ( float coeff ) const { return Point3d(x * coeff, y * coeff, z * coeff); }
	Point3d  operator /  ( float coeff ) const;
};
class MCInterpolation2D 
{
protected:
  double CalcNorma ();
  int PointNum;
public:
  MCInterpolation2D  (const int, const Point3d*);
  MCInterpolation2D  ();
  ~MCInterpolation2D () { }

  void Initialize (const int, const Point3d*);
 
  MCVector ptCoeff;
  MCAPtr <Point3d>  ptData; // Nodes and Function Values

  double m_Min;
  double m_Max; 
  double m_A;
  double Norm;
};

class MCInterpolation_RBF : public MCInterpolation2D
{
public:
  MCInterpolation_RBF(const int, const Point3d*,MCRadialBasis *);
  MCInterpolation_RBF();
  ~MCInterpolation_RBF() { }

  MCRadialBasis *ptRBF;

    //  for regulized surface
  //int m_NumPolinom_x;  //Order of Polinom
  //int m_NumPolinom_y;  //Order of Polinom
  //int isXY; // without xy x2y xy2 ...

  void Initialization(const int, const Point3d*,MCRadialBasis *);
  int Interpolate();
  double CalcValue(double x, double y);
  MCVector CalcSurface(int, int, double, double, double, double, double);
};

class MCRadialBasis 
{
  double m_smooth;
    int TypeRBF;    //enum RBF_FUNCTIONS
public:
  MCRadialBasis()
    { TypeRBF = 0;  m_smooth = 1.; }
  MCRadialBasis(const int _Type, const double _eps)
    { TypeRBF = _Type; m_smooth = _eps; }
  MCRadialBasis(const MCRadialBasis &);
  ~MCRadialBasis(){};

  void SetType(const int _Type)
    { TypeRBF = _Type; }
  void SetEps(const int _eps)
    { m_smooth = _eps; }

  double RBFuction(double);  // calculate RBF, determing by TypeRBF 
    double GaussRBF(double);
    double MultiquadricRBF(double);
    double InverseMultiquadricRBF(double);
  double InverseQuadricRBF(double);
  double TPS_RBF(double);
  double DistanceRBF(double);
    //double MaternLinerRBF(double);
    //double MaternSquareRBF(double);
    //double MacdonaldRBF(double);
    //double BiharmonicLaplaceRBF(double);
    //double BesselRBF(double,int);
};
