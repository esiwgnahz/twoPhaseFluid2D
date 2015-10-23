#include <time.h>
#include "mcmath.h"


//#include "math/deriv.h"
//#include "math/Bessel_func.h"
//
//#include "mclib/mcfile.h"
//#include "RVWnd/rvwnd.rh"
//
//#include "winutil/interact.h"

#define DUMMY_VALUE     -999 // А есть ведь уже INCATIVE_VALUE!!

// h_cur = x_{i} - x_{i-1}
// dy_cur = y_{i} - y_{i-1}
// calculate derivative in x1

#define APPR_EPS  1.E-6

/****************************************************   Derivatives   ***********************************************/
double DerivativeC (double y0 , double y1, double y2, double x0, double x1, double x2)
{
  double h_cur   = x1 - x0;
  double dy_cur  = y1 - y0;
  double h_next  = x2 - x1;
  double dy_next = y2 - y1;
  double dif = (dy_next*h_cur*h_cur + dy_cur*h_next*h_next)/(h_cur*h_next*(h_cur+h_next));
  return dif;
}

double DerivativeC(double y1, double y2, double h )
{
  double dP = y2 - y1;                     // delta p
  return dP/h;
}


double DerivativeL(double y0, double y1, double h)
{
  double dP = y1 - y0;                     // delta p
  return dP/h;
}

double DerivativeR(double y1, double y2, double h)
{
  double dP = y2 - y1;                     // delta p
  return dP/h;
}


//int porabol_approximation (double  *a, double *b, double *c, APPR_STRUCT * ptStruct)
//{   
//   double x1 = ptStruct->x1;
//   double x2 = ptStruct->x2;
//   double x3 = ptStruct->x3;
//
//   double y1 = ptStruct->y1;
//   double y2 = ptStruct->y2;
//   double y3 = ptStruct->y3;
//
//   double A=0;
//   double B=0;
//   double C=0;
//   
//   double dx12 = x1 - x2;
//   double dx23 = x2 - x3;
//
//  if (fabs(dx12) < APPR_EPS || fabs(dx23) < APPR_EPS ) 
//  return 0; 
//
//  A = ((y1-y2)/dx12 - (y1-y3)/(x1-x3))/(dx23);
//  B = (y1-y2)/dx12  - A*(x1+x2);
//  C = y1 - A*x1*x1 - B*x1;
//   
//  *a = A;
//  *b = B;
//  *c = C;
// 
//  return 1;
//}



/****************************************************   MCGaussLinearSystem   ***********************************************/
MCGaussLinearSystem :: MCGaussLinearSystem (int num)
            : m_nDim(num)
            , m_matrix (m_nDim +1, m_nDim+2)
            , ptResult (new double[m_nDim+1])
{
//  memset(&ptResult[0], 0, Ndim+1* sizeof(ptResult[0]));
}

MCGaussLinearSystem :: ~MCGaussLinearSystem ()
{
}

int MCGaussLinearSystem :: Setup (MCMatrix * ptMatrix, double * ptVector)
{
  MCAPtr<double> * A = m_matrix.pointer();  // internal matrix 
  MCAPtr<double> * B = ptMatrix->pointer();

  for (int i=1; i < m_matrix.row(); i++) {
     for (int j=1; j < m_matrix.column()-1; j++) {
       A[i][j] = B[i-1][j-1]; 
     }
  }
    for (int i=1; i < m_matrix.row(); i++) {
    A[i][m_matrix.column()-1] = ptVector[i-1];
  }
  return 1;
}

int MCGaussLinearSystem :: Solve ( )
{
  //double ** AA = m_matrix.pointer();
  MCAPtr<double> * AA = m_matrix.pointer();
  double * Res = ptResult;
  int   i,j,k,l;
  double cz[3];


  for (l=1; l <= m_nDim; l++) {
    cz[1] = 0.0;
    for (i=l; i <= m_nDim; i++) {
      if ( fabs(AA[i][l])>fabs(cz[1]))
      {
        k = i;
        cz[1] = AA[i][l];
      }
    }

    for (j=l; j <= m_nDim+1; j++) {
      if (fabs(cz[1]) < 1.e-8 ) {
        return 0;   // martix is degenerated
      }

      AA[k][j] = AA[k][j]/cz[1];
      if (k != l) {
        cz[2] = AA[k][j];
        AA[k][j] = AA[l][j];
        AA[l][j] = cz[2];
      }
    }

    for (i=l; i<= m_nDim; i++) {
      if (i != l) {
        cz[1] = AA[i][l];
        for (j=l+1;j<= m_nDim+1;j++)
          AA[i][j] = AA[i][j]-AA[l][j]*cz[1];
      }
    }

  }
  Res[m_nDim] = AA[m_nDim][m_nDim+1];
  for (l = m_nDim-1; l>0; l--) {
    Res[l] = AA[l][m_nDim+1];
    for (j = l+1; j < m_nDim+1; j++) {
      Res[l] = Res[l] - AA[l][j]*Res[j];
    }

  }
  return 1;
}

void MCGaussLinearSystem :: GetSolution (double * ptResultArray)
{
  for(int i=1; i < m_matrix.row(); i++) {
    ptResultArray[i-1] = ptResult[i];
  }
}

void MCGaussLinearSystem :: Print (const char * szFile)
{
  /*MCFile file(szFile,_T("wt"));

  for (int i=1; i < m_matrix.row(); i++) {
    _ftprintf(file,_T("%f8.3\n"),ptResult[i]);
  }
  _ftprintf(file,_T("\n"));*/
}
/****************************************************   MCLinearSystem   ***********************************************/
MCLinearSystem :: MCLinearSystem(const int _ni,const int _nj)
  : ni(_ni),
  nj(_nj),
  ptMatrix(_ni,_nj),
  ptImpact(_ni),
  ptResult(_nj)
{ }
MCLinearSystem :: MCLinearSystem(const MCMatrix &_Matrix,const MCVector &_Vector)
  : ni(_Matrix.row()),
  nj(_Matrix.column()),
  ptMatrix(_Matrix),
  ptImpact(_Vector),
  ptResult(_Vector.size())
{ }
MCLinearSystem :: MCLinearSystem(const MCMatrix &ptMatrix)
  : ni(ptMatrix.row()),
  nj(ptMatrix.column()),
  ptMatrix(ptMatrix),
  ptResult(ptMatrix.column())
{ }
MCLinearSystem :: MCLinearSystem(const MCMatrix &A,const MCMatrix &B)
  : ni(A.row()),
  nj(A.column()),
  ptMatrix(A),
  ImpactMatrix(B),
  ResultMatrix(A.column(),B.row())
{ }
void MCLinearSystem :: Initialize (const MCMatrix &_Matrix,const MCVector &_Vector){
    ptMatrix.Initialize(_Matrix);
    ptImpact.Initialize(_Vector);
}
void MCLinearSystem :: Initialize (const MCVector &_Vector){
    ptImpact.Initialize(_Vector);
}
void MCLinearSystem :: Initialize (const MCMatrix &_Matrix){
    ptMatrix.Initialize(_Matrix);
}
void MCLinearSystem :: Print_Matrix (const char * szFile)
{
    ptMatrix.Print(szFile);
}

void MCLinearSystem :: Print_Impact (const char* szFile)
{
    ptImpact.Print(szFile);
}

void MCLinearSystem :: Print_Result (const char * szFile)
{
    ptResult.Print(szFile);
}

/****************************************************   MC_LU_Decomposition   ***********************************************/
MCMatrix MC_LU::permute_copy(const MCMatrix &A,const MCVector &piv,
               const int j0,const int j1){
  MCMatrix X(piv.size(), j1-j0+1);
     for (int i = 0; i < piv.size(); i++) 
        for (int j = j0; j <= j1; j++) 
           X[i][j-j0] = A[piv[i]][j];
  return X;
}
MCVector MC_LU::permute_copy(const MCVector &A,const MCVector &piv){
  //Assert (piv.size() == A.size());
  MCVector x(piv.size());
     for (int i = 0; i < piv.size(); i++) 
           x[i] = A[piv[i]];
  return x;
}
MC_LU :: MC_LU (const int _ni, const int _nj)
  : LU_Matrix(_ni,_nj),
  Pivot_Vector(_ni),
  flag(0),
  MCLinearSystem(_ni,_nj)
{ }
MC_LU :: MC_LU(const MCMatrix& Matr,const MCVector& Vect)
  : LU_Matrix(Matr),
  Pivot_Vector(Matr.row()),
  flag(0),
  MCLinearSystem(Matr,Vect)
{ }
MC_LU :: MC_LU(const MCMatrix& Matr)
  : LU_Matrix(Matr),
  Pivot_Vector(Matr.row()),
  flag(0),
  MCLinearSystem(Matr)
{ }
MC_LU :: MC_LU(const MCMatrix &A,const MCMatrix &B)
  : LU_Matrix(A),
  Pivot_Vector(A.row()),
  flag(0),
  MCLinearSystem(A,B)
{ }
void MC_LU :: LU_Decomposition (const int _infoID)
{
  
  for (int i = 0; i < ni; i++) {
     Pivot_Vector[i] = i;
  }
  pivsign = 1;
  double *LUrowi = 0;
  MCVector LUcolj(ni);

  

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
        int kmax = std::min(i,j);
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
 
}

void MC_LU :: LU_Solve()
{
  /* Dimensions: A is mxn, X is nxk, B is mxk */
  //Assert(!Err);
 
 

  LU_Decomposition();

  //for (int i = 0; i < Pivot_Vector.size(); i++){ 
  //     ptResult[i] = ptImpact[Pivot_Vector[i]];
  //}
  ptResult = permute_copy(ptImpact,Pivot_Vector);


  // solve l*y = b(piv)
  for (int k = 0; k < nj; k++) {
     for (int i = k+1; i < nj; i++) {
           ptResult[i] -= ptResult[k]*LU_Matrix[i][k];
        }
     }
  

  // solve u*x = y;
  for (int k = nj-1; k >= 0; k--) {
        ptResult[k] /= LU_Matrix[k][k];
        for (int i = 0; i < k; i++){
          ptResult[i] -= ptResult[k]*LU_Matrix[i][k];
        }
  }
 
}
void MC_LU :: LU_SolveMatrix(const int _infoID)
{
  /* Dimensions: A is mxn, X is nxk, B is mxk */
  //Assert(!Err);
 // Assert(ImpactMatrix.row() == ni);
 // Assert(!isSingular());

  
  LU_Decomposition();

  // Copy right hand side with pivoting
  int nx = ImpactMatrix.column();
  ResultMatrix = permute_copy(ImpactMatrix,Pivot_Vector,0,nx-1);


  // Solve L*Y = B(piv,:)
  for (int k = 0; k < nj; k++) {
     for (int i = k+1; i < nj; i++) {
        for (int j = 0; j < nx; j++) {
           ResultMatrix[i][j] -= ResultMatrix[k][j]*LU_Matrix[i][k];
        }
     }
  }

  // Solve U*X = Y;
  for (int k = nj-1; k >= 0; k--) {
     for (int j = 0; j < nx; j++) {
        ResultMatrix[k][j] /= LU_Matrix[k][k];
     }
     for (int i = 0; i < k; i++) {
        for (int j = 0; j < nx; j++) {
           ResultMatrix[i][j] -= ResultMatrix[k][j]*LU_Matrix[i][k];
        }
     }
  }
 
}
MCMatrix MC_LU :: Get_U()
{
  //Assert(!Err);
  MCMatrix Temp(nj,nj);
  for (int i = 0; i < nj; i++) {
     for (int j = 0; j < nj; j++) {
        if (i <= j) {
           Temp[i][j] = LU_Matrix[i][j];
        } else {
           Temp[i][j] = 0.0;
        }
     }
  }
  return Temp;
}

MCMatrix MC_LU :: Get_L()
{
  //Assert(!Err);
  MCMatrix Temp(ni,nj);
  for (int i = 0; i < ni; i++) {
     for (int j = 0; j < nj; j++) {
        if (i > j) {
           Temp[i][j] = LU_Matrix[i][j];
        } else if (i == j) {
           Temp[i][j] = 1.0;
        } else {
           Temp[i][j] = 0.0;
        }
     }
  }
  return Temp;
}

double MC_LU :: LU_Determinate()
{
  //  Assert(ni==nj);
  double d = double(pivsign);
  for (int j = 0; j < ni; j++) {
     d *= LU_Matrix[j][j];
  }
  return d;
}
int MC_LU :: isSingular()
{
  //Assert(!Err);
  for (int j = 0; j < nj; j++) {
     if (LU_Matrix[j][j] == 0)
        return 1;
  }
  return 0;
}

/****************************************************   MC_QR_Decomposition   ***********************************************/
MC_QR :: MC_QR(const int _ni, const int _nj)
  : QR_Matrix(_ni,_nj),
  R_Diag_Vector(_nj),
  flag(0),
  MCLinearSystem(_ni,_nj)
{ }
MC_QR :: MC_QR(const MCMatrix& Matr,const MCVector& Vect)
  : QR_Matrix(Matr),
  R_Diag_Vector(Matr.column()),
  flag(0),
  MCLinearSystem(Matr,Vect)
{ }
MC_QR :: MC_QR(const MCMatrix& Matr)
  : QR_Matrix(Matr),
  R_Diag_Vector(Matr.column()),
  flag(0),
  MCLinearSystem(Matr)
{ }
MC_QR :: MC_QR(const MCMatrix &A,const MCMatrix &B)
  : QR_Matrix(A),
  R_Diag_Vector(A.column()),
  flag(0),
  MCLinearSystem(A,B)
{ }
void MC_QR :: QR_Decomposition (const int _infoID)
{    
  //Create a QR factorization object for QR_Matrix.
  //Assert(!Err);
 

  int i=0, j=0, k=0;
  // Main loop.
  for (k = 0; k < nj; k++) {
     // Compute 2-norm of k-th column without under/overflow.
     double nrm = 0;
     for (i = k; i < ni; i++) {
        nrm = hypot(nrm,QR_Matrix[i][k]);
     }

   if (nrm != 0.0) {
        // Form k-th Householder vector.
        if (QR_Matrix[k][k] < 0) {
           nrm = -nrm;
        }
        for (i = k; i < ni; i++) {
           QR_Matrix[i][k] /= nrm;
        }
        QR_Matrix[k][k] += 1.0;

        // Apply transformation to remaining columns.
        for (j = k+1; j < nj; j++) {
           double s = 0.0; 
           for (i = k; i < ni; i++) {
              s += QR_Matrix[i][k]*QR_Matrix[i][j];
           }
           s = -s/QR_Matrix[k][k];
           for (i = k; i < ni; i++) {
              QR_Matrix[i][j] += s*QR_Matrix[i][k];
           }
        }
     }
   
   R_Diag_Vector[k] = -nrm;
  }
 
}

void MC_QR :: QR_Solve(const int _infoID)
{
   /* Least squares solution of A*x = b
    B - m-length array (vector).
    x - n-length array (vector) that minimizes the two norm of Q*R*X-B.
  If B is non-conformant, or if QR.isFullRank() is false,
            the routine returns a null (0-length) vector.  */
  //Assert(!Err);
 // Assert(ptImpact.size() == ni);
  //Assert(isFullRank());

  
  MCVector x(ptImpact);
  QR_Decomposition ();

  // Compute Y = transpose(Q)*b
  for (int k = 0; k < nj; k++) 
  {
        double s = 0.0; 
        for (int i = k; i < ni; i++) 
    {
           s += QR_Matrix[i][k]*x[i];
        }
        s = -s/QR_Matrix[k][k];
        for (int i = k; i < ni; i++) 
    {
           x[i] += s*QR_Matrix[i][k];
        }
  }
 

  // Solve R*X = Y;
  for (int k = nj-1; k >= 0; k--){
     x[k] /= R_Diag_Vector[k];
     for (int i = 0; i < k; i++) {
           x[i] -= x[k]*QR_Matrix[i][k];
     }
  }
 

  for(int i=0; i<nj; i++)
    ptResult[i] = x[i];

 
}
void MC_QR :: QR_SolveMatrix(const int _infoID)
{
   /* Least squares solution of A*X = B
    B - m x k Array (must conform).
    X - n x k Array that minimizes the two norm of Q*R*X-B. If
       B is non-conformant, or if QR.isFullRank() is false,
            the routine returns a null (0x0) array.   */
  //Assert(!Err);
  //Assert(ImpactMatrix.row() == ni);
  //Assert(isFullRank());


  MCMatrix X(ImpactMatrix);
  int nx = ImpactMatrix.column(); 
  int i=0, j=0, k=0;

  QR_Decomposition ();
  // Compute Y = transpose(Q)*B
  for (k = 0; k < nj; k++) {
     for (j = 0; j < nx; j++) {
        double s = 0.0; 
        for (i = k; i < ni; i++) {
           s += QR_Matrix[i][k]*X[i][j];
        }
        s = -s/QR_Matrix[k][k];
        for (i = k; i < ni; i++) {
           X[i][j] += s*QR_Matrix[i][k];
        }
     }
  }


  // Solve R*X = Y;
  for (k = nj-1; k >= 0; k--) {
     for (j = 0; j < nx; j++) {
        X[k][j] /= R_Diag_Vector[k];
     }
     for (i = 0; i < k; i++) {
        for (j = 0; j < nx; j++) {
           X[i][j] -= X[k][j]*QR_Matrix[i][k];
        }
     }
  }
  

  for (i=0; i<nj; i++)
    for (j=0; j<nx; j++)
    ResultMatrix[i][j] = X[i][j];

 
}
MCMatrix MC_QR :: Get_R()
{
  MCMatrix Temp(nj,nj);
  for (int i = 0; i < nj; i++) {
     for (int j = 0; j < nj; j++) {
        if (i < j) {
           Temp[i][j] = QR_Matrix[i][j];
        } else if (i == j) {
           Temp[i][j] = R_Diag_Vector[i];
        } else {
           Temp[i][j] = 0.0;
        }
     }
  }
  return Temp;
}
MCMatrix MC_QR :: Get_Q()
{
  int i=0, j=0, k=0;

  MCMatrix Temp(ni,nj);
  for (k = nj-1; k >= 0; k--) {
     for (i = 0; i < ni; i++) {
        Temp[i][k] = 0.0;
     }
     Temp[k][k] = 1.0;
     for (j = k; j < nj; j++) {
        if (QR_Matrix[k][k] != 0) {
           double s = 0.0;
           for (i = k; i < ni; i++) {
              s += QR_Matrix[i][k]*Temp[i][j];
           }
           s = -s/QR_Matrix[k][k];
           for (i = k; i < ni; i++) {
              Temp[i][j] += s*QR_Matrix[i][k];
           }
        }
     }
  }
  return Temp;
}

MCMatrix MC_QR :: Get_Householder()
{
  MCMatrix Temp(ni,nj);
  for (int i = 0; i < ni; i++) 
  {
     for (int j = 0; j < nj; j++) 
   {
        if (i >= j) {
           Temp[i][j] = QR_Matrix[i][j];
        } else {
           Temp[i][j] = 0.0;
        }
     }
  }
  return Temp;
}
int MC_QR :: isFullRank()
{
  for (int j = 0; j < nj; j++) 
     if (R_Diag_Vector[j] == 0)
        return 0;
  return 1;
}
/****************************************************   MC_Cholesky   ***********************************************/
MC_Cholesky :: MC_Cholesky(const int _ni,const int _nj,const int _spd)
  : Cholesky_Matrix(_ni, _nj),
  isspd(_spd),
  flag(0),
  MCLinearSystem(_ni,_nj){
  //  Assert(_ni == _nj);  // is symmetric?
}
MC_Cholesky :: MC_Cholesky(const MCMatrix& Matr,const MCVector& Vect,const int _spd)
  : Cholesky_Matrix(Matr.row(),Matr.column()),
  MCLinearSystem(Matr,Vect),
  flag(0),
  isspd(_spd){
   // Assert(Matr.row() == Matr.column());  // is symmetric?
  //  Assert(Vect.size()==Matr.row());
}
MC_Cholesky :: MC_Cholesky(const MCMatrix& Matr,const int _spd)
  : Cholesky_Matrix(Matr),
  isspd(_spd),
  flag(0),
  MCLinearSystem(Matr){
  //  Assert(Matr.row() == Matr.column());  // is symmetric?
}
MC_Cholesky :: MC_Cholesky(const MCMatrix &A,const MCMatrix &B,const int _spd)
  : Cholesky_Matrix(A),
  isspd(_spd),
  flag(0),
  MCLinearSystem(A,B){
  //  Assert(A.row() == A.column());  // is symmetric?
}
void MC_Cholesky :: Cholesky_Decomposition(const int _infoID)
{
/*  Constructs a lower triangular matrix L, such that L*L'= A.
  If A is not symmetric positive-definite (SPD), only a
  partial factorization is performed.  If is_spd()
  evalutate true (1) then the factorizaiton was successful.  */
  //Assert(!Err);
  isspd = (ni == nj);
 //   Assert(ni == nj);

 

  {  double d,s;
      // Main loop.
     for (int j = 0; j < ni; j++){
         d= 0.0;
         for (int k = 0; k < j; k++){
            s= 0.0;
            for (int i = 0; i < k; i++){
               s += Cholesky_Matrix[k][i]*Cholesky_Matrix[j][i];
            }
            Cholesky_Matrix[j][k] = s = (ptMatrix[j][k] - s)/Cholesky_Matrix[k][k];
            d = d + s*s;
            isspd = isspd && (ptMatrix[k][j] == ptMatrix[j][k]); 
         }


         d = ptMatrix[j][j] - d;
         isspd = isspd && (d > 0.0);
         Cholesky_Matrix[j][j] = sqrt(d > 0.0 ? d : 0.0);
         for (int k = j+1; k < ni; k++){
            Cholesky_Matrix[j][k] = 0.0;
         }
     
  }
  }


}
void MC_Cholesky :: Cholesky_Solve(const int _infoID)
{
/*  Solve a linear system A*x = b, using the previously computed
  cholesky factorization of A: L*L'.
     B - A Matrix with as many rows as A and any number of columns.
     x - so that L*L'*x = b.
    If b is nonconformat, or if A was not symmetric posidtive
    definite, a null (0x0) array is returned.    */
  //Assert(!Err);
  //Assert(ptImpact.size() == ptMatrix.row());

 

  ptResult.Initialize(ptImpact);
  Cholesky_Decomposition();

  // Solve L*y = b;
  for (int k = 0; k < ni; k++) 
  {
     for (int i = 0; i < k; i++) 
           ptResult[k] -= ptResult[i]*Cholesky_Matrix[k][i];
   ptResult[k] /= Cholesky_Matrix[k][k];
  }
  

  // Solve L'*X = Y;
  for (int k = ni-1; k >= 0; k--) 
  {
     for (int i = k+1; i < ni; i++) 
           ptResult[k] -= ptResult[i]*Cholesky_Matrix[i][k];
     ptResult[k] /= Cholesky_Matrix[k][k];
  }

  
}

void MC_Cholesky :: Cholesky_SolveMatrix(const int _infoID)
{
/*  Solve a linear system A*X = B, using the previously computed
  cholesky factorization of A: L*L'.

     B - A Matrix with as many rows as A and any number of columns.
     X - so that L*L'*X = B.
      If B is nonconformat, or if A was not symmetric posidtive
      definite, a null (0x0) array is returned.    */

  //Assert(!Err);
  //Assert(ImpactMatrix.row() == ni);

  

  //MCMatrix X(ImpactMatrix);
  int nx = ImpactMatrix.column(); 
  Cholesky_Decomposition();
  ResultMatrix = ImpactMatrix;

#if 0
  // Solve L*Y = B;
  for (int k = 0; k < ni; k++) {
     for (int i = k+1; i < ni; i++) {
        for (int j = 0; j < nx; j++) {
           ResultMatrix[i][j] -= ResultMatrix[k][j]*Cholesky_Matrix[k][i];
        }
     }
     for (int j = 0; j < nx; j++) {
        ResultMatrix[k][j] /= Cholesky_Matrix[k][k];
   } 
  }
  // Solve L'*X = Y;
  for (int k = ni-1; k >= 0; k--) {
     for (int j = 0; j < nx; j++) {
        ResultMatrix[k][j] /= Cholesky_Matrix[k][k];
   }
     for (int i = 0; i < k; i++) {
        for (int j = 0; j < nx; j++) {
           ResultMatrix[i][j] -= ResultMatrix[k][j]*Cholesky_Matrix[k][i];
        }
     }
  }
#endif

#if 1
  // Solve L*y = b;
  for (int j=0; j< nx; j++)
  {
    for (int k = 0; k < ni; k++) 
  {
    for (int i = 0; i < k; i++) 
           ResultMatrix[k][j] -= ResultMatrix[i][j]*Cholesky_Matrix[k][i];
      ResultMatrix[k][j] /= Cholesky_Matrix[k][k];
   }
  }
 

  // Solve L'*X = Y;
 for (int j=0; j<nx; j++)
 {
    for (int k = ni-1; k >= 0; k--) 
    {
       for (int i = k+1; i < ni; i++) 
           ResultMatrix[k][j] -= ResultMatrix[i][j]*Cholesky_Matrix[i][k];
       ResultMatrix[k][j] /= Cholesky_Matrix[k][k];
  }
  }
#endif

 
}
/****************************************************   MC_SVD_Decomposition   ***********************************************/
MC_SVD :: MC_SVD(const int _ni, const int _nj)
  : ni(_ni),
  nj(_nj),
  ptMatrix(_ni,_nj),
  U_Matrix(_ni, _nj),
  V_Matrix(_ni, _nj),
  S_Vector(std::min(_ni+1,_nj))
{ }
MC_SVD :: MC_SVD(const MCMatrix& Matr)
  : ni(Matr.row()),
  nj(Matr.column()),
  ptMatrix(Matr),
  U_Matrix(Matr.row(),Matr.column()),
  V_Matrix(Matr.row(),Matr.column()),
  S_Vector(std::min(Matr.row()+1,Matr.column()))
{ }
void MC_SVD :: SVD_Decomposition()
{
  int nij = std::min(ni,nj);
  MCVector e(nj);
  MCVector work(ni);
  MCMatrix A(ptMatrix);
  int wantu = 1;
  int wantv = 1;
  int i=0, j=0, k=0;

  // Reduce A to bidiagonal form, storing the diagonal elements
  // in S_Vector and the super-diagonal elements in e.

  int nct = std::min(ni-1,nj);
  int nrt = std::max(0,std::min(nj-2,ni));
  for (k = 0; k < std::max(nct,nrt); k++) {
     if (k < nct) {

        // Compute the transformation for the k-th column and
        // place the k-th diagonal in S_Vector[k].
        // Compute 2-norm of k-th column without under/overflow.
        S_Vector[k] = 0;
        for (i = k; i < ni; i++) {
           S_Vector[k] = sqrt(S_Vector[k]*S_Vector[k] + A[i][k]*A[i][k]);
        }
        if (S_Vector[k] != 0.0) {
           if (A[k][k] < 0.0) {
              S_Vector[k] = -S_Vector[k];
           }
           for (i = k; i < ni; i++) {
              A[i][k] /= S_Vector[k];
           }
           A[k][k] += 1.0;
        }
        S_Vector[k] = -S_Vector[k];
     }
     for (j = k+1; j < nj; j++) {
        if ((k < nct) && (S_Vector[k] != 0.0))  {

        // Apply the transformation.

           double t(0.0);
           for (i = k; i < ni; i++) {
              t += A[i][k]*A[i][j];
           }
           t = -t/A[k][k];
           for (i = k; i < ni; i++) {
              A[i][j] += t*A[i][k];
           }
        }

        // Place the k-th row of A into e for the
        // subsequent calculation of the row transformation.

        e[j] = A[k][j];
     }
     if (wantu & (k < nct)) {

        // Place the transformation in U_Matrix for subsequent back
        // multiplication.

        for (i = k; i < ni; i++) {
           U_Matrix[i][k] = A[i][k];
        }
     }
     if (k < nrt) {

        // Compute the k-th row transformation and place the
        // k-th super-diagonal in e[k].
        // Compute 2-norm without under/overflow.
        e[k] = 0;
        for (i = k+1; i < nj; i++) {
           e[k] = sqrt(e[k]*e[k] + e[i]*e[i]);
        }
        if (e[k] != 0.0) {
           if (e[k+1] < 0.0) {
              e[k] = -e[k];
           }
           for (i = k+1; i < nj; i++) {
              e[i] /= e[k];
           }
           e[k+1] += 1.0;
        }
        e[k] = -e[k];
        if ((k+1 < ni) & (e[k] != 0.0)) {

        // Apply the transformation.

           for (i = k+1; i < ni; i++) {
              work[i] = 0.0;
           }
           for (j = k+1; j < nj; j++) {
              for (i = k+1; i < ni; i++) {
                 work[i] += e[j]*A[i][j];
              }
           }
           for (j = k+1; j < nj; j++) {
              double t(-e[j]/e[k+1]);
              for (i = k+1; i < ni; i++) {
                 A[i][j] += t*work[i];
              }
           }
        }
        if (wantv) {

        // Place the transformation in V_Matrix for subsequent
        // back multiplication.

           for (i = k+1; i < nj; i++) {
              V_Matrix[i][k] = e[i];
           }
        }
     }
  }

  // Set up the final bidiagonal matrix or order p.

  int p = std::min(nj,ni+1);
  if (nct < nj) {
     S_Vector[nct] = A[nct][nct];
  }
  if (ni < p) {
     S_Vector[p-1] = 0.0;
  }
  if (nrt+1 < p) {
     e[nrt] = A[nrt][p-1];
  }
  e[p-1] = 0.0;

  // If required, generate U_Matrix.

  if (wantu) {
     for (j = nct; j < nij; j++) {
        for (i = 0; i < ni; i++) {
           U_Matrix[i][j] = 0.0;
        }
        U_Matrix[j][j] = 1.0;
     }
     for (k = nct-1; k >= 0; k--) {
        if (S_Vector[k] != 0.0) {
           for (j = k+1; j < nij; j++) {
              double t(0.0);
              for (i = k; i < ni; i++) {
                 t += U_Matrix[i][k]*U_Matrix[i][j];
              }
              t = -t/U_Matrix[k][k];
              for (i = k; i < ni; i++) {
                 U_Matrix[i][j] += t*U_Matrix[i][k];
              }
           }
           for (i = k; i < ni; i++ ) {
              U_Matrix[i][k] = -U_Matrix[i][k];
           }
           U_Matrix[k][k] = 1.0 + U_Matrix[k][k];
           for (i = 0; i < k-1; i++) {
              U_Matrix[i][k] = 0.0;
           }
        } else {
           for (i = 0; i < ni; i++) {
              U_Matrix[i][k] = 0.0;
           }
           U_Matrix[k][k] = 1.0;
        }
     }
  }

  // If required, generate V_Matrix.

  if (wantv) {
     for (k = nj-1; k >= 0; k--) {
        if ((k < nrt) & (e[k] != 0.0)) {
           for (j = k+1; j < nij; j++) {
              double t(0.0);
              for (i = k+1; i < nj; i++) {
                 t += V_Matrix[i][k]*V_Matrix[i][j];
              }
              t = -t/V_Matrix[k+1][k];
              for (i = k+1; i < nj; i++) {
                 V_Matrix[i][j] += t*V_Matrix[i][k];
              }
           }
        }
        for (i = 0; i < nj; i++) {
           V_Matrix[i][k] = 0.0;
        }
        V_Matrix[k][k] = 1.0;
     }
  }

  // Main iteration loop for the singular values.

  int pp = p-1;
  int iter = 0;
  double eps(pow(2.0,-52.0));
  while (p > 0) {
     int k=0;
   int kase=0;

     // Here is where a test for too many iterations would go.

     // This section of the program inspects for
     // negligible elements in the S_Vector and e arrays.  On
     // completion the variables kase and k are set as follows.

     // kase = 1     if S_Vector(p) and e[k-1] are negligible and k<p
     // kase = 2     if S_Vector(k) is negligible and k<p
     // kase = 3     if e[k-1] is negligible, k<p, and
     //              S_Vector(k), ..., S_Vector(p) are not negligible (qr step).
     // kase = 4     if e(p-1) is negligible (convergence).

     for (k = p-2; k >= -1; k--) {
        if (k == -1) {
           break;
        }
        if (abs(e[k]) <= eps*(abs(S_Vector[k]) + abs(S_Vector[k+1]))) {
           e[k] = 0.0;
           break;
        }
     }
     if (k == p-2) {
        kase = 4;
     } else {
        int ks;
        for (ks = p-1; ks >= k; ks--) {
           if (ks == k) {
              break;
           }
           double t( (ks != p ? abs(e[ks]) : 0.) + 
                      (ks != k+1 ? abs(e[ks-1]) : 0.));
           if (abs(S_Vector[ks]) <= eps*t)  {
              S_Vector[ks] = 0.0;
              break;
           }
        }
        if (ks == k) {
           kase = 3;
        } else if (ks == p-1) {
           kase = 1;
        } else {
           kase = 2;
           k = ks;
        }
     }
     k++;

     // Perform the task indicated by kase.

     switch (kase) {

        // Deflate negligible S_Vector(p).

        case 1: {
           double f(e[p-2]);
           e[p-2] = 0.0;
           for (j = p-2; j >= k; j--) {
              double t = sqrt(S_Vector[j]*S_Vector[j] + f*f);
              double cs(S_Vector[j]/t);
              double sn(f/t);
              S_Vector[j] = t;
              if (j != k) {
                 f = -sn*e[j-1];
                 e[j-1] = cs*e[j-1];
              }
              if (wantv) {
                 for (i = 0; i < nj; i++) {
                    t = cs*V_Matrix[i][j] + sn*V_Matrix[i][p-1];
                    V_Matrix[i][p-1] = -sn*V_Matrix[i][j] + cs*V_Matrix[i][p-1];
                    V_Matrix[i][j] = t;
                 }
              }
           }
        }
        break;

        // Split at negligible S_Vector(k).

        case 2: {
           double f(e[k-1]);
           e[k-1] = 0.0;
           for (j = k; j < p; j++) {
              double t = sqrt(S_Vector[j]*S_Vector[j] + f*f);
              double cs( S_Vector[j]/t);
              double sn(f/t);
              S_Vector[j] = t;
              f = -sn*e[j];
              e[j] = cs*e[j];
              if (wantu) {
                 for (i = 0; i < ni; i++) {
                    t = cs*U_Matrix[i][j] + sn*U_Matrix[i][k-1];
                    U_Matrix[i][k-1] = -sn*U_Matrix[i][j] + cs*U_Matrix[i][k-1];
                    U_Matrix[i][j] = t;
                 }
              }
           }
        }
        break;

        // Perform one qr step.

        case 3: {

           // Calculate the shift.

           double scale = std::max(std::max(std::max(std::max(
                   abs(S_Vector[p-1]),abs(S_Vector[p-2])),abs(e[p-2])), 
                   abs(S_Vector[k])),abs(e[k]));
           double sp = S_Vector[p-1]/scale;
           double spm1 = S_Vector[p-2]/scale;
           double epm1 = e[p-2]/scale;
           double sk = S_Vector[k]/scale;
           double ek = e[k]/scale;
           double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
           double c = (sp*epm1)*(sp*epm1);
           double shift = 0.0;
           if ((b != 0.0) || (c != 0.0)) {
              shift = sqrt(b*b + c);
              if (b < 0.0) {
                 shift = -shift;
              }
              shift = c/(b + shift);
           }
           double f = (sk + sp)*(sk - sp) + shift;
           double g = sk*ek;

           // Chase zeros.

           for (j = k; j < p-1; j++) {
              double t = sqrt(g*g + f*f);
              double cs = f/t;
              double sn = g/t;
              if (j != k) {
                 e[j-1] = t;
              }
              f = cs*S_Vector[j] + sn*e[j];
              e[j] = cs*e[j] - sn*S_Vector[j];
              g = sn*S_Vector[j+1];
              S_Vector[j+1] = cs*S_Vector[j+1];
              if (wantv) {
                 for (i = 0; i < nj; i++) {
                    t = cs*V_Matrix[i][j] + sn*V_Matrix[i][j+1];
                    V_Matrix[i][j+1] = -sn*V_Matrix[i][j] + cs*V_Matrix[i][j+1];
                    V_Matrix[i][j] = t;
                 }
              }
              t = sqrt(g*g + f*f);
              cs = f/t;
              sn = g/t;
              S_Vector[j] = t;
              f = cs*e[j] + sn*S_Vector[j+1];
              S_Vector[j+1] = -sn*e[j] + cs*S_Vector[j+1];
              g = sn*e[j+1];
              e[j+1] = cs*e[j+1];
              if (wantu && (j < ni-1)) {
                 for (i = 0; i < ni; i++) {
                    t = cs*U_Matrix[i][j] + sn*U_Matrix[i][j+1];
                    U_Matrix[i][j+1] = -sn*U_Matrix[i][j] + cs*U_Matrix[i][j+1];
                    U_Matrix[i][j] = t;
                 }
              }
           }
           e[p-2] = f;
           iter = iter + 1;
        }
        break;

        // Convergence.

        case 4: {

           // Make the singular values positive.

           if (S_Vector[k] <= 0.0) {
              S_Vector[k] = (S_Vector[k] < 0.0 ? -S_Vector[k] : 0.0);
              if (wantv) {
                 for (i = 0; i <= pp; i++) {
                    V_Matrix[i][k] = -V_Matrix[i][k];
                 }
              }
           }

           // Order the singular values.

           while (k < pp) {
              if (S_Vector[k] >= S_Vector[k+1]) {
                 break;
              }
              double t = S_Vector[k];
              S_Vector[k] = S_Vector[k+1];
              S_Vector[k+1] = t;
              if (wantv && (k < nj-1)) {
                 for (i = 0; i < nj; i++) {
                    t = V_Matrix[i][k+1]; V_Matrix[i][k+1] = V_Matrix[i][k]; V_Matrix[i][k] = t;
                 }
              }
              if (wantu && (k < ni-1)) {
                 for (i = 0; i < ni; i++) {
                    t = U_Matrix[i][k+1]; U_Matrix[i][k+1] = U_Matrix[i][k]; U_Matrix[i][k] = t;
                 }
              }
              k++;
           }
           iter = 0;
           p--;
        }
        break;
     }
  }
}
MCMatrix MC_SVD :: Get_S()
{
  //Return the diagonal matrix of singular values
  MCMatrix Temp(ni,nj);
  for (int i = 0; i < ni; i++) {
     for (int j = 0; j < nj; j++) {
        Temp[i][j] = 0.0;
     }
     Temp[i][i] = S_Vector[i];
  }
  return Temp;
}
int MC_SVD::Rank(){
    //  Effective numerical matrix rank
  double eps = pow(2.0,-52.0);
  double tol = std::max(ni,nj)*S_Vector[0]*eps;
  int r = 0;
  for (int i = 0; i < S_Vector.size(); i++) {
     if (S_Vector[i] > tol) {
        r++;
     }
  }
  return r;
}
double MC_SVD::Cond(){
	
    //  Two norm of condition number (max(S)/min(S))
  return S_Vector[0]/S_Vector[std::min(ni,nj)-1];
}
void MC_SVD::Initialize(const MCMatrix& Matr){
  ni = Matr.row();
  nj = Matr.column();
  ptMatrix.Initialize(Matr);
  U_Matrix.Initialize(Matr.row(),Matr.column());
  V_Matrix.Initialize(Matr.row(),Matr.column());
  S_Vector.Initialize(std::min(Matr.row()+1,Matr.column()));
}
//****************************************************   MC_Eigenvalue    ****************************************************/
MC_Eigenvalue::MC_Eigenvalue(const int _nn)
  : size(_nn),
  Err(1),
  ptMatrix(_nn,_nn),
  EigenMatrix(_nn,_nn),
  Eigendouble(_nn),
  EigenImg(_nn)
{ }
MC_Eigenvalue::MC_Eigenvalue(MCMatrix &Matr)
  : size(Matr.row()),
  Err(0),
  ptMatrix(Matr),
  EigenMatrix(Matr.row(),Matr.row()),
  Eigendouble(Matr.row()),
  EigenImg(Matr.row()),
  issymmetric(Matr.IsSymmetric())
{ }
void MC_Eigenvalue::Initialize(MCMatrix &Matr){
  size = Matr.row();
  Err = 0;
  ptMatrix.Initialize(Matr);
  EigenMatrix.Initialize(Matr.row(),Matr.row());
  Eigendouble.Initialize(Matr.row());
  EigenImg.Initialize(Matr.row());
  issymmetric = Matr.IsSymmetric();
}
void MC_Eigenvalue::EigenDecomposition(){
 // Assert(!Err);  
  if (issymmetric) {
   EigenMatrix = ptMatrix;
     tred2(); // Tridiagonalize
     tql2(); // Diagonalize
  } else {
   Hess.Initialize(ptMatrix);
   ort.Initialize(ptMatrix.row());
     orthes(); // Reduce to Hessenberg form
     hqr2(); // Reduce Hessenberg to double Schur form
  }
}
MCMatrix MC_Eigenvalue::Get_D(){
/*  Computes the block diagonal eigenvalue matrix.
    If the original matrix A is not symmetric, then the eigenvalue 
  matrix D is block diagonal with the double eigenvalues in 1-by-1 
  blocks and any complex eigenvalues,
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
    cases, and A*V = V*D.    */

  MCMatrix Temp(size,size);
  for (int i = 0; i < size; i++) {
     for (int j = 0; j < size; j++) {
        Temp[i][j] = 0.0;
     }
     Temp[i][i] = Eigendouble[i];
     if (EigenImg[i] > 0) {
        Temp[i][i+1] = EigenImg[i];
     } else if (EigenImg[i] < 0) {
        Temp[i][i-1] = EigenImg[i];
     }
  }
  return Temp;
}
void MC_Eigenvalue::tred2(){
   //  This is derived from the Algol procedures tred2 by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.

  for (int j = 0; j < size; j++) {
     Eigendouble[j] = EigenMatrix[size-1][j];
  }

  // Householder reduction to tridiagonal form.
  for (int i = size-1; i > 0; i--) {
     // Scale to avoid under/overflow.
     double scale = 0.0;
     double h = 0.0;
     for (int k = 0; k < i; k++) {
        scale = scale + abs(Eigendouble[k]);
     }
     if (scale == 0.0) {
        EigenImg[i] = Eigendouble[i-1];
        for (int j = 0; j < i; j++) {
           Eigendouble[j] = EigenMatrix[i-1][j];
           EigenMatrix[i][j] = 0.0;
           EigenMatrix[j][i] = 0.0;
        }
     } else {
        // Generate Householder vector.
        for (int k = 0; k < i; k++) {
           Eigendouble[k] /= scale;
           h += Eigendouble[k] * Eigendouble[k];
        }
        double f = Eigendouble[i-1];
        double g = sqrt(h);
        if (f > 0) {
           g = -g;
        }
        EigenImg[i] = scale * g;
        h = h - f * g;
        Eigendouble[i-1] = f - g;
        for (int j = 0; j < i; j++) {
           EigenImg[j] = 0.0;
        }

        // Apply similarity transformation to remaining columns.
        for (int j = 0; j < i; j++) {
           f = Eigendouble[j];
           EigenMatrix[j][i] = f;
           g = EigenImg[j] + EigenMatrix[j][j] * f;
           for (int k = j+1; k <= i-1; k++) {
              g += EigenMatrix[k][j] * Eigendouble[k];
              EigenImg[k] += EigenMatrix[k][j] * f;
           }
           EigenImg[j] = g;
        }
        f = 0.0;
        for (int j = 0; j < i; j++) {
           EigenImg[j] /= h;
           f += EigenImg[j] * Eigendouble[j];
        }
        double hh = f / (h + h);
        for (int j = 0; j < i; j++) {
           EigenImg[j] -= hh * Eigendouble[j];
        }
        for (int j = 0; j < i; j++) {
           f = Eigendouble[j];
           g = EigenImg[j];
           for (int k = j; k <= i-1; k++) {
              EigenMatrix[k][j] -= (f * EigenImg[k] + g * Eigendouble[k]);
           }
           Eigendouble[j] = EigenMatrix[i-1][j];
           EigenMatrix[i][j] = 0.0;
        }
     }
     Eigendouble[i] = h;
  }

  // Accumulate transformations.
  for (int i = 0; i < size-1; i++) {
     EigenMatrix[size-1][i] = EigenMatrix[i][i];
     EigenMatrix[i][i] = 1.0;
     double h = Eigendouble[i+1];
     if (h != 0.0) {
        for (int k = 0; k <= i; k++) {
           Eigendouble[k] = EigenMatrix[k][i+1] / h;
        }
        for (int j = 0; j <= i; j++) {
           double g = 0.0;
           for (int k = 0; k <= i; k++) {
              g += EigenMatrix[k][i+1] * EigenMatrix[k][j];
           }
           for (int k = 0; k <= i; k++) {
              EigenMatrix[k][j] -= g * Eigendouble[k];
           }
        }
     }
     for (int k = 0; k <= i; k++) {
        EigenMatrix[k][i+1] = 0.0;
     }
  }
  for (int j = 0; j < size; j++) {
     Eigendouble[j] = EigenMatrix[size-1][j];
     EigenMatrix[size-1][j] = 0.0;
  }
  EigenMatrix[size-1][size-1] = 1.0;
  EigenImg[0] = 0.0;
}
void MC_Eigenvalue::tql2(){
   //  This is derived from the Algol procedures tql2, by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.

  for (int i = 1; i < size; i++) {
     EigenImg[i-1] = EigenImg[i];
  }
  EigenImg[size-1] = 0.0;

  double f = 0.0;
  double tst1 = 0.0;
  double eps = pow(2.0,-52.0);
  for (int l = 0; l < size; l++) {
     // Find small subdiagonal element
	  tst1 = std::max(tst1,abs(Eigendouble[l]) + abs(EigenImg[l]));
     int m = l;

    // Original while-loop from Java code
     while (m < size) {
        if (abs(EigenImg[m]) <= eps*tst1) {
           break;
        }
        m++;
     }

     // If m == l, Eigendouble[l] is an eigenvalue,
     // otherwise, iterate.
     if (m > l) {
        int iter = 0;
        do {
           iter = iter + 1;  // (Could check iteration count here.)

           // Compute implicit shift
           double g = Eigendouble[l];
           double p = (Eigendouble[l+1] - g) / (2.0 * EigenImg[l]);
           double r = hypot(p,1.0);
           if (p < 0) {
              r = -r;
           }
           Eigendouble[l] = EigenImg[l] / (p + r);
           Eigendouble[l+1] = EigenImg[l] * (p + r);
           double dl1 = Eigendouble[l+1];
           double h = g - Eigendouble[l];
           for (int i = l+2; i < size; i++) {
              Eigendouble[i] -= h;
           }
           f = f + h;

           // Implicit QL transformation.
           p = Eigendouble[m];
           double c = 1.0;
           double c2 = c;
           double c3 = c;
           double el1 = EigenImg[l+1];
           double s = 0.0;
           double s2 = 0.0;
           for (int i = m-1; i >= l; i--) {
              c3 = c2;
              c2 = c;
              s2 = s;
              g = c * EigenImg[i];
              h = c * p;
              r = hypot(p,EigenImg[i]);
              EigenImg[i+1] = s * r;
              s = EigenImg[i] / r;
              c = p / r;
              p = c * Eigendouble[i] - s * g;
              Eigendouble[i+1] = h + s * (c * g + s * Eigendouble[i]);

              // Accumulate transformation.
              for (int k = 0; k < size; k++) {
                 h = EigenMatrix[k][i+1];
                 EigenMatrix[k][i+1] = s * EigenMatrix[k][i] + c * h;
                 EigenMatrix[k][i] = c * EigenMatrix[k][i] - s * h;
              }
           }
           p = -s * s2 * c3 * el1 * EigenImg[l] / dl1;
           EigenImg[l] = s * p;
           Eigendouble[l] = c * p;

           // Check for convergence.
        } while (abs(EigenImg[l]) > eps*tst1);
     }
     Eigendouble[l] = Eigendouble[l] + f;
     EigenImg[l] = 0.0;
  }
 
  // Sort eigenvalues and corresponding vectors.
  for (int i = 0; i < size-1; i++) {
     int k = i;
     double p = Eigendouble[i];
     for (int j = i+1; j < size; j++) {
        if (Eigendouble[j] < p) {
           k = j;
           p = Eigendouble[j];
        }
     }
     if (k != i) {
        Eigendouble[k] = Eigendouble[i];
        Eigendouble[i] = p;
        for (int j = 0; j < size; j++) {
           p = EigenMatrix[j][i];
           EigenMatrix[j][i] = EigenMatrix[j][k];
           EigenMatrix[j][k] = p;
        }
     }
  }
}
void MC_Eigenvalue::orthes(){
  //  This is derived from the Algol procedures orthes and ortran,
  //  by Martin and Wilkinson, Handbook for Auto. Comp.,
  //  Vol.ii-Linear Algebra, and the corresponding
  //  Fortran subroutines in EISPACK.

  int low = 0;
  int high = size-1;

  for (int m = low+1; m <= high-1; m++) {
     // Scale column.
     double scale = 0.0;
     for (int i = m; i <= high; i++) {
        scale = scale + abs(Hess[i][m-1]);
     }
     if (scale != 0.0) {
        // Compute Householder transformation.
        double h = 0.0;
        for (int i = high; i >= m; i--) {
           ort[i] = Hess[i][m-1]/scale;
           h += ort[i] * ort[i];
        }
        double g = sqrt(h);
        if (ort[m] > 0) {
           g = -g;
        }
        h = h - ort[m] * g;
        ort[m] = ort[m] - g;

        // Apply Householder similarity transformation
        // Hess = (I-u*u'/h)*Hess*(I-u*u')/h)
        for (int j = m; j < size; j++) {
           double f = 0.0;
           for (int i = high; i >= m; i--) {
              f += ort[i]*Hess[i][j];
           }
           f = f/h;
           for (int i = m; i <= high; i++) {
              Hess[i][j] -= f*ort[i];
           }
       }

       for (int i = 0; i <= high; i++) {
           double f = 0.0;
           for (int j = high; j >= m; j--) {
              f += ort[j]*Hess[i][j];
           }
           f = f/h;
           for (int j = m; j <= high; j++) {
              Hess[i][j] -= f*ort[j];
           }
        }
        ort[m] = scale*ort[m];
        Hess[m][m-1] = scale*g;
     }
  }

  // Accumulate transformations (Algol's ortran).
  for (int i = 0; i < size; i++) {
     for (int j = 0; j < size; j++) {
        EigenMatrix[i][j] = (i == j ? 1.0 : 0.0);
     }
  }

  for (int m = high-1; m >= low+1; m--) {
     if (Hess[m][m-1] != 0.0) {
        for (int i = m+1; i <= high; i++) {
           ort[i] = Hess[i][m-1];
        }
        for (int j = m; j <= high; j++) {
           double g = 0.0;
           for (int i = m; i <= high; i++) {
              g += ort[i] * EigenMatrix[i][j];
           }
           // Double division avoids possible underflow
           g = (g / ort[m]) / Hess[m][m-1];
           for (int i = m; i <= high; i++) {
              EigenMatrix[i][j] += g * ort[i];
           }
        }
     }
  }
}
void MC_Eigenvalue::cdiv(double xr, double xi, double yr, double yi){
  double r,d;
  if (abs(yr) > abs(yi)) {
     r = yi/yr;
     d = yr + r*yi;
     cdivr = (xr + r*xi)/d;
     cdivi = (xi - r*xr)/d;
  } else {
     r = yr/yi;
     d = yi + r*yr;
     cdivr = (r*xr + xi)/d;
     cdivi = (r*xi - xr)/d;
  }
}
void MC_Eigenvalue::hqr2(){
  //  This is derived from the Algol procedure hqr2,
  //  by Martin and Wilkinson, Handbook for Auto. Comp.,
  //  Vol.ii-Linear Algebra, and the corresponding
  //  Fortran subroutine in EISPACK.
   
      // Initialize
  int nn = size;
  int n = nn-1;
  int low = 0;
  int high = nn-1;
  double eps = pow(2.0,-52.0);
  double exshift = 0.0;
  double p=0,q=0,r=0,s=0,z=0,t,w,x,y;

  // Store roots isolated by balanc and compute matrix norm
  double norm = 0.0;
  for (int i = 0; i < nn; i++) {
     if ((i < low) || (i > high)) {
        Eigendouble[i] = Hess[i][i];
        EigenImg[i] = 0.0;
     }
     for (int j = std::max(i-1,0); j < nn; j++) {
        norm += abs(Hess[i][j]);
     }
  }

  // Outer loop over eigenvalue index
  int iter = 0;
  while (n >= low) {

     // Look for single small sub-diagonal element
     int l = n;
     while (l > low) {
        s = abs(Hess[l-1][l-1]) + abs(Hess[l][l]);
        if (s == 0.0) {
           s = norm;
        }
        if (abs(Hess[l][l-1]) < eps * s) {
           break;
        }
        l--;
     }
   
     // Check for convergence
     // One root found
     if (l == n) {
        Hess[n][n] = Hess[n][n] + exshift;
        Eigendouble[n] = Hess[n][n];
        EigenImg[n] = 0.0;
        n--;
        iter = 0;

     // Two roots found
     } else if (l == n-1) {
        w = Hess[n][n-1] * Hess[n-1][n];
        p = (Hess[n-1][n-1] - Hess[n][n]) / 2.0;
        q = p * p + w;
        z = sqrt(abs(q));
        Hess[n][n] = Hess[n][n] + exshift;
        Hess[n-1][n-1] = Hess[n-1][n-1] + exshift;
        x = Hess[n][n];

        // double pair
        if (q >= 0) {
           if (p >= 0) {
              z = p + z;
           } else {
              z = p - z;
           }
           Eigendouble[n-1] = x + z;
           Eigendouble[n] = Eigendouble[n-1];
           if (z != 0.0) {
              Eigendouble[n] = x - w / z;
           }
           EigenImg[n-1] = 0.0;
           EigenImg[n] = 0.0;
           x = Hess[n][n-1];
           s = abs(x) + abs(z);
           p = x / s;
           q = z / s;
           r = sqrt(p * p+q * q);
           p = p / r;
           q = q / r;

           // Row modification
           for (int j = n-1; j < nn; j++) {
              z = Hess[n-1][j];
              Hess[n-1][j] = q * z + p * Hess[n][j];
              Hess[n][j] = q * Hess[n][j] - p * z;
           }

           // Column modification
           for (int i = 0; i <= n; i++) {
              z = Hess[i][n-1];
              Hess[i][n-1] = q * z + p * Hess[i][n];
              Hess[i][n] = q * Hess[i][n] - p * z;
           }

           // Accumulate transformations
           for (int i = low; i <= high; i++) {
              z = EigenMatrix[i][n-1];
              EigenMatrix[i][n-1] = q * z + p * EigenMatrix[i][n];
              EigenMatrix[i][n] = q * EigenMatrix[i][n] - p * z;
           }

        // Complex pair
        } else {
           Eigendouble[n-1] = x + p;
           Eigendouble[n] = x + p;
           EigenImg[n-1] = z;
           EigenImg[n] = -z;
        }
        n = n - 2;
        iter = 0;

     // No convergence yet
     } else {

        // Form shift
        x = Hess[n][n];
        y = 0.0;
        w = 0.0;
        if (l < n) {
           y = Hess[n-1][n-1];
           w = Hess[n][n-1] * Hess[n-1][n];
        }

        // Wilkinson's original ad hoc shift
        if (iter == 10) {
           exshift += x;
           for (int i = low; i <= n; i++) {
              Hess[i][i] -= x;
           }
           s = abs(Hess[n][n-1]) + abs(Hess[n-1][n-2]);
           x = y = 0.75 * s;
           w = -0.4375 * s * s;
        }

        // MATLAB's new ad hoc shift
        if (iter == 30) {
            s = (y - x) / 2.0;
            s = s * s + w;
            if (s > 0) {
				s = sqrt(s);
                if (y < x) {
                   s = -s;
                }
                s = x - w / ((y - x) / 2.0 + s);
                for (int i = low; i <= n; i++) {
                   Hess[i][i] -= s;
                }
                exshift += s;
                x = y = w = 0.964;
            }
        }
        iter = iter + 1;   // (Could check iteration count here.)

        // Look for two consecutive small sub-diagonal elements
        int m = n-2;
        while (m >= l) {
           z = Hess[m][m];
           r = x - z;
           s = y - z;
           p = (r * s - w) / Hess[m+1][m] + Hess[m][m+1];
           q = Hess[m+1][m+1] - z - r - s;
           r = Hess[m+2][m+1];
           s = abs(p) + abs(q) + abs(r);
           p = p / s;
           q = q / s;
           r = r / s;
           if (m == l) {
              break;
           }
           if (abs(Hess[m][m-1]) * (abs(q) + abs(r)) <
              eps * (abs(p) * (abs(Hess[m-1][m-1]) + abs(z) +
              abs(Hess[m+1][m+1])))) {
                 break;
           }
           m--;
        }

        for (int i = m+2; i <= n; i++) {
           Hess[i][i-2] = 0.0;
           if (i > m+2) {
              Hess[i][i-3] = 0.0;
           }
        }

        // Double QR step involving rows l:n and columns m:n
        for (int k = m; k <= n-1; k++) {
           int notlast = (k != n-1);
           if (k != m) {
              p = Hess[k][k-1];
              q = Hess[k+1][k-1];
              r = (notlast ? Hess[k+2][k-1] : 0.0);
              x = abs(p) + abs(q) + abs(r);
              if (x != 0.0) {
                 p = p / x;
                 q = q / x;
                 r = r / x;
              }
           }
           if (x == 0.0) {
              break;
           }
           s = sqrt(p * p + q * q + r * r);
           if (p < 0) {
              s = -s;
           }
           if (s != 0) {
              if (k != m) {
                 Hess[k][k-1] = -s * x;
              } else if (l != m) {
                 Hess[k][k-1] = -Hess[k][k-1];
              }
              p = p + s;
              x = p / s;
              y = q / s;
              z = r / s;
              q = q / p;
              r = r / p;

              // Row modification
              for (int j = k; j < nn; j++) {
                 p = Hess[k][j] + q * Hess[k+1][j];
                 if (notlast) {
                    p = p + r * Hess[k+2][j];
                    Hess[k+2][j] = Hess[k+2][j] - p * z;
                 }
                 Hess[k][j] = Hess[k][j] - p * x;
                 Hess[k+1][j] = Hess[k+1][j] - p * y;
              }

              // Column modification
              for (int i = 0; i <= std::min(n,k+3); i++) {
                 p = x * Hess[i][k] + y * Hess[i][k+1];
                 if (notlast) {
                    p = p + z * Hess[i][k+2];
                    Hess[i][k+2] = Hess[i][k+2] - p * r;
                 }
                 Hess[i][k] = Hess[i][k] - p;
                 Hess[i][k+1] = Hess[i][k+1] - p * q;
              }

              // Accumulate transformations
              for (int i = low; i <= high; i++) {
                 p = x * EigenMatrix[i][k] + y * EigenMatrix[i][k+1];
                 if (notlast) {
                    p = p + z * EigenMatrix[i][k+2];
                    EigenMatrix[i][k+2] = EigenMatrix[i][k+2] - p * r;
                 }
                 EigenMatrix[i][k] = EigenMatrix[i][k] - p;
                 EigenMatrix[i][k+1] = EigenMatrix[i][k+1] - p * q;
              }
           }  // (s != 0)
        }  // k loop
     }  // check convergence
  }  // while (n >= low)
  
  // Backsubstitute to find vectors of upper triangular form
  if (norm == 0.0) {
     return;
  }

  for (n = nn-1; n >= 0; n--) {
     p = Eigendouble[n];
     q = EigenImg[n];

     // double vector
     if (q == 0) {
        int l = n;
        Hess[n][n] = 1.0;
        for (int i = n-1; i >= 0; i--) {
           w = Hess[i][i] - p;
           r = 0.0;
           for (int j = l; j <= n; j++) {
              r = r + Hess[i][j] * Hess[j][n];
           }
           if (EigenImg[i] < 0.0) {
              z = w;
              s = r;
           } else {
              l = i;
              if (EigenImg[i] == 0.0) {
                 if (w != 0.0) {
                    Hess[i][n] = -r / w;
                 } else {
                    Hess[i][n] = -r / (eps * norm);
                 }

              // Solve double equations
              } else {
                 x = Hess[i][i+1];
                 y = Hess[i+1][i];
                 q = (Eigendouble[i] - p) * (Eigendouble[i] - p) + EigenImg[i] * EigenImg[i];
                 t = (x * s - z * r) / q;
                 Hess[i][n] = t;
                 if (abs(x) > abs(z)) {
                    Hess[i+1][n] = (-r - w * t) / x;
                 } else {
                    Hess[i+1][n] = (-s - y * t) / z;
                 }
              }

              // Overflow control
              t = abs(Hess[i][n]);
              if ((eps * t) * t > 1) {
                 for (int j = i; j <= n; j++) {
                    Hess[j][n] = Hess[j][n] / t;
                 }
              }
           }
        }

     // Complex vector
     } else if (q < 0) {
        int l = n-1;

        // Last vector component imaginary so matrix is triangular
        if (abs(Hess[n][n-1]) > abs(Hess[n-1][n])) {
           Hess[n-1][n-1] = q / Hess[n][n-1];
           Hess[n-1][n] = -(Hess[n][n] - p) / Hess[n][n-1];
        } else {
           cdiv(0.0,-Hess[n-1][n],Hess[n-1][n-1]-p,q);
           Hess[n-1][n-1] = cdivr;
           Hess[n-1][n] = cdivi;
        }
        Hess[n][n-1] = 0.0;
        Hess[n][n] = 1.0;
        for (int i = n-2; i >= 0; i--) {
           double ra,sa,vr,vi;
           ra = 0.0;
           sa = 0.0;
           for (int j = l; j <= n; j++) {
              ra = ra + Hess[i][j] * Hess[j][n-1];
              sa = sa + Hess[i][j] * Hess[j][n];
           }
           w = Hess[i][i] - p;

           if (EigenImg[i] < 0.0) {
              z = w;
              r = ra;
              s = sa;
           } else {
              l = i;
              if (EigenImg[i] == 0) {
                 cdiv(-ra,-sa,w,q);
                 Hess[i][n-1] = cdivr;
                 Hess[i][n] = cdivi;
              } else {

                 // Solve complex equations
                 x = Hess[i][i+1];
                 y = Hess[i+1][i];
                 vr = (Eigendouble[i] - p) * (Eigendouble[i] - p) + EigenImg[i] * EigenImg[i] - q * q;
                 vi = (Eigendouble[i] - p) * 2.0 * q;
                 if ((vr == 0.0) && (vi == 0.0)) {
                    vr = eps * norm * (abs(w) + abs(q) +
                    abs(x) + abs(y) + abs(z));
                 }
                 cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
                 Hess[i][n-1] = cdivr;
                 Hess[i][n] = cdivi;
                 if (abs(x) > (abs(z) + abs(q))) {
                    Hess[i+1][n-1] = (-ra - w * Hess[i][n-1] + q * Hess[i][n]) / x;
                    Hess[i+1][n] = (-sa - w * Hess[i][n] - q * Hess[i][n-1]) / x;
                 } else {
                    cdiv(-r-y*Hess[i][n-1],-s-y*Hess[i][n],z,q);
                    Hess[i+1][n-1] = cdivr;
                    Hess[i+1][n] = cdivi;
                 }
              }

              // Overflow control
              t = std::max(abs(Hess[i][n-1]),abs(Hess[i][n]));
              if ((eps * t) * t > 1) {
                 for (int j = i; j <= n; j++) {
                    Hess[j][n-1] = Hess[j][n-1] / t;
                    Hess[j][n] = Hess[j][n] / t;
                 }
              }
           }
        }
     }
  }

  // Vectors of isolated roots
  for (int i = 0; i < nn; i++) {
     if (i < low || i > high) {
        for (int j = i; j < nn; j++) {
           EigenMatrix[i][j] = Hess[i][j];
        }
     }
  }

  // Back transformation to get eigenvectors of original matrix
  for (int j = nn-1; j >= low; j--) {
     for (int i = low; i <= high; i++) {
        z = 0.0;
        for (int k = low; k <= std::min(j,high); k++) {
           z = z + EigenMatrix[i][k] * Hess[k][j];
        }
        EigenMatrix[i][j] = z;
     }
  }
}
///****************************************************   MCInterpolation2D    ****************************************************/
MCInterpolation2D :: MCInterpolation2D()
  : PointNum(0),
  m_Min (1.E+20),
  m_Max (-1.E+20),
  m_A(0)
{ }
MCInterpolation2D :: MCInterpolation2D(const int Num, const Point3d *ptPoints){
    Initialize(Num,ptPoints);
}

void MCInterpolation2D :: Initialize (const int Num, const Point3d *ptPoints){
  PointNum = Num; 
//  Assert(PointNum);
  ptData = new Point3d [PointNum];
  ptCoeff.Initialize(PointNum);
  m_Min=ptData[0].z; m_Max=ptData[0].z; m_A=0;

  for (int i=0; i<PointNum; i++) {
  //  Assert(ptPoints[i].z!=DUMMY_VALUE);
    ptData[i].x = ptPoints[i].x;
    ptData[i].y = ptPoints[i].y;
    ptData[i].z = ptPoints[i].z;
    if (ptData[i].z> m_Max ) m_Max = ptData[i].z;
    if (ptData[i].z> m_Max ) m_Min = ptData[i].z;
    m_A += ptData[i].z;
  }
  m_A /= PointNum;
  Norm = CalcNorma();
}

double MCInterpolation2D :: CalcNorma (){
//  Assert(PointNum);
  double minx = ptData[0].x;
  double maxx = ptData[0].x;
  double miny = ptData[0].y;
  double maxy = ptData[0].y;

  for (int i=1; i<PointNum; i++) {
    if (ptData[i].x < minx) minx = ptData[i].x;
    if (ptData[i].x > maxx) maxx = ptData[i].x;
    if (ptData[i].y < miny) miny = ptData[i].y;
    if (ptData[i].y > maxy) maxy = ptData[i].y;
  }
  double d=hypot(maxx-minx, maxy-miny);
  if(!d) d=1.;
  return d;
}

///****************************************************   MCInterpolation_RBF    ****************************************************/
MCInterpolation_RBF :: MCInterpolation_RBF()
  :MCInterpolation2D(),
  ptRBF(NULL)
{ }
MCInterpolation_RBF :: MCInterpolation_RBF(const int Num, const Point3d *ptPoints,
                       MCRadialBasis *RBF)
  : MCInterpolation2D(Num,ptPoints),
  ptRBF(RBF)
{ }

int MCInterpolation_RBF :: Interpolate (){
   //Assert(ptRBF==NULL);
   if(ptRBF==NULL)return 1;
   int ni = PointNum;
   int nj = PointNum;
   MCMatrix AA(ni,ni);
 
  for (int i=0; i< ni; i++) {
      for (int j=0; j< nj; j++) {
           double r =   hypot (ptData[i].x -ptData[j].x,ptData[i].y -ptData[j].y)/Norm; 
             AA[i][j]=  ptRBF->RBFuction(r); 
           if(i!=j)AA[j][i]=AA[i][j];
      }
   }

   MCVector B(nj);
   for (int j=0; j< nj; j++) {
      B[j] = ptData[j].z - m_A;
   }

  MC_LU LUlinearEq(AA,B);

  //----------  debug  -------------
  //AA.Print("d://Matrix.dat");
  //---------------------------------

  LUlinearEq.LU_Solve();
  ptCoeff = LUlinearEq.ptResult;

  //----------  debug  -------------
  //ptCoeff.Print("d://ptCoeff.dat");
  //---------------------------------
   return 0;
}

double MCInterpolation_RBF :: CalcValue (double x, double y)
{
  double result = 0;
  double r,tmp;

  for (int i=0; i< PointNum; i++) {
    r = (hypot(ptData[i].x-x, ptData[i].y-y))/Norm;
    tmp = ptRBF->RBFuction(r);
    result += ptCoeff[i]*tmp;
  }
  result += m_A;
  return result;
}

MCVector MCInterpolation_RBF :: CalcSurface(int nx, int ny, double minx, double miny,
                      double dx, double dy, double ang){
    int index = 0;
    MCVector Temp(ny*nx);

    for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++, index++) {
            double y =  miny + i*dx*sin(ang) + j*dy*cos(ang);
            double x =  minx + i*dx*cos(ang) - j*dy*sin(ang);
            Temp[index] = CalcValue(x,y);
        }
    }
    return Temp;
}
void MCInterpolation_RBF :: Initialization(const int Num, const Point3d *ptPoints,
                       MCRadialBasis *RBF){
    Initialize(Num,ptPoints);
    ptRBF = RBF;
}

///****************************************************   MCRadialBasis    ****************************************************/
MCRadialBasis :: MCRadialBasis(const MCRadialBasis &ptRBF)
{
    TypeRBF = ptRBF.TypeRBF;
    m_smooth = ptRBF.m_smooth;
}
double MCRadialBasis :: RBFuction(double Dist)
{
  switch (TypeRBF) {
    case DISTANCE_RBF:
      return DistanceRBF(Dist);
    case GAUSS_RBF: 
      return GaussRBF(Dist);
    case MULTIQUADRIC_RBF: 
      return MultiquadricRBF(Dist);
    case INVERSE_MULTIQUADRIC_RBF:
      return InverseMultiquadricRBF(Dist);
    case INVERT_QUADRIC_RBF:
      return InverseQuadricRBF(Dist);
    case THIN_PLATE_SPLINE_RBF:
      return TPS_RBF(Dist);
    //case MATERN_LINER_RBF:
     // return MaternLinerRBF(Dist);
    //case MATERN_SQUARE_RBF:
     // return MaternSquareRBF(Dist);
    //case BIGARMONIC_LAPLACE_RBF:
     // return BiharmonicLaplaceRBF(Dist);
    //case MACDONALD_RBF:
     // return MacdonaldRBF(Dist);
    //case BESSEL_RBF:
     // return BesselRBF(Dist,2);  // Attention! 'n=2'
      default:
     //     Assert(0);
      return InverseMultiquadricRBF(Dist);
  }
}
double MCRadialBasis :: DistanceRBF (double r) 
{
  double d=r/m_smooth;
  double fun = pow (1 + d*d,1./10.);
  return fun;
}

double MCRadialBasis :: GaussRBF (double r) 
{
  double d=r/m_smooth;
  return exp(-d*d);
}

double MCRadialBasis :: MultiquadricRBF (double r) 
{
  double d=r/m_smooth;
  double fun = sqrt (1 + d*d);
  return fun;
}

double MCRadialBasis :: InverseMultiquadricRBF (double r) 
{
  double d=r/m_smooth;
  double fun = 1/sqrt (1 + d*d);
  return fun;
}

double MCRadialBasis :: InverseQuadricRBF (double r) 
{
  double d=r/m_smooth;
  double fun = 1/(1 + d*d);
  return fun;
}

double MCRadialBasis :: TPS_RBF (double r) 
{
    if(r==0.)r=0.000001;
  double d=r;    //double d=r*m_smooth;
  return r*r*log(r);

  //if(r<=1.0){
   // return pow(d,n)*log(pow(d,d));
  //} else {
   // return pow(d,n)*log(d);
  //}
}

//double MCRadialBasis :: MaternLinerRBF (double r) 
//{
//  double d=r*m_smooth;
//  double fun = exp(-d)*(1 + d);
//  return fun;
//}
//
//double MCRadialBasis :: MaternSquareRBF (double r) 
//{
//  double d=r*m_smooth;
//  double fun = exp(-d)*(3 + 3*d + d*d);
//  return fun;
//}
//
//double MCRadialBasis :: MacdonaldRBF (double r) 
//{
//    if(r==0.)return 0.;      //r=0.01 ;  return 0.01;
//  double d=r*m_smooth;
//  double fun=d*besselk1(d);
//  return fun;
//}
//
//double MCRadialBasis :: BiharmonicLaplaceRBF (double r) 
//{
//    if(r==0.)return 0.;
//  double d=r/m_smooth;
//  //double C=0.577215665;
//  //double tau=0.1;    //0<tau<1.0
//  //double p=sqrt(tau/(1.0-tau));
//  double ui,op;
//  ui=_K0(d);
//  op=log(d);
//  double fun = ui + op;
//  //double fun = _K0(d)+log(d)-log(2.)+C;
//  return fun;
//}
//
//double MCRadialBasis :: BesselRBF (double r, int n) 
//{
//    if(r==0.)return 0.;      //r=0.01 ;  return 0.01;
//  double d=r*m_smooth;
// // Assert(n%2!=0);
//  double fun=besseljn(n/2-1,d)/pow(d,n/2-1);
//  return fun;
//}

// //***************************************CorrectTimeStep***************************************
//void CorrectTimeStep::CorrectStep(double &Step)
//{
//   Step = Step/2;
//}
//
//CorrectTimeStepGeomProgression::CorrectTimeStepGeomProgression()
//  : CorrectTimeStep(),
//  N(10),
//  denominator(2),
//  first_term(1),
//  current_step(1),
//  SplitStep(true)
//{
//}
//
//void CorrectTimeStepGeomProgression::CorrectStep(double &Step)
//{
//  if (SplitStep){
//    Step = first_term*pow(denominator,current_step-1);
//    current_step++;
//    if (current_step > N) SplitStep = false;
//  }
//}
//
//void CorrectTimeStepGeomProgression :: SetData (int _N, double q, double Step)
//{
//  N = _N;
//  denominator = q;
//  first_term = Step*(1.0-denominator)/(1.0-pow(denominator,N));
//}
//
//void CorrectTimeStepGeomProgression :: ReSetSplit(double Step)
//{
//  first_term = Step*(1.0-denominator)/(1.0-pow(denominator,N));
//  SplitStep = true;
//  current_step = 0;
//
//}
//
//CorrectTimeStepforNewton::CorrectTimeStepforNewton() 
//  : CorrectTimeStepGeomProgression(),
//  Newton_fast(0),
//  Newton_slow(0),
//  IterForStepDecrease(6),
//  IterForStepIncrease(3),
//  LimitIter(3),
//  TimeStepMax(1), 
//  TimeStepMin(1)
//{
//}
//
//void CorrectTimeStepforNewton::SetData(double Max, double Min, int N_first, double q, double Step)
//{
//  TimeStepMin = Min;
//  TimeStepMax = Max;
//  N = N_first;
//  denominator = q;
//  first_term = Step*(1.0-denominator)/(1.0-pow(denominator,N));
//}
//
//void CorrectTimeStepforNewton:: CorrectStepforNewton(double &Step, int Newton_iter)
//{
//  if (!SplitStep){ 
//    if (Newton_iter <= IterForStepIncrease){
//      Newton_fast++;
//      Newton_slow = 0;
//    }
//    else if (Newton_iter >= IterForStepDecrease){
//      Newton_slow++;
//      Newton_fast = 0;
//    }
//    else{
//      Newton_fast = 0;
//      Newton_slow = 0;
//    }
//  
//    if (Newton_fast >= LimitIter) {
//      Step = Step*2;
//      Newton_fast =0;
//      if ( Step > TimeStepMax)  Step = TimeStepMax;
//    }
//    else if(Newton_slow >=  LimitIter || Newton_iter > 10){
//      Step = Step/2;
//      if (Step < TimeStepMin)   Step = TimeStepMin;
//      ReSetSplit(Step);
//      Newton_slow =0;
//    }    
//  }
//}
//
//double CorrectTimeStepforNewton::GetMinStep()
//{
//  return TimeStepMin;
//}
//
//bool CorrectTimeStepforNewton::DoesSplitStep()
//{
//  return SplitStep;
//}