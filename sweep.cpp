#include "sweep.h"

double ** auxAlloc2Darray (int ni, int nj)
{
  int _ni = ni+1;
  int _nj = nj+1;
  double ** ptArray = (double **) calloc (_ni,sizeof(double *));
  for (int i=0; i<_ni; i++) {
    ptArray[i]= (double *)calloc (_nj,sizeof(double));
  }
  return ptArray;
}

void auxFree2Darray (double ** ptArray, int ni, int /*nj*/)
{
  int _ni = ni +1;
  //int _nj = nj +1;

  for (int i=0;i<_ni;i++){
    free (ptArray[i]);
  }
  free (ptArray);
}


double *** auxAlloc3Darray (int ni,int nj)
{
  int i,j;
  int _ni = ni +1;
  int _nj = nj +1;
  double *** ptArray = (double ***) calloc(_ni,sizeof(double **));
  if (ptArray == NULL)  return NULL;

  for(i=0; i<_ni;i++) {
    ptArray[i] = (double **) calloc(_nj,sizeof(double*));
    if (ptArray[i] == NULL) return NULL;
    for(j=0;j<_nj;j++)  {
      ptArray[i][j] = (double *) calloc(_nj,sizeof(double));
      if (ptArray[i][j] == NULL) return NULL;
    }
  }
  return ptArray;
}

void auxFree3Darray (double *** ptArray, int ni, int nj)
{
  int _ni = ni +1;
  int _nj = nj +1;

  for (int i=0; i<_ni; i++)  {
    for(int j=0; j<_nj; j++)
      free (ptArray[i][j]);
    free (ptArray[i]);
  }
  free (ptArray);
}



//#define BZ_I_DIM 2
//#define BZ_J_DIM 4

Sweep :: Sweep ()
{
  bz = NULL;
  fz = NULL;
}

Sweep :: Sweep (int PointNum, int EqNum)
{

   BZ_I_DIM = 2*EqNum; 
   BZ_J_DIM = 1+ 3*EqNum;

   bz = auxAlloc2Darray (BZ_I_DIM, BZ_J_DIM);
   fz = auxAlloc2Darray (PointNum, EqNum);

   A.Initialize(EqNum,EqNum);
   B.Initialize(EqNum,EqNum);
   C.Initialize(EqNum,EqNum);
   RightSide.Initialize(EqNum,1);

 
   MZ = EqNum;
   NX = PointNum;
}

Sweep :: ~Sweep ()
{
  Clear ();
}

void Sweep :: Clear ()
{
  if (bz) auxFree2Darray (bz, BZ_I_DIM, BZ_J_DIM);
  if (fz) auxFree2Darray (fz, NX, MZ);
}


void Sweep :: Initialize (int PointNum, int EqNum)
{
   MZ = EqNum;
   NX = PointNum;

   BZ_I_DIM = 2*EqNum; 
   BZ_J_DIM = 1+ 3*EqNum;

   bz = auxAlloc2Darray (BZ_I_DIM, BZ_J_DIM);
   fz = auxAlloc2Darray (NX, MZ);

   A.Initialize(EqNum,EqNum);
   B.Initialize(EqNum,EqNum);
   C.Initialize(EqNum,EqNum);
   RightSide.Initialize(EqNum,1);


}
void Sweep :: LeftBoundAppr (int /*EqNum*/,int )
{
}
void Sweep :: RightBoundAppr (int /*EqNum*/, int )
{
}
void Sweep :: MiddleAppr (int /*m*/, int /*EqNum*/,int )
{
}


void Sweep :: Solve(int NZ,int MZ, int key)
{
  int i,j,k,l,m,n;
  double cz[3];

  double  ***pz;
  double  ***qz;
  pz = NULL;
  qz = NULL;
  pz = auxAlloc3Darray(NZ,MZ);
  qz = auxAlloc3Darray(NZ,MZ);

  for (i=1;i<= MZ;i++) {
    for (j=1;j<= 3*MZ+1;j++) {
      bz[i][j] = 0.0;
    }
  }
  RightBoundAppr(MZ,key);
  goto Pr102;
Pr102:
  for (n=NZ;n>=0;n--)
  {
   if (n==1)  goto Pr104;
   if (n==0)
    {
   for (i=1+MZ;i<= 2*MZ;i++)
   for (j=1;j<= 3*MZ+1;j++)
    bz[i][j] = 0.0;
    goto Pr101;
     }
   m = n-1;
   goto Pr105;

Pr104:    for (i=MZ+1;i<= 2*MZ;i++)
    for (j=1;j<= 3*MZ+1;j++)
    bz[i][j] = 0.0;
   LeftBoundAppr(MZ,key);
  for (k=1;k<= MZ;k++)
  for (j=MZ+1;j<= 2*MZ;j++)
  {
   if ( bz[j][k] != 0.0)
   {
     for (i=1;i<= MZ;i++)
     {
       bz[j][i+MZ] = bz[j][i+MZ]-pz[2][k][i]*bz[j][k];
       bz[j][i+2*MZ] = bz[j][i+2*MZ]-qz[2][k][i]*bz[j][k];
      }
       bz[j][3*MZ+1] = bz[j][3*MZ+1]-(fz[2][k])*bz[j][k];
    }
   }
   for (i=1;i<=MZ;i++)
   for (j=MZ+1;j<= 2*MZ;j++)
     {
     bz[j][i] = bz[j][i+MZ];
     bz[j][i+MZ] = bz[j][i+2*MZ];
     bz[j][i+2*MZ] = 0.0;
    }
     goto Pr101;
Pr105:   for (i=MZ+1;i<= 2*MZ;i++)
   for (j=1;j<= 3*MZ+1;j++)
     bz[i][j] = 0.0;
     MiddleAppr(m,MZ,key);
       goto Pr101;

Pr101:   for (l=1;l<= MZ;l++)
    {
       cz[1] = 0.0;
      // Maximum searching
	   for (i=l;i<= 2*MZ;i++)
		   if ( fabs(bz[i][l])>fabs(cz[1]))
     {
       k = i;
        cz[1] = bz[i][l];
      }
   // Maximum dividing & swapping
        for (j=l;j<= 3*MZ+1;j++)
    {
      bz[k][j] = bz[k][j]/cz[1];
        if (k != l)
          { cz[2] = bz[k][j];
      bz[k][j] = bz[l][j];
      bz[l][j] = cz[2];
           }
     }
       for (i=1;i<= 2*MZ;i++)
   { if (i != l)
     { cz[1] = bz[i][l];
        for (j=l+1;j<= 3*MZ+1;j++)
        bz[i][j] = bz[i][j]-bz[l][j]*cz[1];
      }
	}
	}
    for (i=1;i<= MZ;i++)
  { (fz[n][i]) = bz[i][3*MZ+1];
    if (n != 0)
    { for (j=1;j<= MZ;j++)
      { pz[n][i][j] = bz[i][j+MZ] ;
        if (n != 1)
      qz[n][i][j] = bz[i][j+2*MZ];
       }
     }
   }
    if (n==0)  goto Pr106;
    for (i=1;i<= MZ;i++)
    { bz[i][3*MZ+1] = bz[i+MZ][3*MZ+1];
      for (j=1;j<= 2*MZ;j++)  bz[i][j] = bz[i+MZ][j+MZ];
      for (j=2*MZ+1;j<= 3*MZ;j++)  bz[i][j] = 0.0;
     }
    }
Pr106:
     for (j=1;j<=MZ;j++)
     for (i=1;i<=MZ;i++)
     fz[1][j] = fz[1][j]-pz[1][j][i]*(fz[0][i]);
     for (k=2;k<=NZ;k++)
     for (j=1;j<=MZ;j++)
     for (i=1;i<=MZ;i++)
     fz[k][j] = fz[k][j]-pz[k][j][i]*(fz[k-1][i])-qz[k][j][i]*(fz[k-2][i]);

   
    
   auxFree3Darray (pz,NZ,MZ);
     auxFree3Darray (qz,NZ,MZ);
}

#define top_floor 1
#define ground_floor 2
void Sweep :: construction_bz(int num_eq, int flag)
{
    int i,j;
    switch(flag){
      case top_floor:
        for(i=1;i< num_eq+1;i++){
          bz[i][3*num_eq+1]=RightSide[i-1][0];
          for(j=1;j<num_eq+1;j++){
            bz[i][j] = A[i-1][j-1];
            bz[i][num_eq+j] = B[i-1][j-1];
            bz[i][2*num_eq+j] = C[i-1][j-1];
          }
        }
      break;
      case ground_floor:
        for(i=1;i< num_eq+1;i++){
          bz[num_eq+i][3*num_eq+1]=RightSide[i-1][0];
          for(j=1;j<num_eq+1;j++){
            bz[num_eq+i][j] = A[i-1][j-1];
            bz[num_eq+i][num_eq+j] = B[i-1][j-1];
            bz[num_eq+i][2*num_eq+j] = C[i-1][j-1];
          }
        }          
      break;
    }
}

void MCFlowEq :: construction_from_fz(int N, int m, int key)
{
  int i,j;
  for (i=0;i < N;i++)
  {
		if(key == 1)     pNodeArray[i].NewPress[1] = fz[i][j+1];
		else			 pNodeArray[i].NewPress[0] = fz[i][j+1];
  }


}

/*************************************************************************************************************/

FlowInputData  :: FlowInputData ()
                : cellNum (500)
                , eqNum (1)
                , R (500)
                , r (0.1)
                , timeStart (0)
                , timeStop  (24)               // one day
                , timeStep  (0.5)              // half an hour
        , rateNum   (2) 
{               
};

FlowInputData  :: ~FlowInputData () 
{
}
/*************************************************************************************************************/
MCFlowEq :: MCFlowEq ()
          : Sweep ()
{
  m_iLinear = 0;
}

MCFlowEq :: ~MCFlowEq ()
{
}

void MCFlowEq :: Setup (FlowInputData * ptStruct)
{
   m_iNodeNum  = ptStruct->cellNum +1;                       
   pNodeArray  = (GRID_NODE *)calloc (m_iNodeNum, sizeof(GRID_NODE));
   m_eqNum     = ptStruct->eqNum;
   m_R         = ptStruct->R;
   m_r         = ptStruct->r;

   m_Err       = 1.E-5;
   m_MinLevel  = 1.E-6;
   m_iIterMax  = 10;

 //  BuildGrid    ();
   InitialState ();
}

void MCFlowEq :: InitialState ()
{
    for (int i=0; i < m_iNodeNum; i++) {
     double r =  (pNodeArray + i)->xnew;
     double y =  1;

     for (int ieq =0; ieq < m_eqNum; ieq++) {
       (pNodeArray + i)->OldPress[ieq] = y;          // initial pressure distribution is stationary
       (pNodeArray + i)->NewPress[ieq] = y;
       (pNodeArray + i)->ItrPress[ieq] = y;
     }
  }
}


void MCFlowEq :: Presilka_Sloev ()
{
  m_tPrevStep = m_tStep;
  for (int i=0; i < m_iNodeNum; i++)  {
    GRID_NODE  * pCur = &pNodeArray[i];
    for (int k=0; k < m_eqNum; k++) {
      pCur->PrevPress[k] = pCur->OldPress[k];
      pCur->OldPress[k] = pCur->NewPress[k];
      pCur->ItrPress[k] = pCur->NewPress[k];
    }
  }
}

void MCFlowEq::Presilka_Sloev_Back(){
  m_tStep = m_tPrevStep;
  for (int i=0; i < m_iNodeNum; i++)  {
    GRID_NODE  * pCur = &pNodeArray[i];
    for (int k=0; k < m_eqNum; k++) {
      pCur->OldPress[k] = pCur->PrevPress[k];
      pCur->NewPress[k] = pCur->OldPress[k];
      pCur->ItrPress[k] = pCur->OldPress[k];
    }
  }
}
double  MCFlowEq :: Convergance ()
{
  double relErr = 0;
  for (int i=0; i < m_iNodeNum; i++)  {
    GRID_NODE  * pCur = &pNodeArray[i];
     for (int k=0; k < m_eqNum; k++) {
        double var_new  = pCur->NewPress [k];
        double var_iter = pCur->ItrPress[k];
        if (fabs(var_new) > m_MinLevel)  {
           double cur_relErr = fabs((var_new-var_iter)/var_new);
           if (cur_relErr > relErr) {
               relErr = cur_relErr;
		   }
        }
  }
      }
  return relErr;
}


double  MCFlowEq :: Convergance (int &bParam, int &bNode, double &bVal, int indBeg)
{
  double relErr = 0;
  for (int i=indBeg; i < m_iNodeNum; i++)  {
    GRID_NODE  * pCur = &pNodeArray[i];
     for (int k=0; k < m_eqNum; k++) {
        double var_new  = pCur->NewPress [k];
        double var_iter = pCur->ItrPress[k];
        if (fabs(var_new) > m_MinLevel)  {
           double cur_relErr = fabs((var_new-var_iter)/var_new);
           if (cur_relErr > relErr) {
               relErr = cur_relErr;
			   bNode = i;
			   bParam = k; 
			   bVal = var_new;
           }
        }
  }
      }
  return relErr;
}




void MCFlowEq :: Peresilka_Iter () 
{
  for (int i=0; i < m_iNodeNum; i++)  {
    GRID_NODE  * pCur = &pNodeArray[i];
    for (int k=0; k < m_eqNum; k++) {
      pCur->ItrPress[k] = pCur->NewPress[k];
    }
  }
}
double MCFlowEq::GetStartTime()
{
  return m_tStart;
}
double MCFlowEq::GetStepTime()
{
  return m_tStep;
}
double MCFlowEq::GetStopTime()
{
  return m_tStop;
}
int MCFlowEq::GetIterMax()
{
  return m_iIterMax;
}
double MCFlowEq::GetErr()
{
  return m_Err;
}
int MCFlowEq::GetEqNum()
{
  return m_eqNum;
}

int MCFlowEq::GetNumCell()
{
  return m_iNodeNum;
}

double MCFlowEq::GetCenterOfMass(int i)
{
  return pNodeArray[i].cm;
}
void MCFlowEq :: GetSolution ()
{
  for (int i=0; i < m_iNodeNum; i++){
    for (int k=1; k< m_eqNum+1; k++){
     (pNodeArray + i)->NewPress[k-1] = fz[i][k];
//     (pNodeArray + i)->ItrPress[k-1] = fz[i][k];
   }
  }
}



/****************************************************************************************************************************/

