#include "mcmath.h"
#include <math.h>
class Sweep
{
private: 
  int MZ;                                  // number of equations
  int NX;                                  // number of points in grid
  int BZ_I_DIM;
  int BZ_J_DIM; 

public:
  MCMatrix A;                 // Block matrices
  MCMatrix B;
  MCMatrix C;
  MCMatrix RightSide;
  double ** bz;                             // aux array
  double ** fz; 
  
  Sweep        ();                           // default constructor
  Sweep        (int PointNum, int EqNum);
  ~Sweep       ();

  void Initialize (int PointNum, int EqNum);
  virtual void Solve   (int, int, int );                               
  void construction_bz(int num_eq, int flag);
  void construction_from_fz(int N, int m,int key);

 
  // array for results
private:
  virtual void LeftBoundAppr  (int,int);         // left boundary approximation
  virtual void MiddleAppr     (int, int,int);   // approximation in the middle node
  virtual void RightBoundAppr (int,int);        // right boundary approximation

private:
  Sweep(const Sweep&);
  void operator=(const Sweep&);
  
protected:
  void Clear ();

};
//class MCMatrix;


#define EQ_MZN  7//предусматривает количество уравнений EQ_MZN+1
typedef struct  {

  // ¬ектор неизвестных , нелинейные итерации
  double NewPress[EQ_MZN+1];
  double ItrPress[EQ_MZN+1];
  double OldPress[EQ_MZN+1];
  double PrevPress[EQ_MZN+1];
  MCMatrix alpha;
  MCMatrix beta;
  double xnew;// ≈сли задача одномерна€, координа начала отрезка
  double h;  //длина €чейки
  double hz;// высота €чейки 
  double cm;// центр масс €чейки
  double V; // объем €чейки
} GRID_NODE;
typedef GRID_NODE  * PTRGRID_NODE;


struct  FlowInputData {
  FlowInputData ();
  ~FlowInputData ();
  int  cellNum;    
  int  eqNum; 
  double R;
  double r;

  double timeStart;
  double timeStop;
  double timeStep;

  int  rateNum;
};

// class for solution set of equations 
class MCFlowEq :public Sweep 
{

 protected:
   int   m_iNodeNum;
   int   m_eqNum;

   int   m_iLinear;
   int   m_iIterMax;
   double  m_Err;
   double  m_MinLevel;   

   double  m_tStart;
   double  m_tStop;
   double  m_tStep;
   double  m_tPrevStep;
   
  
           
   double  m_R;
   double  m_r;
   virtual void InitialState   ();
   void   PrintIterArrays    ( const char  * szFile);
   void   PrintDescrepancy   ( const char  * szFile);


public:
   MCPtr <GRID_NODE>  pNodeArray;   


//private:
 //  virtual void BuildGrid ();

public:
  MCFlowEq   ();
  ~MCFlowEq  ();
  virtual void Setup (FlowInputData * ptStruct);
  virtual void Presilka_Sloev ();
  virtual void Presilka_Sloev_Back();
  void construction_from_fz(int N,int m, int key);
  void GetSolution ();
  double  Convergance ();
  double  Convergance  (int& bParam,int& bNode, double& bVal,int indBeg);
  void Peresilka_Iter (); 
  double GetStartTime();
  double GetStopTime();
  double GetStepTime();
  int GetIterMax();
  double GetErr(); 
  int GetEqNum();
  int GetNumCell();
  virtual double GetCenterOfMass(int i);
};





