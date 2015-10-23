#ifndef TWOPHASEMODEL_H_
#define TWOPHASEMODEL_H_

#include <vector>
#include <new>
#include <string>
#include <cassert>
#include <algorithm>
#include <utility>
#include <time.h>

#include "sweep.h"
#include "Interpolate.h"
#include "CylindricalCell.h"
#include "Snapshotter.h"
#include "utils.h"

#define BAR_TO_PA 1.E5
#define EQUALITY_TOLERANCE 1.E-6

#define TEMP 0
#define PRES 1
#define SAT 2

struct Properties
{
	// Vector of start times of periods [sec]
	std::vector<double> timePeriods;
	// Vector of rates [m3/day]
	std::vector<double> rates;
	// Vector of skins
	std::vector<double> skins;
	// Radius of damaged zone [m]
	std::vector<double> radius;
	
	// Time step limits
	// Initial time step [sec]
	double ht;
	// Minimal time step [sec]
	double ht_min;
	// Maximum time step [sec]
	double ht_max;
	// During the time flow rate decreases 'e' times in well test [sec] 
	double alpha;

	// Inner radius of well [m]
	double r_w;
	// Radius of formation [m]
	double r_e;
	
	// Number of cells in radial direction
	double cellsNum_r;
	// Number of cells in vertical direction
	double cellsNum_z;
	
	// Porosity
	double m; 
	// Height of formation [m]
	double height;
	double h1, h2;
	// BHP will be converted to the depth
	double depth_point;

	// Compessibilities [1/Pa]
	double beta_oil;
	double beta_sk;

	// Permeability along radial direction [mD]
	double perm_r;
	// Permeability along vertical direction [mD]
	double perm_z;

	// Data set (saturation, relative oil permeability)
	std::vector< std::pair<double,double> > kr_oil;
	// Data set (saturation, relative gas permeability)
	std::vector< std::pair<double,double> > kr_gas;
	
	// Viscosity of oil [cP]
	double visc_oil;
	// Viscosity of gas [cP]
	double visc_gas;

	// Data set (pressure, oil volume factor) ([Pa], [m3/m3])
	std::vector< std::pair<double,double> > B_oil;
	// Data set (pressure, gas volume factor) ([Pa], [m3/m3])
	std::vector< std::pair<double,double> > B_gas;

	// Data set (pressure, gas content in oil) ([Pa], [m3/m3])
	std::vector< std::pair<double,double> > Rs;
	
	// Densities in STC [kg/m3]
	double dens_oil_stc;
	double dens_gas_stc;
	double dens_sk_stc;

	// Initial formation pressure, right boundary pressure [Pa]
	double p_init;
	// Pressure of fully saturated oil [Pa]
	double p_sat;

	// Initial oil saturation, right boundary of saturation
	double s_init;
	// Initial temperature, right boundary temperature [K]
	double T_init;

	// Thermal properties
	// Joule-thompson coefficients [K/Pa]
	double jt_oil;
	double jt_gas;

	// Adiabatic coefficients [K/Pa]
	double ad_oil;
	double ad_gas;

	// Mass heat capacities [J/kg/K]
	double c_oil;
	double c_gas;
	double c_sk;

	// Heat of phase transition [J/kg]
	double L;

	// Thermal diffusivity coefficient [m2/sec]
	double kappa_eff;

	// Thermal conductivity coefficients [W/m/K]
	double lambda_sk;
	double lambda_oil;
	double lambda_gas;
};

struct Skeleton_Props
{
	// Porosity in STC
	double m; 
	// Permeability along radial direction [mD]
	double perm_r;
	// Permeability along vertical direction [mD]
	double perm_z;
	// Permeability of colmatage zone [mD]
	std::vector<double> perm_eff;
	// Radius of colmatage zone [m]
	std::vector<double> r_eff;
	// Vector of skins
	std::vector<double> skin;
	// Top and bottom depth of perforation
	double h1, h2;

	// Height of formation [m]
	double height;
	// Density of skeleton matter in STC [kg/m3]
	double dens_stc;
	// Compessibility [1/Pa]
	double beta;

	// Thermal properties

	// Mass heat capacity [J/kg/K]
	double c;
	// Thermal conductivity coefficient [W/m/K]
	double lambda;
};

struct Fluid_Props
{
	// Viscosity [cP]
	double visc;
	// Density of fluid in STC [kg/m3]
	double dens_stc;
	// Compessibility [1/Pa]
	double beta;
	// Relative fluid permeability
	Interpolate* kr;
	// Fluid volume factor
	Interpolate* b;

	// Thermal properties

	// Mass heat capacity [J/kg/K]
	double c;
	// Thermal conductivity coefficient [W/m/K]
	double lambda;
	// Joule-thompson coefficient [K/Pa]
	double jt;
	// Adiabatic coefficient [K/Pa]
	double ad;
};

class TwoPhaseModel : public Sweep {
	friend class Snapshotter;
	friend class TimeManager;
private:
	// Fill matrixes belong to inner cell
	void MiddleAppr(int current, int MZ, int key);
	// Fill matrix at left cell
	void LeftBoundAppr(int MZ, int key);
	// Fill matrix at right cell
	void RightBoundAppr (int MZ, int key);

protected:
	// Continuum properties
	Skeleton_Props props_sk;
	Fluid_Props props_oil;
	Fluid_Props props_gas;

	// Rate of the well
	double Q_sum;
	// Vector of <cell number, rate in the cell> for left border cells
	std::vector< std::pair<int,double> > Qcell;
	// The summary height of perforated interval
	double height_perf;

	double r_eff;
	double Perm_eff;
	double skin;

	// Temporal properties

	// Initial time step [sec]
	double ht;
	// Minimal time step [sec]
	double ht_min;
	// Maximum time step [sec]
	double ht_max;
	
	// Initial values of variables
	double T_init;
	double P_init;
	double S_init;
	double P_bub;

	// Number of periods
	int periodsNum;
	// End times of periods [sec]
	std::vector<double> period;
	// Oil rates [m3/day]
	std::vector<double> rate;

	// Grid properties
	
	// Inner radius of well [m]
	double r_w;
	// Radius of formation [m]
	double r_e;
	// Number of cells in radial direction
	double cellsNum_r;
	// Number of cells in vertical direction
	double cellsNum_z;

	// Other properties

	// Gas content in oil
	Interpolate* Rs;
	Interpolate* Prs;

	// Total reservoir volume
	double Volume;
	// Thermal diffusivity coefficient [m2/sec]
	double kappa_eff;
	// Heat of phase transition [J/kg]
	double L;
	// BHP will be converted to the depth
	double depth_point;
	// During the time flow rate decreases 'e' times in well test [sec] 
	double alpha;
	// Dimensions
	double t_dim;
	double R_dim;
	double P_dim;
	double T_dim;
	double Q_dim;

	// Grid
	std::vector<CylCell> cells;
	// Snapshotter
	Snapshotter* snapshotter;

	// Set all properties
	void setProps(Properties& props);
	// Set perforated cells 
	void setPerforated();
	// Set initial conditions
	void setInitialState();
	// Make all properties dimensionless
	void makeDimLess();
	// Build grid
	void buildGridLog();

	// Set data sets
	Interpolate* setDataset(std::vector< std::pair<double,double> >& vec, const double xDim, const double yDim);
	Interpolate* setInvDataset(std::vector< std::pair<double,double> >& vec, const double xDim, const double yDim);

	// Calculate boundary conditions each time step
	void calculateBounds();
	// Calculate perforated ghost cells
	void calculatePerforated();

	// Service functions
	inline double upwindIsCur(int cur, int beta)
	{
		return 0.0 ? cells[cur].u_next.pres < cells[beta].u_next.pres : 1.0;
	};
	inline int getUpwindIdx(int cur, int beta)
	{
		return beta ? cells[cur].u_next.pres < cells[beta].u_next.pres : cur;
	};
	inline void getNeighborIdx(int cur, int* const neighbor)
	{
		neighbor[0] = cur - cellsNum_z - 2; 
		neighbor[1] = cur + cellsNum_z + 2;
		neighbor[2] = cur - 1;
		neighbor[3] = cur + 1;
	};

	// Solving coefficients
	inline double getPoro(double p)
	{
		return props_sk.m * (1.0 + props_sk.beta * (P_init - p) );
	};
	inline double getPoro_dp(double p)
	{
		return -props_sk.m * props_sk.beta;
	};
	inline double getTrans(int idx1, int idx2)
	{
		CylCell& cell1 = cells[idx1];
		CylCell& cell2 = cells[idx2];

		if( abs(idx1 - idx2) == 1)
			return 2.0 * props_sk.perm_z * 2.0 * M_PI * cell1.r * cell1.hr / (cell1.hz + cell2.hz);
		else {
			double k1, k2, S;	
			k1 = (cells[idx1].r > r_eff ? props_sk.perm_r : Perm_eff);
			k2 = (cells[idx2].r > r_eff ? props_sk.perm_r : Perm_eff);
			S = 2.0 * M_PI * cell1.hz * (cell1.r + sign(idx2 - idx1) * cell1.hr / 2.0);
			return 2.0 * k1 * k2 * S / (k1 * cell2.hr + k2 * cell1.hr);
		}
	};
	inline double getKr_oil(double sat_oil)
	{
		if(sat_oil > 1.0)
			return 1.0;
		else
			return props_oil.kr->Solve(sat_oil);
	};
	inline double getKr_oil_ds(double sat_oil)
	{
		if(sat_oil > 1.0)
			return props_oil.kr->DSolve(1.0);
		else
			return props_oil.kr->DSolve(sat_oil);
	};
	inline double getKr_gas(double sat_oil)
	{
		if(sat_oil > 1.0)
			return 0.0;
		else
			return props_gas.kr->Solve(sat_oil);
	};
	inline double getKr_gas_ds(double sat_oil)
	{
		if(sat_oil > 1.0)
			return props_gas.kr->DSolve(1.0);
		else
			return props_gas.kr->DSolve(sat_oil);
	};
	inline double getB_oil(double p, double p_bub, bool SATUR)
	{
		if(SATUR)
			return props_oil.b->Solve(p);
		else
			return props_oil.b->Solve(p_bub) * (1.0 + props_oil.beta * (p_bub - p));
	};
	inline double getB_oil_dp(double p, double p_bub, bool SATUR)
	{
		if(SATUR)
			return props_oil.b->DSolve(p);
		else
			return -props_oil.b->Solve(p_bub) * props_oil.beta;
	};
	inline double getB_gas(double p)
	{
		return props_gas.b->Solve(p);
	};
	inline double getB_gas_dp(double p)
	{
		return props_gas.b->DSolve(p);
	};
	inline double getRs(double p, double p_bub, bool SATUR)
	{
		if(SATUR)
			return Rs->Solve(p);
		else
			return Rs->Solve(p_bub);
	};
	inline double getRs_dp(double p, double p_bub, bool SATUR)
	{
		if(SATUR)
			return Rs->DSolve(p);
		else
			return 0.0;
	};
	inline double getPresFromRs(double rs)
	{
		return Prs->Solve(rs);
	};
	inline void solveP_bub()
	{
		int idx;
		double factRs, dissGas, p_sat_new;

		for(int i = 0; i < cellsNum_r; i++)
			for(int j = 0; j < cellsNum_z; j++)
			{
				idx = (i + 1) * (cellsNum_z + 2) + j + 1;

				Variables& next = cells[idx].u_next;
				Variables& prev = cells[idx].u_prev;

				if(next.sat_oil > 1.0)
					next.sat_oil = 1.0;

				dissGas = (1.0 - next.sat_oil) * getB_oil(next.pres, next.pres_bub, next.SATUR) / ( (1.0 - next.sat_oil) * getB_oil(next.pres, next.pres_bub, next.SATUR) + next.sat_oil * getB_gas(next.pres));
				factRs = getRs(prev.pres, prev.pres_bub, next.SATUR) + dissGas;

				if(getRs(next.pres, next.pres, next.SATUR) > factRs)
				{
					next.pres_bub = getPresFromRs(factRs);
					next.SATUR = false;
				} else {
					next.pres_bub = next.pres;
					next.SATUR = true;
				}
			}
	};

	void construction_from_fz(int N, int n, int key);
	void copyIterLayer();
	double averValue(int var_idx);
	double convergance(int key, int& ind, int& var_idx);

	// First eqn
	double solve_eq1(int cur);
	double solve_eq1_dp(int cur);
	double solve_eq1_ds(int cur);
	double solve_eq1_dp_beta(int cur, int beta);
	double solve_eq1_ds_beta(int cur, int beta);

	// Second eqn
	double solve_eq2(int cur);
	double solve_eq2_dp(int cur);
	double solve_eq2_ds(int cur);
	double solve_eq2_dp_beta(int cur, int beta);
	double solve_eq2_ds_beta(int cur, int beta);

	// Horizontal border approximators
	void TopAppr(int i);
	void BottomAppr(int i);

public:
	TwoPhaseModel();
	~TwoPhaseModel();

	void load(Properties& props);
	// Set time period
	void setPeriod(int period);
	
	void snapshot(int i);
	void snapshot_all(int i);

	int doNextStep();
	void copyNewLayer();
};

#endif /* TWOPHASEMODEL_H_ */