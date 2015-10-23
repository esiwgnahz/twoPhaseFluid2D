#include "TwoPhaseModel.h"

using namespace std;

TwoPhaseModel::TwoPhaseModel()
{
	snapshotter = new Snapshotter();
	snapshotter->setModel(this);
}

TwoPhaseModel::~TwoPhaseModel()
{
}

void TwoPhaseModel::load(Properties& props)
{
	setProps(props);

	buildGridLog();
	setPerforated();
	setInitialState();

	Initialize(cellsNum_r, 2*cellsNum_z);
}

void TwoPhaseModel::setProps(Properties& props)
{
	// Setting initial condition
	T_init = props.T_init;
	P_init = props.p_init;
	S_init = props.s_init;
	P_bub = props.p_sat;

	// Setting grid properties
	r_w = props.r_w;
	r_e = props.r_e;
	cellsNum_r = props.cellsNum_r;
	cellsNum_z = props.cellsNum_z;

	// Setting skeleton properties
	props_sk.m = props.m;
	props_sk.perm_r = MilliDarcyToM2( props.perm_r );
	props_sk.perm_z = MilliDarcyToM2( props.perm_z );
	props_sk.dens_stc = props.dens_sk_stc;
	props_sk.beta = props.beta_sk;
	props_sk.height = props.height;
	props_sk.h1 = props.h1;
	props_sk.h2 = props.h2;

	props_sk.c = props.c_sk;
	props_sk.lambda = props.lambda_sk;

	periodsNum = props.timePeriods.size();
	for(int i = 0; i < periodsNum; i++)
	{
		period.push_back( props.timePeriods[i] );
		rate.push_back( props.rates[i] / 86400.0 );
		props_sk.skin.push_back( props.skins[i] );
		props_sk.r_eff.push_back( props.radius[i] );
		if(props.radius[i] > props.r_w)
			props_sk.perm_eff.push_back( MilliDarcyToM2(props.perm_r * log(props.radius[i] / props.r_w) / (log(props.radius[i] / r_w) + props.skins[i]) ) );
		else
			props_sk.perm_eff.push_back( MilliDarcyToM2(props.perm_r) );
	}

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	// Oil properties
	props_oil.visc = cPToPaSec( props.visc_oil );
	props_oil.dens_stc = props.dens_oil_stc;
	props_oil.beta = props.beta_oil;
	
	props_oil.c = props.c_oil;
	props_oil.lambda = props.lambda_oil;
	props_oil.jt = props.jt_oil;
	props_oil.ad = props.ad_oil;

	// Gas properties
	props_gas.visc = cPToPaSec( props.visc_gas );
	props_gas.dens_stc = props.dens_gas_stc;
	
	props_gas.c = props.c_gas;
	props_gas.lambda = props.lambda_gas;
	props_gas.jt = props.jt_gas;
	props_gas.ad = props.ad_gas;

	alpha = props.alpha;
	depth_point = props.depth_point;
	kappa_eff = props.kappa_eff;
	L = props.L;

	makeDimLess();
	
	// Data sets
	props_oil.kr = setDataset(props.kr_oil, 1.0, 1.0);
	props_gas.kr = setDataset(props.kr_gas, 1.0, 1.0);

	props_oil.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_gas.b = setDataset(props.B_gas, P_dim / BAR_TO_PA, 1.0);
	Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	Prs = setInvDataset(props.Rs, 1.0, P_dim / BAR_TO_PA);
}

void TwoPhaseModel::makeDimLess()
{
	// Main units
	R_dim = r_w;
	t_dim = ht;
	P_dim = BAR_TO_PA;
	T_dim = T_init;

	// Temporal properties
	ht = ht / t_dim;
	ht_min = ht_min / t_dim;
	ht_max = ht_max / t_dim;

	// Initial condition
	T_init = T_init / T_dim;
	P_init = P_init / P_dim;
	P_bub = P_bub / P_dim;

	// Grid properties
	r_w = r_w / R_dim;
	r_e = r_e / R_dim;

	// Skeleton properties
	props_sk.perm_r = props_sk.perm_r / R_dim / R_dim;
	props_sk.perm_z = props_sk.perm_z / R_dim / R_dim;
	props_sk.dens_stc = props_sk.dens_stc / P_dim / t_dim / t_dim * R_dim * R_dim;
	props_sk.beta = props_sk.beta * P_dim;
	props_sk.height = props_sk.height / R_dim;
	props_sk.h1 = (props_sk.h1 - depth_point) / R_dim;
	props_sk.h2 = (props_sk.h2 - depth_point) / R_dim;

	props_sk.c = props_sk.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_sk.lambda = props_sk.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for(int i = 0; i < periodsNum; i++)
	{
		period[i] = period[i] / t_dim;
		rate[i] = rate[i] / Q_dim;
		props_sk.r_eff[i] = props_sk.r_eff[i] / R_dim;
		props_sk.perm_eff[i] = props_sk.perm_eff[i] / R_dim / R_dim;
	}

	// Oil properties
	props_oil.visc = props_oil.visc / (P_dim * t_dim);
	props_oil.dens_stc = props_oil.dens_stc / P_dim / t_dim / t_dim * R_dim * R_dim;
	props_oil.beta = props_oil.beta * P_dim;
	
	props_oil.c = props_oil.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_oil.lambda = props_oil.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_oil.jt = props_oil.jt * P_dim / T_dim;
	props_oil.ad = props_oil.ad * P_dim / T_dim;

	// Gas properties
	props_gas.visc = props_gas.visc / (P_dim * t_dim);
	props_gas.dens_stc = props_gas.dens_stc / P_dim / t_dim / t_dim * R_dim * R_dim;
	
	props_gas.c = props_gas.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_gas.lambda = props_gas.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;;
	props_gas.jt = props_gas.jt * P_dim / T_dim;;
	props_gas.ad = props_gas.ad * P_dim / T_dim;;

	// Rest properties
	alpha = alpha / t_dim;
	//depth_point = 0.0;
	kappa_eff = kappa_eff * t_dim / (R_dim * R_dim);
	L = L / R_dim / R_dim * t_dim * t_dim;
}

void TwoPhaseModel::buildGridLog()
{
	cells.reserve( (cellsNum_r + 2) * (cellsNum_z + 2) );

	Volume = 0.0;
	int counter = 0;

	double r_prev = r_w;
	double logMax = log(r_e / r_w);
	double logStep = logMax / (double)cellsNum_r;
	
	double hz = (props_sk.h2 - props_sk.h1) / (double)cellsNum_z;
	double cm_z = props_sk.h1 - hz / 2.0;
	double hr = r_prev * (exp(logStep) - 1.0);
	double cm_r = r_w - hr / 2.0;

	// Left border
	cells.push_back( CylCell(counter++, cm_r, cm_z, hr, hz) );
	//Volume += cells[cells.size()-1].V;
	cm_z += hz;
	for(int i = 0; i < cellsNum_z; i++)
	{
		cells.push_back( CylCell(counter++, cm_r, cm_z, hr, hz) );
		Volume += cells[cells.size()-1].V;
		cm_z += hz;
	}
	cells.push_back( CylCell(counter++, cm_r, cm_z, hr, hz) );
	//Volume += cells[cells.size()-1].V;

	// Middle cells
	for(int j = 0; j < cellsNum_r; j++)
	{
		cm_z = props_sk.h1 - hz / 2.0;
		cm_r = r_prev * (1.0 + exp(logStep)) / 2.0;
		hr = r_prev * (exp(logStep) - 1.0);

		cells.push_back( CylCell(counter++, cm_r, cm_z, hr, hz, false) );
		//Volume += cells[cells.size()-1].V;
		cm_z += hz;
		for(int i = 0; i < cellsNum_z; i++)
		{
			cells.push_back( CylCell(counter++, cm_r, cm_z, hr, hz, false) );
			Volume += cells[cells.size()-1].V;
			cm_z += hz;
		}
		cells.push_back( CylCell(counter++, cm_r, cm_z, hr, hz, false) );
		//Volume += cells[cells.size()-1].V;

		r_prev = r_prev * exp(logStep);
	}

	// Right border
	cm_z = props_sk.h1 - hz / 2.0;
	cm_r = r_e + hr / 2.0;
	
	cells.push_back( CylCell(counter++, cm_r, cm_z, hr, hz) );
	//Volume += cells[cells.size()-1].V;
	cm_z += hz;

	for(int i = 0; i < cellsNum_z; i++)
	{
		cells.push_back( CylCell(counter++, cm_r, cm_z, hr, hz) );
		//Volume += cells[cells.size()-1].V;
		cm_z += hz;
	}
	cells.push_back( CylCell(counter++, cm_r, cm_z, hr, hz) );
	//Volume += cells[cells.size()-1].V;
}

Interpolate* TwoPhaseModel::setDataset(vector< pair<double,double> >& vec, const double xDim, const double yDim)
{
	sort(vec.begin(), vec.end(), sort_pair_first());

	const int N = vec.size();
	double* x = new double [N];
	double* y  = new double [N];

	for (int i = 0; i < N; i++)
	{
		x[i] = vec[i].first / xDim;
		y[i] = vec[i].second / yDim;
	}

	return new Interpolate(x, y, N);
}

Interpolate* TwoPhaseModel::setInvDataset(vector< pair<double,double> >& vec, const double xDim, const double yDim)
{
	sort(vec.begin(), vec.end(), sort_pair_second());

	const int N = vec.size();
	double* x = new double [N];
	double* y  = new double [N];

	for (int i = 0; i < N; i++)
	{
		x[i] = vec[i].second / xDim;
		y[i] = vec[i].first / yDim;
	}

	return new Interpolate(x, y, N);
}

void TwoPhaseModel::setPeriod(int period)
{
	Q_sum = rate[period];
	for(int i = 0; i < Qcell.size(); i++)
		Qcell[i].second = Q_sum * cells[ Qcell[i].first ].hz / height_perf;
	
	r_eff = props_sk.r_eff[period];
	Perm_eff = props_sk.perm_eff[period];
	skin = props_sk.skin[period];
}

void TwoPhaseModel::setPerforated()
{
	height_perf = 0.0;
	for(int i = 1; i <= cellsNum_z; i++)
	{
		Qcell.push_back( make_pair<int,double>(i, 0.0) );
		height_perf += cells[i].hz;
	}
}

void TwoPhaseModel::setInitialState()
{
	vector<CylCell>::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
	{
		it->u_prev.temp = it->u_iter.temp = it->u_next.temp = T_init;
		it->u_prev.pres = it->u_iter.pres = it->u_next.pres = P_init;
		it->u_prev.sat_oil = it->u_iter.sat_oil = it->u_next.sat_oil = S_init;
		it->u_prev.pres_bub = it->u_iter.pres_bub = it->u_next.pres_bub = P_bub;
		if(P_bub >= P_init)
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = true;
		else
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = false;
	}
}

void TwoPhaseModel::calculateBounds()
{
	// Left border
	calculatePerforated();

	// Top & bottom borders
	int idx_top, idx_bot;
	for(int i = 0; i < cellsNum_r; i++)
	{
		idx_top = (i+1) * (cellsNum_z+2);
		idx_bot = idx_top + cellsNum_z + 1;

		CylCell& top = cells[idx_top];
		CylCell& bot = cells[idx_bot];

		top.u_next = cells[idx_top+1].u_next;
		bot.u_next = cells[idx_bot-1].u_next;
	}

	// Right border
	int idx;
	for(int i = 0; i < cellsNum_z+2; i++)
	{
		idx = (cellsNum_r+1) * (cellsNum_z+2) + i;
		CylCell& right = cells[idx];

		right.u_next.temp = 2.0 * T_init - cells[idx-cellsNum_z-2].u_next.temp;
		right.u_next.pres = 2.0 * P_init - cells[idx-cellsNum_z-2].u_next.pres;
		right.u_next.sat_oil = 2.0 * S_init - cells[idx-cellsNum_z-2].u_next.sat_oil;
	}
}

void TwoPhaseModel::calculatePerforated()
{
	for(int i = 0; i < Qcell.size(); i++)
	{
		Variables& nebr = cells[Qcell[i].first + cellsNum_z + 2].u_next;
		Variables& next = cells[Qcell[i].first].u_next;
		Variables& upwd = cells[getUpwindIdx(Qcell[i].first, Qcell[i].first + cellsNum_z + 2)].u_next;

		next.temp = nebr.temp;
		next.pres = nebr.pres - Qcell[i].second * props_oil.visc * getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR) / getKr_oil(upwd.sat_oil) / getTrans(Qcell[i].first, Qcell[i].first + cellsNum_z + 2);
		next.sat_oil = nebr.sat_oil;
	}
}

void TwoPhaseModel::construction_from_fz(int N, int n, int key)
{
	vector<CylCell>::iterator it;
	if (key == PRES)
	{
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < cellsNum_z; j++)
			{
				Variables& var = cells[(i+1)*(cellsNum_z+2) + (j+1)].u_next;

				var.pres = fz[i][2*j+1];
 				var.sat_oil = fz[i][2*j+2];
			}
		}
	}
	else
	{
		//for(int i = 0; i < N; i++)
		//	cells[i].u_next.temp = fz[i][1];
	}
}

int TwoPhaseModel::doNextStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;	int newton_itr = 0;
	double averPresPrev = averValue(1);
	double averSatPrev = averValue(2);
	double averPres, averSat;
	double dAverPres = 1.0, dAverSat = 1.0;

	while( err_newton > 1.e-4 && ( dAverSat > 1.e-8 || dAverPres > 1.e-4) && newton_itr < 8 )
	{	

		calculateBounds();
		Solve(cellsNum_r-1, 2*cellsNum_z, PRES);
		copyIterLayer();
		construction_from_fz(cellsNum_r, 2*cellsNum_z, PRES);
		solveP_bub();
		calculateBounds();

		err_newton = convergance(PRES, cellIdx, varIdx);

		averPres = averValue(1);					averSat = averValue(2);
		dAverPres = fabs(averPres - averPresPrev);	dAverSat = fabs(averSat - averSatPrev);
		averPresPrev = averPres;					averSatPrev = averSat;

		if(varIdx == PRES)
			cout << "BadPresValue[" << cellIdx  << "]: " << cells[cellIdx].u_next.pres << endl;
		else if(varIdx == SAT)
			cout << "BadSatValue[" << cellIdx  << "]: " << cells[cellIdx].u_next.sat_oil << endl;

		newton_itr++;
	}
	cout << "Newton Iterations = " << newton_itr << endl;

	return newton_itr;
}

/*
void TwoPhaseModel::LeftBoundAppr(int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for (int j = 0 ; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		B[i][i] = -1.0;
		A[i][i] = 1.0;
		RightSide[i][0] = 0.0;
	}

	if(key == PRES)
	{
		int idx;
		vector< pair<int,double> >::iterator it;
		for(it = Qcell.begin(); it != Qcell.end(); ++it)
		{
			CylCell& cell = cells[it->first];
			idx = 2 * it->first;
			
			// First eqn
			B[idx][idx] = -1.0;
			A[idx][idx] = 1.0;
			RightSide[idx][0] = it->second * (cells[it->first+cellsNum_z].r - cell.r) * props_oil.visc * 
								getB_oil(cell.u_next.pres, cell.u_next.pres_bub, cell.u_next.SATUR) / 2.0 / M_PI / cell.r / cell.hz /
								Perm_eff / getKr_oil(cell.u_next.sat_oil);

			// Second eqn
			B[idx+1][idx+1] = -1.0;
			A[idx+1][idx+1] = 1.0;
			RightSide[idx+1][0] = 0.0;


		}
	}
	construction_bz(MZ, 2);
}

void TwoPhaseModel::RightBoundAppr(int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for (int j = 0 ; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		RightSide[i][0] = 0.0;
	}

	if(key == PRES)
	{
		int idx = 0;
		for(int i = 0; i < cellsNum_z; i++)
		{
			idx = 2 * i;

			// First eqn
			B[idx][idx] = 1.0;
			RightSide[idx][0] = P_init;

			// Second eqn
			C[idx+1][idx+1] = 1.0;
			B[idx+1][idx+1] = -1.0;
		}
	}
	construction_bz(MZ,1);
}

*/

void TwoPhaseModel::LeftBoundAppr(int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for (int j = 0 ; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		B[i][i] = 0.0;
		A[i][i] = 0.0;
		RightSide[i][0] = 0.0;
	}

	if(key == PRES)
	{
		int idx = 0;
		int i = (cellsNum_z+2) + 1;
		
		// Top cell
		// First eqn
		//C[idx][idx] = solve_eq1_dp_beta(i, i - cellsNum_z - 2);
		//C[idx][idx+1] = solve_eq1_ds_beta(i, i - cellsNum_z - 2);
		//B[idx][idx-2] = solve_eq1_dp_beta(i, i-1);
		//B[idx][idx-1] = solve_eq1_ds_beta(i, i-1);
		C[idx][idx] = solve_eq1_dp(i);
		C[idx][idx+1] = solve_eq1_ds(i);
		C[idx][idx+2] = solve_eq1_dp_beta(i, i+1);
		C[idx][idx+3] = solve_eq1_ds_beta(i, i+1);
		B[idx][idx] = solve_eq1_dp_beta(i, i + cellsNum_z + 2);
		B[idx][idx+1] = solve_eq1_ds_beta(i, i + cellsNum_z + 2);
		RightSide[idx][0] = -solve_eq1(i) + 
							solve_eq1_dp_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.pres - cells[i-cellsNum_z-2].u_iter.pres) +
							solve_eq1_ds_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.sat_oil - cells[i-cellsNum_z-2].u_iter.sat_oil) + 
							solve_eq1_dp_beta(i, i-1) * (cells[i-1].u_iter.pres - cells[i-1].u_next.pres) + 
							solve_eq1_ds_beta(i, i-1) * (cells[i-1].u_iter.sat_oil - cells[i-1].u_next.sat_oil) + 
							C[idx][idx] * cells[i].u_next.pres + C[idx][idx+1] * cells[i].u_next.sat_oil + 
							C[idx][idx+2] * cells[i+1].u_next.pres + C[idx][idx+3] * cells[i+1].u_next.sat_oil + 
							B[idx][idx] * cells[i+cellsNum_z+2].u_next.pres + B[idx][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			
		// Second eqn
		//C[idx+1][idx] = solve_eq2_dp_beta(i, i - cellsNum_z - 2);
		//C[idx+1][idx+1] = solve_eq2_ds_beta(i, i - cellsNum_z - 2);
		//B[idx+1][idx-2] = solve_eq2_dp_beta(i, i-1);
		//B[idx+1][idx-1] = solve_eq2_ds_beta(i, i-1);
		C[idx+1][idx] = solve_eq2_dp(i);
		C[idx+1][idx+1] = solve_eq2_ds(i);
		C[idx+1][idx+2] = solve_eq2_dp_beta(i, i+1);
		C[idx+1][idx+3] = solve_eq2_ds_beta(i, i+1);
		B[idx+1][idx] = solve_eq2_dp_beta(i, i + cellsNum_z + 2);
		B[idx+1][idx+1] = solve_eq2_ds_beta(i, i + cellsNum_z + 2);
		RightSide[idx+1][0] = -solve_eq2(i) + 
							solve_eq2_dp_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.pres - cells[i-cellsNum_z-2].u_iter.pres) +
							solve_eq2_ds_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.sat_oil - cells[i-cellsNum_z-2].u_iter.sat_oil) + 
							solve_eq2_dp_beta(i, i-1) * (cells[i-1].u_iter.pres - cells[i-1].u_next.pres) + 
							solve_eq2_ds_beta(i, i-1) * (cells[i-1].u_iter.sat_oil - cells[i-1].u_next.sat_oil) + 
							C[idx+1][idx] * cells[i].u_next.pres + C[idx+1][idx+1] * cells[i].u_next.sat_oil + 
							C[idx+1][idx+2] * cells[i+1].u_next.pres + C[idx+1][idx+3] * cells[i+1].u_next.sat_oil + 
							B[idx+1][idx] * cells[i+cellsNum_z+2].u_next.pres + B[idx+1][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			
		idx+=2;

		// Middle cells
		for(i = (cellsNum_z+2) + 2; i < 2 * (cellsNum_z+2) - 2; i++)
		{
			Variables& next = cells[i].u_next;

			// First eqn
			//C[idx][idx] = solve_eq1_dp_beta(i, i - cellsNum_z - 2);
			//C[idx][idx+1] = solve_eq1_ds_beta(i, i - cellsNum_z - 2);
			C[idx][idx-2] = solve_eq1_dp_beta(i, i-1);
			C[idx][idx-1] = solve_eq1_ds_beta(i, i-1);
			C[idx][idx] = solve_eq1_dp(i);
			C[idx][idx+1] = solve_eq1_ds(i);
			C[idx][idx+2] = solve_eq1_dp_beta(i, i+1);
			C[idx][idx+3] = solve_eq1_ds_beta(i, i+1);
			B[idx][idx] = solve_eq1_dp_beta(i, i + cellsNum_z + 2);
			B[idx][idx+1] = solve_eq1_ds_beta(i, i + cellsNum_z + 2);
			RightSide[idx][0] = -solve_eq1(i) + 
								solve_eq1_dp_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.pres - cells[i-cellsNum_z-2].u_iter.pres) +
								solve_eq1_ds_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.sat_oil - cells[i-cellsNum_z-2].u_iter.sat_oil) + 
								C[idx][idx-2] * cells[i-1].u_next.pres + C[idx][idx-1] * cells[i-1].u_next.sat_oil +
								C[idx][idx] * cells[i].u_next.pres + C[idx][idx+1] * cells[i].u_next.sat_oil + 
								C[idx][idx+2] * cells[i+1].u_next.pres + C[idx][idx+3] * cells[i+1].u_next.sat_oil + 
								B[idx][idx] * cells[i+cellsNum_z+2].u_next.pres + B[idx][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			
			// Second eqn
			//C[idx+1][idx] = solve_eq2_dp_beta(i, i - cellsNum_z - 2);
			//C[idx+1][idx+1] = solve_eq2_ds_beta(i, i - cellsNum_z - 2);
			C[idx+1][idx-2] = solve_eq2_dp_beta(i, i-1);
			C[idx+1][idx-1] = solve_eq2_ds_beta(i, i-1);
			C[idx+1][idx] = solve_eq2_dp(i);
			C[idx+1][idx+1] = solve_eq2_ds(i);
			C[idx+1][idx+2] = solve_eq2_dp_beta(i, i+1);
			C[idx+1][idx+3] = solve_eq2_ds_beta(i, i+1);
			B[idx+1][idx] = solve_eq2_dp_beta(i, i + cellsNum_z + 2);
			B[idx+1][idx+1] = solve_eq2_ds_beta(i, i + cellsNum_z + 2);
			RightSide[idx+1][0] = -solve_eq2(i) + 
								solve_eq2_dp_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.pres - cells[i-cellsNum_z-2].u_iter.pres) +
								solve_eq2_ds_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.sat_oil - cells[i-cellsNum_z-2].u_iter.sat_oil) +
								C[idx+1][idx-2] * cells[i-1].u_next.pres + C[idx+1][idx-1] * cells[i-1].u_next.sat_oil +
								C[idx+1][idx] * cells[i].u_next.pres + C[idx+1][idx+1] * cells[i].u_next.sat_oil + 
								C[idx+1][idx+2] * cells[i+1].u_next.pres + C[idx+1][idx+3] * cells[i+1].u_next.sat_oil + 
								B[idx+1][idx] * cells[i+cellsNum_z+2].u_next.pres + B[idx+1][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			
			idx += 2;
		}

		// Bottom cell
		// First eqn
		//C[idx][idx] = solve_eq1_dp_beta(i, i - cellsNum_z - 2);
		//C[idx][idx+1] = solve_eq1_ds_beta(i, i - cellsNum_z - 2);
		C[idx][idx-2] = solve_eq1_dp_beta(i, i-1);
		C[idx][idx-1] = solve_eq1_ds_beta(i, i-1);
		C[idx][idx] = solve_eq1_dp(i);
		C[idx][idx+1] = solve_eq1_ds(i);
		//B[idx][idx+2] = solve_eq1_dp_beta(i, i+1);
		//B[idx][idx+3] = solve_eq1_ds_beta(i, i+1);
		B[idx][idx] = solve_eq1_dp_beta(i, i + cellsNum_z + 2);
		B[idx][idx+1] = solve_eq1_ds_beta(i, i + cellsNum_z + 2);
		RightSide[idx][0] = -solve_eq1(i) + 
							solve_eq1_dp_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.pres - cells[i-cellsNum_z-2].u_iter.pres) +
							solve_eq1_ds_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.sat_oil - cells[i-cellsNum_z-2].u_iter.sat_oil) + 
							C[idx][idx-2] * cells[i-1].u_next.pres + C[idx][idx-1] * cells[i-1].u_next.sat_oil +
							C[idx][idx] * cells[i].u_next.pres + C[idx][idx+1] * cells[i].u_next.sat_oil + 
							solve_eq1_dp_beta(i, i+1) * (cells[i+1].u_iter.pres - cells[i+1].u_next.pres) + 
							solve_eq1_ds_beta(i, i+1) * (cells[i+1].u_iter.sat_oil - cells[i+1].u_next.sat_oil) + 
							B[idx][idx] * cells[i+cellsNum_z+2].u_next.pres + B[idx][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			
		// Second eqn
		//C[idx+1][idx] = solve_eq2_dp_beta(i, i - cellsNum_z - 2);
		//C[idx+1][idx+1] = solve_eq2_ds_beta(i, i - cellsNum_z - 2);
		C[idx+1][idx-2] = solve_eq2_dp_beta(i, i-1);
		C[idx+1][idx-1] = solve_eq2_ds_beta(i, i-1);
		C[idx+1][idx] = solve_eq2_dp(i);
		C[idx+1][idx+1] = solve_eq2_ds(i);
		//B[idx+1][idx+2] = solve_eq2_dp_beta(i, i+1);
		//B[idx+1][idx+3] = solve_eq2_ds_beta(i, i+1);
		B[idx+1][idx] = solve_eq2_dp_beta(i, i + cellsNum_z + 2);
		B[idx+1][idx+1] = solve_eq2_ds_beta(i, i + cellsNum_z + 2);
		RightSide[idx+1][0] = -solve_eq2(i) + 
							solve_eq2_dp_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.pres - cells[i-cellsNum_z-2].u_iter.pres) +
							solve_eq2_ds_beta(i, i - cellsNum_z - 2) * (cells[i-cellsNum_z-2].u_iter.sat_oil - cells[i-cellsNum_z-2].u_iter.sat_oil) + 
							C[idx+1][idx-2] * cells[i-1].u_next.pres + C[idx+1][idx-1] * cells[i-1].u_next.sat_oil +
							C[idx+1][idx] * cells[i].u_next.pres + C[idx+1][idx+1] * cells[i].u_next.sat_oil + 
							solve_eq2_dp_beta(i, i+1) * (cells[i+1].u_iter.pres - cells[i+1].u_next.pres) + 
							solve_eq2_ds_beta(i, i+1) * (cells[i+1].u_iter.sat_oil - cells[i+1].u_next.sat_oil) + 
							B[idx+1][idx] * cells[i+cellsNum_z+2].u_next.pres + B[idx+1][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;
	}

	construction_bz(MZ, 2);
}

void TwoPhaseModel::RightBoundAppr(int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for (int j = 0 ; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		RightSide[i][0] = 0.0;
	}

	if(key == PRES)
	{
		int idx = 0;
		int i = cellsNum_r * (cellsNum_z+2) + 1;
		
		// Top cell
		// First eqn
		B[idx][idx] = solve_eq1_dp_beta(i, i - cellsNum_z - 2);
		B[idx][idx+1] = solve_eq1_ds_beta(i, i - cellsNum_z - 2);
		//B[idx][idx-2] = solve_eq1_dp_beta(i, i-1);
		//B[idx][idx-1] = solve_eq1_ds_beta(i, i-1);
		A[idx][idx] = solve_eq1_dp(i);
		A[idx][idx+1] = solve_eq1_ds(i);
		A[idx][idx+2] = solve_eq1_dp_beta(i, i+1);
		A[idx][idx+3] = solve_eq1_ds_beta(i, i+1);
		//A[idx][idx] = solve_eq1_dp_beta(i, i + cellsNum_z + 2);
		//A[idx][idx+1] = solve_eq1_ds_beta(i, i + cellsNum_z + 2);
		RightSide[idx][0] = -solve_eq1(i) + 
							B[idx][idx] * cells[i-cellsNum_z-2].u_next.pres + B[idx][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
							solve_eq1_dp_beta(i, i-1) * (cells[i-1].u_iter.pres - cells[i-1].u_next.pres) + 
							solve_eq1_ds_beta(i, i-1) * (cells[i-1].u_iter.sat_oil - cells[i-1].u_next.sat_oil) + 
							A[idx][idx] * cells[i].u_next.pres + A[idx][idx+1] * cells[i].u_next.sat_oil + 
							A[idx][idx+2] * cells[i+1].u_next.pres + A[idx][idx+3] * cells[i+1].u_next.sat_oil + 
							solve_eq1_dp_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.pres - cells[i+cellsNum_z+2].u_next.pres ) +
							solve_eq1_ds_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.sat_oil - cells[i+cellsNum_z+2].u_next.sat_oil );
		// Second eqn
		B[idx+1][idx] = solve_eq2_dp_beta(i, i - cellsNum_z - 2);
		B[idx+1][idx+1] = solve_eq2_ds_beta(i, i - cellsNum_z - 2);
		//B[idx+1][idx-2] = solve_eq2_dp_beta(i, i-1);
		//B[idx+1][idx-1] = solve_eq2_ds_beta(i, i-1);
		A[idx+1][idx] = solve_eq2_dp(i);
		A[idx+1][idx+1] = solve_eq2_ds(i);
		A[idx+1][idx+2] = solve_eq2_dp_beta(i, i+1);
		A[idx+1][idx+3] = solve_eq2_ds_beta(i, i+1);
		//A[idx+1][idx] = solve_eq2_dp_beta(i, i + cellsNum_z + 2);
		//A[idx+1][idx+1] = solve_eq2_ds_beta(i, i + cellsNum_z + 2);
		RightSide[idx+1][0] = -solve_eq2(i) + 
							B[idx+1][idx] * cells[i-cellsNum_z-2].u_next.pres + B[idx+1][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
							solve_eq2_dp_beta(i, i-1) * (cells[i-1].u_iter.pres - cells[i-1].u_next.pres) + 
							solve_eq2_ds_beta(i, i-1) * (cells[i-1].u_iter.sat_oil - cells[i-1].u_next.sat_oil) + 
							A[idx+1][idx] * cells[i].u_next.pres + A[idx+1][idx+1] * cells[i].u_next.sat_oil + 
							A[idx+1][idx+2] * cells[i+1].u_next.pres + A[idx+1][idx+3] * cells[i+1].u_next.sat_oil + 
							solve_eq2_dp_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.pres - cells[i+cellsNum_z+2].u_next.pres ) +
							solve_eq2_ds_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.sat_oil - cells[i+cellsNum_z+2].u_next.sat_oil );

		idx+=2;

		// Middle cells
		for(i = cellsNum_r * (cellsNum_z+2) + 2; i < (cellsNum_r+1) * (cellsNum_z+2) - 2; i++)
		{
			Variables& next = cells[i].u_next;

			// First eqn
			B[idx][idx] = solve_eq1_dp_beta(i, i - cellsNum_z - 2);
			B[idx][idx+1] = solve_eq1_ds_beta(i, i - cellsNum_z - 2);
			A[idx][idx-2] = solve_eq1_dp_beta(i, i-1);
			A[idx][idx-1] = solve_eq1_ds_beta(i, i-1);
			A[idx][idx] = solve_eq1_dp(i);
			A[idx][idx+1] = solve_eq1_ds(i);
			A[idx][idx+2] = solve_eq1_dp_beta(i, i+1);
			A[idx][idx+3] = solve_eq1_ds_beta(i, i+1);
			//A[idx][idx] = solve_eq1_dp_beta(i, i + cellsNum_z + 2);
			//A[idx][idx+1] = solve_eq1_ds_beta(i, i + cellsNum_z + 2);
			RightSide[idx][0] = -solve_eq1(i) + 
								B[idx][idx] * cells[i-cellsNum_z-2].u_next.pres + B[idx][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
								A[idx][idx-2] * cells[i-1].u_next.pres + A[idx][idx-1] * cells[i-1].u_next.sat_oil +
								A[idx][idx] * cells[i].u_next.pres + A[idx][idx+1] * cells[i].u_next.sat_oil + 
								A[idx][idx+2] * cells[i+1].u_next.pres + A[idx][idx+3] * cells[i+1].u_next.sat_oil + 
								solve_eq1_dp_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.pres - cells[i+cellsNum_z+2].u_next.pres ) +
								solve_eq1_ds_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.sat_oil - cells[i+cellsNum_z+2].u_next.sat_oil );
			// Second eqn
			B[idx+1][idx] = solve_eq2_dp_beta(i, i - cellsNum_z - 2);
			B[idx+1][idx+1] = solve_eq2_ds_beta(i, i - cellsNum_z - 2);
			A[idx+1][idx-2] = solve_eq2_dp_beta(i, i-1);
			A[idx+1][idx-1] = solve_eq2_ds_beta(i, i-1);
			A[idx+1][idx] = solve_eq2_dp(i);
			A[idx+1][idx+1] = solve_eq2_ds(i);
			A[idx+1][idx+2] = solve_eq2_dp_beta(i, i+1);
			A[idx+1][idx+3] = solve_eq2_ds_beta(i, i+1);
			//A[idx+1][idx] = solve_eq2_dp_beta(i, i + cellsNum_z + 2);
			//A[idx+1][idx+1] = solve_eq2_ds_beta(i, i + cellsNum_z + 2);
			RightSide[idx+1][0] = -solve_eq2(i) + 
								B[idx+1][idx] * cells[i-cellsNum_z-2].u_next.pres + B[idx+1][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
								A[idx+1][idx-2] * cells[i-1].u_next.pres + A[idx+1][idx-1] * cells[i-1].u_next.sat_oil +
								A[idx+1][idx] * cells[i].u_next.pres + A[idx+1][idx+1] * cells[i].u_next.sat_oil + 
								A[idx+1][idx+2] * cells[i+1].u_next.pres + A[idx+1][idx+3] * cells[i+1].u_next.sat_oil + 
								solve_eq2_dp_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.pres - cells[i+cellsNum_z+2].u_next.pres ) +
								solve_eq2_ds_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.sat_oil - cells[i+cellsNum_z+2].u_next.sat_oil );
			idx += 2;
		}

		// Bottom cell
		// First eqn
		B[idx][idx] = solve_eq1_dp_beta(i, i - cellsNum_z - 2);
		B[idx][idx+1] = solve_eq1_ds_beta(i, i - cellsNum_z - 2);
		A[idx][idx-2] = solve_eq1_dp_beta(i, i-1);
		A[idx][idx-1] = solve_eq1_ds_beta(i, i-1);
		A[idx][idx] = solve_eq1_dp(i);
		A[idx][idx+1] = solve_eq1_ds(i);
		//B[idx][idx+2] = solve_eq1_dp_beta(i, i+1);
		//B[idx][idx+3] = solve_eq1_ds_beta(i, i+1);
		//A[idx][idx] = solve_eq1_dp_beta(i, i + cellsNum_z + 2);
		//A[idx][idx+1] = solve_eq1_ds_beta(i, i + cellsNum_z + 2);
		RightSide[idx][0] = -solve_eq1(i) + 
							B[idx][idx] * cells[i-cellsNum_z-2].u_next.pres + B[idx][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
							A[idx][idx-2] * cells[i-1].u_next.pres + A[idx][idx-1] * cells[i-1].u_next.sat_oil +
							A[idx][idx] * cells[i].u_next.pres + A[idx][idx+1] * cells[i].u_next.sat_oil + 
							solve_eq1_dp_beta(i, i+1) * (cells[i+1].u_iter.pres - cells[i+1].u_next.pres) + 
							solve_eq1_ds_beta(i, i+1) * (cells[i+1].u_iter.sat_oil - cells[i+1].u_next.sat_oil) + 
							solve_eq1_dp_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.pres - cells[i+cellsNum_z+2].u_next.pres ) +
							solve_eq1_ds_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.sat_oil - cells[i+cellsNum_z+2].u_next.sat_oil );
		// Second eqn
		B[idx+1][idx] = solve_eq2_dp_beta(i, i - cellsNum_z - 2);
		B[idx+1][idx+1] = solve_eq2_ds_beta(i, i - cellsNum_z - 2);
		A[idx+1][idx-2] = solve_eq2_dp_beta(i, i-1);
		A[idx+1][idx-1] = solve_eq2_ds_beta(i, i-1);
		A[idx+1][idx] = solve_eq2_dp(i);
		A[idx+1][idx+1] = solve_eq2_ds(i);
		//B[idx+1][idx+2] = solve_eq2_dp_beta(i, i+1);
		//B[idx+1][idx+3] = solve_eq2_ds_beta(i, i+1);
		//A[idx+1][idx] = solve_eq2_dp_beta(i, i + cellsNum_z + 2);
		//A[idx+1][idx+1] = solve_eq2_ds_beta(i, i + cellsNum_z + 2);
		RightSide[idx+1][0] = -solve_eq2(i) + 
							B[idx+1][idx] * cells[i-cellsNum_z-2].u_next.pres + B[idx+1][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
							A[idx+1][idx-2] * cells[i-1].u_next.pres + A[idx+1][idx-1] * cells[i-1].u_next.sat_oil +
							A[idx+1][idx] * cells[i].u_next.pres + A[idx+1][idx+1] * cells[i].u_next.sat_oil + 
							solve_eq2_dp_beta(i, i+1) * (cells[i+1].u_iter.pres - cells[i+1].u_next.pres) + 
							solve_eq2_ds_beta(i, i+1) * (cells[i+1].u_iter.sat_oil - cells[i+1].u_next.sat_oil) + 
							solve_eq2_dp_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.pres - cells[i+cellsNum_z+2].u_next.pres ) +
							solve_eq2_ds_beta(i, i + cellsNum_z + 2) * (cells[i+cellsNum_z+2].u_iter.sat_oil - cells[i+cellsNum_z+2].u_next.sat_oil );
	}

/*	for(int i = 0; i < MZ; i++)
	{
		for(int j = 0; j < MZ; j++)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
		//cout << RightSide[i][0] << " ";
	}
	cout << endl;*/

	construction_bz(MZ,1);
}

void TwoPhaseModel::MiddleAppr(int current, int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for(int j = 0; j < MZ; j++)
		{
			A[i][j] = 0.0;
			B[i][j] = 0.0;
			C[i][j] = 0.0;
		}
		RightSide[i][0] = 0.0;
	}

	if(key == PRES)
	{
		int idx = 0;
		int i = (current+1) * (cellsNum_z+2) + 1;
		
		// Top cell
		TopAppr(i);
		idx+=2;

		// Middle cells
		for(i = (current+1) * (cellsNum_z+2) + 2; i < (current+2) * (cellsNum_z+2) - 2; i++)
		{
			Variables& next = cells[i].u_next;

			// First eqn
			C[idx][idx] = solve_eq1_dp_beta(i, i - cellsNum_z - 2);
			C[idx][idx+1] = solve_eq1_ds_beta(i, i - cellsNum_z - 2);
			B[idx][idx-2] = solve_eq1_dp_beta(i, i-1);
			B[idx][idx-1] = solve_eq1_ds_beta(i, i-1);
			B[idx][idx] = solve_eq1_dp(i);
			B[idx][idx+1] = solve_eq1_ds(i);
			B[idx][idx+2] = solve_eq1_dp_beta(i, i+1);
			B[idx][idx+3] = solve_eq1_ds_beta(i, i+1);
			A[idx][idx] = solve_eq1_dp_beta(i, i + cellsNum_z + 2);
			A[idx][idx+1] = solve_eq1_ds_beta(i, i + cellsNum_z + 2);
			RightSide[idx][0] = -solve_eq1(i) + 
								C[idx][idx] * cells[i-cellsNum_z-2].u_next.pres + C[idx][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
								B[idx][idx-2] * cells[i-1].u_next.pres + B[idx][idx-1] * cells[i-1].u_next.sat_oil +
								B[idx][idx] * cells[i].u_next.pres + B[idx][idx+1] * cells[i].u_next.sat_oil + 
								B[idx][idx+2] * cells[i+1].u_next.pres + B[idx][idx+3] * cells[i+1].u_next.sat_oil + 
								A[idx][idx] * cells[i+cellsNum_z+2].u_next.pres + A[idx][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			
			// Second eqn
			C[idx+1][idx] = solve_eq2_dp_beta(i, i - cellsNum_z - 2);
			C[idx+1][idx+1] = solve_eq2_ds_beta(i, i - cellsNum_z - 2);
			B[idx+1][idx-2] = solve_eq2_dp_beta(i, i-1);
			B[idx+1][idx-1] = solve_eq2_ds_beta(i, i-1);
			B[idx+1][idx] = solve_eq2_dp(i);
			B[idx+1][idx+1] = solve_eq2_ds(i);
			B[idx+1][idx+2] = solve_eq2_dp_beta(i, i+1);
			B[idx+1][idx+3] = solve_eq2_ds_beta(i, i+1);
			A[idx+1][idx] = solve_eq2_dp_beta(i, i + cellsNum_z + 2);
			A[idx+1][idx+1] = solve_eq2_ds_beta(i, i + cellsNum_z + 2);
			RightSide[idx+1][0] = -solve_eq2(i) + 
								C[idx+1][idx] * cells[i-cellsNum_z-2].u_next.pres + C[idx+1][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
								B[idx+1][idx-2] * cells[i-1].u_next.pres + B[idx+1][idx-1] * cells[i-1].u_next.sat_oil +
								B[idx+1][idx] * cells[i].u_next.pres + B[idx+1][idx+1] * cells[i].u_next.sat_oil + 
								B[idx+1][idx+2] * cells[i+1].u_next.pres + B[idx+1][idx+3] * cells[i+1].u_next.sat_oil + 
								A[idx+1][idx] * cells[i+cellsNum_z+2].u_next.pres + A[idx+1][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			
			idx += 2;
		}

		// Bottom cell
		BottomAppr(i);
	}

	/*for(int i = 0; i < MZ; i++)
	{
		for(int j = 0; j < MZ; j++)
		{
			cout << B[i][j] << " ";
		}
		//cout << RightSide[i][0];
		cout << endl;
	}
	cout << endl;*/
	/*for(int i = 0; i < MZ; i++)
	{
		for(int j = 0; j < MZ; j++)
		{
			cout << B[i][j] << " ";
		}
		cout << endl;
	}*/

	construction_bz(MZ,2);
}

void TwoPhaseModel::TopAppr(int i)
{
	int idx = 0;

	// First eqn
	C[idx][idx] = solve_eq1_dp_beta(i, i - cellsNum_z - 2);
	C[idx][idx+1] = solve_eq1_ds_beta(i, i - cellsNum_z - 2);
	//B[idx][idx-2] = solve_eq1_dp_beta(i, i-1);
	//B[idx][idx-1] = solve_eq1_ds_beta(i, i-1);
	B[idx][idx] = solve_eq1_dp(i);
	B[idx][idx+1] = solve_eq1_ds(i);
	B[idx][idx+2] = solve_eq1_dp_beta(i, i+1);
	B[idx][idx+3] = solve_eq1_ds_beta(i, i+1);
	A[idx][idx] = solve_eq1_dp_beta(i, i + cellsNum_z + 2);
	A[idx][idx+1] = solve_eq1_ds_beta(i, i + cellsNum_z + 2);
	RightSide[idx][0] = -solve_eq1(i) + 
						C[idx][idx] * cells[i-cellsNum_z-2].u_next.pres + C[idx][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
						solve_eq1_dp_beta(i, i-1) * (cells[i-1].u_iter.pres - cells[i-1].u_next.pres) + 
						solve_eq1_ds_beta(i, i-1) * (cells[i-1].u_iter.sat_oil - cells[i-1].u_next.sat_oil) + 
						B[idx][idx] * cells[i].u_next.pres + B[idx][idx+1] * cells[i].u_next.sat_oil + 
						B[idx][idx+2] * cells[i+1].u_next.pres + B[idx][idx+3] * cells[i+1].u_next.sat_oil + 
						A[idx][idx] * cells[i+cellsNum_z+2].u_next.pres + A[idx][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			
	// Second eqn
	C[idx+1][idx] = solve_eq2_dp_beta(i, i - cellsNum_z - 2);
	C[idx+1][idx+1] = solve_eq2_ds_beta(i, i - cellsNum_z - 2);
	//B[idx+1][idx-2] = solve_eq2_dp_beta(i, i-1);
	//B[idx+1][idx-1] = solve_eq2_ds_beta(i, i-1);
	B[idx+1][idx] = solve_eq2_dp(i);
	B[idx+1][idx+1] = solve_eq2_ds(i);
	B[idx+1][idx+2] = solve_eq2_dp_beta(i, i+1);
	B[idx+1][idx+3] = solve_eq2_ds_beta(i, i+1);
	A[idx+1][idx] = solve_eq2_dp_beta(i, i + cellsNum_z + 2);
	A[idx+1][idx+1] = solve_eq2_ds_beta(i, i + cellsNum_z + 2);
	RightSide[idx+1][0] = -solve_eq2(i) + 
						C[idx+1][idx] * cells[i-cellsNum_z-2].u_next.pres + C[idx+1][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
						solve_eq2_dp_beta(i, i-1) * (cells[i-1].u_iter.pres - cells[i-1].u_next.pres) + 
						solve_eq2_ds_beta(i, i-1) * (cells[i-1].u_iter.sat_oil - cells[i-1].u_next.sat_oil) + 
						B[idx+1][idx] * cells[i].u_next.pres + B[idx+1][idx+1] * cells[i].u_next.sat_oil + 
						B[idx+1][idx+2] * cells[i+1].u_next.pres + B[idx+1][idx+3] * cells[i+1].u_next.sat_oil + 
						A[idx+1][idx] * cells[i+cellsNum_z+2].u_next.pres + A[idx+1][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			

	/*// First eqn
	B[0][0] = 1.0;
	B[0][2] = -1.0;

	// Second eqn
	B[1][1] = 1.0;
	B[1][3] = -1.0;*/
}

void TwoPhaseModel::BottomAppr(int i)
{
	int idx = 2 * cellsNum_z - 2;
	
	C[idx][idx] = solve_eq1_dp_beta(i, i - cellsNum_z - 2);
	C[idx][idx+1] = solve_eq1_ds_beta(i, i - cellsNum_z - 2);
	B[idx][idx-2] = solve_eq1_dp_beta(i, i-1);
	B[idx][idx-1] = solve_eq1_ds_beta(i, i-1);
	B[idx][idx] = solve_eq1_dp(i);
	B[idx][idx+1] = solve_eq1_ds(i);
	//B[idx][idx+2] = solve_eq1_dp_beta(i, i+1);
	//B[idx][idx+3] = solve_eq1_ds_beta(i, i+1);
	A[idx][idx] = solve_eq1_dp_beta(i, i + cellsNum_z + 2);
	A[idx][idx+1] = solve_eq1_ds_beta(i, i + cellsNum_z + 2);
	RightSide[idx][0] = -solve_eq1(i) + 
						C[idx][idx] * cells[i-cellsNum_z-2].u_next.pres + C[idx][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
						B[idx][idx-2] * cells[i-1].u_next.pres + B[idx][idx-1] * cells[i-1].u_next.sat_oil +
						B[idx][idx] * cells[i].u_next.pres + B[idx][idx+1] * cells[i].u_next.sat_oil + 
						solve_eq1_dp_beta(i, i+1) * (cells[i+1].u_iter.pres - cells[i+1].u_next.pres) + 
						solve_eq1_ds_beta(i, i+1) * (cells[i+1].u_iter.sat_oil - cells[i+1].u_next.sat_oil) + 
						A[idx][idx] * cells[i+cellsNum_z+2].u_next.pres + A[idx][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			
	
	// Second eqn
	C[idx+1][idx] = solve_eq2_dp_beta(i, i - cellsNum_z - 2);
	C[idx+1][idx+1] = solve_eq2_ds_beta(i, i - cellsNum_z - 2);
	B[idx+1][idx-2] = solve_eq2_dp_beta(i, i-1);
	B[idx+1][idx-1] = solve_eq2_ds_beta(i, i-1);
	B[idx+1][idx] = solve_eq2_dp(i);
	B[idx+1][idx+1] = solve_eq2_ds(i);
	//B[idx+1][idx+2] = solve_eq2_dp_beta(i, i+1);
	//B[idx+1][idx+3] = solve_eq2_ds_beta(i, i+1);
	A[idx+1][idx] = solve_eq2_dp_beta(i, i + cellsNum_z + 2);
	A[idx+1][idx+1] = solve_eq2_ds_beta(i, i + cellsNum_z + 2);
	RightSide[idx+1][0] = -solve_eq2(i) + 
						C[idx+1][idx] * cells[i-cellsNum_z-2].u_next.pres + C[idx+1][idx+1] * cells[i-cellsNum_z-2].u_next.sat_oil + 
						B[idx+1][idx-2] * cells[i-1].u_next.pres + B[idx+1][idx-1] * cells[i-1].u_next.sat_oil +
						B[idx+1][idx] * cells[i].u_next.pres + B[idx+1][idx+1] * cells[i].u_next.sat_oil + 
						solve_eq2_dp_beta(i, i+1) * (cells[i+1].u_iter.pres - cells[i+1].u_next.pres) + 
						solve_eq2_ds_beta(i, i+1) * (cells[i+1].u_iter.sat_oil - cells[i+1].u_next.sat_oil) + 
						A[idx+1][idx] * cells[i+cellsNum_z+2].u_next.pres + A[idx+1][idx+1] * cells[i+cellsNum_z+2].u_next.sat_oil;			

	/*// First eqn
	B[idx][idx] = 1.0;
	B[idx][idx-2] = -1.0;

	// Second eqn
	B[idx+1][idx+1] = 1.0;
	B[idx+1][idx-1] = -1.0;*/
}

double TwoPhaseModel::solve_eq1(int cur)
{
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	CylCell& cell = cells[cur];
	Variables& next = cell.u_next;
	Variables& prev = cell.u_prev;
	double H = 0.0;
	
	H = ( getPoro(next.pres) * next.sat_oil / getB_oil(next.pres, next.pres_bub, next.SATUR) - 
				getPoro(prev.pres) * prev.sat_oil / getB_oil(prev.pres, prev.pres_bub, prev.SATUR) );

	for(int i = 0; i < 4; i++)
	{
		Variables& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;

		H += ht / cell.V * getTrans(cur, neighbor[i]) * (next.pres - cells[ neighbor[i] ].u_next.pres) *
			getKr_oil(upwd.sat_oil) / props_oil.visc / getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR);
	}

	return H;
}

double TwoPhaseModel::solve_eq1_dp(int cur)
{
	double upwind;
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	CylCell& cell = cells[cur];
	Variables& next = cell.u_next;
	double Boil_upwd;
	double Boil = getB_oil(next.pres, next.pres_bub, next.SATUR);
	double H = 0.0;

	H = (next.sat_oil * getPoro_dp(next.pres) - 
		getPoro(next.pres) * next.sat_oil * getB_oil_dp(next.pres, next.pres_bub, next.SATUR) / Boil ) / Boil;

	for(int i = 0; i < 4; i++)
	{
		upwind = upwindIsCur(cur, neighbor[i]);
		Variables& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;
		Boil_upwd = getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR);

		H += ht / cell.V * getTrans(cur, neighbor[i]) * 
			( getKr_oil(upwd.sat_oil) / props_oil.visc / Boil_upwd - 
			upwind * (next.pres - cells[ neighbor[i] ].u_next.pres) * getKr_oil(upwd.sat_oil) / props_oil.visc / Boil_upwd / Boil_upwd * getB_oil_dp(upwd.pres, upwd.pres_bub, upwd.SATUR) );
	}

	return H;
}

double TwoPhaseModel::solve_eq1_ds(int cur)
{
	double upwind;	
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	CylCell& cell = cells[cur];
	Variables& next = cell.u_next;
	double B_oil_upwd;
	double H = 0.0;

	H = getPoro(next.pres) / getB_oil(next.pres, next.pres_bub, next.SATUR);

	for(int i = 0; i < 4; i++)
	{
		upwind = upwindIsCur(cur, neighbor[i]);
		Variables& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;

		H += ht / cell.V * getTrans(cur, neighbor[i]) * 
			upwind * (next.pres - cells[ neighbor[i] ].u_next.pres) * getKr_oil_ds(upwd.sat_oil) / props_oil.visc / getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR);
	}

	return H;
}

double TwoPhaseModel::solve_eq1_dp_beta(int cur, int beta)
{
	CylCell& cell = cells[cur];

	double upwind = upwindIsCur(cur, beta);
	Variables& upwd = cells[ getUpwindIdx(cur, beta) ].u_next;
	double Boil_upwd = getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR);

	return -ht / cell.V * getTrans(cur, beta) * 
			( getKr_oil(upwd.sat_oil) / props_oil.visc / Boil_upwd + 
			(1.0 - upwind) * (cell.u_next.pres - cells[beta].u_next.pres) * getKr_oil(upwd.sat_oil) / props_oil.visc / Boil_upwd / Boil_upwd * getB_oil_dp(upwd.pres, upwd.pres_bub, upwd.SATUR) );
}

double TwoPhaseModel::solve_eq1_ds_beta(int cur, int beta)
{
	CylCell& cell = cells[cur];

	double upwind = upwindIsCur(cur, beta);
	Variables& upwd = cells[ getUpwindIdx(cur, beta) ].u_next;

	return ht / cell.V * getTrans(cur, beta) * (1.0 - upwind) * (cell.u_next.pres - cells[beta].u_next.pres) *
			getKr_oil_ds(upwd.sat_oil) / props_oil.visc / getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR);
}

double TwoPhaseModel::solve_eq2(int cur)
{
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	CylCell& cell = cells[cur];
	Variables& next = cell.u_next;
	Variables& prev = cell.u_prev;
	double H = 0.0;

	H = (getPoro(next.pres) * ( (1.0 - next.sat_oil) / getB_gas(next.pres) + next.sat_oil / getB_oil(next.pres, next.pres_bub, next.SATUR) ) -
				getPoro(prev.pres) * ( (1.0 - prev.sat_oil) / getB_gas(prev.pres) + prev.sat_oil / getB_oil(prev.pres, prev.pres_bub, prev.SATUR) ) );

	for(int i = 0; i < 4; i++)
	{
		Variables& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;

		H += ht / cell.V * getTrans(cur, neighbor[i]) * (next.pres - cells[ neighbor[i] ].u_next.pres) * 
			( getKr_oil(upwd.sat_oil) * getRs(upwd.pres, upwd.pres_bub, upwd.SATUR) / props_oil.visc / getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR) +
			getKr_gas(upwd.sat_oil) / props_gas.visc / getB_gas(upwd.pres) );
	}

	return H;
}

double TwoPhaseModel::solve_eq2_dp(int cur)
{
	double upwind;	
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	CylCell& cell = cells[cur];
	Variables& next = cell.u_next;
	double Boil_upwd, Bgas_upwd, rs_upwd;
	double Boil = getB_oil(next.pres, next.pres_bub, next.SATUR);
	double Bgas = getB_gas(next.pres);
	double rs = getRs(next.pres, next.pres_bub, next.SATUR);
	double H = 0.0;

	H = ( (next.sat_oil * rs / Boil + (1.0 - next.sat_oil) / Bgas) * getPoro_dp(next.pres) - 
		getPoro(next.pres) * ( (1.0 - next.sat_oil) / Bgas / Bgas * getB_gas_dp(next.pres) + 
		next.sat_oil * rs / Boil / Boil * getB_oil_dp(next.pres, next.pres_bub, next.SATUR) - 
		next.sat_oil / Boil * getRs_dp(next.pres, next.pres_bub, next.SATUR) ) );

	for(int i = 0; i < 4; i++)
	{
		upwind = upwindIsCur(cur, neighbor[i]);
		Variables& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;
		Boil_upwd = getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR);
		Bgas_upwd = getB_gas(upwd.pres);
		rs_upwd = getRs(upwd.pres, upwd.pres_bub, upwd.SATUR);

		H += ht / cell.V * getTrans(cur, neighbor[i]) * 
			( getKr_oil(upwd.sat_oil) * rs_upwd / props_oil.visc / Boil_upwd + getKr_gas(upwd.sat_oil) / props_gas.visc / Bgas_upwd + 
			upwind * (next.pres - cells[ neighbor[i] ].u_next.pres) * 
			( getKr_oil(upwd.sat_oil) / props_oil.visc / Boil_upwd * (getRs_dp(upwd.pres, upwd.pres_bub, upwd.SATUR) - getB_oil_dp(upwd.pres, upwd.pres_bub, upwd.SATUR) / Boil_upwd) - 
			getKr_gas(upwd.sat_oil) / props_gas.visc / Bgas_upwd / Bgas_upwd * getB_gas_dp(upwd.pres) ));
	}

	return H;
}

double TwoPhaseModel::solve_eq2_ds(int cur)
{
	double upwind;	
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	CylCell& cell = cells[cur];
	Variables& next = cell.u_next;
	double B_oil_upwd;
	double H = 0.0;

	H = getPoro(next.pres) * ( getRs(next.pres, next.pres_bub, next.SATUR) / getB_oil(next.pres, next.pres_bub, next.SATUR) - 1.0 / getB_gas(next.pres) );
	double tempp = getB_gas(next.pres) ;
	for(int i = 0; i < 4; i++)
	{
		upwind = upwindIsCur(cur, neighbor[i]);
		Variables& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;

		H += ht / cell.V * getTrans(cur, neighbor[i]) * upwind * (next.pres - cells[ neighbor[i] ].u_next.pres) * 
			( getRs(upwd.pres, upwd.pres_bub, upwd.SATUR) * getKr_oil_ds(upwd.sat_oil) / props_oil.visc / getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR) + 
			getKr_gas_ds(upwd.sat_oil) / props_gas.visc / getB_gas(upwd.pres) );
	}

	return H;
}

double TwoPhaseModel::solve_eq2_dp_beta(int cur, int beta)
{
	CylCell& cell = cells[cur];

	double upwind = upwindIsCur(cur, beta);
	Variables& upwd = cells[ getUpwindIdx(cur, beta) ].u_next;
	double Boil_upwd = getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR);
	double Bgas_upwd = getB_gas(upwd.pres);
	double rs_upwd = getRs(upwd.pres, upwd.pres_bub, upwd.SATUR);

	return -ht / cell.V * getTrans(cur, beta) * 
			( getKr_oil(upwd.sat_oil) * rs_upwd / props_oil.visc / Boil_upwd + getKr_gas(upwd.sat_oil) / props_gas.visc / Bgas_upwd + 
			(1.0 - upwind) * (cells[cur].u_next.pres - cells[beta].u_next.pres) * 
			( getKr_oil(upwd.sat_oil) / props_oil.visc / Boil_upwd * (getRs_dp(upwd.pres, upwd.pres_bub, upwd.SATUR) - getB_oil_dp(upwd.pres, upwd.pres_bub, upwd.SATUR) / Boil_upwd) - 
			getKr_gas(upwd.sat_oil) / props_gas.visc / Bgas_upwd / Bgas_upwd * getB_gas_dp(upwd.pres) ));
}

double TwoPhaseModel::solve_eq2_ds_beta(int cur, int beta)
{
	CylCell& cell = cells[cur];

	double upwind = upwindIsCur(cur, beta);
	Variables& upwd = cells[ getUpwindIdx(cur, beta) ].u_next;

	return ht / cell.V * getTrans(cur, beta) * (1.0 - upwind) * (cell.u_next.pres - cells[beta].u_next.pres) *
		( getRs(upwd.pres, upwd.pres_bub, upwd.SATUR) * getKr_oil_ds(upwd.sat_oil) / props_oil.visc / getB_oil(upwd.pres, upwd.pres_bub, upwd.SATUR) +
		getKr_gas_ds(upwd.sat_oil) / props_gas.visc / getB_gas(upwd.pres) );
}

double TwoPhaseModel::averValue(int var_idx)
{
	int idx;
	double val = 0.0;

	for(int i = 0; i < cellsNum_r; i++)
		for(int j = 0; j < cellsNum_z; j++)
		{
			idx = (i + 1) * (cellsNum_z + 2) + j + 1;
			val += cells[idx].V * cells[idx].u_next.values[var_idx];
		}

	return val / Volume;
}

double TwoPhaseModel::convergance(int key, int& ind, int& var_idx)
{
	int idx;
	double relErr = 0.0;
	double cur_relErr = 0.0;
	ind = 0;	var_idx = 0;

	if(key == PRES)
	{
		double var_next, var_iter;
		
		for(int i = 0; i < cellsNum_r; i++)
			for(int j = 0; j < cellsNum_z; j++)
			{
				idx = (i + 1) * (cellsNum_z + 2) + j + 1;
				CylCell& cell = cells[idx];

				// Pressure
				var_next = cell.u_next.pres;	var_iter = cell.u_iter.pres;
				if(fabs(var_next) > EQUALITY_TOLERANCE)
				{
					cur_relErr = fabs( (var_next - var_iter) / var_next );
					if(cur_relErr > relErr)
					{
						relErr = cur_relErr;
						ind  = idx;
						var_idx = 1;
					}
				}

				// Saturation
				var_next = cell.u_next.sat_oil;	var_iter = cell.u_iter.sat_oil;
				if(fabs(var_next) > EQUALITY_TOLERANCE)
				{
					cur_relErr = fabs( (var_next - var_iter) / var_next );
					if(cur_relErr > relErr)
					{
						relErr = cur_relErr;
						ind  = idx;
						var_idx = 2;
					}
				}
			}
	}
	return relErr;
}

void TwoPhaseModel::copyNewLayer()
{
	vector<CylCell>::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
	{
		it->u_prev = it->u_next;
		it->u_iter = it->u_next;
	}
}

void TwoPhaseModel::copyIterLayer()
{
	vector<CylCell>::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
	{
		it->u_iter = it->u_next;
	}
}

void TwoPhaseModel::snapshot(int i)
{
	snapshotter->snapshot(i);
}

void TwoPhaseModel::snapshot_all(int i)
{
	snapshotter->snapshot_all(i);
}
