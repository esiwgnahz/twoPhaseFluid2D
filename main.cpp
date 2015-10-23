#include <new>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <iostream>

#include "TwoPhaseModel.h"
#include "TimeManager.h"

using namespace std;

void setDataFromFile(vector< pair<double,double> >& vec, string fileName)
{
	ifstream file;
	file.open(fileName.c_str(), ifstream::in);
	
	double temp1, temp2;
	while( !file.eof() )
	{
		file >> temp1;
		if( file.eof() )
			break;
		file >> temp2;
		vec.push_back(make_pair(temp1, temp2));
	}

	file.close();
}

Properties* getProps()
{
	Properties* props = new Properties();

	props->cellsNum_r = 500;
	props->cellsNum_z = 10;

	props->timePeriods.push_back(1296000.0);
	props->timePeriods.push_back(2000000.0);
	
	props->rates.push_back(10.0);
	props->rates.push_back(0.0);
 
	props->skins.push_back(0.0); 
	props->skins.push_back(0.0);

	props->radius.push_back(0.1524);
	props->radius.push_back(0.1524);

	props->ht = 3600.0;
	props->ht_min = 100.0;
	props->ht_max  = 100000.0;
	
	props->alpha = 7200.0;

	props->r_w = 0.1524;
	props->r_e = 3000.0;

	props->p_init = 160.0 * 1.E5;
	props->p_sat = 160.0 * 1.E5;
	props->T_init = 121.111 + 273.15;
	props->s_init = 0.999;
	props->height = 10.0;
	props->h1 = 1500.0;
	props->h2 = 1510.0;
	props->depth_point = 1500.0;
	props->perm_r = 50.0;
	props->perm_z = 0.0;
	props->m = 0.01;
	props->visc_gas = 0.02833;
	props->visc_oil = 0.25137;
	props->dens_oil_stc = 800.026;
	props->dens_gas_stc = 0.72275;
	props->dens_sk_stc = 2000.0;
	props->beta_oil = 4.0*1.E-9;
	props->beta_sk = 4.0*1.E-9;

	// Thermal properties
	props->jt_oil = 4.5*1.e-6;
	props->jt_gas = -1.4*1.e-5;
	props->ad_oil = 2.e-6;
	props->ad_gas = 2.e-6;
	props->c_oil = 1600.0;
	props->c_gas = 2000.0;
	props->c_sk = 800.0;
	props->kappa_eff = 5.5*1e-6;
	props->L = -50.0*1.e3;

	props->lambda_sk = 7.0;
	props->lambda_oil = 0.4;
	props->lambda_gas = 0.05;
	
	// Defining relative permeabilities
	setDataFromFile(props->kr_oil, "koil.txt");
	setDataFromFile(props->kr_gas, "kgas.txt");
	
	// Defining volume factors
	//props->byDefault.B_oil = true;
	setDataFromFile(props->B_oil, "Boil.txt");
	//props->byDefault.B_gas = false;
	setDataFromFile(props->B_gas, "Bgas.txt");

	//props->byDefault.Rs = true;
	setDataFromFile(props->Rs, "Rs.txt");

	return props;
}

int main(int argc, char* argv[])
{
	Properties* props = getProps();
	TwoPhaseModel model;
	model.load(*props);

	TimeManager manager(&model);
	manager.start();

	return 0;
}