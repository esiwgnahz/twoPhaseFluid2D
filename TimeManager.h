#ifndef TIMEMANAGER_H_
#define TIMEMANAGER_H_

#include "TwoPhaseModel.h"
#include <iostream>

using namespace std;

class TwoPhaseModel;

class TimeManager
{
protected:
	TwoPhaseModel* model;

	// Time dimension
	double t_dim;
	// Current time & time since well is stopped
	double cur_t, cur_t_log;
	// Index of current time period
	int curTimePeriod;
	// Rate, pressure and temperature before well test
	double q0;
	double p0;
	double t0;

	bool first;

	// Cells indexes
	int idx1, idx2;

	// Output file streams
	ofstream plot_Sdyn;
	ofstream plot_Pdyn;
	ofstream plot_Tdyn;
	ofstream plot_Psat;
	ofstream plot_Rs;

public:
	TimeManager(TwoPhaseModel* _model) : model(_model)
	{
		first = true;
		cur_t = cur_t_log = 0.0;
		curTimePeriod = 0;

		idx1 = int(model->cellsNum_z / 2) + 1;
		idx2 = idx1 + model->cellsNum_z + 2;

		t_dim = model->t_dim;

		plot_Tdyn.open("model/T_dyn.dat", ofstream::out);
		plot_Pdyn.open("model/P_dyn.dat", ofstream::out);
		plot_Sdyn.open("model/S_dyn.dat", ofstream::out);
	};

	~TimeManager()
	{
		plot_Tdyn.close();
		plot_Pdyn.close();
		plot_Sdyn.close();
	};

	// Time step & rate control
	inline void control(int iteration)
	{
		writeData();

		if(cur_t >= model->period[curTimePeriod])
		{
			curTimePeriod++;
			model->ht = model->ht_min;
			model->setPeriod(curTimePeriod);
		}

		if(model->ht <= model->ht_max && iteration < 4)
			model->ht = model->ht * 1.5;
		else if(iteration > 6 && model->ht > model->ht_min)
			model->ht = model->ht / 1.5;

		if(cur_t + model->ht > model->period[curTimePeriod])
			model->ht = model->period[curTimePeriod] - cur_t;

		cur_t += model->ht;
	};

	// Rate after stopping
	/*inline void rateAfterStop()
	{
		if(cur_t > model->getTimePeriod(model->getNumberOfPeriods() - 2))
		{
			if(first) {
				q0 = model->getQ();
				p0 = model->getPress(0);
				t0 = model->getTemperature(0);
				model->SetbKvdPart(true);
				
				plot_Plog.open("model/P_log.dat", ofstream::out);
				plot_Tlog.open("model/T_log.dat", ofstream::out);
				plot_Slog.open("model/S_log.dat", ofstream::out);
	
				first = false;
			}

			//model->setQ(q0 * exp( -cur_t_log/model->getAlpha() ));
			//std::cout << "q0=" << q0 << "\tq=" << q0 * exp( -cur_t_log/model->getAlpha() ) << "\tcur_t_log=" << cur_t_log << std::endl;
			//model->writeData();

			writeLogData();
	
			cur_t_log += model->getHt();
			model->setCurTime(cur_t_log);
		}
	};*/

	inline double getCurTime() const
	{
		return cur_t;
	};

	inline void writeData()
	{
		//plot_Tdyn << cur_t * SEC_TO_HOUR << "\t" << model->getTemperature(0) << endl;
		plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << (model->cells[idx1].u_next.pres + model->cells[idx2].u_next.pres) / 2.0 << endl;
		plot_Sdyn << cur_t * t_dim / 3600.0 << "\t" << (model->cells[idx1].u_next.sat_oil + model->cells[idx2].u_next.sat_oil) / 2.0 << endl; 
		//plot_Psat << cur_t * SEC_TO_HOUR << "\t" << model->getPress_sat(0) * PA_TO_BAR << endl;
		//plot_Rs << cur_t * SEC_TO_HOUR << "\t" << model->Solve_Rs(model->getPress(0) * PA_TO_BAR, model->getPress_sat(0) * PA_TO_BAR, 0) << endl;
	};

	/*inline void writeLogData()
	{
		plot_Plog << cur_t_log * SEC_TO_HOUR << " " << (model->getPress(0) - p0 ) * PA_TO_BAR << endl;
		plot_Tlog << cur_t_log * SEC_TO_HOUR << " " << model->getTemperature(0) - t0 << endl;
		plot_Slog << cur_t_log * SEC_TO_HOUR << " " << model->getSoil(0) << endl;
	};*/

	inline void start()
	{
		int counter = 0;
		int newton_itr = 1;
		
		model->setPeriod(curTimePeriod);
		//double Tt = model->period[model->period.size()-1];
		
		model->snapshot_all(counter++);
		
		//while(cur_t <= Tt)
		//{
			control(newton_itr);
			newton_itr = model->doNextStep();
			model->snapshot_all(counter++);
		//}

		
	};
};

#endif /* TIMEMANAGER_H_ */