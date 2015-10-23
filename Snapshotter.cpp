#include "Snapshotter.h"
#include "utils.h"

using namespace std;

Snapshotter::Snapshotter()
{
	pattern = "model/snap_%{STEP}.grdecl";
}

Snapshotter::~Snapshotter()
{
}

void Snapshotter::setModel(TwoPhaseModel* _model)
{
	model = _model;
}

string Snapshotter::replace(string filename, string from, string to)
{
	size_t start_pos = 0;
    while((start_pos = filename.find(from, start_pos)) != string::npos) 
	{
		filename.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
	return filename;
}

string Snapshotter::getFileName(int i)
{
	string filename = pattern;
	return replace(filename, "%{STEP}", to_string(i));
}

void Snapshotter::snapshot(int i)
{
	string gridName = getFileName(i);

	ofstream grid;
	grid.open(gridName.c_str(), ofstream::out);

	int n = model->cellsNum_r * model->cellsNum_z;

	grid << "SPECGRID" << endl;
	grid << "  " << model->cellsNum_r << "  " << model->cellsNum_z << "  " << 1;
	grid << "  " << 1 << "  " << "F /" << endl;

	grid << "COORDSYS" << endl;
	grid << "  " << "1  34  'INCOMP  '  /" << endl;

	grid << "COORD" << endl;
	int idx = 0;
	double R_dim = model->R_dim;
	double z;
	int counter = 0;
	for(int k = 0; k < model->cellsNum_z+1; k++)
	{
		idx = k+1;
		z = R_dim * (model->cells[idx].z - model->cells[idx].hz / 2.0);
		for(int j = 0; j < model->cellsNum_r+1; j++)
		{
			idx = (j+1)*(model->cellsNum_z+2) + (k+1);
			grid << R_dim * (model->cells[idx].r - model->cells[idx].hr / 2.0) << "  " \
					<< z << " " \
					<< 1.0;
			counter++;
			if(counter % 2 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	if(counter % 2 == 0)
		grid << endl;
	else
		grid << endl << endl;

	grid << "ACTNUM" << endl;
	counter = 0;
	while(counter < n)
	{
		grid << counter << "  ";
		if(counter % 6 == 0)
			grid << endl;
		counter++;
	}

	grid << endl;

	// Pressure
	grid << "PRES" << endl;
	counter = 0;
	for(int k = 0; k < model->cellsNum_z; k++)
	{
		for(int j = 0; j < model->cellsNum_r; j++)
		{
			idx = (j+1)*(model->cellsNum_z+2) + (k+1);
			grid << model->cells[idx].u_next.pres;
			
			counter++;
			if(counter % 2 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid << endl;

	// Saturation
	grid << "SATUR" << endl;
	counter = 0;
	for(int k = 0; k < model->cellsNum_z; k++)
	{
		for(int j = 0; j < model->cellsNum_r; j++)
		{
			idx = (j+1)*(model->cellsNum_z+2) + (k+1);
			grid << model->cells[idx].u_next.sat_oil;
			
			counter++;
			if(counter % 2 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid.close();
}

void Snapshotter::snapshot_all(int i)
{
	string gridName = getFileName(i);

	ofstream grid;
	grid.open(gridName.c_str(), ofstream::out);

	int n = (model->cellsNum_r + 2) * (model->cellsNum_z + 2);

	grid << "SPECGRID" << endl;
	grid << "  " << model->cellsNum_r+2 << "  " << model->cellsNum_z+2 << "  " << 1;
	grid << "  " << 1 << "  " << "F /" << endl;

	grid << "COORDSYS" << endl;
	grid << "  " << "1  34  'INCOMP  '  /" << endl;

	grid << "COORD" << endl;
	int idx = 0;
	double R_dim = model->R_dim;
	double z;
	int counter = 0;

	idx = 0;
	z = R_dim * (model->cells[idx].z - model->cells[idx].hz / 2.0);
	grid << R_dim * (model->cells[idx].r - model->cells[idx].hr / 2.0) << "  " \
			<< z << " " \
			<< 1.0 << " ";
	counter++;
	for(int k = 0; k < model->cellsNum_r+2; k++)
	{
		idx = k * (model->cellsNum_z+2);
		grid << R_dim * (model->cells[idx].r + model->cells[idx].hr / 2.0) << "  " \
			<< z << " " \
			<< 1.0;
		counter++;

		if(counter % 2 == 0)
			grid << endl;
		else
			grid << "  ";
	}

	for(int k = 0; k < model->cellsNum_z+2; k++)
	{
		idx = k;
		z = R_dim * (model->cells[idx].z + model->cells[idx].hz / 2.0);
		grid << R_dim * (model->cells[idx].r - model->cells[idx].hr / 2.0) << "  " \
			<< z << " " \
			<< 1.0;
		counter++;
		
		if(counter % 2 == 0)
			grid << endl;
		else
			grid << "  ";
		
		for(int j = 0; j < model->cellsNum_r+2; j++)
		{
			idx = j*(model->cellsNum_z+2) + k;
			grid << R_dim * (model->cells[idx].r + model->cells[idx].hr / 2.0) << "  " \
					<< z << " " \
					<< 1.0;
			counter++;

			if(counter % 2 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid << endl;

	grid << "ACTNUM" << endl;
	counter = 0;
	while(counter < n)
	{
		grid << counter << "  ";
		if(counter % 6 == 0)
			grid << endl;
		counter++;
	}

	grid << endl;

	// Pressure
	grid << "PRES" << endl;
	counter = 0;
	for(int k = 0; k < model->cellsNum_z+2; k++)
	{
		for(int j = 0; j < model->cellsNum_r+2; j++)
		{
			idx = j * (model->cellsNum_z+2) + k;
			grid << model->cells[idx].u_next.pres;
			
			counter++;
			if(counter % 6 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid << endl;

	// Saturation
	grid << "SATUR" << endl;
	counter = 0;
	for(int k = 0; k < model->cellsNum_z+2; k++)
	{
		for(int j = 0; j < model->cellsNum_r+2; j++)
		{
			idx = j * (model->cellsNum_z+2) + k;
			grid << model->cells[idx].u_next.sat_oil;
			
			counter++;
			if(counter % 6 == 0)
				grid << endl;
			else
				grid << "  ";
		}
	}

	grid.close();
}