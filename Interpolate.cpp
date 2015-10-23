#include "Interpolate.h"
#include <iostream>

Interpolate::Interpolate()
{
}

Interpolate::~Interpolate()
{
}

Interpolate::Interpolate(double *ptx, double *pty, int _N)
{
	Nf = 1000000;
	Flag = new int[Nf+1];
	N  = _N;
	x = new double[N];
	y = new double[N];

	for (int i= 0; i < N; i++){
		x[i] = ptx[i];
		y[i] = pty[i];
	}
	
	xmin = x[0];
	xmax = x[N-1];
	int j;
	temp = (xmax - xmin)/Nf;
	
	for (int i = 0; i < N-1; i++)
		for (j = floor((x[i] - xmin)/temp); j <= floor((x[i+1] - xmin)/temp); j++)
			Flag[j] = i;
	
	IsSetDiff = false;

}

Interpolate::Interpolate(double *ptx, double *pty, double *dpty, int N)
{
	Nf = 10000;
	Flag = new int[Nf+1];
	x = new double[N];
	y = new double[N];
	dy = new double[N];
	for (int i= 0; i < N; i++){
		x[i] = ptx[i];
		y[i] = pty[i];
		dy[i] = dpty[i];
	}
	xmin = x[0];
	xmax = x[N-1];
	temp = (xmax - xmin)/Nf;
	for (int i = 0; i < N-1; i++){
		for (int j = floor((x[i] - xmin)/temp); j <= floor((x[i+1] - xmin)/temp); j++)
		Flag[j] = i;
	}
	IsSetDiff = true;

}

Interpolate::Interpolate(double *ptx, double *pty, double *dpty, double *d2pty, int N)
{
	Nf = 10000;
	Flag = new int[Nf+1];
	x = new double[N];
	y = new double[N];
	dy = new double[N];
	d2y = new double[N];

	for (int i= 0; i < N; i++){
		x[i] = ptx[i];
		y[i] = pty[i];
		dy[i] = dpty[i];
		d2y[i]  = d2pty[i];
	}
	xmin = x[0];
	xmax = x[N-1];
	temp = (xmax - xmin)/Nf;
	for (int i = 0; i < N-1; i++){
		for (int j = floor((x[i] - xmin)/temp); j <= floor((x[i+1] - xmin)/temp); j++)
		Flag[j] = i;
	}
	IsSetDiff = true;
	IsSetDiff2 = true;
}

double Interpolate::Solve(double arg)
{
	if (arg < xmin ) 
		arg = xmin;
	else if(arg > xmax) 
		arg = xmax;

	int index = floor( (arg - xmin) / temp );
	double dx = x[Flag[index]+1] - x[Flag[index]];
	double dy = y[Flag[index]+1] - y[Flag[index]];
	
	return (arg - x[Flag[index]])/dx*dy + y[Flag[index]];
}

double Interpolate::DSolve(double arg)
{
	double result;
	double dx;
	double dcy;
	
	if(arg < xmin)
		return 0;
	else if(arg > xmax)
		return 0;

	int i = floor((arg - xmin)/temp);
	int index = Flag[i];
	
	if (IsSetDiff)
	{
		 dx = x[index+1] - x[index];
		 dcy = dy[index+1] - dy[index];
		 result = (arg - x[index])/dx*dcy + dy[index];
	} else {
		dx = (x[index] - x[index + 1]);
		dcy = (y[index] - y[index + 1]);
		result = dcy/dx;
	}

	return result;
}

double Interpolate::D2Solve(double arg)
{
	double result;
	double dx;
	double dcy;
	
	int i = floor((arg - xmin)/temp);
	if (i > Nf)
		i= Nf;
	if (i < 0)
		i = 0;
    int index = Flag[i];
	
	if (IsSetDiff2)
	{
		 dx = x[index+1] - x[index];
		 dcy = d2y[index+1] - d2y[index];
		 result = (arg - x[index])/dx*dcy + d2y[index];
	} else 
		result = 0;
	
	return result;

}