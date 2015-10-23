#ifndef CYLINDRICALCELL_H_
#define CYLINDRICALCELL_H_

#define VAR_NUM 4

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

struct Variables {
	union {
		double values[VAR_NUM];
		
		struct {
			double temp;
			double pres;
			double sat_oil;
			double pres_bub;
		};
	};
	bool SATUR;

	Variables& operator=(const Variables& a)
	{
		for(int i = 0; i < VAR_NUM; i++)
			values[i] = a.values[i];
	
		SATUR = a.SATUR;

		return *this;
	}
};

class CylCell {
public:
	int number;
	union {
		double cm [2];
		struct {
			double r;
			double z;
		};
	};
	
	double hr, hz;
	double V;

	Variables u_next;
	Variables u_iter;
	Variables u_prev;

	bool isGhost;

	CylCell();
	CylCell(int _number, double _r, double _z);
	CylCell(int _number, double _r, double _z, double _hr, double _hz);
	CylCell(int _number, double _r, double _z, double _hr, double _hz, bool _isGhost);
	~CylCell();
};

namespace std {
	inline std::ostream& operator<<(std::ostream& os, const CylCell& cell)
	{
		os << "\nCell number:\t" << cell.number << std::endl;
		os << "Mass center:\t" << cell.r << "\t" << cell.z << std::endl;
		os << "Values:\t";
		for(int i = 0; i < VAR_NUM; i++)
			os << cell.u_next.values[i] << "\t";
		os << std::endl;

		return os;
	}
}

#endif /* CYLINDRICALCELL_H_ */