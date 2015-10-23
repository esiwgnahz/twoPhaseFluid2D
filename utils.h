#ifndef UTILS_H_
#define UTILS_H_

#include <sstream>
#include <string>

inline bool IsNan(double a)
{
	if (a!=a)  return true;
	return false;
};

inline double sign(int a)
{
	if (a > 0) return 1;
	else if (a < 0) return -1;
	else return 0;
};

inline double MilliDarcyToM2(double perm)
{
	return perm * 0.986923 * 1.E-15;
};

inline double M2toMilliDarcy(double perm)
{
	return perm * 1.E15 / 0.986923;
};

inline double cPToPaSec(double visc)
{
 return visc/1000.0;
};

template <typename T>
inline std::string to_string(T value)
{
	std::ostringstream os ;
	os << value;
	return os.str();
};

struct sort_pair_first {
    bool operator() (const std::pair<double,double> &left, const std::pair<double,double> &right) 
	{
        return left.first < right.first;
    }
};

struct sort_pair_second {
    bool operator() (const std::pair<double,double> &left, const std::pair<double,double> &right) 
	{
        return left.second < right.second;
    }
};

#endif /* UTILS_H_ */