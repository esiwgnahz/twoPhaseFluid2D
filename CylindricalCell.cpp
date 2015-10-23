#include "CylindricalCell.h"

CylCell::CylCell()
{
	number = -1;
	isGhost = true;
	V = hr = hz = 0.0;
}

CylCell::CylCell(int _number, double _r, double _z) : number(_number), r(_r), z(_z)
{
	isGhost = true;
	V = hr = hz = 0.0;
}

CylCell::CylCell(int _number, double _r, double _z, double _hr, double _hz) : number(_number), r(_r), z(_z), hr(_hr), hz(_hz)
{
	isGhost = true;
	V = 2 * M_PI * r * hr * hz;
}

CylCell::CylCell(int _number, double _r, double _z, double _hr, double _hz, bool _isGhost) : number(_number), r(_r), z(_z), hr(_hr), hz(_hz), isGhost(_isGhost)
{
	V = 2 * M_PI * r * hr * hz;
}

CylCell::~CylCell()
{
}