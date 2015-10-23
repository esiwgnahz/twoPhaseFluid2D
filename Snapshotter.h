#ifndef SNAPSHOTTER_H_
#define SNAPSHOTTER_H_

#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

#include "CylindricalCell.h"
#include "TwoPhaseModel.h"

class TwoPhaseModel;

class Snapshotter {
protected:
	TwoPhaseModel* model;
	std::string pattern;

	std::string replace(std::string filename, std::string from, std::string to);
	std::string getFileName(int i);
public:
	Snapshotter();
	~Snapshotter();

	void setModel(TwoPhaseModel* _model);
	void snapshot(int i);
	void snapshot_all(int i);

};

#endif /* SNAPSHOTTER_H_ */