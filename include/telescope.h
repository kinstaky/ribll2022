#ifndef __TELESCOPE_H__
#define __TELESCOPE_H__

#include <vector>

#include "include/detector/detector.h"

namespace ribll {

const size_t max_layers = 8;

class Telescope {
public:

	Telescope(size_t layer);

	~Telescope();

	int AttachDetector(size_t layer, Detector *detector);

	int Correlate();

private:
	// layer number
	size_t layers_;
	// list of detectors
	std::vector<Detector*> detectors_[max_layers];
};

};

#endif