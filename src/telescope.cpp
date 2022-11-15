#include "include/telescope.h"

namespace ribll {

Telescope::Telescope(size_t layer)
: layers_(layer) {
}


Telescope::~Telescope() {
}


int Telescope::AttachDetector(size_t layer, Detector *detector) {
	if (layer >= layers_ || layer >= max_layers) {
		std::cerr << "Layer should be less than "
			<< (layers_ < max_layers ? layers_ : max_layers) << ".\n"; 
		return -1;
	}
	detectors_[layer].push_back(detector);
	return 0;
}


int Telescope::Correlate() {
	for (size_t i = 0; i < layers_; ++i) {
		for (Detector *detector : detectors_[i]) {
			detector->
		}
	}

}

