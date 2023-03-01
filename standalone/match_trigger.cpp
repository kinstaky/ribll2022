#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "include/detectors.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " run detector(..)\n"
		<< "  run                       run number\n"
		<< "  detector(..)              detector(s) name\n";
}

const std::map<std::string, std::pair<double, double>> window_edge_map = {
	{"tafcsi", {-1300.0, 800.0}},
	{"tabcsi", {-1000.0, 800.0}},
	{"t0csi", {-1600.0, 800.0}},
	{"t1csi", {-1000.0, 800.0}}
};

int main(int argc, char **argv) {
	// check parameters
	if (argc < 3) {
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[1]);
	std::vector<std::string> detector_names;
	for (int i = 2; i < argc; ++i) {
		detector_names.push_back(std::string(argv[i]));
	}

	for (auto detector_name : detector_names) {
		std::shared_ptr<Detector> detector = CreateDetector(detector_name, run);
		if (!detector) continue;

		std::pair<double, double> window_edge = {-1000.0, 1000.0};
		// search for special window left and right edge
		auto search = window_edge_map.find(detector_name);
		// found special edge and assign to window_edge
		if (search != window_edge_map.end()) {
			window_edge = search->second;
		}
		if (detector->MatchTrigger(window_edge.first, window_edge.second)) {
			std::cerr << "Error: Match trigger for "
				<< detector_name << " failed.\n";
			continue;
		}
	}
	return 0;
}