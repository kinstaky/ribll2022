#include <iostream>
#include <string>
#include <map>

#include "include/detectors.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " run detector\n"
		<< "  run                       run number\n"
		<< "  detector                  detector name\n";
}

const std::map<std::string, std::pair<double, double>> window_edge_map = {
	{"tafcsi", {-1300.0, 800.0}},
	{"tabcsi", {-1000.0, 800.0}},
	{"t0csi", {-1600.0, 800.0}},
	{"t1csi", {-1000.0, 800.0}}
};

int main(int argc, char **argv) {
	// check parameters
	if (argc != 3) {
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[1]);
	std::string detector_name = std::string(argv[2]);

	std::shared_ptr<Detector> detector = CreateDetector(detector_name, run);
	if (!detector) return -1;

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
		return -1;
	}
	return 0;
}