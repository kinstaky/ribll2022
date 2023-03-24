#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "include/detectors.h"

using namespace ribll;


/// @brief print usage of the program
/// @param[in] name program name
///
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run detector(..)\n"
		<< "  run               Set run number.\n"
		<< "  detector(..)      Choose detector(s) name.\n"
		<< "Options:\n"
		<< "  -h                Print this help information.\n"
		<< "  -t  tag           Set trigger tag.\n"
		<< "  -e                Extract trigger.\n";
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @param[out] extract extract or not
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag,
	bool &extract
) {
	// initialize
	help = false;
	trigger_tag.clear();
	extract = false;
	// start index of positional arugments
	int result = 0;
	for (result = 1; result < argc; ++result) {
		// assumed that all options have read
		if (argv[result][0] != '-') break;
		// short option contains only one letter
		if (argv[result][2] != 0) return -result;
		if (argv[result][1] == 'h') {
			help = true;
			return result;
		} else if (argv[result][1] == 't') {
			// option of trigger tag
			// get tag in next argument
			++result;
			// miss arguemnt behind option
			if (result == argc) return -argc;
			trigger_tag = argv[result];
		} else if (argv[result][1] == 'e') {
			// option of extract tag
			extract = true;
		} else {
			return -result;
		}
	}
	return result;
}


const std::map<std::string, std::pair<double, double>> window_edge_map = {
	{"tof", {400.0, 800.0}},
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

	// print usage
	bool help = false;
	// trigger tag
	std::string trigger_tag;
	// trigger extract tag
	bool extract = false;
	// start index of postional paramters
	int pos_start = ParseArguments(argc, argv, help, trigger_tag, extract);

	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}

	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invalid option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option needs parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}

	if (pos_start+1 == argc) {
		// positional arguments less than 2
		std::cerr << "Error: Miss detector argument(s).\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// run number
	int run = atoi(argv[pos_start]);
	// list of detector names
	std::vector<std::string> detector_names;
	for (int i = pos_start+1; i < argc; ++i) {
		detector_names.push_back(std::string(argv[i]));
	}

	for (auto detector_name : detector_names) {
		// merge ADSSD trigger
		if (detector_name == "ta") {
			int result = MergeAdssdTrigger(trigger_tag, run);
			if (result) {
				std::cerr << "Error: Merge ADSSD trigger failed.\n";
			}
			continue;
		}

		std::shared_ptr<Detector> detector = CreateDetector(detector_name, run);
		if (!detector) continue;

		std::pair<double, double> window_edge = {-1000.0, 1000.0};
		// search for special window left and right edge
		auto search = window_edge_map.find(detector_name);
		// found special edge and assign to window_edge
		if (search != window_edge_map.end()) {
			window_edge = search->second;
		}

		if (!extract) {
			// extract tag empty, just match trigger
			int result = detector->MatchTrigger(
				trigger_tag,
				window_edge.first,
				window_edge.second
			);
			if (result) {
				// failed
				std::cerr << "Error: Match " << trigger_tag << " trigger for "
					<< detector_name << " failed.\n";
			}
		} else {
			int result = detector->ExtractTrigger(
				trigger_tag,
				window_edge.first,
				window_edge.second
			);
			if (result) {
				// failed
				std::cerr << "Error: Extract from " << trigger_tag
					<< " trigger with " << detector_name << " failed.\n";
			}
		}
	}
	return 0;
}