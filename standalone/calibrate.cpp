#include <iostream>
#include <vector>
#include <string>

#include "include/telescopes.h"
#include "include/detectors.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run name\n"
		"  run               Set run number.\n"
		"  name              Set telescope/detector name.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n"
		"  -a                Alpha calibrate.\n"
		"  -c                CsI calibrate.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @param[out] alpha_calibrate alpha calibration
/// @param[out] csi_calibrate csi calibration
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag,
	bool &alpha_calibrate,
	bool &csi_calibrate
) {
	// initialize
	help = false;
	trigger_tag.clear();
	alpha_calibrate = false;
	csi_calibrate = false;
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
		} else if (argv[result][1] == 'a') {
			alpha_calibrate = true;
		} else if (argv[result][1] == 'c') {
			csi_calibrate = true;
		} else {
			return -result;
		}
	}
	return result;
}

int main(int argc, char **argv) {
	if (argc < 2) {
		PrintUsage(argv[0]);
		return -1;
	}

	// help flag
	bool help = false;
	// trigger tag
	std::string tag;
	// alpha calibration
	bool alpha = false;
	// csi calibration
	bool csi = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag, alpha, csi);

	// need help
	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}

	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}

	if (pos_start + 1 >= argc) {
		// positional arguments less than 3
		std::cerr << "Error: Miss telescope argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}


	int run = atoi(argv[pos_start]);
	std::string telescope_name = argv[pos_start+1];

	std::shared_ptr<Telescope> telescope =
		CreateTelescope(telescope_name, run, tag);
	if (!telescope) {
		std::shared_ptr<Detector> detector =
			CreateDetector(telescope_name, run, tag);
		if (!detector) {
			std::cerr << "Error: Telescope/Detector "
				<< telescope_name << " not found.\n";
			return -1;
		}
		if (detector->Calibrate()) {
			std::cerr << "Error: Calibrate "
				<< telescope_name<< " failed.\n";
			return -1;
		}
		return 0;
	}

	if (alpha) {
		if (telescope->AlphaCalibrate()) {
			std::cerr << "Error: Alpha calibrate "
				<< telescope_name << " failed.\n";
			return -1;
		}
	} else if (csi) {
		if (telescope->CsiCalibrate()) {
			std::cerr << "Error: Calibrate CsI(Tl) in "
				<< telescope_name << " failed.\n";
			return -1;
		}
	} else {
		if (telescope->Calibrate()) {
			std::cerr << "Error: Calibrate "
				<< telescope_name << " failed.\n";
			return -1;
		}
	}
	return 0;
}