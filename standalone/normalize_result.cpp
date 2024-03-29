#include "include/detectors.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run detector[...]\n"
		"  run               Set run number.\n"
		"  detector          Set detector name.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n"
		"  -i num            Set iteration mode.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @param[out] iteartion iteration mode
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag,
	int &iteration
) {
	// initialize
	help = false;
	trigger_tag.clear();
	iteration = 0;
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
		} else if (argv[result][1] == 'i') {
			// option of iteration flag
			// get number in next argument
			++result;
			// miss arguemnt behind option
			if (result == argc) return -argc;
			iteration = atoi(argv[result]);
		} else {
			return -result;
		}
	}
	return result;
}

int main(int argc, char **argv) {
	if (argc < 3) {
		PrintUsage(argv[0]);
		return -1;
	}

	// help flag
	bool help = false;
	// trigger tag
	std::string tag;
	// iteration flag
	int iteration = 0;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag, iteration);

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
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}


	int run = atoi(argv[pos_start]);
	// list of detector names
	std::vector<std::string> dssd_names;
	for (int i = pos_start+1; i < argc; ++i) {
		dssd_names.push_back(std::string(argv[i]));
	}

	for (auto dssd_name : dssd_names) {
		std::shared_ptr<Dssd> dssd = CreateDssd(dssd_name, run, tag);
		if (!dssd) continue;

		if (dssd->NormalizeResult(iteration)) {
			std::cerr << "Error: Normalize result from "
				<< dssd_name << " failed.\n";
			continue;
		}
	}

	return 0;
}