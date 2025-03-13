#include <TFile.h>
#include <TString.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1F.h>

#include "include/defs.h"

using namespace ribll;

constexpr int csi_index = 0;

double CsiEnergy(int index, double channel) {
	double a0 = power_csi_param[index][0];
	double a1 = power_csi_param[index][1];
	double a2 = power_csi_param[index][2];
	return pow((channel - a2) / a0, 1.0 / a1);
}


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] taf_index\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -s                Use simulated data.\n";
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] sim use simulated data
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	bool &sim
) {
	// initialize
	help = false;
	sim = false;
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
		} else if (argv[result][1] == 's') {
			sim = true;
		} else {
			return -result;
		}
	}
	return result;
}


int main(int argc, char **argv) {
	if (argc > 3) {
		PrintUsage(argv[0]);
		return -1;
	}
	// help flag
	bool help = false;
	// simulated data
	bool sim = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, sim);
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
	if (pos_start >= argc) {
		// positional arguments less than 1
		std::cerr << "Error: Parameter case not found.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// TAF index
	int taf_index = atoi(argv[pos_start]);


	TString file_name;
	if (sim) {
		file_name.Form(
			"%s%staf%dcsi-calibrate-sim-0002.root",
			kGenerateDataPath,
			kShowDir,
			taf_index
		);
	} else {
		file_name.Form(
			"%s%staf%dcsi-calibrate.root",
			kGenerateDataPath,
			kShowDir,
			taf_index
		);
	}
	TFile ipf(file_name, "read");
	TGraph *gdpa, *gdpb;
	gdpa = (TGraph*)ipf.Get("gdpa");
	gdpb = (TGraph*)ipf.Get("gdpb");

	constexpr double range[][2] = {
		{5.0, 15.0},
		{5.0, 25.0},
		{10.0, 30.0},
		{15.0, 35.0},
		{20.0, 40.0},
		{20.0, 50.0},
		{25.0, 55.0},
		{25.0, 65.0}
	};

	TString output_file_name = TString::Format(
		"%s%staf%dcsi-error-range%s.root",
		kGenerateDataPath,
		kShowDir,
		taf_index,
		sim ? "-sim" : ""
	);
	TFile opf(output_file_name, "recreate");
	std::vector<TH1F> hist_range_a, hist_range_b;
	for (int i = 0; i < 8; ++i) {
		hist_range_a.emplace_back(
            TString::Format("ha%d", 10+i*5),
			TString::Format("correlated range of %d to %d MeV", 9+i*5, 11+i*5),
			100, range[i][0], range[i][1]
		);
		hist_range_b.emplace_back(
            TString::Format("hb%d", 10+i*5),
			TString::Format("correlated range of %d to %d MeV", 9+i*5, 11+i*5),
			100, range[i][0], range[i][1]
		);
	}

	int n = gdpa->GetN();
	double *x = gdpa->GetX();
	double *y = gdpa->GetY();
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < 8; ++j) {
			double mid = 10.0 + 5.0 * j;
			if (x[i] > mid - 1.0 && x[i] < mid + 1.0) {
				hist_range_a[j].Fill(CsiEnergy(taf_index*2, y[i]));
				break;
			}
		}
	}

	constexpr double fit_range[][2] = {
		{7.0, 13.0},
		{10.0, 20.0},
		{15.0, 25.0},
		{20.0, 30.0},
		{25.0, 35.0},
		{30.0, 40.0},
		{30.0, 50.0},
		{35.0, 55.0}
	};

	// fit
	double fit_result_a[8][2];
	for (int i = 0; i < 8; ++i) {
		TF1 *f1 = new TF1(
			TString::Format("fa%d", i), "gaus",
			fit_range[i][0], fit_range[i][1]
		);
		hist_range_a[i].Fit(f1, "RQ+");
		fit_result_a[i][0] = f1->GetParameter(1);
		fit_result_a[i][1] = f1->GetParameter(2);
	}
	for (size_t i = 0; i < 8; ++i) {
		std::cout << fit_result_a[i][0] << ", " << fit_result_a[i][1] << "\n";
	}



	n = gdpb->GetN();
	x = gdpb->GetX();
	y = gdpb->GetY();
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < 8; ++j) {
			double mid = 10.0 + 5.0 * j;
			if (x[i] > mid - 1.0 && x[i] < mid + 1.0) {
				hist_range_b[j].Fill(CsiEnergy(taf_index*2+1, y[i]));
				break;
			}
		}
	}

	// fit
	double fit_result_b[8][2];
	for (int i = 0; i < 8; ++i) {
		TF1 *f1 = new TF1(
			TString::Format("fb%d", i), "gaus",
			fit_range[i][0], fit_range[i][1]
		);
		hist_range_b[i].Fit(f1, "RQ+");
		fit_result_b[i][0] = f1->GetParameter(1);
		fit_result_b[i][1] = f1->GetParameter(2);
	}
	for (size_t i = 0; i < 8; ++i) {
		std::cout << fit_result_b[i][0] << ", " << fit_result_b[i][1] << "\n";
	}

	opf.cd();
	for (TH1F &hist : hist_range_a) hist.Write();
	for (TH1F &hist : hist_range_b) hist.Write();
	opf.Close();
	ipf.Close();

}