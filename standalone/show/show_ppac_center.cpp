#include <iostream>
#include <iomanip>
#include <fstream>

#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TRandom3.h>

#include "include/event/particle_event.h"
#include "include/statistics/center_statistics.h"

using namespace ribll;

constexpr double d1z = 100.0;
constexpr double ppac_xz[3] = {-695.2, -454.2, -275.2};
constexpr double ppac_yz[3] = {-689.2, -448.2, -269.2};

double SimpleFit(const double *x, const double *y, double &k, double &b) {
	int n = 3;
	double sumx = 0.0;
	double sumy = 0.0;
	double sumxy = 0.0;
	double sumx2 = 0.0;
	for (int i = 0; i < n; ++i) {
		sumx += x[i];
		sumy += y[i];
		sumxy += x[i] * y[i];
		sumx2 += x[i] * x[i];
	}
	k = (sumxy - sumx*sumy/double(n)) / (sumx2 - sumx*sumx/double(n));
	b = (sumy - k*sumx) / double(n);
	double chi2 = 0.0;
	for (int i = 0; i < n; ++i) {
		double t = y[i] - k*x[i] - b;
		chi2 += t * t;
	}
	return chi2;
}


void FitAndFill(
	const double *x,
	const double *y,
	const double x0,
	const double y0,
	TH1F *h
) {
	// // fix T0D1 and fit 3 PPAC points
	// double a[3];
	// for (size_t i = 0; i < 3; ++i) {
	// 	a[i] = x[i] - x0;
	// }
	// double b[3];
	// for (size_t i = 0; i < 3; ++i) {
	// 	b[i] = y0 - y[i];
	// }
	// double k = -(a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
	// 	/ (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	// for (size_t i = 0; i < 3; ++i) {
	// 	h[i].Fill(a[i]*k + b[i]);
	// }

	// fit 3 PPAC points and calculate T0D1 offset
	// double k, b;
	// SimpleFit(x, y, k, b);
	// h[0].Fill(y[0] - b - k * x[0]);
	// h[1].Fill(y[1] - b - k * x[1]);
	// h[2].Fill(y[2] - b - k * x[2]);


	// // fit 3 PPAC points and calculate T0D1 offset
	// double k, b;
	// SimpleFit(x, y, k, b);
	// h[0].Fill(b+x0*k-y0);

	// fix first PPAC and T0D1 points and calculate offset of other two PPACs
	// hjx's method
	double k = (y0 - y[0]) / (x0 - x[0]);
	double b = y0 - k * x0;
	h[0].Fill(y[1] - b - k * x[1]);
	h[1].Fill(y[2] - b - k * x[2]);
}


double GausPol0(double *x, double *par) {
	return par[3] + par[0]*exp((x[0]-par[1])*(x[0]-par[1])/2.0/par[2]/par[2]);
}



void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run [end_run]\n"
		"  run               Set run number.\n"
		"  end_run           Set the last run, inclusive.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n"
		"Examples:\n"
		"  " << name << " 600        Show center of run 600.\n"
		"  " << name << " 600 700    Show center from run 600 to 700.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag
) {
	// initialize
	help = false;
	trigger_tag.clear();
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
		} else {
			return -result;
		}
	}
	return result;
}


/// @brief search for offset in one run
/// @param[in] run run number
/// @param[in] tag trigger tag
/// @returns 0 if success, -1 otherwise
///
int CalculateOffset(unsigned int run, const std::string &tag) {
	// t0 file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-particle-%s%04u.root",
		kGenerateDataPath,
		kParticleDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// t0 file
	TFile t0_file(t0_file_name, "read");
	// t0 tree
	TTree *tree = (TTree*)t0_file.Get("tree");
	if (!tree) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// xppac file name
	TString xppac_file_name;
	xppac_file_name.Form(
		"%s%sxppac-particle-%s%04u.root",
		kGenerateDataPath,
		kParticleDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	tree->AddFriend("xppac=tree", xppac_file_name);
	// input t0 event
	ParticleEvent t0_event;
	// input particle type event
	ParticleEvent xppac_event;
	unsigned short xflag;
	unsigned short yflag;
	// setup branches
	t0_event.SetupInput(tree);
	xppac_event.SetupInput(tree, "xppac.");
	tree->SetBranchAddress("xppac.xflag", &xflag);
	tree->SetBranchAddress("xppac.yflag", &yflag);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sxppac-center-%04u.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// delta X of PPAC
	TH1F hdx[3]{
		TH1F("hdx0", "#Deltax0", 1000, -10, 10),
		TH1F("hdx1", "#Deltax1", 1000, -10, 10),
		TH1F("hdx2", "#Deltax2", 1000, -10, 10)
	};
	// delta Y of PPAC
	TH1F hdy[3]{
		TH1F("hdy0", "#Deltay0", 1000, -10, 10),
		TH1F("hdy1", "#Deltay1", 1000, -10, 10),
		TH1F("hdy2", "#Deltay2", 1000, -10, 10)
	};

	// random number generator
	TRandom3 generator(tree->GetEntries());

	// total number of entries
	long long entries = tree->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling residual of run %u   0%%", run);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		tree->GetEntry(entry);
		if (
			t0_event.num != 1 || t0_event.mass[0] != 14
			|| xflag != 0x7 || yflag != 0x7
		) continue;

		// continuous d1x
		double d1x = t0_event.x[0] + generator.Rndm() - 0.5;
		double d1y = t0_event.y[0] + generator.Rndm() - 0.5;

		FitAndFill(ppac_xz, xppac_event.x, d1z, d1x, hdx);
		FitAndFill(ppac_yz, xppac_event.y, d1z, d1y, hdy);
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	CenterStatistics statistics(run, "xppac", tag);

	// fit offset
	for (int i = 0; i < 2; ++i) {
		TF1 fx(TString::Format("fx%d", i), "gaus", -5, 5);
		fx.SetParameter(0, 20);
		fx.SetParameter(1, 0.0);
		fx.SetParameter(2, 1.0);
		// fx.SetParameter(3, 1.0);
		hdx[i].Fit(&fx, "QR+");
		statistics.x_offset[i] = fx.GetParameter(1);

		TF1 fy(TString::Format("fy%d", i), "gaus", -5, 5);
		fy.SetParameter(0, 20);
		fy.SetParameter(1, 0.0);
		fy.SetParameter(2, 1.0);
		// fy.SetParameter(3, 1.0);
		hdy[i].Fit(&fy, "QR+");
		statistics.y_offset[i] = fy.GetParameter(1);
	}

	for (size_t i = 0; i < 3; ++i) {
		hdx[i].Write();
		hdy[i].Write();
	}
	// close files
	opf.Close();
	t0_file.Close();

	statistics.Write();
	statistics.Print();
	return 0;
}


/// @brief show offsets of multiple runs in graph
/// @param[in] run the start run number
/// @param[in] end_run the last run number
/// @param[in] tag trigger tag
/// @returns 0 if success, -1 otherwise
///
int ShowOffsets(
	unsigned int run,
	unsigned int end_run,
	const std::string &tag
) {
	// input csv file name
	TString input_file_name;
	input_file_name.Form(
		"%sstatistics/center.csv", kGenerateDataPath
	);
	// input file
	std::ifstream fin(input_file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open file " << input_file_name << " failed.\n";
		return -1;
	}

	// output root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sxppac-center-%s%04u-%04u.root",
		kGenerateDataPath,
		kShowDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run,
		end_run
	);
	// output root file
	TFile opf(output_file_name, "recreate");
	// PPAC x offset
	TGraph gx[3];
	// PPAC y offset
	TGraph gy[3];
	// PPAC x offsets
	TMultiGraph mgx;
	// PPAC y offsets
	TMultiGraph mgy;
	// PPAC offsets
	TMultiGraph mg;

	// add graph to multigraph
	for (size_t i = 0; i < 3; ++i) {
		mgx.Add(gx+i);
		mgy.Add(gy+i);
		mg.Add(gx+i);
		mg.Add(gy+i);
	}
	// colors
	const int colors[3] = {1, 600, 632};
	for (size_t i = 0; i < 3; ++i) {
		// set line color
		gx[i].SetLineColor(colors[i]);
		gy[i].SetLineColor(colors[i]);
		// set marker style
		gx[i].SetMarkerStyle(20);
		gy[i].SetMarkerStyle(21);
		// set marker color
		gx[i].SetMarkerColor(colors[i]);
		gy[i].SetMarkerColor(colors[i]);
	}
	// make x legend
	TLegend *legend_x = new TLegend(0.8, 0.8, 0.95, 0.95);
	for (size_t i = 0; i < 3; ++i) {
		legend_x->AddEntry(gx+i, TString::Format("x%ld", i));
	}
	mgx.GetListOfFunctions()->Add(legend_x);
	// make y legend
	TLegend *legend_y = new TLegend(0.8, 0.8, 0.95, 0.95);
	for (size_t i = 0; i < 3; ++i) {
		legend_y->AddEntry(gy+i, TString::Format("y%ld", i));
	}
	mgy.GetListOfFunctions()->Add(legend_y);
	// make x and y legend
	TLegend *legend_xy = new TLegend(0.8, 0.75, 0.95, 0.95);
	for (size_t i = 0; i < 3; ++i) {
		legend_xy->AddEntry(gx+i, TString::Format("x%ld", i));
		legend_xy->AddEntry(gy+i, TString::Format("y%ld", i));
	}
	mg.GetListOfFunctions()->Add(legend_xy);

	// buffer to read lines
	std::string buffer;
	// read first title line
	std::getline(fin, buffer);
	// statistics entry to read
	CenterStatistics statistics;
	// read lines
	fin >> statistics;
	while (fin.good()) {
		if (
			statistics.Run() >= run
			&& statistics.Run() <= end_run
			&& (
				(statistics.Tag() == "-" && tag.empty())
				|| (statistics.Tag() ==  tag)
			)
			&& statistics.Telescope() == "xppac"
		) {
			for (size_t i = 0; i < 3; ++i) {
				gx[i].AddPoint(statistics.Run(), statistics.x_offset[i]);
				gy[i].AddPoint(statistics.Run(), statistics.y_offset[i]);
			}
		}
		fin >> statistics;
	}

	// save graphs
	for (size_t i = 0; i < 3; ++i) {
		gx[i].Write(TString::Format("gx%ld", i));
		gy[i].Write(TString::Format("gy%ld", i));
	}
	mgx.Write("mgx");
	mgy.Write("mgy");
	mg.Write("mg");
	// close files
	opf.Close();
	return 0;
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
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag);

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
		std::cerr << "Error: Miss run argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	// run number
	unsigned int run = atoi(argv[pos_start]);
	unsigned int end_run = run;
	if (pos_start + 1 < argc) {
		end_run = atoi(argv[pos_start+1]);
	}

	if (end_run == run) {
		if (CalculateOffset(run, tag)) return -1;
	} else {
		if (ShowOffsets(run, end_run, tag)) return -1;
	}

	return 0;
}
