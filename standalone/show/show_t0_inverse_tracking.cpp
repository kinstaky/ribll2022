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

#include "include/event/t0_event.h"
#include "include/event/particle_event.h"
#include "include/event/threebody_info_event.h"
#include "include/statistics/center_statistics.h"
#include "include/optimize_utilities.h"

using namespace ribll;

constexpr double d1z = 100.0;
constexpr double d2z = 111.76;
const double t0z[3] = {100.0, 111.76, 123.52};
// constexpr double ppac_xz[4] = {-695.2, -633.7, -454.2, -275.2};
// constexpr double ppac_yz[4] = {-689.2, -627.7, -448.2, -269.2};
constexpr double ppac_correctx[3] = {0.0, -2.23, -3.40};
constexpr double ppac_correcty[3] = {0.0, 0.84, 1.78};
constexpr int change_run = 717;

double SimpleFit(const double *x, const double *y, const int n, double &k, double &b) {
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
	const double *z,
	const double *x,
	const double z1,
	const double x1,
	const double z2,
	const double x2,
	TH1F *ppactx,
	TH1F *t0tx,
	TH1F *diffx,
	TH1F *t0dx
) {
	// fit 3 PPAC points and calculate T0D1 offset
	double ppac_k, ppac_b;
	SimpleFit(z, x, ppac_k, ppac_b);
	t0dx->Fill(ppac_b+z1*ppac_k-x1);
	ppactx->Fill(ppac_b);
	double t0b = (z1*x2 - z2*x1) / (z1 - z2);
	t0tx->Fill(t0b);
	diffx->Fill(ppac_b-t0b);
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
		"  -i                Use threebody info, 4He, 10Be.\n"
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
/// @param[out] info use threebody information
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag,
	bool &info
) {
	// initialize
	help = false;
	trigger_tag.clear();
	info = false;
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
			info = true;
		} else {
			return -result;
		}
	}
	return result;
}


int InverseTrack(unsigned int run, const std::string &tag) {
	// t0 file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
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
	T0Event t0_event;
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
		"%s%st0-inverse-%04u.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// reaction point X from PPAC
	TH1F txppac("ptx", "target x from PPAC", 300, -30, 30);
	// reaction point Y from PPAC
	TH1F typpac("pty", "target y from PPAC", 300, -30, 30);
	// reaction point X from T0
	TH1F txt0("ttx", "target x from T0", 300, -30, 30);
	// reaction point Y from T0
	TH1F tyt0("tty", "target y from T0", 300, -30, 30);
	// dX from PPAC and T0 at T0D1
	TH1F diffx("dx", "target x difference", 300, -30, 30);
	// dY from PPAC and T0 at T0D1
	TH1F diffy("dy", "target y difference", 300, -30, 30);
	// reaction point dX from PPAC and T0
	TH1F t0dx("t0dx", "x position difference", 300, -30, 30);
	// reaction point dY from PPAC and T0
	TH1F t0dy("t0dy", "y position difference", 300, -30, 30);

	// random number generator
	TRandom3 generator(tree->GetEntries());

	double using_xz[3] = {ppac_xz[0], ppac_xz[2], ppac_xz[3]};
	double using_yz[3] = {ppac_yz[0], ppac_yz[2], ppac_yz[3]};
	if (run >= change_run) {
		using_xz[0] = ppac_xz[1];
		using_yz[0] = ppac_yz[1];
	}

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

		// continuous d1 position
		double d1x = t0_event.x[0][0] + generator.Rndm() - 0.5;
		double d1y = t0_event.y[0][0] + generator.Rndm() - 0.5;

		// continuous d2 position
		double d2x = t0_event.x[0][1] + generator.Rndm()*2.0 - 1.0;
		double d2y = t0_event.y[0][1] + generator.Rndm()*2.0 - 1.0;

		// correct PPAC
		for (int i = 0; i < 3; ++i) {
			xppac_event.x[i] += ppac_correctx[i];
			xppac_event.y[i] += ppac_correcty[i];
		}

		FitAndFill(
			using_xz, xppac_event.x,
			d1z, d1x, d2z, d2x,
			&txppac, &txt0, &diffx, &t0dx
		);
		FitAndFill(
			using_yz, xppac_event.y,
			d1z, d1y, d2z, d2y,
			&typpac, &tyt0, &diffy, &t0dy
		);
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	txppac.Write();
	typpac.Write();
	txt0.Write();
	tyt0.Write();
	diffx.Write();
	diffy.Write();
	t0dx.Write();
	t0dy.Write();
	// close files
	opf.Close();
	t0_file.Close();
	return 0;
}


int InverseTrackInfo() {
	// input file name
	TString info_file_name;
	info_file_name.Form(
		"%s%sthreebody.root",
		kGenerateDataPath,
		kInformationDir
	);
	// info file
	TFile info_file(info_file_name, "read");
	// info tree
	TTree *tree = (TTree*)info_file.Get("tree");
	if (!tree) {
		std::cerr << "Error: Get tree from "
			<< info_file_name << " failed.\n";
		return -1;
	}
	// input info event
	ThreeBodyInfoEvent info;
	// setup branches
	info.SetupInput(tree);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-inverse-tb.root",
		kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// T0 10Be and 4He target dx
	TH1F t0_dx12("htdx12", "T0 10Be and 4He target dx", 100, -25, 25);
	// T0 10Be and 4He target dy
	TH1F t0_dy12("htdy12", "T0 10Be and 4He target dy", 100, -25, 25);
	// PPAC and T0 10Be tracking target dx
	TH1F ppac_t0_dx1(
		"hptdx1", "PPAC and T0 10Be tracking target dx", 100, -25, 25
	);
	// PPAC and T0 10Be tracking target dy
	TH1F ppac_t0_dy1(
		"hptdy1", "PPAC and T0 10Be tracking target dy", 100, -25, 25
	);
	// PPAC and T0 4He tracking target dx
	TH1F ppac_t0_dx2(
		"hptdx2", "PPAC and T0 4He tracking target dx", 100, -25, 25
	);
	// PPAC and T0 4He tracking target dy
	TH1F ppac_t0_dy2(
		"hptdy2", "PPAC and T0 4He tracking target dy", 100, -25, 25
	);
	TTree opt("tree", "T0 inverse tracking");
	// output data
	// PPAC reaction ponit
	double ptx, pty;
	// T0 10Be reaction point
	double betx, bety;
	// T0 4He reaction point
	double hetx, hety;
	// setup output branches
	opt.Branch("ptx", &ptx, "ptx/D");
	opt.Branch("pty", &pty, "pty/D");
	opt.Branch("betx", &betx, "betx/D");
	opt.Branch("bety", &bety, "bety/D");
	opt.Branch("hetx", &hetx, "hetx/D");
	opt.Branch("hety", &hety, "hety/D");

	// random number generator
	TRandom3 generator(tree->GetEntries());

	constexpr double using_xz[3] = {
		-695.2, -454.2, -275.2
	};
	constexpr double using_yz[3] = {
		-689.2, -448.2, -269.2
	};

	// total number of entries
	long long entries = tree->GetEntries();
	for (long long entry = 0; entry < entries; ++entry) {
		tree->GetEntry(entry);

		// correct PPAC
		for (int i = 0; i < 3; ++i) {
			info.ppac_x[i] += ppac_correctx[i];
			info.ppac_y[i] += ppac_correcty[i];
		}

		double xk, xb, yk, yb;
		TrackPpac(info.ppac_xflag, using_xz, info.ppac_x, xk, xb);
		TrackPpac(info.ppac_yflag, using_yz, info.ppac_y, yk, yb);

		double xk0, xb0, yk0, yb0;
		int len0 = info.layer[0] == 1 ? 2 : 3;
		SimpleFit(t0z, info.be_x, len0, xk0, xb0);
		SimpleFit(t0z, info.be_y, len0, yk0, yb0);

		double xk1, xb1, yk1, yb1;
		int len1 = info.layer[1] == 1 ? 2 : 3;
		SimpleFit(t0z, info.he_x, len1, xk1, xb1);
		SimpleFit(t0z, info.he_y, len1, yk1, yb1);

		// fill to histogram
		t0_dx12.Fill(xb1 - xb0);
		t0_dy12.Fill(yb1 - yb0);
		ppac_t0_dx1.Fill(xb0 - xb);
		ppac_t0_dy1.Fill(yb0- yb);
		ppac_t0_dx2.Fill(xb1 - xb);
		ppac_t0_dy2.Fill(yb1 - yb);

		// fill tree
		ptx = xb;
		pty = yb;
		betx = xb0;
		bety = yb0;
		hetx = xb1;
		hety = yb1;
		opt.Fill();
	}

	// save histograms
	t0_dx12.Write();
	t0_dy12.Write();
	ppac_t0_dx1.Write();
	ppac_t0_dy1.Write();
	ppac_t0_dx2.Write();
	ppac_t0_dy2.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	info_file.Close();
	return 0;
}


int main(int argc, char **argv) {
	if (argc < 1) {
		PrintUsage(argv[0]);
		return -1;
	}

	// help flag
	bool help = false;
	// trigger tag
	std::string tag;
	// info flag
	bool info = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag, info);

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

	if (info) {
		if (InverseTrackInfo()) {
			std::cerr << "Error: Inverse track threebody info error.\n";
			return -1;
		}
	} else {
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
		if (run == end_run) {
			if (InverseTrack(run, tag)) {
				std::cerr << "Error: Inverse track error.\n";
				return -1;
			}
		}
	}

	return 0;
}
