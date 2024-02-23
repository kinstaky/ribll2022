#include <cmath>
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/ppac_track.h"
#include "include/event/threebody_info_event.h"


using namespace ribll;


double ThreeBodyProcess(
	const ThreeBodyInfoEvent &event,
	const double a0,
	const double a1,
	const double a2,
	double &be_kinetic,
	double &he_kinetic,
	double &d_kinetic,
	double &c_kinetic
) {
	// 10Be kinetic energy
	// T0D1 energy
	be_kinetic = t0_param[0][0] + t0_param[0][1] * event.be_channel[0];
	// T0D2 energy
	be_kinetic += t0_param[1][0] + t0_param[1][1] * event.be_channel[1];
	// T0D3 energy
	if (event.layer[0] > 1) {
		be_kinetic += t0_param[2][0] + t0_param[2][1] * event.be_channel[2];
	}
	// T0S1 energy
	if (event.layer[0] > 2) {
		be_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
	}
	// T0S2 energy
	if (event.layer[0] > 3) {
		be_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
	}
	// T0S3 energy
	if (event.layer[0] > 4) {
		be_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
	}

	// 4He kinetic energy
	// T0D1 energy
	he_kinetic = t0_param[0][0] + t0_param[0][1] * event.he_channel[0];
	// T0D2 energy
	he_kinetic += t0_param[1][0] + t0_param[1][1] * event.he_channel[1];
	// T0D3 energy
	if (event.layer[1] > 1) {
		he_kinetic += t0_param[2][0] + t0_param[2][1] * event.he_channel[2];
	}
	// T0S1 energy
	if (event.layer[1] > 2) {
		he_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
	}
	// T0S2 energy
	if (event.layer[1] > 3) {
		he_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
	}
	// T0S3 energy
	if (event.layer[1] > 4) {
		he_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
	}

	// 2H kinetic energy
	d_kinetic = event.tafd_energy + pow(
		(event.csi_channel - a2) / a0,
		1.0 / a1
	);

	double ppac_correct[2][3] = {
		{0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0}
	};
	PpacOffsetX(ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]);
	PpacOffsetY(ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]);
	// calculate reaction point
	double ppac_cx[3], ppac_cy[3];
	for (int i = 0; i < 3; ++i) {
		ppac_cx[i] = event.xppac_x[i] - ppac_correct[0][i];
		ppac_cy[i] = event.xppac_y[i] - ppac_correct[1][i];
	}
	// slope and intercept
	double xk, yk, xb, yb;
	TrackMultiplePpac(event.xppac_xflag, ppac_xz, ppac_cx, xk, xb);
	TrackMultiplePpac(event.xppac_yflag, ppac_yz, ppac_cy, yk, yb);

	// 10Be momentum
	double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - xb,
		event.be_y[0] - yb,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - xb,
		event.he_y[0] - yb,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		event.d_x - xb,
		event.d_y - yb,
		135.0
	);
	p_d = p_d.Unit() * d_momentum;

	// beam 14C momentum vector
	ROOT::Math::XYZVector p_c = p_be + p_he + p_d;
// std::cout << p_be.X() << ", " << p_be.Y() << ", " << p_be.Z() << "\n"
// 	<< p_he.X() << ", " << p_he.Y() << ", " << p_he.Z() << "\n"
// 	<< p_d.X() << ", " << p_d.Y() << ", " << p_d.Z() << "\n"
// 	<< p_c.X() << ", " << p_c.Y() << ", " << p_c.Z() << "\n";

	// 14C momentum
	double c_momentum = p_c.R();
	// 14C kinetic energy
	c_kinetic =
		sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

	// three-fold Q value
	double q = be_kinetic + he_kinetic + d_kinetic - c_kinetic;

	return q;
}




void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options]\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -c                Set calculate T0D2 flags.\n"
		"  -i                Iteration mode.\n";
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help
) {
	// initialize
	help = false;
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
		}else {
			return -result;
		}
	}
	return result;
}


int main(int argc, char **argv) {

	// help flag
	bool help = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help);

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

	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody.root", kGenerateDataPath, kInformationDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input data
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file
	TFile opf(
		TString::Format(
			"%s%sthreebody-group.root",
			kGenerateDataPath,
			kOptimizeDir
		),
		"recreate"
	);
	// group parameters Q value spectrum
	std::vector<TH1F> group_q_spectrum;
	std::vector<TGraph> group_g;
	// output tree
	TTree opt("tree", "group Q mean and sigma");
	// output data
	double q_mean;
	double q_sigma;
	// setup output branches
	opt.Branch("q_mean", &q_mean, "mean/D");
	opt.Branch("q_sigma", &q_sigma, "sigma/D");


	std::vector<double> csi_parameters[12][3];

	// read CsI calibration parameters from file
	for (int i = 0; i < 12; ++i) {
		std::ifstream fin(TString::Format(
			"%s%scsi-group-parameters-%d.txt",
			kGenerateDataPath,
			kOptimizeDir,
			i
		).Data());
		double a0, a1, a2;
		while (fin.good()) {
			fin >> a0 >> a1 >> a2;
			csi_parameters[i][0].push_back(a0);
			csi_parameters[i][1].push_back(a1);
			csi_parameters[i][2].push_back(a2);
		}
		fin.close();
	}
	int groups = csi_parameters[0][0].size();


	// kinetic energy
	double be_kinetic, he_kinetic, d_kinetic, c_kinetic;
	// loop to calculate Q value
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		group_q_spectrum.emplace_back(
			TString::Format("hgq%lld", entry),
			"group Q",
			300, -23, -8
		);
		group_g.emplace_back();

		for (int i = 0; i < groups; ++i) {
			double q_value = ThreeBodyProcess(
				event,
				csi_parameters[event.csi_index][0][i],
				csi_parameters[event.csi_index][1][i],
				csi_parameters[event.csi_index][2][i],
				be_kinetic, he_kinetic, d_kinetic, c_kinetic
			);

			group_q_spectrum[entry].Fill(q_value);
			group_g[entry].AddPoint(csi_parameters[event.csi_index][1][i], q_value);
		}

		q_mean = group_q_spectrum[entry].GetMean();
		q_sigma = group_q_spectrum[entry].GetStdDev(1);
		// fillt to tree
		opt.Fill();
	}


	// save histograms
	opf.cd();
	for (size_t i = 0; i < group_q_spectrum.size(); ++i) {
		group_q_spectrum[i].Write();
		group_g[i].Write(TString::Format("gg%ld", i));
	}
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}