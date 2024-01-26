#include <iostream>

#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/detect_event.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options]\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set detected data type.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] tag tag get from arguments
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &tag
) {
	// initialize
	help = false;
	tag.clear();
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
			tag = argv[result];
		} else {
			return -result;
		}
	}
	return result;
}

int main(int argc, char **argv) {
	bool help = false;
	std::string tag = "";
	ParseArguments(argc, argv, help, tag);

	// need help
	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}

	// input file name
	TString input_file_name = TString::Format(
		"%s%sdetect.root",
		kGenerateDataPath, kSimulateDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	// input event
	DetectEvent detect;
	// setup input branches
	detect.SetupInput(ipt);

	// single PPAC tracking data
	bool single_ppac_track = false;
	double ppac_tx[3], ppac_ty[3];
	if (tag.empty()) {
		// do nothing
	} else if (
		tag.length() >= 4
		&& tag.substr(0, 3) == "spt"
	) {
		single_ppac_track = true;
		ipt->AddFriend("spt=tree", TString::Format(
			"%s%ssingle-ppac-track-%s.root",
			kGenerateDataPath, kSimulateDir, tag.substr(3).c_str()
		));
		// setup output branches
		ipt->SetBranchAddress("spt.ppac_tx", ppac_tx);
		ipt->SetBranchAddress("spt.ppac_ty", ppac_ty);
	} else {
		std::cerr << "Error: Invalid tag: " << tag << "\n";
		return -1;
	}

	// output file name
	TString output_file_name = TString::Format(
		"%s%srebuild%s.root",
		kGenerateDataPath,
		kSimulateDir,
		tag.empty() ? "" : ("-"+tag).c_str()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of Q value
	TH1F hist_q("hq", "Q", 300, -23, -8);
	// output tree
	TTree opt("tree", "simulated rebuild");
	// output data
	// kinetic energy of particles
	double be_kinetic, he_kinetic, d_kinetic, c_kinetic;
	// rebuild Q value
	double rebuild_q;
	// rebuild 10Be excited, 0-groud state, 1-3.5MeV, 2-6MeV
	int rebuild_state;
	// excited energy of 14C
	double c_excited;
	// setup output branches
	opt.Branch("be_kinetic", &be_kinetic, "bek/D");
	opt.Branch("he_kinetic", &he_kinetic, "hek/D");
	opt.Branch("d_kinetic", &d_kinetic, "dk/D");
	opt.Branch("c_kinetic", &c_kinetic, "ck/D");
	opt.Branch("q", &rebuild_q, "q/D");
	opt.Branch("state", &rebuild_state, "state/I");
	opt.Branch("c_excited", &c_excited, "cex/D");
	opt.Branch("vaild", &detect.valid, "valid/I");

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Rebuilding   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get data
		ipt->GetEntry(entry);

		// reaction point
		double tx = detect.tx;
		double ty = detect.ty;
		if (single_ppac_track) {
			tx = ppac_tx[0];
			ty = ppac_ty[0];
		}

		// rebuild Q value spectrum
		// momentum vector of fragments
		ROOT::Math::XYZVector fp[2];
		for (size_t i = 0; i < 2; ++i) {
			fp[i] = ROOT::Math::XYZVector(
				detect.t0x[i][0] - tx,
				detect.t0y[i][0] - ty,
				detect.t0z[i][0]
			);
			// fragment mass
			double mass = i == 0 ? mass_10be : mass_4he;
			// fragment kinetic energy
			double kinetic = i == 0 ?
				detect.be_kinetic : detect.he_kinetic;
			// fragment momentum value
			double momentum = sqrt(
				pow(kinetic, 2.0) + 2.0 * kinetic * mass
			);
			fp[i] = fp[i].Unit() * momentum;
		}
		// recoil 2H momentum vector
		ROOT::Math::XYZVector rp(
			detect.tafx - tx,
			detect.tafy - ty,
			detect.tafz
		);
		// recoild 2H momentum value
		double recoil_momentum = sqrt(
			pow(detect.d_kinetic, 2.0) + 2.0 * detect.d_kinetic * mass_2h
		);
		rp = rp.Unit() * recoil_momentum;
		// rebuild beam momentum vector
		ROOT::Math::XYZVector bp = fp[0] + fp[1] + rp;
		// rebuild beam kinetic energy
		detect.c_kinetic = sqrt(bp.Dot(bp) + pow(mass_14c, 2.0)) - mass_14c;
		// Q value
		rebuild_q = detect.be_kinetic + detect.he_kinetic
			+ detect.d_kinetic - detect.c_kinetic;

		// rebuild 14C excited energy
		// rebuild state
		rebuild_state = -1;
		double rebuild_be_excited = 0.0;
		if (rebuild_q > -13 && rebuild_q < -10) {
			rebuild_state = 0;
			rebuild_be_excited = 0.0;
		} else if (rebuild_q > -16 && rebuild_q < -14) {
			rebuild_state = 1;
			rebuild_be_excited = 3.368;
		} else if (rebuild_q > -19 && rebuild_q < -16.7) {
			rebuild_state = 2;
			rebuild_be_excited = 6.179;
		}
		// calcualted excited 14C momentum vecotr
		ROOT::Math::XYZVector cbp = fp[0] + fp[1];
		// excited 14C momentum value
		double c_momentum = cbp.R();
		// excited 14C total energy
		double c_energy =
			(detect.be_kinetic + mass_10be + rebuild_be_excited)
			+ (detect.he_kinetic + mass_4he);
		// excited 14C mass
		double excited_c_mass = sqrt(
			pow(c_energy, 2.0) - pow(c_momentum, 2.0)
		);
		// excited energy of 14C
		c_excited = excited_c_mass - mass_14c;

		// fill histpogram
		if (detect.valid == 7) hist_q.Fill(rebuild_q);
		// fill tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms
	hist_q.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}