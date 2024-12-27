#include <iostream>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <Math/Vector3D.h>

#include "include/event/channel_v2_event.h"
#include "include/calculator/csi_energy_calculator.h"

using namespace ribll;


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] tag [recoil]\n"
		<< "  tag           Channel file tag.\n"
		<< "  recoil        Recoil particle mass, default is 2(2H).\n"
		<< "Options:\n"
		<< "  -h            Print this help information.\n"
		<< "  -c            Calculate recoil energy from TAFD only.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] calculate calculate recoil energy from TAFD only
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	bool &calculate
) {
	// initialize
	help = false;
	calculate = false;
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
		} else if (argv[result][1] == 'c') {
			// option of binding events
			calculate = true;
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
	// calculate flag
	bool calculate = false;
	// parse arguments
	int pos_start = ParseArguments(argc, argv, help, calculate);

	// need help
	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}
	// check arguments
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
		std::cerr << "Error: Miss tag argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	// get tag
	std::string tag(argv[pos_start]);
	// get image recoil mass
	int image_recoil_mass =
		pos_start + 1 == argc ?  2 : atoi(argv[pos_start+1]);

	if (image_recoil_mass != 1 && image_recoil_mass != 2) {
		std::cerr << "Error: Invalid image recoil mass number, "
			<< image_recoil_mass << "\n";
		return -1;
	}
	std::cout << "Image recoil particle is " << image_recoil_mass << "H\n";

	// input file name
	TString input_file_name = TString::Format(
		"%s%sC14-10Be-4He-%s-v2.root",
		kGenerateDataPath,
		kChannelDir,
		tag.c_str()
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
	// input event
	ChannelV2Event channel;
	// setup input branches
	channel.SetupInput(ipt);

	std::cout << "Read from " << input_file_name << "\n";

	// output file name
	TString output_file_name = TString::Format(
		"%s%sC14-10Be-4He-%s-image-%dH-q%s.root",
		kGenerateDataPath,
		kSpectrumDir,
		tag.c_str(),
		image_recoil_mass,
		calculate ? "-calc-taf" : ""
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "image recoil Q");
	// output data
	int taf_index, csi_index, tafd_front_strip;
	bool tafcsi_valid;
	double beam_kinetic, q;
	opt.Branch("taf_index", &taf_index, "tafi/I");
	opt.Branch("csi_index", &csi_index, "ci/I");
	opt.Branch("tafd_front_strip", &tafd_front_strip, "tafdfs/I");
	opt.Branch("tafcsi_valid", &tafcsi_valid, "cv/O");
	opt.Branch("beam_kinetic", &beam_kinetic, "bk/D");
	opt.Branch("q", &q, "q/D");

	// get mass
	const double beam_mass = mass_14c;
	const double recoil_mass = image_recoil_mass == 1 ? mass_1h : mass_2h;

	// calculators
	elc::CsiEnergyCalculator recoil_calculator(
		TString::Format("%dH", image_recoil_mass).Data()
	);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		if (channel.valid != 0) continue;

		taf_index = channel.taf_index;
		csi_index = channel.csi_index;
		tafd_front_strip = channel.tafd_front_strip;
		tafcsi_valid = channel.tafcsi_valid;

		// get recoil direction
		ROOT::Math::XYZVector recoil_direction(
			channel.recoil_x - channel.tx,
			channel.recoil_y - channel.ty,
			channel.recoil_z
		);
		recoil_direction = recoil_direction.Unit();
		// get recoil energy
		double recoil_kinetic = channel.recoil_kinetic;
		if (calculate) {
			recoil_kinetic = recoil_calculator.Energy(
				recoil_direction.Theta(),
				channel.tafd_energy,
				tafd_thickness[channel.taf_index]
			);
		}
		double recoil_momentum = MomentumFromKinetic(
			recoil_mass, recoil_kinetic
		);
		ROOT::Math::XYZVector recoil_p = recoil_direction * recoil_momentum;

		// get fragment momentum vector
		ROOT::Math::XYZVector fragment_p[2];
		for (int i = 0; i < 2; ++i) {
			fragment_p[i] = ROOT::Math::XYZVector(
				channel.fragment_x[i] - channel.tx,
				channel.fragment_y[i] - channel.ty,
				channel.fragment_z[i]
			);
			fragment_p[i] =
				fragment_p[i].Unit() * channel.fragment_momentum[i];
		}

		// get beam momentum and energy
		ROOT::Math::XYZVector beam_p =
			fragment_p[0] + fragment_p[1] + recoil_p;
		beam_kinetic = KineticFromMomentum(
			beam_mass, beam_p.R()
		);

		q = channel.fragment_kinetic[0] + channel.fragment_kinetic[1]
			+ recoil_kinetic - beam_kinetic;

		opt.Fill();
	}

	// save
	opf.cd();
	opt.Write();
	opf.Close();

	return 0;
}