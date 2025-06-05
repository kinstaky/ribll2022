#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <Math/Vector3D.h>
#include <TRandom3.h>

#include "include/defs.h"

using namespace ribll;

constexpr double be_excited_energy[4] = {
	0.0, 3.368, 6.179, 7.542
};

int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v3.root",
		kGenerateDataPath,
		kSpectrumDir
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
	int valid, be_state;
	double be_kinetic, he_kinetic;
	double be_dx, be_dy, be_dz, he_dx, he_dy, he_dz;
	// setup input branches
	ipt->SetBranchAddress("valid", &valid);
	ipt->SetBranchAddress("be_state", &be_state);
	ipt->SetBranchAddress("be_kinetic", &be_kinetic);
	ipt->SetBranchAddress("he_kinetic", &he_kinetic);
	ipt->SetBranchAddress("be_dx", &be_dx);
	ipt->SetBranchAddress("be_dy", &be_dy);
	ipt->SetBranchAddress("be_dz", &be_dz);
	ipt->SetBranchAddress("he_dx", &he_dx);
	ipt->SetBranchAddress("he_dy", &he_dy);
	ipt->SetBranchAddress("he_dz", &he_dz);

	// // output file name
	TString output_file_name = TString::Format(
		"%s%sevent-mixing.root",
		kGenerateDataPath,
		kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// event mixing histogram
	TH1F he[3] = {
		TH1F("he0", "excited energy", 50, 11.0, 26.0),
		TH1F("he1", "excited energy", 50, 11.0, 26.0),
		TH1F("he2", "excited energy", 50, 11.0, 26.0)
	};

	TRandom3 generator(0);

	for (int state = 0; state < 3; ++state) {
		std::vector<double> beks;
		std::vector<double> heks;
		std::vector<ROOT::Math::XYZVector> d_bes;
		std::vector<ROOT::Math::XYZVector> d_hes;


		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);

			if (valid != 0) continue;
			if (be_state != state) continue;

			beks.push_back(be_kinetic);
			heks.push_back(he_kinetic);
			d_bes.push_back(ROOT::Math::XYZVector(
				be_dx, be_dy, be_dz
			));
			d_hes.push_back(ROOT::Math::XYZVector(
				he_dx, he_dy, he_dz
			));
		}

		// loop
		int count = 0;
		while (count < 100000) {
			auto be_index = generator.Rndm() * beks.size();
			auto he_index = generator.Rndm() * heks.size();

			ROOT::Math::XYZVector d_be = d_bes[be_index];
			double be_momentum = MomentumFromKinetic(mass_10be, beks[be_index]);
			ROOT::Math::XYZVector p_be = d_be * be_momentum;
			ROOT::Math::XYZVector d_he = d_hes[he_index];
			double he_momentum = MomentumFromKinetic(mass_4he, heks[he_index]);
			ROOT::Math::XYZVector p_he = d_he * he_momentum;
			ROOT::Math::XYZVector p_c = p_be + p_he;
			double mass_excited_10be = mass_10be
				+ (state >= 0 ? be_excited_energy[state] : 0.0);
			double c_energy = (beks[be_index] + mass_excited_10be) + (heks[he_index] + mass_4he);
			double c_mass = sqrt(
				pow(c_energy, 2.0) - p_c.Mag2()
			);
			double excited_energy = c_mass - mass_14c;

			he[state].Fill(excited_energy);
			++count;
		}

	}


	opf.cd();
	he[0].Write();
	he[1].Write();
	he[2].Write();
	opf.Close();
	ipf.Close();

	return 0;
}