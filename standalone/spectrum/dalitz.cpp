#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <Math/Vector3D.h>

#include "include/defs.h"
#include "include/event/threebody_info_event.h"

using namespace ribll;

constexpr double mass_6li = 6.0151223 * u;
constexpr double mass_12b = 12.011609 * u;

struct DalitzEvent {
	double c14_ex;
	double li6_ex;
	double b12_ex;
};

int FillV2(
	TH1F &hist_li_ex,
	TH1F &hist_b_ex,
	TTree &tree,
	DalitzEvent &event
) {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody-10Be-dv2-2.root", kGenerateDataPath, kSpectrumDir
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
	int valid;
	double be_kinetic[4], he_kinetic[4], d_kinetic, c_kinetic[4], q[4];
	double c_ex[3][4];
	double be_dx, be_dy, be_dz;
	double he_dx, he_dy, he_dz;
	double d_dx, d_dy, d_dz;
	// setup input branches
	ipt->SetBranchAddress("valid", &valid);
	ipt->SetBranchAddress("be_kinetic_target", be_kinetic);
	ipt->SetBranchAddress("he_kinetic_target", he_kinetic);
	ipt->SetBranchAddress("d_kinetic_target", &d_kinetic);
	ipt->SetBranchAddress("c_kinetic_target", c_kinetic);
	ipt->SetBranchAddress("q", q);
	ipt->SetBranchAddress("stateless_excited_energy", c_ex);
	ipt->SetBranchAddress("be_dx", &be_dx);
	ipt->SetBranchAddress("be_dy", &be_dy);
	ipt->SetBranchAddress("be_dz", &be_dz);
	ipt->SetBranchAddress("he_dx", &he_dx);
	ipt->SetBranchAddress("he_dy", &he_dy);
	ipt->SetBranchAddress("he_dz", &he_dz);
	ipt->SetBranchAddress("d_dx", &d_dx);
	ipt->SetBranchAddress("d_dy", &d_dy);
	ipt->SetBranchAddress("d_dz", &d_dz);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		if (valid != 0) continue;
		double be10_ex = 0.0;
		int be_state = 0;
		if (q[0] < -11 && q[0] > -13) {
			be10_ex = 0.0;
			be_state = 0;
		} else if (q[0] < -14.5 && q[0] > -16.0) {
			be10_ex = 3.368;
			be_state = 1;
		} else if (q[0] < -17.0 && q[0] > -20.0) {
			be10_ex = 6.179;
			be_state = 2;
		} else {
			continue;
		}
		// 10Be direction vector
		ROOT::Math::XYZVector d_be(be_dx, be_dy, be_dz);
		// 4He direction vector
		ROOT::Math::XYZVector d_he(he_dx, he_dy, he_dz);
		// 2H direction vector
		ROOT::Math::XYZVector d_d(d_dx, d_dy, d_dz);
		// 10Be momentum
		double be_momentum = MomentumFromKinetic(
			mass_10be+be10_ex, be_kinetic[0]
		);
		ROOT::Math::XYZVector p_be = d_be * be_momentum;
		// 4He momentum
		double he_momentum = MomentumFromKinetic(
			mass_4he, he_kinetic[0]
		);
		ROOT::Math::XYZVector p_he = d_he * he_momentum;
		// 2H momentum
		double d_momentum = MomentumFromKinetic(
            mass_2h, d_kinetic
        );
        ROOT::Math::XYZVector p_d = d_d * d_momentum;

		// excited 6Li momentum vector
		ROOT::Math::XYZVector p_li = p_d + p_he;
		// excited 6Li momentum
		double li_momentum = p_li.R();
		// excited 6Li total energy
		double li_energy = (he_kinetic[0] + mass_4he) + (d_kinetic + mass_2h);
		// excited 6Li mass
		double li_mass = sqrt(
			pow(li_energy, 2.0) - pow(li_momentum, 2.0)
		);
		// excited energy of 6Li
		double li_ex = li_mass - mass_6li;

		// excited 12B momentum vector
		ROOT::Math::XYZVector p_b = p_be + p_d;
		// excited 12B momentum
		double b_momentum = p_b.R();
		// excited 12B total energy
		double b_energy = (be_kinetic[0] + be10_ex + mass_10be)
			+ (d_kinetic + mass_2h);
		// excited 12B mass
		double b_mass = sqrt(
			pow(b_energy, 2.0)- pow(b_momentum, 2.0)
		);
		// excited energy of 12B
		double b_ex = b_mass - mass_12b;

		hist_li_ex.Fill(li_ex);
		hist_b_ex.Fill(b_ex);
		event.b12_ex = b_ex;
		event.li6_ex = li_ex;
		event.c14_ex = c_ex[be_state][0];
		tree.Fill();
	}

	// close file
	ipf.Close();
	return 0;
}


int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%sdalitz.root",
		kGenerateDataPath, kSpectrumDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// excitation spectrum
	TH1F hist_li6_ex("li6ex", "6Li excitation energy", 100, -100, 100);
	TH1F hist_b12_ex("b12ex", "12B excitation energy", 100, -100, 100);
	// output tree
	TTree opt("tree", "dalitz");
	// output event
	DalitzEvent event;
	// setup output branches
	opt.Branch("c14_ex", &event.c14_ex, "cex/D");
	opt.Branch("li6_ex", &event.li6_ex, "liex/D");
	opt.Branch("b12_ex", &event.b12_ex, "bex/D");

	FillV2(hist_li6_ex, hist_b12_ex, opt, event);

	// save
	opf.cd();
	hist_li6_ex.Write();
	hist_b12_ex.Write();
	opt.Write();
	// close file
	opf.Close();
	return 0;
}