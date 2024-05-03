#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Math/Vector3D.h>
#include <TH1F.h>

#include "include/event/threebody_info_event.h"

using namespace ribll;

int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody-9Be.root",
		kGenerateDataPath,
		kInformationDir
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

	// output file name
	TString output_file_name = TString::Format(
		"%s%sthreebody-9Be.root",
		kGenerateDataPath,
		kSpectrumDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "9Be+4He+n + 2H spectrum");
	// output data
	int ppac_flag, taf_flag, bind, hole;
	double q;
	// setup output branches
	opt.Branch("ppac_flag", &ppac_flag, "pflag/I");
	opt.Branch("taf_flag", &taf_flag, "tflag/I");
	opt.Branch("bind", &bind, "bind/I");
	opt.Branch("hole", &hole, "hole/I");
	opt.Branch("q", &q, "q/D");


	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		// filter events
		if (event.taf_flag != 0) continue;
		if ((event.ppac_flag & 1) == 0) continue;
		if (event.xppac_track[0] < 2 || event.xppac_track[1] < 2) continue;
		if (event.bind != 0) continue;
		if (event.hole[0] || event.hole[1]) continue;

		// target position
		const double tx = event.xptx;
		const double ty = event.xpty;


		// 9Be momentum
		double be_momentum = MomentumFromKinetic(mass_9be, event.t0_energy[0]);
		// 10Be momentum vector
		ROOT::Math::XYZVector p_be(
			event.be_x[0] - tx,
			event.be_y[0] - ty,
			100.0
		);
		p_be = p_be.Unit() * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, event.t0_energy[1]);
		// 4He momentum vector
		ROOT::Math::XYZVector p_he(
			event.he_x[0] - tx,
			event.he_y[0] - ty,
			100.0
		);
		p_he = p_he.Unit() * he_momentum;

		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, event.taf_energy);
		// 2H momentum vector
		ROOT::Math::XYZVector p_d(
			event.d_x - tx,
			event.d_y - ty,
			135.0
		);
		p_d = p_d.Unit() * d_momentum;

		// beam 14C momentum vector
		ROOT::Math::XYZVector p_beam(
			event.beam_px, event.beam_py, event.beam_pz
		);

		// neutron momentum
		ROOT::Math::XYZVector p_n = p_beam - p_be - p_he - p_d;
		// neutron kinetic energy
		double neutron_kinetic = sqrt(p_n.Mag2() + pow(mass_n, 2.0)) - mass_n;

		q = event.t0_energy[0] + event.t0_energy[1]
			+ event.taf_energy + neutron_kinetic - 389.5;

		opt.Fill();
	}

	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}