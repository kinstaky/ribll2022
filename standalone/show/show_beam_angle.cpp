// get beam angle from PPAC and rebuilt 10Be+4He+2H respectively
// compare their different to find whether is possible to distinct
// 9Be and 10Be

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/ppac_track.h"

using namespace ribll;



int main(int argc, char **argv) {
	if (argc > 2) {
		std::cout << "Usage: " << argv[0] << " mass\n"
			<< "  mass            Be mass number, 9 or 10, default is 10\n";
		return -1;
	}
	// Be mass number
	int mass = 10;;
	if (argc == 2) {
		mass = atoi(argv[1]);
	}


	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody-calc-%dBe.root",
		kGenerateDataPath, kInformationDir, mass
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
	// threebody information event
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);


	// output file name
	TString output_file_name = TString::Format(
		"%s%sthreebody-beam-angle-%dBe.root",
		kGenerateDataPath, kShowDir, mass
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "beam angle");
	// output data
	double px[2], py[2], pz[2];
	bool valid;
	// setup output branches
	opt.Branch("px", px, "px[2]/D");
	opt.Branch("py", py, "py[2]/D");
	opt.Branch("pz", pz, "pz[2]/D");
	opt.Branch("valid", &valid, "valid/O");


	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		// jump not so good events
		if (
			event.taf_flag != 0
			|| event.xppac_track[0] < 2
			|| event.xppac_track[1] < 2
			|| event.bind != 0
			|| event.hole[0]
			|| event.hole[1] 
		) {
			valid = false;
			opt.Fill();
			continue;
		}

		// set valid
		valid = true;

		// get beam angle from XIA PPAC
		// fitted trace parameters
		double xk, xb, yk, yb;
		// using xz and yz
		double using_ppac_xz[3] = {ppac_xz[0], ppac_xz[1], ppac_xz[2]};
		double using_ppac_yz[3] = {ppac_yz[0], ppac_yz[1], ppac_yz[2]};
		if (event.run >= ppac_change_run) {
			using_ppac_xz[0] = all_ppac_xz[1];
			using_ppac_yz[0] = all_ppac_yz[1];
		}
		// track
		TrackMultiplePpac(event.xppac_xflag, using_ppac_xz, event.xppac_x, xk, xb);
		TrackMultiplePpac(event.xppac_yflag, using_ppac_yz, event.xppac_y, yk, yb);
		// get vector
		ROOT::Math::XYZVector p1(xk, yk, 1.0);
		// unify	
		p1 = p1.Unit();
		// fill data
		px[0] = p1.X();
		py[0] = p1.Y();
		pz[0] = p1.Z();

		// get beam angle from 10Be, 4He and 2H
		// 10Be momentum
		double be_momentum = MomentumFromKinetic(mass_10be, event.t0_energy[0]);
		// 10Be momentum vector
		ROOT::Math::XYZVector p_be(
			event.be_x[0] - event.xptx,
			event.be_y[0] - event.xpty,
			100.0
		);
		p_be = p_be.Unit() * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, event.t0_energy[1]);
		// 4He momentum vector
		ROOT::Math::XYZVector  p_he(
			event.he_x[0] - event.xptx,
			event.he_y[0] - event.xpty,
			100.0
		);
		p_he = p_he.Unit() * he_momentum;

		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, event.taf_energy);
		// 2H momentum vector
		ROOT::Math::XYZVector p_d(
			event.d_x - event.xptx,
			event.d_y - event.xpty,
			135.0
		);
		p_d = p_d.Unit() * d_momentum;

		// get beam momentum from momentum conversion
		ROOT::Math::XYZVector p_beam = p_be + p_he + p_d;
		// unify
		ROOT::Math::XYZVector p2 = p_beam.Unit();
		// fill data
		px[1] = p2.X();
		py[1] = p2.Y();
		pz[1] = p2.Z();

		opt.Fill();
	}

	// save tre
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}