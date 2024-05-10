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



int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%sthreebody-beam-angle.root",
		kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// 10Be beam angle
	TH1F hist_be10_angle("ha10", "10Be beam angle", 100, 0.999, 1);
	// 9Be beam angle
	TH1F hist_be9_angle("ha9", "9Be beam angle", 100, 0.999, 1);


	// input file name
	TString be10_info_file_name = TString::Format(
		"%s%sthreebody-10Be.root",
		kGenerateDataPath, kInformationDir
	);
	// input file
	TFile be10_info_file(be10_info_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)be10_info_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< be10_info_file_name << " failed.\n";
		return -1;
	}
	// threebody information event
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// add friend
	ipt->AddFriend("s=tree", TString::Format(
		"%s%sthreebody-10Be-2.root", kGenerateDataPath, kSpectrumDir
	));
	// input data
	double be_kinetic[4], he_kinetic[4], d_kinetic;
	ipt->SetBranchAddress("s.be_kinetic", be_kinetic);
	ipt->SetBranchAddress("s.he_kinetic", he_kinetic);
	ipt->SetBranchAddress("s.d_kinetic", &d_kinetic);

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
			continue;
		}

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

		// get beam angle from 10Be, 4He and 2H
		// 10Be momentum
		double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic[0]);
		// 10Be momentum vector
		ROOT::Math::XYZVector p_be(
			event.be_x[0] - event.xptx,
			event.be_y[0] - event.xpty,
			100.0
		);
		p_be = p_be.Unit() * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic[0]);
		// 4He momentum vector
		ROOT::Math::XYZVector  p_he(
			event.he_x[0] - event.xptx,
			event.he_y[0] - event.xpty,
			100.0
		);
		p_he = p_he.Unit() * he_momentum;

		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
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

		hist_be10_angle.Fill(p1.Dot(p2));
	}
	be10_info_file.Close();




	// input file name
	TString be9_info_file_name = TString::Format(
		"%s%sthreebody-9Be.root",
		kGenerateDataPath, kInformationDir
	);
	// input file
	TFile be9_info_file(be9_info_file_name, "read");
	// input tree
	ipt = (TTree*)be9_info_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< be9_info_file_name << " failed.\n";
		return -1;
	}
	// setup input branches
	event.SetupInput(ipt);

	// add friend
	ipt->AddFriend("s=tree", TString::Format(
		"%s%sthreebody-9Be-2.root", kGenerateDataPath, kSpectrumDir
	));
	ipt->SetBranchAddress("s.be_kinetic", be_kinetic);
	ipt->SetBranchAddress("s.he_kinetic", he_kinetic);
	ipt->SetBranchAddress("s.d_kinetic", &d_kinetic);

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
			continue;
		}


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

		// get beam angle from 10Be, 4He and 2H
		// 10Be momentum
		double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic[3]);
		// 10Be momentum vector
		ROOT::Math::XYZVector p_be(
			event.be_x[0] - event.xptx,
			event.be_y[0] - event.xpty,
			100.0
		);
		p_be = p_be.Unit() * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic[3]);
		// 4He momentum vector
		ROOT::Math::XYZVector  p_he(
			event.he_x[0] - event.xptx,
			event.he_y[0] - event.xpty,
			100.0
		);
		p_he = p_he.Unit() * he_momentum;

		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
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

		hist_be9_angle.Fill(p1.Dot(p2));
	}
	be9_info_file.Close();

	// save histograms
	opf.cd();
	hist_be10_angle.Write();
	hist_be9_angle.Write();
	// close files
	opf.Close();
	return 0;
}