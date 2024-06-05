#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <Math/Vector3D.h>
#include <TRandom3.h>

#include "include/event/threebody_info_event.h"
#include "include/ppac_track.h"

using namespace ribll;

int main(int argc, char **argv) {
	std::string suffix = "";
	if (argc > 1) {
		suffix = std::string(argv[1]);
	}
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody%s.root",
		kGenerateDataPath,
		kInformationDir,
		suffix.empty() ? "" : ("-"+suffix).c_str()
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

	// spectrum v2 file name
	TString spectrum_file_name = TString::Format(
		"%s%sthreebody%s-2.root",
		kGenerateDataPath,
		kSpectrumDir,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	// add friend
	ipt->AddFriend("s=tree", spectrum_file_name);
	// additional data
	int be_state[4];
	// setup input branches
	ipt->SetBranchAddress("s.be_state", be_state);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sppac-effect%s-2.root",
		kGenerateDataPath,
		kShowDir,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// target position, Q value, and exicted energy
	// in different target position case
	// 0: multiple PPAC track
	// 1: single PPAC track
	// 2: gaussian random
	// 3: fixed zero point 
	double tx[4], ty[4];
	double q[4], ex[4];
	// output tree
	TTree opt("tree", "PPAC effect 2");
	// setup output branches
	opt.Branch("target_x", tx, "tx[4]/D");
	opt.Branch("target_y", ty, "ty[4]/D");
	opt.Branch("q", q, "q[4]/D");
	opt.Branch("ex", ex, "ex[4]/D");

	// random number generator
	TRandom3 generator(639234);

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Processing   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// showing process
		if (entry / entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get data
		ipt->GetEntry(entry);

		if (event.xppac_track[0] < 2 || event.xppac_track[1] < 2) continue;
		if (event.taf_flag != 0) continue;
		if (event.bind != 0) continue;
		if (event.hole[0] || event.hole[1]) continue;
		if (be_state[0] < 0 || be_state[0] > 2) continue;

		// actually used PPAC
		double using_ppac_xz[3] = {ppac_xz[0], ppac_xz[1], ppac_xz[2]};
		double using_ppac_yz[3] = {ppac_yz[0], ppac_yz[1], ppac_yz[2]};
		if (event.run >= ppac_change_run) {
			using_ppac_xz[0] = all_ppac_xz[1];
			using_ppac_yz[0] = all_ppac_yz[1];
		}
		double be10_excited_energy = 0.0;
		if (be_state[0] == 1) be10_excited_energy = 3.368;
		else if (be_state[0] == 2) be10_excited_energy = 6.179;

		for (int i = 0; i < 4; ++i) {
			if (i == 0) {
				double xk, xb, yk, yb;
				TrackMultiplePpac(
					event.xppac_xflag, using_ppac_xz, event.xppac_x, xk, xb
				);
				TrackMultiplePpac(
					event.xppac_yflag, using_ppac_yz, event.xppac_y, yk, yb
				);
				tx[i] = xb;
				ty[i] = yb;
			} else if (i == 1) {
				if (event.xppac_xflag & 0x1) {
					tx[i] = DeutronRelativeApproximateTrack(
						event.t0_energy[0], event.t0_energy[1],
						event.taf_energy, using_ppac_xz[0],
						event.be_x[0], event.he_x[0], event.d_x, event.d_y,
						event.xppac_x[0]
					);
				} else if (event.xppac_xflag & 0x2) {
					tx[i] = DeutronRelativeApproximateTrack(
						event.t0_energy[0], event.t0_energy[1],
						event.taf_energy, using_ppac_xz[1],
						event.be_x[0], event.he_x[0], event.d_x, event.d_y,
						event.xppac_x[1]
					);					
				}
				if (event.xppac_yflag & 0x1) {
					ty[i] = DeutronRelativeApproximateTrack(
						event.t0_energy[0], event.t0_energy[1],
						event.taf_energy, using_ppac_yz[0],
						event.be_y[0], event.he_y[0], event.d_y, event.d_x,
						event.xppac_y[0]
					);
				} else if (event.xppac_yflag & 0x2) {
					ty[i] = DeutronRelativeApproximateTrack(
						event.t0_energy[0], event.t0_energy[1],
						event.taf_energy, using_ppac_yz[1],
						event.be_y[0], event.he_y[0], event.d_y, event.d_x,
						event.xppac_y[1]
					);					
				}
			} else if (i == 2) {
				tx[i] = generator.Gaus(0.0, 6.0);
				ty[i] = generator.Gaus(0.0, 6.0);
			} else if (i == 3) {
				tx[i] = 0.0;
				ty[i] = 0.0;
			}

			// 10Be momentum
			double be_momentum =
				MomentumFromKinetic(mass_10be, event.t0_energy[0]);
			// 10Be momentum vector
			ROOT::Math::XYZVector p_be(
				event.be_x[0] - tx[i],
				event.be_y[0] - ty[i],
				100.0
			);
			p_be = p_be.Unit() * be_momentum;

			// 4He momentum
			double he_momentum =
				MomentumFromKinetic(mass_4he, event.t0_energy[1]);
			// 4He momentum vector
			ROOT::Math::XYZVector p_he(
				event.he_x[0] - tx[i],
				event.he_y[0] - ty[i],
				100.0
			);
			p_he = p_he.Unit() * he_momentum;

			// 2H momentum
			double d_momentum =
				MomentumFromKinetic(mass_2h, event.taf_energy);
			// 2H momentum vector
			ROOT::Math::XYZVector p_d(
				event.d_x - tx[i],
				event.d_y - ty[i],
				135.0
			);
			p_d = p_d.Unit() * d_momentum;

			// beam 14C momentum vector
			ROOT::Math::XYZVector p_beam = p_be + p_he + p_d;

			// 14C momentum
			double beam_momentum = p_beam.R();
			// 14C kinematic energy
			event.c14_kinetic =
				sqrt(pow(beam_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

			// three-fold Q value
			q[i] = event.t0_energy[0] + event.t0_energy[1]
				+ event.taf_energy - event.c14_kinetic;

			// excited 14C momentum vector
			ROOT::Math::XYZVector p_excited_c = p_be + p_he;
			// excited 14C momentum
			double excited_c_momentum = p_excited_c.R();
			// excited 14C total energy
			double excited_c_energy =
				(event.t0_energy[0] + mass_10be + be10_excited_energy)
				+ (event.t0_energy[1] + mass_4he);
			// excited 14C mass
			double excited_c_mass = sqrt(
				pow(excited_c_energy, 2.0) - pow(excited_c_momentum, 2.0)
			);
			// excited energy of 14C
			ex[i] = excited_c_mass - mass_14c;
		}

		

		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	// save tree
	opf.cd();
	opt.Write();
	// close file
	opf.Close();
	ipf.Close();	
	return 0;
}