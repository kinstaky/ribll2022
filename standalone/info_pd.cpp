#include <iostream>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/channel_event.h"
#include "include/event/ta_event.h"
#include "include/event/t0_event.h"
#include "include/event/particle_event.h"

using namespace ribll;

int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%sinfo/pd.root",
		kGenerateDataPath
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "information of 15C pd reaction");
	// output data
	// energy
	double tafd_energy;
	int csi_index;
	double csi_channel;
	double t0d1_channel;
	double t0d2_channel;
	// position
	double p0x;
	double p0y;
	double p0z;
	double p1x;
	double p1y;
	double p1z;
	// state
	int state;
	// setup output branches
	opt.Branch("tafd_energy", &tafd_energy, "tafde/D");
	opt.Branch("csi_index", &csi_index, "ci/I");
	opt.Branch("csi_channel", &csi_channel, "csic/D");
	opt.Branch("t0d1_channel", &t0d1_channel, "d1c/D");
	opt.Branch("t0d2_channel", &t0d2_channel, "d2c/D");
	opt.Branch("p0x", &p0x, "p0x/D");
	opt.Branch("p0y", &p0y, "p0y/D");
	opt.Branch("p0z", &p0z, "p0z/D");
	opt.Branch("p1x", &p1x, "p1x/D");
	opt.Branch("p1y", &p1y, "p1y/D");
	opt.Branch("p1z", &p1z, "p1z/D");
	opt.Branch("state", &state, "state/I");

	for (unsigned int run = 433; run <= 442; ++run) {
		std::cout << "Processing run " << run << "\n";
		// input file name
		TString input_file_name = TString::Format(
			"%s%sC15-14C-2H-%04d.root",
			kGenerateDataPath,
			kChannelDir,
			run
		);
		// input file
		TFile input_file(input_file_name, "read");
		// input tree
		TTree *ipt = (TTree*)input_file.Get("tree");
		// input channel event
		ChannelEvent channel;
		// t0 index
		int t0_index;
		// setup input branches
		channel.SetupInput(ipt);
		ipt->SetBranchAddress("t0_index", &t0_index);

		std::vector<long long> valid_entries;
		std::vector<int> taf_indexes;
		std::vector<int> t0_indexes;
		std::vector<short> states;

		// loop to record information in channel
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);
			valid_entries.push_back(channel.entry);
			taf_indexes.push_back(channel.taf_index);
			t0_indexes.push_back(t0_index);
			double q = channel.parent_energy + 1.006
				- channel.daughter_energy[0] - channel.daughter_energy[1];
			if (q > -3.0 && q < 2.5) states.push_back(0);
			else if (q > 3.0 && q < 9.0) states.push_back(1);
			else states.push_back(-1);
		}
		// close file
		input_file.Close();

		// open telescope file to get more information
		// t0 telescope
		TString t0_file_name = TString::Format(
			"%s%st0-telescope-ta-%04u.root",
			kGenerateDataPath,
			kTelescopeDir,
			run
		);
		// t0 file
		TFile t0_file(t0_file_name, "read");
		// t0 tree
		TTree *t0_tree = (TTree*)t0_file.Get("tree");
		// t0 event
		T0Event t0;
		// setup input branches
		t0.SetupInput(t0_tree);

		// taf events
		TaEvent taf[6];
		// add taf telescopes as friends
		for (int i = 0; i < 6; ++i) {
			t0_tree->AddFriend(
				TString::Format("taf%d=tree", i),
				TString::Format(
					"%s%staf%d-telescope-ta-%04u.root",
					kGenerateDataPath,
					kTelescopeDir,
					i,
					run
				)
			);
			// setup input branches
			taf[i].SetupInput(t0_tree, TString::Format("taf%d.", i).Data());
		}

		// XPPAC events
		ParticleEvent xppac;
		// add xppac as friend
		t0_tree->AddFriend(
			"xppac=tree",
			TString::Format(
				"%s%sxppac-particle-ta-%04u.root",
				kGenerateDataPath,
				kParticleDir,
				run
			)
		);
		// setup input branches
		xppac.SetupInput(t0_tree, "xppac.");

		for (size_t i = 0; i < valid_entries.size(); ++i) {
			t0_tree->GetEntry(valid_entries[i]);
			tafd_energy = taf[taf_indexes[i]].energy[0][0];
			if (taf[taf_indexes[i]].flag[0] == 0x3) {
				csi_index = taf_indexes[i]*2;
			} else if (taf[taf_indexes[i]].flag[0] == 0x5) {
				csi_index = taf_indexes[i]*2 + 1;
			} else {
				csi_index = -1;
			}
			csi_channel = taf[taf_indexes[i]].energy[0][1];
			t0d1_channel = t0.energy[t0_indexes[i]][0];
			t0d2_channel = t0.energy[t0_indexes[i]][1];

			// p0 from T0
			ROOT::Math::XYZVector p0(
				t0.x[t0_indexes[i]][0] - xppac.x[3],
				t0.y[t0_indexes[i]][0] - xppac.y[3],
				t0.z[t0_indexes[i]][0] - xppac.z[3]
			);
			p0 = p0.Unit();
			p0x = p0.X();
			p0y = p0.Y();
			p0z = p0.Z();

			// p1 from TAF
			ROOT::Math::Polar3DVector p1p(
				taf[taf_indexes[i]].radius[0][0],
				taf[taf_indexes[i]].theta[0][0],
				taf[taf_indexes[i]].phi[0][0]
			);
			ROOT::Math::XYZVector p1(
				p1p.X() - xppac.x[3],
				p1p.Y() - xppac.y[3],
				p1p.Z() - xppac.z[3]
			);
			p1 = p1.Unit();
			p1x = p1.X();
			p1y = p1.Y();
			p1z = p1.Z();

			state = states[i];
			opt.Fill();
		}
		t0_file.Close();
	}
	opf.cd();
	// save trees
	opt.Write();
	// close files
	opf.Close();
	return 0;
}