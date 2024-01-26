#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/generate_event.h"
#include "include/event/detect_event.h"
#include "include/ppac_track.h"

using namespace ribll;

int main() {
	// input generate file name
	TString generate_file_name = TString::Format(
		"%s%sgenerate.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// generate file
	TFile generate_file(generate_file_name, "read");
	// generate tree
	TTree *ipt = (TTree*)generate_file.Get("tree");

	// detect file name
	TString detect_file_name = TString::Format(
		"%s%sdetect.root",
		kGenerateDataPath,
		kSimulateDir
	);
	ipt->AddFriend("detect=tree", detect_file_name);

	// input data
	GenerateEvent generate;
	DetectEvent detect;
	// setup input branches
	generate.SetupInput(ipt);
	detect.SetupInput(ipt, "detect.");

	// ouptut file name
	TString track_file_name = TString::Format(
		"%s%strack.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// output file
	TFile opf(track_file_name, "recreate");
	// PPAC tracking target dx
	TH1F ppac_dx("hpdx", "PPAC tracking target dx", 500, -25, 25);
	// PPAC tracking target dy
	TH1F ppac_dy("hpdy", "PPAC tracking target dy", 500, -25, 25);
	// T0 10Be tracking target dx
	TH1F t0_dx1("htdx1", "T0 10Be tracking target dx", 100, -25, 25);
	// T0 10Be tracking target dy
	TH1F t0_dy1("htdy1", "T0 10Be tracking target dy", 100, -25, 25);
	// T0 4He tracking target dx
	TH1F t0_dx2("htdx2", "T0 4He tracking target dx", 500, -25, 25);
	// T0 4He tracking target dy
	TH1F t0_dy2("htdy2", "T0 4He tracking target dy", 500, -25, 25);
	// T0 10Be and 4He target dx
	TH1F t0_dx12("htdx12", "T0 10Be and 4He target dx", 100, -25, 25);
	// T0 10Be and 4He target dy
	TH1F t0_dy12("htdy12", "T0 10Be and 4He target dy", 100, -25, 25);
	// PPAC and T0 10Be tracking target dx
	TH1F ppac_t0_dx1(
		"hptdx1", "PPAC and T0 10Be tracking target dx", 100, -25, 25
	);
	// PPAC and T0 10Be tracking target dy
	TH1F ppac_t0_dy1(
		"hptdy1", "PPAC and T0 10Be tracking target dy", 100, -25, 25
	);
	// PPAC and T0 4He tracking target dx
	TH1F ppac_t0_dx2(
		"hptdx2", "PPAC and T0 4He tracking target dx", 100, -25, 25
	);
	// PPAC and T0 4He tracking target dy
	TH1F ppac_t0_dy2(
		"hptdy2", "PPAC and T0 4He tracking target dy", 100, -25, 25
	);
	// output tree
	TTree opt("tree", "track simulate data");
	// output data
	bool valid[2];
	double t0tx[2], t0ty[2];
	// setup output branches
	opt.Branch("vaild", valid, "valid[2]/O");
	opt.Branch("tx", &(generate.target_x), "tx/D");
	opt.Branch("ty", &(generate.target_y), "ty/D");
	opt.Branch("ppactx", &(detect.tx), "ptx/D");
	opt.Branch("ppacty", &(detect.ty), "pty/D");
	opt.Branch("t0tx", t0tx, "ttx[2]/D");
	opt.Branch("t0ty", t0ty, "tty[2]/D");

	// total entries
	long long entries = ipt->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Tracking   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		ipt->GetEntry(entry);

		if (detect.valid != 7) {
			valid[0] = valid[1] = false;
			opt.Fill();
			continue;
		}

		for (int i = 0; i < 2; ++i) {
			valid[i] = true;
			if (detect.t0_layer[i] == 0) {
				valid[i] = false;
			} else {
				// fitted k and b
				double k, b;
				int len = detect.t0_layer[i] == 1 ? 2 : 3;
				// fit T0 x
				SimpleFit(t0z, detect.t0x[i], len, k, b);
				t0tx[i] = b;
				// fit T0 y
				SimpleFit(t0z, detect.t0y[i], len, k, b);
				t0ty[i] = b;
			}
		}
		opt.Fill();

		// fill histograms
		ppac_dx.Fill(detect.tx - generate.target_x);
		ppac_dy.Fill(detect.ty - generate.target_y);
		if (valid[0]) {
			t0_dx1.Fill(t0tx[0] - generate.target_x);
			t0_dy1.Fill(t0ty[0] - generate.target_y);
			ppac_t0_dx1.Fill(t0tx[0] - detect.tx);
			ppac_t0_dy1.Fill(t0ty[0] - detect.ty);
		}
		if (valid[1]) {
			t0_dx2.Fill(t0tx[1] - generate.target_x);
			t0_dy2.Fill(t0ty[1] - generate.target_y);
			ppac_t0_dx2.Fill(t0tx[1] - detect.tx);
			ppac_t0_dy2.Fill(t0ty[1] - detect.ty);
		}
		if (valid[0] && valid[1]) {
			t0_dx12.Fill(t0tx[1] - t0tx[0]);
			t0_dy12.Fill(t0ty[1] - t0ty[0]);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms
	ppac_dx.Write();
	ppac_dy.Write();
	t0_dx1.Write();
	t0_dy1.Write();
	t0_dx2.Write();
	t0_dy2.Write();
	t0_dx12.Write();
	t0_dy12.Write();
	ppac_t0_dx1.Write();
	ppac_t0_dy1.Write();
	ppac_t0_dx2.Write();
	ppac_t0_dy2.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	generate_file.Close();
	return 0;
}