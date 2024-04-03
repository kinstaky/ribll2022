#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <THStack.h>
#include <TChain.h>
#include <Math/Vector3D.h>

#include "include/event/t0_event.h"
#include "include/event/particle_event.h"
#include "include/ppac_track.h"

using namespace ribll;

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " run end_run\n";
        return -1;
    }
    int run = atoi(argv[1]);
    int end_run = atoi(argv[2]);

    // input t0 chain
	TChain t0_chain("t0", "t0 telescope");
	for (int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		t0_chain.AddFile(TString::Format(
			"%s%st0-telescope-ta-%04d.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			i
		));
	}
	// input event
	T0Event t0_event;
	// setup input branches
	t0_event.SetupInput(&t0_chain);

    // input XPPA chain
    TChain ppac_chain("ppac", "PPAC particle");
    for (int i = run; i <= end_run; ++i) {
        if (i == 628) continue;
        ppac_chain.AddFile(TString::Format(
            "%s%sxppac-particle-ta-%04d.root/tree",
            kGenerateDataPath,
            kParticleDir,
            i
        ));
    }
    // add friend
    t0_chain.AddFriend(&ppac_chain, "ppac");
    // input event
    ParticleEvent ppac_event;
    unsigned short xflag, yflag;
    // setup input branches
    ppac_event.SetupInput(&t0_chain, "ppac.");
    t0_chain.SetBranchAddress("ppac.xflag", &xflag);
    t0_chain.SetBranchAddress("ppac.yflag", &yflag);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-pid-angle-ta-%04d-%04d.root",
		kGenerateDataPath,
		kShowDir,
		run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
    // 1D histogram of 10Be angle
    TH1F hist_be_angle("hbea", "10Be cos(#theta)", 100, 0.9, 1);
    // 1D histogram of 4He angle
    TH1F hist_he_angle("hhea", "4He cos(#theta)", 100, 0.9, 1);
    // 4He PID in D1D2 slice and cos(theta) in [0.97, 0.98], [0.98, 0.99], [0.99, 1.0]
    TH2F pid_d1d2_he[3] = {
        TH2F("hepid1", "4He PID in D1D2 slice", 1000, 0, 30000, 1000, 0, 30000),
        TH2F("hepid2", "4He PID in D1D2 slice", 1000, 0, 30000, 1000, 0, 30000),
        TH2F("hepid3", "4He PID in D1D2 slice", 1000, 0, 30000, 1000, 0, 30000)
    };
    // set color
    pid_d1d2_he[0].SetMarkerColor(kBlack);
    pid_d1d2_he[1].SetMarkerColor(kRed);
    pid_d1d2_he[2].SetMarkerColor(kBlue);
    THStack stack_pid_d1d2_he("hepid", "4He PID in D1D2 slice");
    stack_pid_d1d2_he.Add(pid_d1d2_he);
    stack_pid_d1d2_he.Add(pid_d1d2_he+1);
    stack_pid_d1d2_he.Add(pid_d1d2_he+2);


    // total number of entries
	long long entries = t0_chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling T0 pid   0%%");
	fflush(stdout);
	// loop to fill pid histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		t0_chain.GetEntry(entry);

        // corrected PPAC position
        double ppac_x[3] = {
            ppac_event.x[0]-ppac_correct[0][0],
            ppac_event.x[1]-ppac_correct[0][1],
            ppac_event.x[2]-ppac_correct[0][2]
        };
        double ppac_y[3] = {
            ppac_event.y[0]-ppac_correct[1][0],
            ppac_event.y[1]-ppac_correct[1][1],
            ppac_event.y[2]-ppac_correct[1][2]
        };
        // fitted parameters
        double xk, yk, tx, ty;
        // track PPAC and get target X and Y
        TrackMultiplePpac(xflag, ppac_xz, ppac_x, xk, tx);
        TrackMultiplePpac(yflag, ppac_yz, ppac_y, yk, ty);

        // loop T0 particles
        for (int i = 0; i < t0_event.num; ++i) {
            ROOT::Math::XYZVector direction(
                t0_event.x[i][0] - tx,
                t0_event.y[i][0] - ty,
                100.0
            );
            if (t0_event.charge[i] == 2 && t0_event.mass[i] == 4) {
                double angle = cos(direction.Theta());
                hist_he_angle.Fill(angle);
                if (angle >= 0.97 && angle < 1.0 && (t0_event.flag[i]&0x3) == 0x3 ) {
                    int index = int((angle - 0.97) * 100.0);
                    pid_d1d2_he[index].Fill(t0_event.energy[i][1], t0_event.energy[i][0]);
                }
            } else if (t0_event.charge[i] == 4 && t0_event.mass[i] == 10) {
                hist_be_angle.Fill(cos(direction.Theta()));
            }
        }
	}
	// show finish
	printf("\b\b\b\b100%%\n");

    // save histograms
    hist_be_angle.Write();
    hist_he_angle.Write();
    for (auto &hist : pid_d1d2_he) hist.Write();
    stack_pid_d1d2_he.Write();
    // close file
    opf.Close();
    return 0;
}