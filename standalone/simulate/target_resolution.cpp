#include <iostream>
#include <vector>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <Math/Vector3D.h>

#include "include/event/channel_v2_event.h"
#include "include/event/generate_event.h"
#include "include/event/particle_event.h"
#include "include/ppac_track.h"

// #define ISOLATED

using namespace ribll;

int main() {
	// input file name
	TString resolution_file_name = TString::Format(
		"%s%sfull-resolution.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// channel file
	TFile resolution_file(resolution_file_name, "read");
	// input tree
	TTree *resolution_tree = (TTree*)resolution_file.Get("tree");
	if (!resolution_tree) {
		std::cerr << "Error: Get tree from "
			<< resolution_file_name << " failed.\n";
        return -1;
	}
	// input data
	double tx, ty, bex, bey, hex, hey, dx, dy;
	double bek, hek, dk, ck;
	unsigned short ppac_xflag, ppac_yflag;
	double ppac_x[3], ppac_y[3];
	double stx, sty, mptx, mpty, spatx[3], spaty[3], spftx[9], spfty[9];
	double t0bex, t0bey, t0hex, t0hey, tafdx, tafdy;
	double t0bek, t0hek, tafdk;
	// setup input branches
	resolution_tree->SetBranchAddress("tx", &tx);
    resolution_tree->SetBranchAddress("ty", &ty);
    resolution_tree->SetBranchAddress("bex", &bex);
    resolution_tree->SetBranchAddress("bey", &bey);
    resolution_tree->SetBranchAddress("hex", &hex);
    resolution_tree->SetBranchAddress("hey", &hey);
    resolution_tree->SetBranchAddress("dx", &dx);
    resolution_tree->SetBranchAddress("dy", &dy);
    resolution_tree->SetBranchAddress("bek", &bek);
    resolution_tree->SetBranchAddress("hek", &hek);
    resolution_tree->SetBranchAddress("dk", &dk);
    resolution_tree->SetBranchAddress("ck", &ck);
	resolution_tree->SetBranchAddress("ppac_xflag", &ppac_xflag);
	resolution_tree->SetBranchAddress("ppac_yflag", &ppac_yflag);
	resolution_tree->SetBranchAddress("ppac_x", ppac_x);
	resolution_tree->SetBranchAddress("ppac_y", ppac_y);
	resolution_tree->SetBranchAddress("stx", &stx);
	resolution_tree->SetBranchAddress("sty", &sty);
	resolution_tree->SetBranchAddress("mptx", &mptx);
	resolution_tree->SetBranchAddress("mpty", &mpty);
	resolution_tree->SetBranchAddress("t0bex", &t0bex);
	resolution_tree->SetBranchAddress("t0bey", &t0bey);
	resolution_tree->SetBranchAddress("t0hex", &t0hex);
	resolution_tree->SetBranchAddress("t0hey", &t0hey);
	resolution_tree->SetBranchAddress("tafdx", &tafdx);
	resolution_tree->SetBranchAddress("tafdy", &tafdy);
	resolution_tree->SetBranchAddress("t0bek", &t0bek);
	resolution_tree->SetBranchAddress("t0hek", &t0hek);
	resolution_tree->SetBranchAddress("tafdk", &tafdk);

	// fsolve file name
	TString fsolve_file_name = TString::Format(
#ifdef ISOLATED
        "%s%ssingle-ppac-generate-fsolve-sim.root",
        kGenerateDataPath,
        kSimulateDir
#else
		"%s%ssingle-ppac-channel-fsolve-sim.root",
		kGenerateDataPath,
        kChannelDir
#endif
    );
	resolution_tree->AddFriend("fsolve=tree", fsolve_file_name);
	resolution_tree->SetBranchAddress("fsolve.ftx", spftx);
	resolution_tree->SetBranchAddress("fsolve.fty", spfty);

    // output file name
    TString output_file_name = TString::Format(
        "%s%starget-resolution.root",
        kGenerateDataPath,
        kSimulateDir
    );
    // output file
    TFile opf(output_file_name, "recreate");
	// target position
	TH1F target_3ppac[2]  = {
		TH1F("h3ptx", "target X from 3 PPACs", 100, -3, 3),
		TH1F("h3pty", "target Y from 3 PPACs", 100, -3, 3)
	};
	TH1F target_2ppac[6] = {
		TH1F("h2ptx0", "target X from PPAC1 and PPAC2", 100, -3, 3),
		TH1F("h2pty0", "target Y from PPAC1 and PPAC2", 100, -3, 3),
		TH1F("h2ptx1", "target X from PPAC0 and PPAC2", 100, -3, 3),
		TH1F("h2pty1", "target Y from PPAC0 and PPAC2", 100, -3, 3),
		TH1F("h2ptx2", "target X from PPAC0 and PPAC1", 100, -3, 3),
		TH1F("h2pty2", "target Y from PPAC0 and PPAC1", 100, -3, 3)
	};
	TH1F target_1ppac_approx[6] = {
		TH1F("h1patx0", "target X from PPAC0 approx", 100, -3, 3),
		TH1F("h1paty0", "target Y from PPAC0 approx", 100, -3, 3),
		TH1F("h1patx1", "target X from PPAC1 approx", 100, -3, 3),
		TH1F("h1paty1", "target Y from PPAC1 approx", 100, -3, 3),
		TH1F("h1patx2", "target X from PPAC2 approx", 100, -3, 3),
        TH1F("h1paty2", "target Y from PPAC2 approx", 100, -3, 3)
	};
	std::vector<TH1F> target_1ppac_fsolve;
	for (size_t i = 0; i < 9; ++i) {
		target_1ppac_fsolve.emplace_back(
			TString::Format("h1pftx%ld%ld", i/3, i%3),
			TString::Format("target X from PPAC%ldX-PPAC%ldY fsolve", i/3, i%3),
#ifdef ISOLATED
			100, -0.5, 0.5
#else
			100, -3, 3
#endif
		);
		target_1ppac_fsolve.emplace_back(
			TString::Format("h1pfty%ld%ld", i/3, i%3),
			TString::Format("target Y from PPAC%ldX-PPAC%ldY fsolve", i/3, i%3),
#ifdef ISOLATED
			100, -0.5, 0.5
#else
			100, -3, 3
#endif
		);
	}

	long long entries = resolution_tree->GetEntries();
	long long entry100 = entries / 100 + 1;
	printf("Processing   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		resolution_tree->GetEntry(entry);


#ifdef ISOLATED
		for (int i = 0; i < 3; ++i) {
			if ((ppac_xflag & (1 << i)) != 0) {
				spatx[i] = DeutronRelativeApproximateTrack(
                    bek, hek, dk, ppac_xz[i],
                    bex, hex, dx, dy, ppac_x[i]
                );
			}
			if ((ppac_yflag & (1 << i)) != 0) {
				spaty[i] = DeutronRelativeApproximateTrack(
                    bek, hek, dk, ppac_yz[i],
                    bey, hey, dy, dx, ppac_y[i]
                );
			}
		}

#else
		for (int i = 0; i < 3; ++i) {
			if ((ppac_xflag & (1 << i)) != 0) {
				spatx[i] = DeutronRelativeApproximateTrack(
					t0bek, t0hek, tafdk, ppac_xz[i],
					t0bex, t0hex, tafdx, tafdy, ppac_x[i]
				);
			}
			if ((ppac_yflag & (1 << i)) != 0) {
				spaty[i] = DeutronRelativeApproximateTrack(
					t0bek, t0hek, tafdk, ppac_yz[i],
					t0bey, t0hey, tafdy, tafdx, ppac_y[i]
				);
			}
		}
#endif

		if ((ppac_xflag & 1) == 1) {
			target_1ppac_approx[0].Fill(spatx[0] - tx);
		}
		if ((ppac_xflag & 2) == 2) {
			target_1ppac_approx[2].Fill(spatx[1] - tx);
		}
		if ((ppac_xflag & 4) == 4) {
			target_1ppac_approx[4].Fill(spatx[2] - tx);
		}
		if (ppac_xflag == 3) {
			target_2ppac[4].Fill(mptx - tx);
		} else if (ppac_xflag == 5) {
			target_2ppac[2].Fill(mptx - tx);
		} else if (ppac_xflag == 6) {
			target_2ppac[0].Fill(mptx - tx);
		} else if (ppac_xflag == 7) {
			target_3ppac[0].Fill(mptx - tx);
		}

		if ((ppac_yflag & 1) == 1) {
			target_1ppac_approx[1].Fill(spaty[0] - ty);
		}
		if ((ppac_yflag & 2)== 2) {
			target_1ppac_approx[3].Fill(spaty[1] - ty);
		}
		if ((ppac_yflag & 4) == 4) {
			target_1ppac_approx[5].Fill(spaty[2] - ty);
		}
		if (ppac_yflag == 3) {
			target_2ppac[5].Fill(mpty - ty);
		} else if (ppac_yflag == 5) {
			target_2ppac[3].Fill(mpty - ty);
		} else if (ppac_yflag == 6) {
			target_2ppac[1].Fill(mpty - ty);
		} else if (ppac_yflag == 7) {
			target_3ppac[1].Fill(mpty - ty);
		}

		for (int x = 0; x < 3; ++x) {
			if ((ppac_xflag & (1 << x)) == 0) continue;
			for (int y = 0; y < 3; ++y) {
				if ((ppac_yflag & (1 << y)) == 0) continue;
				target_1ppac_fsolve[2*(x*3+y)].Fill(spftx[x*3+y] - tx);
				target_1ppac_fsolve[2*(x*3+y)+1].Fill(spfty[x*3+y] - ty);
			}
		}
	}
	printf("\b\b\b\b100%%\n");


	// fit target points
	// 3 PPACs
	for (size_t i = 0; i < 2; ++i) {
		TF1 *f1 = new TF1(TString::Format("f3p%c", "xy"[i]), "gaus", -5, 5);
		f1->SetNpx(10000);
		target_3ppac[i].Fit(f1, "QR+");
		std::cout << "3 PPAC " << "XY"[i] << ": "
			<< f1->GetParameter(1) << ", " << f1->GetParameter(2) << "\n";
	}
	// 2 PPACs
	for (size_t i = 0; i < 6; ++i) {
		TF1 *f1 = new TF1(TString::Format("f2p%c%ld", "xy"[i%2], i/2), "gaus", -5, 5);
		f1->SetNpx(10000);
		target_2ppac[i].Fit(f1, "QR+");
		std::cout << "2 PPAC " << "XY"[i%2] << i/2 << ": "
			<< f1->GetParameter(1) << ", " << f1->GetParameter(2) << "\n";
	}
	// 1 PPAC approximate
	for (size_t i = 0; i < 6; ++i) {
		TF1 *f1 = new TF1(TString::Format("f1ap%c%ld", "xy"[i%2], i/2), "gaus", -3, 3);
		f1->SetNpx(10000);
		target_1ppac_approx[i].Fit(f1, "QR+");
		std::cout << "1 PPAC " << "XY"[i%2] << i/2 << ": "
			<< f1->GetParameter(1) << ", " << f1->GetParameter(2) << "\n";
	}
	// 1 PPAC fsolve
	for (size_t i = 0; i < target_1ppac_fsolve.size(); ++i) {
		TF1 *f1 = new TF1(TString::Format("f1fp%c%ld", "xy"[i%2], i/2), "gaus", -2, 2);
		f1->SetNpx(10000);
		target_1ppac_fsolve[i].Fit(f1, "QR+");
		std::cout << "1 PPAC " << "XY"[i%2] << " X" << (i/2)/3 << "Y" << (i/2)%3
		    << ": " << f1->GetParameter(1) << ", " << f1->GetParameter(2) << "\n";
	}


	opf.cd();
	for (size_t i = 0; i < 2; ++i) target_3ppac[i].Write();
	for (size_t i = 0; i < 6; ++i) target_2ppac[i].Write();
	for (size_t i = 0; i < 6; ++i) target_1ppac_approx[i].Write();
	for (TH1F &hist : target_1ppac_fsolve) hist.Write();
	opf.Close();

	resolution_file.Close();

	return 0;
}