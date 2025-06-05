#include <iostream>
#include <vector>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <Math/Vector3D.h>

#include "include/event/channel_v2_event.h"
#include "include/event/generate_event.h"
#include "include/calculator/target_energy_calculator.h"
#include "include/event/particle_event.h"
#include "include/ppac_track.h"

// #define ISOLATED
// #define BE_EXCITED_MASS
// #define DEPTH
// #define DEPTH_RAW
// #define THREE_PPAC
// #define TWO_PPAC
// #define SINGLE_PPAC_APPROX
// #define SINGLE_PPAC_FSOLVE
// #define GAMMA_DECAY
// #define BE_POS
// #define HE_POS
// #define H_POS
// #define BE_ENERGY
// #define HE_ENERGY
// #define H_ENERGY
#define FULL

constexpr int PPAC_INDEX = 0;


#if defined(BE_EXCITED_MASS)
	constexpr double bins = 500;
	constexpr double border = 0.5;
#elif defined(DEPTH)
	constexpr double bins = 900;
	constexpr double border = 3;
#elif defined(DEPTH_RAW)
	constexpr double bins = 900;
	constexpr double border = 3;
#elif defined(GAMMA_DECAY)
	constexpr double bins = 500;
	constexpr double border = 0.5;
#elif defined(FULL)
	constexpr double bins = 500;
	constexpr double border = 5;
#else
	constexpr double bins = 400;
	constexpr double border = 2;
#endif

#if defined(BE_EXCITED_MASS)
	constexpr double exbins = 400;
	constexpr double exborder = 0.2;
#elif defined(DEPTH)
	constexpr double exbins = 400;
	constexpr double exborder = 0.4;
#elif defined(SINGLE_PPAC_APPROX) && defined(ISOLATED)
	constexpr double exbins = 1000;
	constexpr double exborder = 0.05;
#elif defined(SINGLE_PPAC_FSOLVE) && defined(ISOLATED)
	constexpr double exbins = 1000;
	constexpr double exborder = 0.05;
#else
	constexpr double exbins = 1000;
	constexpr double exborder = 1;
#endif

using namespace ribll;

int main() {
	// input file name
	TString channel_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v2-sim.root",
		kGenerateDataPath,
		kChannelDir
	);
	// channel file
	TFile channel_file(channel_file_name, "read");
	// input tree
	TTree *channel_tree = (TTree*)channel_file.Get("tree");
	if (!channel_tree) {
		std::cerr << "Error: Get tree from "
			<< channel_file_name << " failed.\n";
        return -1;
	}
	// setup input branches
	ChannelV2Event channel;
    channel.SetupInput(channel_tree);

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
	TFile fsolve_file(fsolve_file_name, "read");
	TTree *fsolve_tree = (TTree*)fsolve_file.Get("tree");
	if (!fsolve_tree) {
		std::cerr << "Error: Get tree from "
			<< fsolve_file_name << "\n";
		return -1;
	}
	double ftx[9], fty[9];
	fsolve_tree->SetBranchAddress("ftx", ftx);
	fsolve_tree->SetBranchAddress("fty", fty);


	// generate file name
	TString generate_file_name = TString::Format(
        "%s%sgenerate-0002.root",
        kGenerateDataPath,
        kSimulateDir
    );
    // generate file
    TFile generate_file(generate_file_name, "read");
    // input tree
    TTree *generate_tree = (TTree*)generate_file.Get("tree");
	if (!generate_tree) {
		std::cerr << "Error: Get tree from "
			<< generate_file_name << " failed.\n";
		return -1;
	}
	GenerateEvent generate;
	// setup input branches
	generate.SetupInput(generate_tree);

	// PPAC file name
	TString ppac_file_name = TString::Format(
		"%s%sxppac-particle-sim-ta-0002.root",
		kGenerateDataPath,
		kParticleDir
	);
	// add friend
	generate_tree->AddFriend("ppac=tree", ppac_file_name);
	// PPAC event
	ParticleEvent ppac;
	unsigned short ppac_xflag, ppac_yflag;
	// setup input branches
	ppac.SetupInput(generate_tree, "ppac.");
	generate_tree->SetBranchAddress("ppac.xflag", &ppac_xflag);
	generate_tree->SetBranchAddress("ppac.yflag", &ppac_yflag);


    // output file name
    TString output_file_name = TString::Format(
        "%s%sfull-resolution.root",
        kGenerateDataPath,
        kSimulateDir
    );
    // output file
    TFile opf(output_file_name, "recreate");
	// Q difference
	TH1F hdq[3] = {
		TH1F("hdq0", "Q difference of g.s.", bins, -border, border),
		TH1F("hdq1", "Q difference of 2_1^+", bins, -border, border),
		TH1F("hdq2", "Q difference of 6MeV", bins, -border, border)
	};
	// excitation energy difference
	TH1F hdex[3] = {
		TH1F(
			"hdex0", "excitation energy difference of g.s.",
			exbins, -exborder, exborder
		),
		TH1F(
			"hdex1", "excitation energy difference of 2_1^+",
			exbins, -exborder, exborder
		),
		TH1F(
			"hdex2", "excitation energy difference of 6MeV",
			exbins, -exborder, exborder
		)
	};
	// TH1F hdmex[3] = {
	// 	TH1F(
	// 		"hdmex0", "22MeV excitation energy difference of g.s.",
	// 		exbins, -exborder, exborder
	// 	),
	// 	TH1F(
	// 		"hdmex1", "22MeV excitation energy difference of 2_1^+",
	// 		exbins, -exborder, exborder
	// 	),
	// 	TH1F(
	// 		"hdmex2", "22MeV excitation energy difference of 6MeV",
	// 		exbins, -exborder, exborder
	// 	)
	// };
	TH2F hdex2[3] = {
		TH2F(
			"hdexe0", "excitation energy difference of g.s. VS ex energy",
			100, 12, 32,
			exbins, -exborder, exborder
		),
		TH2F(
			"hdexe1", "excitation energy difference of 2_1^+ VS ex energy",
			100, 15, 35,
			exbins, -exborder, exborder
		),
		TH2F(
			"hdexe2", "excitation energy difference of 6MeV VS ex energy",
			100, 18, 38,
			exbins, -exborder, exborder
		)
	};
	// output tree
    TTree opt("tree", "tree");
    // ouutput data
	double tx, ty, bex, bey, hex, hey, dx, dy;
	double bek, hek, dk, ck;
	double ppac_x[3], ppac_y[3];
	double stx, sty, mptx, mpty, spatx[3], spaty[3], spftx[3], spfty[3];
	double t0bex, t0bey, t0hex, t0hey, tafdx, tafdy;
	double t0bek, t0hek, tafdk;
	double q, sq, xce, sxce;
	// setup output branches
	opt.Branch("tx", &tx, "tx/D");
	opt.Branch("ty", &ty, "ty/D");
    opt.Branch("bex", &bex, "bex/D");
    opt.Branch("bey", &bey, "bey/D");
    opt.Branch("hex", &hex, "hex/D");
    opt.Branch("hey", &hey, "hey/D");
    opt.Branch("dx", &dx, "dx/D");
    opt.Branch("dy", &dy, "dy/D");
    opt.Branch("bek", &bek, "bek/D");
    opt.Branch("hek", &hek, "hek/D");
    opt.Branch("dk", &dk, "dk/D");
	opt.Branch("ck", &ck, "ck/D");
	opt.Branch("ppac_xflag", &ppac_xflag, "pxflag/s");
	opt.Branch("ppac_yflag", &ppac_yflag, "pyflag/s");
	opt.Branch("ppac_x", ppac_x, "px[3]/D");
	opt.Branch("ppac_y", ppac_y, "py[3]/D");
	opt.Branch("stx", &stx, "stx/D");
	opt.Branch("sty", &sty, "sty/D");
    opt.Branch("mptx", &mptx, "mptx/D");
    opt.Branch("mpty", &mpty, "mpty/D");
    opt.Branch("spatx", spatx, "spatx[3]/D");
    opt.Branch("spaty", spaty, "spaty[3]/D");
	opt.Branch("spftx", spftx, "spftx[3]/D");
	opt.Branch("spfty", spfty, "spfty[3]/D");
	opt.Branch("t0bex", &t0bex, "t0bex/D");
    opt.Branch("t0bey", &t0bey, "t0bey/D");
    opt.Branch("t0hex", &t0hex, "t0hex/D");
    opt.Branch("t0hey", &t0hey, "t0hey/D");
    opt.Branch("tafdx", &tafdx, "tafdx/D");
    opt.Branch("tafdy", &tafdy, "tafdy/D");
    opt.Branch("t0bek", &t0bek, "t0bek/D");
    opt.Branch("t0hek", &t0hek, "t0hek/D");
    opt.Branch("tafdk", &tafdk, "tafdk/D");
    opt.Branch("q", &q, "q/D");
    opt.Branch("sq", &sq, "sq/D");
	opt.Branch("xce", &xce, "xce/D");
	opt.Branch("sxce", &sxce, "sxce/D");

	elc::TargetEnergyCalculator be10_target("10Be", "CD2", 9.53);
	elc::TargetEnergyCalculator he4_target("4He", "CD2", 9.53);
	elc::TargetEnergyCalculator h2_target("2H", "CD2", 9.53);


	long long entries = channel_tree->GetEntries();
	long long entry100 = entries / 100 + 1;
	long long fsolve_entry = 0;
	printf("Processing   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < channel_tree->GetEntriesFast(); ++entry) {
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		channel_tree->GetEntry(entry);
		if (channel.valid != 0) continue;
		generate_tree->GetEntry(channel.entry);
#ifdef SINGLE_PPAC_FSOLVE
		fsolve_tree->GetEntry(fsolve_entry);
#endif

		tx = generate.target_x;
		ty = generate.target_y;
		bex = generate.excited_fragment0_x;
		bey = generate.excited_fragment0_y;
		hex = generate.fragment_x[1];
		hey = generate.fragment_y[1];
		dx = generate.rx;
		dy = generate.ry;
		bek = generate.excited_fragment0_kinetic_in_target;
		hek = generate.fragment_kinetic_in_target[1];
		dk = generate.recoil_kinetic_in_target;

		// calculate Q
		ROOT::Math::XYZVector be_direction(bex-tx, bey-ty, 100.0);
		be_direction = be_direction.Unit();
		double mass_excited_10be = mass_10be + generate.fragment_excited_energy;
		ROOT::Math::XYZVector bep = be_direction * MomentumFromKinetic(
			mass_excited_10be, bek
		);

		ROOT::Math::XYZVector he_direction(hex-tx, hey-ty, 100.0);
		he_direction = he_direction.Unit();
		ROOT::Math::XYZVector hep = he_direction * MomentumFromKinetic(mass_4he, hek);

		ROOT::Math::XYZVector d_direction(dx-tx, dy-ty, 135.0);
		d_direction = d_direction.Unit();
		ROOT::Math::XYZVector dp = d_direction * MomentumFromKinetic(mass_2h, dk);

		ROOT::Math::XYZVector cp = bep + hep + dp;
		ck = KineticFromMomentum(mass_14c, cp.R());
		q = bek + hek + dk - ck;

		for (size_t i = 0; i < 3; ++i) {
			ppac_x[i] = ppac.x[i];
			ppac_y[i] = ppac.y[i];
		}

		// calculate excited energy
		ROOT::Math::XYZVector xcp = bep + hep;
		double xc_energy = bek + mass_excited_10be + hek + mass_4he;
		double xc_mass = sqrt(pow(xc_energy, 2.0) - xcp.Mag2());
		xce = xc_mass - mass_14c;
		xce = generate.beam_excited_energy;

		stx = channel.tx;
		sty = channel.ty;
		t0bex = channel.fragment_x[0];
		t0bey = channel.fragment_y[0];
		t0hex = channel.fragment_x[1];
		t0hey = channel.fragment_y[1];
		tafdx = channel.recoil_x;
		tafdy = channel.recoil_y;
		t0bek = channel.fragment_kinetic[0];
		t0hek = channel.fragment_kinetic[1];
		tafdk = channel.recoil_kinetic;
		sq = q;
		double tmp;
		if (ppac_xflag >= 3 && ppac_xflag != 4) {
			TrackMultiplePpac(ppac_xflag, ppac_xz, ppac.x, tmp, mptx);
		}
		if (ppac_yflag >= 3 && ppac_yflag != 4) {
			TrackMultiplePpac(ppac_yflag, ppac_yz, ppac.y, tmp, mpty);
		}
		for (int i = 0; i < 3; ++i) {
#ifdef ISOLATED
			if ((ppac_xflag & (1 << i)) != 0) {
				spatx[i] = DeutronRelativeApproximateTrack(
                    bek, hek, dk, ppac_xz[i],
                    bex, hex, dx, dy, ppac.x[i]
                );
			}
			if ((ppac_yflag & (1 << i)) != 0) {
				spaty[i] = DeutronRelativeApproximateTrack(
                    bek, hek, dk, ppac_yz[i],
                    bey, hey, dy, dx, ppac.y[i]
                );
			}
#else
			if ((ppac_xflag & (1 << i)) != 0) {
				spatx[i] = DeutronRelativeApproximateTrack(
					t0bek, t0hek, tafdk, ppac_xz[i],
					t0bex, t0hex, tafdx, tafdy, ppac.x[i]
				);
			}
			if ((ppac_yflag & (1 << i)) != 0) {
				spaty[i] = DeutronRelativeApproximateTrack(
					t0bek, t0hek, tafdk, ppac_yz[i],
					t0bey, t0hey, tafdy, tafdx, ppac.y[i]
				);
			}
#endif
		}
		for (int i = 0; i < 3; ++i) {
			spftx[i] = ftx[i*3+i];
			spfty[i] = fty[i*3+i];
		}
		double sbek = bek;
		double shek = hek;
		double sdk = dk;

#if defined(BE_EXCITED_MASS)
		bep = be_direction * MomentumFromKinetic(mass_10be, bek);
#elif defined(DEPTH)
		sbek = be10_target.Energy(
			-0.5 / cos(be_direction.Theta()),
			be10_target.Energy(
				(1.0 - generate.depth) / cos(generate.excited_fragment0_theta),
				generate.excited_fragment0_kinetic_in_target
			)
		);
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, sbek);
		shek = he4_target.Energy(
			-0.5 / cos(he_direction.Theta()),
			generate.fragment_kinetic_after_target[1]
		);
		hep = he_direction * MomentumFromKinetic(mass_4he, shek);
		sdk = h2_target.Energy(
			-0.5 / cos(he_direction.Theta()),
			generate.recoil_kinetic_after_target
		);
		dp = d_direction * MomentumFromKinetic(mass_2h, sdk);
#elif defined(DEPTH_RAW)
		sbek = be10_target.Energy(
			(1.0 - generate.depth) / cos(generate.excited_fragment0_theta),
			generate.excited_fragment0_kinetic_in_target
		);
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, sbek);
		shek = generate.fragment_kinetic_after_target[1];
		hep = he_direction * MomentumFromKinetic(mass_4he, shek);
		sdk = generate.recoil_kinetic_after_target;
		dp = d_direction * MomentumFromKinetic(mass_2h, sdk);
#elif defined(THREE_PPAC)
		if (ppac_xflag != 7) continue;
		if (ppac_yflag != 7) continue;
		be_direction = ROOT::Math::XYZVector(
			bex - mptx, bey - mpty, 100.0
		).Unit();
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, bek);
		he_direction = ROOT::Math::XYZVector(
			hex - mptx, hey - mpty, 100.0
		).Unit();
		hep = he_direction * MomentumFromKinetic(mass_4he, hek);
		d_direction = ROOT::Math::XYZVector(
            dx - mptx, dy - mpty, 135.0
        ).Unit();
		dp = d_direction * MomentumFromKinetic(mass_2h, dk);
#elif defined(TWO_PPAC)
		int flag = (~(1 << PPAC_INDEX)) & 0x7;
		if ((ppac_xflag & flag) != flag) continue;
		if ((ppac_yflag & flag) != flag) continue;
		be_direction = ROOT::Math::XYZVector(
			bex - mptx, bey - mpty, 100.0
		).Unit();
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, bek);
		he_direction = ROOT::Math::XYZVector(
			hex - mptx, hey - mpty, 100.0
		).Unit();
		hep = he_direction * MomentumFromKinetic(mass_4he, hek);
		d_direction = ROOT::Math::XYZVector(
            dx - mptx, dy - mpty, 135.0
        ).Unit();
		dp = d_direction * MomentumFromKinetic(mass_2h, dk);

#elif defined(SINGLE_PPAC_APPROX)
	#ifdef ISOLATED
		if ((ppac_xflag & (1 << PPAC_INDEX)) == 0) continue;
		if ((ppac_yflag & (1 << PPAC_INDEX)) == 0) continue;
		be_direction = ROOT::Math::XYZVector(
			bex - spatx[PPAC_INDEX], bey - spaty[PPAC_INDEX], 100.0
		).Unit();
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, bek);
		he_direction = ROOT::Math::XYZVector(
			hex - spatx[PPAC_INDEX], hey - spaty[PPAC_INDEX], 100.0
		).Unit();
		hep = he_direction * MomentumFromKinetic(mass_4he, hek);
		d_direction = ROOT::Math::XYZVector(
            dx - spatx[PPAC_INDEX], dy - spaty[PPAC_INDEX], 135.0
        ).Unit();
		dp = d_direction * MomentumFromKinetic(mass_2h, dk);
	#else
		if ((ppac_xflag & (1 << PPAC_INDEX)) == 0) continue;
		if ((ppac_yflag & (1 << PPAC_INDEX)) == 0) continue;
		be_direction = ROOT::Math::XYZVector(
			t0bex - spatx[PPAC_INDEX], t0bey - spaty[PPAC_INDEX], 100.0
		).Unit();
		bep = be_direction * MomentumFromKinetic(mass_10be, t0bek);
		he_direction = ROOT::Math::XYZVector(
			t0hex - spatx[PPAC_INDEX], t0hey - spaty[PPAC_INDEX], 100.0
		).Unit();
		hep = he_direction * MomentumFromKinetic(mass_4he, t0hek);
		d_direction = ROOT::Math::XYZVector(
			tafdx - spatx[PPAC_INDEX], tafdy - spaty[PPAC_INDEX], 135.0
		).Unit();
		dp = d_direction * MomentumFromKinetic(mass_2h, tafdk);
		sbek = t0bek;
		shek = t0hek;
		sdk = tafdk;
	#endif

#elif defined(SINGLE_PPAC_FSOLVE)
	#ifdef ISOLATED
		if ((ppac_xflag & (1 << PPAC_INDEX)) == 0) {
			++fsolve_entry;
			continue;
		}
		if ((ppac_yflag & (1 << PPAC_INDEX)) == 0) {
			++fsolve_entry;
			continue;
		}
		be_direction = ROOT::Math::XYZVector(
			bex - spftx[PPAC_INDEX], bey - spfty[PPAC_INDEX], 100.0
		).Unit();
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, bek);
		he_direction = ROOT::Math::XYZVector(
			hex - spftx[PPAC_INDEX], hey - spfty[PPAC_INDEX], 100.0
		).Unit();
		hep = he_direction * MomentumFromKinetic(mass_4he, hek);
		d_direction = ROOT::Math::XYZVector(
            dx - spftx[PPAC_INDEX], dy - spfty[PPAC_INDEX], 135.0
        ).Unit();
		dp = d_direction * MomentumFromKinetic(mass_2h, dk);
	#else
		if ((ppac_xflag & (1 << PPAC_INDEX)) == 0) {
			++fsolve_entry;
			continue;
		}
		if ((ppac_yflag & (1 << PPAC_INDEX)) == 0) {
			++fsolve_entry;
			continue;
		}
		be_direction = ROOT::Math::XYZVector(
			t0bex - spatx[PPAC_INDEX], t0bey - spaty[PPAC_INDEX], 100.0
		).Unit();
		bep = be_direction * MomentumFromKinetic(mass_10be, t0bek);
		he_direction = ROOT::Math::XYZVector(
			t0hex - spatx[PPAC_INDEX], t0hey - spaty[PPAC_INDEX], 100.0
		).Unit();
		hep = he_direction * MomentumFromKinetic(mass_4he, t0hek);
		d_direction = ROOT::Math::XYZVector(
			tafdx - spatx[PPAC_INDEX], tafdy - spaty[PPAC_INDEX], 135.0
		).Unit();
		dp = d_direction * MomentumFromKinetic(mass_2h, tafdk);
		sbek = t0bek;
		shek = t0hek;
		sdk = tafdk;
	#endif

#elif defined(GAMMA_DECAY)
		sbek = generate.fragment_kinetic_in_target[0];
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, sbek);

#elif defined(BE_POS)
		be_direction = ROOT::Math::XYZVector(
			channel.fragment_x[0] - tx,
			channel.fragment_y[0] - ty,
			100.0
		).Unit();
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, bek);

#elif defined(HE_POS)
		he_direction = ROOT::Math::XYZVector(
			channel.fragment_x[1] - tx,
			channel.fragment_y[1] - ty,
			100.0
		).Unit();
		hep = he_direction * MomentumFromKinetic(mass_4he, hek);

#elif defined(H_POS)
		d_direction = ROOT::Math::XYZVector(
			channel.recoil_x - tx,
			channel.recoil_y - ty,
			135.0
		).Unit();
		dp = d_direction * MomentumFromKinetic(mass_2h, dk);


#elif defined(BE_ENERGY)
		sbek = t0bek;
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, sbek);

#elif defined(HE_ENERGY)
		shek = t0hek;
		hep = he_direction * MomentumFromKinetic(mass_4he, shek);

#elif defined(H_ENERGY)
		sdk = tafdk;
		dp = d_direction * MomentumFromKinetic(mass_2h, sdk);

#elif defined(FULL)
		sbek = t0bek;
		shek = t0hek;
		sdk = tafdk;
		be_direction = ROOT::Math::XYZVector(
			t0bex - stx, t0bey - sty, 100.0
		).Unit();
		// bep = be_direction * MomentumFromKinetic(mass_10be, sbek);
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, sbek);
		he_direction = ROOT::Math::XYZVector(
            t0hex - stx, t0hey - sty, 100.0
        ).Unit();
		hep = he_direction * MomentumFromKinetic(mass_4he, shek);
		d_direction = ROOT::Math::XYZVector(
            tafdx - stx, tafdy - sty, 135.0
        ).Unit();
		dp = d_direction * MomentumFromKinetic(mass_2h, sdk);
#endif
		cp = bep + hep + dp;
		double sck = KineticFromMomentum(mass_14c, cp.R());
		sq = sbek + shek + sdk - sck;
#if defined(SINGLE_PPAC_APPROX) || defined(SINGLE_PPAC_FSOLVE)
		bep = be_direction * MomentumFromKinetic(mass_excited_10be, t0bek);
#endif
		xcp = bep + hep;
		double sxc_energy = sbek + mass_excited_10be + shek + mass_4he;
		// double sxc_energy = sbek + mass_10be + shek + mass_4he;
		double sxc_mass = sqrt(pow(sxc_energy, 2.0) - xcp.Mag2());
		sxce = sxc_mass - mass_14c;


		hdq[generate.fragment_state].Fill(sq-q);
		hdex[generate.fragment_state].Fill(sxce-xce);
		// if (xce > 21.5 && xce < 22.5) {
		// 	hdmex[generate.fragment_state].Fill(sxce-xce);
		// }
		hdex2[generate.fragment_state].Fill(xce, sxce-xce);
		opt.Fill();
		++fsolve_entry;
	}
	printf("\b\b\b\b100%%\n");

	// fit
	for (int i = 0; i < 3; ++i) {
		TF1 *f1 = new TF1(TString::Format("f%d", i), "gaus", -border, border);
		f1->SetNpx(10000);
		hdq[i].Fit(f1, "QR+");
		std::cout << "dQ" << i << ": "
			<< f1->GetParameter(1) << ", " << f1->GetParameter(2) << "\n";
	}

	for (int i = 0; i < 3; ++i) {
		TF1 *f1 = new TF1(TString::Format("fx%d", i), "gaus", -exborder, exborder);
		f1->SetNpx(10000);
		hdex[i].Fit(f1, "QR+");
		std::cout << "dEx " << i << ": "
			<< f1->GetParameter(1)*1e3 << ", " << f1->GetParameter(2)*1e3 << "\n";
	}

	// for (int i = 0; i < 3; ++i) {
	// 	TF1 *f1 = new TF1(TString::Format("fmx%d", i), "gaus", -exborder, exborder);
	// 	f1->SetNpx(10000);
	// 	hdmex[i].Fit(f1, "QR+");
	// 	std::cout << "22 MeV dEx " << i << ": "
	// 		<< f1->GetParameter(1)*1e3 << ", " << f1->GetParameter(2)*1e3 << "\n";
	// }

	opf.cd();
	for (size_t i = 0; i < 3; ++i) hdq[i].Write();
	for (size_t i = 0; i < 3; ++i) hdex[i].Write();
	// for (size_t i = 0; i < 3; ++i) hdmex[i].Write();
	for (size_t i = 0; i < 3; ++i) hdex2[i].Write();
	opt.Write();
	opf.Close();

	generate_file.Close();
	channel_file.Close();

	return 0;
}