// spectrum 程序版本2，版本1太臃肿了，重新写一个
#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <Math/Vector3D.h>

#include "include/event/channel_v2_event.h"
#include "include/calculator/target_energy_calculator.h"
#include "include/calculator/csi_energy_calculator.h"
#include "include/ppac_track.h"

using namespace ribll;

constexpr double be_excited_energy[4] = {
	0.0, 3.368, 6.179, 7.542
};

constexpr double correct_tafd_thickness[6] = {
	166.0, 160.0, 150.0, 160.0, 164.0, 170.0
};

// Q value correct
constexpr double q_correct[6] = {
	0.3, -0.3, -0.4, 0.0, 0.4, 0.4
};


// straight parameters
constexpr double sa12 = 0.47;
constexpr double sb12 = -0.042;
constexpr double sa23 = 0.59;
constexpr double sb23 = -0.042;


inline double Straight(double de, double e, double a, double b) {
	return sqrt(de*e + a*de*de) + b*e;
}


int main(int argc, char **argv) {
	if (argc > 2) {
		std::cout << "Usage: " << argv[0] << " [recoil].\n"
			<< "  recoil         recoil particle, 1H or 2H (default)";
		return -1;
	}
	std::string recoil_name = "2H";
	if (argc != 1) recoil_name = std::string(argv[1]);
	if (recoil_name != "1H" && recoil_name != "2H") {
		std::cerr << "Error: Invalid recoil " << recoil_name << "\n";
		return -1;
	}

	// input file name
	TString input_file_name = TString::Format(
		"%s%sC14-10Be-4He-%s-v2.root",
		kGenerateDataPath,
		kChannelDir,
		recoil_name.c_str()
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
	ChannelV2Event channel;
	// setup input branches
	channel.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sC14-10Be-4He-%s-v3.root",
		kGenerateDataPath,
		kSpectrumDir,
		recoil_name.c_str()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// Q spectrum for single TAF
	TH1F hq_taf[6];
	for (int i = 0; i < 6; ++i) {
		hq_taf[i] = TH1F(
			TString::Format("hq%d", i),
			TString::Format("Q value for TAF%d", i),
			30, -23, -8
		);
	}
	// Q value spectrum
	TH1F hq("hq", "Q value", 70, -23, -8);
	// output tree
	TTree opt("tree", "spectrum v3");
	// output data
	int valid, ppac_flag, target_flag, taf_flag;
	// PID in correct range? index 0 for origin, 1 for calculated
	int straight[2];
	// kinetic energy with target energy lost
	double be_kinetic_target, he_kinetic_target, d_kinetic_target;
	double c_momentum_target, c_kinetic_target;
	// Q value
	double q;
	// 10Be state, under different T0D2 energy method
	int be_state;
	// excited energy with target energy lost
	double excited_energy_target;

	// setup output branches
	opt.Branch("valid", &valid, "valid/I");
	opt.Branch("ppac_flag", &ppac_flag, "pflag/I");
	opt.Branch("taf_flag", &taf_flag, "tflag/I");
	opt.Branch("target_flag", &target_flag, "tarflag/I");
	opt.Branch("straight", &straight, "straight[2]/I");
	opt.Branch("be_kinetic_target", &be_kinetic_target, "bekt/D");
	opt.Branch("he_kinetic_target", &he_kinetic_target, "hekt/D");
	opt.Branch("d_kinetic_target", &d_kinetic_target, "dkt/D");
	opt.Branch("c_momentum_target", &c_momentum_target, "cpt/D");
	opt.Branch("c_kinetic_target", &c_kinetic_target, "ckt/D");
	opt.Branch("q", &q, "q/D");
	opt.Branch("be_state", &be_state, "bes/I");
	opt.Branch("excited_energy_target", &excited_energy_target, "ext/D");

	// target energy calculator
	elc::TargetEnergyCalculator be10_target("10Be", "CD2", 9.53);
	elc::TargetEnergyCalculator he4_target("4He", "CD2", 9.53);
	elc::TargetEnergyCalculator h2_target("2H", "CD2", 9.53);
	// CsI energy calculator
	elc::CsiEnergyCalculator h2_csi("2H");


	// loop to process
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		// get data
		ipt->GetEntry(entry);

		if (channel.hole != 0) continue;
		// get flags
		valid = 0;
		taf_flag = 0;
		if (recoil_name == "2H") {
			if (!channel.tafcsi_valid && channel.tafd_front_strip >= 13) {
				taf_flag = 1;
			} else if (channel.tafcsi_valid && channel.tafd_front_strip <= 13) {
				taf_flag = 2;
			}
		} else if (recoil_name == "1H") {
			if (channel.tafcsi_valid && channel.tafd_front_strip >= 12) {
				taf_flag = 1;
			}
		}
		ppac_flag = 0;
		if (channel.ppac_xnum >= 2 && channel.ppac_ynum >= 2) {
			ppac_flag = 2;
		} else if (channel.ppac_xnum >= 1 && channel.ppac_ynum >= 1) {
			ppac_flag = 1;
		}
		if (pow(channel.tx-2.5, 2.0) + pow(channel.ty+2.0, 2.0) < 196.0) {
			target_flag = 1;
		} else {
			target_flag = 0;
		}
		straight[0] = straight[1] = 0;

		double &tx = channel.tx;
		double &ty = channel.ty;

		// 10Be direction vector
		ROOT::Math::XYZVector d_be(
			channel.fragment_x[0] - tx,
			channel.fragment_y[0] - ty,
			100.0
		);
		d_be = d_be.Unit();
		// 4He direction vector
		ROOT::Math::XYZVector d_he(
			channel.fragment_x[1] - tx,
			channel.fragment_y[1] - ty,
			100.0
		);
		d_he = d_he.Unit();
		// 2H direction vector
		ROOT::Math::XYZVector d_d(
			channel.recoil_x - tx,
			channel.recoil_y - ty,
			135.0
		);
		d_d = d_d.Unit();

		// sum up
		double be_kinetic = channel.fragment_kinetic[0];
		// consider energy loss in target
		be_kinetic_target =
			be10_target.Energy(-0.5/cos(d_be.Theta()), be_kinetic);


		// 4He kinetic energy
		double he_kinetic = channel.fragment_kinetic[1];
		// consider energy loss in target
		he_kinetic_target =
			he4_target.Energy(-0.5/cos(d_he.Theta()), he_kinetic);

		// deutron kinetic energy
		double d_kinetic = 0.0;
		if (taf_flag == 1) {
			double tafd_energy = channel.tafd_energy;
			double csi_energy = h2_csi.Energy(
				d_d.Theta(), tafd_energy, correct_tafd_thickness[channel.taf_index]
			);
			d_kinetic = tafd_energy + csi_energy;
		} else if (taf_flag == 2) {
			d_kinetic = channel.recoil_energy;
		}
		// consider energy lost in target
		d_kinetic_target =
			h2_target.Energy(-0.5/cos(d_d.Theta()), d_kinetic);

		// calculate energy
		double using_ppac_xz[3] = {ppac_xz[0], ppac_xz[1], ppac_xz[2]};
		double using_ppac_yz[3] = {ppac_yz[0], ppac_yz[1], ppac_yz[2]};
		if (channel.run >= ppac_change_run) {
			using_ppac_xz[0] = all_ppac_xz[1];
			using_ppac_yz[0] = all_ppac_yz[1];
		}
		double iter_tx, iter_ty;
		for (int i = 0; i < 2; ++i) {
			// calculate iterate reaction point
			TrackPpac(
				channel.ppac_xflag, using_ppac_xz, channel.ppac_x,
				be_kinetic, he_kinetic, d_kinetic,
				channel.fragment_x[0], channel.fragment_x[1],
				channel.recoil_x, channel.recoil_y,
				iter_tx
			);
			TrackPpac(
				channel.ppac_yflag, using_ppac_yz, channel.ppac_y,
				be_kinetic, he_kinetic, d_kinetic,
				channel.fragment_y[0], channel.fragment_y[1],
				channel.recoil_y, channel.recoil_x,
				iter_ty
			);

			// calculate again
			// 10Be direction vector
			d_be = ROOT::Math::XYZVector(
				channel.fragment_x[0] - iter_tx,
				channel.fragment_y[0] - iter_ty,
				100.0
			).Unit();
			// 4He direction vector
			d_he = ROOT::Math::XYZVector(
				channel.fragment_x[1] - iter_tx,
				channel.fragment_y[1] - iter_ty,
				100.0
			).Unit();
			// 2H direction vector
			d_d = ROOT::Math::XYZVector(
				channel.recoil_x - iter_tx,
				channel.recoil_y - iter_ty,
				135.0
			).Unit();

			// calculated recoil kinetic
			d_kinetic = channel.tafd_energy + h2_csi.Energy(
				d_d.Theta(), channel.tafd_energy,
				correct_tafd_thickness[channel.taf_index]
			);
			d_kinetic_target = 
				h2_target.Energy(-0.5/cos(d_d.Theta()), d_kinetic);
		}
		if (pow(iter_tx-2.5, 2.0) + pow(iter_ty+2.0, 2.0) < 196.0) {
			target_flag = 1;
		} else {
			target_flag = 0;
		}

		// 10Be momentum
		double be_momentum= MomentumFromKinetic(mass_10be, be_kinetic_target);
		ROOT::Math::XYZVector p_be = d_be * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic_target);
		ROOT::Math::XYZVector p_he = d_he * he_momentum;

		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic_target);
		ROOT::Math::XYZVector p_d = d_d * d_momentum;

		// beam 14C momentum vector
		ROOT::Math::XYZVector p_c = p_be + p_he + p_d;

		// 14C momentum
		c_momentum_target = p_c.R();
		// 14C kinematic energy
		c_kinetic_target =
			sqrt(pow(c_momentum_target, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

		// three body Q value
		q = be_kinetic_target + he_kinetic_target + d_kinetic_target - c_kinetic_target;

		// correct Q
		q += q_correct[channel.taf_index];
		// get 10Be state from Q value
		if (q < -11.0 && q > -13.3) be_state = 0;
		else if (q < -14.5 && q > -16.5) be_state = 1;
		else if (q < -17.3 && q > -19) be_state = 2;
		// else if (q < -19 && q > -20.5) be_state = 3;
		else be_state = -1;
		// get 10Be excited energy form state
		double be10_excited_energy = be_state >= 0
			? be_excited_energy[be_state]
			: 0.0;

		// excited 14C momentum vector
		ROOT::Math::XYZVector p_excited_c = p_be + p_he;
		// excited 14C momentum
		double excited_c_momentum = p_excited_c.R();
		// excited 14C total energy
		double excited_c_energy =
			(be_kinetic_target + mass_10be + be10_excited_energy)
			+ (he_kinetic_target + mass_4he);
		// excited 14C mass
		double excited_c_mass = sqrt(
			pow(excited_c_energy, 2.0) - pow(excited_c_momentum, 2.0)
		);
		// excited energy of 14C
		excited_energy_target = excited_c_mass - mass_14c;

		if (ppac_flag == 0) valid |= 1;
		if (taf_flag == 0) valid |= 2;
		if (!channel.t0_valid) valid |= 4;
		if (target_flag == 0) valid |= 8;
		if (c_kinetic_target < 360.0) valid |= 32;

		if (valid == 0 && taf_flag == 1) {
			// fill Q spectrums
			hq_taf[channel.taf_index].Fill(q - q_correct[channel.taf_index]);
			hq.Fill(q);
		}

		// fill to tree
		opt.Fill();
	}


	constexpr double q_init_param[12] = {
		10, -10, 1,
		40, -12, 1,
		50, -15.5, 1,
		10, -18, 1
	};
	// fit TAF Q spectrum
	TF1 *fq = new TF1("fq", "gaus(0)+gaus(3)+gaus(6)+gaus(9)", -23, -8);
	fq->SetParameters(q_init_param);
	fq->SetParLimits(4, -14, -11);
	fq->SetParLimits(7, -17, -14);
	fq->SetParLimits(10, -19, -17);
	hq.Fit(fq, "QR+");
	std::cout << fq->GetParameter(1) << ", "
		<< fq->GetParameter(4) << ", "
		<< fq->GetParameter(7) << ", "
		<< fq->GetParameter(10) << "\n";
	for (size_t i = 0; i < 4; ++i) {
		TF1 *f1 = new TF1(TString::Format("fq%ld", i), "gaus", -23, -8);
		if (i == 0) f1->SetLineColor(kOrange);
		else f1->SetLineColor(kBlue);
		f1->SetParameters(fq->GetParameters()+3*i);
		hq.GetListOfFunctions()->Add(f1);
	}

	// save
	opf.cd();
	// fill spectrums
	for (size_t i = 0; i < 6; ++i) hq_taf[i].Write();
	hq.Write();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}