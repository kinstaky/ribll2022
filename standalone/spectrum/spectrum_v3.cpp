// spectrum 程序版本3，仅从 channel_v2 读取
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
#include "include/ppac_track.h"

using namespace ribll;

constexpr double be_excited_energy[4] = {
	0.0, 3.368, 6.179, 7.542
};

constexpr double correct_tafd_thickness[6] = {
	166.0, 160.0, 150.0, 160.0, 164.0, 170.0
};

// Q value correct
constexpr double q_tafd_correct[6] = {
	0.3, -0.3, -0.4, 0.0, 0.4, 0.4
};
constexpr double q_csi_correct[12] = {
	-0.46, -0.25, 0.18, 0.49,
	0.49, 0.55, 0.28, 0.30,
	-0.04, -0.10, -0.12, -0.64
};

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
	// kinetic energy with target energy lost
	double be_kinetic, he_kinetic, d_kinetic;
	double c_momentum, c_kinetic;
	// particle direction
	double be_dx, be_dy, be_dz;
	double he_dx, he_dy, he_dz;
	double d_dx, d_dy, d_dz;
	// Q value
	double q;
	// 10Be state, under different T0D2 energy method
	int be_state;
	// excited energy with target energy lost
	double excited_energy;
	// image 1H part
	double image_q;

	// setup output branches
	opt.Branch("valid", &valid, "valid/I");
	opt.Branch("ppac_flag", &ppac_flag, "pflag/I");
	opt.Branch("taf_flag", &taf_flag, "tflag/I");
	opt.Branch("target_flag", &target_flag, "tarflag/I");
	opt.Branch("be_kinetic", &be_kinetic, "bek/D");
	opt.Branch("he_kinetic", &he_kinetic, "hek/D");
	opt.Branch("d_kinetic", &d_kinetic, "dk/D");
	opt.Branch("c_momentum", &c_momentum, "cp/D");
	opt.Branch("c_kinetic", &c_kinetic, "ck/D");
	opt.Branch("be_dx", &be_dx, "bedx/D");
	opt.Branch("be_dy", &be_dy, "bedy/D");
	opt.Branch("be_dz", &be_dz, "bedz/D");
	opt.Branch("he_dx", &he_dx, "hedx/D");
	opt.Branch("he_dy", &he_dy, "hedy/D");
	opt.Branch("he_dz", &he_dz, "hedz/D");
	opt.Branch("d_dx", &d_dx, "ddx/D");
	opt.Branch("d_dy", &d_dy, "ddy/D");
	opt.Branch("d_dz", &d_dz, "ddz/D");
	opt.Branch("q", &q, "q/D");
	opt.Branch("be_state", &be_state, "bes/I");
	opt.Branch("excited_energy", &excited_energy, "ex/D");
	opt.Branch("image_q", &image_q, "imq/D");

	// loop to process
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		// get data
		ipt->GetEntry(entry);

		if (channel.hole != 0) continue;
		// get flags
		valid = 0;
		if (!channel.t0_valid) valid |= 1;
		taf_flag = 0;
		if (recoil_name == "2H") {
			if (!channel.tafcsi_valid && channel.tafd_front_strip >= 13) {
				taf_flag = 1;
			} else if (channel.tafcsi_valid && channel.tafd_front_strip <= 13) {
				taf_flag = 2;
			} else {
				valid |= 2;
			}
		} else if (recoil_name == "1H") {
			if (channel.tafcsi_valid && channel.tafd_front_strip >= 12) {
				taf_flag = 1;
			} else {
				valid |= 2;
			}
		}

		if (!channel.ppac_valid) valid |= 4;
		if (pow(channel.tx-2.5, 2.0) + pow(channel.ty+2.0, 2.0) < 196.0) {
			target_flag = 1;
		} else {
			target_flag = 0;
			valid |= 8;
		}

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
		be_kinetic = channel.fragment_kinetic[0];
		// 4He kinetic energy
		he_kinetic = channel.fragment_kinetic[1];
		// deutron kinetic energy
		d_kinetic = channel.recoil_kinetic;

		// 10Be momentum
		double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
		ROOT::Math::XYZVector p_be = d_be * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic);
		ROOT::Math::XYZVector p_he = d_he * he_momentum;

		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
		ROOT::Math::XYZVector p_d = d_d * d_momentum;

		// beam 14C momentum vector
		ROOT::Math::XYZVector p_c = p_be + p_he + p_d;

		// 14C momentum
		c_momentum = p_c.R();
		// 14C kinematic energy
		c_kinetic =
			sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

		// three body Q value
		q = be_kinetic + he_kinetic + d_kinetic - c_kinetic;

		// correct Q
		if (taf_flag == 1) {
			q += q_tafd_correct[channel.taf_index];
		} else if (taf_flag == 2) {
			q -= q_csi_correct[channel.csi_index];
		}
		// get 10Be state from Q value
		if (q < -11.0 && q > -13.0) be_state = 0;
		else if (q < -15.0 && q > -16.0) be_state = 1;
		else if (q < -18.0 && q > -20.0) be_state = 2;
		// else if (q < -19 && q > -20.5) be_state = 3;
		else be_state = -1;

		double mass_excited_10be = mass_10be
			+ (be_state >= 0 ? be_excited_energy[be_state] : 0.0);
 		be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
		p_be = d_be * be_momentum;
		// excited 14C momentum vector
		ROOT::Math::XYZVector p_excited_c = p_be + p_he;
		// excited 14C momentum
		double excited_c_momentum = p_excited_c.R();
		// excited 14C total energy
		double excited_c_energy = (be_kinetic + mass_excited_10be)
			+ (he_kinetic + mass_4he);
		// excited 14C mass
		double excited_c_mass = sqrt(
			pow(excited_c_energy, 2.0) - pow(excited_c_momentum, 2.0)
		);
		// excited energy of 14C
		excited_energy = excited_c_mass - mass_14c;


		if (c_kinetic < 360.0) valid |= 16;
		if (valid == 0 && taf_flag == 1) {
			// fill Q spectrums
			hq_taf[channel.taf_index].Fill(
				q - q_tafd_correct[channel.taf_index]
			);
			hq.Fill(q);
		}

		// // image 1H part
		// // calculated recoil kinetic
		// double image_p_kinetic = channel.tafd_energy + h1_csi.Energy(
		// 	d_d.Theta(), channel.tafd_energy,
		// 	correct_tafd_thickness[channel.taf_index]
		// );
		// double image_p_kinetic_target =
		// 	h2_target.Energy(-0.5/cos(d_d.Theta()), image_p_kinetic);
		// // 2H momentum
		// double image_p_momentum =
		// 	MomentumFromKinetic(mass_1h, image_p_kinetic_target);
		// ROOT::Math::XYZVector image_p_p = d_d * image_p_momentum;
		// // beam 14C momentum vector
		// ROOT::Math::XYZVector image_p_c = p_be + p_he + image_p_p;
		// // 14C momentum
		// double image_c_momentum_target = image_p_c.R();
		// // 14C kinematic energy
		// double image_c_kinetic_target =
		// 	sqrt(pow(image_c_momentum_target, 2.0)
		// 	+ pow(mass_14c, 2.0)) - mass_14c;
		// // three body Q value
		// image_q = be_kinetic_target + he_kinetic_target
		// 	+ image_p_kinetic_target - image_c_kinetic_target;
		// if (image_q > 1) valid |= 64;

		// fill particle direction
		be_dx = d_be.X();
		be_dy = d_be.Y();
		be_dz = d_be.Z();
		he_dx = d_he.X();
		he_dy = d_he.Y();
		he_dz = d_he.Z();
		d_dx = d_d.X();
		d_dy = d_d.Y();
		d_dz = d_d.Z();

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