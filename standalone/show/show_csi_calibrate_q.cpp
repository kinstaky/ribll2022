#include <iostream>
#include <vector>
#include <fstream>

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <Math/Vector3D.h>
#include <TRandom3.h>
#include <TF1.h>

#include "include/event/threebody_info_event.h"

using namespace ribll;


/// @brief rebuild threebody reaction process
/// @param[in] event input event
/// @param[in] csi_energy CsI energy
/// @returns Q value
///
double ThreeBodyProcess(const ThreeBodyInfoEvent &event, double csi_energy) {
	double tx = event.xptx;
	double ty = event.xpty;

	// 10Be momentum
	double be_momentum = MomentumFromKinetic(mass_10be, event.t0_energy[0]);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - tx,
		event.be_y[0] - ty,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, event.t0_energy[1]);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - tx,
		event.he_y[0] - ty,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	double taf_energy = event.tafd_energy + csi_energy;
	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, taf_energy);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		event.d_x - tx,
		event.d_y - ty,
		135.0
	);
	p_d = p_d.Unit() * d_momentum;

	// beam 14C momentum vector
	ROOT::Math::XYZVector p_beam = p_be + p_he + p_d;

	// 14C momentum
	double beam_momentum = p_beam.R();
	// 14C kinematic energy
	double c14_kinetic =
		sqrt(pow(beam_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

	// three-fold Q value
	double q = event.t0_energy[0] + event.t0_energy[1]
		+ taf_energy - c14_kinetic;

	return q;
}


double QFit(double *x, double *par) {
	return par[0] * exp(-0.5*pow((x[0]-par[1])/par[2], 2.0))
		+ par[3] * exp(-0.5*pow((x[0]-par[4])/par[5], 2.0))
		+ par[6] * exp(-0.5*pow((x[0]-par[7])/par[8], 2.0));
}


double QFit3(double *x, double *par) {
	return par[0] * exp(-0.5*pow((x[0]-par[1])/par[2], 2.0))
		+ par[3] * exp(-0.5*pow((x[0]-par[1]-2.5)/par[4], 2.0))
		+ par[5] * exp(-0.5*pow((x[0]-par[6])/par[7], 2.0));
}


int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody.root", kGenerateDataPath, kInformationDir
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
	// input event
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);


	// output file name
	TString output_file_name = TString::Format(
		"%s%stafcsi-calibrate-q.root", kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histograms
	// histogram-linear-calibration-Q-CsI
	std::vector<TH1F> hlqc;
	for (int i = 0; i < 12; ++i) {
		hlqc.emplace_back(
			TString::Format("hlqc%d", i),
			TString::Format("linear calibrated 1 Q value of CsI %d", i),
			90, -23, -8
		);
	}
	// histogram-linear-Q-aligned
	TH1F hlqa("hlqa", "aligned linear calibrated CsI Q value", 90, -23, -8);
	// histogram-power-Q-aligned
	TH1F hpqa("hpqa", "aligned power calibrated CsI Q value", 90, -23, -8);
	// histogram-optimized-2-Q-aligned
	TH1F hoqa("hoqa", "aligned optimized CsI Q value", 90, -23, -8);
	// output tree
	TTree opt("tree", "different calibrated Q");
	// output data, different Q
	double linear_q, power_q, optimized_q;
	bool good;
	// setup output branches
	opt.Branch("linear_q", &linear_q, "lq/D");
	opt.Branch("power_q", &power_q, "pq/D");
	opt.Branch("opt_q", &optimized_q, "oq/D");
	opt.Branch("good", &good, "good/O");

	// Q value offset
	double q_offset[12] = {
		-0.39, -0.13, 0.21, 0.56,
		0.49, 0.50, 0.19, 0.29,
		-0.12, -0.20, -0.39, -0.66
	};

	// calibrate parameters
	double linear_param[12][2];
	double power_param[12][3];


	// read linear calibrate parameters
	// loop TAFs
	for (int index = 0; index < 6; ++index) {
		// file name   
		TString file_name = TString::Format(
			"%s%staf%dcsi-cali-param-2H-linear.txt",
			kGenerateDataPath,
			kCalibrationDir,
			index
		);
		// input stream
		std::ifstream linear_fin(file_name.Data());
		if (!linear_fin.good()) {
			std::cerr << "Error: Get parameters from file "
				<< file_name << " failed.\n";
			return -1;
		}
		linear_fin >> linear_param[index*2][0] >> linear_param[index*2][1]
			>> linear_param[index*2+1][0] >> linear_param[index*2+1][1];
		// close files
		linear_fin.close();
	}
	// show read parameteres
	std::cout << "Read linear calibrate parameters\n";
	for (int i = 0; i < 12; ++i) {
		std::cout << linear_param[i][0] << ", "
			<< linear_param[i][1] << "\n";
	}

	// read power calibrate parameters
	// loop TAFs
	for (int index = 0; index < 6; ++index) {
		// file name   
		TString file_name = TString::Format(
			"%s%staf%dcsi-cali-param-2H.txt",
			kGenerateDataPath,
			kCalibrationDir,
			index
		);
		// input stream
		std::ifstream power_fin(file_name.Data());
		if (!power_fin.good()) {
			std::cerr << "Error: Get parameters from file "
				<< file_name << " failed.\n";
			return -1;
		}
		for (int i = 0; i < 2; ++i) {
			power_fin >> power_param[index*2+i][0]
				>> power_param[index*2+i][1]
				>> power_param[index*2+i][2];
		}
		// close files
		power_fin.close();
	}
	// show read parameteres
	std::cout << "Read power calibrate parameters 2\n";
	for (int i = 0; i < 12; ++i) {
		std::cout << power_param[i][0] << ", "
			<< power_param[i][1] << ", "
			<< power_param[i][2] << "\n";
	}


	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Processing   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);

		if (
			event.taf_flag == 0
			&& (event.ppac_flag & 1) == 1
			&& event.bind == 0
			&& !event.hole[0]
			&& !event.hole[1]
			&& event.csi_channel > 1400.0
		) {
			good = true;
		} else {
			good = false;
		}
		// linear-calibrated 1 CsI energy
		double linear_csi_energy =
			linear_param[event.csi_index][0]
			+ linear_param[event.csi_index][1] * event.csi_channel;
		// linear-calibrated 1 Q value
		linear_q = ThreeBodyProcess(event, linear_csi_energy);
		// fill to histogram
		hlqc[event.csi_index].Fill(linear_q);
		hlqa.Fill(linear_q - q_offset[event.csi_index]);
	
		// power-calibrated CsI energy
		double power_csi_energy = pow(
			(event.csi_channel - power_param[event.csi_index][2])
				/ power_param[event.csi_index][0],
			1.0 / power_param[event.csi_index][1]
		);
		// power-calibrated CsI Q value
		power_q = ThreeBodyProcess(event, power_csi_energy);
		// fill to histogram
		hpqa.Fill(power_q - q_offset[event.csi_index]);

		// optimized CsI energy
		double opt_csi_energy = pow(
			(event.csi_channel - csi_param[event.csi_index][2])
				/ csi_param[event.csi_index][0],
			1.0 / csi_param[event.csi_index][1]
		);
		// optimized Q value
		optimized_q = ThreeBodyProcess(event, opt_csi_energy);
		// fill to histogram
		hoqa.Fill(optimized_q - q_offset[event.csi_index]);

		// fill tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit parameters
	double offset[12][2];
	double sigma[12][3];

	// fit linear-calibrated CsI Q (hlqc)
	// loop CsI index
	for (int i = 0; i < 12; ++i) {
		TF1 *flcqc = new TF1(
			TString::Format("flcqc%d", i), QFit3, -23, -8, 8
		);
		// set initial value of parameters
		double flcqc_init_value[8] = {
			20.0, -18.0, 1.0,
			20.0, 1.0,
			10.0, -12.0, 1.0
		};
		flcqc->SetParameters(flcqc_init_value);
		// fit
		hlqc[i].Fit(flcqc, "QR+");
		// get fitted result
		offset[i][0] = flcqc->GetParameter(1) + 18.0;
		offset[i][1] = flcqc->GetParameter(6) + 12.0;
		sigma[i][0] = flcqc->GetParameter(2);
		sigma[i][1] = flcqc->GetParameter(4);
		sigma[i][2] = flcqc->GetParameter(7);
	}

	std::cout << "CsI offset1 offset2 sigma1 sigma2 sigma3\n";
	for (int i = 0; i < 12; ++i) {
		std::cout << i << " "
			<< offset[i][0] << " "
			<< offset[i][1] << " "
			<< sigma[i][0] << " "
			<< sigma[i][1] << " "
			<< sigma[i][2] << "\n";
	}


	// 4 histograms to fit 
	TH1F *fit_hist[3] = {&hlqa, &hpqa, &hoqa};
	// title to print
	std::string title[3] = {
		"\nAligned linear calibrated Q value:\n",
		"\nAligned power calibrated Q value:\n",
		"\nAligned optimized Q value:\n"
	};
	// initial value
	double fit_init_value[9] = {
		100.0, -18.0, 1.0,
		120.0, -15.5, 1.0,
		40.0, -12.0, 1.0
	};
	// loop and fit 4 histograms
	for (int i = 0; i < 3; ++i) {
		TF1 *f = new TF1(TString::Format("f%d", i), QFit, -23, -8, 9);
		// set initial value
		f->SetParameters(fit_init_value);
		// fit
		fit_hist[i]->Fit(f, "QR+");
		// print
		std::cout << title[i]
			<< f->GetParameter(1) << ", "
			<< f->GetParameter(4) << ", "
			<< f->GetParameter(7) << "\n"
			<< f->GetParameter(2) << ", "
			<< f->GetParameter(5) << ", "
			<< f->GetParameter(8) << "\n";
	}

	// save histograms
	opf.cd();
	for (auto &hist : hlqc) hist.Write();
	hlqa.Write();
	hpqa.Write();
	hoqa.Write();
	opt.Write();
	// close output file
	opf.Close();
	// close input file
	ipf.Close();
	return 0;
}
