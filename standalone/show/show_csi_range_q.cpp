#include <iostream>
#include <vector>

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


double QFit2(double *x, double *par) {
	return par[0] * exp(-0.5*pow((x[0]-par[1])/par[2], 2.0))
		+ par[3] * exp(-0.5*pow((x[0]-par[1]-2.5)/par[4], 2.0))
		+ par[5] * exp(-0.5*pow((x[0]-par[1]-6.0)/par[6], 2.0));
}


double QFit3(double *x, double *par) {
	return par[0] * exp(-0.5*pow((x[0]-par[1])/par[2], 2.0))
		+ par[3] * exp(-0.5*pow((x[0]-par[1]-2.811)/par[4], 2.0))
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

	// input CsI energy range file name
	TString csi_range_file_name = TString::Format(
		"%s%stafcsi-energy-range-0618-0746.root", kGenerateDataPath, kShowDir
	);
	// input file
	TFile csi_range_file(csi_range_file_name, "read");
	// CsI energy range histograms
	TH1F *hist_range[12][16][48];
	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 16; ++j) {
			for (int k = 0; k < 48; ++k) {
				hist_range[i][j][k] = (TH1F*)csi_range_file.Get(
					TString::Format("he%dt%dr%d", i, 2*j+140, k)
				);
			}
		}
	}

	// output file name
	TString output_file_name = TString::Format(
		"%s%stafcsi-range-q.root", kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histograms
	// histogram-Q-entry, single entry's Q values
	std::vector<TH1F> hqe;
	for (long long i = 0; i < ipt->GetEntriesFast(); ++i) {
		hqe.emplace_back(
			TString::Format("hqe%lld", i),
			TString::Format("Q value of entry %lld, selected thickness", i),
			90, -23, -8
		);
	}
	// histogram-Q-CsI, single CsI's Q values of all entries
	std::vector<TH1F> hqc[12];
	for (int i = 0; i < 12; ++i) {
		for (int t = 0; t < 16; ++t) {
			hqc[i].emplace_back(
				TString::Format("hqc%dt%d", i, t*2+140),
				TString::Format("Q value of CsI %d, thickness %d", i, 2*t+140),
				150, -23, -8
			);

		}
	}
	// histogram-mean-Q-CsI, single CsI's mean Q value, at selected thickness
	std::vector<TH1F> hmqc;
	for (int i = 0; i < 12; ++i) {
		hmqc.emplace_back(
			TString::Format("hmqc%d", i),
			TString::Format("mean Q value of CsI %d", i),
			90, -23, -8
		);
	}
	// histogram-mean-Q, at selected thickness respectively
	TH1F hmq("hmq", "mean Q value", 90, -23, -8);
	// histogram-Q-aligned
	TH1F hqa("hqa", "align Q value", 150, -23, -8);

	// random CsI energy event tree
	TTree opt("tree", "random CsI energy events");
	// output event
	ThreeBodyInfoEvent random_event;
	// setup output branches
	random_event.SetupOutput(&opt);

	int taf_selected_thickness[6] = {
		166, 154, 158, 162, 150, 164
	};

	double q_offset[12] = {
		-0.46, -0.25, 0.18, 0.49,
		0.49, 0.55, 0.28, 0.30,
		-0.04, -0.10, -0.12, -0.64
	};

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
			&& event.csi_channel < 11000.0
		) {
			int range_index = int((event.csi_channel - 1400.0) / 200.0);
			for (int i = 0; i < 1000; ++i) {
				for (int t = 0; t < 16; ++t) {
					double csi_energy =
						hist_range[event.csi_index][t][range_index]->GetRandom();
					double q = ThreeBodyProcess(event, csi_energy);
					// fill hqc
					hqc[event.csi_index][t].Fill(q);
					// fill hqe and hqa at selected thickness
					if (2*t+140 == taf_selected_thickness[event.csi_index/2]) {
						hqe[entry].Fill(q);
						hqa.Fill(q-q_offset[event.csi_index]);
					}

					random_event = event;
					random_event.q = q;
					opt.Fill();
				}
				
			}

			// fit and get mean Q
			TF1 *f1 = new TF1(
				TString::Format("fqe%lld", entry), "gaus", -23.0, -8.0
			);
			f1->SetParameter(0, 100.0);
			f1->SetParameter(1, -15.0);
			f1->SetParameter(2, 1.0);
			if (hqe[entry].GetMean() < 0.0) {
				hqe[entry].Fit(f1, "RQ+");
				double mean_q = f1->GetParameter(1);
				hmqc[event.csi_index].Fill(mean_q);
				hmq.Fill(mean_q-q_offset[event.csi_index]);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit parameters
	double offset[12][16][2];
	double sigma[12][16][3];

	// fit single CsI's Q value, with different thickness
	// loop CsI index
	for (int i = 0; i < 12; ++i) {
		// loop TAFD thickness
		for (int j = 0; j < 16; ++j) {
			TF1 *f3p = new TF1(
				TString::Format("fs%dt%d", i, 2*j+140), QFit3, -23, -8, 8
			);
			// set initial value of parameters
			f3p->SetParameter(0, 50000.0);
			f3p->SetParameter(1, -18.0);
			f3p->SetParameter(2, 1.0);
			f3p->SetParameter(3, 80000.0);
			f3p->SetParameter(4, 1.0);
			f3p->SetParameter(5, 20000.0);
			f3p->SetParameter(6, -12.0);
			f3p->SetParameter(7, 1.0);
			// fit
			hqc[i][j].Fit(f3p, "QR+");
			// get fitted result
			offset[i][j][0] = f3p->GetParameter(1) + 18.1915;
			offset[i][j][1] = f3p->GetParameter(6) + 12.0125;
			sigma[i][j][0] = f3p->GetParameter(2);
			sigma[i][j][1] = f3p->GetParameter(4);
			sigma[i][j][2] = f3p->GetParameter(7);
		}
	}

	std::cout << "CsI thickness(um) offset1 offset2 sigma1 sigma2 sigma3\n";
	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 16; ++j) {
			std::cout << i << " " << 2*j+140 << " "
				<< offset[i][j][0] << " " << offset[i][j][1] << " "
				<< sigma[i][j][0] << " " << sigma[i][j][1] << " "
				<< sigma[i][j][2] << "\n";
		}
	}


	// fit mean Q value in selected thickness (hmq)
	TF1 *fmq = new TF1("fmq", QFit, -23, -8, 9);
	// set initial value
	double fmq_init_value[9] = {
		100.0, -18.0, 1.0,
		120.0, -15.5, 1.0,
		40.0, -12.0, 1.0
	};
	fmq->SetParameters(fmq_init_value);
	// fit
	hmq.Fit(fmq, "QR+");
	std::cout << "\nMean Q value:\n"
		<< fmq->GetParameter(1) << ", "
		<< fmq->GetParameter(4) << ", "
		<< fmq->GetParameter(7) << "\n"
		<< fmq->GetParameter(2) << ", "
		<< fmq->GetParameter(5) << ", "
		<< fmq->GetParameter(8) << "\n";


	// fit aligned Q value (hqa)
	TF1 *fqa = new TF1("fqa", QFit, -23, -8, 9);
	// set initial value
	double fqa_init_value[9] = {
		60000.0, -18.0, 1.0,
		70000.0, -15.5, 1.0,
		20000.0, -12.0, 1.0
	};
	fqa->SetParameters(fqa_init_value);
	// fit
	hqa.Fit(fqa, "QR+");
	std::cout << "\nAligned Q value:\n"
		<< fqa->GetParameter(1) << ", "
		<< fqa->GetParameter(4) << ", "
		<< fqa->GetParameter(7) << "\n"
		<< fqa->GetParameter(2) << ", "
		<< fqa->GetParameter(5) << ", "
		<< fqa->GetParameter(8) << "\n";



	opf.cd();
	for (auto &hist : hqe) hist.Write();
	for (int i = 0; i < 12; ++i) {
		for (auto &hist : hqc[i]) hist.Write();
	}
	for (auto &hist : hmqc) hist.Write();
	hmq.Write();
	hqa.Write();
	opt.Write();
	opf.Close();

	// close input files
	csi_range_file.Close();
	ipf.Close();
	return 0;
}
