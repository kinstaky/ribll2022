#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TF1.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/ppac_track.h"

using namespace ribll;

constexpr int strips = 8;
constexpr int param_counts[strips] = {
	3, 3, 3, 3, 3, 3, 3, 3
};

constexpr double tafd_phi_start[6] = {
	117.6, 57.6, -2.4, -62.4, -122.4, 177.6
};


constexpr double theory_q[4] = {-12.0, -15.5, -18.0, -19.5};


bool NextGroup(int *param) {
	param[0]++;
	for (int i = 0; i < strips; ++i) {
		if (param[i] == param_counts[i]) {
			param[i] = 0;
			if (i == strips-1) return false;
			param[i+1]++;
		} else {
			break;
		}
	}
	return true;
}


int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody.root",
		kGenerateDataPath, kInformationDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	// detect event
	ThreeBodyInfoEvent info;
	// setup input branches
	info.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sresist-strip.root", kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of Q value
	std::vector<TH1F> hist_q;
	// histogram of beam kinetic
	std::vector<TH1F> hist_c_kinetic;
	// output tree
	TTree opt("tree", "resist strip");
	// output data
	bool valid;
	double q_mean, q_sigma;
	double represent_q[4];
	// represent state, the last index is the state of closest one
	int represent_state[4];
	// Q difference with theory, the last index is the closest one
	double q_diff[4];
	// the most same state
	int most_same_state;

	// setup output branches
	opt.Branch("vaild", &valid, "valid/O");
	opt.Branch("q_mean", &q_mean, "qmean/D");
	opt.Branch("q_sigma", &q_sigma, "qsigma/D");
	opt.Branch("represent_q", represent_q, "rq[4]/D");
	opt.Branch("represent_state", represent_state, "rs[4]/I");
	opt.Branch("q_diff", q_diff, "qdiff[4]/D");
	opt.Branch("same_state", &most_same_state, "ss/I");

	// parameters
	int strip_params[strips];

	// statistics
	int total_events = 0;
	int same_state_num[3];
	for (int i = 0; i < 3; ++i) same_state_num[i] = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Looping   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);

		// add histograms
		hist_q.emplace_back(
			TString::Format("hq%lld", entry), "Q values",
			500, -40, 10
		);
		hist_c_kinetic.emplace_back(
			TString::Format("hck%lld", entry), "beam kinetic energy",
			100, 350, 400
		);

		// check TAF
		if (info.taf_flag != 0) {
			valid = false;
			opt.Fill();
			continue;
		}

		// check PPAC
		if (info.xppac_xflag == 0 || info.xppac_yflag == 0) {
			valid = false;
			opt.Fill();
			continue;
		}
		valid = true;
		int ppac_x_index = 0;
		int ppac_y_index = 0;
		for (; ppac_x_index < 3; ++ppac_x_index) {
			if ((info.xppac_xflag & (1 << ppac_x_index)) != 0) break;
		}
		for (; ppac_y_index < 3; ++ppac_y_index) {
			if ((info.xppac_yflag & (1 << ppac_y_index)) != 0) break;
		}

		// initialize parameters
		for (int i = 0; i < strips; ++i) strip_params[i] = 0;
		// loop
		do {
			// change positions
			double ct0x[2], ct0y[2], ctafx, ctafy, cppacx, cppacy;
			// change 10Be position
			ct0x[0] = info.be_x[0] + (strip_params[0]-1)*0.5;
			ct0y[0] = info.be_y[0] + (strip_params[4]-1)*0.5;
			// change 4He position
			ct0x[1] = info.he_x[0] + (strip_params[1]-1)*0.5;
			ct0y[1] = info.he_y[0] + (strip_params[5]-1)*0.5;
			// change 2H position
			double ctaf_r =
				102.5/16.0 * (info.d_x_strip + strip_params[2]*0.5) + 32.6;
			double ctaf_phi =
				-55.2/8.0 * (info.d_y_strip + strip_params[6]*0.5)
				+ tafd_phi_start[info.csi_index/2];
			ctaf_phi *= TMath::DegToRad();
			double mid_phi =
				(tafd_phi_start[info.csi_index/2] - 27.6) * TMath::DegToRad();
			ctafx = ctaf_r * cos(ctaf_phi) + 34.4*cos(mid_phi);
			ctafy = ctaf_r * sin(ctaf_phi) + 34.4*sin(mid_phi);
			// change 14C position
			cppacx = info.xppac_x[ppac_x_index]
				+ ppac_correct[0][ppac_x_index]
				+ (strip_params[3]-1)*0.5;
			cppacy = info.xppac_y[ppac_y_index]
				+ ppac_correct[1][ppac_y_index]
				+ (strip_params[7]-1)*0.5;

			double d_kinetic = info.tafd_energy + pow(
				(info.csi_channel - csi_param[info.csi_index][2])
					/ csi_param[info.csi_index][0],
				1.0 / csi_param[info.csi_index][1]
			);

			// reaction point from single PPAC track
			double tx = DeutronRelativeApproximateTrack(
				info.t0_energy[0], info.t0_energy[1], d_kinetic,
				ppac_xz[ppac_x_index], ct0x[0], ct0x[1], ctafx, ctafy, cppacx
			);
			double ty = DeutronRelativeApproximateTrack(
				info.t0_energy[0], info.t0_energy[1], d_kinetic,
				ppac_yz[ppac_y_index], ct0y[0], ct0y[1], ctafy, ctafx, cppacy
			);

			// momentum vector of fragments
			ROOT::Math::XYZVector fp[2];
			for (size_t i = 0; i < 2; ++i) {
				fp[i] = ROOT::Math::XYZVector(
					ct0x[i] - tx,
					ct0y[i] - ty,
					100.0
				);
				// fragment mass
				double mass = i == 0 ? mass_10be : mass_4he;
				// fragment kinetic energy
				double kinetic = info.t0_energy[i];
				// fragment momentum value
				double momentum = sqrt(
					pow(kinetic, 2.0) + 2.0 * kinetic * mass
				);
				fp[i] = fp[i].Unit() * momentum;
			}
			// recoil 2H momentum vector
			ROOT::Math::XYZVector rp(
				ctafx - tx,
				ctafy - ty,
				135.0
			);
			// recoild 2H momentum value
			double recoil_momentum = sqrt(
				pow(d_kinetic, 2.0) + 2.0 * d_kinetic * mass_2h
			);
			rp = rp.Unit() * recoil_momentum;
			// rebuild beam momentum vector
			ROOT::Math::XYZVector bp = fp[0] + fp[1] + rp;
			// rebuild beam kinetic energy
			double c_kinetic = sqrt(bp.Dot(bp) + pow(mass_14c, 2.0)) - mass_14c;
			// Q value
			double rebuild_q = info.t0_energy[0] + info.t0_energy[1]
				+ d_kinetic - c_kinetic;

			hist_q[entry].Fill(rebuild_q);
			hist_c_kinetic[entry].Fill(c_kinetic);

		} while (NextGroup(strip_params));

		// fit Q value histogram
		TF1 *q_fit = new TF1(
			TString::Format("qf%lld", entry), "gaus", -40, 10
		);
		// set initial parameters
		q_fit->SetParameter(0, 100);
		q_fit->SetParameter(1, info.q);
		q_fit->SetParameter(2, 1.0);
		// fit
		hist_q[entry].Fit(q_fit, "RQ+");
		// get fitted parameter
		q_mean = q_fit->GetParameter(1);
		q_sigma = q_fit->GetParameter(2);

		// represent Q
		for (int i = 0; i < 3; ++i) {
			// change positions
			double ct0x[2], ct0y[2], ctafx, ctafy, cppacx, cppacy;
			// keep 10Be position
			ct0x[0] = info.be_x[0];
			ct0y[0] = info.be_y[0];
			// keep 4He position
			ct0x[1] = info.he_x[0];
			ct0y[1] = info.he_y[0];
			// change 2H position
			double ctaf_r =
				102.5/16.0 * (info.d_x_strip + i*0.5) + 32.6;
			double ctaf_phi =
				-55.2/8.0 * (info.d_y_strip + 0.5)
				+ tafd_phi_start[info.csi_index/2];
			ctaf_phi *= TMath::DegToRad();
			double mid_phi =
				(tafd_phi_start[info.csi_index/2] - 27.6) * TMath::DegToRad();
			ctafx = ctaf_r * cos(ctaf_phi) + 34.4*cos(mid_phi);
			ctafy = ctaf_r * sin(ctaf_phi) + 34.4*sin(mid_phi);
			// change 14C position
			cppacx =
				info.xppac_x[ppac_x_index] + ppac_correct[0][ppac_x_index];
			cppacy =
				info.xppac_y[ppac_y_index] + ppac_correct[1][ppac_y_index];

			double d_kinetic = info.tafd_energy + pow(
				(info.csi_channel - csi_param[info.csi_index][2])
					/ csi_param[info.csi_index][0],
				1.0 / csi_param[info.csi_index][1]
			);

			// reaction point from single PPAC track
			double tx = DeutronRelativeApproximateTrack(
				info.t0_energy[0], info.t0_energy[1], d_kinetic,
				ppac_xz[ppac_x_index], ct0x[0], ct0x[1], ctafx, ctafy, cppacx
			);
			double ty = DeutronRelativeApproximateTrack(
				info.t0_energy[0], info.t0_energy[1], d_kinetic,
				ppac_yz[ppac_y_index], ct0y[0], ct0y[1], ctafy, ctafx, cppacy
			);

			// momentum vector of fragments
			ROOT::Math::XYZVector fp[2];
			for (size_t i = 0; i < 2; ++i) {
				fp[i] = ROOT::Math::XYZVector(
					ct0x[i] - tx,
					ct0y[i] - ty,
					100.0
				);
				// fragment mass
				double mass = i == 0 ? mass_10be : mass_4he;
				// fragment kinetic energy
				double kinetic = info.t0_energy[i];
				// fragment momentum value
				double momentum = sqrt(
					pow(kinetic, 2.0) + 2.0 * kinetic * mass
				);
				fp[i] = fp[i].Unit() * momentum;
			}
			// recoil 2H momentum vector
			ROOT::Math::XYZVector rp(
				ctafx - tx,
				ctafy - ty,
				135.0
			);
			// recoild 2H momentum value
			double recoil_momentum = sqrt(
				pow(d_kinetic, 2.0) + 2.0 * d_kinetic * mass_2h
			);
			rp = rp.Unit() * recoil_momentum;
			// rebuild beam momentum vector
			ROOT::Math::XYZVector bp = fp[0] + fp[1] + rp;
			// rebuild beam kinetic energy
			double c_kinetic = sqrt(bp.Dot(bp) + pow(mass_14c, 2.0)) - mass_14c;
			// Q value
			represent_q[i] = info.t0_energy[0] + info.t0_energy[1]
				+ d_kinetic - c_kinetic;
		}

		// find the state of represent Q, reference closest Q value
		for (int i = 0; i < 3; ++i) {
			q_diff[i] = represent_q[i] - theory_q[0];
			represent_state[i] = 0;
			for (int j = 1; j < 4; ++j) {
				double diff = represent_q[i] - theory_q[j];
				if (fabs(diff) < fabs(q_diff[i])) {
					q_diff[i] = diff;
					represent_state[i] = j;
				}
			}
		}
		// search closest state
		represent_q[3] = represent_q[0];
		represent_state[3] = represent_state[0];
		q_diff[3] = q_diff[0];
		if (q_diff[1] < q_diff[3]) {
			q_diff[3] = q_diff[1];
			represent_q[3] = represent_q[1];
			represent_state[3] = represent_state[1];
		}
		if (q_diff[2] < q_diff[3]) {
			q_diff[3] = q_diff[2];
			represent_q[3] = represent_q[2];
			represent_state[3] = represent_state[2];
		}

		// represent Q in same state
		int same_num = 0;
		if (represent_state[0] == represent_state[1]) {
			++same_num;
			most_same_state = represent_state[0];
		}
		if (represent_state[0] == represent_state[2]) {
			++same_num;
			most_same_state = represent_state[0];
		}
		if (represent_state[1] == represent_state[2]) {
			++same_num;
			most_same_state = represent_state[1];
		}
		if (same_num == 0) ++same_state_num[0];
		else if (same_num == 1) ++same_state_num[1];
		else if (same_num == 3) ++same_state_num[2];
		++total_events;
		// search for the most same state
		if (same_num == 0) {
			most_same_state = -1;
		}
		// fill tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "total events " << total_events << "\n"
		<< "0 same state " << same_state_num[0] << "\n"
		<< "2 same state " << same_state_num[1] << "\n"
		<< "3 same state " << same_state_num[2] << "\n";

	// save histograms
	for (auto &hist : hist_q) hist.Write();
	for (auto &hist : hist_c_kinetic) hist.Write();
	// save tree
	opt.Write();
	// close files
	ipf.Close();
	opf.Close();
	return 0;
}