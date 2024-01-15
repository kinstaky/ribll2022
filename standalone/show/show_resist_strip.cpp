#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TF1.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"

using namespace ribll;

constexpr double pi = 3.1415926;
constexpr double u = 931.494;
constexpr double c14_mass = 13.9999505089 * u;
constexpr double be10_mass = 10.0113403769 * u;
constexpr double he4_mass = 4.0015060943 * u;
constexpr double h2_mass = 2.0135531980 * u;

constexpr double ppac_xz[3] = {-695.2, -454.2, -275.2};
constexpr double ppac_yz[3] = {-689.2, -448.2, -269.2};

constexpr double ppac_correctx[3] = {0.0, -2.23, -3.40};
constexpr double ppac_correcty[3] = {0.0, 0.84, 1.78};

constexpr int strips = 8;
constexpr int param_counts[strips] = {
	3, 3, 3, 3, 3, 3, 3, 3
};

constexpr double tafd_phi_start[6] = {
	117.6, 57.6, -2.4, -62.4, -122.4, 177.6
};


constexpr double csi_param[12][3] = {
	{273.278, 0.9159, -350.988},
	{280.876, 0.8937, -584.531},
	{386.048, 0.8575, -616.000},
	{319.845, 0.9122, -368.000},
	{346.592, 0.9312, -843.328},
	{319.096, 0.9146, -377.708},
	{289.534, 0.9638, -516.000},
	{381.000, 0.8550, -632.000},
	{510.623, 0.8433, -916.000},
	{386.525, 0.8638, -610.244},
	{285.668, 0.9571, -413.028},
	{300.415, 0.9672, -165.075}
};


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


/// @brief single PPAC tracking, optimize deutron relative approximate method
/// @param[in] be_kinetic kinetic energy of 10Be
/// @param[in] he_kinetic kinetic energy of 4He
/// @param[in] d_kinetic kinetic energy of 2H
/// @param[in] ppac_z distance of PPAC from target
/// @param[in] be_position 10Be position on T0
/// @param[in] he_position 4He position on T0
/// @param[in] d_position 2H position on TAFD
/// @param[in] d_position_v 2H position on TAFD other direction
/// @param[in] c_position beam position on PPAC
/// @returns reaction point
///
double DeutronRelativeApproximateTrack(
	double be_kinetic,
	double he_kinetic,
	double d_kinetic,
	double ppac_z,
	double be_position,
	double he_position,
	double d_position,
	double d_position_v,
	double c_position
) {
	// 10Be parameter
	double a_be = sqrt(
		(2.0 * be10_mass + be_kinetic) * be_kinetic
	) / 100.0;
	// 4He parameter
	double a_he = sqrt(
		(2.0 * he4_mass + he_kinetic) * he_kinetic
	) / 100.0;
	// 2H parameter
	double a_d = sqrt(
		(2.0 * h2_mass + d_kinetic) * d_kinetic
	) / sqrt(
		135.0 * 135.0
		+ d_position * d_position
		+ d_position_v * d_position_v
	);
	// 14C parameter
	double a_c = -sqrt(
		(2.0 * c14_mass + 386.0) * 386.0
	) / ppac_z;
	// calculate
	double numerator = a_be * be_position;
	numerator += a_he * he_position;
	numerator += a_d * d_position;
	numerator += a_c * c_position;
	double denominator = a_be + a_he + a_d + a_c;
	return numerator / denominator;
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
		"%s%sresist_strip.root", kGenerateDataPath, kShowDir
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
	// setup output branches
	opt.Branch("vaild", &valid, "valid/O");
	opt.Branch("q_mean", &q_mean, "qmean/D");
	opt.Branch("q_sigma", &q_sigma, "qsigma/D");

	// parameters
	int strip_params[strips];

	// total number of entries
	long long entries = ipt->GetEntries();
	for (long long entry = 0; entry < entries; ++entry) {
		ipt->GetEntry(entry);
		// std::cout << "entry: " << entry << "\n";

		// add histograms
		hist_q.emplace_back(
			TString::Format("hq%lld", entry), "Q values",
			500, -40, 10
		);
		hist_c_kinetic.emplace_back(
			TString::Format("hck%lld", entry), "beam kinetic energy",
			100, 350, 400
		);

		// check PPAC
		if (info.ppac_xflag == 0 || info.ppac_yflag == 0) {
			valid = false;
			opt.Fill();
			continue;
		}
		valid = true;
		int ppac_x_index = 0;
		int ppac_y_index = 0;
		for (; ppac_x_index < 3; ++ppac_x_index) {
			if ((info.ppac_xflag & (1 << ppac_x_index)) != 0) break;
		}
		for (; ppac_y_index < 3; ++ppac_y_index) {
			if ((info.ppac_yflag & (1 << ppac_y_index)) != 0) break;
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
			cppacx = info.ppac_x[ppac_x_index]
				+ ppac_correctx[ppac_x_index]
				+ (strip_params[3]-1)*0.5;
			cppacy = info.ppac_y[ppac_y_index]
				+ ppac_correcty[ppac_y_index]
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
				double mass = i == 0 ? be10_mass : he4_mass;
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
				pow(d_kinetic, 2.0) + 2.0 * d_kinetic * h2_mass
			);
			rp = rp.Unit() * recoil_momentum;
			// rebuild beam momentum vector
			ROOT::Math::XYZVector bp = fp[0] + fp[1] + rp;
			// rebuild beam kinetic energy
			double c_kinetic = sqrt(bp.Dot(bp) + pow(c14_mass, 2.0)) - c14_mass;
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

		// fill tree
		opt.Fill();
	}

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