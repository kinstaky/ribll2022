#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TF1.h>
#include <Math/Vector3D.h>

#include "include/event/detect_event.h"
#include "include/simulate_defs.h"

using namespace ribll;


constexpr int strips = 8;
constexpr int param_counts[strips] = {
	3, 3, 3, 3, 3, 3, 3, 3
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
		"%s%sdetect.root",
		kGenerateDataPath, kSimulateDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	// detect event
	DetectEvent detect;
	// setup input branches
	detect.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sresist_strip.root", kGenerateDataPath, kSimulateDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of Q value
	TH1F hist_q("hq", "Q", 300, -23, -8);
	// histogram of beam kinetic
	TH1F hist_c_kinetic("hck", "beam kinetic energy", 500, 350, 400);

	// rebuild Q value
	double rebuild_q;
	// parameters
	int strip_params[strips];


	ipt->GetEntry(1);
	// initialize parameters
	for (int i = 0; i < strips; ++i) strip_params[i] = 0;
	// loop
	do {
		// for (int i = 0; i < strips; ++i) {
		// 	std::cout << strip_params[i] << " \n"[i==strips-1];
		// }

		// change positions
		double ct0x[2], ct0y[2], ctafx, ctafy, cppacx, cppacy;
		ct0x[0] = detect.t0x[0][0] + (strip_params[0]-1)*0.7;
		ct0x[1] = detect.t0x[0][1] + (strip_params[1]-1)*0.7;
		ctafx = detect.tafx + (strip_params[2]-1);
		cppacx = detect.ppacx[0] + (strip_params[3]-1) * 0.7;
		ct0y[0] = detect.t0y[0][0] + (strip_params[4]-1)*0.7;
		ct0y[1] = detect.t0y[0][1] + (strip_params[5]-1)*0.7;
		ctafy = detect.tafy + (strip_params[6]-1);
		cppacy = detect.ppacy[0] + (strip_params[7]-1)*0.7;

		// reaction point from single PPAC track
		double tx = DeutronRelativeApproximateTrack(
			detect.be_kinetic, detect.he_kinetic, detect.d_kinetic,
			ppac_xz[0], ct0x[0], ct0x[1], ctafx, ctafy, cppacx
		);
		double ty = DeutronRelativeApproximateTrack(
			detect.be_kinetic, detect.he_kinetic, detect.d_kinetic,
			ppac_yz[0], ct0y[0], ct0y[1], ctafy, ctafx, cppacy
		);

		// momentum vector of fragments
		ROOT::Math::XYZVector fp[2];
		for (size_t i = 0; i < 2; ++i) {
			fp[i] = ROOT::Math::XYZVector(
				ct0x[i] - tx,
				ct0y[i] - ty,
				detect.t0z[i][0]
			);
			// fragment mass
			double mass = i == 0 ? be10_mass : he4_mass;
			// fragment kinetic energy
			double kinetic = i == 0 ?
				detect.be_kinetic : detect.he_kinetic;
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
			detect.tafz
		);
		// recoild 2H momentum value
		double recoil_momentum = sqrt(
			pow(detect.d_kinetic, 2.0) + 2.0 * detect.d_kinetic * h2_mass
		);
		rp = rp.Unit() * recoil_momentum;
		// rebuild beam momentum vector
		ROOT::Math::XYZVector bp = fp[0] + fp[1] + rp;
		// rebuild beam kinetic energy
		double c_kinetic = sqrt(bp.Dot(bp) + pow(c14_mass, 2.0)) - c14_mass;
		// Q value
		rebuild_q = detect.be_kinetic + detect.he_kinetic
			+ detect.d_kinetic - c_kinetic;

		hist_q.Fill(rebuild_q);
		hist_c_kinetic.Fill(c_kinetic);

	} while (NextGroup(strip_params));

	// save histograms
	hist_q.Write();
	hist_c_kinetic.Write();
	// close files
	ipf.Close();
	opf.Close();
	return 0;
}