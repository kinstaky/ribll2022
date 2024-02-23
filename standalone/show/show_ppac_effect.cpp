#include <fstream>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/ppac_track.h"

using namespace ribll;

double ThreeBodyProcess(
	const ThreeBodyInfoEvent &event,
	const double target_x,
	const double target_y,
	double &be_kinetic,
	double &he_kinetic,
	double &d_kinetic,
	double &c_kinetic
) {
	// 10Be kinetic energy
	// T0D1 energy
	be_kinetic = t0_param[0][0] + t0_param[0][1] * event.be_channel[0];
	// T0D2 energy
	be_kinetic += t0_param[1][0] + t0_param[1][1] * event.be_channel[1];
	// T0D3 energy
	if (event.layer[0] > 1) {
		be_kinetic += t0_param[2][0] + t0_param[2][1] * event.be_channel[2];
	}
	// T0S1 energy
	if (event.layer[0] > 2) {
		be_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
	}
	// T0S2 energy
	if (event.layer[0] > 3) {
		be_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
	}
	// T0S3 energy
	if (event.layer[0] > 4) {
		be_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
	}

	// 4He kinetic energy
	// T0D1 energy
	he_kinetic = t0_param[0][0] + t0_param[0][1] * event.he_channel[0];
	// T0D2 energy
	he_kinetic += t0_param[1][0] + t0_param[1][1] * event.he_channel[1];
	// T0D3 energy
	if (event.layer[1] > 1) {
		he_kinetic += t0_param[2][0] + t0_param[2][1] * event.he_channel[2];
	}
	// T0S1 energy
	if (event.layer[1] > 2) {
		he_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
	}
	// T0S2 energy
	if (event.layer[1] > 3) {
		he_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
	}
	// T0S3 energy
	if (event.layer[1] > 4) {
		he_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
	}

	double a0 = csi_param[event.csi_index][0];
	double a1 = csi_param[event.csi_index][1];
	double a2 = csi_param[event.csi_index][2];
	// 2H kinetic energy
	d_kinetic = event.tafd_energy + pow(
		(event.csi_channel - a2) / a0,
		1.0 / a1
	);

	// 10Be momentum
	double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - target_x,
		event.be_y[0] - target_y,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - target_x,
		event.he_y[0] - target_y,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		event.d_x - target_x,
		event.d_y - target_y,
		135.0
	);
	p_d = p_d.Unit() * d_momentum;

	// beam 14C momentum vector
	ROOT::Math::XYZVector p_c = p_be + p_he + p_d;
// std::cout << p_be.X() << ", " << p_be.Y() << ", " << p_be.Z() << "\n"
// 	<< p_he.X() << ", " << p_he.Y() << ", " << p_he.Z() << "\n"
// 	<< p_d.X() << ", " << p_d.Y() << ", " << p_d.Z() << "\n"
// 	<< p_c.X() << ", " << p_c.Y() << ", " << p_c.Z() << "\n";

	// 14C momentum
	double c_momentum = p_c.R();
	// 14C kinetic energy
	c_kinetic =
		sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

	// three-fold Q value
	double q = be_kinetic + he_kinetic + d_kinetic - c_kinetic;

	return q;
}


/// @brief single PPAC tracking, approximate method
/// @param[in] be_kinetic kinetic energy of 10Be
/// @param[in] he_kinetic kinetic energy of 4He
/// @param[in] d_kinetic kinetic energy of 2H
/// @param[in] ppac_z distance of PPAC from target
/// @param[in] be_position 10Be position on T0
/// @param[in] he_position 4He position on T0
/// @param[in] d_position 2H position on TAFD
/// @param[in] c_position beam position on PPAC
/// @returns reaction point
///
double ApproximateTrack(
	double be_kinetic,
	double he_kinetic,
	double d_kinetic,
	double ppac_z,
	double be_position,
	double he_position,
	double d_position,
	double c_position
) {
	// 10Be parameter
	double a_be = sqrt(10.0 * be_kinetic) / 100.0;
	// 4He parameter
	double a_he = sqrt(4.0 * he_kinetic) / 100.0;
	// 2H parameter
	double a_d = sqrt(2.0 * d_kinetic) / 135.0;
	// 14C parameter
	double a_c = -sqrt(14.0 * 385.0) / ppac_z;
	// calculate
	double numerator = a_be * be_position;
	numerator += a_he * he_position;
	numerator += a_d * d_position;
	numerator += a_c * c_position;
	double denominator = a_be + a_he + a_d + a_c;
	return numerator / denominator;
}


/// @brief single PPAC tracking, optimize approximate method
/// @param[in] param optimized parameters
/// @param[in] be_kinetic kinetic energy of 10Be
/// @param[in] he_kinetic kinetic energy of 4He
/// @param[in] d_kinetic kinetic energy of 2H
/// @param[in] be_position 10Be position on T0
/// @param[in] he_position 4He position on T0
/// @param[in] d_position 2H position on TAFD
/// @param[in] c_position beam position on PPAC
/// @returns reaction point
///
double OptimizeApproximateTrack(
	double *param,
	double be_kinetic,
	double he_kinetic,
	double d_kinetic,
	double be_position,
	double he_position,
	double d_position,
	double c_position
) {
	// 10Be parameter
	double a_be = sqrt(be_kinetic) * param[0];
	// 4He parameter
	double a_he = sqrt(he_kinetic) * param[1];
	// 2H parameter
	double a_d = sqrt(d_kinetic) * param[2];
	// 14C parameter
	double a_c = -sqrt(385.0) * param[3];
	// calculate
	double numerator = a_be * be_position;
	numerator += a_he * he_position;
	numerator += a_d * d_position;
	numerator += a_c * c_position;
	double denominator = a_be + a_he + a_d + a_c;
	return numerator / denominator;
}


/// @brief single PPAC tracking, optimize deutron approximate method
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
double DeutronApproximateTrack(
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
	double a_be = sqrt(10.0 * be_kinetic) / 100.0;
	// 4He parameter
	double a_he = sqrt(4.0 * he_kinetic) / 100.0;
	// 2H parameter
	double a_d =
		sqrt(2.0 * d_kinetic)
		/ sqrt(
			135.0 * 135.0
			+ d_position * d_position
			+ d_position_v * d_position_v
		);
	// 14C parameter
	double a_c = -sqrt(14.0 * 385.0) / ppac_z;
	// calculate
	double numerator = a_be * be_position;
	numerator += a_he * he_position;
	numerator += a_d * d_position;
	numerator += a_c * c_position;
	double denominator = a_be + a_he + a_d + a_c;
	return numerator / denominator;
}


/// @brief single PPAC approximate tracking with deutron optimized
///		and relative effect in iteration mode
/// @param[in] be_kinetic kinetic energy of 10Be
/// @param[in] he_kinetic kinetic energy of 4He
/// @param[in] d_kinetic kinetic energy of 2H
/// @param[in] ppac_xz distance of PPAC from target
/// @param[in] ppac_yz distance of PPAC from target
/// @param[in] be_position_x 10Be x position on T0
/// @param[in] be_position_y 10Be y position on T0
/// @param[in] he_position_x 4He x position on T0
/// @param[in] he_position_y 4He y position on T0
/// @param[in] d_position_x 2H x position on TAFD
/// @param[in] d_position_y 2H y position on TAFD
/// @param[in] c_position_x beam x position on PPAC
/// @param[in] c_position_y beam y position on PPAC
/// @param[out] sptx reaction point x
/// @param[out] spty reaction point y
/// @param[out] iteration iteration times
///
void ApproximateTrackIteration(
	double be_kinetic,
	double he_kinetic,
	double d_kinetic,
	double ppac_xz,
	double ppac_yz,
	double be_position_x,
	double be_position_y,
	double he_position_x,
	double he_position_y,
	double d_position_x,
	double d_position_y,
	double c_position_x,
	double c_position_y,
	double &sptx,
	double &spty,
	int &iteration
) {
	constexpr double ck = 388.0;
	// x direction
	// 10Be parameter
	double a_be = sqrt(
		(2.0 * mass_10be + be_kinetic) * be_kinetic
	) / 100.0;
	// 4He x parameter
	double a_he = sqrt(
		(2.0 * mass_4he + he_kinetic) * he_kinetic
	) / 100.0;
	// 2H x parameter
	double a_d = sqrt(
		(2.0 * mass_2h + d_kinetic) * d_kinetic
	) / sqrt(
		135.0 * 135.0
		+ d_position_x * d_position_x
		+ d_position_y * d_position_y
	);
	// 14C x parameter
	double a_c_x = -sqrt(
		(2.0 * mass_14c + ck) * ck
	) / ppac_xz;
	// calculate
	double numerator_x = a_be * be_position_x;
	numerator_x += a_he * he_position_x;
	numerator_x += a_d * d_position_x;
	numerator_x += a_c_x * c_position_x;
	double denominator_x = a_be + a_he + a_d + a_c_x;
	sptx = numerator_x / denominator_x;

	// y direction
	// 14C y parameter
	double a_c_y = -sqrt(
		(2.0 * mass_14c + ck) * ck
	) / ppac_yz;
	// calculate
	double numerator_y = a_be * be_position_y;
	numerator_y += a_he * he_position_y;
	numerator_y += a_d * d_position_y;
	numerator_y += a_c_y * c_position_y;
	double denominator_y = a_be + a_he + a_d + a_c_y;
	spty = numerator_y / denominator_y;

	// iteration
	// reaction point correct
	double cx, cy;
	cx = cy = 10.0;
	iteration = 0;
	while ((fabs(cx) > 0.01 || fabs(cy) > 0.01) && iteration < 50) {
		++iteration;
		// x direction
		// 10Be parameter
		double a_be = sqrt(
			(2.0 * mass_10be + be_kinetic) * be_kinetic
		) / sqrt(
			100.0 * 100.0
			+ pow(be_position_x - sptx, 2.0)
			+ pow(be_position_y - spty, 2.0)
		);
		// 4He x parameter
		double a_he = sqrt(
			(2.0 * mass_4he + he_kinetic) * he_kinetic
		) / sqrt(
			100.0 * 100.0
			+ pow(he_position_x - sptx, 2.0)
			+ pow(he_position_y - spty, 2.0)
		);
		// 2H x parameter
		double a_d = sqrt(
			(2.0 * mass_2h + d_kinetic) * d_kinetic
		) / sqrt(
			135.0 * 135.0
			+ pow(d_position_x - sptx, 2.0)
			+ pow(d_position_y - spty, 2.0)
		);
		// 14C x parameter
		double a_c_x = sqrt(
			(2.0 * mass_14c + ck) * ck
		) / sqrt(
			ppac_xz * ppac_xz
			+ pow(c_position_x - sptx, 2.0)
			+ pow(c_position_y - spty, 2.0)
		);
		// calculate
		double numerator_x = a_be * (be_position_x - sptx);
		numerator_x += a_he * (he_position_x - sptx);
		numerator_x += a_d * (d_position_x - sptx);
		numerator_x += a_c_x * (c_position_x - sptx);
		double denominator_x = a_be + a_he + a_d + a_c_x;
		cx = numerator_x / denominator_x;

		// y direction
		// 14C y parameter
		double a_c_y = sqrt(
			(2.0 * mass_14c + ck) * ck
		) / sqrt(
			ppac_yz * ppac_yz
			+ pow(c_position_x - sptx, 2.0)
			+ pow(c_position_y - spty, 2.0)
		);
		// calculate
		double numerator_y = a_be * (be_position_y - spty);
		numerator_y += a_he * (he_position_y - spty);
		numerator_y += a_d * (d_position_y - spty);
		numerator_y += a_c_y * (c_position_y - spty);
		double denominator_y = a_be + a_he + a_d + a_c_y;
		cy = numerator_y / denominator_y;

		// correct
		sptx += cx;
		spty += cy;
	}
}


int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody-only-xppac.root", kGenerateDataPath, kInformationDir
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
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sppac-effect.root", kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of Q value with reaction point from PPAC
	TH1F hist_q_ppac("hqp", "Q with PPAC", 90, -23, -8);
	// histogram of Q value with reaction point at origin
	TH1F hist_q_origin("hqo", "Q with origin reaction point", 90, -23, -8);
	// histogram of Q value with reaction point from T0
	TH1F hist_q_t0("hqt", "Q with T0 reaction point", 60, -23, -8);
	// histograms of Q value with single PPAC tracking approximation (spta)
	TH1F *hist_q_spta[3];
	for (int i = 0; i < 3; ++i) {
		hist_q_spta[i] = new TH1F(
			TString::Format("hqspta%d", i),
			TString::Format("Q with single PPAC(%d) track reaction point", i),
			90, -23, -8
		);
	}
	// histograms of Q value with single PPAC tracking
	// in optimized approximation (sptao)
	TH1F *hist_q_sptao[3];
	for (int i = 0; i < 3; ++i) {
		hist_q_sptao[i] = new TH1F(
			TString::Format("hqsptao%d", i),
			TString::Format("Q with sptao %d method reaction point", i),
			90, -23, -8
		);
	}
	// histograms of Q value with single PPAC trcking
	// in approximation method with deutron optimized
	TH1F *hist_q_sptad[3];
	for (int i = 0; i < 3; ++i) {
		hist_q_sptad[i] = new TH1F(
			TString::Format("hqsptad%d", i),
			TString::Format("Q with sptad %d method reaction point", i),
			90, -23, -8
		);
	}
	// histogram of rebuilt 14C kinetic
	TH1F hist_c_kinetic("hck", "rebuilt 14C kinetic", 100, 350, 400);
	// output tree
	TTree opt("tree", "Q in different condition");
	// output data
	// multiple PPAC tracking target x and y
	double mptx, mpty;
	// T0 inverse tracking target x and y
	double t0tx, t0ty;
	// single PPAC tracking approximate target x and y
	double sptax[3], sptay[3];
	// single PPAC tracking optimized approximate target x and y
	double sptaox[3], sptaoy[3];
	// single PPAC tracking approximate deutron optimized target x and y
	double sptadx[3], sptady[3];
	// single PPAC tracking approximate deutron optimized and relative
	// target x and y
	double sptadrx[3], sptadry[3];
	// single PPAC tracking approximate with iterations
	double sptaix[3], sptaiy[3];
	int iteration[3];
	// flag of single PPAC tracking
	int spt_flag;
	// q value from different tracking condition
	double q_ppac, q_origin, q_t0;
	// q value from single PPAC tracking
	double q_spta[3], q_sptao[3], q_sptad[3], q_sptadr[3], q_sptai[3];
	// setup output branches
	opt.Branch("mptx", &mptx, "mptx/D");
	opt.Branch("mpty", &mpty, "mpty/D");
	opt.Branch("t0tx", &t0tx, "t0tx/D");
	opt.Branch("t0ty", &t0ty, "t0ty/D");
	opt.Branch("spt_flag", &spt_flag, "spt_flag/I");
	opt.Branch("sptax", sptax, "sptax[3]/D");
	opt.Branch("sptay", sptay, "sptay[3]/D");
	opt.Branch("sptaox", sptaox, "sptaox[3]/D");
	opt.Branch("sptaoy", sptaoy, "sptaoy[3]/D");
	opt.Branch("sptadx", sptadx, "sptadx[3]/D");
	opt.Branch("sptady", sptady, "sptady[3]/D");
	opt.Branch("sptadrx", sptadrx, "sptadrx[3]/D");
	opt.Branch("sptadry", sptadry, "sptadry[3]/D");
	opt.Branch("sptaix", sptaix, "sptaix[3]/D");
	opt.Branch("sptaiy", sptaiy, "sptaiy[3]/D");
	opt.Branch("iteration", iteration, "iter[3]/I");
	opt.Branch("q_ppac", &q_ppac, "qp/D");
	opt.Branch("q_origin", &q_origin, "qo/D");
	opt.Branch("q_t0", &q_t0, "qt/D");
	opt.Branch("q_spta", q_spta, "qspta[3]/D");
	opt.Branch("q_sptao", q_sptao, "qsptao[3]/D");
	opt.Branch("q_sptad", q_sptad, "qsptad[3]/D");
	opt.Branch("q_sptadr", q_sptadr, "qsptadr[3]/D");
	opt.Branch("q_sptai", q_sptai, "qsptai[3]/D");

	double be_kinetic, he_kinetic, d_kinetic, c_kinetic;

	// optimized approximation parameters
	double opt_approx_params[6][4];
	// read from file
	// file name
	TString param_file_name = TString::Format(
		"%s%ssingle-ppac-track-a.txt", kGenerateDataPath, kOptimizeDir
	);
	// input file stream
	std::ifstream fin(param_file_name.Data());
	// check file
	if (!fin.good()) {
		std::cerr << "Error: Open file " << param_file_name << " failed.\n";
		return -1;
	}
	// read parameters
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 4; ++j) {
			fin >> opt_approx_params[i][j];
		}
	}
	// close file
	fin.close();

	// total number of entries
	long long entries = ipt->GetEntries();
	for (long long entry = 0; entry < entries; ++entry) {
		// get data
		ipt->GetEntry(entry);

		// calculate reaction point from PPAC data
		// initialize correct parameters
		double ppac_correct[2][3] = {
			{0.0, 0.0, 0.0},
			{0.0, 0.0, 0.0}
		};
		// calculate correct offsets
		PpacOffsetX(ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]);
		PpacOffsetY(ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]);
		// calculate reaction point
		double ppac_cx[3], ppac_cy[3];
		// correct PPAC
		for (int i = 0; i < 3; ++i) {
			ppac_cx[i] = event.xppac_x[i] - ppac_correct[0][i];
			ppac_cy[i] = event.xppac_y[i] - ppac_correct[1][i];
		}
		// slope and intercept
		double xk, yk;
		// fit and get reaction point
		TrackMultiplePpac(event.xppac_xflag, ppac_xz, ppac_cx, xk, mptx);
		TrackMultiplePpac(event.xppac_yflag, ppac_yz, ppac_cy, yk, mpty);
		q_ppac = ThreeBodyProcess(
			event, mptx, mpty,
			be_kinetic, he_kinetic, d_kinetic, c_kinetic
		);
		hist_q_ppac.Fill(q_ppac);

		// calculate reaction point from origin
		q_origin = ThreeBodyProcess(
			event, 0.0, 0.0,
			be_kinetic, he_kinetic, d_kinetic, c_kinetic
		);
		// if (q_ppac < -17.0) hist_q_origin.Fill(q_origin);
		hist_q_origin.Fill(q_origin);

		// calculate reaction point from T0
		double xk1, yk1;
		// int len1 = event.layer[1] == 1 ? 2 : 3;
		if (event.layer[1] >= 2) {
			SimpleFit(t0z, event.he_x, 3, xk1, t0tx);
			SimpleFit(t0z, event.he_y, 3, yk1, t0ty);
			q_t0 = ThreeBodyProcess(
				event, t0tx, t0ty,
				be_kinetic, he_kinetic, d_kinetic, c_kinetic
			);
			// if (q_ppac < -17.0) hist_q_t0.Fill(q_t0);
			hist_q_t0.Fill(q_t0);
		}

		// single PPAC tracking
		spt_flag = 0;
		for (int i = 0; i < 3; ++i) {
			unsigned short flag = 1 << i;
			if ((event.xppac_xflag & event.xppac_yflag & flag) != flag) {
				continue;
			}
			spt_flag |= flag;
			// apprximate method
			sptax[i] = ApproximateTrack(
				be_kinetic, he_kinetic, d_kinetic, ppac_xz[i],
				event.be_x[0], event.he_x[0], event.d_x, ppac_cx[i]
			);
			sptay[i] = ApproximateTrack(
				be_kinetic, he_kinetic, d_kinetic, ppac_yz[i],
				event.be_y[0], event.he_y[0], event.d_y, ppac_cy[i]
			);
			q_spta[i] = ThreeBodyProcess(
				event, sptax[i], sptay[i],
				be_kinetic, he_kinetic, d_kinetic, c_kinetic
			);
			hist_q_spta[i]->Fill(q_spta[i]);

			// optimized approximate method
			sptaox[i] = OptimizeApproximateTrack(
				opt_approx_params[2*i],
				be_kinetic, he_kinetic, d_kinetic,
				event.be_x[0], event.he_x[0], event.d_x, ppac_cx[i]
			);
			sptaoy[i] = OptimizeApproximateTrack(
				opt_approx_params[2*i+1],
				be_kinetic, he_kinetic, d_kinetic,
				event.be_y[0], event.he_y[0], event.d_y, ppac_cy[i]
			);
			q_sptao[i] = ThreeBodyProcess(
				event, sptaox[i], sptaoy[i],
				be_kinetic, he_kinetic, d_kinetic, c_kinetic
			);
			hist_q_sptao[i]->Fill(q_sptao[i]);

			// approximate deutron optmized method
			sptadx[i] = DeutronApproximateTrack(
				be_kinetic, he_kinetic, d_kinetic, ppac_xz[i],
				event.be_x[0], event.he_x[0], event.d_x, event.d_y, ppac_cx[i]
			);
			sptady[i] = DeutronApproximateTrack(
				be_kinetic, he_kinetic, d_kinetic, ppac_yz[i],
				event.be_y[0], event.he_y[0], event.d_y, event.d_x, ppac_cy[i]
			);
			q_sptad[i] = ThreeBodyProcess(
				event, sptadx[i], sptady[i],
				be_kinetic, he_kinetic, d_kinetic, c_kinetic
			);
			hist_q_sptad[i]->Fill(q_sptad[i]);

			// approximate deutron optimized and relative method
			sptadrx[i] = DeutronRelativeApproximateTrack(
				be_kinetic, he_kinetic, d_kinetic, ppac_xz[i],
				event.be_x[0], event.he_x[0], event.d_x, event.d_y, ppac_cx[i]
			);
			sptadry[i] = DeutronRelativeApproximateTrack(
				be_kinetic, he_kinetic, d_kinetic, ppac_yz[i],
				event.be_y[0], event.he_y[0], event.d_y, event.d_x, ppac_cy[i]
			);
			q_sptadr[i] = ThreeBodyProcess(
				event, sptadrx[i], sptadry[i],
				be_kinetic, he_kinetic, d_kinetic, c_kinetic
			);

			// apprximate iteration method
			ApproximateTrackIteration(
				be_kinetic, he_kinetic, d_kinetic,
				ppac_xz[i], ppac_yz[i],
				event.be_x[0], event.be_y[0],
				event.he_x[0], event.he_y[0],
				event.d_x, event.d_y,
				event.xppac_x[i], event.xppac_y[i],
				sptaix[i], sptaiy[i], iteration[i]
			);
			q_sptai[i] = ThreeBodyProcess(
				event, sptaix[i], sptaiy[i],
				be_kinetic, he_kinetic, d_kinetic, c_kinetic
			);
			if (i == 0) hist_c_kinetic.Fill(c_kinetic);
		}

		// fill
		opt.Fill();
	}

	// save histograms
	hist_q_ppac.Write();
	hist_q_origin.Write();
	hist_q_t0.Write();
	for (int i = 0; i < 3; ++i) hist_q_spta[i]->Write();
	for (int i = 0; i < 3; ++i) hist_q_sptao[i]->Write();
	for (int i = 0; i < 3; ++i) hist_q_sptad[i]->Write();
	hist_c_kinetic.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}