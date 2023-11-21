#include <cmath>
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMarker.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TTree.h>
#include <TMath.h>
#include <Math/Vector3D.h>

#include <ceres/ceres.h>
#include <glog/logging.h>

#include "include/optimize_utilities.h"
#include "include/defs.h"

using namespace ribll;

constexpr double excitation = 6.0938;
constexpr double beam_kinematic = 430.0;


///	@brief get cos(theta) from energy
/// @param[in] beam_mass mass of beam
/// @param[in] fragment1_mass mass of fragment1
/// @param[in] fragment2_mass mass of fragment2
/// @param[in] q Q value
/// @param[in] beam_kinematic kinematic energy of beam
/// @param[in] fragment1_kinematic kinematic energy of fragment1
/// @returns angle theta of fragment1 and beam
///
double ThetaEnergy(
	double beam_mass,
	double fragment1_mass,
	double fragment2_mass,
	double q,
	double beam_kinematic,
	double fragment1_kinematic
) {
	double fragment2_kinematic = q + beam_kinematic - fragment1_kinematic;

	double momentum_beam = MomentumFromKinematic(beam_mass, beam_kinematic);
	double momentum_fragment1 =
		MomentumFromKinematic(fragment1_mass, fragment1_kinematic);
	double momentum_fragment2 =
		MomentumFromKinematic(fragment2_mass, fragment2_kinematic);

	double cos_theta =
		(
			pow(momentum_fragment1, 2.0)
			+ pow(momentum_beam, 2.0)
			- pow(momentum_fragment2, 2.0)
		)
		/ (2.0 * momentum_fragment1 * momentum_beam);

	return cos_theta;
}


class CostFunctor {
public:

	CostFunctor(
		const int state,
		const double tafd_energy,
		const double csi_channel,
		const double beam_energy,
		const double beam_px,
		const double beam_py,
		const double beam_pz,
		const double tafx,
		const double tafy,
		const unsigned short xflag,
		const unsigned short yflag,
		const double *ppacx,
		const double *ppacy,
		const double csi_power
	)
	: state_(state)
	, tafd_energy_(tafd_energy)
	, csi_channel_(csi_channel)
	, beam_energy_(beam_energy)
	, beam_px_(beam_px)
	, beam_py_(beam_py)
	, beam_pz_(beam_pz)
	, tafx_(tafx)
	, tafy_(tafy)
	, xflag_(xflag)
	, yflag_(yflag)
	, csi_power_(csi_power)
	{
		for (int i = 0; i < 3; ++i) {
			ppacx_[i] = ppacx[i];
			ppacy_[i] = ppacy[i];
		}
	}


	bool operator()(
		const double * const csi_param,
		double *residual
	) const {
		// calculate fragment kinematic energy
		double csi_energy = pow(
			(csi_channel_ - csi_param[1]) / csi_param[0],
			1.0 / csi_power_
		);

		double fragment1_kinematic = tafd_energy_ + csi_energy;

		// calculate kinematic and momentum
		double q = 1.0064550;
		q = state_ == 0 ? q : q - 6.0938;
		// calculated theta
		double calculated_theta = ThetaEnergy(
			mass_15c, mass_2h, mass_14c,
			q, 430.0, fragment1_kinematic
		);

		if (calculated_theta < -1.0 || calculated_theta > 1.0) {
			return false;
		}

		// ppac corrected x and y
		double ppac_x_correct[3], ppac_y_correct[3];
		double ppac_cx[3], ppac_cy[3];

		ppac_x_correct[0] = 0.0;
		ppac_y_correct[0] = 0.0;

		// calculate offset of the other two PPACs
		PpacOffsetX(ppac_x_correct[0], ppac_x_correct[1], ppac_x_correct[2]);
		PpacOffsetY(ppac_y_correct[0], ppac_y_correct[1], ppac_y_correct[2]);

		for (int i = 0; i < 3; ++i) {
			ppac_cx[i] = ppacx_[i] - ppac_x_correct[i];
			ppac_cy[i] = ppacy_[i] - ppac_y_correct[i];
		}

		// calculate the reaction point
		double xk, yk, xb, yb;
		TrackPpac(xflag_, ppac_xz, ppac_cx, xk, xb);
		TrackPpac(yflag_, ppac_yz, ppac_cy, yk, yb);

		// fragment1 direction
		double p1x = tafx_ - xb;
		double p1y = tafy_ - yb;
		double p1z = 135.0;
		double p1 = sqrt(pow(p1x, 2.0) + pow(p1y, 2.0) + pow(p1z, 2.0));
		p1x /= p1;
		p1y /= p1;
		p1z /= p1;
		double theta = p1x * beam_px_ + p1y * beam_py_ + p1z * beam_pz_;

		residual[0] = calculated_theta - theta;

		return true;
	}


private:
	int state_;
	double tafd_energy_;
	double csi_channel_;
	double beam_energy_;
	double beam_px_;
	double beam_py_;
	double beam_pz_;
	double tafx_;
	double tafy_;
	unsigned short xflag_;
	unsigned short yflag_;
	double ppacx_[3];
	double ppacy_[3];
	double csi_power_;
};


int main(int, char **argv) {
	google::InitGoogleLogging(argv[0]);

	// input file name
	TString input_file_name = TString::Format(
		"%s%spd.root", kGenerateDataPath, kInformationDir
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
	// indexes
	int csi_index;
	int state;
	// energy
	double tafd_energy;
	// channel
	double csi_channel;
	double d1_channel;
	double d2_channel;
	// position
	double d1x;
	double d1y;
	double tafx;
	double tafy;
	// target position
	double tx;
	double ty;
	// beam direction
	double bpx, bpy, bpz;
	// strip
	unsigned short tafd_pi_strip;
	// ppac flag
	unsigned short ppac_xflag, ppac_yflag;
	// ppac position
	double ppac_x[3], ppac_y[3];

	// setup input branches
	ipt->SetBranchAddress("csi_index", &csi_index);
	ipt->SetBranchAddress("state", &state);
	ipt->SetBranchAddress("tafd_energy", &tafd_energy);
	ipt->SetBranchAddress("csi_channel", &csi_channel);
	ipt->SetBranchAddress("d1_channel", &d1_channel);
	ipt->SetBranchAddress("d2_channel", &d2_channel);
	ipt->SetBranchAddress("d1_x", &d1x);
	ipt->SetBranchAddress("d1_y", &d1y);
	ipt->SetBranchAddress("rx", &tafx);
	ipt->SetBranchAddress("ry", &tafy);
	ipt->SetBranchAddress("tx", &tx);
	ipt->SetBranchAddress("ty", &ty);
	ipt->SetBranchAddress("beam_px", &bpx);
	ipt->SetBranchAddress("beam_py", &bpy);
	ipt->SetBranchAddress("beam_pz", &bpz);
	ipt->SetBranchAddress("tafd_pi_strip", &tafd_pi_strip);
	ipt->SetBranchAddress("ppac_xflag", &ppac_xflag);
	ipt->SetBranchAddress("ppac_yflag", &ppac_yflag);
	ipt->SetBranchAddress("ppac_x", ppac_x);
	ipt->SetBranchAddress("ppac_y", ppac_y);

	// output file
	TFile opf(
		TString::Format(
			"%s%spd-group.root",
			kGenerateDataPath,
			kOptimizeDir
		),
		"recreate"
	);
	// E-theta graph
	std::vector<TGraph> energy_theta[12];
	// theorical E-theta graph
	TGraph theory_energy_theta[2];
	for (int i = 0; i < 2; ++i) {
		theory_energy_theta[i].SetLineColor(kRed);
	}
	// multigraph
	std::vector<TMultiGraph*> mg[12];
	// missing mass Q value spectrum sepertaed by CsI index
	std::vector<TH1F> sep_miss_q[12];

	// output tree
	// TTree opt("tree", "pd optimization tree");
	// output data
	double taf_energy;
	double csi_energy;
	// angle of beam and fragment 2H
	double angle;
	// angle in lab coordinate
	// double taf_theta;
	// Q value in missing mass method
	double miss_q_value;
	// setup output branches
	// opt.Branch("csi_index", &csi_index, "ci/I");
	// opt.Branch("taf_energy", &taf_energy, "tafe/D");
	// opt.Branch("tafd_energy", &tafd_energy, "tafde/D");
	// opt.Branch("csi_energy", &csi_energy, "csie/D");
	// opt.Branch("csi_channel", &csi_channel, "csic/D");
	// opt.Branch("angle", &angle, "theta/D");
	// opt.Branch("taf_theta", &taf_theta, "taft/D");
	// opt.Branch("q", &miss_q_value, "q/D");

	// double t0_param[4] = {
	// 	0.0876983, 0.00524694,
	// 	0.159542, 0.00644286
	// };

	double csi_param[12][3] = {
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0},
		{1000.0, -2000.0}
	};

	std::vector<double> a0_group[12];
	std::vector<double> a1_group[12];
	std::vector<double> a2_group[12];

	for (double csi_power = 0.64; csi_power < 1.2; csi_power += 0.01) {
		std::cout << "----------------------\nTrying power " << csi_power << "\n";
		// initialize
		for (int i = 0; i < 12; ++i) {
			sep_miss_q[i].emplace_back(
				TString::Format("hq%dg%ld", i, sep_miss_q[i].size()),
				"Q",
				100, -5, 15
			);
			energy_theta[i].emplace_back();
			mg[i].push_back(new TMultiGraph);
		}
		ceres::Problem problem[12];
		// loop to add observations points to problem model
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);

			ceres::CostFunction *cost_function;

			cost_function = new ceres::NumericDiffCostFunction<
				CostFunctor, ceres::CENTRAL, 1, 2
			> (
				new CostFunctor(
					state,
					tafd_energy, csi_channel,
					beam_kinematic,
					bpx, bpy, bpz,
					tafx, tafy,
					ppac_xflag, ppac_yflag,
					ppac_x, ppac_y,
					csi_power
				)
			);
			problem[csi_index].AddResidualBlock(
				cost_function,
				nullptr,
				csi_param[csi_index]
			);
		}

		for (int i = 0; i < 12; ++i) {
			// slove optimization problem
			ceres::Solver::Options options;
			options.max_num_iterations = 100;
			options.linear_solver_type = ceres::DENSE_QR;
			// options.minimizer_progress_to_stdout = true;
			ceres::Solver::Summary summary;
			ceres::Solve(options, problem+i, &summary);
			// std::cout << summary.BriefReport() << "\n";
		}

		// print optimization results
		// std::cout << "Optimization calibration results:\n";
		// for (int i = 0; i < 12; ++i) {
		// 	std::cout << i << " "
		// 		<< csi_param[i][0] << " "
		// 		<< csi_power << " "
		// 		<< csi_param[i][1] << "\n";
		// }

		// store
		for (int i = 0; i < 12; ++i) {
			a0_group[i].push_back(csi_param[i][0]);
			a1_group[i].push_back(csi_power);
			a2_group[i].push_back(csi_param[i][1]);
		}
	}

	for (int i = 0; i < 12; ++i) {
		// save calibration parameters
		std::ofstream fout_calibrate(
			TString::Format(
				"%s%scsi-group-parameters-%d.txt",
				kGenerateDataPath,
				kOptimizeDir,
				i
			).Data()
		);
		for (size_t j = 0; j < a0_group[i].size(); ++j) {
			fout_calibrate
				<< a0_group[i][j] << " "
				<< a1_group[i][j] << " "
				<< a2_group[i][j] << "\n";
		}
		fout_calibrate.close();
	}


	double ppac_correct[2][3] = {
		{0.0, -2.26, -3.41},
		{0.0, 0.95, 1.89}
	};
	PpacOffsetX(ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]);
	PpacOffsetY(ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]);

	// loop to calculate Q value
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		double ppac_cx[3], ppac_cy[3];
		for (int i = 0; i < 3; ++i) {
			ppac_cx[i] = ppac_x[i] - ppac_correct[0][i];
			ppac_cy[i] = ppac_y[i] - ppac_correct[1][i];
		}
		// calcuate target point from PPAC
		double xk, yk, xb, yb;
		TrackPpac(ppac_xflag, ppac_xz, ppac_cx, xk, xb);
		TrackPpac(ppac_yflag, ppac_yz, ppac_cy, yk, yb);
		tx = xb;
		ty = yb;

		// beam momentum
		double beam_momentum =
			MomentumFromKinematic(mass_15c, beam_kinematic);
		// beam momentum direction
		ROOT::Math::XYZVector p2(bpx, bpy, bpz);
		// beam momentum vector
		ROOT::Math::XYZVector bp = p2.Unit() * beam_momentum;

		// calculate angle theta of beam direction and fragment1 direction
		ROOT::Math::XYZVector p1(
			tafx - tx,
			tafy - ty,
			135.0
		);
		p1 = p1.Unit();
		angle = acos(p1.Dot(p2));
		// taf_theta = atan(sqrt(pow(tafx, 2.0) + pow(tafy, 2.0)) / 135.0);


		for (size_t i = 0; i < a0_group[0].size(); ++i) {
			// calculate fragment1 kinematic energy
			csi_energy = pow(
				(csi_channel-a2_group[csi_index][i])/a0_group[csi_index][i],
				1.0/a1_group[csi_index][i]
			);
			taf_energy = tafd_energy + csi_energy;

			// calculate fragment1 momentum
			double momentum1 = MomentumFromKinematic(mass_2h, taf_energy);
			// fragment1 momentum vector
			ROOT::Math::XYZVector f1p = p1 * momentum1;

			// calculated fragment2 momentum
			ROOT::Math::XYZVector cf2p = bp - f1p;
			// calculated fragment2 energy
			double cf2e = sqrt(
				pow(cf2p.R(), 2.0) + pow(mass_14c, 2.0)
			) - mass_14c;
			// Q value from missing mass method
			miss_q_value =
				beam_kinematic + 1.0064550 - taf_energy - cf2e;

			// fill to histograms
			energy_theta[csi_index][i].AddPoint(angle, taf_energy);
			sep_miss_q[csi_index][i].Fill(miss_q_value);
		}


		// opt.Fill();
	}


	for (int i = 0; i < 2; ++i) {
		double q_value = i == 0 ? 1.0064550 : 1.0064550 - 6.0938;
		for (double energy = 0.0; energy < 70.0; energy += 0.5) {
			double cos_theta = ThetaEnergy(
				mass_15c, mass_2h, mass_14c,
				q_value, 430.0, energy
			);
			if (cos_theta > -1.0 && cos_theta < 1.0) {
				theory_energy_theta[i].AddPoint(acos(cos_theta), energy);
			}
		}
	}

	for (int i = 0; i < 12; ++i) {
		for (size_t j = 0; j < mg[i].size(); ++j) {
			mg[i][j]->Add(&(energy_theta[i][j]), "AP*");
			mg[i][j]->Add(theory_energy_theta, "C");
			mg[i][j]->Add(theory_energy_theta+1, "C");
		}
	}

	// save histograms
	opf.cd();
	theory_energy_theta[0].Write("gt0");
	theory_energy_theta[1].Write("gt1");
	// for (int i = 0; i < 12; ++i) {
	// 	for (size_t j = 0; j < energy_theta[i].size(); ++j) {
	// 		energy_theta[i][j].Write(TString::Format("g%dg%ld", i, j));
	// 	}
	// }
	for (int i = 0; i < 12; ++i) {
		for (size_t j = 0; j < mg[i].size(); ++j) {
			mg[i][j]->Write(TString::Format("mg%dg%ld", i, j));
		}
	}
	for (int i = 0; i < 12; ++i) {
		for (auto &hist : sep_miss_q[i]) hist.Write();
	}
	// save tree
	// opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}
