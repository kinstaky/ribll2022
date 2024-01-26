#include <iostream>

#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/ppac_track.h"


using namespace ribll;

// double pos_param[6][2] = {
// 	{-3.00068, -0.0174848},
// 	{0.530552, -4.95597},
// 	{-0.453959, 1.39616},
// 	{0.817438, 1.99788},
// 	{1.37519, -1.14332},
// 	{2.08862, 3.909}
// };

double pos_param[6][2] = {
	{-0.42, 0.05},
	{0.59, 0.46},
	{1.22, 0.43},
	{0.92, -0.46},
	{0.69, -3.48},
	{-0.95, 2.11}
};

// double pos_param[6][2] = {
// 	{1.77919, 1.34676},
// 	{1.78847, -2.74042},
// 	{1.2214, 0.427209},
// 	{0.921257, -0.460719},
// 	{0.692692, -3.47688},
// 	{-0.95645, 2.10895}
// };


/// @brief rebuild threebody reaction process
/// @param[in] event input event
/// @param[out] tx reaction point x
/// @param[out] ty reaction point y
/// @param[out] be_kinetic 10Be kinetic energy
/// @param[out] he_kinetic 4He kinetic energy
/// @param[out] d_kinetic 2H kinetic energy
/// @param[out] c_kinetic 14C kinetic energy
/// @returns Q value
///
double ThreeBodyProcess(
	const ThreeBodyInfoEvent &event,
	double &tx,
	double &ty,
	double &be_kinetic,
	double &he_kinetic,
	double &d_kinetic,
	double &c_kinetic
) {
	// 10Be kinematic energy
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

	// 4He kinematic energy
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

	// 2H kinematic energy
	d_kinetic = event.tafd_energy + pow(
		(event.csi_channel - csi_param[event.csi_index][2]) / csi_param[event.csi_index][0],
		1.0 / csi_param[event.csi_index][1]
	);


	// calculate reaction point
	double ppac_cx[3], ppac_cy[3];
	for (int i = 0; i < 3; ++i) {
		ppac_cx[i] = event.ppac_x[i] - ppac_correct[0][i];
		ppac_cy[i] = event.ppac_y[i] - ppac_correct[1][i];
	}
	TrackPpac(
		event.ppac_xflag, ppac_xz, ppac_cx,
		be_kinetic, he_kinetic, d_kinetic,
		event.be_x[0], event.he_x[0], event.d_x, event.d_y,
		tx
	);
	TrackPpac(
		event.ppac_yflag, ppac_yz, ppac_cy,
		be_kinetic, he_kinetic, d_kinetic,
		event.be_y[0], event.he_y[0], event.d_y, event.d_x,
		ty
	);

	// 10Be momentum
	double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - tx,
		event.be_y[0] - ty,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - tx,
		event.he_y[0] - ty,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		event.d_x - tx + pos_param[event.csi_index/2][0],
		event.d_y - ty + pos_param[event.csi_index/2][1],
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
	// 14C kinematic energy
	c_kinetic =
		sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

	// three-fold Q value
	double q = be_kinetic + he_kinetic + d_kinetic - c_kinetic;

	return q;
}


/// @brief calculate 14C exicted energy
/// @param[in] event input event
/// @param[in] be_kinetic 10Be kinetic energy
/// @param[in] he_kinetic 4He kinetic energy
/// @param[in] tx reaction point x
/// @param[in] ty reaction point y
/// @param[in] excited_10be excited energy of 10Be
/// @returns excited energy of 14C
///
double TwoBodyProcess(
	const ThreeBodyInfoEvent &event,
	const double be_kinetic,
	const double he_kinetic,
	const double tx,
	const double ty,
	const double excited_10be
) {
	// 10Be momentum
	double be_momentum = MomentumFromKinetic(
		mass_10be + excited_10be, be_kinetic
	);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - tx,
		event.be_y[0] - ty,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - tx,
		event.he_y[0] - ty,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	// excited 14C momentum vector
	ROOT::Math::XYZVector p_c = p_be + p_he;
	// excited 14C momentum
	double c_momentum = p_c.R();
	// excited 14C total energy
	double c_energy =
		(be_kinetic + mass_10be + excited_10be)
		+ (he_kinetic + mass_4he);
	// excited 14C mass
	double excited_c_mass = sqrt(
		pow(c_energy, 2.0) - pow(c_momentum, 2.0)
	);
	// excited energy of 14C
	double excited_14c = excited_c_mass - mass_14c;

	return excited_14c;
}


class Spectrum {
public:

	Spectrum(int n): n_(n) {}

	double operator()(double *x, double *par) {
		double result = 0.0;
		for (int i = 0; i < n_; ++i) {
			result += par[3*i] * exp(-0.5 * pow((x[0] - par[3*i+1]) / par[3*i+2], 2.0));
		}
		return result;
	}

private:
	int n_;
};


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
	// input data
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// ipt->AddFriend(
	// 	"group=tree",
	// 	TString::Format(
	// 		"%s%sthreebody-group.root",
	// 		kGenerateDataPath,
	// 		kOptimizeDir
	// 	)
	// );
	// double mean;
	// double sigma;
	// ipt->SetBranchAddress("group.q_mean", &mean);
	// ipt->SetBranchAddress("group.q_sigma", &sigma);

	// simulated file to get efficiency
	TString simulate_file_name = TString::Format(
		"%s%sdetect.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// simulate file
	TFile simulate_file(simulate_file_name, "read");
	TGraph *eff0 = (TGraph*)simulate_file.Get("gef0");
	TGraph *eff1 = (TGraph*)simulate_file.Get("gef1");
	TGraph *eff2 = (TGraph*)simulate_file.Get("gef2");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sthreebody.root",
		kGenerateDataPath,
		kSpectrumDir
	);
	// output file
	TFile output_file(output_file_name, "recreate");
	// histogram of Q value
	TH1F hist_q("hq", "Q value", 60, -23, -8);
	hist_q.SetLineColor(kBlack);
	// histogram of Q value seperated by CsI index
	std::vector<TH1F> sep_hist_q;
	for (int i = 0; i < 12; ++i) sep_hist_q.emplace_back(
		TString::Format("hq%d", i), "Q", 60, -23, -8
	);
	// histogram of 14C decay to all state of 10Be
	TH1F c_spec_all("hsca", "spectrum of 14C", 100, 10, 40);
	// spectrum of 14C decay to 10Be ground state
	TH1F c_spec_0("hsc0", "spectrum of 14C to 10Be ground state", 100, 12, 32);
	// spectrum of 14C decay to 10Be 3.5 MeV state
	TH1F c_spec_1("hsc1", "spectrum of 14C to 10Be 3.3MeV state", 100, 12, 32);
	// spectrum of 14C decay to 10Be 6 MeV state
	TH1F c_spec_2("hsc2", "spectrum of 14C to 10Be 6MeV state", 100, 12, 32);
	// spectrum of 14C decay to 10Be 7.5 MeV state
	TH1F c_spec_3("hsc3", "spectrum of 14C to 10Be 7.5MeV state", 100, 12, 32);
	// spectrum in the range of hjx's article
	TH1F hjx_spectrum0("hjx0", "spectrum", 35, 12, 19);
	// sepectrum of first excited state in hjx's article range
	TH1F hjx_spectrum1("hjx1", "spectrum", 45, 16, 25);
	// spectrum of 6MeV excited state in Baba's article predicted range
	TH1F sigma_predicted_spectrum2("baba2", "spectrum", 55, 19, 30);
	hjx_spectrum0.SetLineColor(kBlack);
	hjx_spectrum1.SetLineColor(kBlack);
	sigma_predicted_spectrum2.SetLineColor(kBlack);
	// corrected efficiency
	TGraph gefc[3];
	// output tree
	TTree opt("tree", "spectrum");
	// output data
	double be_kinetic, he_kinetic, d_kinetic, c_kinetic;
	double threebody_q;
	int be_state;
	double c_excited;
	bool ppac_valid;
	// setup output branches
	opt.Branch("be_kinetic", &be_kinetic, "bek/D");
	opt.Branch("he_kinetic", &he_kinetic, "hek/D");
	opt.Branch("d_kinetic", &d_kinetic, "dk/D");
	opt.Branch("c_kinetic", &c_kinetic, "ck/D");
	opt.Branch("q", &threebody_q, "q/D");
	opt.Branch("be_state",  &be_state, "state/I");
	opt.Branch("c_excited", &c_excited, "cex/D");
	// opt.Branch("mean", &mean, "mean/D");
	// opt.Branch("sigma", &sigma, "sigma/D");
	opt.Branch("ppac_valid", &ppac_valid, "pvalid/O");

	constexpr double q_correct[12] = {
		0.5, 0.0, -0.6, -0.6, -0.8, -0.6,
		-0.2, -0.5, 0.0, 0.0, 0.0, -0.3
	};

	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		// check PPAC flags
		ppac_valid = true;
		if (event.ppac_xflag == 0 || event.ppac_yflag == 0) {
			ppac_valid = false;
			opt.Fill();
			continue;
		}

		// reaction point
		double tx, ty;
		threebody_q = ThreeBodyProcess(
			event, tx, ty,
			be_kinetic, he_kinetic, d_kinetic, c_kinetic
		);

		// q value correct
		threebody_q += q_correct[event.csi_index];

		// if (threebody_q > -13.5 && threebody_q < -10) be_state = 0;
		// else if (threebody_q > -16 && threebody_q < -14.5) be_state = 1;
		// else if (threebody_q > -19 && threebody_q < -16.5) be_state = 2;
		// else be_state = -1;

		if (threebody_q > -14 && threebody_q < -10) be_state = 0;
		else if (threebody_q > -16.5 && threebody_q < -14.5) be_state = 1;
		else if (threebody_q > -19 && threebody_q < -17) be_state = 2;
		else be_state = -1;

		double be_excited = 0.0;
		if (be_state == 1) be_excited = 3.368;
		else if (be_state == 2) be_excited = 6.179;

		c_excited = TwoBodyProcess(
			event, be_kinetic, he_kinetic,
			tx, ty, be_excited
		);

		// if (sigma < 0.2) {
			hist_q.Fill(threebody_q);
			sep_hist_q[event.csi_index].Fill(threebody_q);
			c_spec_all.Fill(c_excited);
			if (be_state == 0){
				c_spec_0.Fill(c_excited);
				hjx_spectrum0.Fill(c_excited);
			} else if (be_state == 1) {
				c_spec_1.Fill(c_excited);
				hjx_spectrum1.Fill(c_excited);
			} else if (be_state == 2) {
				c_spec_2.Fill(c_excited);
				sigma_predicted_spectrum2.Fill(c_excited);
			} else if (be_state == 3) {
				c_spec_3.Fill(c_excited);
			}
		// }

		opt.Fill();
	}


	// fit total Q spectrum
	Spectrum total_q_spectrum(4);
	// fit histograms
	TF1 fit_q_spectrum("fq", total_q_spectrum, -23, -8, 12);
	double fit_q_initial_parameters[12] = {
		20.0, -11.5, 1.0,
		60.0, -15.0, 1.0,
		50.0, -17.5, 1.0,
		10.0, -19.0, 0.5
	};
	fit_q_spectrum.SetParameters(fit_q_initial_parameters);
	// set limits to A
	fit_q_spectrum.SetParLimits(0, 0.1, 100.0);
	fit_q_spectrum.SetParLimits(3, 0.1, 100.0);
	fit_q_spectrum.SetParLimits(6, 0.1, 100.0);
	fit_q_spectrum.SetParLimits(9, 1.0, 100.0);
	// set limits to mean
	fit_q_spectrum.SetParLimits(1, -13.0, -11.0);
	fit_q_spectrum.SetParLimits(4, -16.0, -13.5);
	fit_q_spectrum.SetParLimits(7, -18.5, -16.5);
	fit_q_spectrum.SetParLimits(10, -21.0, -19.0);
	// set limits to sigma
	// fit_q_spectrum.SetParLimits(11, 0.1, 10.0);
	fit_q_spectrum.SetNpx(1000);
	hist_q.Fit(&fit_q_spectrum, "R+");

	double q_parameters[12];
	fit_q_spectrum.GetParameters(q_parameters);
	TF1 *fit_q_spectrums[4];
	for (int i = 0; i < 4; ++i) {
		fit_q_spectrums[i] = new TF1(
			TString::Format("fq%d", i), "gaus", -23, -8
		);
		fit_q_spectrums[i]->SetParameters(q_parameters+3*i);
		fit_q_spectrums[i]->SetLineColor(kBlue);
		fit_q_spectrums[i]->SetNpx(200);
		hist_q.GetListOfFunctions()->Add(fit_q_spectrums[i]);
	}

	std::cout << "Fit Q spectrum result:\n";
	for (int i = 0; i < 4; ++i) {
		std::cout << q_parameters[i*3] << " " << q_parameters[i*3+1]
			<< " " << q_parameters[i*3+2] << "\n";
	}


	// fit hjx range spectrum
	Spectrum hjx_peaks(7);
	// fitting function
	TF1 fit_hjx("fhjx", hjx_peaks, 12, 19, 21);
	double fit_hjx_initial_parameters[21] = {
		2.0, 13.6, 0.05,
		5.8, 14.3, 0.08,
		6.0, 14.9, 0.05,
		10.0, 15.7, 0.05,
		5.0, 16.0, 0.1,
		3.0, 16.9, 0.1,
		5.0, 18.0, 0.4
	};

	fit_hjx.SetParameters(fit_hjx_initial_parameters);

	fit_hjx.SetParLimits(0, 1.0, 3.0);
	fit_hjx.SetParLimits(1, 13.4, 14.0);
	fit_hjx.SetParLimits(2, 0.0, 0.2);

	fit_hjx.SetParLimits(3, 2.0, 10.0);
	fit_hjx.SetParLimits(4, 14.1, 14.7);
	fit_hjx.SetParLimits(5, 0.0, 0.2);

	fit_hjx.SetParLimits(6, 2.0, 10.0);
	fit_hjx.SetParLimits(7, 14.8, 15.2);
	fit_hjx.SetParLimits(8, 0.0, 0.4);

	fit_hjx.SetParLimits(9, 2.0, 30.0);
	fit_hjx.SetParLimits(10, 15.3, 15.9);
	fit_hjx.SetParLimits(11, 0.0, 0.2);

	fit_hjx.SetParLimits(12, 0.0, 10.0);
	fit_hjx.SetParLimits(13, 15.9, 16.4);
	fit_hjx.SetParLimits(14, 0.0, 0.25);

	fit_hjx.SetParLimits(15, 0.0, 10.0);
	fit_hjx.SetParLimits(16, 16.5, 17.2);
	fit_hjx.SetParLimits(17, 0.0, 0.15);

	fit_hjx.SetParLimits(18, 0.0, 10.0);
	fit_hjx.SetParLimits(19, 17.6, 18.6);
	fit_hjx.SetParLimits(20, 0.0, 1.0);

	fit_hjx.SetNpx(1000);
	hjx_spectrum0.Fit(&fit_hjx, "R+");
	double hjx_final_parameters[21];
	fit_hjx.GetParameters(hjx_final_parameters);
	TF1 *fit_hjxs[7];
	for (int i = 0; i < 7; ++i) {
		fit_hjxs[i] = new TF1(
			TString::Format("fhjx0%d", i), "gaus", 12, 19
		);
		fit_hjxs[i]->SetParameters(hjx_final_parameters+3*i);
		fit_hjxs[i]->SetLineColor(kBlue);
		fit_hjxs[i]->SetNpx(200);
		hjx_spectrum0.GetListOfFunctions()->Add(fit_hjxs[i]);
	}

	std::cout << "Fit hjx spectrum result:\n";
	for (int i = 0; i < 7; ++i) {
		std::cout << hjx_final_parameters[i*3] << " " << hjx_final_parameters[i*3+1]
			<< " " << hjx_final_parameters[i*3+2] << "\n";
	}



	// first excited state
	Spectrum first_peaks(3);
	// fitting function
	TF1 fit_hjx1("fhjx1", first_peaks, 16, 19, 9);
	double fit_hjx1_initial_parameters[9] = {
		10.0, 17.3, 0.08,
		5.0, 17.9, 0.1,
		15.0, 18.5, 0.1,
	};

	fit_hjx1.SetParameters(fit_hjx1_initial_parameters);

	fit_hjx1.SetParLimits(0, 2.0, 20.0);
	fit_hjx1.SetParLimits(1, 17.0, 17.5);
	fit_hjx1.SetParLimits(2, 0.0, 0.2);

	fit_hjx1.SetParLimits(3, 2.0, 20.0);
	fit_hjx1.SetParLimits(4, 17.2, 18.0);
	fit_hjx1.SetParLimits(5, 0.0, 0.4);

	fit_hjx1.SetParLimits(6, 2.0, 30.0);
	fit_hjx1.SetParLimits(7, 18.0, 19.0);
	fit_hjx1.SetParLimits(8, 0.0, 0.4);


	fit_hjx1.SetNpx(1000);
	hjx_spectrum1.Fit(&fit_hjx1, "R+");
	double hjx1_final_parameters[9];
	fit_hjx1.GetParameters(hjx1_final_parameters);
	TF1 *fit_hjxs1[3];
	for (int i = 0; i < 3; ++i) {
		fit_hjxs1[i] = new TF1(
			TString::Format("fhjx1%d", i), "gaus", 16, 19
		);
		fit_hjxs1[i]->SetParameters(hjx1_final_parameters+3*i);
		fit_hjxs1[i]->SetLineColor(kBlue);
		fit_hjxs1[i]->SetNpx(200);
		hjx_spectrum1.GetListOfFunctions()->Add(fit_hjxs1[i]);
	}

	std::cout << "Fit hjx spectrum first excited result:\n";
	for (int i = 0; i < 3; ++i) {
		std::cout << hjx1_final_parameters[i*3] << " " << hjx1_final_parameters[i*3+1]
			<< " " << hjx1_final_parameters[i*3+2] << "\n";
	}



	// 6MeV excited state
	Spectrum second_peaks(3);
	// fitting function
	TF1 fit_baba2("fbaba2", second_peaks, 21, 30, 9);
	double fit_baba2_initial_parameters[9] = {
		20.0, 21.4, 0.05,
		20.0, 22.2, 0.05,
		10.0, 23.6, 0.2
	};

	fit_baba2.SetParameters(fit_baba2_initial_parameters);

	fit_baba2.SetParLimits(0, 2.0, 50.0);
	fit_baba2.SetParLimits(1, 21.0, 22.0);
	fit_baba2.SetParLimits(2, 0.0, 1.0);

	fit_baba2.SetParLimits(3, 2.0, 50.0);
	fit_baba2.SetParLimits(4, 21.8, 23.0);
	fit_baba2.SetParLimits(5, 0.0, 1.0);

	fit_baba2.SetParLimits(6, 0.0, 50.0);
	fit_baba2.SetParLimits(7, 23.0, 25.0);
	fit_baba2.SetParLimits(8, 0.0, 1.0);


	fit_baba2.SetNpx(1000);
	sigma_predicted_spectrum2.Fit(&fit_baba2, "R+");
	double baba2_final_parameters[9];
	fit_baba2.GetParameters(baba2_final_parameters);
	TF1 *fit_babas2[3];
	for (int i = 0; i < 3; ++i) {
		fit_babas2[i] = new TF1(
			TString::Format("fbaba2%d", i), "gaus", 21, 30
		);
		fit_babas2[i]->SetParameters(baba2_final_parameters+3*i);
		fit_babas2[i]->SetLineColor(kBlue);
		fit_babas2[i]->SetNpx(200);
		sigma_predicted_spectrum2.GetListOfFunctions()->Add(fit_babas2[i]);
	}

	std::cout << "Fit predicted sigma bond spectrum result:\n";
	for (int i = 0; i < 3; ++i) {
		std::cout << baba2_final_parameters[i*3] << " " << baba2_final_parameters[i*3+1]
			<< " " << baba2_final_parameters[i*3+2] << "\n";
	}

	gefc[0].AddPoint(12.0, 0.0);
	gefc[1].AddPoint(15.0, 0.0);
	gefc[2].AddPoint(18.0, 0.0);
	for (int i = 0; i < 41; ++i) {
		gefc[0].AddPoint(12.0+i*0.5, 25.0*eff0->Eval(12.0+i*0.5));
	}
	for (int i = 0; i < 35; ++i) {
		gefc[1].AddPoint(15.0+i*0.5, 40.0*eff1->Eval(15.0+i*0.5));
	}
	for (int i = 0; i < 29; ++i) {
		gefc[2].AddPoint(18.0+i*0.5, 40.0*eff2->Eval(18.0+i*0.5));
	}

	// save
	hist_q.Write();
	for (TH1F &hist : sep_hist_q) hist.Write();
	c_spec_all.Write();
	c_spec_0.Write();
	c_spec_1.Write();
	c_spec_2.Write();
	c_spec_3.Write();
	hjx_spectrum0.Write();
	hjx_spectrum1.Write();
	sigma_predicted_spectrum2.Write();
	gefc[0].Write("gef0");
	gefc[1].Write("gef1");
	gefc[2].Write("gef2");
	// save tree
	opt.Write();
	// close files
	output_file.Close();
	return 0;
}