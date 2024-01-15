#include <iostream>

#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>
#include <Math/Vector3D.h>

#include "include/simulate_defs.h"
#include "include/event/generate_event.h"
#include "include/event/detect_event.h"

using namespace ribll;

void SimpleFit(const double *x, double *y, double &k, double &b) {
	int n = 3;
	double sumx = 0.0;
	double sumy = 0.0;
	double sumxy = 0.0;
	double sumx2 = 0.0;
	for (int i = 0; i < n; ++i) {
		sumx += x[i];
		sumy += y[i];
		sumxy += x[i] * y[i];
		sumx2 += x[i] * x[i];
	}
	k = (sumxy - sumx*sumy/double(n)) / (sumx2 - sumx*sumx/double(n));
	b = (sumy - k*sumx) / double(n);
}


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] method\n"
		"  method       Set tracking method:\n"
		"                 a -- momentum approximation.\n"
		"                 p -- average beam direction.\n"
		"                 ad -- momentum approximation with 2H optimized.\n"
		"                 adr -- momentum approximation with 2H optimized\n"
		"                    and relative effect.\n"
		"                 ai -- momentum approximation iteration.\n"
		"Options:\n"
		"  -h             Print this help information.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help
) {
	// initialize
	help = false;
	// start index of positional arugments
	int result = 0;
	for (result = 1; result < argc; ++result) {
		// assumed that all options have read
		if (argv[result][0] != '-') break;
		// short option contains only one letter
		if (argv[result][2] != 0) return -result;
		if (argv[result][1] == 'h') {
			help = true;
			return result;
		} else {
			return -result;
		}
	}
	return result;
}


/// @brief single PPAC approximate tracking
/// @param[in] detect simulated detect event
/// @param[out] sptx reaction point x result
/// @param[out] spty reaction point y result
///
void AppproximateTrack(
	const DetectEvent &detect,
	double *sptx,
	double *spty
) {
	for (int i = 0; i < 3; ++i) {
		// x direction
		// 10Be parameter
		double a_be = sqrt(10.0 * detect.be_kinetic) / 100.0;
		// 4He x parameter
		double a_he = sqrt(4.0 * detect.he_kinetic) / 100.0;
		// 2H x parameter
		double a_d = sqrt(2.0 * detect.d_kinetic) / 135.0;
		// 14C x parameter
		double a_c_x = -sqrt(14.0 * 385.0) / ppac_xz[i];
		// calculate
		double numerator_x = a_be * detect.t0x[0][0];
		numerator_x += a_he * detect.t0x[1][0];
		numerator_x += a_d * detect.tafx;
		numerator_x += a_c_x * detect.ppacx[i];
		double denominator_x = a_be + a_he + a_d + a_c_x;
		sptx[i] = numerator_x / denominator_x;

		// y direction
		// 14C y parameter
		double a_c_y = -sqrt(14.0 * 385.0) / ppac_yz[i];
		// calculate
		double numerator_y = a_be * detect.t0y[0][0];
		numerator_y += a_he * detect.t0y[1][0];
		numerator_y += a_d * detect.tafy;
		numerator_y += a_c_y * detect.ppacy[i];
		double denominator_y = a_be + a_he + a_d + a_c_y;
		spty[i] = numerator_y / denominator_y;
	}
}


/// @brief single PPAC approximate tracking with deutron optimized
/// @param[in] detect simulated detect event
/// @param[out] sptx reaction point x result
/// @param[out] spty reaction point y result
///
void AppproximateTrackDeutron(
	const DetectEvent &detect,
	double *sptx,
	double *spty
) {
	for (int i = 0; i < 3; ++i) {
		// x direction
		// 10Be parameter
		double a_be = sqrt(10.0 * detect.be_kinetic) / 100.0;
		// 4He x parameter
		double a_he = sqrt(4.0 * detect.he_kinetic) / 100.0;
		// 2H x parameter
		double a_d =
			sqrt(2.0 * detect.d_kinetic)
			/ sqrt(
				135.0 * 135.0
				+ detect.tafx * detect.tafx
				+ detect.tafy * detect.tafy
			);
		// 14C x parameter
		double a_c_x = -sqrt(14.0 * 385.0) / ppac_xz[i];
		// calculate
		double numerator_x = a_be * detect.t0x[0][0];
		numerator_x += a_he * detect.t0x[1][0];
		numerator_x += a_d * detect.tafx;
		numerator_x += a_c_x * detect.ppacx[i];
		double denominator_x = a_be + a_he + a_d + a_c_x;
		sptx[i] = numerator_x / denominator_x;

		// y direction
		// 14C y parameter
		double a_c_y = -sqrt(14.0 * 385.0) / ppac_yz[i];
		// calculate
		double numerator_y = a_be * detect.t0y[0][0];
		numerator_y += a_he * detect.t0y[1][0];
		numerator_y += a_d * detect.tafy;
		numerator_y += a_c_y * detect.ppacy[i];
		double denominator_y = a_be + a_he + a_d + a_c_y;
		spty[i] = numerator_y / denominator_y;
	}
}


/// @brief single PPAC approximate tracking with deutron optimized
///		and relative effect
/// @param[in] detect simulated detect event
/// @param[out] sptx reaction point x result
/// @param[out] spty reaction point y result
///
void AppproximateTrackDeutronRelative(
	const DetectEvent &detect,
	double *sptx,
	double *spty
) {
	for (int i = 0; i < 3; ++i) {
		// x direction
		// 10Be parameter
		double a_be = sqrt(
			(2.0 * be10_mass + detect.be_kinetic) * detect.be_kinetic
		) / 100.0;
		// 4He x parameter
		double a_he = sqrt(
			(2.0 * he4_mass + detect.he_kinetic) * detect.he_kinetic
		) / 100.0;
		// 2H x parameter
		double a_d = sqrt(
			(2.0 * h2_mass + detect.d_kinetic) * detect.d_kinetic
		) / sqrt(
			135.0 * 135.0
			+ detect.tafx * detect.tafx
			+ detect.tafy * detect.tafy
		);
		// 14C x parameter
		double a_c_x = -sqrt(
			(2.0 * c14_mass + 385.0) * 385.0
		) / ppac_xz[i];
		// calculate
		double numerator_x = a_be * detect.t0x[0][0];
		numerator_x += a_he * detect.t0x[1][0];
		numerator_x += a_d * detect.tafx;
		numerator_x += a_c_x * detect.ppacx[i];
		double denominator_x = a_be + a_he + a_d + a_c_x;
		sptx[i] = numerator_x / denominator_x;

		// y direction
		// 14C y parameter
		double a_c_y = -sqrt(
			(2.0 * c14_mass + 385.0) * 385.0
		) / ppac_yz[i];
		// calculate
		double numerator_y = a_be * detect.t0y[0][0];
		numerator_y += a_he * detect.t0y[1][0];
		numerator_y += a_d * detect.tafy;
		numerator_y += a_c_y * detect.ppacy[i];
		double denominator_y = a_be + a_he + a_d + a_c_y;
		spty[i] = numerator_y / denominator_y;
	}
}


/// @brief single PPAC approximate tracking with deutron optimized
///		and relative effect in itertaion mode
/// @param[in] detect simulated detect event
/// @param[out] sptx reaction point x result
/// @param[out] spty reaction point y result
/// @param[out] iterations iteration times
///
void AppproximateTrackIteration(
	const DetectEvent &detect,
	double *sptx,
	double *spty,
	int *iterations
) {
	constexpr double ck = 388.0;
	for (int i = 0; i < 3; ++i) iterations[i] = 0;
	for (int i = 0; i < 3; ++i) {
		// x direction
		// 10Be parameter
		double a_be = sqrt(
			(2.0 * be10_mass + detect.be_kinetic) * detect.be_kinetic
		) / 100.0;
		// 4He x parameter
		double a_he = sqrt(
			(2.0 * he4_mass + detect.he_kinetic) * detect.he_kinetic
		) / 100.0;
		// 2H x parameter
		double a_d = sqrt(
			(2.0 * h2_mass + detect.d_kinetic) * detect.d_kinetic
		) / sqrt(
			135.0 * 135.0
			+ detect.tafx * detect.tafx
			+ detect.tafy * detect.tafy
		);
		// 14C x parameter
		double a_c_x = -sqrt(
			(2.0 * c14_mass + ck) * ck
		) / ppac_xz[i];
		// calculate
		double numerator_x = a_be * detect.t0x[0][0];
		numerator_x += a_he * detect.t0x[1][0];
		numerator_x += a_d * detect.tafx;
		numerator_x += a_c_x * detect.ppacx[i];
		double denominator_x = a_be + a_he + a_d + a_c_x;
		sptx[i] = numerator_x / denominator_x;

		// y direction
		// 14C y parameter
		double a_c_y = -sqrt(
			(2.0 * c14_mass + ck) * ck
		) / ppac_yz[i];
		// calculate
		double numerator_y = a_be * detect.t0y[0][0];
		numerator_y += a_he * detect.t0y[1][0];
		numerator_y += a_d * detect.tafy;
		numerator_y += a_c_y * detect.ppacy[i];
		double denominator_y = a_be + a_he + a_d + a_c_y;
		spty[i] = numerator_y / denominator_y;

		// iteration
		// reaction point correct
		double cx, cy;
		cx = cy = 10.0;
		while ((fabs(cx) > 0.01 || fabs(cy) > 0.01) && iterations[i] < 50) {
			++iterations[i];
			// x direction
			// 10Be parameter
			double a_be = sqrt(
				(2.0 * be10_mass + detect.be_kinetic) * detect.be_kinetic
			) / sqrt(
				100.0 * 100.0
				+ pow(detect.t0x[0][0] - sptx[i], 2.0)
				+ pow(detect.t0y[0][0] - spty[i], 2.0)
			);
			// 4He x parameter
			double a_he = sqrt(
				(2.0 * he4_mass + detect.he_kinetic) * detect.he_kinetic
			) / sqrt(
				100.0 * 100.0
				+ pow(detect.t0x[1][0] - sptx[i], 2.0)
				+ pow(detect.t0y[1][0] - spty[i], 2.0)
			);
			// 2H x parameter
			double a_d = sqrt(
				(2.0 * h2_mass + detect.d_kinetic) * detect.d_kinetic
			) / sqrt(
				135.0 * 135.0
				+ pow(detect.tafx - sptx[i], 2.0)
				+ pow(detect.tafy - spty[i], 2.0)
			);
			// 14C x parameter
			double a_c_x = sqrt(
				(2.0 * c14_mass + ck) * ck
			) / sqrt(
				ppac_xz[i] * ppac_xz[i]
				+ pow(detect.ppacx[i] - sptx[i], 2.0)
				+ pow(detect.ppacy[i] - spty[i], 2.0)
			);
			// calculate
			double numerator_x = a_be * (detect.t0x[0][0] - sptx[i]);
			numerator_x += a_he * (detect.t0x[1][0] - sptx[i]);
			numerator_x += a_d * (detect.tafx - sptx[i]);
			numerator_x += a_c_x * (detect.ppacx[i] - sptx[i]);
			double denominator_x = a_be + a_he + a_d + a_c_x;
			cx = numerator_x / denominator_x;

			// y direction
			// 14C y parameter
			double a_c_y = sqrt(
				(2.0 * c14_mass + ck) * ck
			) / sqrt(
				ppac_yz[i] * ppac_yz[i]
				+ pow(detect.ppacx[i] - sptx[i], 2.0)
				+ pow(detect.ppacy[i] - spty[i], 2.0)
			);
			// calculate
			double numerator_y = a_be * (detect.t0y[0][0] - spty[i]);
			numerator_y += a_he * (detect.t0y[1][0] - spty[i]);
			numerator_y += a_d * (detect.tafy - spty[i]);
			numerator_y += a_c_y * (detect.ppacy[i] - spty[i]);
			double denominator_y = a_be + a_he + a_d + a_c_y;
			cy = numerator_y / denominator_y;

			// correct
			sptx[i] += cx;
			spty[i] += cy;
		}
	}
}


constexpr double project_theta = 0.00790883;
constexpr double project_phi = 0.189434;

void ProjectTrack(
	const DetectEvent &detect,
	double *sptx,
	double *spty
) {
	double xk = tan(project_theta) * cos(project_phi);
	double yk = tan(project_theta) * sin(project_phi);
	for (int i = 0; i < 3; ++i) {
		// x direction
		sptx[i] = detect.ppacx[i] + xk * -ppac_xz[i];
		// y direction
		spty[i] = detect.ppacy[i] + yk * -ppac_yz[i];
	}
}

const std::string methods[] = {
	"a", "p", "ad", "adr", "ai"
};

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cerr << "Error: Too few arguments!\n";
		PrintUsage(argv[0]);
		return -1;
	}

	bool help = false;
	int pos_start = ParseArguments(argc, argv, help);

	// need help
	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}
	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}
	if (pos_start >= argc) {
		// positional arguments less than 1
		std::cerr << "Error: Parameter method not found.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// method
	std::string method = std::string(argv[pos_start]);
	// method valid flag
	bool valid_method = false;
	// check method
	for (const auto &m : methods) {
		if (m == method) {
			valid_method = true;
			break;
		}
	}
	if (!valid_method) {
		std::cerr << "Error: Invalid method " << method << "\n";
		return -1;
	}

	// detected data file name
	TString detect_file_name = TString::Format(
		"%s%sdetect.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// detected data file
	TFile detect_file(detect_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)detect_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< detect_file_name << " failed.\n";
		return -1;
	}
	// add generate friend
	ipt->AddFriend("gen=tree", TString::Format(
		"%s%sgenerate.root", kGenerateDataPath, kSimulateDir
	));
	// detect event
	DetectEvent detect;
	// generate event
	GenerateEvent generate;
	// setup input branches
	detect.SetupInput(ipt);
	generate.SetupInput(ipt, "gen.");

	// output file name
	TString output_file_name = TString::Format(
		"%s%ssingle-ppac-track-%s.root",
		kGenerateDataPath, kSimulateDir, method.c_str()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of dX of reaction point with generated data
	TH1F *hist_gen_dx[3];
	for (int i = 0; i < 3; ++i) {
		hist_gen_dx[i] = new TH1F(
			TString::Format("hgdx%d", i),
			TString::Format("dX with generated data"),
			400, -20, 20
		);
	}
	// histogram of dY of reaction point with generated data
	TH1F *hist_gen_dy[3];
	for (int i = 0; i < 3; ++i) {
		hist_gen_dy[i] = new TH1F(
			TString::Format("hgdy%d", i),
			TString::Format("dY with generated data"),
			400, -20, 20
		);
	}
	// histogram of dX of reaction point with detected data
	TH1F *hist_det_dx[3];
	for (int i = 0; i < 3; ++i) {
		hist_det_dx[i] = new TH1F(
			TString::Format("hddx%d", i),
			TString::Format("dX with detected data"),
			400, -20, 20
		);
	}
	// histogram of dY of reaction point with detected data
	TH1F *hist_det_dy[3];
	for (int i = 0; i < 3; ++i) {
		hist_det_dy[i] = new TH1F(
			TString::Format("hddy%d", i),
			TString::Format("dY with detected data"),
			400, -20, 20
		);
	}
	// graph of double PPAC relation
	TGraph gdptx[3], gdpty[3];
	// histogram of beam direction theta
	TH1F hist_beam_theta("hbtheta", "beam theta", 100, 0, 10.0);
	// histogram of beam direction phi
	TH1F hist_beam_phi("hbphi", "beam phi", 1000, -180.0, 180.0);
	// output tree
	TTree opt("tree", "single ppac track");
	// output data
	double sptx[3], spty[3];
	int iterations[3];
	// setup output branches
	opt.Branch("generate_tx", &generate.target_x, "gtx/D");
	opt.Branch("generate_ty", &generate.target_y, "gty/D");
	opt.Branch("detect_tx", &detect.tx, "dtx/D");
	opt.Branch("detect_ty", &detect.ty, "dty/D");
	opt.Branch("sptx", sptx, "sptx[3]/D");
	opt.Branch("spty", spty, "spty[3]/D");
	opt.Branch("iteration", iterations, "iter[3]/I");

	// sum of beam direction angle(theta)
	double sum_theta = 0.0;
	// sun of beam direction angle(phi)
	double sum_phi = 0.0;
	// total number of valid entires
	int valid_num = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Tracking with single PPAC   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// showing process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get data
		ipt->GetEntry(entry);

		if (method == "a") {
			AppproximateTrack(detect, sptx, spty);
		} else if (method == "ad") {
			AppproximateTrackDeutron(detect, sptx, spty);
		} else if (method == "adr") {
			AppproximateTrackDeutronRelative(detect, sptx, spty);
		} else if (method == "ai") {
			AppproximateTrackIteration(detect, sptx, spty, iterations);
		} else if (method == "p") {
			ProjectTrack(detect, sptx, spty);
		}

		// loop PPAC
		for (int i = 0; i < 3; ++i) {
			// fill to histogram
			hist_gen_dx[i]->Fill(sptx[i] - generate.target_x);
			hist_gen_dy[i]->Fill(spty[i] - generate.target_y);
			hist_det_dx[i]->Fill(sptx[i] - detect.tx);
			hist_det_dy[i]->Fill(spty[i] - detect.ty);
		}

		// fill to tree
		opt.Fill();

		// calculate beam direction
		if (detect.valid != 7) continue;
		// PPAC beam parameters
		double xk, xb, yk, yb;
		// fit and get direction parameters
		SimpleFit(ppac_xz, detect.ppacx, xk, xb);
		SimpleFit(ppac_yz, detect.ppacy, yk, yb);
		// beam direction momentum
		ROOT::Math::XYZVector dir(xk, yk, 1.0);
		dir = dir.Unit();
		// accumulate
		sum_theta += dir.Theta();
		sum_phi += dir.Phi();
		++valid_num;

		// fill angle distribution
		hist_beam_theta.Fill(dir.Theta() / pi * 180.0);
		hist_beam_phi.Fill(dir.Phi() / pi * 180.0);

		// fill double ppac correlation graph
		gdptx[0].AddPoint(
			sptx[1]-generate.target_x, sptx[2]-generate.target_x
		);
		gdptx[1].AddPoint(
			sptx[2]-generate.target_x, sptx[0]-generate.target_x
		);
		gdptx[2].AddPoint(
			sptx[0]-generate.target_x, sptx[1]-generate.target_x
		);
		gdpty[0].AddPoint(
			spty[1]-generate.target_y, spty[2]-generate.target_y
		);
		gdpty[1].AddPoint(
			spty[2]-generate.target_y, spty[0]-generate.target_y
		);
		gdpty[2].AddPoint(
			spty[0]-generate.target_y, spty[1]-generate.target_y
		);
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit double ppac track x
	for (int i = 0; i < 3; ++i) {
		TF1 *f = new TF1(TString::Format("dptxf%d", i), "pol1", -1, 1);
		f->SetParameter(0, 0);
		f->SetParameter(1, 1.0);
		gdptx[i].Fit(f, "RQ+");
		std::cout << "DPT x" << i << ": " << f->GetParameter(1)
			<< " " << f->GetParameter(0) << "\n";
	}
	// fit double ppac track y
	for (int i = 0; i < 3; ++i) {
		TF1 *f = new TF1(TString::Format("dptyf%d", i), "pol1", -1, 1);
		f->SetParameter(0, 0);
		f->SetParameter(1, 1.0);
		gdpty[i].Fit(f, "RQ+");
		std::cout << "DPT y" << i << ": " << f->GetParameter(1)
			<< " " << f->GetParameter(0) << "\n";
	}


	// average direction
	double avg_theta = sum_theta / double(valid_num);
	double avg_phi = sum_phi / double(valid_num);

	std::cout << "Average: "
		<< avg_theta << ", " << avg_phi << "\n";

	// save histograms
	for (int i = 0; i < 3; ++i) hist_gen_dx[i]->Write();
	for (int i = 0; i < 3; ++i) hist_gen_dy[i]->Write();
	for (int i = 0; i < 3; ++i) hist_det_dx[i]->Write();
	for (int i = 0; i < 3; ++i) hist_det_dy[i]->Write();
	for (int i = 0; i < 3; ++i) {
		gdptx[i].Write(TString::Format("gdptx%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		gdpty[i].Write(TString::Format("gdpty%d", i));
	}
	hist_beam_theta.Write();
	hist_beam_phi.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	detect_file.Close();
	return 0;
}