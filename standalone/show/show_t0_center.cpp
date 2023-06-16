#include <iostream>
#include <iomanip>

#include <Math/Vector3D.h>
#include <TDirectoryFile.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>

#include "include/event/t0_event.h"
#include "include/event/particle_type_event.h"
#include "include/event/particle_event.h"
#include "include/statistics/center_statistics.h"

using namespace ribll;

const ROOT::Math::XYZVector d1_center{0.0, 0.0, 100.0};
const ROOT::Math::XYZVector d2_center{0.0, 0.0, 111.76};
const ROOT::Math::XYZVector d3_center{0.0, 0.0, 123.52};

/// @brief print usage of this program
/// @param[in] name program name
///
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run [end_run]\n"
		"  run               Set run number.\n"
		"  end_run           Set the last run, inclusive.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n"
		"Examples:\n"
		"  " << name << " 600        Show center of run 600.\n"
		"  " << name << " 600 700    Show center from run 600 to 700.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag
) {
	// initialize
	help = false;
	trigger_tag.clear();
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
		} else if (argv[result][1] == 't') {
			// option of trigger tag
			// get tag in next argument
			++result;
			// miss arguemnt behind option
			if (result == argc) return -argc;
			trigger_tag = argv[result];
		} else {
			return -result;
		}
	}
	return result;
}

/// @brief search for offset in one run
/// @param[in] run run number
/// @param[in] tag trigger tag
/// @returns 0 if success, -1 otherwise
///
int CalculateOffset(
	unsigned int run,
	const std::string &tag
) {
	// t0 file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// t0 file
	TFile t0_file(t0_file_name, "read");
	// t0 tree
	TTree *tree = (TTree*)t0_file.Get("tree");
	if (!tree) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// particle type file name
	TString type_file_name;
	type_file_name.Form(
		"%s%st0-particle-type-%s%04u.root",
		kGenerateDataPath,
		kParticleIdentifyDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	tree->AddFriend("type=tree", type_file_name);
	// ppac particle file name
	TString ppac_file_name;
	ppac_file_name.Form(
		"%s%sxppac-particle-%s%04u.root",
		kGenerateDataPath,
		kParticleDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	tree->AddFriend("ppac=tree", ppac_file_name);
	// input t0 event
	T0Event t0_event;
	// input particle type event
	ParticleTypeEvent type_event;
	// input ppac event
	ParticleEvent ppac_event;
	// setup branches
	t0_event.SetupInput(tree);
	type_event.SetupInput(tree, "type.");
	ppac_event.SetupInput(tree, "ppac.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-center-%s%04u.root",
		kGenerateDataPath,
		kShowDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// offset of detectors
	TH1F hd[4]{
		TH1F("d2dx", "T0D2 #Deltax", 1000, -10, 10),
		TH1F("d2dy", "T0D2 #Deltay", 1000, -10, 10),
		TH1F("d3dx", "T0D3 #Deltax", 200, -10, 10),
		TH1F("d3dy", "T0D3 #Deltay", 200, -10, 10)
	};
	TH1F hist_cos_theta[3]{
		TH1F("d1d2t", "T0 D1D2 cos#theta", 100, 0.99, 1),
		TH1F("d1d3t", "T0 D1D3 cos#theta", 100, 0.99, 1),
		TH1F("d2d3t", "T0 D2D3 cos#theta", 100, 0.99, 1)
	};
	TH2F hd2dxx("d2dxx", "T0D2 #Deltax VS x", 800, -40, 40, 1000, -10, 10);
	TH2F hd2dxy("d2dxy", "T0D2 #Deltax VS y", 800, -40, 40, 1000, -10, 10);
	TH2F hd2dyx("d2dyx", "T0D2 #Deltay VS x", 800, -40, 40, 1000, -10, 10);
	TH2F hd2dyy("d2dyy", "T0D2 #Deltay VS y", 800, -40, 40, 1000, -10, 10);
	TH2F hd3dxx("d3dxx", "T0D3 #Deltax VS x", 800, -40, 40, 1000, -10, 10);
	TH2F hd3dxy("d3dxy", "T0D3 #Deltax VS y", 800, -40, 40, 1000, -10, 10);
	TH2F hd3dyx("d3dyx", "T0D3 #Deltay VS x", 800, -40, 40, 1000, -10, 10);
	TH2F hd3dyy("d3dyy", "T0D3 #Deltay VS y", 800, -40, 40, 1000, -10, 10);
	TTree opt("tree", "delta x and y of d2d3");
	// recoreded x
	double x[2];
	// recorded y
	double y[2];
	// calculated x
	double cx[2];
	// calculated y
	double cy[2];
	opt.Branch("x", x, "x[2]/D");
	opt.Branch("y", y, "y[2]/D");
	opt.Branch("cx", cx, "cx[2]/D");
	opt.Branch("cy", cy, "cy[2]/D");

	// random number generator
	TRandom3 generator(tree->GetEntries());

	// total number of entries
	long long entries = tree->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling offset of run %u   0%%", run);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		tree->GetEntry(entry);
		// initialize
		for (size_t i = 0; i < 2; ++i) {
			x[i] = y[i] = cx[i] = cy[i] = -1e5;
		}
		for (unsigned short i = 0; i < t0_event.num; ++i) {
			// check t0 flag
			if ((t0_event.flag[i] & 0x3) != 0x3) continue;
			// check type
			if (type_event.mass[i] <= 0 || type_event.charge[i] <= 0) continue;
			// check PPAC tracking
			if (ppac_event.num != 4) continue;
			// in target condition
			if (pow(ppac_event.x[3]+3.3, 2.0)+pow(ppac_event.y[3]-1.0, 2.0) > 225.0) {
				continue;
			}

			// target point position
			ROOT::Math::XYZVector target(
				ppac_event.x[3], ppac_event.y[3], 0.0
			);

			// get T0D1 position
			// T0D1 x status flag
			int t0d1_xflag = (t0_event.status[i]>>2) & 0x3;
			// jump for special case
			if (t0d1_xflag != 0 && t0d1_xflag != 1) continue;
			// T0D1 x
			double d1x = t0_event.x[i][0];
			// continuous single strip event
			if (t0d1_xflag == 0) {
				d1x += generator.Rndm() - 0.5;
			}
			// T0D1 y status flag
			int t0d1_yflag = t0_event.status[i] & 0x3;
			// jump for special case
			if (t0d1_yflag != 0 && t0d1_yflag != 1) continue;
			// T0D1 y
			double d1y = t0_event.y[i][0];
			// continous single strip event
			if (t0d1_yflag == 0) {
				d1y += generator.Rndm() - 0.5;
			}
			// d1 position in lab coordinate
			ROOT::Math::XYZVector d1_pos(d1x, d1y, d1_center.Z());

			// get T0D2 position
			// T0D2 x status flag
			int t0d2_xflag = (t0_event.status[i]>>4) & 0x3;
			// jump for special case
			if (t0d2_xflag != 0 && t0d2_xflag != 1) continue;
			// T0D2 x
			double d2x = t0_event.x[i][1];
			// continuous single strip event
			if (t0d2_xflag == 0) {
				d2x += generator.Rndm() * 2.0 - 1.0;
			}
			// T0D2 y status flag
			int t0d2_yflag = (t0_event.status[i]>>6) & 0x3;
			// jump for special case
			if (t0d2_yflag != 0 && t0d2_yflag != 1) continue;
			// T0D2 y
			double d2y = t0_event.y[i][1];
			//  continouts single strip event
			if (t0d2_yflag == 0) {
				d2y += generator.Rndm() * 2.0 - 1.0;
			}
			// d2 position in lab coordinate
			ROOT::Math::XYZVector d2_pos(d2x, d2y, d2_center.Z());
			// d1 position relate to target point
			ROOT::Math::XYZVector d1_rel_pos = d1_pos - target;
			// d2 posiition relate to target point
			ROOT::Math::XYZVector d2_rel_pos = d2_pos - target;
			// cos(theta) of d1 and d2 relative position
			double d1d2_cos_theta = d1_rel_pos.Dot(d2_rel_pos)
				/ (d1_rel_pos.R() * d2_rel_pos.R());
			// fill to histogram
			hist_cos_theta[0].Fill(d1d2_cos_theta);
			// calculated d2 position from target point and d1 position
			ROOT::Math::XYZVector d2_cal_pos = target;
			d2_cal_pos += d1_rel_pos * (d2_center.Z()/d1_center.Z());
			// fill offsets
			hd[0].Fill((d2_cal_pos - d2_pos).X());
			hd[1].Fill((d2_cal_pos - d2_pos).Y());
			// fill offset correlation
			hd2dxx.Fill(d2_cal_pos.X(), (d2_cal_pos - d2_pos).X());
			hd2dxy.Fill(d2_cal_pos.Y(), (d2_cal_pos - d2_pos).X());
			hd2dyx.Fill(d2_cal_pos.X(), (d2_cal_pos - d2_pos).Y());
			hd2dyy.Fill(d2_cal_pos.Y(), (d2_cal_pos - d2_pos).Y());
			// fill branches
			if (i == 0) {
				x[0] = d2_pos.X();
				y[0] = d2_pos.Y();
				cx[0] = d2_cal_pos.X();
				cy[0] = d2_cal_pos.Y();
			}

			if (t0_event.flag[0] == 0x7) {
				// get T0D3 position
				// T0D3 x status flag
				int t0d3_xflag = (t0_event.status[i]>>8) & 0x3;
				// jump for special case
				if (t0d3_xflag != 0 && t0d3_xflag != 1) continue;
				// T0D3 x
				double d3x = t0_event.x[i][2];
				// continuous single strip event
				if (t0d3_xflag == 0) {
					d3x += generator.Rndm() * 2.0 - 1.0;
				}
				// T0D3 y status flag
				int t0d3_yflag = (t0_event.status[i]>>10) & 0x3;
				// jump for special case
				if (t0d3_yflag != 0 && t0d3_yflag != 1) continue;
				// T0D3 y
				double d3y = t0_event.y[i][2];
				//  continouts single strip event
				if (t0d3_yflag == 0) {
					d3y += generator.Rndm() * 2.0 - 1.0;
				}
				// d3 particle position
				ROOT::Math::XYZVector d3_pos(d3x, d3y, d3_center.Z());
				// d3 position relate to target point
				ROOT::Math::XYZVector d3_rel_pos = d3_pos - target;
				// cos(theta) of d1 and d3 relative position
				double d1d3_cos_theta = d1_rel_pos.Dot(d3_rel_pos)
					/ (d1_rel_pos.R() * d3_rel_pos.R());
				// cso(theta) of d2 and d3 relative position
				double d2d3_cos_theta = d2_rel_pos.Dot(d3_rel_pos)
					/ (d2_rel_pos.R() * d3_rel_pos.R());
				// fill to histgram
				hist_cos_theta[1].Fill(d1d3_cos_theta);
				hist_cos_theta[2].Fill(d2d3_cos_theta);
				// calculated d3 position from target point and d1 position
				ROOT::Math::XYZVector d3_cal_pos = target;
				d3_cal_pos += d1_rel_pos * (d3_center.Z()/d1_center.Z());
				// fill offsets
				hd[2].Fill((d3_cal_pos - d3_pos).X());
				hd[3].Fill((d3_cal_pos - d3_pos).Y());
				// fill offset correlation
				hd3dxx.Fill(d3_cal_pos.X(), (d3_cal_pos - d3_pos).X());
				hd3dxy.Fill(d3_cal_pos.Y(), (d3_cal_pos - d3_pos).X());
				hd3dyx.Fill(d3_cal_pos.X(), (d3_cal_pos - d3_pos).Y());
				hd3dyy.Fill(d3_cal_pos.Y(), (d3_cal_pos - d3_pos).Y());
				// fill branches
				if (i == 0) {
					x[1] = d3_pos.X();
					y[1] = d3_pos.Y();
					cx[1] = d3_cal_pos.X();
					cy[1] = d3_cal_pos.Y();
				}
			}
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// statistics
	CenterStatistics statistics(run, "t0", tag);

	// fit offset
	for (int i = 0; i < 2; ++i) {
		TF1 fx(TString::Format("fx%d", i), "gaus", -2, 2);
		fx.SetParameter(0, 20);
		fx.SetParameter(1, 0.0);
		fx.SetParameter(2, 1.0);
		hd[i*2].Fit(&fx, "QR+");
		statistics.x_offset[i] = fx.GetParameter(1);

		TF1 fy(TString::Format("fy%d", i), "gaus", -2, 2);
		fy.SetParameter(0, 20);
		fy.SetParameter(1, 0.0);
		fy.SetParameter(2, 1.0);
		hd[i*2+1].Fit(&fy, "QR+");
		statistics.y_offset[i] = fy.GetParameter(1);
	}

	// save histograms
	for (size_t i = 0; i < 4; ++i) {
		hd[i].Write();
	}
	for (size_t i = 0; i < 3; ++i) {
		hist_cos_theta[i].Write();
	}
	hd2dxx.Write();
	hd2dxy.Write();
	hd2dyx.Write();
	hd2dyy.Write();
	hd3dxx.Write();
	hd3dxy.Write();
	hd3dyx.Write();
	hd3dyy.Write();
	opt.Write();
	// close files
	opf.Close();
	t0_file.Close();

	statistics.Write();
	statistics.Print();

	return 0;
}


/// @brief show offsets of multiple runs in graph
/// @param[in] run the start run number
/// @param[in] end_run the last run number
/// @param[in] tag trigger tag
/// @returns 0 if success, -1 otherwise
///
int ShowOffsets(
	unsigned int run,
	unsigned int end_run,
	const std::string &tag
) {
	// input csv file name
	TString input_file_name;
	input_file_name.Form(
		"%sstatistics/center.csv", kGenerateDataPath
	);
	// input file
	std::ifstream fin(input_file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open file " << input_file_name << " failed.\n";
		return -1;
	}

	// output root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-center-%s%04u-%04u.root",
		kGenerateDataPath,
		kShowDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run,
		end_run
	);
	// output root file
	TFile opf(output_file_name, "recreate");
	// T0 x offset
	TGraph gx[2];
	// T0 y offset
	TGraph gy[2];
	// T0D2 offsets
	TMultiGraph gd2;
	// T0D3 offsets
	TMultiGraph gd3;
	// T0 offsets
	TMultiGraph gd2d3;

	// add graph to multigraph
	gd2.Add(gx);
	gd2.Add(gy);
	gd3.Add(gx+1);
	gd3.Add(gy+1);
	gd2d3.Add(gx);
	gd2d3.Add(gy);
	gd2d3.Add(gx+1);
	gd2d3.Add(gy+1);

	// colors
	const int colors[2] = {1, 600};
	for (size_t i = 0; i < 2; ++i) {
		// set line color
		gx[i].SetLineColor(colors[0]);
		gy[i].SetLineColor(colors[1]);
		// set marker style
		gx[i].SetMarkerStyle(20+i);
		gy[i].SetMarkerStyle(20+i);
		// set marker color
		gx[i].SetMarkerColor(colors[0]);
		gy[i].SetMarkerColor(colors[1]);
	}
	// make d2 legend
	TLegend *legend_d2 = new TLegend(0.8, 0.8, 0.95, 0.95);
	legend_d2->AddEntry(gx, "x");
	legend_d2->AddEntry(gy, "y");
	gd2.GetListOfFunctions()->Add(legend_d2);
	// make d3 legend
	TLegend *legend_d3 = new TLegend(0.8, 0.8, 0.95, 0.95);
	legend_d3->AddEntry(gx+1, "x");
	legend_d3->AddEntry(gy+1, "y");
	gd3.GetListOfFunctions()->Add(legend_d3);
	// make d2d3 legend
	TLegend *legend_d2d3 = new TLegend(0.8, 0.75, 0.95, 0.95);
	for (size_t i = 0; i < 2; ++i) {
		legend_d2d3->AddEntry(gx+i, TString::Format("d%ldx", i+2));
		legend_d2d3->AddEntry(gy+i, TString::Format("d%ldy", i+2));
	}
	gd2d3.GetListOfFunctions()->Add(legend_d2d3);

	// buffer to read lines
	std::string buffer;
	// read first title line
	std::getline(fin, buffer);
	// statistics entry to read
	CenterStatistics statistics;
	// read lines
	fin >> statistics;
	while (fin.good()) {
		if (
			statistics.Run() >= run
			&& statistics.Run() <= end_run
			&& (
				(statistics.Tag() == "-" && tag.empty())
				|| (statistics.Tag() ==  tag)
			)
			&& statistics.Telescope() == "t0"
		) {
			for (size_t i = 0; i < 2; ++i) {
				gx[i].AddPoint(statistics.Run(), statistics.x_offset[i]);
				gy[i].AddPoint(statistics.Run(), statistics.y_offset[i]);
			}
		}
		fin >> statistics;
	}

	// save graphs
	for (size_t i = 0; i < 2; ++i) {
		gx[i].Write(TString::Format("d%ldx", i+2));
		gy[i].Write(TString::Format("d%ldy", i+2));
	}
	gd2.Write("d2");
	gd3.Write("d3");
	gd2d3.Write("d2d3");
	// close files
	opf.Close();
	return 0;
}


int main(int argc, char **argv) {
	if (argc < 2) {
		PrintUsage(argv[0]);
		return -1;
	}

	// help flag
	bool help = false;
	// trigger tag
	std::string tag;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag);

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
		std::cerr << "Error: Miss run argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	// run number
	unsigned int run = atoi(argv[pos_start]);
	unsigned int end_run = run;
	if (pos_start + 1 < argc) {
		end_run = atoi(argv[pos_start+1]);
	}

	if (end_run == run) {
		if (CalculateOffset(run, tag)) return -1;
	} else {
		if (ShowOffsets(run, end_run, tag)) return -1;
	}

	return 0;
}
