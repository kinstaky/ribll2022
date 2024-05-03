#include <iostream>
#include <fstream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>

#include "include/event/threebody_info_event.h"

using namespace ribll;

const double d2_time_threshold = 40.0;
const double d3_time_threshold = 40.0;
const double d3_adj_time_threshold = 100.0;

int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody.root", kGenerateDataPath, kInformationDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	// input event
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output ROOT file name
	TString output_file_name = TString::Format(
		"%s%scheck-threebody-time.root", kGenerateDataPath, kInformationDir 
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of T0D2 time difference
	TH1F hist_d2_diff(
		"hd2dt", "T0D2 time difference", 100, 0, 500
	);
	// hitogram of T0D2 adjacent strip time difference
	TH1F hist_d2_adj_diff(
		"hd2adt", "T0D2 ajdacent strip time difference", 100, 0, 500
	);
	// histogram of T0D3 time difference
	TH1F hist_d3_diff(
		"hd3dt", "T0D3 time difference", 100, 0, 500
	);
	// histogram of T0D3 adjacent strip time difference
	TH1F hist_d3_adj_diff(
		"hd3adt", "T0D3 adjacent strip time difference", 100, 0, 500
	);
	// histogram of T0D1 time difference, different particle and side
	// 0: BeX, 1: BeY(except module 2), 2: HeX, 3: HeY(except module 2)
	TH1F *hist_d1_diff[4];
	hist_d1_diff[0] = new TH1F(
		"hd1dt0", "T0D1 BeX time difference", 100, -200, 200
	);
	hist_d1_diff[1] = new TH1F(
		"hd1dt1", "T0D1 BeY time difference(without module 2)", 100, -200, 200
	);
	hist_d1_diff[2] = new TH1F(
		"hd1dt2", "T0D1 HeX time difference", 100, -200, 200
	);
	hist_d1_diff[3] = new TH1F(
		"hd1dt3", "T0D1 HeY time difference(without module 2)", 100, -200, 200
	);
	// histogram of T0D1 module 2 time difference, different particle
	// 0: BeY, 1: HeY
	TH1F *hist_d1_diff_m2[2];
	hist_d1_diff_m2[0] = new TH1F(
		"hd1dt1m2", "T0D1 BeY time difference(module 2)", 100, -200, 200
	);
	hist_d1_diff_m2[1] = new TH1F(
		"hd1dt3m2", "T0D1 HeY time difference(module 2)", 100, -200, 200
	);
	// histogram of T0D1 adjacent strip time difference, different particle and side
	// 0: BeX, 1: BeY(except module 2), 2: HeX, 3: HeY(except module 2)
	TH1F *hist_d1_adj_diff[4];
	hist_d1_adj_diff[0] = new TH1F(
		"hd1adt0", "T0D1 BeX adjacent strip time difference", 100, -200, 200
	);
	hist_d1_adj_diff[1] = new TH1F(
		"hd1adt1", "T0D1 BeY adjacent strip time difference(without module 2)",
		100, -200, 200
	);
	hist_d1_adj_diff[2] = new TH1F(
		"hd1adt2", "T0D1 HeX adjacent strip time difference", 100, -200, 200
	);
	hist_d1_adj_diff[3] = new TH1F(
		"hd1adt3", "T0D1 HeY adjacent strip time difference(without module 2)",
		100, -200, 200
	);
	// histogram of T0D1 module 2 adjacent strip time difference, different particle 
	// 0: BeY, 1: HeY
	TH1F *hist_d1_adj_diff_m2[2];
	hist_d1_adj_diff_m2[0] = new TH1F(
		"hd1adt1m2", "T0D1 BeX adjacent strip time difference(module 2)",
		100, -200, 200
	);
	hist_d1_adj_diff_m2[1] = new TH1F(
		"hd1adt3m2", "T0D1 HeX adjacent strip time difference(module 2)",
		100, -200, 200
	);

	// output tree
	TTree opt("tree", "check threebody time");
	// output data
	bool valid;
	double d2_avg_time;
	double d2_diff_time[4];
	double d2_adj_diff_time[4];
	double d2_adj_ratio[4];
	bool d2_main_valid, d2_adj_valid, d2_valid;
	double d3_diff_time[4];
	double d3_adj_diff_time[4];
	double d3_adj_ratio[4];
	bool d3_main_valid, d3_adj_valid, d3_valid;
	double d1_diff_time[4];
	int d1_module[4];
	double d1_adj_diff_time[4];
	double d1_adj_ratio[4];
	int d1_adj_module[4];
	bool d1_main_valid, d1_adj_valid, d1_valid;
	double ssd_diff_time[3];
	bool ssd_valid;
	int taf_index;
	double tafd_diff_time[2];
	bool tafd_valid;
	// setup output branches
	opt.Branch("d2_avg_time", &d2_avg_time, "d2avgt/D");
	opt.Branch("d2_diff_time", d2_diff_time, "d2dt[4]/D");
	opt.Branch("d2_adj_diff_time", d2_adj_diff_time, "d2adt[4]/D");
	opt.Branch("d2_ratio", d2_adj_ratio, "d2r[4]/D");
	opt.Branch("d2_main_valid", &d2_main_valid, "d2mv/O");
	opt.Branch("d2_adj_valid", &d2_adj_valid, "d2av/O");
	opt.Branch("d2_valid", &d2_valid, "d2v/O");
	opt.Branch("d3_diff_time", d3_diff_time, "d3dt[4]/D");
	opt.Branch("d3_adj_diff_time", d3_adj_diff_time, "d3adt[4]/D");
	opt.Branch("d3_ratio", d3_adj_ratio, "d3r[4]/D");
	opt.Branch("d3_main_valid", &d3_main_valid, "d3mv/O");
	opt.Branch("d3_adj_valid", &d3_adj_valid, "d3av/O");
	opt.Branch("d3_valid", &d3_valid, "d3v/O");
	opt.Branch("d1_module", d1_module, "d1m[4]/I");
	opt.Branch("d1_diff_time", d1_diff_time, "d1dt[4]/D");
	opt.Branch("d1_adj_diff_time", d1_adj_diff_time, "d1adt[4]/D");
	opt.Branch("d1_ratio", d1_adj_ratio, "d1r[4]/D");
	opt.Branch("d1_adj_module", d1_adj_module, "d1am[4]/I");
	opt.Branch("d1_main_valid", &d1_main_valid, "d1mv/O");
	opt.Branch("d1_adj_valid", &d1_adj_valid, "d1av/O");
	opt.Branch("d1_valid", &d1_valid, "d1v/O");
	opt.Branch("valid", &valid, "v/O");
	opt.Branch("ssd_diff_time", ssd_diff_time, "sdt[3]/D");
	opt.Branch("ssd_valid", &ssd_valid, "sv/O");
	opt.Branch("taf_index", &taf_index, "ti/I");
	opt.Branch("tafd_diff_time", tafd_diff_time, "tafddt[2]/D");
	opt.Branch("tafd_valid", &tafd_valid, "tafdv/O");

	// statistics
	int total_count = int(ipt->GetEntries());
	int d1_valid_count = 0;
	int d2_valid_count = 0;
	int d3_valid_count = 0;
	int ssd_valid_count = 0;
	int tafd_valid_count = 0;
	int valid_count = 0;

	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		// time
		double d2_times[4] {
			event.be_x_time[1][0],
			event.be_y_time[1][0],
			event.he_x_time[1][0],
			event.he_y_time[1][0]
		};
		double d2_adj_times[4] {
			event.be_x_hit[1] == 2 ? event.be_x_time[1][1] : -1e5,
			event.be_y_hit[1] == 2 ? event.be_y_time[1][1] : -1e5,
			event.he_x_hit[1] == 2 ? event.he_x_time[1][1] : -1e5,
			event.he_y_hit[1] == 2 ? event.he_y_time[1][1] : -1e5
		};
		d2_adj_ratio[0] = event.be_x_channel[1][1]
			/ (event.be_x_channel[1][0] + event.be_x_channel[1][1]);
		d2_adj_ratio[1] = event.be_y_channel[1][1]
			/ (event.be_y_channel[1][0] + event.be_y_channel[1][1]);
		d2_adj_ratio[2] = event.he_x_channel[1][1]
			/ (event.he_x_channel[1][0] + event.he_x_channel[1][1]);
		d2_adj_ratio[3] = event.he_y_channel[1][1]
			/ (event.he_y_channel[1][0] + event.he_y_channel[1][1]);

		// sum
		double d2_sum_time = 0.0;
		for (const double &t : d2_times) d2_sum_time += t;
		// average D2 time
		d2_avg_time = d2_sum_time / 4.0;
		for (int i = 0; i < 4; ++i)  {
			d2_diff_time[i] = fabs(d2_times[i] - d2_avg_time);
			hist_d2_diff.Fill(d2_diff_time[i]);
		}
		// check 
		d2_main_valid = true;
		for (int i = 0; i < 4; ++i) {
			if (d2_diff_time[i] > d2_time_threshold) d2_main_valid = false;
		}
		// jump if T0D2 time check failed
		if (!d2_main_valid) {
			d2_adj_valid = false;
			d2_valid = false;
			opt.Fill();
			continue;
		}
		// check adjacent strips
		d2_adj_valid = true;
		for (int i = 0; i < 4; ++i) {
			if (d2_adj_times[i] > -9e4) {
				d2_adj_diff_time[i] = fabs(d2_adj_times[i] - d2_avg_time);
				if (d2_adj_diff_time[i] > d2_time_threshold) {
					d2_adj_valid = false;
				}
				hist_d2_adj_diff.Fill(d2_adj_diff_time[i]);
			} else {
				d2_adj_diff_time[i] = -1e5;
				d2_adj_ratio[i] = 1.0;
			}
		}

		d2_valid = d2_main_valid && d2_adj_valid;
		if (!d2_valid) {
			opt.Fill();
			continue;
		}


		// initialize T0D3 time
		for (int i = 0; i < 4; ++i) {
			d3_diff_time[i] = -1e5;
			d3_adj_diff_time[i] = -1e5;
			d3_adj_ratio[i] = 1.0;
		}
		// get 10Be T0D3 time
		if (event.layer[0] == 2) {
			d3_diff_time[0] = fabs(event.be_x_time[2][0] - d2_avg_time);
			d3_diff_time[1] = fabs(event.be_y_time[2][0] - d2_avg_time);
			hist_d3_diff.Fill(d3_diff_time[0]);
			hist_d3_diff.Fill(d3_diff_time[1]);
			// Be D3X adjacent strip
			if (event.be_x_hit[2] == 2) {
				d3_adj_diff_time[0] = fabs(event.be_x_time[2][1]-d2_avg_time);
				d3_adj_ratio[0] = event.be_x_channel[2][1]
					/ (event.be_x_channel[2][0] + event.be_x_channel[2][1]);
				hist_d3_adj_diff.Fill(d3_adj_diff_time[0]);
			}
			// Be D3Y adjacent strip
			if (event.be_y_hit[2] == 2) {
				d3_adj_diff_time[1] = fabs(event.be_y_time[2][1]-d2_avg_time);
				d3_adj_ratio[1] = event.be_y_channel[2][1]
					/ (event.be_y_channel[2][0] + event.be_y_channel[2][1]);
				hist_d3_adj_diff.Fill(d3_adj_diff_time[1]);
			}
		}
		// get 4He T0D3 time
		if (event.layer[1] >= 2) {
			d3_diff_time[2] = fabs(event.he_x_time[2][0] - d2_avg_time);
			d3_diff_time[3] = fabs(event.he_y_time[2][0] - d2_avg_time);
			hist_d3_diff.Fill(d3_diff_time[0]);
			hist_d3_diff.Fill(d3_diff_time[1]);
			// He D3X adjacent strip
			if (event.he_x_hit[2] == 2) {
				d3_adj_diff_time[2] = fabs(event.he_x_time[2][1]-d2_avg_time);
				d3_adj_ratio[2] = event.he_x_channel[2][1]
					/ (event.he_x_channel[2][0] + event.he_x_channel[2][1]);
				hist_d3_adj_diff.Fill(d3_adj_diff_time[2]);
			}
			// He D3Y adjacent strip
			if (event.he_y_hit[2] == 2) {
				d3_adj_diff_time[3] = fabs(event.he_y_time[2][1]-d2_avg_time);
				d3_adj_ratio[3] = event.he_y_channel[2][1]
					/ (event.he_y_channel[2][0] + event.he_y_channel[2][1]);
				hist_d3_adj_diff.Fill(d3_adj_diff_time[3]);
			}
		}
		// check T0D3 time now
		d3_main_valid = d3_adj_valid = d3_valid = true;
		for (int i = 0; i < 4; ++i) {
			if (d3_diff_time[i] > d3_time_threshold) d3_main_valid = false; 
			if (d3_adj_diff_time[i] > d3_adj_time_threshold) {
				d3_adj_valid = false;
			}
		}
		d3_valid = d3_main_valid && d3_adj_valid;


		// initialize T0D1 time
		for (int i = 0; i < 4; ++i) {
			d1_diff_time[i] = -1e5;
			d1_adj_diff_time[i] = -1e5;
			d1_adj_ratio[i] = 1.0;
		}
		// get 10Be T0D1 time
		d1_module[0] = event.be_x_strip[0][0] / 16;
		d1_module[1] = event.be_y_strip[0][0] / 16;
		d1_diff_time[0] = event.be_x_time[0][0] - d2_avg_time;
		d1_diff_time[1] = event.be_y_time[0][0] - d2_avg_time;
		hist_d1_diff[0]->Fill(d1_diff_time[0]);
		if (d1_module[1] != 2) {
			hist_d1_diff[1]->Fill(d1_diff_time[1]);
		} else {
			hist_d1_diff_m2[0]->Fill(d1_diff_time[1]);
		}
		// Be D1X adjacent strip
		if (event.be_x_hit[0] == 2) {
			d1_adj_diff_time[0] = event.be_x_time[0][1] - d2_avg_time;
			d1_adj_ratio[0] = event.be_x_channel[0][1]
				/ (event.be_x_channel[0][0] + event.be_x_channel[0][1]);
			d1_adj_module[0] = event.be_x_strip[0][1] / 16;
			hist_d1_adj_diff[0]->Fill(d1_adj_diff_time[0]);
		}
		// Be D1Y adjacent strip
		if (event.be_y_hit[0] == 2) {
			d1_adj_diff_time[1] = event.be_y_time[0][1] - d2_avg_time;
			d1_adj_ratio[1] = event.be_y_channel[0][1]
				/ (event.be_y_channel[0][0] + event.be_y_channel[0][1]);
			d1_adj_module[1] = event.be_y_strip[0][1] / 16;
			if (d1_adj_module[1] != 2) {
				hist_d1_adj_diff[1]->Fill(d1_adj_diff_time[1]);
			} else {
				hist_d1_adj_diff_m2[0]->Fill(d1_adj_diff_time[1]);
			}
		}
		// get 4He T0D1 time
		d1_module[2] = event.he_x_strip[0][0] / 16;
		d1_module[3] = event.he_y_strip[0][0] / 16;
		d1_diff_time[2] = event.he_x_time[0][0] - d2_avg_time;
		d1_diff_time[3] = event.he_y_time[0][0] - d2_avg_time;
		hist_d1_diff[2]->Fill(d1_diff_time[2]);
		if (d1_module[3] != 2) {
			hist_d1_diff[3]->Fill(d1_diff_time[3]);
		} else {
			hist_d1_diff_m2[1]->Fill(d1_diff_time[3]);
		}
		// He D1X adjacent strip
		if (event.he_x_hit[0] == 2) {
			d1_adj_diff_time[2] = event.he_x_time[0][1] - d2_avg_time;
			d1_adj_ratio[2] = event.he_x_channel[0][1]
				/ (event.he_x_channel[0][0] + event.he_x_channel[0][1]);
			d1_adj_module[2] = event.he_x_strip[0][1] / 16;
			hist_d1_adj_diff[2]->Fill(d1_adj_diff_time[2]);
		}
		// He D1Y adjacent strip
		if (event.he_y_hit[0] == 2) {
			d1_adj_diff_time[3] = event.he_y_time[0][1] - d2_avg_time;
			d1_adj_ratio[3] = event.he_y_channel[0][1]
				/ (event.he_y_channel[0][0] + event.he_y_channel[0][1]);
			d1_adj_module[3] = event.he_y_strip[0][1] / 16;
			if (d1_adj_module[3] != 2) {
				hist_d1_adj_diff[3]->Fill(d1_adj_diff_time[3]);
			} else {
				hist_d1_adj_diff_m2[1]->Fill(d1_adj_diff_time[3]);
			}

		}

		// check T0D1 time now
		d1_main_valid = d1_adj_valid = d1_valid = true;
		// main strip
		if (d1_diff_time[0] < 50.0 || d1_diff_time[0] > 110.0) {
			d1_main_valid = false;
		}
		if (d1_module[1] != 2) {
			if (d1_diff_time[1] < 40.0 || d1_diff_time[1] > 110.0) {
				d1_main_valid = false;
			}
		} else {
			if (d1_diff_time[1] < -160.0 || d1_diff_time[1] > -70.0) {
				d1_main_valid = false;
			}
		}
		if (d1_diff_time[2] < 50.0 || d1_diff_time[2] > 115.0) {
			d1_main_valid = false;
		}
		if (d1_module[3] != 2) {
			if (d1_diff_time[3] < 40.0 || d1_diff_time[3] > 105.0) {
				d1_main_valid = false;
			}
		} else {
			if (d1_diff_time[3] < -110.0 || d1_diff_time[3] > -10.0) {
				d1_main_valid = false;
			}
		}
		// adjacent strip
		if (
			event.be_x_hit[0] == 2
			&& (d1_adj_diff_time[0] < 60.0 || d1_adj_diff_time[0] > 200.0)
		) {
			d1_adj_valid = false;
		}
		if (event.be_y_hit[0] == 2 && d1_adj_module[1] != 2) {
			if (d1_adj_diff_time[1] < 0.0 || d1_adj_diff_time[1] > 130.0) {
				d1_adj_valid = false;
			}
		} else if (event.be_y_hit[0] == 2 && d1_adj_module[1] == 2) {
			if (d1_adj_diff_time[1] < -150.0 || d1_adj_diff_time[1] > 100.0) {
				d1_adj_valid = false;
			}
		}
		if (
			event.he_x_hit[0] == 2
			&& (d1_adj_diff_time[2] < 60.0 || d1_adj_diff_time[2] > 140.0)
		) {
			d1_adj_valid = false;
		}
		if (
			event.he_y_hit[0] == 2
			&& (d1_adj_diff_time[3] < -50.0 || d1_adj_diff_time[3] > 150.0)
		) {
			d1_adj_valid = false;
		}
		d1_valid = d1_main_valid && d1_adj_valid;


		// get SSD time difference
		for (int i = 0; i < 3; ++i) ssd_diff_time[i] = -1e5;
		for (int i = 0; i < event.layer[1]-2; ++i) {
			// jump CsI layer
			if (i == 3) break;
			ssd_diff_time[i] = event.ssd_time[i] - d2_avg_time;
		}
		// check SSD time
		ssd_valid = true;
		if (event.layer[1] >= 3) {
			if (ssd_diff_time[0] < -260.0 || ssd_diff_time[0] > -120.0) {
				ssd_valid = false;
			}
		}
		if (event.layer[1] >= 4) {
			if (ssd_diff_time[1] < -220.0 || ssd_diff_time[1] > -110.0) {
				ssd_valid = false;
			}
		}
		if (event.layer[1] >= 5) {
			if (ssd_diff_time[2] < -190.0 || ssd_diff_time[1] > -120.0) {
				ssd_valid = false;
			}
		}

		
		// get TAFD time difference
		taf_index = event.csi_index / 2;
		tafd_diff_time[0] = event.d_x_time - d2_avg_time;
		tafd_diff_time[1] = event.d_y_time - d2_avg_time;
		// check time
		tafd_valid = true;
		if (taf_index < 2) {
			if (tafd_diff_time[0] < 820.0 || tafd_diff_time[0] > 900.0) {
				tafd_valid = false;
			}
		} else {
			if (tafd_diff_time[0] < -20 || tafd_diff_time[0] > 60.0) {
				tafd_valid = false;
			}
			if (tafd_diff_time[1] < -20 || tafd_diff_time[1] > 60.0) {
				tafd_valid = false;
			}
		}


		valid = d1_valid && d2_valid && d3_valid && ssd_valid && tafd_valid;
		d1_valid_count += d1_valid ? 1 : 0;
		d2_valid_count += d2_valid ? 1 : 0;
		d3_valid_count += d3_valid ? 1 : 0;
		ssd_valid_count += ssd_valid ? 1 : 0;
		tafd_valid_count += tafd_valid ? 1 : 0;
		valid_count += valid ? 1 : 0;
		
		opt.Fill();
	}

	// save histograms
	hist_d2_diff.Write();
	hist_d2_adj_diff.Write();
	hist_d3_diff.Write();
	hist_d3_adj_diff.Write();
	for (auto &hist : hist_d1_diff) hist->Write();
	for (auto &hist : hist_d1_diff_m2) hist->Write();
	for (auto &hist : hist_d1_adj_diff) hist->Write();
	for (auto &hist : hist_d1_adj_diff_m2) hist->Write();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	// show statistics
	std::cout << "D1 valid " << d1_valid_count
		<< "  " << double(d1_valid_count) / double(total_count) << "\n"
		<< "D2 valid " << d2_valid_count
		<< "  " << double(d2_valid_count) / double(total_count) << "\n"
		<< "D3 valid " << d3_valid_count
		<< "  " << double(d3_valid_count) / double(total_count) << "\n"
		<< "SSD valid " << ssd_valid_count
		<< "  " << double(ssd_valid_count) / double(total_count) << "\n"
		<< "TAFD valid " << tafd_valid_count
		<< "  " << double(tafd_valid_count) / double(total_count) << "\n"
		<< "Valid " << valid_count
		<< "  " << double(valid_count) / double(total_count) << "\n";

	return 0;
}