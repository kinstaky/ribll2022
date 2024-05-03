#include <iostream>
#include <fstream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>

#include "include/event/threebody_info_event.h"

using namespace ribll;

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


	// output csv file name
	TString output_csv_file_name = TString::Format(
		"%s%sthreebody-info.csv", kGenerateDataPath, kInformationDir
	);
	// output file
	std::ofstream fout(output_csv_file_name.Data());


	fout << "index, run, entry, PPAC, TAF, CsI, BeLayer, HeLayer, "
		<< "tx, ty, D2Bedx, D2Bedy, D2Hedx, D2Hedy, D3Bedx, D3Bedy, D3Hedx, D3Hedy,"
		<< "D1BedE, D1HedE, D2BedE, D2HedE, D3BedE, D3HedE,"
		<< "adjacent, ratio, ch1, ch2";
	fout << "\n";
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		// position
		// calculated Be D2 X
		double calc_be_d2x =
			t0z[1] / t0z[0] * event.be_x[0]
			- (t0z[1] - t0z[0]) / t0z[0] * event.xptx;
		// calculated Be D2 Y
		double calc_be_d2y =
			t0z[1] / t0z[0] * event.be_y[0]
			- (t0z[1] - t0z[0]) / t0z[0] * event.xpty;
		// calculated Be D3 X
		double calc_be_d3x =
			t0z[2] / t0z[0] * event.be_x[0]
			- (t0z[2] - t0z[0]) / t0z[0] * event.xptx;
		// calculated Be D3 Y
		double calc_be_d3y =
			t0z[2] / t0z[0] * event.be_y[0]
			- (t0z[2] - t0z[0]) / t0z[0] * event.xpty;

		// calculated He D2 X
		double calc_he_d2x =
			t0z[1] / t0z[0] * event.he_x[0]
			- (t0z[1] - t0z[0]) / t0z[0] * event.xptx;
		// calculated He D2 Y
		double calc_he_d2y =
			t0z[1] / t0z[0] * event.he_y[0]
			- (t0z[1] - t0z[0]) / t0z[0] * event.xpty;
		// calculated He D3 X
		double calc_he_d3x =
			t0z[2] / t0z[0] * event.he_x[0]
			- (t0z[2] - t0z[0]) / t0z[0] * event.xptx;
		// calculated He D3 Y
		double calc_he_d3y =
			t0z[2] / t0z[0] * event.he_y[0]
			- (t0z[2] - t0z[0]) / t0z[0] * event.xpty;

		
		// channel
		double channel_diff[2][3];
		for (int i = 0; i < 3; ++i) {
			// Be
			channel_diff[0][i] =
				event.be_x_channel[i][0] - event.be_y_channel[i][0];
			if (event.be_x_hit[i] == 2) {
				channel_diff[0][i] += event.be_x_channel[i][1];
			}
			if (event.be_y_hit[i] == 2) {
				channel_diff[0][i] -= event.be_y_channel[i][1];
			}
			channel_diff[0][i] = fabs(channel_diff[0][i]);

			// He
			channel_diff[1][i] =
				event.he_x_channel[i][0] - event.he_y_channel[i][0];
			if (event.he_x_hit[i] == 2) {
				channel_diff[1][i] += event.he_x_channel[i][1];
			}
			if (event.he_y_hit[i] == 2) {
				channel_diff[1][i] -= event.he_y_channel[i][1];
			}
			channel_diff[1][i] = fabs(channel_diff[1][i]);
		}
		// D1X bind strips
		if (event.bind == 6) {
			channel_diff[0][0] = fabs(
				event.be_x_channel[0][0]
				- event.be_y_channel[0][0] - event.he_y_channel[0][0]
			);
			channel_diff[1][0] = channel_diff[0][0];
		}
		// D1Y bind strips
		if (event.bind == 9) {
			channel_diff[0][0] = fabs(
				event.be_x_channel[0][0] + event.he_x_channel[0][0]
				- event.be_y_channel[0][0]
			);
			channel_diff[1][0] = channel_diff[0][0];
		}
		// D2X bind strips
		if ((event.bind & 4) != 0) {
			channel_diff[0][1] = fabs(
				event.be_x_channel[1][0]
				- event.be_y_channel[1][0] - event.he_y_channel[1][0]
			);
			channel_diff[1][1] = channel_diff[0][1];
		}
		// D2Y bind strips
		if ((event.bind & 8) != 0) {
			channel_diff[0][1] = fabs(
				event.be_x_channel[1][0] + event.he_x_channel[1][0]
				- event.be_y_channel[1][0]
			);
			channel_diff[1][1] = channel_diff[0][1];
		}
		// search for adjacent strip and get ratio
		std::string min_name;
		double min_ratio = 1.0;
		double min_channel[2]{0.0, 0.0};
		for (int i = 0; i < 3; ++i) {
			if (event.be_x_hit[i] == 2 && event.layer[0] >= i) {
				double sum =
					event.be_x_channel[i][0] + event.be_x_channel[i][1];
				double ratio = event.be_x_channel[i][1] / sum;
				if (ratio < min_ratio) {
					min_ratio = ratio;
					min_name = "BeD" + std::to_string(i+1) + "X";
					min_channel[0] = event.be_x_channel[i][0];
					min_channel[1] = event.be_x_channel[i][1];
				}
			}
			if (event.be_y_hit[i] == 2 && event.layer[0] >= i) {
				double sum =
					event.be_y_channel[i][0] + event.be_y_channel[i][1];
				double ratio = event.be_y_channel[i][1] / sum;
				if (ratio < min_ratio) {
					min_ratio = ratio;
					min_name = "BeD" + std::to_string(i+1) + "Y";
					min_channel[0] = event.be_y_channel[i][0];
					min_channel[1] = event.be_y_channel[i][1];
				}
			}
			if (event.he_x_hit[i] == 2 && event.layer[1] >= i) {
				double sum =
					event.he_x_channel[i][0] + event.he_x_channel[i][1];
				double ratio = event.he_x_channel[i][1] / sum;
				if (ratio < min_ratio) {
					min_ratio = ratio;
					min_name = "HeD" + std::to_string(i+1) + "X";
					min_channel[0] = event.he_x_channel[i][0];
					min_channel[1] = event.he_x_channel[i][1];
				}
			}
			if (event.he_y_hit[i] == 2 && event.layer[1] >= i) {
				double sum =
					event.he_y_channel[i][0] + event.he_y_channel[i][1];
				double ratio = event.he_y_channel[i][1] / sum;
				if (ratio < min_ratio) {
					min_ratio = ratio;
					min_name = "HeD" + std::to_string(i+1) + "Y";
					min_channel[0] = event.he_y_channel[i][0];
					min_channel[1] = event.he_y_channel[i][1];
				}
			}
		}

		// print information
		fout << entry << ", " << event.run << ", " << event.entry
			<< ", " << event.ppac_flag << ", " << event.taf_flag
			<< ", " << event.csi_index;
		if (event.layer[0] == 1 || event.layer[0] == 2) {
			fout << ", D" << event.layer[0]+1;
		} else if (event.layer[0] < 6) {
			fout << ", S" << event.layer[0]-2;
		} else {
			fout << ", CsI";
		}
		if (event.layer[1] == 1 || event.layer[1] == 2) {
			fout << ", D" << event.layer[1]+1;
		} else if (event.layer[1] < 6) {
			fout << ", S" << event.layer[1]-2;
		} else {
			fout << ", CsI";
		}
		// print position
		fout << std::fixed << std::setprecision(1)
			<< ", " << event.xptx
			<< ", " << event.xpty
			<< ", " << fabs(calc_be_d2x - event.be_x[1])
			<< ", " << fabs(calc_be_d2y - event.be_y[1])
			<< ", " << fabs(calc_he_d2x - event.he_x[1])
			<< ", " << fabs(calc_he_d2y - event.he_y[1]);
		if (event.layer[0] < 2) fout << ", ";
		else fout << ", " << fabs(calc_be_d3x - event.be_x[2]);
		if (event.layer[0] < 2) fout << ", ";
		else fout << ", " << fabs(calc_be_d3y - event.be_y[2]);
		if (event.layer[1] < 2) fout << ", ";
		else fout << ", " << fabs(calc_he_d3x - event.he_x[2]);
		if (event.layer[1] < 2) fout << ", ";
		else fout << ", " << fabs(calc_he_d3y - event.he_y[2]);
		// print channel
		fout << std::fixed << std::setprecision(0)
			<< ", " << channel_diff[0][0] << ", " << channel_diff[1][0]
			<< ", " << channel_diff[0][1] << ", " << channel_diff[1][1];
		if (event.layer[0] >= 2) fout << ", " << channel_diff[0][2];
		else fout << ", ";
		if (event.layer[1] >= 2) fout << ", " << channel_diff[1][2];
		else fout << ", ";
		if (min_ratio == 1.0) {
			fout << ", , , , ";
		} else {
			fout << ", " << min_name
				<< std::fixed << std::setprecision(2)
				<< ", " << min_ratio
				<< std::fixed << std::setprecision(0)
				<< ", " << min_channel[0]
				<< ", " << min_channel[1];
		}

		// finish
		fout << "\n";
	}

	// close files
	fout.close();
	ipf.Close();

	return 0;
}