#include "include/detector/ppac.h"

#include <TH1F.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TString.h>
#include <TH2F.h>
#include <Math/Vector3D.h>

#include "include/event/ppac_event.h"
#include "include/event/particle_event.h"

namespace ribll {

constexpr unsigned int kChangePpacRun = 717;

constexpr double sum_range[][19][2] = {
	{
		// 618 <= run <= 640
		// x0, y0
		{100, 120},		// x0-a0
		{100, 120},		// y0-a0
		{100, 120},		// x0-a1
		{100, 120},		// y0-a1
		{70, 90},		// x0-a2
		{70, 90},		// y0-a2
		// x1, y1
		{85, 105},		// x1-a0
		{100, 120},		// y1-a0
		{90, 110},		// x1-a1
		{100, 120},		// y1-a1
		{70, 90},		// x1-a2
		{70, 90},		// y1-a2
		// x2, y2
		{100, 120},		// x2-a0
		{110, 130},		// y2-a0
		{100, 120},		// x2-a1
		{110, 130},		// y2-a1
		{80, 100},		// x2-a2
		{80, 100}		// y2-a2
	},
	{
		// run == 641, 642, 644, 645, 648, 649, 651
		// x0, y0
		{75, 95},		// x0-a0
		{70, 90},		// y0-a0
		{85, 105},		// x0-a1
		{85, 105},		// y0-a1
		{70, 90},		// x0-a2
		{70, 90},		// y0-a2
		// x1, y1
		{90, 110},		// x1-a0
		{85, 105},		// y1-a0
		{100, 120},		// x1-a1
		{95, 115},		// y1-a1
		{80, 100},		// x1-a2
		{80, 100},		// y1-a2
		// x2, y2
		{90, 110},		// x2-a0
		{85, 105},		// y2-a0
		{100, 120},		// x2-a1
		{100, 120},		// y2-a1
		{85, 105},		// x2-a2
		{85, 105}		// y2-a2
	}
};
constexpr double fit_range[][18][2] = {
	{
		// 618 <= run <= 640
		// x0, y0
		{107, 109},		// x0-a0
		{106, 108},		// y0-a0
		{107, 110},		// x0-a1
		{106, 109},		// y0-a1
		{79, 82},		// x0-a2
		{77, 80},		// y0-a2
		// x1, y1
		{102, 105},		// x1-a0
		{109, 112},		// y1-a0
		{103, 105},		// x1-a1
		{110.5, 113},	// y1-a1
		{74, 77},		// x1-a2
		{81, 84},		// y1-a2
		// x2, y2
		{113, 116},		// x2-a0
		{120, 124},		// y2-a0
		{114, 117},		// x2-a1
		{121, 124},		// y2-a1
		{85.5, 88},		// x2-a2
		{93.5, 96}		// y2-a2
	},
	{
		// run == 641, 642, 644, 645, 648, 649, 651
		// x0, y0
		{83, 86},		// x0-a0
		{81, 84},		// y0-a0
		{95, 98},		// x0-a1
		{94, 97},		// y0-a1
		{79, 82},		// x0-a2
		{77, 80},		// y0-a2
		// x1, y1
		{94, 97},		// x1-a0
		{93, 96},		// y1-a0
		{107, 110},		// x1-a1
		{106, 109},		// y1-a1
		{90, 93},		// x1-a2
		{89, 92},		// y1-a2
		// x2, y2
		{97, 100},		// x2-a0
		{95, 98},		// y2-a0
		{110, 113},		// x2-a1
		{107, 110},		// y2-a1
		{93, 96},		// x2-a2
		{91, 94}		// y2-a2
	}
};
constexpr double run_correct[][3] = {
	// run == 643
	{10, 0, 10},
	// run == 646, 650, 652, 653, 658, 659, 671
	{0, 0, 10},
	// run == 645, 647, 654, 655, 665, 668
	{0, -10, 0},
};


Ppac::Ppac(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: Detector(run, name, tag) {
}


/// @brief convert map event to fundamental event
/// @param[in] trigger_time trigger time to match
/// @param[in] match_map map_events order by trigger time
/// @param[out] fundamental_event converted fundamental event
/// @param[inout] statistics information about statistics
///
void FillEvent(
	double trigger_time,
	const std::multimap<double, PpacMapEvent> &match_map,
	PpacFundamentalEvent &fundamental_event,
	MatchTriggerStatistics &statistics
) {
	// initialize fundamental event
	fundamental_event.flag = 0;
	fundamental_event.cfd_flag = 0;
	fundamental_event.hit = 0;
	fundamental_event.x_hit = 0;
	fundamental_event.y_hit = 0;
	for (size_t i = 0; i < ppac_num; ++i) {
		fundamental_event.x1[i] = -1e5;
		fundamental_event.x2[i] = -1e5;
		fundamental_event.y1[i] = -1e5;
		fundamental_event.y2[i] = -1e5;
		fundamental_event.anode[i] = -1e5;
	}

	auto range = match_map.equal_range(trigger_time);
	for (auto iter = range.first; iter != range.second; ++iter) {
		unsigned int index = iter->second.index * 5 + iter->second.side;
		if ((fundamental_event.flag & (1 << index)) != 0) {
			// conflict event
			++statistics.conflict_events;
			fundamental_event.flag = 0;
			break;
		} else {
			++fundamental_event.hit;
			fundamental_event.flag |= 1 << index;
			fundamental_event.cfd_flag |=
				iter->second.cfd_flag ? (1 << index) : 0;
			// fill x or y or a
			double time = iter->second.time - trigger_time;
			if (iter->second.side == 0) {
				fundamental_event.x1[iter->second.index] = time;
			} else if (iter->second.side == 1) {
				fundamental_event.x2[iter->second.index] = time;
			} else if (iter->second.side == 2) {
				fundamental_event.y1[iter->second.index] = time;
			} else if (iter->second.side == 3) {
				fundamental_event.y2[iter->second.index] = time;
			} else {
				fundamental_event.anode[iter->second.index] = time;
			}
		}
	}

	// fill x hit and y hit
	for (size_t i = 0; i < ppac_num; ++i) {
		unsigned int x_mask = 0x3 << (i * 5);
		fundamental_event.x_hit +=
			(fundamental_event.flag & x_mask) == x_mask ? 1 : 0;
		unsigned int y_mask = 0xc << (i * 5);
		fundamental_event.y_hit +=
			(fundamental_event.flag & y_mask) == y_mask ? 1 : 0;
	}

	// match event is able to track
	if (fundamental_event.x_hit >= 2 && fundamental_event.y_hit >= 2) {
		++statistics.match_events;
	}
	// flag over 0 means match at least one channel
	if (fundamental_event.flag) {
		statistics.used_events += fundamental_event.hit;
	}

	return;
}


int Ppac::MatchTrigger(
	double window_left,
	double window_right
) {
	if (name_ == "xppac") {
		return Detector::MatchTrigger<PpacMapEvent, PpacFundamentalEvent>(
			window_left,
			window_right,
			FillEvent
		);
	} else {
		// vppac
		return Detector::VmeMatchTrigger<PpacFundamentalEvent>();
	}
}


int Ppac::GetSumRange(
	PpacFundamentalEvent &fundamental,
	TTree *ipt,
	double *range
) {
	if (tag_.empty() || name_ == "vppac") {
		// 1D histogram of sum of time to find the range in large range
		std::vector<TH1F> large_sum_time;
		for (int i = 0; i < 18; ++i) {
			large_sum_time.emplace_back(
				TString::Format("hls%c%da%d", "xy"[i%2], i/6, (i/2)%3),
				TString::Format(
					"sum time of %c%d reference a%d",
					"xy"[i%2], i/6, (i/2)%3
				),
				5000, 50, 150
			);
		}

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filling sum time histograms   0%%");
		fflush(stdout);
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			ipt->GetEntry(entry);
			// fill sum time x
			for (int i = 0; i < 3; ++i) {
				unsigned int xflag = 0x3 << (5*i);
				if ((fundamental.flag & xflag) != xflag) continue;
				for (int j = 0; j < 3; ++j) {
					unsigned aflag = 0x10 << (5*j);
					if ((fundamental.flag & aflag) != aflag) continue;
					large_sum_time[(i*3+j)*2].Fill(
						fundamental.x1[i] + fundamental.x2[i]
						- fundamental.anode[j] - fundamental.anode[j]
					);
				}
			}
			// fill sum time y
			for (int i = 0; i < 3; ++i) {
				unsigned int yflag = 0xc << (5*i);
				if ((fundamental.flag & yflag) != yflag) continue;
				for (int j = 0; j < 3; ++j) {
					unsigned aflag = 0x10 << (5*j);
					if ((fundamental.flag & aflag) != aflag) continue;
					large_sum_time[(i*3+j)*2+1].Fill(
						fundamental.y1[i] + fundamental.y2[i]
						- fundamental.anode[j] - fundamental.anode[j]
					);
				}
			}

		}
		// show finish
		printf("\b\b\b\b100%%\n");

		// search for maximal bin
		// max bin
		int max_bin[18];
		for (int i = 0; i < 18; ++i) {
			max_bin[i] = 1;
			double max_content = large_sum_time[i].GetBinContent(1);
			// search for maximum value
			for (int bin = 2; bin <= large_sum_time[i].GetNbinsX(); ++bin) {
				if (large_sum_time[i].GetBinContent(bin) > max_content) {
					max_bin[i] = bin;
					max_content = large_sum_time[i].GetBinContent(bin);
				}
			}
			// std::cout << "max " << max_bin[i] << " " << max_content << std::endl;
		}

		// 1D histogram of sum of time to find the range
		std::vector<TH1F> sum_time;
		for (int i = 0; i < 18; ++i) {
			sum_time.emplace_back(
				TString::Format("hs%c%da%d", "xy"[i%2], i/6, (i/2)%3),
				TString::Format(
					"sum time of %c%d reference a%d",
					"xy"[i%2], i/6, (i/2)%3
				),
				1000,
				large_sum_time[i].GetBinCenter(max_bin[i] - 500) - 0.01,
				large_sum_time[i].GetBinCenter(max_bin[i] + 499) + 0.01
			);
			for (int bin = 1; bin <= 1000; ++bin) {
				sum_time[i].SetBinContent(
					bin,
					large_sum_time[i].GetBinContent(
						max_bin[i] - 501 + bin
					)
				);
			}
		}

		// fit and get range
		for (int i = 0; i < 18; ++i) {
			// left border of fit range
			double fit_range_left =
				large_sum_time[i].GetBinCenter(max_bin[i]) - 1.01;
			// right border of fit range
			double fit_range_right =
				large_sum_time[i].GetBinCenter(max_bin[i]) + 1.99;
			// fit function
			TF1 fit(
				TString::Format("f%c%da%d", "xy"[i%2], i/6, (i%6)/2),
				"gaus", fit_range_left, fit_range_right
			);
			sum_time[i].Fit(&fit, "QR+");
			double mean = fit.GetParameter(1);
			double sigma = fit.GetParameter(2);
			range[i*2] = mean - sigma * 3.0;
			range[i*2+1] = mean + sigma * 3.0;
		}
		// save histograms
		for (TH1F &hist : sum_time) {
			hist.Write();
		}
		for (TH1F &hist : large_sum_time) {
			hist.Write();
		}

		// write ranges to file
		// range file name
		TString range_file_name = TString::Format(
			"%s%s%s-range-%04u.txt",
			kGenerateDataPath,
			kTimeDir,
			name_.c_str(),
			run_
		);
		// output range text file
		std::ofstream fout(range_file_name.Data());
		if (!fout.good()) {
			std::cerr << "Error: Open range file "
				<< range_file_name << " failed.\n";
			return -1;
		}
		// write range
		for (size_t i = 0; i < 18; ++i) {
			fout << range[i*2] << " " << range[i*2+1] << "\n";
		}
		// close file
		fout.close();
	} else {
		// read ranges from file
		// range file name
		TString range_file_name = TString::Format(
			"%s%s%s-range-%04u.txt",
			kGenerateDataPath,
			kTimeDir,
			name_.c_str(),
			run_
		);
		// input range file
		std::ifstream fin(range_file_name.Data());
		if (!fin.good()) {
			std::cerr << "Error: Read range file "
				<< range_file_name << " failed.\n";
			return -1;
		}
		// read range
		for (size_t i = 0; i < 18; ++i) {
			fin >> range[i*2] >> range[i*2+1];
		}
		// close file
		fin.close();
	}
	return 0;
}


int Ppac::Merge(double) {
	bool vppac = name_[0] == 'v';
	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	// input event
	PpacFundamentalEvent fundamental;
	// setup input branches
	fundamental.SetupInput(ipt);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-merge-%s%04u.root",
		kGenerateDataPath,
		kMergeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output time difference histogram
	std::vector<TH1F> hist_diff_time;
	for (int i = 0; i < 3; ++i) {
		hist_diff_time.emplace_back(
			TString::Format("hdx%d", i),
			TString::Format("diff time of x%d", i),
			vppac ? 400 : 1000, -100, 100
		);
		hist_diff_time.emplace_back(
			TString::Format("hdy%d", i),
			TString::Format("diff time of y%d", i),
			vppac ? 400 : 1000, -100, 100
		);
	}
	// output tree
	TTree opt("tree", "ppac merge tree");
	// output merge event
	PpacMergeEvent merge;
	// setup output branches
	merge.SetupOutput(&opt);

	// range get from fitting
	double range[18][2];

	if (GetSumRange(fundamental, ipt, (double*)range)) {
		std::cerr << "Error: GetSumRange failed.\n";
		return -1;
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling merge events   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		// initialize merge event
		merge.xflag = 0;
		merge.yflag = 0;
		// check anode
		int valid_anode = 0;
		for (unsigned int i = 0; i < 3; ++i) {
			unsigned int aflag = 0x10 << (i * 5);
			if ((fundamental.flag & aflag) == aflag) {
				++valid_anode;
			}
		}
		// fill x
		for (int i = 0; i < 3; ++i) {
			// sum of x in range
			int x_in_range = 0;
			// x flag
			unsigned int xflag = 0x3 << (i * 5);
			// check x flag
			if ((fundamental.flag & xflag) != xflag) continue;
			// check x sum
			for (int j = 0; j < 3; ++j) {
				unsigned int aflag = 0x10 << (j*5);
				if ((fundamental.flag & aflag) != aflag) continue;
				double sum = fundamental.x1[i] + fundamental.x2[i]
					- fundamental.anode[j] - fundamental.anode[j];
				size_t index = (i*3 + j) * 2;
				if (sum > range[index][0] && sum < range[index][1]) {
					++x_in_range;
				}
			}
			// pass check
			if (
				(valid_anode == 1 && x_in_range == 1)
				|| (valid_anode >= 2 && x_in_range >= 2)
			) {
				// fill x flag
				merge.xflag |= 0x1 << i;
				// fill x value
				merge.x[i] = fundamental.x1[i] - fundamental.x2[i];
				// fill difference histogram
				hist_diff_time[i*2].Fill(merge.x[i]);
			}
		}
		// fill y
		for (int i = 0; i < 3; ++i) {
			// sum of y in range
			int y_in_range = 0;
			// y flag
			unsigned int yflag = 0xc << (i * 5);
			// check y flag
			if ((fundamental.flag & yflag) != yflag) continue;
			// check sum of y
			for (int j = 0; j < 3; ++j) {
				unsigned int aflag = 0x10 << (j*5);
				if ((fundamental.flag & aflag) != aflag) continue;
				double sum = fundamental.y1[i] + fundamental.y2[i]
					- fundamental.anode[j] - fundamental.anode[j];
				size_t index = ((i*3) + j)*2 + 1;
				if (sum > range[index][0] && sum < range[index][1]) {
					++y_in_range;
				}
			}
			// pass check
			if (
				(valid_anode == 1 && y_in_range == 1)
				|| (valid_anode >= 2 && y_in_range >= 2)
			) {
				// fill y flag
				merge.yflag |= 0x1 << i;
				// fill y value
				merge.y[i] = fundamental.y1[i] - fundamental.y2[i];
				// fill difference histogram
				hist_diff_time[i*2+1].Fill(merge.y[i]);
			}
		}
		for (int i = 0; i < 3; ++i) merge.time[i] = fundamental.anode[i];
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	opf.cd();
	// save histograms
	for (auto &hist : hist_diff_time) hist.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}


int Calibrate() {
	std::cerr << "PPAC::Calibrate is not implemented yet.\n";
	return -1;
}


double PositionXFromXIA(unsigned int run, double time, int index) {
	double result = 0.0;
	double offset[3];
	if (run <= 431) {
		offset[0] = 9.9;
		offset[1] = 9.3;
		offset[2] = -1.6;
		if (run == 422 || run == 423 || run == 427 || run == 431) {
			offset[0] -= 10.0;
		}
	} else if (run < 618) {
		offset[0] = -2.1;
		offset[1] = -2.7;
		offset[2] = -3.9;
	} else if (run <= 640) {
		offset[0] = -0.1;
		offset[1] = 1.3;
		offset[2] = -9.7;
	} else if (run <= 716) {
		offset[0] = 15.9;
		offset[1] = 1.3;
		offset[2] = 2.3;
		if (run == 643 || run == 691 || run == 696) {
			for (size_t i = 0; i < 3; ++i) offset[i] -= run_correct[0][i];
		} else if (
			run == 646 || run == 650 || run == 652 || run == 653
			|| run == 658 || run == 659 || run == 671 || run == 676
			|| run == 681 || run == 699 || run == 703 || run == 708 || run == 709
			|| run == 713 || run == 716
		) {
			for (size_t i = 0; i < 3; ++i) offset[i] -= run_correct[1][i];
		} else if (
			run == 645 || run == 647 || run == 654 || run == 655
			|| run == 665 || run == 668 || run == 679 || run == 680
			|| run == 683 || run == 685 || run == 687 || run == 701
			|| run == 702 || run == 704 || run == 705 || run == 712
			|| run == 715
		) {
			for (size_t i = 0; i < 3; ++i) offset[i] -= run_correct[2][i];
		}
	} else {
		offset[0] = 17.1;
		offset[1] = 1.3;
		offset[2] = 2.3;
		if (run ==742 || run == 743 || run == 753 || run == 755) {
			for (size_t i = 0; i < 3; ++i) offset[i] -= run_correct[1][i];
		} else if (run == 756 || run == 758 || run == 759 || run == 760 || run == 761) {
			for (size_t i = 0; i < 3; ++i) offset[i] -= run_correct[2][i];
		}
	}
	// // double correct[3] = {1.87, -0.97, -2.52};
	// double correct[3] = {0.0, -2.23, -3.40};
	// if (run <= 452) {
	// 	correct[1] = -2.26;
	// 	correct[2] = -3.41;
	// } else if (run >= 675 && run <= 716) {
	// 	correct[1] = -2.27 + 0.12;
	// 	correct[2] = -3.32 - 0.12;
	// } else if (run >= 739) {
	// 	correct[1] = -2.23 + 0.81 - 0.23;
	// 	correct[2] = -3.40 + 0.46 - 0.15;
	// }
	// result = (time - offset[index]) / 4.0 + correct[index];
	result = (time - offset[index]) / 4.0;
	// return round(result);
	return result;
}


double PositionYFromXIA(unsigned int run, double time, size_t index) {
	double result = 0.0;
	double offset[3];
	if (run <= 431) {
		offset[0] = -2.7;
		offset[1] = 5.9;
		offset[2] = 1.1;
	} else if (run < 618) {
		offset[0] = -8.7;
		offset[1] = 1.9;
		offset[2] = -2.9;
	} else if (run <= 640) {
		offset[0] = -10.7;
		offset[1] = 13.7;
		offset[2] = -2.9;
	} else if (run <= 716) {
		offset[0] = -14.7;
		offset[1] = 1.7;
		offset[2] = 3.1;
	} else {
		offset[0] = -14.9;
		offset[1] = 1.9;
		offset[2] = 3.1;
	}
	// // double correct[3] = {-0.85, 0.25, 1.42};
	// double correct[3] = {0.0, 0.84, 1.78};
	// if (run <= 452) {
	// 	correct[1] = 0.95;
	// 	correct[2] = 1.89;
	// } else if (run >= 675 && run <= 716) {
	// 	correct[1] = 0.92 + 0.06;
	// 	correct[2] = 1.92 - 0.05;
	// } else if (run >= 739) {
	// 	correct[1] = 0.84 - 0.03 - 0.09;
	// 	correct[2] = 1.78 - 0.01 - 0.03;
	// }
	// result = (time - offset[index]) / -4.0 + correct[index];
	result = (time - offset[index]) / -4.0;
	// return round(result);
	return result;
}


double SimpleFit(const double *x, double *y, double &k, double &b) {
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
	double chi2 = 0.0;
	for (int i = 0; i < n; ++i) {
		double t = y[i] - k*x[i] - b;
		chi2 += t * t;
	}
	return chi2;
}


void AddTrace(TH2F *h, double k, double b, int min, int max) {
	for (int i = min; i < max; ++i) {
		h->Fill(i, double(i)*k+b);
	}
	return;
}

constexpr double ppac_xz[4] = {-695.2, -633.7, -454.2, -275.2};
constexpr double ppac_yz[4] = {-689.2, -627.7, -448.2, -269.2};

int Ppac::Track() {
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-merge-%s%04u.root",
		kGenerateDataPath,
		kMergeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	TFile ipf(input_file_name, "read");
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input data
	PpacMergeEvent merge;
	// setup input branches
	merge.SetupInput(ipt);

	// output file
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-particle-%s%04u.root",
		kGenerateDataPath,
		kParticleDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	TFile opf(output_file_name, "recreate");
	// 2D histogram of xz tracking
	TH2F hxz("hxz", "xz tracking", 900, -800, 100, 1000, -30, 30);
	// 2D histogram of yz tracking
	TH2F hyz("hyz", "yz tracking", 900, -800, 100, 1000, -30, 30);
	// 2D histogram of reaction point distribution
	TH2F htarget("ht", "target", 300, -30, 30, 300, -30, 30);
	// 1D histogram of chi2/ndf of xz tracking
	TH1F hc2nx("hc2nx", "chi2/ndf of xz tracking", 1000, 0, 10);
	// 1D histogram of chi2/ndf of yz tracking
	TH1F hc2ny("hc2ny", "chi2/ndf of yz tracking", 1000, 0, 10);
	// 1D histogram of (x3-x1)/(x2-x1)
	TH1F hrdx("hrdx", "(x3-x1)/(x2-x1)", 1000, -10, 10);
	// 1D histogram of (y3-y1)/(y2-y1)
	TH1F hrdy("hrdy", "(y3-y1)/(y2-y1)", 1000, -10, 10);
	// output particle(beam) tree
	TTree opt("tree", "beam");
	// output particle event
	ParticleEvent particle;
	// output ppac flag
	unsigned short xflag;
	unsigned short yflag;
	// setup output branches
	particle.SetupOutput(&opt);
	opt.Branch("xflag", &xflag, "xflag/s");
	opt.Branch("yflag", &yflag, "yflag/s");


	double xz[3]{ppac_xz[0], ppac_xz[2], ppac_xz[3]};
	double yz[3]{ppac_yz[0], ppac_yz[2], ppac_yz[3]};
	if (run_ >= kChangePpacRun) {
		xz[0] = ppac_xz[1];
		yz[0] = ppac_yz[1];
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Tracking ppac events   0%%");
	fflush(stdout);
	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		// x hit
		size_t xhit = 0;
		// y hit
		size_t yhit = 0;
		// x position
		double x_position[3];
		// valid x index
		size_t x_index[3];
		// y position
		double y_position[3];
		// valid y index
		size_t y_index[3];
		// check flag and fill hit and position
		for (size_t i = 0; i < 3; ++i) {
			if ((merge.xflag & (1 << i)) != 0) {
				x_position[xhit] = PositionXFromXIA(run_, merge.x[i], i);
				x_index[xhit] = i;
				++xhit;
			}
			if ((merge.yflag & (1 << i)) != 0) {
				y_position[yhit] = PositionYFromXIA(run_, merge.y[i], i);
				y_index[yhit] = i;
				++yhit;
			}
		}
		// track
		// x tracking parameters
		double xk, xb;
		// y tracking parameters
		double yk, yb;
		// track x
		if (xhit == 2) {
			xk = (x_position[1] - x_position[0])
				/ (xz[x_index[1]] - xz[x_index[0]]);
			xb = x_position[0] - xk * xz[x_index[0]];
		} else if (xhit == 3) {
			double c2nx = SimpleFit(xz, x_position, xk, xb);
			hc2nx.Fill(c2nx);
			hrdx.Fill(
				(x_position[2]-x_position[0]) / (x_position[1]-x_position[0])
			);
		}
		// track y
		if (yhit == 2) {
			yk = (y_position[1] - y_position[0])
				/ (yz[y_index[1]] - yz[y_index[0]]);
			yb = y_position[0] - yk * yz[y_index[0]];
		} else if (yhit == 3) {
			double c2ny = SimpleFit(yz, y_position, yk, yb);
			hc2ny.Fill(c2ny);
			hrdy.Fill(
				(y_position[2]-y_position[0]) / (y_position[1]-y_position[0])
			);
		}
		// record histogram
		if (xhit >= 2) {
			AddTrace(&hxz, xk, xb, -800, 100);
		}
		if (yhit >= 2) {
			AddTrace(&hyz, yk, yb, -800, 100);
		}
		if (xhit >= 2 && yhit >= 2) htarget.Fill(xb, yb);

		// rebuild beam particle
		particle.num = 3;
		for (size_t i = 0; i < 3; ++i) {
			particle.time[i] = merge.time[i];
			particle.x[i] = PositionXFromXIA(run_, merge.x[i], i);
			particle.y[i] = PositionYFromXIA(run_, merge.y[i], i);
		}
		xflag = merge.xflag;
		yflag = merge.yflag;
		if (xhit >= 2 && yhit >= 2) {
			particle.x[3] = xb;
			particle.y[3] = yb;
			particle.z[3] = 0.0;
			ROOT::Math::XYZVector p(xk, yk, 1.0);
			p = p.Unit();
			particle.px[3] = p.X();
			particle.py[3] = p.Y();
			particle.pz[3] = p.Z();
			if (merge.time[2] > -9e4) {
				particle.time[3] = merge.time[2];
			} else if (merge.time[0] > -9e4) {
				particle.time[3] = merge.time[0];
			} else if (merge.time[1] > -9e4) {
				particle.time[3] = merge.time[1];
			} else {
				particle.time[3] = -1e5;
			}
			++particle.num;
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms
	hxz.Write();
	hyz.Write();
	htarget.Write();
	hc2nx.Write();
	hc2ny.Write();
	hrdx.Write();
	hrdy.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}


}