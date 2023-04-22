#include "include/detector/ppac.h"

#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TH2F.h>
#include <Math/Vector3D.h>

#include "include/event/ppac_event.h"
#include "include/event/particle_event.h"

namespace ribll {

constexpr unsigned int kChangePpacRun = 717;

constexpr double sum_x_range[9][2] = {
	{100, 120},		// x0-a0
	{100, 120},		// x0-a1
	{70, 90},		// x0-a2
	{90, 110},		// x1-a0
	{90, 110},		// x1-a1
	{70, 90},		// x1-a2
	{100, 120},		// x2-a0
	{100, 120},		// x2-a1
	{80, 100}		// x3-a0
};
constexpr double sum_y_range[9][2] = {
	{100, 120},		// y0-a0
	{100, 120},		// y0-a1
	{70, 90},		// y0-a2
	{100, 120},		// y1-a0
	{100, 120},		// y1-a1
	{70, 90},		// y1-a2
	{110, 130},		// y2-a0
	{110, 130},		// y2-a1
	{80, 100}		// y2-a2
};
constexpr double fit_range[18][2] = {
	// x0 or y0
	{107, 109},
	{106, 108},
	{107, 110},
	{106, 109},
	{79, 82},
	{77, 80},
	// x1 or y1
	{102, 105},
	{109, 112},
	{103, 105},
	{110.5, 113},
	{74, 77},
	{81, 84},
	// x2 or y2
	{113, 116},
	{120, 124},
	{114, 117},
	{121, 124},
	{85.5, 88},
	{93.5, 96}
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


int Ppac::Merge(double) {
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
	// output sum time histgrams
	std::vector<TH1F> hist_sum_time;
	for (int i = 0; i < 9; ++i) {
		hist_sum_time.emplace_back(
			TString::Format("hsx%da%d", i/3, i%3),
			TString::Format("sum time of x%d reference a%d", i/3, i%3),
			1000, sum_x_range[i][0], sum_x_range[i][1]
		);
		hist_sum_time.emplace_back(
			TString::Format("hsy%da%d", i/3, i%3),
			TString::Format("sum time of y%d reference a%d", i/3, i%3),
			1000, sum_y_range[i][0], sum_y_range[i][1]
		);
	}
	// output tree
	TTree opt("tree", "ppac merge tree");
	// output merge event
	PpacMergeEvent merge;
	// setup output branches
	merge.SetupOutput(&opt);


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
				hist_sum_time[(i*3+j)*2].Fill(
					fundamental.x1[i] + fundamental.x2[i]
					-fundamental.anode[j] - fundamental.anode[j]
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
				hist_sum_time[(i*3+j)*2+1].Fill(
					fundamental.y1[i] + fundamental.y2[i]
					-fundamental.anode[j] - fundamental.anode[j]
				);
			}
		}

	}
	// show finish
	printf("\b\b\b\b100%%\n");

	double range[18][2];

	// fit sum time
	for (int i = 0; i < 18; ++i) {
		// fit function
		TF1 fit(
			TString::Format("f%c%da%d", "xy"[i%2], i/6, (i%6)/2),
			"gaus",
			fit_range[i][0],
			fit_range[i][1]
		);
		hist_sum_time[i].Fit(&fit, "QR+");
		double mean = fit.GetParameter(1);
		double sigma = fit.GetParameter(2);
		range[i][0] = mean - sigma * 3.0;
		range[i][1] = mean + sigma * 3.0;
	}


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
			int x_in_range = 0;
			unsigned xflag = 0x3 << (i * 5);
			if ((fundamental.flag & xflag) != xflag) continue;
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
			if (
				(valid_anode == 1 && x_in_range == 1)
				|| (valid_anode >= 2 && x_in_range >= 2)
			) {
				merge.xflag |= 0x1 << i;
				merge.x[i] = fundamental.x1[i] - fundamental.x2[i];
			}
		}
		// fill y
		for (int i = 0; i < 3; ++i) {
			int y_in_range = 0;
			unsigned yflag = 0xc << (i * 5);
			if ((fundamental.flag & yflag) != yflag) continue;
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
			if (
				(valid_anode == 1 && y_in_range == 1)
				|| (valid_anode >= 2 && y_in_range >= 2)
			) {
				merge.yflag |= 0x1 << i;
				merge.y[i] = fundamental.y1[i] - fundamental.y2[i];
			}
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	opf.cd();
	// save histograms
	for (auto &hist : hist_sum_time) hist.Write();
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


double PositionXFromXIA(double time, int index) {
	double result = 0.0;
	if (index == 0) {
		result = (time + 0.1) / 4.0;
	} else if (index == 1) {
		result = (time - 1.1) / 4.0;
	} else {
		result = (time + 9.7) / 4.0;
	}
	// return round(result);
	return result;
}


double PositionYFromXIA(double time, size_t index) {
	double result = 0.0;
	if (index == 0) {
		result = round((time + 6.3) / 4.0);
	} else if (index == 1) {
		result = round((time - 17.5) / 4.0);
	} else {
		result = round((time - 0.9) / 4.0);
	}
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
				x_position[xhit] = PositionXFromXIA(merge.x[i], i);
				x_index[xhit] = i;
				++xhit;
			}
			if ((merge.yflag & (1 << i)) != 0) {
				y_position[yhit] = PositionYFromXIA(merge.y[i], i);
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
				/ (yz[x_index[1]] - yz[y_index[0]]);
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
		xflag = merge.xflag;
		yflag = merge.yflag;
		if (xhit >= 2 && yhit >= 2) {
			particle.num = 1;
			particle.x[0] = xb;
			particle.y[0] = yb;
			particle.z[0] = 0.0;
			ROOT::Math::XYZVector p(xk, yk, 1.0);
			p = p.Unit();
			particle.px[0] = p.X();
			particle.py[0] = p.Y();
			particle.pz[0] = p.Z();
		} else {
			particle.num = 0;
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