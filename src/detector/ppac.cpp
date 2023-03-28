#include "include/detector/ppac.h"

#include "include/event/ppac_event.h"

namespace ribll {


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

// int PPAC::ReadEvents() {
// 	events_.clear();
// 	// setup input file and tree
// 	TString file_name;
// 	file_name.Form("%s%sxppac-map-%04d.root", kGenerateDataPath, kMappingDir, run_);
// 	TFile *ipf = new TFile(file_name, "read");
// 	if (!ipf) {
// 		std::cerr << "Error: open file " << file_name << " failed.\n";
// 		return -1;
// 	}
// 	TTree *ipt = (TTree*)ipf->Get("tree");
// 	if (!ipt) {
// 		std::cerr << "Error: get tree from " << file_name << " failed.\n";
// 		ipf->Close();
// 		return -1;
// 	}

// 	// input data
// 	Long64_t timestamp;
// 	// setup input branch
// 	ipt->SetBranchAddress("timestamp", &timestamp);
// 	ipt->SetBranchAddress("index", &event_.index);
// 	ipt->SetBranchAddress("side", &event_.side);
// 	ipt->SetBranchAddress("time", &event_.time);
// 	event_.used = false;
	
// 	// show process
// 	printf("Reading ppac events   0%%");
// 	fflush(stdout);
// 	Long64_t entry100 = ipt->GetEntries() / 100;
// 	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
// 		// show process
// 		if (entry % entry100 == 0) {
// 			printf("\b\b\b\b%3lld%%", entry / entry100);
// 			fflush(stdout);
// 		}

// 		ipt->GetEntry(entry);
// 		events_.insert(std::make_pair(timestamp, event_));
// 	}
// 	// show finish
// 	printf("\b\b\b\b100%%\n");
// 	return 0;
// }

// int PPAC::Correlate() {
// 	// read ppac events
// 	if (ReadEvents()) {
// 		std::cerr << "Error: read ppac events failed.\n";
// 		return -1;
// 	}

// 	// setup input file
// 	TString input_file_name;
// 	input_file_name.Form("%s%sxt-map-%04d.root", kGenerateDataPath, kMappingDir, run_);
// 	TFile *ipf = new TFile(input_file_name, "read");
// 	if (!ipf) {
// 		std::cerr << "Error: open file " << input_file_name << " failed.\n";
// 		return -1;
// 	}
// 	TTree *ipt = (TTree*)ipf->Get("tree");
// 	if (!ipt) {
// 		std::cerr << "Error: get tree from " << input_file_name << " failed.\n";
// 		ipf->Close();
// 		return -1;
// 	}
// 	// input data and branch
// 	Long64_t timestamp;
// 	// setup input trigger branch
// 	ipt->SetBranchAddress("timestamp", &timestamp);            


// 	// setup output file
// 	TString output_file_name;
// 	output_file_name.Form("%s%sppac-corr-%04d.root", kGenerateDataPath, kCorrelationDir, run_); 
// 	TFile *opf = new TFile(output_file_name, "recreate");
// 	if (!opf) {
// 		std::cerr << "Error: create file " << output_file_name << " failed.\n";
// 		return -1;
// 	}
// 	TH1F *hist_look_window = new TH1F(
// 		"ht", "time window", 1000, -look_window, look_window
// 	);
// 	TTree *opt = new TTree("tree", "correlation tree of ppac");
// 	// output data
// 	long long ppac_timestamp;
// 	int ppac_flag;
// 	unsigned short ppac_hit;
// 	unsigned short ppac_xhit;
// 	unsigned short ppac_yhit;
// 	double ppac_x[ppac_num];
// 	double ppac_y[ppac_num];
// 	// setup output branches
// 	opt->Branch("timestamp", &ppac_timestamp, "ts/L");
// 	opt->Branch("flag", &ppac_flag, "flag/I");
// 	opt->Branch("hit", &ppac_hit, "hit/s");
// 	opt->Branch("xhit", &ppac_xhit, "xhit/s");
// 	opt->Branch("yhit", &ppac_yhit, "yhit/s");
// 	opt->Branch("x", ppac_x, TString::Format("x[%llu]/D", ppac_num).Data());
// 	opt->Branch("y", ppac_y, TString::Format("y[%llu]/D", ppac_num).Data());


// 	// show process
// 	printf("Correlating ppac events   0%%");
// 	fflush(stdout);
// 	Long64_t entry100 = ipt->GetEntries() / 100;
// 	long long correlation_num = 0;
// 	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
// 		// show process
// 		if (entry % entry100 == 0) {
// 			printf("\b\b\b\b%3lld%%", entry / entry100);
// 			fflush(stdout);
// 		}

// 		// get trigger timestamp
// 		ipt->GetEntry(entry);

// 		std::vector<Event> filling_events;
		
// 		// search events
// 		ppac_flag = 0;
// 		ppac_hit = ppac_xhit = ppac_yhit = 0;
// 		for (size_t i = 0; i < ppac_num; ++i) {
// 			ppac_x[i] = 0.0;
// 			ppac_y[i] = 0.0;
// 		}
// 		for (
// 			auto iter = events_.lower_bound(timestamp - look_window);
// 			iter != events_.upper_bound(timestamp + look_window);
// 			++iter
// 		) {
// 			hist_look_window->Fill(iter->first - timestamp);
			
// 			if (iter->first <= timestamp - search_window) continue;
// 			if (iter->first >= timestamp + search_window) continue;
// 			if (iter->second.used) continue;
// 			int flag_bit = 1 << (iter->second.index * 5 + iter->second.side);
// 			if ((flag_bit & ppac_flag) != 0) continue;

// 			iter->second.used = true;
// 			if (!ppac_flag) {
// 				// first event
// 				ppac_timestamp = iter->first;
// 			}
// 			ppac_flag |= flag_bit;
// 			++ppac_hit;
// 			filling_events.push_back(iter->second);
// 		}

// 		// sort events
// 		std::sort(filling_events.begin(), filling_events.end(), [](const Event &x, const Event &y) {
// 			return (x.index < y.index) || (x.index == y.index && x.side < y.side);
// 		});

// 		// fill to tree
// 		if (filling_events.empty()) continue;
// 		for (size_t i = 0; i < filling_events.size()-1; ++i) {
// 			if (filling_events[i].side == 0 && filling_events[i+1].side == 1) {
// 				int bit_flag
// 					= 3 << (filling_events[i].index * 5 + filling_events[i].side);
// 				if ((bit_flag & ppac_flag) == bit_flag) {
// 					++ppac_xhit;
// 					ppac_x[filling_events[i].index]
// 						= filling_events[i].time - filling_events[i+1].time;
// 				}
// 			} else if (filling_events[i].side == 2 && filling_events[i+1].side == 3) {
// 				int bit_flag
// 					= 3 << (filling_events[i].index * 5 + filling_events[i].side);
// 				if ((bit_flag & ppac_flag) == bit_flag) {
// 					++ppac_yhit;
// 					ppac_y[filling_events[i].index]
// 						= filling_events[i].time - filling_events[i+1].time;
// 				}
// 			}
// 		}
// 		if (ppac_xhit && ppac_yhit) {
// 			opt->Fill();
// 			correlation_num++;
// 		}
// 	}
// 	// show finish
// 	printf("\b\b\b\b100%%\n");
// 	hist_look_window->Write();
// 	opt->Write();
// 	opf->Close();

// 	// show correlation rate
// 	std::cout << "ppac correlation rate " << correlation_num << " / " << ipt->GetEntries()
// 		<< " " << double(correlation_num) / double(ipt->GetEntries()) << "\n";


// 	ipf->Close();

	
// 	return 0;
// }



// void SimpleFit(const double *x, double *y, double &k, double &b) {
// 	int n = 3;
// 	double sumx = 0.0;
// 	double sumy = 0.0;
// 	double sumxy = 0.0;
// 	double sumx2 = 0.0;
// 	for (int i = 0; i < n; ++i) {
// 		sumx += x[i];
// 		sumy += y[i];
// 		sumxy += x[i] * y[i];
// 		sumx2 += x[i] * x[i];
// 	}
// 	k = (sumxy - sumx*sumy/double(n)) / (sumx2 - sumx*sumx/double(n));
// 	b = (sumy - k*sumx) / double(n);
// 	// double chi2 = 0.0;
// 	// for (int i = 0; i < n; ++i) {
// 	// 	double t = y[i] - k*x[i] - b;
// 	// 	chi2 += t * t;
// 	// }
// 	// return chi2;
// }


// double PositionX(double time, int index) {
// 	double result = 0.0;
// 	if (index == 0) {
// 		result = time / 4.0 - 2.007985;
// 	} else if (index == 1) {
// 		result = time / 4.0 - 1.805155;
// 	} else {
// 		result = time / 4.0 + 0.9156165;
// 	}
// 	return result;
// }


// double PositionY(double time, int index) {
// 	double result = 0.0;
// 	if (index == 0) {
// 		result = time / 4.0 + 0.133592;
// 	} else if (index == 1) {
// 		result = time / 4.0 - 1.92654;
// 	} else {
// 		result = time / 4.0 -0.7477595;
// 	}
// 	return result;
// }


// void AddTrace(TH2D *h, double k, double b, int min, int max) {
// 	for (int i = min; i < max; ++i) {
// 		h->Fill(i, double(i)*k+b);
// 	}
// 	return;
// }

// const double ppac_xz[3] = {-695.2, -454.2, -275.2};
// const double ppac_yz[3] = {-689.2, -448.2, -269.2};

// int PPAC::Tracking() {
// 	TString input_file_name;
// 	input_file_name.Form("%s%sppac-corr-%04d.root", kGenerateDataPath, kCorrelationDir, run_); 
// 	TFile *ipf = new TFile(input_file_name, "read");
// 	TTree *ipt = (TTree*)ipf->Get("tree");
// 	if (!ipt) {
// 		std::cerr << "Error: get tree from " << input_file_name << " failed.\n";
// 		return -1;
// 	}


// 	// input data
// 	long long ppac_timestamp;
// 	int ppac_flag;
// 	unsigned short ppac_hit;
// 	unsigned short ppac_xhit;
// 	unsigned short ppac_yhit;
// 	double ppac_x[ppac_num];
// 	double ppac_y[ppac_num];
// 	// setup output branches
// 	ipt->SetBranchAddress("timestamp", &ppac_timestamp);
// 	ipt->SetBranchAddress("flag", &ppac_flag);
// 	ipt->SetBranchAddress("hit", &ppac_hit);
// 	ipt->SetBranchAddress("xhit", &ppac_xhit);
// 	ipt->SetBranchAddress("yhit", &ppac_yhit);
// 	ipt->SetBranchAddress("x", ppac_x);
// 	ipt->SetBranchAddress("y", ppac_y);


// 	// output file
// 	TString output_file_name;
// 	output_file_name.Form("%s%sppac-track-%04d.root", kGenerateDataPath, kTelescopeDir, run_);
// 	TFile *opf = new TFile(output_file_name, "recreate");
	
// 	TH2D *hxz = new TH2D("hxz", "xz tracking", 900, -800, 100, 1000, -30, 30);
// 	TH2D *hyz = new TH2D("hyz", "yz tracking", 900, -800, 100, 1000, -30, 30);
// 	TH2D *target = new TH2D("target", "target", 300, -30, 30, 300, -30, 30);

// 	// show process
// 	printf("Tracking ppac events   0%%");
// 	fflush(stdout);
// 	Long64_t entry100 = ipt->GetEntries() / 100;
// 	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
// 		// show process
// 		if (entry % entry100 == 0) {
// 			printf("\b\b\b\b%3lld%%", entry / entry100);
// 			fflush(stdout);
// 		}


// 		ipt->GetEntry(entry);
	
// 		if (ppac_hit != 15) continue;

// 		double pos_x[3];
// 		for (int i = 0; i < 3; ++i) {
// 			pos_x[i] = PositionX(ppac_x[i], i);
// 		}
// 		double xk, xb;
// 		SimpleFit(ppac_xz, pos_x, xk, xb);
// 		AddTrace(hxz, xk, xb, -800, 100);


// 		double pos_y[3];
// 		for (int i = 0; i < 3; ++i) {
// 			pos_y[i] = PositionY(ppac_y[i], i);
// 		}
// 		double yk, yb;
// 		SimpleFit(ppac_yz, pos_y, yk, yb);
// 		AddTrace(hyz, yk, yb, -800, 100);

// 		target->Fill(xb, yb);
// 	}
// 	// show finish
// 	printf("\b\b\b\b100%%\n");
// 	hxz->Write();
// 	hyz->Write();
// 	target->Write();

// 	opf->Close();

// 	ipf->Close();

// 	return 0;
// }


}