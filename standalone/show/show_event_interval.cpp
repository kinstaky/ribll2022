#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TH1F.h>

#include "include/defs.h"
#include "include/statistics/xt_period_statistics.h"

using namespace ribll;


int DssdInterval(int run, const std::string &detector) {
	int strips = 32;
	if (detector == "t0d1") strips = 64;

	// file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-map-nc-%04u.root",
		kGenerateDataPath, kMappingDir, detector.c_str(), run
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
	unsigned short side;
	unsigned short strip;
	double time;
	// setup branches
	ipt->SetBranchAddress("side", &side);
	ipt->SetBranchAddress("strip", &strip);
	ipt->SetBranchAddress("time", &time);

	TString output_file_name;
	output_file_name.Form(
		"%s%sshow-interval-%s-%04u.root",
		kGenerateDataPath, kShowDir, detector.c_str(), run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histogram
	std::vector<TH1F> hist_time;
	for (int i = 0; i < strips; ++i) {
		hist_time.emplace_back(
			TString::Format("hf%d", i), "time interval of front strip",
			1000, 0, 10'000'000
		);
	}
	for (int i = 0; i < strips; ++i) {
		hist_time.emplace_back(
			TString::Format("hb%d", i), "time interval of back strip",
			1000, 0, 10'000'000
		);
	}

	// last entry's trigger time
	double last_time[2][64];
	double min_time[2][64];
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < strips; ++j) {
			last_time[i][j] = -1e5;
			min_time[i][j] = 1e20;
		}
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100;
	// show start
	std::cout << "Filling histogram   0%";
	fflush(stdout);
	// loop to fill histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// showing process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		ipt->GetEntry(entry);
		// check if it's the first event
		if (last_time[side][strip] < -9e4) {
			last_time[side][strip] = time;
			continue;
		}
		// interval
		double interval = time - last_time[side][strip];
		// fill period
		hist_time[side*strips+strip].Fill(interval);
		// check if shorter than min time
		min_time[side][strip] = interval < min_time[side][strip] ?
			interval : min_time[side][strip];
		// update trigger time of last entry
		last_time[side][strip] = time;
	}
	// show end
	std::cout << "\b\b\b\b100%\n";

	// fit
	// for (int i = 0; i < 128; ++i) {
	// 	TF1 *f1 = new TF1(TString::Format("fit%d", i), "[0]*exp([1]*x)", 0, 20'000, 90'000);
	// 	f1->SetParameter(0, hist_time[i].GetMaximum());
	// 	f1->SetParameter(1, 1e-6);
	// 	hist_time[i].Fit(f1, "RQ+");
	// 	std::cout << "fb"[i>=64] << i%64 << " " << f1->GetParameter(1)*1e9 << "\n";
	// }

	for (int i = 0; i < strips; ++i) {
		std::cout << "f" << i << " " << min_time[0][i] << "\n";
	}
	for (int i = 0; i < strips; ++i) {
		std::cout << "b" << i << " " << min_time[1][i] << "\n";
	}

	// write histogram to output file
	for (auto &hist : hist_time) {
		hist.Write();
	}
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}


int TriggerInterval(int run, const std::string &detector) {
	// file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-map-%04u.root",
		kGenerateDataPath, kMappingDir, detector.c_str(), run
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
	// input trigger time
	double trigger_time;
	// setup branches
	ipt->SetBranchAddress("time", &trigger_time);

	TString output_file_name;
	output_file_name.Form(
		"%s%sshow-interval-%s-%04u.root",
		kGenerateDataPath, kShowDir, detector.c_str(), run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histogram
	TH1F hist_time("ht", "time interval of trigger", 1000, 0, 100'000);

	// last entry's trigger time
	double last_time;
	// get the first entry and fill to last time variable
	ipt->GetEntry(0);
	last_time = trigger_time;

	// longer than one century
	double min_time = 1e20;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100;
	// show start
	std::cout << "Filling histogram   0%";
	fflush(stdout);
	// loop to fill histogram
	for (long long entry = 1; entry < entries; ++entry) {
		// showing process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		ipt->GetEntry(entry);
		// interval
		double interval = trigger_time - last_time;
		// fill interval
		hist_time.Fill(interval);
		// check if shorter than min time
		min_time = interval < min_time ? interval : min_time;
		// update trigger time of last entry
		last_time = trigger_time;
	}
	// show end
	std::cout << "\b\b\b\b100%\n";

	// fit
	TF1 *f1 = new TF1("f1", "[0]*exp(-[1] * x)", 20'000, 90'000);
	f1->SetParameter(0, hist_time.GetMaximum());
	f1->SetParameter(1, 1e-5);
	hist_time.Fit(f1, "RQ+");
	std::cout << "Rate " << f1->GetParameter(1)*1e9 << "\n";

	// write histogram to output file
	hist_time.Write();
	// close files
	opf.Close();
	ipf.Close();

	XiaTriggerPeriodStatistics statistics(run, min_time);
	statistics.Write();
	statistics.Print();

	return 0;
}


int TofInterval(int run) {
	const int bins = 1'000;
	// file name
	TString input_file_name;
	input_file_name.Form(
		"%s%stof-map-%04u.root",
		kGenerateDataPath, kMappingDir, run
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
	// input ToF time
	double time;
	// input ToF index
	unsigned short index;
	// setup branches
	ipt->SetBranchAddress("index", &index);
	ipt->SetBranchAddress("time", &time);

	TString output_file_name;
	output_file_name.Form(
		"%s%sshow-interval-tof-%04u.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histogram
	// original time interval histogram
	TH1F hist_time[2];
	hist_time[0] = TH1F("ht0", "time interval of ToF1", bins, 0, 100'000);
	hist_time[1] = TH1F("ht1", "time interval of ToF2", bins, 0, 100'000);
	// time interval histogram without expo distribution
	TH1F hist_time_peak[2];
	hist_time_peak[0] =
		TH1F("hp0", "time interval of ToF1 peaks", bins, 0, 100'000);
	hist_time_peak[1] =
		TH1F("hp1", "time interval of ToF2 peaks", bins, 0, 100'000);

	// last entry's trigger time
	double last_time[2];
	last_time[0] = last_time[1] = -1e5;

	// longer than one century
	double min_time[2];
	min_time[0] = min_time[1] = 1e20;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100;
	// show start
	std::cout << "Filling histogram   0%";
	fflush(stdout);
	// loop to fill histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// showing process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		ipt->GetEntry(entry);
		// check if it's the first event
		if (last_time[index] < -9e4) {
			last_time[index] = time;
			continue;
		}
		// interval
		double interval = time - last_time[index];
		// fill interval
		hist_time[index].Fill(interval);
		// check if shorter than min time
		min_time[index] = interval < min_time[index] ?
			interval : min_time[index];
		// update last time
		last_time[index] = time;
	}
	// show end
	std::cout << "\b\b\b\b100%\n";

	TF1 *fit_func[2];
	for (int i = 0; i < 2; ++i) {
		fit_func[i] = new TF1(TString::Format("f%d", i), "expo", 20'000, 90'000);
		fit_func[i]->SetParameter(0, hist_time[i].GetMaximum());
		fit_func[i]->SetParameter(1, 1e-5);
		hist_time[i].Fit(fit_func[i], "RQ+");
		std::cout << "Rate " << fit_func[i]->GetParameter(1)*-1e9
			<< ", min time " << min_time[i] << "\n";
	}

	// sub the expo distribution
	for (int i = 0; i < 2; ++i) {
		for (int bin = 1; bin <= bins; ++bin) {
			double content = hist_time[i].GetBinContent(bin);
			content -= fit_func[i]->Eval(hist_time[i].GetBinCenter(bin));
			if (bin < bins/20) hist_time_peak[i].SetBinContent(bin, 0);
			else hist_time_peak[i].SetBinContent(bin, content);
		}
	}

	// write histogram to output file
	hist_time[0].Write();
	hist_time[1].Write();
	hist_time_peak[0].Write();
	hist_time_peak[1].Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}


int MissToFInterval(int run) {
	const int bins = 1'000;
	// file name
	TString input_file_name;
	input_file_name.Form(
		"%s%stof-fundamental-%04d.root",
		kGenerateDataPath, kFundamentalDir, run
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
	// add xia trigger friend
	ipt->AddFriend("xt=tree", TString::Format(
		"%s%sxt-map-%04d.root",
		kGenerateDataPath, kMappingDir, run
	));
	// input ToF time
	double tof_time[2];
	// input XIA trigger time
	double xt_time;
	// setup branches
	ipt->SetBranchAddress("time", tof_time);
	ipt->SetBranchAddress("xt.time", &xt_time);

	TString output_file_name;
	output_file_name.Form(
		"%s%sshow-interval-tof-miss-%04u.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histogram
	// original time interval histogram
	TH1F hist_time[2];
	hist_time[0] =
		TH1F("ht0", "time interval of missing ToF1", bins, 0, 100'000);
	hist_time[1] =
		TH1F("ht1", "time interval of missing ToF2", bins, 0, 100'000);

	ipt->GetEntry(0);
	// last entry's trigger time
	double last_time = xt_time;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100;
	// show start
	std::cout << "Filling histogram   0%";
	fflush(stdout);
	// loop to fill histogram
	for (long long entry = 1; entry < entries; ++entry) {
		// showing process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		ipt->GetEntry(entry);

		for (int i = 0; i < 2; ++i) {
			if (tof_time[i] > -9e4) continue;
			// interval
			double interval = xt_time - last_time;
			// fill interval
			hist_time[i].Fill(interval);
		}
		// update last time
		last_time = xt_time;
	}
	// show end
	std::cout << "\b\b\b\b100%\n";

	// TF1 *fit_func[2];
	// for (int i = 0; i < 2; ++i) {
	// 	fit_func[i] = new TF1(TString::Format("f%d", i), "expo", 20'000, 90'000);
	// 	fit_func[i]->SetParameter(0, hist_time[i].GetMaximum());
	// 	fit_func[i]->SetParameter(1, 1e-5);
	// 	hist_time[i].Fit(fit_func[i], "RQ+");
	// 	std::cout << "Rate " << fit_func[i]->GetParameter(1)*-1e9
	// 		<< ", min time " << min_time[i] << "\n";
	// }

	// write histogram to output file
	hist_time[0].Write();
	hist_time[1].Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}


int PpacInterval(int run) {
	// file name
	TString input_file_name;
	input_file_name.Form(
		"%s%sxppac-map-%04u.root",
		kGenerateDataPath, kMappingDir, run
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
	unsigned short index;
	unsigned short side;
	double time;
	// setup branches
	ipt->SetBranchAddress("index", &index);
	ipt->SetBranchAddress("side", &side);
	ipt->SetBranchAddress("time", &time);

	TString output_file_name;
	output_file_name.Form(
		"%s%sshow-interval-xppac-%04u.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histogram
	std::vector<TH1F> hist_time;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 5; ++j) {
			hist_time.emplace_back(
				TString::Format("h%ds%d", i, j),
				TString::Format("time interval of xppac %d side %d", i, j),
				1000, 0, 100'000
			);
		}
	}

	// last entry's trigger time
	double last_time[3][5];
	double min_time[3][5];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 5; ++j) {
			last_time[i][j] = -1e5;
			min_time[i][j] = 1e20;
		}
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100;
	// show start
	std::cout << "Filling histogram   0%";
	fflush(stdout);
	// loop to fill histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// showing process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		ipt->GetEntry(entry);
		// check if it's the first event
		if (last_time[index][side] < -9e4) {
			last_time[index][side] = time;
			continue;
		}
		// interval
		double interval = time - last_time[index][side];
		// fill period
		hist_time[index*5+side].Fill(interval);
		// check if shorter than min time
		min_time[index][side] = interval < min_time[index][side] ?
			interval : min_time[index][side];
		// update trigger time of last entry
		last_time[index][side] = time;
	}
	// show end
	std::cout << "\b\b\b\b100%\n";

	// fit
	// for (int i = 0; i < 15; ++i) {
	// 	TF1 *f1 = new TF1(TString::Format("fit%d", i), "[0]*exp([1]*x)", 0, 20'000, 90'000);
	// 	f1->SetParameter(0, hist_time[i].GetMaximum());
	// 	f1->SetParameter(1, 1e-6);
	// 	hist_time[i].Fit(f1, "RQ+");
	// 	std::cout << "fb"[i>=64] << i%64 << " " << f1->GetParameter(1)*1e9 << "\n";
	// }

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 5; ++j) {
			std::cout << "xppac" << i << " side " << j
				<< " " << min_time[0][i] << "\n";
		}
	}

	// write histogram to output file
	for (auto &hist : hist_time) {
		hist.Write();
	}
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}


int main(int argc, char **argv) {

	// check arguments
	if (argc != 3 && argc != 4) {
		std::cout << "Usage: " << argv[0] << " [options] run \n"
			<< "  run                 run number.\n"
			<< "Options:\n"
			<< "  -m                  missing events.\n";
		return -1;
	}

	if (argc == 4) {
		// check options
		if (argv[1][0] != '-' || argv[1][1] != 'm') {
			std::cerr << "Error: Invalid option " << argv[1] << "\n";
			return -1;
		}
		// run number
		int run = atoi(argv[2]);
		// detector name
		std::string detector(argv[3]);
		// check detector name
		if (detector != "tof") {
			std::cerr << "Error: Invalid detector " << detector << "\n";
			return -1;
		}
		if (MissToFInterval(run)) return -1;
		return 0;
	}

	// run number
	int run = atoi(argv[1]);
	// detector name
	std::string detector(argv[2]);

	if (detector == "xt" || detector == "vt") {
		if (TriggerInterval(run, detector)) return -1;
	} else if (detector =="tof") {
		if (TofInterval(run)) return -1;
	} else if (detector == "xppac") {
		if (PpacInterval(run)) return -1;
	} else {
		if (DssdInterval(run, detector)) return -1;
	}

	return 0;
}
