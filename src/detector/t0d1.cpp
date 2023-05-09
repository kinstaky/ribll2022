#include "include/detector/t0d1.h"

#include <TF1.h>

#include "include/event/t0_event.h"
#include "include/event/particle_type_event.h"

namespace ribll {

// center of t0d1, in mm
const ROOT::Math::XYZVector t0d1_center{0.0, 0.0, 100.0};
// x range of t0d1, in mm
const std::pair<double, double> t0d1_x_range{-32.0, 32.0};
// y range of t0d1, in mm
const std::pair<double, double> t0d1_y_range{-32.0, 32.0};

T0d1::T0d1(unsigned int run, const std::string &tag)
: Dssd(run, "t0d1", tag) {

	center_ = t0d1_center;
	x_range_ = t0d1_x_range;
	y_range_ = t0d1_y_range;
}

//-----------------------------------------------------------------------------
//								geometry
//-----------------------------------------------------------------------------

ROOT::Math::XYZVector T0d1::CalculatePosition(double fs, double bs) const {
	double x = (x_range_.second - x_range_.first) / BackStrip();
	x = x * (bs + 0.5) + x_range_.first;
	double y = (y_range_.second - y_range_.first) / FrontStrip();
	y = y * (fs + 0.5) + y_range_.first;
	ROOT::Math::XYZVector result(x, y, 0.0);
	result += center_;
	return result;
}


//-----------------------------------------------------------------------------
//								normalize
//-----------------------------------------------------------------------------

int T0d1::NormalizeFilter(int iteration) {
	if (iteration == 1) {
		// telescope file name
		TString tele_file_name;
		tele_file_name.Form(
			"%s%st0-telescope-%s%04u.root",
			kGenerateDataPath,
			kTelescopeDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_
		);
		// telescope file
		TFile tele_file(tele_file_name, "read");
		// telescope tree
		TTree *ipt = (TTree*)tele_file.Get("tree");
		// telescope event
		T0Event event;
		// setup input branches
		event.SetupInput(ipt);

		// output file name
		TString filter_file_name;
		filter_file_name.Form(
			"%s%st0d1-filter-%s%04u-%d.root",
			kGenerateDataPath,
			kNormalizeDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_,
			iteration
		);
		// filter file
		TFile filter_file(filter_file_name, "recreate");
		// filter tree
		TTree opt("tree", "filter tree");
		// filter flag
		unsigned short flag;
		// setup output branch
		opt.Branch("flag", &flag, "f/s");

		// filter event, for statistics
		long long filter_num = 0;

		// get cut
		std::unique_ptr<TCutG> cut =
			ReadCut(kParticleIdentifyDir, "t0-d1d2-norm-1");
		if (!cut) {
			std::cerr << "Error: Read cut from file failed.\n";
			return -1;
		}

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filtering T0D1 events in iteration %d   0%%", iteration);
		fflush(stdout);
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			// get event
			ipt->GetEntry(entry);
			// initialize flag
			flag = 0;
			if (
				event.num == 1
				&& (event.flag[0] & 0x3) == 0x3
				&& cut->IsInside(event.energy[0][1], event.energy[0][0])
			) {
				flag = 0x1;
				++filter_num;
			}
			opt.Fill();
		}
		// show finish
		printf("\b\b\b\b100%%\n");

		// save tree
		opt.Write();
		// close files
		filter_file.Close();
		tele_file.Close();

		// show statistics
		std::cout << "Filter rate " << filter_num << " / " << entries
			<< "  " << double(filter_num) / double(entries) << "\n";

	} else if (iteration == 2) {

		// telescope file name
		TString tele_file_name;
		tele_file_name.Form(
			"%s%st0-telescope-%s%04u.root",
			kGenerateDataPath,
			kTelescopeDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_
		);
		// telescope file
		TFile tele_file(tele_file_name, "read");
		// telescope tree
		TTree *ipt = (TTree*)tele_file.Get("tree");
		// particle type file name
		TString type_file_name;
		type_file_name.Form(
			"%s%st0-particle-type-%s%04u.root",
			kGenerateDataPath,
			kParticleIdentifyDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_
		);
		ipt->AddFriend("type=tree", type_file_name);
		// telescope event
		T0Event event;
		// particle type event
		ParticleTypeEvent type;
		// setup input branches
		event.SetupInput(ipt);
		type.SetupInput(ipt, "type.");

		// output file name
		TString filter_file_name;
		filter_file_name.Form(
			"%s%st0d1-filter-%s%04u-%d.root",
			kGenerateDataPath,
			kNormalizeDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_,
			iteration
		);
		// filter file
		TFile filter_file(filter_file_name, "recreate");
		// filter tree
		TTree opt("tree", "filter tree");
		// filter flag
		unsigned short flag;
		// setup output branch
		opt.Branch("flag", &flag, "f/s");

		// filter event, for statistics
		long long filter_num = 0;

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filtering T0D1 events in iteration %d   0%%", iteration);
		fflush(stdout);
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			// get event
			ipt->GetEntry(entry);
			// initialize flag
			flag = 0;
			if (
				event.num == 1
				&& type.charge[0] > 0
				&& type.mass[0] > 0
			) {
				flag = 0x1;
				++filter_num;
			}
			opt.Fill();
		}
		// show finish
		printf("\b\b\b\b100%%\n");

		// save tree
		opt.Write();
		// close files
		filter_file.Close();
		tele_file.Close();

		// show statistics
		std::cout << "Filter rate " << filter_num << " / " << entries
			<< "  " << double(filter_num) / double(entries) << "\n";

	} else {
		std::cerr << "Error: T0D1 normalize filter in iteration "
			<< iteration << " is not implemented yet.\n";
		return -1;
	}
	return 0;
}


int T0d1::NormalizeSides(TChain *chain, int iteration) {
	constexpr size_t side[] = {0, 1, 0, 1, 0, 1};
	constexpr std::pair<size_t, size_t> ref_strip[] = {
		{26, 27},
		{35, 36},
		{23, 30},
		{32, 39},
		{16, 48},
		{16, 48}
	};
	constexpr std::pair<size_t, size_t> norm_strip[] = {
		{32, 39},
		{23, 30},
		{16, 48},
		{16, 48},
		{0, 64},
		{0, 64}
	};

	for (size_t i = 0; i < 6; ++i) {
		if (StripsNormalize(
			chain,
			side[i],
			ref_strip[i].first, ref_strip[i].second,
			norm_strip[i].first, norm_strip[i].second,
			iteration
		)) {
			std::cerr << "Error: Normalize strips [" << norm_strip[i].first << ", "
				<< norm_strip[i].second << ") in side " << side[i]
				<< " reference strips [" << ref_strip[i].first << ", "
				<< ref_strip[i].second << ") failed.\n";
			return -1;
		}
	}
	return 0;
}


int T0d1::AnalyzeTime() {
	// input fundamental file name
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
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// reference time file name
	TString ref_file_name;
	ref_file_name.Form(
		"%s%sreftime-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// add reference time tree
	ipt->AddFriend("ref=tree", ref_file_name);
	// input event
	DssdFundamentalEvent event;
	// input reference time
	double ref_time;
	// setup input branches
	event.SetupInput(ipt);
	ipt->SetBranchAddress("ref.time", &ref_time);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-time-%s%04u.root",
		kGenerateDataPath,
		kTimeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "time");
	// output time flag event
	DssdTimeEvent time_event;
	// setup output branches
	time_event.SetupOutput(&opt);

	// read cuts
	std::unique_ptr<TCutG> cutf1 = ReadCut(kTimeDir, "t0d1-f1");
	std::unique_ptr<TCutG> cutf2 = ReadCut(kTimeDir, "t0d1-f2");
	std::unique_ptr<TCutG> cutb = ReadCut(kTimeDir, "t0d1-b");
	if (!cutf1 || !cutf2 || !cutb) {
		std::cerr << "Error: Read cut failed.\n";
		return -1;
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Analyzing time   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		time_event.front_hit = event.front_hit;
		time_event.back_hit = event.back_hit;
		if (ref_time < -9e4) {
			for (unsigned short i = 0; i < event.front_hit; ++i) {
				time_event.front_time_flag[i] = 9;
			}
			for (unsigned short i = 0; i < event.back_hit; ++i) {
				time_event.back_time_flag[i] = 9;
			}
		} else {
			for (unsigned short i = 0; i < event.front_hit; ++i) {
				if (event.front_strip[i] < 32 || event.front_strip[i] >= 48) {
					if (cutf1->IsInside(
						event.front_energy[i], event.front_time[i]-ref_time
					)) {
						time_event.front_time_flag[i] = 0;
					} else {
						time_event.front_time_flag[i] = 8;
					}
				} else {
					if (cutf2->IsInside(
						event.front_energy[i], event.front_time[i]-ref_time
					)) {
						time_event.front_time_flag[i] = 0;
					} else {
						time_event.front_time_flag[i] = 8;
					}
				}
			}
			for (unsigned short i = 0; i < event.back_hit; ++i) {
				if (cutb->IsInside(
					event.back_energy[i], event.back_time[i]-ref_time
				)) {
					time_event.back_time_flag[i] = 0;
				} else {
					time_event.back_time_flag[i] = 8;
				}
			}
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}

}		// namespace ribll

