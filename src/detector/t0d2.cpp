#include "include/detector/t0d2.h"

#include "include/event/t0_event.h"
#include "include/event/particle_type_event.h"

namespace ribll {

// center of t0d2, in mm
const ROOT::Math::XYZVector t0d2_center{-1.03, -0.86, 111.76};
// x range of t0d2, in mm
const std::pair<double, double> t0d2_x_range{-32.0, 32.0};
// y range of t0d2, in mm
const std::pair<double, double> t0d2_y_range{-32.0, 32.0};


T0d2::T0d2(unsigned int run, const std::string &tag)
: Dssd(run, "t0d2", tag) {

	center_ = t0d2_center;
	x_range_ = t0d2_x_range;
	y_range_ = t0d2_y_range;
}


int T0d2::NormalizeFilter(int iteration) {
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
			"%s%st0d2-filter-%s%04u-%d.root",
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

		// get cuts
		std::unique_ptr<TCutG> cut1 =
			ReadCut(kParticleIdentifyDir, "t0-d1d2-norm-1");
		std::unique_ptr<TCutG> cut2 =
			ReadCut(kParticleIdentifyDir, "t0-d2d3-norm-1");
		if (!cut1 || !cut2) {
			std::cerr << "Error: Read cut from file failed.\n";
			return -1;
		}

		// filter event, for statistics
		long long filter_num = 0;

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filtering T0D2 events in iteration %d   0%%", iteration);
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
			if (event.num != 1) {
				opt.Fill();
				continue;
			}
			if (
				(
					(event.flag[0] & 0x3) == 0x3
					&& cut1->IsInside(event.energy[0][1], event.energy[0][0])
				) || (
					(event.flag[0] & 0x6) == 0x6
					&& cut2->IsInside(event.energy[0][2], event.energy[0][1])
				)
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
			"%s%st0d2-filter-%s%04u-%d.root",
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
		printf("Filtering T0D2 events in iteration %d   0%%", iteration);
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
		std::cerr << "Error: T0D2 normalize filter in iteration "
			<< iteration << " is not implemented yet.\n";
		return -1;
	}
	return 0;
}


int T0d2::NormalizeSides(TChain *chain, int iteration) {
	constexpr size_t side[] = {0, 1, 0, 1, 0, 1};
	constexpr std::pair<size_t, size_t> ref_strip[] = {
		{19, 20},
		{14, 15},
		{15, 25},
		{10, 20},
		{5, 25},
		{10, 30}
	};
	constexpr std::pair<size_t, size_t> norm_strip[] = {
		{10, 20},
		{15, 25},
		{5, 25},
		{10, 30},
		{0, 32},
		{0, 32}
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
	return 0;
}



int T0d2::AnalyzeTime() {
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
	std::unique_ptr<TCutG> cutf = ReadCut(kTimeDir, "t0d2-f");
	std::unique_ptr<TCutG> cutf2 = ReadCut(kTimeDir, "t0d2-f2");
	std::unique_ptr<TCutG> cutb = ReadCut(kTimeDir, "t0d2-b");
	std::unique_ptr<TCutG> cutb2 = ReadCut(kTimeDir, "t0d2-b2");
	if (!cutf || !cutf2 || !cutb || !cutb2) {
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
				if (cutf->IsInside(
					event.front_energy[i], event.front_time[i]-ref_time
				)) {
					time_event.front_time_flag[i] = 0;
				} else if (cutf2->IsInside(
					event.front_energy[i], event.front_time[i]-ref_time
				)) {
					time_event.front_time_flag[i] = 0;
				} else {
					time_event.front_time_flag[i] = 8;
				}
			}
			for (unsigned short i = 0; i < event.back_hit; ++i) {
				if (cutb->IsInside(
					event.back_energy[i], event.back_time[i]-ref_time
				)) {
					time_event.back_time_flag[i] = 0;
				} else if (cutb2->IsInside(
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

