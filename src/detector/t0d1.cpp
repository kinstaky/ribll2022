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


int T0d1::FilterTimeCurve() {
	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%st0d1-result-%s%04u-0.root",
		kGenerateDataPath,
		kNormalizeDir,
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
	// add friend
	ipt->AddFriend("ref=tree", ref_file_name);
	// input event
	DssdFundamentalEvent event;
	// reference time
	double ref_time;
	// setup input branches
	event.SetupInput(ipt);
	ipt->SetBranchAddress("ref.time", &ref_time);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0d1-time-curve-%s%04u.root",
		kGenerateDataPath,
		kFilterDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "time curve filter");
	// output flag
	unsigned short front_flag[8];
	unsigned short back_flag[8];
	// setup output branches
	opt.Branch("front_hit", &event.front_hit, "fhit/s");
	opt.Branch("back_hit", &event.back_hit, "bhit/s");
	opt.Branch("front_flag", front_flag, "fflag[fhit]/s");
	opt.Branch("back_flag", back_flag, "bflag[bhit]/s");

	// read fs14 cut
	std::unique_ptr<TCutG> fs14_cut = ReadCut(kTimeDir, "t0d1-fs14-curve");

	// number of filter events
	long long filter_num = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filtering   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// if (entry % 10'000 == 0) {
		// 	std::cout << "entry " << entry << "\n";
		// }
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		for (unsigned short i = 0; i < event.front_hit; ++i) {
			if (
				ref_time > -9e4
				&& event.front_strip[i] == 14
				&& (event.cfd_flag & (1 << i)) == 0
				&&  fs14_cut->IsInside(
					event.front_energy[i],
					event.front_time[i] - ref_time
				)
			) {
				front_flag[i] = 1;
				++filter_num;
			} else {
				front_flag[i] = 0;
			}
		}
		// fill tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	std::cout << "Filter events " << filter_num << "\n";
	return 0;
}


double TimeCurveFunc(double *x, double *par) {
	return par[0] + par[1]*exp(-par[2]*x[0]);
}

int T0d1::FitTimeCurve() {
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
	// add friend
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
		"%s%s%s-fit-curve-%s%04u.root",
		kGenerateDataPath,
		kTimeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// graph of dt VS energy of adjacent smaller strips
	TGraph dt_energy_small[64];
	// graph of dt VS energy of adjacent bigger strips
	TGraph dt_energy_big[64];
	// residual of fitted curve of smaller strips
	std::vector<TH1F> res_small;
	for (size_t i = 0; i < Strip(1); ++i) {
		res_small.emplace_back(
			TString::Format("rs%ld", i), "residual time",
			1000, -500, 500
		);
	}
	// residual of fitted curve of bigger strips
	std::vector<TH1F> res_big;
	for (size_t i = 0; i < Strip(1); ++i) {
		res_big.emplace_back(
			TString::Format("rb%ld", i), "residual time",
			1000, -500, 500
		);
	}

	// read normalized parameters from file
	if (ReadNormalizeParameters()) {
		std::cerr << "Error: Read normalize parameters from file failed.\n";
		return -1;
	}
	// read time normalized parameters from file
	if (ReadNormalizeTimeParameters()) {
		std::cerr <<
			"Error: Read normalize time parameters from file failed.\n";
		return -1;
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// l/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling dt-energy graph   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		if (ref_time < -9e4) continue;
		// fill back events
		if (event.back_hit == 2) {
			size_t small = 0;
			size_t big = 1;
			if (event.back_strip[small] > event.back_strip[big]) {
				small = 1;
				big = 0;
			}
			// small strip
			size_t small_strip = event.back_strip[small];
			// big strip
			size_t big_strip = event.back_strip[big];
			// small strip's energy
			double small_energy = NormEnergy(
				1, small_strip, event.back_energy[small]
			);
			// small strip's time
			double small_time = event.back_time[small]
				- norm_time_params_[1][small_strip]
				- ref_time;
			// big strip's energy
			double big_energy = NormEnergy(
				1, big_strip, event.back_energy[big]
			);
			// big strip's time
			double big_time = event.back_time[big]
				- norm_time_params_[1][big_strip]
				- ref_time;
			if (small_strip+1 == big_strip) {
				if ((event.cfd_flag & (1 << (small + 8))) == 0) {
					dt_energy_small[small_strip].AddPoint(
						small_energy, small_time
					);
				}
				if ((event.cfd_flag & (1 << (big + 8))) == 0) {
					dt_energy_big[big_strip].AddPoint(
						big_energy, big_time
					);
				}
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit
	// fit functions
	TF1 *fit_small[64];
	TF1 *fit_big[64];
	for (size_t i = 0; i < Strip(1)-1; ++i) {
		if (dt_energy_small[i].GetN() < 10) continue;
		fit_small[i] = new TF1(
			TString::Format("fs%ld", i),
			TimeCurveFunc, 0, 60000, 3
		);
		fit_small[i]->FixParameter(0, 0.0);
		fit_small[i]->SetParameter(1, 200.0);
		fit_small[i]->SetParameter(2, 1/5000.0);
		dt_energy_small[i].Fit(fit_small[i], "RQ+ rob=0.7");
	}
	for (size_t i = 1; i < Strip(1); ++i) {
		if (dt_energy_big[i].GetN() < 10) continue;
		fit_big[i] = new TF1(
			TString::Format("fb%ld", i),
			TimeCurveFunc, 0, 60000, 3
		);
		fit_big[i]->FixParameter(0, 0.0);
		fit_big[i]->SetParameter(1, 200.0);
		fit_big[i]->SetParameter(2, 1/5000.0);
		dt_energy_big[i].Fit(fit_big[i], "RQ+ rob=0.7");
	}

	// fill residual graphs
	for (size_t i = 0; i < Strip(1)-1; ++i) {
		int points = dt_energy_small[i].GetN();
		if (points < 10) continue;
		double *energy = dt_energy_small[i].GetX();
		double *dt = dt_energy_small[i].GetY();
		for (int j = 0; j < points; ++j) {
			res_small[i].Fill(dt[j] - fit_small[i]->Eval(energy[j]));
		}
	}
	for (size_t i = 1; i < Strip(1); ++i) {
		int points = dt_energy_big[i].GetN();
		if (points < 10) continue;
		double *energy = dt_energy_big[i].GetX();
		double *dt = dt_energy_big[i].GetY();
		for (int j = 0; j < points; ++j) {
			res_big[i].Fill(dt[j] - fit_big[i]->Eval(energy[j]));
		}
	}

	// save graphs
	for (size_t i = 0; i < Strip(1)-1; ++i) {
		dt_energy_small[i].Write(TString::Format("gs%ld", i));
	}
	for (size_t i = 1; i < Strip(1); ++i) {
		dt_energy_big[i].Write(TString::Format("gb%ld", i));
	}
	for (TH1F &hist : res_small) {
		hist.Write();
	}
	for (TH1F &hist :res_big) {
		hist.Write();
	}
	// close files
	opf.Close();
	ipf.Close();


	// fit
	return 0;
}


int T0d1::ReadTimeCuts() {
	// condition 1 cuts
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-msc0"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-fsc15"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-msc1"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-fsc31"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-msc2"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-fsc47"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-msc3"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-msc4"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-bsc15"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-msc5"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-bsc31"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-msc6"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-bsc47"));
	time_cuts_[1].push_back(ReadCut(kTimeDir, "t0d1-msc7"));
	// condition 2 cuts
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-mbc0"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-fbc16"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-mbc1"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-fbc32"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-mbc2"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-fbc48"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-mbc3"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-mbc4"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-bbc16"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-mbc5"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-bbc32"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-mbc6"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-bbc48"));
	time_cuts_[2].push_back(ReadCut(kTimeDir, "t0d1-mbc7"));
	return 0;
}


bool T0d1::CheckTime(
	int condition,
	size_t side,
	unsigned short strip,
	double energy,
	double time
) {
	if (condition == 0) {
		if (side == 0) {
			if (strip < 16) {
				if (fabs(time) < 20) return true;
			} else if (strip < 32) {
				if (fabs(time) < 10) return true;
			} else if (strip < 48) {
				if (time > -65 && time < 20) return true;
			} else {
				if (fabs(time) < 20) return true;
			}
		} else {
			if (fabs(time) < 10) return true;
		}
	} else if (condition == 1) {
		if (side == 0) {
			if (strip < 15) {
				if (time_cuts_[1][0]->IsInside(energy, time)) return true;
			} else if (strip == 15) {
				if (time_cuts_[1][1]->IsInside(energy, time)) return true;
			} else if (strip < 31) {
				if (time_cuts_[1][2]->IsInside(energy, time)) return true;
			} else if (strip == 31) {
				if (time_cuts_[1][3]->IsInside(energy, time)) return true;
			} else if (strip < 47) {
				if (time_cuts_[1][4]->IsInside(energy, time)) return true;
			} else if (strip == 47) {
				if (time_cuts_[1][5]->IsInside(energy, time)) return true;
			} else {
				if (time_cuts_[1][6]->IsInside(energy, time)) return true;
			}
		} else {
			if (strip < 15) {
				if (time_cuts_[1][7]->IsInside(energy, time)) return true;
			} else if (strip == 15) {
				if (time_cuts_[1][8]->IsInside(energy, time)) return true;
			} else if (strip < 31) {
				if (time_cuts_[1][9]->IsInside(energy, time)) return true;
			} else if (strip == 31) {
				if (time_cuts_[1][10]->IsInside(energy, time)) return true;
			} else if (strip < 47) {
				if (time_cuts_[1][11]->IsInside(energy, time)) return true;
			} else if (strip == 47) {
				if (time_cuts_[1][12]->IsInside(energy, time)) return true;
			} else {
				if (time_cuts_[1][13]->IsInside(energy, time)) return true;
			}
		}
	} else if (condition == 2) {
		if (side == 0) {
			if (strip < 16) {
				if (time_cuts_[2][0]->IsInside(energy, time)) return true;
			} else if (strip == 16) {
				if (time_cuts_[2][1]->IsInside(energy, time)) return true;
			} else if (strip < 32) {
				if (time_cuts_[2][2]->IsInside(energy, time)) return true;
			} else if (strip == 32) {
				if (time_cuts_[2][3]->IsInside(energy, time)) return true;
			} else if (strip < 48) {
				if (time_cuts_[2][4]->IsInside(energy, time)) return true;
			} else if (strip == 48) {
				if (time_cuts_[2][5]->IsInside(energy, time)) return true;
			} else {
				if (time_cuts_[2][6]->IsInside(energy, time)) return true;
			}
		} else {
			if (strip < 16) {
				if (time_cuts_[2][7]->IsInside(energy, time)) return true;
			} else if (strip == 16) {
				if (time_cuts_[2][8]->IsInside(energy, time)) return true;
			} else if (strip < 32) {
				if (time_cuts_[2][9]->IsInside(energy, time)) return true;
			} else if (strip == 32) {
				if (time_cuts_[2][10]->IsInside(energy, time)) return true;
			} else if (strip < 48) {
				if (time_cuts_[2][11]->IsInside(energy, time)) return true;
			} else if (strip == 48) {
				if (time_cuts_[2][12]->IsInside(energy, time)) return true;
			} else {
				if (time_cuts_[2][13]->IsInside(energy, time)) return true;
			}
		}
	}
	return false;
}


}		// namespace ribll

