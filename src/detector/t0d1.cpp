#include "include/detector/t0d1.h"

#include <TF1.h>

#include "include/event/filter_event.h"
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

int T0d1::NormalizeSides(TChain *chain, int iteration) {
	constexpr size_t side[] = {0, 1, 0, 1, 0, 1};
	constexpr std::pair<size_t, size_t> ref_strip[] = {
		{27, 28},
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


int T0d1::Normalize(unsigned int end_run, int iteration) {
	if (iteration == 0) return Dssd::Normalize(end_run, iteration);

	constexpr NormalizeInfo normalize_info[] = {
		{0, 23, 24, 27, 31},
		{1, 27, 31, 20, 27},
		{0, 20, 27, 14, 48},
		{1, 14, 48, 5, 55},
		{0, 5, 55, 0, 64},
		{0, 0, 64, 0, 64}
	};

	// input filter file name
	TString filter_file_name = TString::Format(
		"%s%st0d1-normalize-filter-%d-%s%04d.root",
		kGenerateDataPath,
		kFilterDir,
		iteration,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// input file
	TFile *ipf = new TFile(filter_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< filter_file_name << " failed\n";
		return -1;
	}
	// input event
	FilterEvent filter_event;
	// setup input branches
	filter_event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%st0d1-normalize-fit-%s%s%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		iteration == 0 ? "" : (std::to_string(iteration)+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");


	// initialize normalized parameters
	for (size_t i = 0; i < FrontStrip(); ++i) {
		norm_params_[0][i][0] = 0.0;
		norm_params_[0][i][1] = 0.0;
		norm_params_[0][i][2] = 1.0;
		norm_params_[0][i][3] = 0.0;
	}
	for (size_t i = 0; i < BackStrip(); ++i) {
		norm_params_[1][i][0] = 0.0;
		norm_params_[1][i][1] = 0.0;
		norm_params_[1][i][2] = 1.0;
		norm_params_[1][i][3] = 0.0;
	}

	bool has_normalized[2][64];
	// initialize
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < 64; ++j) {
			has_normalized[i][j] = false;
		}
	}
	int first_side = 1 - normalize_info[0].side;
	int first_strip = normalize_info[0].ref_start;
	has_normalized[first_side][first_strip] = true;

	for (const auto &info : normalize_info) {
		// energy graph
		std::vector<TGraph> ge;
		ge.resize(Strip(info.side));

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of total entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf(
			"Filling for side %d [%d, %d) reference strips [%d, %d)   0%%",
			info.side,
			info.norm_start, info.norm_end,
			info.ref_start, info.ref_end
		);
		fflush(stdout);
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			ipt->GetEntry(entry);

			for (int i = 0; i < filter_event.num; ++i) {
				unsigned short fs = filter_event.front_strip[i];
				unsigned short bs = filter_event.back_strip[i];
				double fe = filter_event.front_energy[i];
				double be = filter_event.back_energy[i];
				if (info.side == 0) {
					// jump if not reference strips
					if (bs < info.ref_start || bs >= info.ref_end) continue;
					// jump if not normalize strips
					if (fs < info.norm_start || fs >= info.norm_end) continue;
					// jump if has normalized
					if (has_normalized[0][fs]) continue;
					// jump if referece has not normalized
					if (!has_normalized[1][bs]) continue;
					// fill to graph
					ge[fs].AddPoint(fe, NormEnergy(1, bs, be));
				} else {
					// jump if not reference strips
					if (fs < info.ref_start || fs >= info.ref_end) continue;
					// jump if not normalize strips
					if (bs < info.norm_start || bs >= info.norm_end) continue;
					// jump if has normalized
					if (has_normalized[1][bs]) continue;
					// jump if referece has not normalized
					if (!has_normalized[0][fs]) continue;
					// fill to graph
					ge[bs].AddPoint(be, NormEnergy(0, fs, fe));
				}
			}
		}
		// show finish
		printf("\b\b\b\b100%%\n");

		// fitting
		for (size_t i = info.norm_start; i < info.norm_end; ++i) {
			if (has_normalized[info.side][i]) continue;
			// only fits when over 10 points
			if (ge[i].GetN() > 5) {
				// fitting function
				TF1 energy_fit("efit", "pol1", 0, 60000);
				// set initial value
				energy_fit.SetParameter(0, 0.0);
				energy_fit.SetParameter(1, 1.0);
				// fit
				ge[i].Fit(&energy_fit, "QR+ ROB=0.8");
				// store the normalized parameters
				norm_params_[info.side][i][1] = energy_fit.GetParameter(0);
				norm_params_[info.side][i][2] = energy_fit.GetParameter(1);
				// set as normalized
				has_normalized[info.side][i] = true;
			}
			// store the graph
			ge[i].Write(TString::Format("g%c%ld", "fb"[info.side], i));
			// print normalized paramters on screen
			std::cout << i
				<< ", " << norm_params_[info.side][i][1]
				<< ", " << norm_params_[info.side][i][2] << "\n";
		}
	}

	// close files
	opf.Close();
	ipf->Close();

	// write parameters
	if (WriteNormalizeParameters(iteration)) {
		std::cerr << "Error: write normalize paramters to file failed.\n";
		return -1;
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
			for (int i = 0; i < event.front_hit; ++i) {
				time_event.front_time_flag[i] = 9;
			}
			for (int i = 0; i < event.back_hit; ++i) {
				time_event.back_time_flag[i] = 9;
			}
		} else {
			for (int i = 0; i < event.front_hit; ++i) {
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
			for (int i = 0; i < event.back_hit; ++i) {
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
	return true;
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


bool Cut11Beam(
	DssdFundamentalEvent &event,
	TH1F &hist_de,
	TH1F &hist_dt
) {
	// for convenience
	int &fhit = event.front_hit;
	int &bhit = event.back_hit;
	double *fe = event.front_energy;
	double *be = event.back_energy;
	double *ft = event.front_time;
	double *bt = event.back_time;

	for (int i = 0; i < fhit; ++i) {
		// cut single strip beam
		if (fe[i] < 31'000) continue;
		for (unsigned short j = 0; j < bhit; ++j) {
			hist_de.Fill(fe[i] - be[j]);
			hist_dt.Fill(ft[i] - bt[j]);

			if (
				fabs(ft[i] - bt[j]) < 80.0
				&& fabs(fe[i] - be[j]) < 500.0
			) {
				// // move front events
				// for (unsigned short k = i+1; k < fhit; ++k) {
				// 	fs[k-1] = fs[k];
				// 	fe[k-1] = fe[k];
				// 	ft[k-1] = ft[k];
				// 	cfd |= (cfd >> 1) & (1<<(k-1));
				// }
				// // move back events
				// for (unsigned short k = j+1; k < bhit; ++k) {
				// 	bs[k-1] = bs[k];
				// 	be[k-1] = be[k];
				// 	bt[k-1] = bt[k];
				// 	cfd |= (cfd >> 1) & (1 <<(k+8-1));
				// }
				// --fhit;
				// --bhit;
				event.Erase(0, i);
				event.Erase(1, j);
				return true;
			}
		}
	}
	return false;
}


bool Cut12Beam(
	DssdFundamentalEvent &event,
	TH1F &hist_de,
	TH1F &hist_dt,
	TH1F &hist_adt
) {
	// for convenience
	int &fhit = event.front_hit;
	int &bhit = event.back_hit;
	unsigned short *bs = event.back_strip;
	double *fe = event.front_energy;
	double *be = event.back_energy;
	double *ft = event.front_time;
	double *bt = event.back_time;

	for (int i = 0; i < fhit; ++i) {
		if (fe[i] < 31'000) continue;
		// search for adjacent strips in back side
		for (int j = 0; j < bhit-1; ++j) {
			for (int k = j+1; k < bhit; ++k) {
				if (abs(bs[j]-bs[k]) != 1) continue;
				hist_de.Fill(fe[i]-be[j]-be[k]);
				hist_dt.Fill(ft[i]-bt[j]);
				// hist_dt.Fill(ft[i]-bt[k]);
				hist_adt.Fill(bt[j]-bt[k]);
				if (fabs(fe[i]-be[j]-be[k]) > 500) continue;
				if (fabs(ft[i]-bt[j]) > 80) continue;
				if (bt[j]-bt[k] < -200 || bt[j]-bt[k] > 20) continue;
				event.Erase(0, i);
				event.Erase(1, j);
				event.Erase(1, k);
				return true;
			}
		}
	}
	return false;
}


bool Cut21Beam(
	DssdFundamentalEvent &event,
	TH1F &hist_de,
	TH1F &hist_dt,
	TH1F &hist_adt
) {
	// for convenience
	int &fhit = event.front_hit;
	int &bhit = event.back_hit;
	unsigned short *fs = event.front_strip;
	double *fe = event.front_energy;
	double *be = event.back_energy;
	double *ft = event.front_time;
	double *bt = event.back_time;

	for (int i = 0; i < bhit; ++i) {
		if (be[i] < 31'000) continue;
		// search for adjacent strips in back side
		for (int j = 0; j < fhit-1; ++j) {
			for (int k = j+1; k < fhit; ++k) {
				if (abs(fs[j]-fs[k]) != 1) continue;
				hist_de.Fill(be[i]-fe[j]-fe[k]);
				hist_dt.Fill(bt[i]-ft[j]);
				// hist_dt.Fill(bt[i]-ft[k]);
				hist_adt.Fill(ft[j]-ft[k]);
				if (fabs(be[i]-fe[j]-fe[k]) > 500) continue;
				if (fabs(bt[i]-ft[j]) > 80) continue;
				if (ft[j]-ft[k] < -200 || ft[j]-ft[k] > 100) continue;
				event.Erase(0, j);
				event.Erase(0, k);
				event.Erase(1, i);
				return true;
			}
		}
	}
	return false;
}


bool Cut22Beam(
	DssdFundamentalEvent &event,
	TH1F &hist_de,
	TH1F &hist_dt,
	TH1F &hist_adt
) {
	// for convenience
	int &fhit = event.front_hit;
	int &bhit = event.back_hit;
	unsigned short *fs = event.front_strip;
	unsigned short *bs = event.back_strip;
	double *fe = event.front_energy;
	double *be = event.back_energy;
	double *ft = event.front_time;
	double *bt = event.back_time;

	for (int i = 0; i < fhit; ++i) {
		for (int j = i+1; j < fhit; ++j) {
			if (abs(fs[i]-fs[j]) != 1) continue;
			if (fe[i]+fe[j] < 31'000) continue;
			// search for adjacent strips in back side
			for (int k = 0;  k < bhit; ++k) {
				for (int l = k+1; l < bhit; ++l) {
					if (abs(bs[k]-bs[l] != 1)) continue;
					hist_de.Fill(fe[i]+fe[j]-be[k]-be[l]);
					hist_dt.Fill(ft[i]-bt[k]);
					hist_adt.Fill(ft[i]-ft[j]);
					hist_adt.Fill(bt[k]-bt[l]);
					if (fabs(fe[i]+fe[j]-be[k]-be[l]) > 500) continue;
					if (fabs(ft[i]-bt[k]) > 80) continue;
					if (ft[i]-ft[j] < -180 || ft[i]-ft[j] > 80) continue;
					if (bt[k]-bt[l] < -180 || bt[k]-bt[l] > 80) continue;
					event.Erase(0, i);
					event.Erase(0, j);
					event.Erase(1, k);
					event.Erase(1, l);
					return true;
				}
			}
		}
	}
	return false;
}


int T0d1::CutBeamThreshold() {
	// input file name
	TString input_file_name = TString::Format(
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
	// input event
	DssdFundamentalEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%st0d1-result-%sbtc-%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// 1D histogram of energy difference of 11 beam
	TH1F de_11_beam("de11b", "#DeltaE", 1000, -2000, 2000);
	// 1D histogram of time difference of 11 beam
	TH1F dt_11_beam("dt11b", "#Deltat", 1000, -500, 500);
	// 1D histogram of energy difference of 12 beam
	TH1F de_12_beam("de12b", "#DeltaE", 1000, -2000, 2000);
	// 1D histogram of time difference of 12 beam
	TH1F dt_12_beam("dt12b", "#Deltat", 1000, -500, 500);
	// 1D histogram of time difference of 12 beam ajdacent strips
	TH1F adt_12_beam("adt12b", "#Deltat", 1000, -500, 500);
	// 1D histogram of energy difference of 21 beam
	TH1F de_21_beam("de21b", "#DeltaE", 1000, -2000, 2000);
	// 1D histogram of time difference of 21 beam
	TH1F dt_21_beam("dt21b", "#Deltat", 1000, -500, 500);
	// 1D histogram of time difference of 12 beam ajdacent strips
	TH1F adt_21_beam("adt21b", "#Deltat", 1000, -500, 500);
	// 1D histogram of energy difference of 22 beam
	TH1F de_22_beam("de22b", "#DeltaE", 1000, -2000, 2000);
	// 1D histogram of time difference of 22 beam
	TH1F dt_22_beam("dt22b", "#Deltat", 1000, -500, 500);
	// 1D histogram of time difference of 22 beam ajdacent strips
	TH1F adt_22_beam("adt22b", "#Deltat", 1000, -500, 500);
	// output tree
	TTree opt("tree", "normalize result without beam and threshold");
	// setup output branches
	event.SetupOutput(&opt);

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100;
	// show start
	printf("Cutting events of beam or under threshold   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);

		while (
			Cut11Beam(event, de_11_beam, dt_11_beam)
			|| Cut12Beam(event, de_12_beam, dt_12_beam, adt_12_beam)
			|| Cut21Beam(event, de_21_beam, dt_21_beam, adt_21_beam)
			|| Cut22Beam(event, de_22_beam, dt_22_beam, adt_22_beam)
		);
		// cut threshold
		while (
			event.front_hit > 0
			&& event.front_energy[event.front_hit-1] < 800.0
		) {
			--event.front_hit;
		}
		while (
			event.back_hit > 0
			&& event.back_energy[event.back_hit-1] < 800.0
		) {
			--event.back_hit;
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms
	de_11_beam.Write();
	dt_11_beam.Write();
	de_12_beam.Write();
	dt_12_beam.Write();
	adt_12_beam.Write();
	de_21_beam.Write();
	dt_21_beam.Write();
	adt_21_beam.Write();
	de_22_beam.Write();
	dt_22_beam.Write();
	adt_22_beam.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}


}		// namespace ribll

