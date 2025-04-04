#include <TF1.h>
#include <TH2F.h>
#include <TRandom3.h>

#include "include/detector/t0d2.h"

#include "include/event/filter_event.h"

namespace ribll {

// center of t0d2, in mm
const ROOT::Math::XYZVector t0d2_center{-1.03, -0.86, 111.76};
// const ROOT::Math::XYZVector t0d2_center{-1.03, -0.86, 113.78};
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
}


int T0d2::Normalize(unsigned int, int iteration) {
	constexpr NormalizeInfo normalize_info[] = {
		{0, 12, 13, 10, 20},
		{1, 9, 20, 0, 32},
		{0, 0, 32, 0, 32}
	};

	// input file
	TFile *ipf = nullptr;
	// input tree
	TTree *ipt = nullptr;
	// input event
	DssdNormalizeEvent fundamental_event;
	FilterEvent filter_event;

	// setup input file, tree, and branches in different conditions
	if (iteration == 0) {
		// input file name
		TString input_file_name = TString::Format(
			"%s%st0d2-fundamental-%s%04u.root",
			kGenerateDataPath,
			kFundamentalDir,
			tag_.empty() ? "" : (tag_+"-").c_str(),
			run_
		);
		// input file
		ipf = new TFile(input_file_name, "read");
		// input tree
		ipt = (TTree*)ipf->Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< input_file_name << " failed\n";
			return -1;
		}
		// setup input branches
		fundamental_event.SetupInput(ipt);
	} else {
		// input filter file name
		TString filter_file_name = TString::Format(
			"%s%st0d2-normalize-filter-%d-%s%04d.root",
			kGenerateDataPath,
			kFilterDir,
			iteration,
			tag_.empty() ? "" : (tag_+"-").c_str(),
			run_
		);
		// input file
		ipf = new TFile(filter_file_name, "read");
		// input tree
		ipt = (TTree*)ipf->Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< filter_file_name << " failed\n";
			return -1;
		}
		// setup input branches
		filter_event.SetupInput(ipt);
	}

	// pixel resolution file name
	TString pixel_file_name = TString::Format(
		"%s%sshow-t0d2-pixel-%04u.root",
		kGenerateDataPath,
		kShowDir,
		run_
	);
	// pixel resolution file
	TFile pixel_file(pixel_file_name, "read");
	// resolution histogram
	TH2F *hbr = (TH2F*)pixel_file.Get("hbr");

	// output file name
	TString output_file_name = TString::Format(
		"%s%st0d2-normalize-fit-%s%s%04u.root",
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

	bool has_normalized[2][32];
	// initialize
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < 32; ++j) {
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

			if (iteration == 0) {
				// ignore multiple hit events
				if (
					fundamental_event.front_hit != 1
					|| fundamental_event.back_hit != 1
				) continue;

				unsigned short fs = fundamental_event.front_strip[0];
				unsigned short bs = fundamental_event.back_strip[0];
				double fe = fundamental_event.front_energy[0];
				double be = fundamental_event.back_energy[0];
				if (fe > 20000 || be > 20000) continue;
				if (fe < 4000 || be < 4000) continue;

				if (info.side == 0) {
					// jump if not reference strips
					if (bs < info.ref_start || bs >= info.ref_end) continue;
					// jump if not normalize strips
					if (fs < info.norm_start || fs >= info.norm_end) continue;
					// jump if has normalized
					if (has_normalized[info.side][fs]) continue;
					// jump if resolution over 8%
					if (hbr->GetBinContent(fs+1, bs+1) > 0.08) continue;
					// fill to graph
					ge[fs].AddPoint(fe, NormEnergy(1, bs, be));
				} else {
					// jump if not reference strips
					if (fs < info.ref_start || fs >= info.ref_end) continue;
					// jump if not normalize strips
					if (bs < info.norm_start || bs >= info.norm_end) continue;
					// jump if has normalized
					if (has_normalized[info.side][bs]) continue;
					// jump if resolution over 8%
					if (hbr->GetBinContent(fs+1, bs+1) > 0.08) continue;
					// fill to graph
					ge[bs].AddPoint(be, NormEnergy(0, fs, fe));
				}
			} else {
				// iteration > 0
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
	pixel_file.Close();
	ipf->Close();

	// write parameters
	if (WriteNormalizeParameters(iteration)) {
		std::cerr << "Error: write normalize paramters to file failed.\n";
		return -1;
	}

	return 0;
}


/// @brief fill hole flag
/// @param[in] merge_event DSSD merge level event 
/// @param[in] index index to fill 
/// @param[in] hbr histogram of pixel resolution 
/// @param[in] hole hole flag to fill
///  
void FillHole(
	const DssdMergeEvent &merge_event,
	int index,
	TH2F *hbr,
	bool *hole
) {
	// get strip
	int xstrip = int(merge_event.x[index]);
	int ystrip = int(merge_event.y[index]);

	// check the basic pixel
	// get resolution
	double res1 = hbr->GetBinContent(xstrip+1, ystrip+1);
	// compare
	hole[index] = res1 > 0.08;
	if (hole[index]) return;

	// check y+1 pixel
	if ((merge_event.merge_tag[index] & 1) != 0) {
		// get resolution
		double res = hbr->GetBinContent(xstrip+1, ystrip+2);
		// compare
		if (res > 0.08) {
			hole[index] = true;
			return;
		}
	}

	// check x+1 pixel
	if ((merge_event.merge_tag[index] & 2) != 0) {
		// get resolution
		double res = hbr->GetBinContent(xstrip+2, ystrip+1);
		// compare
		if (res > 0.08) {
			hole[index] = true;
			return;
		}
	}

	// check x+1,y+1 pixel
	if (merge_event.merge_tag[index] == 3) {
		// get resolution
		double res = hbr->GetBinContent(xstrip+2, ystrip+2);
		// compare
		if (res > 0.08) {
			hole[index] = true;
			return;
		}
	}

	return;
}



void FillRandomHole(
	const DssdMergeEvent &merge_event,
	int index,
	TH2F *hahp,
	TRandom3 *generator,
	bool *hole
) {
	// get strip
	int xstrip = int(merge_event.x[index]);
	int ystrip = int(merge_event.y[index]);
	if (ystrip == 17) {
		hole[index] = true;
		return;
	}

	// check the basic (x,y) pixel
	// get possiblity
	double possibility_basic = hahp->GetBinContent(xstrip+1, ystrip+1);
	// get random number
	double rand_basic = generator->Rndm();
	// compare
	hole[index] = rand_basic < possibility_basic;
	if (hole[index]) return;


	// check x,y+1 pixel, f1b2 and f2b2
	if ((merge_event.merge_tag[index] & 1) != 0) {
		// get possibility
		double possiblity = hahp->GetBinContent(xstrip+1, ystrip+2);
		// get random number
		double random = generator->Rndm();
		// compare
		if (random < possiblity) {
			hole[index] = true;
			return;
		}
	}

	// check x+1,y pixel, f2b1 and f2b2
	if ((merge_event.merge_tag[index] & 2) != 0) {
		// get possibility
		double possiblity = hahp->GetBinContent(xstrip+2, ystrip+1);
		// get random number
		double random = generator->Rndm();
		// compare
		if (random < possiblity) {
			hole[index] = true;
			return;
		}
	}

	// check x+1,y+1 pixel, f2b2
	if (merge_event.merge_tag[index] == 3) {
		// get possibility
		double possiblity = hahp->GetBinContent(xstrip+2, ystrip+2);
		// get random number
		double random = generator->Rndm();
		// compare
		if (random < possiblity) {
			hole[index] = true;
			return;
		}
	}
	
	return;
}


int T0d2::Merge(double energy_diff, int supplementary) {
	std::string tag = tag_;
	if (supplementary == 1) {
		tag = tag.empty() ? "s1" : tag+"-s1";
	}

	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-result-%s%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
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
	// input event
	DssdFundamentalEvent fundamental_event;
	// setup input branches
	fundamental_event.SetupInput(ipt);
	// for convenient
	int &fhit = fundamental_event.front_hit;
	int &bhit = fundamental_event.back_hit;

	// pixel resolution file name
	TString pixel_file_name = TString::Format(
		"%s%sshow-t0d2-pixel-%04u.root",
		kGenerateDataPath,
		kShowDir,
		run_
	);
	// pixel resolution file
	TFile pixel_file(pixel_file_name, "read");
	// resolution histogram
	TH2F *hbr = (TH2F*)pixel_file.Get("hbr");
	if (!hbr) {
		std::cerr << "Error: Get pixel resolution from "
			<< pixel_file_name << " failed.\n";
		return -1;
	}

	TH2F *hahp = nullptr;
	if (run_ == 2) {
		// average hole possiblity file name
		TString avg_hole_file_name = TString::Format(
			"%s%saverage-hole.root",
			kGenerateDataPath, kHoleDir
		);
		// average hole possiblity file
		TFile *avg_hole_file = new TFile(avg_hole_file_name, "read");
		// average hole possiblity histogram
		hahp = (TH2F*)avg_hole_file->Get("hahp");
		if (!hahp) {
			std::cerr << "Error: Get average hole possibility from "
				<< avg_hole_file_name << " failed.\n";
			return -1;
		}
	}

	TRandom3 generator(0);

	// output file name
	TString merge_file_name;
	merge_file_name.Form(
		"%s%s%s-merge-%s%04u.root",
		kGenerateDataPath,
		kMergeDir,
		name_.c_str(),
		tag.empty() ? "" : (tag+"-").c_str(),
		run_
	);
	// output file
	TFile opf(merge_file_name, "recreate");
	// output tree
	TTree opt("tree", "tree of merged events");
	// output event
	DssdMergeEvent merge_event;
	bool hole[8];
	// setup output branches
	merge_event.SetupOutput(&opt);
	opt.Branch("hole", hole, "hole[hit]/O");

	MergeStatistics statistics(run_, name_, tag);
	long long four_hit = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100;
	// show start
	printf("Writing merged events   0%%");
	fflush(stdout);
	// loop over events
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);

		int merge_num = 0;
		if (supplementary == 0) {
			merge_num = FillMergeEvent(
				fundamental_event, energy_diff, merge_event
			);
		} else if (supplementary == 1) {
			merge_num = FillMergeEventSupplementary(
				fundamental_event, energy_diff, merge_event
			);

		}

		if (fhit > 0 && bhit > 0) ++statistics.total;
		if (merge_num == 1) ++statistics.one_hit;
		else if (merge_num == 2) ++statistics.two_hit;
		else if (merge_num == 3) ++statistics.three_hit;
		else if (merge_num == 4) ++four_hit;
		for (int i = 0; i < merge_event.hit; ++i) {
			if (run_ == 2) {
				FillRandomHole(merge_event, i, hahp, &generator, hole);
			} else {
				FillHole(merge_event, i, hbr, hole);
			}
			auto position =
				CalculatePosition(merge_event.x[i], merge_event.y[i]);
			merge_event.x[i] = position.X();
			merge_event.y[i] = position.Y();
			merge_event.z[i] = position.Z();
		}
		SortMergeEvent(merge_event);
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save trees
	opt.Write();
	// close files
	opf.Close();
	pixel_file.Close();
	ipf.Close();

	statistics.merged =
		statistics.one_hit + statistics.two_hit
		+ statistics.three_hit + four_hit;
	// save and show statistics
	statistics.Write();
	statistics.Print();
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
			for (int i = 0; i < event.front_hit; ++i) {
				time_event.front_time_flag[i] = 9;
			}
			for (int i = 0; i < event.back_hit; ++i) {
				time_event.back_time_flag[i] = 9;
			}
		} else {
			for (int i = 0; i < event.front_hit; ++i) {
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
			for (int i = 0; i < event.back_hit; ++i) {
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

