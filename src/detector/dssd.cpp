#include "include/detector/dssd.h"

#include <array>

#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH2F.h>

#include "include/event/dssd_event.h"
#include "include/statistics/normalize_statistics.h"


namespace ribll {

Dssd::Dssd(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: Detector(run, name, tag) {
}

//-----------------------------------------------------------------------------
//								geometry
//-----------------------------------------------------------------------------

ROOT::Math::XYZVector Dssd::CalculatePosition(double fs, double bs) const {
	double x = (x_range_.second - x_range_.first) / FrontStrip();
	x = x * (fs + 0.5) + x_range_.first;
	double y = (y_range_.second - y_range_.first) / BackStrip();
	y = y * (bs + 0.5) + y_range_.first;
	ROOT::Math::XYZVector result(x, y, 0.0);
	result += center_;
	return result;
}


//-----------------------------------------------------------------------------
//								match trigger
//-----------------------------------------------------------------------------

/// @brief convert map event to fundamental event
/// @param[in] trigger_time trigger time to match
/// @param[in] match_map map_events order by trigger time
/// @param[out] fundamental_event converted fundamental event
/// @param[out] statistics information about statistics
/// @returns 0 if filling event is valid, -1 otherwise
///
int FillEvent(
	double trigger_time,
	const std::multimap<double, DssdMapEvent> &match_map,
	DssdFundamentalEvent &fundamental_event,
	MatchTriggerStatistics &statistics
) {
	// initialize output fundamental event
	fundamental_event.front_hit = 0;
	fundamental_event.back_hit = 0;
	fundamental_event.cfd_flag = 0;

	std::vector<DssdMapEvent> front_events;
	std::vector<DssdMapEvent> back_events;

	// check match events number
	auto range = match_map.equal_range(trigger_time);
	for (auto iter = range.first; iter != range.second; ++iter) {
		if (iter->second.side == 0) {
			front_events.push_back(iter->second);
		} else {
			back_events.push_back(iter->second);
		}
	}

	// sort events by energy
	std::sort(
		front_events.begin(),
		front_events.end(),
		[](const DssdMapEvent &x, const DssdMapEvent &y) {
			return x.energy > y.energy;
		}
	);
	std::sort(
		back_events.begin(),
		back_events.end(),
		[](const DssdMapEvent &x, const DssdMapEvent &y) {
			return x.energy > y.energy;
		}
	);

	// record events
	if (front_events.size() > 8 || back_events.size() > 8) {
		// Front or back events is more than 8, which is out of consideration.
		fundamental_event.front_hit = 0;
		fundamental_event.back_hit = 0;
		++statistics.oversize_events;
	} else {
		fundamental_event.front_hit = front_events.size();
		fundamental_event.back_hit = back_events.size();
	}

	// check multiple hit on the same strip
	bool strip_multiple_hit = false;
	for (size_t i = 1; i < front_events.size(); ++i) {
		if (front_events[i].strip == front_events[i-1].strip) {
			strip_multiple_hit = true;
			fundamental_event.front_hit = 0;
		}
	}
	for (size_t i = 1; i < back_events.size(); ++i) {
		if (back_events[i].strip == back_events[i-1].strip) {
			strip_multiple_hit = true;
			fundamental_event.back_hit = 0;
		}
	}
	statistics.conflict_events += strip_multiple_hit ? 1 : 0;


	for (unsigned short i = 0; i < fundamental_event.front_hit; ++i) {
		fundamental_event.front_strip[i] = front_events[i].strip;
		fundamental_event.front_time[i] = front_events[i].time - trigger_time;
		fundamental_event.front_energy[i] = front_events[i].energy;
		fundamental_event.cfd_flag |=
			front_events[i].cfd_flag ? (1 << i) : 0;
	}
	for (unsigned short i = 0; i < fundamental_event.back_hit; ++i) {
		fundamental_event.back_strip[i] = back_events[i].strip;
		fundamental_event.back_time[i] = back_events[i].time - trigger_time;
		fundamental_event.back_energy[i] = back_events[i].energy;
		fundamental_event.cfd_flag |=
			back_events[i].cfd_flag ? (1 << (i + 8)) : 0;
	}

	if (
		fundamental_event.front_hit > 0
		&& fundamental_event.back_hit > 0
	) {
		++statistics.match_events;
		statistics.used_events +=
			fundamental_event.front_hit + fundamental_event.back_hit;
		return 0;
	}

	return -1;
}


void FillEventInMatch(
	double trigger_time,
	const std::multimap<double, DssdMapEvent> &match_map,
	DssdFundamentalEvent &fundamental_event,
	MatchTriggerStatistics &statistics
) {
	FillEvent(trigger_time, match_map, fundamental_event, statistics);
}


int FillEventInExtract(
	double trigger_time,
	const std::multimap<double, DssdMapEvent> &match_map,
	DssdFundamentalEvent &fundamental_event,
	MatchTriggerStatistics &statistics
) {
	return
		FillEvent(trigger_time, match_map, fundamental_event, statistics);
}


int Dssd::MatchTrigger(
	double window_left,
	double window_right
) {
	return Detector::MatchTrigger<DssdMapEvent, DssdFundamentalEvent>(
		window_left,
		window_right,
		FillEventInMatch
	);
}


int Dssd::ExtractTrigger(
	double window_left,
	double window_right
) {
	return Detector::ExtractTrigger<DssdMapEvent, DssdFundamentalEvent>(
		name_,
		window_left,
		window_right,
		FillEventInExtract
	);
}



//-----------------------------------------------------------------------------
//									normalize
//-----------------------------------------------------------------------------

int Dssd::ReadNormalizeParameters(int iteration) {
	// setup file name
	TString file_name;
	file_name.Form(
		"%s%s%s%s.txt",
		kGenerateDataPath, kNormalizeDir, name_.c_str(),
		iteration == -1 ? "" : ("-"+std::to_string(iteration)).c_str()
	);
	// open file
	std::ifstream fin(file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: open normalize file "
			<< file_name << " failed.\n";
		return -1;
	}
	// setup strip number
	size_t strip_num;
	// read front strip number and normalize parameters
	fin >> strip_num;
	// read normalized paramters for front strips
	for (size_t i = 0; i < strip_num; ++i) {
		fin >> norm_params_[0][i][0] >> norm_params_[0][i][1];
	}
	//  read back strip number
	fin >> strip_num;
	// read normalized parameters for back strips
	for (size_t i = 0; i < strip_num; ++i) {
		fin >> norm_params_[1][i][0] >> norm_params_[1][i][1];
	}
	// close file
	fin.close();

	return 0;
}


int Dssd::WriteNormalizeParameters(int iteration) {
	// setup file name
	TString file_name;
	file_name.Form(
		"%s%s%s%s.txt",
		kGenerateDataPath, kNormalizeDir, name_.c_str(),
		iteration == -1 ? "" : ("-"+std::to_string(iteration)).c_str()
	);
	// open output file
	std::ofstream fout(file_name.Data());
	if (!fout.good()) {
		std::cerr << "Error: open normalize parameters file "
			<< file_name << " failed.\n";
		return -1;
	}
	// write front strip number
	fout << FrontStrip() << "\n";
	// write normalized paramters for front strips
	for (size_t i = 0; i < FrontStrip(); ++i) {
		fout << norm_params_[0][i][0]<< " "
			<< norm_params_[0][i][1] << "\n";
	}
	// write back strip number
	fout << BackStrip() << "\n";
	// write normalized parameters for back strips
	for (size_t i = 0; i < BackStrip(); ++i) {
		fout << norm_params_[1][i][0] << " "
			<< norm_params_[1][i][1] << "\n";
	}
	// close file
	fout.close();

	// write to default parameters file
	if (iteration != -1) {
		return WriteNormalizeParameters(-1);
	}

	return 0;
}


bool Dssd::NormEnergyCheck(size_t, const DssdNormalizeEvent&) const {
	return true;
}


int Dssd::SideNormalize(
	TChain *chain,
	size_t side,
	size_t ref_strip,
	int
) {
	DssdNormalizeEvent events[2];
	events[0].SetupInput(chain);

	// energy graph fe:be or be:fe
	TGraph *ge;
	ge = new TGraph[Strip(side)];

	// total number of entries
	long long entries = chain->GetEntries();
	// 1/100 of total entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling for side %ld refer strip %ld   0%%", side, ref_strip);
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		chain->GetEntry(entry);
		// ignore multiple hit events
		if (events[0].front_hit != 1 || events[0].back_hit != 1) continue;

		unsigned short &fs = events[0].front_strip[0];
		unsigned short &bs = events[0].back_strip[0];
		double &fe = events[0].front_energy[0];
		double &be = events[0].back_energy[0];

		if (side == 0) {
			// jump if not refer strip
			if (bs != ref_strip) continue;
			if (!NormEnergyCheck(side, events[0])) continue;
			ge[fs].AddPoint(fe, NormEnergy(1, bs, be));
		} else {
			if (fs != ref_strip) continue;
			if (!NormEnergyCheck(side, events[0])) continue;
			ge[bs].AddPoint(be, NormEnergy(0, fs, fe));
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fitting
	std::cout << "side " << side << " normalize parameters.\n";
	for (size_t i = 0; i < Strip(side); ++i) {
		// only fits when over 10 points
		if (ge[i].GetN() > 10) {
			// fitting function
			TF1 energy_fit("efit", "pol1", 0, 60000);
			// set initial value
			energy_fit.SetParameter(0, 0.0);
			energy_fit.SetParameter(1, 1.0);
			// fit
			ge[i].Fit(&energy_fit, "QR+ ROB=0.9");
			// store the normalized parameters
			norm_params_[side][i][0] = energy_fit.GetParameter(0);
			norm_params_[side][i][1] = energy_fit.GetParameter(1);
		}
		// store the graph
		ge[i].Write(TString::Format("g%c%ld", "fb"[side], i));
		// print normalized paramters on screen
		std::cout << i << " " << norm_params_[side][i][0]
			<< ", " << norm_params_[side][i][1] << "\n";
	}

	// residual
	TGraph *res;
	res = new TGraph[Strip(side)];
	for (size_t i = 0; i < Strip(side); ++i) {
		int point = ge[i].GetN();
		double *gex = ge[i].GetX();
		double *gey = ge[i].GetY();
		for (int j = 0; j < point; ++j) {
			res[i].AddPoint(gex[j], NormEnergy(side, i, gex[j])-gey[j]);
		}
		res[i].Write(TString::Format("res%c%ld", "fb"[side], i));
	}

	// free memory
	delete[] ge;
	delete[] res;

	return 0;
}


int Dssd::NormalizeSides(TChain*, int) {
	std::cerr << "Error: Sides Normalize is not implemented yet.\n";
	return -1;
}


int Dssd::NormalizeResult(int iteration) {
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
	// input event
	DssdNormalizeEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-normalize-result-%s%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "normalized energy tree");
	// setup branches
	event.SetupOutput(&opt);

	// read normalize parameters from file
	if (ReadNormalizeParameters(iteration)) {
		std::cerr << "Error: Read normalize parameters from file failed.\n";
		return -1;
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// l/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf(
		"Writing normalized energy for %s in run %u   0%%",
		name_.c_str(), run_
	);
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		for (unsigned short i = 0; i < event.front_hit; ++i) {
			event.front_energy[i] =
				NormEnergy(0, event.front_strip[i], event.front_energy[i]);
		}
		for (unsigned short i = 0; i < event.back_hit; ++i) {
			event.back_energy[i] =
				NormEnergy(1, event.back_strip[i], event.back_energy[i]);
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// write output tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}


int Dssd::Normalize(
	unsigned int end_run,
	int iteration
) {
	// setup input chain
	TChain chain("tree", "chain");
	for (unsigned int i = run_; i <= end_run; ++i) {
		if (i == 628) continue;
		chain.AddFile(TString::Format(
			"%s%s%s-fundamental-%s%04u.root/tree",
			kGenerateDataPath,
			kFundamentalDir,
			name_.c_str(),
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			i
		));
	}

	// initialize normalized parameters
	for (size_t i = 0; i < FrontStrip(); ++i) {
		norm_params_[0][i][0] = 0.0;
		norm_params_[0][i][1] = 1.0;
	}
	for (size_t i = 0; i < BackStrip(); ++i) {
		norm_params_[1][i][0] = 0.0;
		norm_params_[1][i][1] = 1.0;
	}

	// read normalized parameters from file if in iteration mode
	if (iteration > 0 && ReadNormalizeParameters()) {
		std::cerr << "Error: read normalize parameters from file failed.\n";
		return -1;
	}

	// setup normalize root file
	TString normalize_file_name;
	normalize_file_name.Form(
		"%s%s%s-normalize-fit-%s%04u-%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_,
		end_run
	);
	// output file
	TFile opf(normalize_file_name, "recreate");

	if (NormalizeSides(&chain, iteration)) return -1;

	// close files
	opf.Close();

	// write parameters
	if (WriteNormalizeParameters(iteration)) {
		std::cerr << "Error: write normalize paramters to file failed.\n";
		return -1;
	}

	NormalizeStatistics statistics(run_, name_, tag_, end_run, iteration);
	statistics.Write();
	statistics.Print();
	return 0;
}


int Dssd::ShowNormalize() {
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
			<< input_file_name << "failed.\n";
		return -1;
	}
	// input event
	DssdFundamentalEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-norm-result-%s%04u.root",
		kGenerateDataPath,
		kShowDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// energy difference range
	const double diff_e_range = 5000;
	// 2D hisogram fe VS be
	TH2F fe_vs_be("hfbe", "fe:be", 1000, 0, 60000, 1000, 0, 60000);
	// 2D histogram fe-be VS fs
	TH2F diff_e_vs_fs(
		"hefs", "fe-be:fs",
		FrontStrip(), 0, FrontStrip(), 1000, 0, diff_e_range
	);
	// 2D histogram fe-be VS bs
	TH2F diff_e_vs_bs(
		"hebs", "fe-be:bs",
		BackStrip(), 0, BackStrip(), 1000, -diff_e_range, diff_e_range
	);
	// 2D histogram fe-be VS fe
	TH2F diff_e_vs_fe(
		"hefe", "fe-be:fe",
		1000, 0, 60000, 1000, 0, diff_e_range
	);
	// 2D histogram fe-be VS be
	TH2F diff_e_vs_be(
		"hebe", "fe-be:be",
		1000, 0, 60000, 1000, -diff_e_range, diff_e_range
	);
	// 1D histogram fe-be
	TH1F diff_e("hde", "fe-be", 1000, -diff_e_range, diff_e_range);
	// 1D histogram abs(fe-be)/(fe+be)
	TH1F relative_diff_e("hrde", "(fe-be)/(fe+be)", 1000, -1, 1);
	// 2D histogram (fe-be)/(fe+be) VS fe
	TH2F relative_diff_e_vs_fe(
		"hrdefe", "(fe-be)/(fe+be):fe",
		1000, 0, 60000, 1000, 0,  1
	);
	// output tree
	TTree opt("tree", "normalized energy tree");
	// setup branches
	event.SetupOutput(&opt);

	// read normalize parameters
	if (ReadNormalizeParameters()) {
		std::cerr << "Error: Read normalize parameters failed.\n";
		return -1;
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of total number of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Showing normalize result   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		for (unsigned short i = 0; i < event.front_hit; ++i) {
			event.front_energy[i] = NormEnergy(
				0, event.front_strip[i], event.front_energy[i]
			);
		}
		for (unsigned short i = 0; i < event.back_hit; ++i) {
			event.back_energy[i] = NormEnergy(
				1, event.back_strip[i], event.back_energy[i]
			);
		}
		opt.Fill();

		if (event.front_hit != 1 || event.back_hit != 1) continue;
		double fe = event.front_energy[0];
		double be = event.back_energy[0];
		fe_vs_be.Fill(be, fe);
		diff_e_vs_fs.Fill(event.front_strip[0], fabs(fe-be));
		diff_e_vs_bs.Fill(event.back_strip[0], fe-be);
		diff_e_vs_fe.Fill(fe, fabs(fe-be));
		diff_e_vs_be.Fill(be, fe-be);
		diff_e.Fill(fe-be);
		relative_diff_e.Fill((fe-be)/(fe+be));
		relative_diff_e_vs_fe.Fill(fe, fabs((fe-be)/(fe+be)));
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// save histograms
	fe_vs_be.Write();
	diff_e_vs_fs.Write();
	diff_e_vs_bs.Write();
	diff_e_vs_fe.Write();
	diff_e_vs_be.Write();
	diff_e.Write();
	relative_diff_e.Write();
	relative_diff_e_vs_fe.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}

//-----------------------------------------------------------------------------
//								merge
//-----------------------------------------------------------------------------

/// @brief merge front hit 1 and back hit 1 events
/// @param[in] fundamental input fundamental event
/// @param[in] merge output merged event
/// @param[in] diff_tolerance energy difference tolerance
/// @param[in] hde histogram to record energy difference
/// @param[in] statistics merge case statistics
///
void Merge11Event(
	const DssdFundamentalEvent &fundamental,
	DssdMergeEvent &merge,
	double diff_tolerance,
	TH1F *hde,
	MergeCaseStatistics* statistics
) {
	const unsigned short *fs = fundamental.front_strip;
	const unsigned short *bs = fundamental.back_strip;
	const double *fe = fundamental.front_energy;
	const double *be = fundamental.back_energy;
	// assume that only one particle
	double diff = fe[0] - be[0];
	hde[0].Fill(diff);
	++statistics[0].total;

	// fill events
	if (fabs(diff) < diff_tolerance) {
		merge.hit = 1;
		merge.case_tag = 0;
		merge.energy[0] = fe[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[0];
		merge.z[0] = 0.0;
		++statistics[0].merged;
	}
}


void Merge12Event(
	const DssdFundamentalEvent &fundamental,
	DssdMergeEvent &merge,
	double diff_tolerance,
	TH1F *hde,
	MergeCaseStatistics* statistics
) {
	const unsigned short *fs = fundamental.front_strip;
	const unsigned short *bs = fundamental.back_strip;
	const double *fe = fundamental.front_energy;
	const double *be = fundamental.back_energy;
	// assume that ajacent strips in back side
	double diff_a = fe[0] - be[0] - be[1];
	if (abs(bs[0]-bs[1]) == 1) hde[0].Fill(diff_a);
	++statistics[0].total;
	// assume that one strip can be negelect in back side
	double diff_b = fe[0] - be[0];
	if (abs(bs[0]-bs[1]) != 1) hde[1].Fill(diff_b);
	++statistics[1].total;

	// fill events
	if (fabs(diff_a) < diff_tolerance && abs(bs[0]-bs[1]) == 1) {
		// adjacent back strips event
		merge.hit = 1;
		merge.case_tag = 100;
		merge.energy[0] = fe[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[0] + bs[1]/(be[0]+be[1])*(bs[1]-bs[0]);
		merge.z[0] = 0.0;
		++statistics[0].merged;
	} else if (fabs(diff_b) < diff_tolerance) {
		merge.hit = 1;
		merge.case_tag = 101;
		merge.energy[0] = fe[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[0];
		merge.z[0] = 0.0;
		++statistics[1].merged;
	}
}


void Merge21Event(
	const DssdFundamentalEvent &fundamental,
	DssdMergeEvent &merge,
	double diff_tolerance,
	TH1F *hde,
	MergeCaseStatistics* statistics
) {
	const unsigned short *fs = fundamental.front_strip;
	const unsigned short *bs = fundamental.back_strip;
	const double *fe = fundamental.front_energy;
	const double *be = fundamental.back_energy;
	// assume that ajacent strips in front side
	double diff_a = fe[0] + fe[1] - be[0];
	if (abs(fs[0]-fs[1]) == 1) hde[0].Fill(diff_a);
	++statistics[0].total;
	// assume that one strip can be negelect in front side
	double diff_b = fe[0] - be[0];
	if (abs(fs[0]-fs[1]) != 1) hde[1].Fill(diff_b);
	++statistics[1].total;

	// fill events
	if (fabs(diff_a) < diff_tolerance && abs(fs[0]-fs[1]) == 1) {
		// adjacent event at the front side
		merge.hit = 1;
		merge.case_tag = 200;
		merge.energy[0] = be[0];
		merge.x[0] = fs[0] + fe[1]/(fe[0]+fe[1])*(fs[1]-fs[0]);
		merge.y[0] = bs[0];
		merge.z[0] = 0.0;
		++statistics[0].merged;
	} else if (fabs(diff_b) < diff_tolerance) {
		// small energy in the seoncd front strip can be neglected,
		// or the second back strip was thrown
		merge.hit = 1;
		merge.case_tag = 201;
		merge.energy[0] = fe[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[0];
		merge.z[0] = 0.0;
		++statistics[1].merged;
	}
}

void Merge22Event(
	const DssdFundamentalEvent &fundamental,
	DssdMergeEvent &merge,
	double diff_tolerance,
	TH1F *hde,
	MergeCaseStatistics* statistics
) {
	const unsigned short *fs = fundamental.front_strip;
	const unsigned short *bs = fundamental.back_strip;
	const double *fe = fundamental.front_energy;
	const double *be = fundamental.back_energy;
	// assume that there are two particles
	double diff_a11 = fe[0] - be[0];
	double diff_a12 = fe[1] - be[1];
	double diff_a21 = fe[0] - be[1];
	double diff_a22 = fe[1] - be[0];
	if (fabs(diff_a11)+fabs(diff_a12) < fabs(diff_a21)+fabs(diff_a22)) {
		hde[0].Fill(diff_a11);
		hde[0].Fill(diff_a12);
	} else {
		hde[0].Fill(diff_a21);
		hde[0].Fill(diff_a22);
	}
	++statistics[0].total;
	// assume that only one particle but ajacent strips in both sides
	double diff_b = fe[0] + fe[1] - be[0] - be[1];
	if (abs(fs[0]-fs[1])==1 && abs(bs[0]-bs[1])==1) {
		hde[1].Fill(diff_b);
	}
	++statistics[1].total;

	// fill events
	if (
		fabs(diff_a11) < diff_tolerance
		&& fabs(diff_a12) < diff_tolerance
	) {
		// tow particles event
		merge.hit = 2;
		merge.case_tag = 300;
		merge.energy[0] = fe[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[0];
		merge.z[0] = 0.0;
		merge.energy[1] = fe[1];
		merge.x[1] = fs[1];
		merge.y[1] = bs[1];
		merge.z[1] = 0.0;
		++statistics[0].merged;
	} else if (
		fabs(diff_a21) < diff_tolerance
		&& fabs(diff_a22) < diff_tolerance
	) {
		// two particles event
		merge.hit = 2;
		merge.case_tag = 300;
		if (fe[0] > fe[1]) {
			merge.energy[0] = fe[0];
			merge.x[0] = fs[0];
			merge.y[0] = bs[1];
			merge.z[0] = 0.0;
			merge.energy[1] = fe[1];
			merge.x[1] = fs[1];
			merge.y[1] = bs[0];
			merge.z[1] = 0.0;
		} else {
			merge.energy[0] = fe[1];
			merge.x[0] = fs[1];
			merge.y[0] = bs[0];
			merge.z[0] = 0.0;
			merge.energy[1] = fe[0];
			merge.x[1] = fs[0];
			merge.y[1] = bs[1];
			merge.z[1] = 0.0;
		}
		++statistics[0].merged;
	} else if (
		fabs(diff_b) < 1.414 * diff_tolerance
		&& abs(fs[0] - fs[1]) == 1
		&& abs(bs[0] - bs[1]) == 1
	) {
		// one particle event, two adacent strips
		merge.hit = 1;
		merge.case_tag = 301;
		merge.energy[0] = fe[0] + fe[1];
		merge.x[0] = fs[0] + fe[1]/(fe[0]+fe[1])*(fs[1]-fs[0]);
		merge.y[0] = bs[0] + be[1]/(be[0]+be[1])*(bs[1]-bs[0]);
		merge.z[0] = 0.0;
		++statistics[1].merged;
	}
}


void Merge23Event(
	const DssdFundamentalEvent &fundamental,
	DssdMergeEvent &merge,
	double diff_tolerance,
	TH1F *hde,
	MergeCaseStatistics* statistics
) {
	const unsigned short *fs = fundamental.front_strip;
	const unsigned short *bs = fundamental.back_strip;
	const double *fe = fundamental.front_energy;
	const double *be = fundamental.back_energy;
	// assume that there are two particles
	// 1. bs[0] and bs[1] are adjacent strips
	double diff_a11 = fe[0] - be[0] - be[1];
	double diff_a12 = fe[1] - be[2];
	// 2. bs[0] and bs[2] are adjacent strips
	double diff_a21 = fe[0] - be[0] - be[2];
	double diff_a22 = fe[1] - be[1];
	// 3. bs[1] and bs[2] are adjacent strips
	double diff_a31 = fe[0] - be[0];
	double diff_a32 = fe[1] - be[1] - be[2];
	double diff_a33 = fe[0] - be[1] - be[2];
	double diff_a34 = fe[1] - be[0];
	// fill to histogram
	if (abs(bs[0]-bs[1]) == 1) {
		hde[0].Fill(diff_a11);
		hde[0].Fill(diff_a12);
	} else if (abs(bs[0]-bs[2]) == 1) {
		hde[0].Fill(diff_a21);
		hde[0].Fill(diff_a22);
	} else if (abs(bs[1]-bs[2]) == 1) {
		if (fabs(diff_a31)+fabs(diff_a32) < fabs(diff_a33)+fabs(diff_a34)) {
			hde[0].Fill(diff_a31);
			hde[0].Fill(diff_a32);
		} else {
			hde[0].Fill(diff_a33);
			hde[0].Fill(diff_a34);
		}
	}
	++statistics[0].total;

	// fill events
	if (
		abs(bs[0] - bs[1]) == 1
		&& fabs(diff_a11) < diff_tolerance
		&& fabs(diff_a12) < diff_tolerance
	) {
		// particle 1: f0b0b1,
		// particle 2: f1b2
		merge.hit = 2;
		merge.case_tag = 400;
		merge.energy[0] = fe[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[0] + be[1]/(be[0]+be[1])*(bs[1]-bs[0]);
		merge.z[0] = 0.0;
		merge.energy[1] = fe[1];
		merge.x[1] = fs[1];
		merge.y[1] = bs[2];
		merge.z[1] = 0.0;
		++statistics[0].merged;
	} else if (
		abs(bs[0] - bs[2]) == 1
		&& fabs(diff_a21) < diff_tolerance
		&& fabs(diff_a22) < diff_tolerance
	) {
		// particle 1: f0b0b2
		// particle 2: f1b1
		merge.hit = 2;
		merge.case_tag = 400;
		merge.energy[0] = fe[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[0] + bs[2]/(be[0]+be[2])*(bs[2]-bs[0]);
		merge.z[0] = 0.0;
		merge.energy[1] = fe[1];
		merge.x[1] = fs[1];
		merge.y[1] = bs[1];
		merge.z[1] = 0.0;
		++statistics[0].merged;
	} else if (
		abs(bs[1] - bs[2]) == 1
		&& fabs(diff_a31) < diff_tolerance
		&& fabs(diff_a32) < diff_tolerance
		&& fabs(diff_a31) + fabs(diff_a32) < fabs(diff_a33) + fabs(diff_a34)
	) {
		// particle 1: f0b0
		// particle 2: f1b1b2
		merge.hit = 2;
		merge.case_tag = 400;
		merge.energy[0] = fe[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[0];
		merge.z[0] = 0.0;
		merge.energy[1] = fe[1];
		merge.x[1] = fs[1];
		merge.y[1] = bs[1] + bs[2]/(be[1]+be[2])*(bs[2]-bs[1]);
		merge.z[1] = 0.0;
		++statistics[0].merged;
	} else if (
		abs(bs[1] - bs[2]) == 1
		&& fabs(diff_a33) < diff_tolerance
		&& fabs(diff_a34) < diff_tolerance
		&& fabs(diff_a33) + fabs(diff_a34) < fabs(diff_a31) + fabs(diff_a32)
	) {
		// particle 1: f0b1b2
		// particle 2: f1b0
		merge.hit = 2;
		merge.case_tag = 400;
		merge.energy[0] = fe[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[1] + bs[2]/(be[1]+be[2])*(bs[2]-bs[1]);
		merge.z[0] = 0.0;
		merge.energy[1] = fe[1];
		merge.x[1] = fs[1];
		merge.y[1] = bs[0];
		merge.z[1] = 0.0;
		++statistics[0].merged;
	}
}


void Merge32Event(
	const DssdFundamentalEvent &fundamental,
	DssdMergeEvent &merge,
	double diff_tolerance,
	TH1F *hde,
	MergeCaseStatistics* statistics
) {
	const unsigned short *fs = fundamental.front_strip;
	const unsigned short *bs = fundamental.back_strip;
	const double *fe = fundamental.front_energy;
	const double *be = fundamental.back_energy;
	// assume that there are two particles
	// 1. fs[0] and fs[1] are adjacent strips
	double diff_a11 = fe[0] + fe[1] - be[0];
	double diff_a12 = fe[2] - be[1];
	// 2. fs[0] and fs[2] are adjacent strips
	double diff_a21 = fe[0] + fe[2] - be[0];
	double diff_a22 = fe[1] - be[1];
	// 3. fs[1] and fs[2] are adjacent strips
	double diff_a31 = fe[0] - be[0];
	double diff_a32 = fe[1] + fe[2] - be[1];
	double diff_a33 = fe[0] - be[1];
	double diff_a34 = fe[1] + fe[2] - be[0];
	// fill to histogram
	if (abs(fs[0]-fs[1]) == 1) {
		hde[0].Fill(diff_a11);
		hde[0].Fill(diff_a12);
	} else if (abs(fs[0]-fs[2]) == 1) {
		hde[0].Fill(diff_a21);
		hde[0].Fill(diff_a22);
	} else if (abs(fs[1]-fs[2]) == 1) {
		if (fabs(diff_a31)+fabs(diff_a32) < fabs(diff_a33)+fabs(diff_a34)) {
			hde[0].Fill(diff_a31);
			hde[0].Fill(diff_a32);
		} else {
			hde[0].Fill(diff_a33);
			hde[0].Fill(diff_a34);
		}
	}
	++statistics[0].total;

	// fill events
	if (
		abs(fs[0] - fs[1]) == 1
		&& fabs(diff_a11) < diff_tolerance
		&& fabs(diff_a12) < diff_tolerance
	) {
		// particle 1: f0f1b0,
		// particle 2: f2b1
		merge.hit = 2;
		merge.case_tag = 500;
		merge.energy[0] = be[0];
		merge.x[0] = fs[0] + fs[1]/(fe[0]+fe[1])*(fs[1]-fs[0]);
		merge.y[0] = bs[0];
		merge.z[0] = 0.0;
		merge.energy[1] = be[1];
		merge.x[1] = fs[2];
		merge.y[1] = bs[1];
		merge.z[1] = 0.0;
		++statistics[0].merged;
	} else if (
		abs(fs[0] - fs[2]) == 1
		&& fabs(diff_a21) < diff_tolerance
		&& fabs(diff_a22) < diff_tolerance
	) {
		// particle 1: f0f2b0
		// particle 2: f1b1
		merge.hit = 2;
		merge.case_tag = 500;
		merge.energy[0] = be[0];
		merge.x[0] = fs[0] + fs[2]/(fe[0]+fe[2])*(fs[2]-fs[0]);
		merge.y[0] = bs[0];
		merge.z[0] = 0.0;
		merge.energy[1] = be[1];
		merge.x[1] = fs[1];
		merge.y[1] = bs[1];
		merge.z[1] = 0.0;
		++statistics[0].merged;
	} else if (
		abs(fs[1] - fs[2]) == 1
		&& fabs(diff_a31) < diff_tolerance
		&& fabs(diff_a32) < diff_tolerance
		&& fabs(diff_a31) + fabs(diff_a32) < fabs(diff_a33) + fabs(diff_a34)
	) {
		// particle 1: f0b0
		// particle 2: f1f2b1
		merge.hit = 2;
		merge.case_tag = 500;
		merge.energy[0] = be[0];
		merge.x[0] = fs[0];
		merge.y[0] = bs[0];
		merge.z[0] = 0.0;
		merge.energy[1] = be[1];
		merge.x[1] = fs[1] + fs[2]/(fe[1]+fe[2])*(fs[2]-fs[1]);
		merge.y[1] = bs[1];
		merge.z[1] = 0.0;
		++statistics[0].merged;
	} else if (
		abs(fs[1] - fs[2]) == 1
		&& fabs(diff_a33) < diff_tolerance
		&& fabs(diff_a34) < diff_tolerance
		&& fabs(diff_a33) + fabs(diff_a34) < fabs(diff_a31) + fabs(diff_a32)
	) {
		// particle 1: f0b1
		// particle 2: f1f2b0
		merge.hit = 2;
		merge.case_tag = 500;
		merge.energy[0] = be[1];
		merge.x[0] = fs[0];
		merge.y[0] = bs[1];
		merge.z[0] = 0.0;
		merge.energy[1] = be[0];
		merge.x[1] = fs[1] + fs[2]/(fe[1]+fe[2])*(fs[2]-fs[1]);
		merge.y[1] = bs[0];
		merge.z[1] = 0.0;
		++statistics[0].merged;
	}
}


// void Merge33Event(
// 	const DssdFundamentalEvent &fundamental,
// 	DssdMergeEvent &merge,
// 	double diff,
// 	TH1F *hde,
// 	MergeCaseStatistics* statistics
// ) {

// }


/// @brief recursive selection sort particles in merge event by energy
/// @param[inout] merge merge event to sort
/// @param[in] start start index of particles to sort,
///		assume that particles ahead are sorted
/// @warning this function can only sort particles under 2
void SortMergeEvent(DssdMergeEvent &merge) {
	if (merge.hit <= 1) return;
	if (merge.hit == 2 && merge.energy[0] < merge.energy[1]) {
		double tmp;
		// swap energy
		tmp = merge.energy[0];
		merge.energy[0] = merge.energy[1];
		merge.energy[1] = tmp;
		// swap x
		tmp = merge.x[0];
		merge.x[0] = merge.x[1];
		merge.x[1] = tmp;
		// swap y
		tmp = merge.y[0];
		merge.y[0] = merge.y[1];
		merge.y[1] = tmp;
		// swap z
		tmp = merge.z[0];
		merge.z[0] = merge.z[1];
		merge.z[1] = tmp;
	}
	return;
}



int Dssd::Merge(double energy_diff) {
	// input file name
	TString fundamental_file_name;
	fundamental_file_name.Form(
		"%s%s%s-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// input file
	TFile ipf(fundamental_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< fundamental_file_name << " failed.\n";
	}
	// input event
	DssdFundamentalEvent fundamental_event;
	// setup input branches
	fundamental_event.SetupInput(ipt);
	// for convenient
	unsigned short &fhit = fundamental_event.front_hit;
	unsigned short &bhit = fundamental_event.back_hit;
	unsigned short *fs = fundamental_event.front_strip;
	unsigned short *bs = fundamental_event.back_strip;
	double *fe = fundamental_event.front_energy;
	double *be = fundamental_event.back_energy;

	// output file name
	TString merge_file_name;
	merge_file_name.Form(
		"%s%s%s-merge-%s%04u.root",
		kGenerateDataPath,
		kMergeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(merge_file_name, "recreate");
	// difference of front and back side energy
	TH1F hde[]{
		TH1F("hde11", "fe[0]-be[0]", 400, -2000, 2000),
		TH1F("hde12a", "fe[0]+fe[1]-be[0]", 400, -2000, 2000),
		TH1F("hde12b", "fe[0]-be[0]", 400, -2000, 2000),
		TH1F("hde21a", "fe[0]-be[0]-be[1]", 400, -2000, 2000),
		TH1F("hde21b", "fe[0]-be[0]", 400, -2000, 2000),
		TH1F("hde22a", "fe[a]-be[a]", 400, -2000, 2000),
		TH1F("hde22b", "fe-be", 400, -2000, 2000),
		TH1F("hde32a", "fe[a]-be[a], fe[b]-be[b]-be[c]", 400, -2000, 2000),
		TH1F("hde23a", "fe[a]-be[a], fe[b]+fe[c]-be[b]", 400, -2000, 2000)
	};
	// output tree
	TTree opt("tree", "tree of merged events");
	// output event
	DssdMergeEvent merge_event;
	// setup output branches
	merge_event.SetupOutput(&opt);

	// read normalized parameters
	if (ReadNormalizeParameters()) {
		std::cerr << "Error: Read normalize parameters failed.\n";
		return -1;
	}

	// summary statistics
	MergeStatistics statistics(run_, name_, tag_);
	// statistics for each case
	MergeCaseStatistics case_statistics[]{
		MergeCaseStatistics(run_, name_, tag_, "f1b1", 0, energy_diff),
		MergeCaseStatistics(run_, name_, tag_, "f1b2", 0, energy_diff),
		MergeCaseStatistics(run_, name_, tag_, "f1b2", 1, energy_diff),
		MergeCaseStatistics(run_, name_, tag_, "f2b1", 0, energy_diff),
		MergeCaseStatistics(run_, name_, tag_, "f2b1", 1, energy_diff),
		MergeCaseStatistics(run_, name_, tag_, "f2b2", 0, energy_diff),
		MergeCaseStatistics(run_, name_, tag_, "f2b2", 1, energy_diff),
		MergeCaseStatistics(run_, name_, tag_, "f2b3", 0, energy_diff),
		MergeCaseStatistics(run_, name_, tag_, "f3b2", 0, energy_diff)
	};
	// function pointer to function that merge differenct case events
	void(*MergeFunction[])(
		const DssdFundamentalEvent&,
		DssdMergeEvent&,
		double,
		TH1F*,
		MergeCaseStatistics*
	) = {
		Merge11Event,
		Merge12Event, Merge21Event,
		Merge22Event,
		Merge23Event, Merge32Event,
		// Merge33Event
	};
	// start index of different cases in hde and statistics
	int case_index[]{0, 1, 3, 5, 7, 8};

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
		// initialize
		merge_event.hit = 0;

		// normalize energy
		for (unsigned short i = 0; i < fhit; ++i) {
			fe[i] = NormEnergy(0, fs[i], fe[i]);
		}
		for (unsigned short i = 0; i < bhit; ++i) {
			be[i] = NormEnergy(1, bs[i], be[i]);
		}
		if (fhit > 0 && bhit > 0) ++statistics.total;
		// A trick to get case number by fhit and bhit.
		// This only works with following cases with hits:
		// 0. f1-b1; 1. f1-b2; 2. f2-b1; 3. f2-b2
		// 4. f2-b3; 5. f3-b2; 6. f3-b3
		int case_num = 2*fhit  + bhit - 3;
		// check case is in the 7 cases memtioned
		if (case_num < 0 || case_num >= 6 || abs(fhit-bhit) > 1) {
			opt.Fill();
			continue;
		}
		// call merge function to fill merge_event
		MergeFunction[case_num](
			fundamental_event,
			merge_event,
			energy_diff,
			hde + case_index[case_num],
			case_statistics + case_index[case_num]
		);
		for (unsigned short i = 0; i < merge_event.hit; ++i) {
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

	// save and close files
	for (auto &hist : hde) hist.Write();
	opt.Write();
	opf.Close();
	ipf.Close();

	// summarize one hit merged events
	for (size_t i = 0; i < 5; ++i) {
		statistics.one_hit += case_statistics[i].merged;
	}
	statistics.one_hit += case_statistics[6].merged;
	// summarize two hit merged events
	statistics.two_hit += case_statistics[5].merged;
	statistics.two_hit += case_statistics[7].merged;
	statistics.two_hit += case_statistics[8].merged;
	statistics.merged =
		statistics.one_hit+ statistics.two_hit + statistics.three_hit;
	// save and show statistics
	for (auto &s : case_statistics) {
		s.Write();
		s.Print();
	}
	statistics.Write();
	statistics.Print();

	return 0;
}


}