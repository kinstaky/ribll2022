#include "include/detector/dssd.h"

#include <array>
#include <exception>

#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TMarker.h>

#include "include/event/dssd_event.h"
#include "include/event/filter_event.h"
#include "include/event/particle_event.h"
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


	for (int i = 0; i < fundamental_event.front_hit; ++i) {
		fundamental_event.front_strip[i] = front_events[i].strip;
		fundamental_event.front_time[i] = front_events[i].time - trigger_time;
		fundamental_event.front_energy[i] = front_events[i].energy;
		fundamental_event.cfd_flag |=
			front_events[i].cfd_flag ? (1 << i) : 0;
		fundamental_event.front_decode_entry[i] = front_events[i].decode_entry;
	}
	for (int i = 0; i < fundamental_event.back_hit; ++i) {
		fundamental_event.back_strip[i] = back_events[i].strip;
		fundamental_event.back_time[i] = back_events[i].time - trigger_time;
		fundamental_event.back_energy[i] = back_events[i].energy;
		fundamental_event.cfd_flag |=
			back_events[i].cfd_flag ? (1 << (i + 8)) : 0;
		fundamental_event.back_decode_entry[i] = back_events[i].decode_entry;
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
		// "%s%s%s-wk-b0%s.txt",
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
		for (size_t j = 0; j < 4; ++j) {
			fin >> norm_params_[0][i][j];
		}
	}
	//  read back strip number
	fin >> strip_num;
	// read normalized parameters for back strips
	for (size_t i = 0; i < strip_num; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			fin >> norm_params_[1][i][j];
		}
	}
	std::string title, version;
	fin >> title >> version;
	std::cout << "Read normalize parameters "
		<< title << " " << version << "\n";
	// close file
	fin.close();

	return 0;
}


int Dssd::WriteNormalizeParameters(int iteration) {
	// setup file name
	TString file_names[2];
	file_names[0].Form(
		"%s%s%s.txt",
		kGenerateDataPath, kNormalizeDir, name_.c_str()
	);
	file_names[1].Form(
		"%s%s%s-%d.txt",
		kGenerateDataPath, kNormalizeDir, name_.c_str(), iteration
	);
	for (int i = 0; i < (iteration >= 0 ? 2 : 1); ++i) {
		// open output file
		std::ofstream fout(file_names[i].Data());
		if (!fout.good()) {
			std::cerr << "Error: open normalize parameters file "
				<< file_names[i] << " failed.\n";
			return -1;
		}
		// write front strip number
		fout << FrontStrip() << "\n";
		// write normalized paramters for front strips
		for (size_t i = 0; i < FrontStrip(); ++i) {
			fout << norm_params_[0][i][0]<< " "
				<< norm_params_[0][i][1]<< " "
				<< norm_params_[0][i][2]<< " "
				<< norm_params_[0][i][3] << "\n";
		}
		// write back strip number
		fout << BackStrip() << "\n";
		// write normalized parameters for back strips
		for (size_t i = 0; i < BackStrip(); ++i) {
			fout << norm_params_[1][i][0] << " "
				<< norm_params_[1][i][1] << " "
				<< norm_params_[1][i][2] << " "
				<< norm_params_[1][i][3] << "\n";
		}
		// add some information about the paramters
		fout << "version " << iteration << "\n";
		// close file
		fout.close();
	}
	return 0;
}


int Dssd::StripsNormalize(
	TChain *chain,
	size_t side,
	size_t ref_start,
	size_t ref_end,
	size_t norm_start,
	size_t norm_end,
	int iteration
) {
	// event to normalize
	DssdNormalizeEvent event;
	// filter event
	unsigned short filter;
	event.SetupInput(chain);
	if (iteration > 0) {
		chain->SetBranchAddress("filter.flag", &filter);
	}

	// energy graph fe:be or be:fe
	TGraph *ge = new TGraph[Strip(side)];

	// total number of entries
	long long entries = chain->GetEntries();
	// 1/100 of total entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf(
		"Filling for side %ld [%ld, %ld) reference strips [%ld, %ld)   0%%",
		side, norm_start, norm_end, ref_start, ref_end
	);
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		chain->GetEntry(entry);
		// ignore multiple hit events
		if (event.front_hit != 1 || event.back_hit != 1) continue;
		if (iteration > 0 && filter != 0x1) continue;

		unsigned short &fs = event.front_strip[0];
		unsigned short &bs = event.back_strip[0];
		double &fe = event.front_energy[0];
		double &be = event.back_energy[0];
		if (fe > 20000 || be > 20000) continue;

		if (side == 0) {
			// jump if not reference strips
			if (bs < ref_start || bs >= ref_end) continue;
			// jump if not normalize strips
			if (fs < norm_start || fs >= norm_end) continue;
			// jump if has normalized
			if (has_normalized_[side][fs]) continue;
			// fill to graph
			ge[fs].AddPoint(fe, NormEnergy(1, bs, be));
		} else {
			// jump if not reference strips
			if (fs < ref_start || fs >= ref_end) continue;
			// jump if not normalize strips
			if (bs < norm_start || bs >= norm_end) continue;
			// jump if has normalized
			if (has_normalized_[side][bs]) continue;
			// fill to graph
			ge[bs].AddPoint(be, NormEnergy(0, fs, fe));
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fitting
	std::cout << "side " << side << " normalize parameters.\n";
	for (size_t i = norm_start; i < norm_end; ++i) {
		if (has_normalized_[side][i]) continue;
		// only fits when over 10 points
		if (ge[i].GetN() > 5) {
			// fitting function
			TF1 energy_fit("efit", "pol1", 0, 60000);
			// set initial value
			// energy_fit.FixParameter(0, 0.0);
			energy_fit.SetParameter(0, 0.0);
			energy_fit.SetParameter(1, 1.0);
			// energy_fit.SetParameter(2, 0.0);
			// fit
			ge[i].Fit(&energy_fit, "QR+ ROB=0.8");
			// store the normalized parameters
			// norm_params_[side][i][0] = energy_fit.GetParameter(0);
			norm_params_[side][i][1] = energy_fit.GetParameter(0);
			norm_params_[side][i][2] = energy_fit.GetParameter(1);
			// norm_params_[side][i][3] = energy_fit.GetParameter(2);
		}
		// store the graph
		ge[i].Write(TString::Format("g%c%ld", "fb"[side], i));
		// set as normalized
		has_normalized_[side][i] = true;
		// print normalized paramters on screen
		std::cout << i
			<< " " << norm_params_[side][i][0]
			<< ", " << norm_params_[side][i][1]
			<< ", " << norm_params_[side][i][2]
			<< ", " << norm_params_[side][i][3] << "\n";
	}

	// residual
	TGraph *res = new TGraph[Strip(side)];
	for (size_t i = norm_start; i < norm_end; ++i) {
		if (has_normalized_[side][i]) continue;
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

	// setup filter chain in iteration mode
	// filter chain for filtering events
	TChain filter_chain("filter", "filter chain");
	if (iteration > 0) {
		for (unsigned int i = run_; i <= end_run; ++i) {
			if (i == 628) continue;
			filter_chain.AddFile(TString::Format(
				"%s%s%s-filter-%s%04u-%d.root/tree",
				kGenerateDataPath,
				kNormalizeDir,
				name_.c_str(),
				tag_.empty() ? "" : (tag_ + "-").c_str(),
				i,
				iteration
			));
		}
		chain.AddFriend(&filter_chain);
	}

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

	// setup normalize root file
	TString normalize_file_name;
	normalize_file_name.Form(
		"%s%s%s-normalize-fit-%s%04u-%04u%s.root",
		kGenerateDataPath,
		kNormalizeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_,
		end_run,
		iteration == 0 ? "" : ("-"+std::to_string(iteration)).c_str()
	);
	// output file
	TFile opf(normalize_file_name, "recreate");

	// initialize
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < 64; ++j) {
			has_normalized_[i][j] = false;
		}
	}

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


int Dssd::NormalizeResult(int iteration) {
	// input file
	TFile *ipf = nullptr;
	// input tree
	TTree *ipt = nullptr;
	// input event
	DssdFundamentalEvent event;
	FilterEvent filter_event;

	// setup input file, tree, and branches in different conditions
	if (iteration == 0) {
		// input file name
		TString fundamental_file_name = TString::Format(
			"%s%s%s-fundamental-%s%04u.root",
			kGenerateDataPath,
			kFundamentalDir,
			name_.c_str(),
			tag_.empty() ? "" : (tag_+"-").c_str(),
			run_
		);
		// input file
		ipf = new TFile(fundamental_file_name, "read");
		// input tree
		ipt = (TTree*)ipf->Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< fundamental_file_name << " failed.\n";
			return -1;
		}
		// setup input branches
		event.SetupInput(ipt);
	} else {
		// input filter file name
		TString filter_file_name = TString::Format(
			"%s%s%s-normalize-filter-%d-%s%04u.root",
			kGenerateDataPath,
			kFilterDir,
			name_.c_str(),
			iteration,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_
		);
		// input file
		ipf = new TFile(filter_file_name, "read");
		// input tree
		ipt = (TTree*)ipf->Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< filter_file_name << " failed.\n";
			return -1;
		}
		// setup input branches
		filter_event.SetupInput(ipt);
	}

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-result-%s%s%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
		name_.c_str(),
		iteration == 0 ? "" : (std::to_string(iteration)+"-").c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// energy difference range
	const double diff_e_range = 2000;
	// 2D hisogram fe VS be
	TH2F fe_be("hfbe", "fe:be", 1000, 0, 60000, 1000, 0, 60000);
	// 2D histogram fe-be VS fs
	TH2F de_fs(
		"hefs", "fe-be:fs",
		FrontStrip(), 0, FrontStrip(), 1000, 0, diff_e_range
	);
	// 2D histogram fe-be VS bs
	TH2F de_bs(
		"hebs", "fe-be:bs",
		BackStrip(), 0, BackStrip(), 1000, -diff_e_range, diff_e_range
	);
	// 2D histogram fe-be VS fe
	TH2F de_fe(
		"hefe", "fe-be:fe",
		1000, 0, 60000, 1000, 0, diff_e_range
	);
	// 2D histogram fe-be VS be
	TH2F de_be(
		"hebe", "fe-be:be",
		1000, 0, 60000, 1000, -diff_e_range, diff_e_range
	);
	// 1D histogram fe-be
	TH1F de("hde", "fe-be", 1000, -diff_e_range, diff_e_range);
	// 1D histogram abs(fe-be)/(fe+be)
	TH1F rde("hrde", "(fe-be)/(fe+be)", 1000, -1, 1);
	// 2D histogram (fe-be)/(fe+be) VS fe
	TH2F rde_fe(
		"hrdefe", "(fe-be)/(fe+be):fe",
		1000, 0, 60000, 1000, 0,  1
	);
	// output tree
	TTree opt("tree", "normalized energy tree");
	// setup branches
	event.SetupOutput(&opt);

	// read normalize parameters from file
	if (ReadNormalizeParameters()) {
		std::cerr << "Error: Read normalize parameters from file failed.\n";
		return -1;
	}
	// read time normalize parameters from file
	// if (ReadNormalizeTimeParameters()) {
	// 	std::cerr <<
	// 		"Error: Read normalize time parameters from file failed.\n";
	// 	return -1;
	// }

	// total number of entries
	long long entries = ipt->GetEntries();
	// l/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf(
		"Writing normalized result for %s in run %u   0%%",
		name_.c_str(), run_
	);
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		if (iteration == 0) {
			for (int i = 0; i < event.front_hit; ++i) {
				event.front_energy[i] =
					NormEnergy(0, event.front_strip[i], event.front_energy[i]);
				event.front_time[i] -= norm_time_params_[0][event.front_strip[i]];
				event.front_fundamental_index[i] = i;
			}
			for (int i = 0; i < event.back_hit; ++i) {
				event.back_energy[i] =
					NormEnergy(1, event.back_strip[i], event.back_energy[i]);
				event.back_time[i] -= norm_time_params_[1][event.back_strip[i]];
				event.back_fundamental_index[i] = i;
			}
			event.Sort();

			// energy cut
			// for (int i = event.front_hit-1; i >= 0; --i) {
			// 	if (event.front_energy[i] < 1200.0) --event.front_hit;
			// }
			// for (int i = event.back_hit-1; i >= 0; --i) {
			// 	if (event.back_energy[i] < 1200.0) --event.back_hit;
			// }

			if (event.front_hit == 1 && event.back_hit == 1) {
				double fe = event.front_energy[0];
				double be = event.back_energy[0];
				fe_be.Fill(be, fe);
				de_fs.Fill(event.front_strip[0], fabs(fe-be));
				de_bs.Fill(event.back_strip[0], fe-be);
				de_fe.Fill(fe, fabs(fe-be));
				de_be.Fill(be, fe-be);
				de.Fill(fe-be);
				rde.Fill((fe-be)/(fe+be));
				rde_fe.Fill(fe, fabs((fe-be)/(fe+be)));
			}
		} else {
			event.front_hit = filter_event.num;
			event.back_hit = filter_event.num;
			for (int i = 0; i < filter_event.num; ++i) {
				event.front_strip[i] = filter_event.front_strip[i];
				event.back_strip[i] = filter_event.back_strip[i];
				event.front_energy[i] = NormEnergy(
					0,
					filter_event.front_strip[i],
					filter_event.front_energy[i]
				);
				event.back_energy[i] = NormEnergy(
					1,
					filter_event.back_strip[i],
					filter_event.back_energy[i]
				);
				double fe = event.front_energy[i];
				double be = event.back_energy[i];
				fe_be.Fill(be, fe);
				de_fs.Fill(event.front_strip[i], fabs(fe-be));
				de_bs.Fill(event.back_strip[i], fe-be);
				de_fe.Fill(fe, fabs(fe-be));
				de_be.Fill(be, fe-be);
				de.Fill(fe-be);
				rde.Fill((fe-be)/(fe+be));
				rde_fe.Fill(fe, fabs((fe-be)/(fe+be)));
			}
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms
	fe_be.Write();
	de_fs.Write();
	de_bs.Write();
	de_fe.Write();
	de_be.Write();
	de.Write();
	rde.Write();
	rde_fe.Write();
	// write output tree
	opt.Write();
	// close files
	opf.Close();
	ipf->Close();

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
		for (int i = 0; i < event.front_hit; ++i) {
			event.front_energy[i] = NormEnergy(
				0, event.front_strip[i], event.front_energy[i]
			);
		}
		for (int i = 0; i < event.back_hit; ++i) {
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


std::unique_ptr<TCutG> Dssd::ReadCut(
	const std::string &dir,
	const std::string &name
) const {
	// open cut file
	std::ifstream fin(TString::Format(
		"%s%scut/%s.txt",
		kGenerateDataPath,
		dir.c_str(),
		name.c_str()
	));
	if (!fin.good()) return nullptr;
	// cut
	std::unique_ptr<TCutG> cut = std::make_unique<TCutG>();
	// point index
	int point;
	// point position
	double x, y;
	// loop to read points
	while (fin.good()) {
		fin >> point >> x >> y;
		cut->SetPoint(point, x, y);
	}
	// close file
	fin.close();
	return cut;
}

//-----------------------------------------------------------------------------
//								merge
//-----------------------------------------------------------------------------

int Dssd::CutBeamThreshold() {
	std::cerr << "Error: Dssd::CutBeamThreshold is not implemented yet.\n";
	return -1;
}


/// @brief swap particles in merge event
/// @param[inout] merge merge event
/// @param[in] i first index to swap
/// @param[in] j second index to swap
///
void SwapMergeEvent(DssdMergeEvent &merge, size_t i, size_t j) {
	// swap flag
	unsigned int tf = merge.flag[i];
	merge.flag[i] = merge.flag[j];
	merge.flag[j] = tf;
	// swap merge flag
	unsigned short tmf = merge.merge_tag[i];
	merge.merge_tag[i] = merge.merge_tag[j];
	merge.merge_tag[j] = tmf;
	// swap energy
	double tmp = merge.energy[i];
	merge.energy[i] = merge.energy[j];
	merge.energy[j] = tmp;
	// swap time
	tmp = merge.time[i];
	merge.time[i] = merge.time[j];
	merge.time[j] = tmp;
	// swap x
	tmp = merge.x[i];
	merge.x[i] = merge.x[j];
	merge.x[j] = tmp;
	// swap y
	tmp = merge.y[i];
	merge.y[i] = merge.y[j];
	merge.y[j] = tmp;
	// swap z
	tmp = merge.z[i];
	merge.z[i] = merge.z[j];
	merge.z[j] = tmp;
	// swap time flag
	int ttf = merge.time_flag[i];
	merge.time_flag[i] = merge.time_flag[j];
	merge.time_flag[j] = ttf;
}


/// @brief sort particles in merge event by energy
/// @param[inout] merge merge event to sort
///
void Dssd::SortMergeEvent(DssdMergeEvent &merge) {
	if (merge.hit <= 1) return;
	else if (merge.hit == 2 && merge.energy[0] < merge.energy[1]) {
		SwapMergeEvent(merge, 0, 1);
	} else if (merge.hit == 3) {
		if (
			merge.energy[0] < merge.energy[1]
			&& merge.energy[0] < merge.energy[2]
		) {
			SwapMergeEvent(merge, 0, 2);
		}
		if (merge.energy[1] < merge.energy[2]) {
			SwapMergeEvent(merge, 1, 2);
		}
		if (merge.energy[0] < merge.energy[1]) {
			SwapMergeEvent(merge, 0, 1);
		}
	}
	return;
}


/// @brief search front ajacent strip to specific strip
/// @param[in] event normalize result event
/// @param[in] strip strip to search
/// @param[in] flag used strip flag
/// @returns index of adjacent strip if found, -1 otherwise
///
int Dssd::SearchFrontAdjacentStrips(
	const DssdFundamentalEvent &event,
	unsigned short strip,
	unsigned short flag
) {
	for (int i = 0; i < event.front_hit; ++i) {
		if (event.front_strip[i] == strip) continue;
		if ((flag & (1 << i)) != 0) continue;
		if (abs(event.front_strip[i] - strip) != 1) continue;
		// returns the first adjacent strip, since it has the largest energy
		return int(i);
	}
	return -1;
}


/// @brief search back ajacent strip to specific strip
/// @param[in] event normalize result event
/// @param[in] strip strip to search
/// @param[in] flag used strip flag
/// @returns index of adjacent strip if found, -1 otherwise
///
int Dssd::SearchBackAdjacentStrips(
	const DssdFundamentalEvent &event,
	unsigned short strip,
	unsigned short flag
) {
	for (int i = 0; i < event.back_hit; ++i) {
		if (event.back_strip[i] == strip) continue;
		if ((flag & (1 << (i+8))) != 0) continue;
		if (abs(event.back_strip[i] - strip) != 1) continue;
		// returns the first adjacent strip, since it has the largest energy
		return int(i);
	}
	return -1;
}


/// @brief Fill merge event from normalize result event, version 2
/// @param[in] event normalize result event
/// @param[in] energy_diff energy difference tolerance
/// @param[out] merge merge event to fill
/// @returns merge hit
int Dssd::FillMergeEvent(
	const DssdFundamentalEvent &event,
	double energy_diff,
	DssdMergeEvent &merge
) {
	// initialize merge event
	merge.hit = 0;
	unsigned short used_flag = 0;

	// for convenience
	const int &fhit = event.front_hit;
	const int &bhit = event.back_hit;
	const unsigned short *fs = event.front_strip;
	const unsigned short *bs = event.back_strip;
	const double *fe = event.front_energy;
	const double *be = event.back_energy;
	const double *ft = event.front_time;
	// const double *bt = event.back_time;

	// found new event in the loop, initialize to true to start loop
	bool found = true;
	// loop to find merge event
	while (found && merge.hit < 4) {
		found = false;
		// search f1b1 event
		for (int i = 0; i < fhit; ++i) {
			// jump if it's used event
			if ((used_flag & (1 << i)) != 0) continue;
			// check f1b1 event
			for (int j = 0; j < bhit; ++j) {
				// jump if it's used event
				if ((used_flag & (0x100 << j)) != 0) continue;
				// found, check energy now
				double de = fe[i] - be[j];
				// jump if energy is out of range
				if (fabs(de) < energy_diff) {
					int num = merge.hit;
					merge.flag[num] = (0x1 << i) | (0x100 << j);
					merge.merge_tag[num] = 0;
					merge.energy[num] = fe[i];
					merge.time[num] = ft[i];
					merge.x[num] = fs[i];
					merge.y[num] = bs[j];
					merge.z[num] = 0.0;
					merge.hit += 1;
					used_flag |= (0x1 << i) | (0x100 << j);
					found = true;
					break;
				}
			}
		}
		if (found) continue;
		// try to find f2b2 or f2b1 event
		for (int i = 0; i < fhit; ++i) {
			// jump if it's used event
			if ((used_flag & (1 << i)) != 0) continue;
			// search adjacent strip
			int fi = SearchFrontAdjacentStrips(event, fs[i], used_flag);
			// not found, continue
			if (fi < 0) continue;
			// total energy of adjacent strips in front side
			double total_front_energy = fe[i] + fe[fi];
			// found adjacent strips in front side, check back side now
			for (int j = 0; j < bhit; ++j) {
				// jump if it's used event
				if ((used_flag & (0x100 << j)) != 0) continue;
				// search adjacent strip
				int bi = SearchBackAdjacentStrips(event, bs[j], used_flag);
				// not found, continue
				if (bi < 0) continue;
				// found, check energy now
				double total_back_energy = be[j] + be[bi];
				double de = total_front_energy - total_back_energy;
				// jump if energy is out of range
				if (fabs(de) > energy_diff) continue;
				// check if this is a two particle event
				double de1 = fe[i] - be[j];
				double de2 = fe[fi] - be[bi];
				if (
					fabs(de1) < energy_diff
					&& fabs(de2) < energy_diff
				) {
					// fill f2b2-2 event
					int num = merge.hit;
					merge.flag[num] = (0x1 << i) | (0x100 << j);
					merge.merge_tag[num] = 0;
					merge.energy[num] = fe[i];
					merge.time[num] = ft[i];
					merge.x[num] = fs[i];
					merge.y[num] = bs[i];
					merge.z[num] = 0.0;
					merge.flag[num+1] = (0x1 << fi) | (0x100 << bi);
					merge.merge_tag[num+1] = 0;
					merge.energy[num+1] = fe[fi];
					merge.time[num+1] = ft[fi];
					merge.x[num+1] = fs[fi];
					merge.y[num+1] = bs[bi];
					merge.z[num+1] = 0.0;
					merge.hit += 2;
					used_flag |= (0x1 << i) | (0x100 << j);
					used_flag |= (0x1 << fi) | (0x100 << bi);
					found = true;
					break;
				} else {
					int num = merge.hit;
					merge.flag[num] = (0x1 << i) | (0x1 << fi);
					merge.flag[num] |= (0x100 << j)  | (0x100 << bi);
					merge.merge_tag[num] = 3;
					merge.energy[num] = total_front_energy;
					merge.time[num] = ft[i];
					merge.x[num] =
						fs[i] + fs[fi]/(fe[i]+fe[fi])*(fs[fi]-fs[i]);
					merge.y[num] =
						bs[j] + bs[bi]/(be[j]+be[bi])*(bs[bi]-bs[j]);
					merge.z[num] = 0.0;
					merge.hit += 1;
					used_flag |= (0x1 << i) | (0x100 << j);
					used_flag |= (0x1 << fi) | (0x100 << bi);
					found = true;
					break;
				}
			}
			// found f2b2 event
			if (found) break;
			// adjacent strips in back side not found, check f2b1 event
			for (int j = 0; j < bhit; ++j) {
				// jump if it's used event
				if ((used_flag & (0x100 << j)) != 0) continue;
				// energy difference of sides
				double de = total_front_energy - be[j];
				if (fabs(de) < energy_diff) {
					int num = merge.hit;
					merge.flag[num] = (0x1 << i) | (0x1 << fi) | (0x100 << j);
					merge.merge_tag[num] = 2;
					merge.energy[num] = total_front_energy;
					merge.time[num] = ft[i];
					merge.x[num] =
						fs[i] + fs[fi]/(fe[i]+fe[fi])*(fs[fi]-fs[i]);
					merge.y[num] = bs[j];
					merge.z[num] = 0.0;
					merge.hit += 1;
					used_flag |= (0x1 << i) | (0x1 << fi) | (0x100 << j);
					found = true;
					break;
				}
			}
			// found f2b1 event
			if (found) break;
		}
		// found f2b2 or f2b1 event
		if (found) continue;
		// check f1b2 or f1b1 event now
		for (int i = 0; i < fhit; ++i) {
			// jump if it's used event
			if ((used_flag & (1 << i)) != 0) continue;
			// check back side now
			for (int j = 0; j < bhit; ++j) {
				// jump if it's used event
				if ((used_flag & (1 << (j+8))) != 0) continue;
				// search adjacent strip
				int bi = SearchBackAdjacentStrips(event, bs[j], used_flag);
				// not found, continue
				if (bi < 0) continue;
				// found, check energy now
				double de = fe[i] - be[j] - be[bi];
				// jump if energy is out of range
				if (fabs(de) < energy_diff) {
					int num = merge.hit;
					merge.flag[num] = (0x1 << i) | (0x100 << j) | (0x100 << bi);
					merge.merge_tag[num] = 1;
					merge.energy[num] = fe[i];
					merge.time[num] = ft[i];
					merge.x[num] = fs[i];
					merge.y[num] =
						bs[j] + bs[bi]/(be[j]+be[bi])*(bs[bi]-bs[j]);
					merge.z[num] = 0.0;
					merge.hit += 1;
					used_flag |= (0x1 << i) | (0x100 << j) | (0x100 << bi);
					found = true;
					break;
				}
			}
			// found f1b2 event
			if (found) break;
		}
	}
	if (merge.hit > 4) merge.hit = 4;
	return merge.hit;
}



/// @brief Fill merge event from normalize result event, version 2
/// @param[in] event normalize result event
/// @param[in] energy_diff energy difference tolerance
/// @param[out] merge merge event to fill
/// @returns merge hit
int Dssd::FillMergeEventSupplementary(
	const DssdFundamentalEvent &event,
	double energy_diff,
	DssdMergeEvent &merge
) {
	// initialize merge event
	merge.hit = 0;
	unsigned short used_flag = 0;

	// for convenience
	const int &fhit = event.front_hit;
	const int &bhit = event.back_hit;
	const unsigned short *fs = event.front_strip;
	const unsigned short *bs = event.back_strip;
	const double *fe = event.front_energy;
	const double *be = event.back_energy;
	const double *ft = event.front_time;
	// const double *bt = event.back_time;

	// found new event in the loop, initialize to true to start loop
	bool found = true;
	// loop to find merge event
	while (found && merge.hit < 4) {
		found = false;
		// search f1b1 event
		for (int i = 0; i < fhit; ++i) {
			// jump if it's used event
			if ((used_flag & (1 << i)) != 0) continue;
			// check f1b1 event
			for (int j = 0; j < bhit; ++j) {
				// jump if it's used event
				if ((used_flag & (0x100 << j)) != 0) continue;
				// found, check energy now
				double de = fe[i] - be[j];
				// jump if energy is out of range
				if (fabs(de) < energy_diff) {
					int num = merge.hit;
					merge.flag[num] = (0x1 << i) | (0x100 << j);
					merge.merge_tag[num] = 0;
					merge.energy[num] = fe[i];
					merge.time[num] = ft[i];
					merge.x[num] = fs[i];
					merge.y[num] = bs[j];
					merge.z[num] = 0.0;
					merge.hit += 1;
					used_flag |= (0x1 << i) | (0x100 << j);
					found = true;
					break;
				}
			}
		}
		if (found) continue;
		// search f1b2 bind events
		for (int i = 0; i < fhit; ++i) {
			// jump if it's used event
			if ((used_flag & (1 << i)) != 0) continue;
			// check back side
			for (int j = 0; j < bhit; ++j) {
				// ignore used event
				if ((used_flag & (1 << (j+8))) != 0) continue;
				// search for binding event
				for (int k = j+1; k < bhit; ++k) {
					// ignore used event
					if ((used_flag & (1 << (k+8))) != 0) continue;
					// found, check energy
					double de = fe[i] - be[j] - be[k];
					if (fabs(de) < energy_diff && abs(bs[j]-bs[k]) != 1) {
						int num = merge.hit;
						merge.flag[num] = (0x1 << i) | (0x100 << j);
						merge.merge_tag[num] = 4;
						merge.energy[num] = fe[i] * be[j] / (be[j] + be[k]);
						merge.time[num] = ft[i];
						merge.x[num] = fs[i];
						merge.y[num] = bs[j];
						merge.z[num] = 0.0;
						merge.flag[num+1] = (0x1 << i) | (0x100 << k);
						merge.merge_tag[num+1] = 4;
						merge.energy[num+1] = fe[i] * be[k] / (be[j] + be[k]);
						merge.time[num+1] = ft[i];
						merge.x[num+1] = fs[i];
						merge.y[num+1] = bs[k];
						merge.z[num+1] = 0.0;
						merge.hit += 2;
						used_flag |= (0x1 << i) | (0x100 << j) | (0x100 << k);
						found = true;
						break;
					}
				}
			}
		}
		if (found) continue;
		// search f2b1 bind events
		for (int i = 0; i < bhit; ++i) {
			// jump if it's used event
			if ((used_flag & (0x100 << i)) != 0) continue;
			// check front side now
			for (int j = 0; j < fhit; ++j) {
				// ignore used event
				if ((used_flag & (0x1 << j)) != 0) continue;
				for (int k = j+1; k < fhit; ++k) {
					// ignore used event
					if ((used_flag & (0x1 << k)) != 0) continue;
					// found, check energy
					double de = fe[j] + fe[k] - be[i];
					if (fabs(de) < energy_diff && abs(fs[j]-fs[k]) != 1) {
						int num = merge.hit;
						merge.flag[num] = (0x100 << i) | (0x1 << j);
						merge.merge_tag[num] = 5;
						merge.energy[num] = fe[j];
						merge.time[num] = ft[j];
						merge.x[num] = fs[j];
						merge.y[num] = bs[i];
						merge.z[num] = 0.0;
						merge.flag[num+1] = (0x100 << i) | (0x1 << k);
						merge.merge_tag[num+1] = 5;
						merge.energy[num+1] = fe[k];
						merge.time[num+1] = ft[k];
						merge.x[num+1] = fs[k];
						merge.y[num+1] = bs[i];
						merge.z[num+1] = 0.0;
						merge.hit += 2;
						used_flag |= (0x100 << i) | (0x1 << j) | (0x1 << k);
						found = true;
						break;
					}
				}
			}
		}
	}
	if (merge.hit > 4) merge.hit = 4;
	return merge.hit;
}


int Dssd::Merge(double energy_diff, int supplementary) {
	std::string tag = tag_;
	if (supplementary == 1) {
		if (tag.empty()) tag = "s1";
		else tag += "-s1";
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
	// double ref_time;
	// setup input branches
	fundamental_event.SetupInput(ipt);
	// ipt->SetBranchAddress("ref.time", &ref_time);
	// for convenient
	int &fhit = fundamental_event.front_hit;
	int &bhit = fundamental_event.back_hit;

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
	// setup output branches
	merge_event.SetupOutput(&opt);

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
	ipf.Close();

	statistics.merged =
		statistics.one_hit + statistics.two_hit
		+ statistics.three_hit + four_hit;
	// save and show statistics
	statistics.Write();
	statistics.Print();

	return 0;
}



//-----------------------------------------------------------------------------
//									time
//-----------------------------------------------------------------------------

int Dssd::AnalyzeTime() {
	std::cerr << "Error: Dssd::AnalyzeTime is not implemented yet.\n";
	return -1;
}


int Dssd::ReadNormalizeTimeParameters() {
	// setup file name
	TString file_name;
	file_name.Form(
		"%s%s%s-norm-param-%04u.txt",
		kGenerateDataPath, kTimeDir, name_.c_str(), run_
	);
	// open file
	std::ifstream fin(file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open normalize parameters file "
			<< file_name << " failed.\n";
		return -1;
	}
	// read normalized paramters for front strips
	for (size_t i = 0; i < Strip(0); ++i) {
		fin >> norm_time_params_[0][i];
	}
	// read normalized parameters for back strips
	for (size_t i = 0; i < Strip(1); ++i) {
		fin >> norm_time_params_[1][i];
	}
	// close file
	fin.close();
	return 0;
}


int Dssd::WriteNormalizeTimeParameters() {
	// setup file name
	TString file_name;
	file_name.Form(
		"%s%s%s-norm-param-%04u.txt",
		kGenerateDataPath, kTimeDir, name_.c_str(), run_
	);
	// open file
	std::ofstream fout(file_name.Data());
	if (!fout.good()) {
		std::cerr << "Error: Open normalize parameters file "
			<< file_name << " failed.\n";
		return -1;
	}
	// read normalized paramters for front strips
	for (size_t i = 0; i < Strip(0); ++i) {
		fout << norm_time_params_[0][i] << "\n";
	}
	// read normalized parameters for back strips
	for (size_t i = 0; i < Strip(1); ++i) {
		fout << norm_time_params_[1][i] << "\n";
	}
	// close file
	fout.close();
	return 0;
}


int Dssd::NormalizeTime() {
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
	// input fundamental event
	DssdFundamentalEvent event;
	// input reference time
	double ref_time;
	// setup input branches
	event.SetupInput(ipt);
	ipt->SetBranchAddress("ref.time", &ref_time);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-normalize-%s%04u.root",
		kGenerateDataPath,
		kTimeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// time 1D histogram with CFD
	std::vector<TH1F> hist_time_cfd;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < Strip(side); ++i) {
			if (side == 0 && i >= 32 && i < 48) {
				hist_time_cfd.emplace_back(
					TString::Format("%ctc%ld", "fb"[side], i),
					"time with CFD",
					1000, -200, 800
				);
			} else {
				hist_time_cfd.emplace_back(
					TString::Format("%ctc%ld", "fb"[side], i),
					"time with CFD",
					500, -200, 800
				);
			}
		}
	}
	// time 1D histogram with LE
	std::vector<TH1F> hist_time_le;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < Strip(side); ++i) {
			hist_time_le.emplace_back(
				TString::Format("%ctl%ld", "fb"[side], i),
				"time with LE",
				400, -200, 800
			);
		}
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling time to histogram   0%%");
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
		if (event.front_hit == 1) {
			if ((event.cfd_flag & 0x1) == 0) {
				hist_time_cfd[event.front_strip[0]].Fill(
					event.front_time[0] - ref_time
				);
			} else {
				hist_time_le[event.front_strip[0]].Fill(
					event.front_time[0] - ref_time
				);
			}
		}
		// fill back events
		if (event.back_hit == 1) {
			if ((event.cfd_flag & 0x100) == 0) {
				hist_time_cfd[Strip(0)+event.back_strip[0]].Fill(
					event.back_time[0] - ref_time
				);
			} else {
				hist_time_le[Strip(0)+event.back_strip[0]].Fill(
					event.back_time[0] - ref_time
				);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit front
	for (size_t i = 0; i < Strip(0); ++i) {
		// max value in all bins
		double max_value = 0.0;
		// bin number with max value
		int max_bin = 0;
		// search for max value
		for (int j = 0; j < hist_time_cfd[i].GetNbinsX(); ++j) {
			if (hist_time_cfd[i].GetBinContent(j) > max_value) {
				if (
					i >= 32 &&  i < 48
					&& hist_time_cfd[i].GetBinCenter(j) < 60
				) {
					continue;
				}
				max_value = hist_time_cfd[i].GetBinContent(j);
				max_bin = j;
			}
		}
		// marker to point the max value
		TMarker *marker = new TMarker(
			hist_time_cfd[i].GetBinCenter(max_bin), max_value, 20
		);
		marker->SetMarkerColor(kRed);
		hist_time_cfd[i].GetListOfFunctions()->Add(marker);
		// record paramters
		norm_time_params_[0][i] = hist_time_cfd[i].GetBinCenter(max_bin);
	}
	// fit back
	for (size_t i = Strip(0); i < Strip(0)+Strip(1); ++i) {
		// max value in all bins
		double max_value = 0.0;
		// bin number with max value
		int max_bin = 0;
		// search for max value
		for (int j = 0; j < hist_time_cfd[i].GetNbinsX(); ++j) {
			if (hist_time_cfd[i].GetBinContent(j) > max_value) {
				max_value = hist_time_cfd[i].GetBinContent(j);
				max_bin = j;
			}
		}
		// marker to point the max value
		TMarker *marker = new TMarker(
			hist_time_cfd[i].GetBinCenter(max_bin), max_value, 20
		);
		marker->SetMarkerColor(kRed);
		hist_time_cfd[i].GetListOfFunctions()->Add(marker);
		// record parameters
		norm_time_params_[1][i-Strip(0)]
			= hist_time_cfd[i].GetBinCenter(max_bin);
	}

	// write parameters
	if (WriteNormalizeTimeParameters()) {
		std::cerr << "Error: Write normalize time parameters failed.\n";
		return -1;
	}

	// save histograms
	for (TH1F &hist : hist_time_cfd) {
		hist.Write();
	}
	for (TH1F &hist : hist_time_le) {
		hist.Write();
	}
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}


int Dssd::FilterTimeCurve() {
	std::cerr << "Error: Dssd::FilterTimeCurve is not implemented yet.\n";
	return -1;
}


int Dssd::FitTimeCurve() {
	std::cerr << "Error: Dssd::FitTimeCurve is not implemented yet.\n";
	return -1;
}


bool Dssd::CheckTime(int, size_t, unsigned short, double, double) {
	return true;
}


int Dssd::ReadTimeCuts() {
	// std::cerr << "Error: Dssd::ReadTimeCuts is not implemented yet.\n";
	return 0;
}

}