#include "include/detector/dssd.h"

#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>

#include "include/event/dssd_event.h"


namespace ribll {

Dssd::Dssd(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: Detector(run, name, tag) {
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

	// sort events by strip
	std::sort(
		front_events.begin(),
		front_events.end(),
		[](const DssdMapEvent &x, const DssdMapEvent &y) {
			return x.strip < y.strip;
		}
	);
	std::sort(
		back_events.begin(),
		back_events.end(),
		[](const DssdMapEvent &x, const DssdMapEvent &y) {
			return x.strip < y.strip;
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

int Dssd::ReadNormalizeParameters() {
	// setup file name
	TString file_name;
	file_name.Form(
		"%s%s%s.txt",
		kGenerateDataPath, kNormalizeDir, name_.c_str()
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


int Dssd::WriteNormalizeParameters() {
	// setup file name
	TString file_name;
	file_name.Form(
		"%s%s%s.txt",
		kGenerateDataPath, kNormalizeDir, name_.c_str()
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

	return 0;
}


bool Dssd::NormEnergyCheck(size_t, const DssdFundamentalEvent&) {
	return true;
}


int Dssd::SideNormalize(
	TChain *chain,
	size_t side,
	size_t ref_strip,
	bool
) {
	DssdFundamentalEvent event;
	event.SetupInput(chain);

	// energy graph fe:be or be:fe
	TGraph **ge;
	ge = new TGraph*[Strip(side)];
	for (size_t i = 0; i < Strip(side); ++i) {
		ge[i] = new TGraph;
	}

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
		if (event.front_hit != 1 || event.back_hit != 1) continue;

		unsigned short &fs = event.front_strip[0];
		unsigned short &bs = event.back_strip[0];
		double &fe = event.front_energy[0];
		double &be = event.back_energy[0];

		if (side == 0) {
			// jump if not refer strip
			if (bs != ref_strip) continue;
			if (!NormEnergyCheck(side, event)) continue;
			ge[fs]->AddPoint(fe, NormEnergy(1, bs, be));
		} else {
			if (fs != ref_strip) continue;
			if (!NormEnergyCheck(side, event)) continue;
			ge[bs]->AddPoint(be, NormEnergy(0, fs, fe));
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fitting
	std::cout << "side " << side << " normalize parameters.\n";
	for (size_t i = 0; i < Strip(side); ++i) {
		// only fits when over 10 points
		if (ge[i]->GetN() > 10) {
			// fitting function
			TF1 *energy_fit = new TF1("efit", "pol1", 0, 60000);
			// set initial value
			energy_fit->SetParameter(0, 0.0);
			energy_fit->SetParameter(1, 1.0);
			// fit
			ge[i]->Fit(energy_fit, "QR+ ROB=0.9");
			// store the normalized parameters
			norm_params_[side][i][0] = energy_fit->GetParameter(0);
			norm_params_[side][i][1] = energy_fit->GetParameter(1);
		}
		// store the graph
		ge[i]->Write(TString::Format("g%c%ld", (side ? 'b' : 'f'),i));
		// print normalized paramters on screen
		std::cout << i << " " << norm_params_[side][i][0]
			<< ", " << norm_params_[side][i][1] << "\n";
	}

	// residual
	TGraph **res;
	res = new TGraph*[Strip(side)];
	for (size_t i = 0; i < Strip(side); i++) res[i] = new TGraph;
	for (size_t i = 0; i < Strip(side); ++i) {
		int point = ge[i]->GetN();
		double *gex = ge[i]->GetX();
		double *gey = ge[i]->GetY();
		for (int j = 0; j < point; ++j) {
			res[i]->AddPoint(gex[j], NormEnergy(side, i, gex[j])-gey[j]);
		}
		res[i]->Write(TString::Format("res%c%ld", side ? 'b' : 'f', i));
	}

	// free memory
	for (size_t i = 0; i < Strip(side); ++i) {
		delete ge[i];
	}
	delete[] ge;

	return 0;
}


int Dssd::NormalizeSides(TChain*, bool) {
	std::cerr << "Error: Sides Normalize is not implemented yet.\n";
	return -1;
}


int Dssd::WriteNormalizeFiles(unsigned int run, const std::string &tag) {
	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		name_.c_str(),
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// input file
	TFile *ipf = new TFile(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
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
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-normalize-%s%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
		name_.c_str(),
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// output tree
	TTree *opt = new TTree("tree", "normalized energy tree");
	// setup branches
	opt->Branch("front_hit", &event.front_hit, "fhit/s");
	opt->Branch("back_hit", &event.back_hit, "bhit/s");
	opt->Branch("front_energy", event.front_energy, "fe[fhit]/D");
	opt->Branch("back_energy", event.back_energy, "be[bhit]/D");

	// total number of entries
	long long entries = ipt->GetEntries();
	// l/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf(
		"Writing normalized energy for %s in run %u   0%%",
		name_.c_str(), run
	);
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
		for (unsigned short i = 0; i < event.front_hit; ++i) {
			event.front_energy[i] =
				NormEnergy(0, event.front_strip[i], event.front_energy[i]);
		}
		for (unsigned short i = 0; i < event.back_hit; ++i) {
			event.back_energy[i] =
				NormEnergy(1, event.back_strip[i], event.back_energy[i]);
		}
		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// write output tree
	opt->Write();
	// close files
	opf->Close();
	ipf->Close();

	return 0;
}


int Dssd::Normalize(
	unsigned int length,
	bool iteration
) {
	// setup input chain
	TChain *chain = new TChain("tree", "chain");
	for (unsigned int i = 0; i < length; ++i) {
		TString file_name;
		file_name.Form(
			"%s%s%s-fundamental-%s%04u.root",
			kGenerateDataPath,
			kFundamentalDir,
			name_.c_str(),
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			i+run_
		);
		if (!chain->AddFile(file_name)) {
			std::cerr << "Error: Read file "
				<< file_name << "failed.\n";
			return -1;
		}
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
	if (iteration && ReadNormalizeParameters()) {
		std::cerr << "Error: read normalize parameters from file failed.\n";
		return -1;
	}

	// setup normalize record root file
	TString normalize_file_name;
	normalize_file_name.Form(
		"%s%s%s-normalize-fit-%s%04u-%u.root",
		kGenerateDataPath,
		kNormalizeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_,
		length
	);
	// output file
	TFile *opf = new TFile(normalize_file_name, "recreate");

	if (NormalizeSides(chain, iteration)) return -1;

	// close files
	opf->Close();

	// write parameters
	if (WriteNormalizeParameters()) {
		std::cerr << "Error: write normalize paramters to file failed.\n";
		return -1;
	}

	// write trees
	std::cout << "---------------------------------------------------------\n"; 
	for (unsigned int i = run_; i < run_ + length; i++) {
		if (WriteNormalizeFiles(i, tag_)) {
			std::cerr << "Error: Write normalized energy to file in run "
				<< i << " failed.\n";
			return -1;
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
//									merge
//-----------------------------------------------------------------------------

int Dssd::Merge(double) {
	std::cerr << "Error: Dssd::Merge not implemented yet.\n";
	return -1;
}


}