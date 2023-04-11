#include "include/detector/t0d1.h"

namespace ribll {

T0d1::T0d1(unsigned int run, const std::string &tag)
: Dssd(run, "t0d1", tag) {
}


int T0d1::NormalizeSides(TChain *chain, bool iteration) {
	if (SideNormalize(chain, 0, 26, iteration)) {
		std::cerr << "Error: Normalize first side failed.\n";
		return -1;
	}
	if (SideNormalize(chain, 1, 35, iteration)) {
		std::cerr << "Error: Normalize second side failed.\n";
		return -1;
	}
	return 0;
}


bool T0d1::NormEnergyCheck(
	size_t,
	const DssdFundamentalEvent&
) const {
	return true;
}


int T0d1::Merge(double) {
	// // input file name
	// TString fundamental_file_name;
	// fundamental_file_name.Form(
	// 	"%s%s%s-fundamental-%s%04u.root",
	// 	kGenerateDataPath,
	// 	kFundamentalDir,
	// 	name_.c_str(),
	// 	tag_.empty() ? "" : (tag_+"-").c_str(),
	// 	run_
	// );
	// // input file
	// TFile *ipf = new TFile(fundamental_file_name, "read");
	// // input tree
	// TTree *ipt = (TTree*)ipf->Get("tree");
	// if (!ipt) {
	// 	std::cerr << "Error: Get tree from "
	// 		<< fundamental_file_name << " failed.\n";
	// }
	// // input event
	// DssdFundamentalEvent fundamental_event;
	// // setup input branches
	// fundamental_event.SetupInput(ipt);
	// // for convenient
	// unsigned short &fhit = fundamental_event.front_hit;
	// unsigned short &bhit = fundamental_event.back_hit;
	// unsigned short *fs = fundamental_event.front_strip;
	// unsigned short *bs = fundamental_event.back_strip;
	// double *fe = fundamental_event.front_energy;
	// double *be = fundamental_event.back_energy;

	// // output file name
	// TString merge_file_name;
	// merge_file_name.Form(
	// 	"%s%s%s-merge-%s%04u.root",
	// 	kGenerateDataPath,
	// 	kMergeDir,
	// 	name_.c_str(),
	// 	tag_.empty() ? "" : (tag_+"-").c_str(),
	// 	run_
	// );
	// // output file
	// TFile *opf = new TFile(merge_file_name, "recreate");
	// // relatetive difference of front and back side energy
	// TH1F *hrd = new TH1F(
	// 	"hrd", "relateive difference of front and back side energy",
	// 	1000, 0.0, 1.0
	// );
	// // output tree
	// TTree *opt = new TTree("tree", "tree of merged events");
	// // output event
	// DssdMergeEvent merge_event;
	// // setup output branches
	// merge_event.SetupOutput(opt);

	// // read normalized parameters
	// if (ReadNormalizeParameters()) {
	// 	std::cerr << "Error: Read normalize parameters failed.\n";
	// 	return -1;
	// }

	// // total number of entries
	// long long entries = ipt->GetEntries();
	// // 1/100 of entries
	// long long entry100 = entries / 100;
	// // show start
	// printf("Writing merged events   0%%");
	// fflush(stdout);
	// // loop over events
	// for (long long entry = 0; entry < entries; ++entry) {
	// 	// show process
	// 	if (entry % entry100 == 0) {
	// 		printf("\b\b\b\b%3lld%%", entry / entry100);
	// 		fflush(stdout);
	// 	}

	// 	ipt->GetEntry(entry);
	// 	// initialize
	// 	merge_event.hit = 0;

	// 	if (fhit == 1 && bhit == 1) {
	// 		if (fe[0] < 1e4 && be[0] < 1e4) {
	// 			fe[0] = NormEnergy(0, fs[0], fe[0]);
	// 			be[0] = NormEnergy(1, bs[1], be[0]);
	// 			double diff = RelativeDifference(fe[0], be[0]);
	// 			hrd->Fill(diff);
	// 			if (diff < energy_diff) {
	// 				merge_event.energy[0] = fe[0];
	// 				auto position = CalculatePosition(
	// 					fundamental_event.front_strip[0],
	// 					fundamental_event.back_strip[0]
	// 				);
	// 				merge_event.radius[0] = position.R();
	// 				merge_event.phi[0] = position.Phi();
	// 				merge_event.theta[0] = position.Theta();
	// 				merge_event.hit = 1;
	// 			}
	// 		}
	// 	}

	// 	opt->Fill();
	// }
	// // show finish
	// printf("\b\b\b\b100%%\n");

	// // save and close files
	// hrd->Write();
	// opt->Write();
	// opf->Close();
	// ipf->Close();

	return 0;
}

// bool T0D1::NormalizeFrontEnergyCheck(
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
// 	)
// 		return false;


// 	if (!iteration)
// 	{
// 		if (correlation.front_strip[0] < 16)
// 		{
// 		}
// 		else if (correlation.front_strip[0] <= 32)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 18000
// 					&& correlation.front_energy[0] < 30000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 30000
// 					&& correlation.back_energy[0] < 44000
// 				)
// 			)
// 			{
// 				return false;
// 			}

// 		}
// 		else if (correlation.front_strip[0] < 48)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 23000
// 					&& correlation.front_energy[0] < 27000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 34000
// 					&& correlation.back_energy[0] < 41000
// 				)
// 			)
// 			{
// 				return false;
// 			}

// 			// if (
// 			// 	abs(correlation.front_energy[0] - norm_back_energy)
// 			// 	> 0.6 * norm_back_energy
// 			// ) return false;

// 		}
// 		else
// 		{

// 		}

// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs
// 			(
// 				NormalizeEnergy
// 				(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- correlation.back_energy[0]
// 			)
// 			> 1000
// 		) return false;

// 		if (correlation.front_strip[0] < 16)
// 		{


// 		}
// 		else if (correlation.front_strip[0] <= 32)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 23000
// 					&& correlation.front_energy[0] < 26000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 34000
// 					&& correlation.back_energy[0] < 40000
// 				)
// 			)
// 			{
// 				return false;
// 			}


// 		}
// 		else if (correlation.front_strip[0] < 48)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 17000
// 					&& correlation.front_energy[0] < 35000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 25000
// 					&& correlation.back_energy[0] < 49000
// 				)
// 			)
// 			{
// 				return false;
// 			}

// 		}
// 		else
// 		{
// 		}
// 	}

// 	return true;
// }



// bool T0D1::NormalizeBackEnergyCheck(const CorrelatedEvent &correlation, bool iteration) {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
// 	)
// 		return false;

// 	if (!iteration)
// 	{
// 		// if (correlation.back_strip[0] < 16) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 22000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 34000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 32) {
// 		// 	if (
// 		// 		abs(correlation.front_energy[0] - correlation.back_energy[0])
// 		// 		> 0.3 * correlation.back_energy[0]
// 		// 	) return false;

// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 24000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 38000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 48) {
// 		// 	if (
// 		// 		abs(correlation.front_energy[0] - correlation.back_energy[0])
// 		// 		> 0.5 * correlation.back_energy[0]
// 		// 	) return false;

// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 19000 && correlation.back_energy[0] < 28000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 29000 && correlation.back_energy[0] > 42000
// 		// 		)
// 		// 	)) return false;

// 		// } else {

// 		// }

// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs(
// 				NormalizeEnergy(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- NormalizeEnergy(
// 					1,
// 					correlation.back_strip[0],
// 					correlation.back_energy[0]
// 				)
// 			)
// 			> 1000
// 		)
// 			return false;


// 		// if (correlation.back_strip[0] < 16) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 22000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 34000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 32) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 24000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 38000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 48) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 19000 && correlation.back_energy[0] < 28000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 29000 && correlation.back_energy[0] > 42000
// 		// 		)
// 		// 	)) return false;

// 		// } else {

// 		// }

// 	}

// 	return true;
// }



}		// namespace ribll

