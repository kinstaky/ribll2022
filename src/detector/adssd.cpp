#include "include/detector/adssd.h"

#include <iostream>

#include "TH2F.h"
#include "TMarker.h"

namespace ribll {

Adssd::Adssd(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: Dssd(run, name, tag) {
}


ROOT::Math::Polar3DVector Adssd::CalculatePosition(
	unsigned short front_strip,
	unsigned short back_strip
) const {
	double radius = (radius_range_.second-radius_range_.first) / FrontStrip();
	radius  = radius * (front_strip + 0.5) + radius_range_.first;
	double phi = (phi_range_.second - phi_range_.first) / BackStrip();
	phi = phi * (back_strip + 0.5) + phi_range_.first;
	ROOT::Math::Polar3DVector result(radius, 90.0*TMath::DegToRad(), phi);
	result += center_;
	return result;
}


int Adssd::Merge(double) {
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
	// output tree
	TTree opt("tree", "tree of merged events");
	// output event
	AdssdMergeEvent merge_event;
	// setup output branches
	merge_event.SetupOutput(&opt);

	// calibrated parameters
	double cali_params[16][2];
	// calibration parameters file name
	TString param_file_name;
	param_file_name.Form(
		"%s%s%s-alpha-cali-param.txt",
		kGenerateDataPath,
		kCalibrationDir,
		name_.c_str()
	);
	// parameters file
	std::ifstream fin(param_file_name);
	if (!fin.good()) {
		std::cerr << "Error: Open file "
			<< param_file_name << " failed.\n";
		return -1;
	}
	// read parameters from file
	for (size_t i = 0; i < 16; ++i) {
		fin >> cali_params[i][0] >> cali_params[i][1];
		double tmp;
		fin >> tmp >> tmp >> tmp;
	}
	// close parameters file
	fin.close();


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

		if (fhit == 1 && bhit == 1) {
			fe[0] = cali_params[fs[0]][0] + cali_params[fs[0]][1] * fe[0];
			merge_event.energy[0] = fe[0];
			auto position = CalculatePosition(fs[0], bs[0]);
			merge_event.radius[0] = position.R();
			merge_event.phi[0] = position.Phi();
			merge_event.theta[0] = position.Theta();
			merge_event.decode_entry[0] =
				fundamental_event.front_decode_entry[0];
			merge_event.hit = 1;
		}

		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save and close files
	opt.Write();
	opf.Close();
	ipf.Close();

	return 0;
}


int Adssd::AnalyzeTrace() {
	// input trace file name
	TString trace_file_name;
	trace_file_name.Form(
		"%s%s%s-trace-ta-%04u.root",
		kGenerateDataPath,
		kTraceDir,
		name_.c_str(),
		run_
	);
	// input trace file
	TFile ipf(trace_file_name, "read");
	// input trace tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< trace_file_name << " failed.\n";
		return -1;
	}
	// trace points
	unsigned short points;
	// trace data
	unsigned short raw_trace[1024];
	// setup input branches
	ipt->SetBranchAddress("point", &points);
	ipt->SetBranchAddress("trace", raw_trace);

	// trace in floating numbers
	double trace[1024];

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-rise-ta-%04u.root",
		kGenerateDataPath,
		kTraceDir,
		name_.c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// all trace in histogram
	TH2F hist_trace(
		"ht", "histogram of all trace",
		1000, 0, 1000, 1000, -1000, 9000
	);
	// graph for checking trace
	TGraph graph[10];
	// filled graph number
	size_t graph_num = 0;
	// output tree
	TTree opt("tree", "rise time");
	// output data
	double rise_time;
	double magnitude;
	// setup output branches
	opt.Branch("rise_time", &rise_time, "rt/D");
	opt.Branch("magnitude", &magnitude, "mag/D");

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100;
	// show start
	printf("Analyzing trace   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get trace
		ipt->GetEntry(entry);

		if (points == 0) {
			rise_time = -1.0;
			opt.Fill();
			continue;
		}

		// smooth trace
		constexpr unsigned short smooth = 20;
		// sum of several points
		double sum = 0.0;
		for (unsigned short i = 0; i < smooth; ++i) {
			sum += raw_trace[i];
		}
		for (unsigned short i = 0; i < points-smooth; ++i) {
			double tmp = sum;
			sum -= raw_trace[i];
			sum += raw_trace[i+smooth];
			trace[i] = tmp / double(smooth);
		}
		trace[points-smooth] = sum / double(smooth);

		// correct points
		points -= smooth - 1;

		// find baseline
		// base line length
		constexpr unsigned short baseline_length = 40;
		double baseline = 0.0;
		for (unsigned short i = 0; i < baseline_length; ++i) {
			baseline += trace[i];
		}
		baseline /= double(baseline_length);
		for (unsigned short i = 0; i < points; ++i) {
			trace[i] -= baseline;
		}

		// search for first peak over threshold
		// threshold
		constexpr double threshold = 200.0;
		// search for first over threshold point
		// over threshold point
		unsigned short over_point = points;
		for (unsigned short i = 0; i < points; ++i) {
			if (trace[i] > threshold) {
				over_point = i;
				break;
			}
		}
		// trace is under the threshold
		if (over_point == points) {
			rise_time = -2.0;
			opt.Fill();
			continue;
		}
		// search for peak
		// peak point
		unsigned short peak_point = points;
		for (unsigned short i = over_point; i < points; ++i) {
			if (trace[i] > trace[i+1]) {
				peak_point = i;
				break;
			}
		}
		// peak not found
		if (peak_point == points) {
			rise_time = -3.0;
			opt.Fill();
			continue;
		}

		// search for 90% and 10% point and calculate the rise time
		// point with 0.9*peak
		unsigned short point90 = points;
		for (unsigned short i = over_point; i < peak_point; ++i) {
			if (trace[i] > 0.9 * trace[peak_point]) {
				point90 = i;
				break;
			}
		}
		// point with 0.1*peak
		unsigned short point10 = points;
		for (unsigned short i = over_point; i > 0; --i) {
			if (trace[i] < 0.1 * trace[peak_point]) {
				point10 = i;
				break;
			}
		}
		// point90 or point10 not found
		if (point90 == points || point10 == points) {
			rise_time = -4.0;
			opt.Fill();
			continue;
		}

		// calculate rise time
		// rise_time = double(point90 - point10);
		double point90_linear = double(point90)
			- (trace[point90] - 0.9*trace[peak_point])
			/ (trace[point90] - trace[point90-1]);
		double point10_linear = double(point10+1)
			- (trace[point10+1] - 0.1*trace[peak_point])
			/ (trace[point10+1] - trace[point10]);
		rise_time = point90_linear - point10_linear;
		magnitude = trace[peak_point];
		// fill to tree
		opt.Fill();

		// fill all trace to histogram
		for (unsigned short i = 0; i < points; ++i) {
			hist_trace.Fill(i, trace[i]);
		}
		// fill first 10 trace to graph for checking
		if (graph_num < 10) {
			for (unsigned short i = 0; i < points; ++i) {
				graph[graph_num].AddPoint(i, trace[i]);
			}
			// add peak marker
			TMarker *marker_peak = new TMarker(
				peak_point, trace[peak_point], 20
			);
			marker_peak->SetMarkerColor(kRed);
			graph[graph_num].GetListOfFunctions()->Add(marker_peak);
			// add 90% peak point
			TMarker *marker_90 = new TMarker(
				point90, trace[point90], 20
			);
			marker_90->SetMarkerColor(kGreen);
			graph[graph_num].GetListOfFunctions()->Add(marker_90);
			// add 10% peak point
			TMarker *marker_10 = new TMarker(
				point10, trace[point10], 20
			);
			marker_10->SetMarkerColor(kGreen);
			graph[graph_num].GetListOfFunctions()->Add(marker_10);

			++graph_num;
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// save graphs
	hist_trace.Write();
	for (size_t i = 0; i < graph_num; ++i) {
		graph[i].Write(TString::Format("g%ld", i));
	}
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}

}		// namespace ribll