#include "include/detector/adssd.h"

#include <iostream>

#include <TF1.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TMarker.h>
#include <TMultiGraph.h>

#include "include/event/ta_event.h"

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


int Adssd::Merge(double, int) {
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
			// merge_event.theta[0] = position.Theta();
			merge_event.theta[0] = fs[0];
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


/// @brief function to fit expotional decay tail
/// @param[in] x time
/// @param[in] par parameters, par0-A, par1-tau
/// @returns A*exp(-t/tau)
///
double ExpDecay(double *x, double *par) {
	double t = x[0];
	double A = par[0];
	double tau = par[1];
	return A * exp(-t / tau);
}

/// @brief function to fit linear rise head
/// @param[in] x time
/// @param[in] par parameters, par0-k, par1-t0
/// @returns 0 if t < t0, k(t-t0) if t >= t0
///
double LinearRise(double *x, double *par) {
	double t = x[0];
	double k = par[0];
	double t0 = par[1];

	if (t < t0) return 0.0;
	else return k * (t - t0);
}

double ExpRiseExpDecay(double *x, double *par) {
	double t = x[0];
	double t0 = par[3];
	if (t < t0) return 0.0;
	double A = par[0];
	double tau = par[1];
	double theta = par[2];
	t -= t0;
	double ret = A * (exp(-t/tau) - exp(-t/theta));
	return ret;
}

double LinearRiseExpDecay(double *x, double *par) {
	double t = x[0];
	double t0 = par[3];
	if (t < t0) return 0.0;
	double A = par[0];
	double tau = par[1];
	double p1 = par[2];
	double t1 = A / p1;

	t -= t0;
	double ret;
	if (t > t1) ret = A * exp(-t/tau);
	else ret = p1 * t + A * exp(-t/tau) - A;
	return ret;
}

double CosineVibrate(double *x, double *par) {
	double t = x[0];
	double t0 = par[4];
	if (t < t0) return 0.0;

	double B = par[0];
	double omega = par[1];
	double phi = par[2];
	double kappa = par[3];

	t -= t0;
	double ret = B * cos(omega*t+phi) * exp(-t/kappa);
	return ret;
}

double LinearRiseExpDecayCosineVibrate(double *x, double *par) {
	double *f1Par, f2Par[5];
	f1Par = par;
	f2Par[0] = par[4];
	f2Par[1] = par[5];
	f2Par[2] = par[6];
	f2Par[3] = par[7];
	f2Par[4] = par[0] / par[2] + par[3];
	return LinearRiseExpDecay(x, f1Par) + CosineVibrate(x, f2Par);
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
	// input taf telescope file name
	TString tele_file_name;
	tele_file_name.Form(
		"%s%staf%d-telescope-ta-%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		name_[4] - '0',
		run_
	);
	// add friend
	ipt->AddFriend("tele=tree", tele_file_name);
	// trace points
	unsigned short points;
	// trace data
	unsigned short raw_trace[1024];
	// input tele event
	TaEvent tele_event;
	// setup input branches
	ipt->SetBranchAddress("point", &points);
	ipt->SetBranchAddress("trace", raw_trace);
	tele_event.SetupInput(ipt, "tele.");

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
		1000, 0, 1000, 2000, 0, 2
	);
	// relative value at point
	TH1F hist_relative(
		"hr", "relative value at one point", 1000, 0, 1
	);
	// trace in the energy range
	TMultiGraph energy_trace[20];
	// filled graph number
	size_t graph_num = 0;
	// output tree
	TTree opt("tree", "rise time");
	// output data
	// trace rise time
	double rise_time;
	// trace magnitude
	double magnitude;
	// trace fitted a
	double fit_a;
	// trace fitted tau
	double fit_tau;
	// trace fitted t0
	double fit_t0;
	// trace fitted k
	double fit_k;
	// trace fitted b
	double fit_b;
	// trace fitted omega
	double fit_omega;
	// trace fitted phi
	double fit_phi;
	// trace fitted kappa
	double fit_kappa;
	// setup output branches
	opt.Branch("rise_time", &rise_time, "rt/D");
	opt.Branch("magnitude", &magnitude, "mag/D");
	opt.Branch("fit_a", &fit_a, "fa/D");
	opt.Branch("fit_tau", &fit_tau, "ftau/D");
	opt.Branch("fit_k", &fit_k, "fk/D");
	opt.Branch("fit_t0", &fit_t0, "ft0/D");
	opt.Branch("fit_b", &fit_b, "fb/D");
	opt.Branch("fit_omega", &fit_omega, "fomega/D");
	opt.Branch("fit_phi", &fit_phi, "fphi/D");
	opt.Branch("fit_kappa", &fit_kappa, "fkappa/D");

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

		// fill all trace to histogram
		for (unsigned short i = 0; i < points; ++i) {
			hist_trace.Fill(i, trace[i]/trace[peak_point]);
		}
		// fill relative value at point 500
		hist_relative.Fill(trace[900]/trace[peak_point]);
		// fill to energy range trace histogram
		constexpr double energy_range = 0.3;
		if (tele_event.num == 1 && tele_event.flag[0] == 0x1) {
			for (int i = 0; i < 20; ++i) {
				if (
					tele_event.energy[0][0] < double(i)+energy_range
					&& tele_event.energy[0][0] > double(i)-energy_range
				) {
					TGraph *g = new TGraph;
					for (unsigned short j = 0; j < points; ++j) {
						// energy_trace[i].Fill(j, trace[j]);
						g->AddPoint(j, trace[j]);
					}
					energy_trace[i].Add(g);
					break;
				}
			}
		}

		// graph of trace to fit (and save)
		TGraph g;
		// add graph points
		for (unsigned short i = 0; i < points; ++i) {
			g.AddPoint(i, trace[i]);
		}

		// fit trace
		// 1. fit the exponential decay tail
		TF1 *exp_decay = new TF1(
			TString::Format("tail%lld", entry),
			ExpDecay, 400, 800, 2
		);
		exp_decay->SetParName(0, "A");
		exp_decay->SetParName(1, "tau");
		exp_decay->SetParameter(0, trace[peak_point]*1.6);
		// this initial value of parameter is edited manually
		exp_decay->SetParameter(1, 22000);
		exp_decay->SetLineColor(kGreen);
		g.Fit(exp_decay, "RQ");
		fit_a = exp_decay->GetParameter(0);
		fit_tau = exp_decay->GetParameter(1);

 		// 2. fit the linear rise head
		TF1 *linear_rise = new TF1(
			TString::Format("rise%lld", entry),
			LinearRise, 60, (peak_point-2), 2
		);
		linear_rise->SetParName(0, "k");
		linear_rise->SetParName(1, "t0");
		// this initial value of parameter is edited manually
		linear_rise->SetParameter(0, trace[peak_point]/20.0);
		// this initial value of parameter is edited manually
		linear_rise->SetParameter(1, 75.0);
		linear_rise->SetLineColor(kGreen);
		g.Fit(linear_rise, "RQ");
		fit_k = linear_rise->GetParameter(0);
		fit_t0 = linear_rise->GetParameter(1);

		// exp rise exp decay
		// TF1 *exp_rise_exp_decay = new TF1(
		// 	TString::Format("exp_rise_decay%lld", entry),
		// 	ExpRiseExpDecay, 60, 800, 4
		// );
		// exp_rise_exp_decay->SetParameter(0, fit_a);
		// exp_rise_exp_decay->SetParameter(1, fit_tau);
		// exp_rise_exp_decay->SetParameter(2, 1.0/fit_k);
		// exp_rise_exp_decay->SetParameter(3, fit_t0);
		// exp_rise_exp_decay->SetLineColor(kBlue);
		// exp_rise_exp_decay->SetNpx(5000);
		// g.Fit(exp_rise_exp_decay, "RQ");


		// parameter correction
		fit_k += 1.0/fit_tau;
		fit_a /= exp(fit_t0 / fit_tau);

		// 3. combinse the linear rise and exponential decay
		TF1 *linear_rise_exp_decay = new TF1(
			TString::Format("linear_rise_exp_decay%lld", entry),
			LinearRiseExpDecay, 60, 900, 4
		);
		linear_rise_exp_decay->SetParName(0, "A");
		linear_rise_exp_decay->SetParName(1, "tau");
		linear_rise_exp_decay->SetParName(2, "k");
		linear_rise_exp_decay->SetParName(3, "t0");
		double lred_par[4] = {fit_a, fit_tau, fit_k, fit_t0};
		linear_rise_exp_decay->SetParameters(lred_par);
		linear_rise_exp_decay->SetNpx(5000);
		linear_rise_exp_decay->SetLineColor(kBlue);
		g.Fit(linear_rise_exp_decay, "RQ");
		fit_a = linear_rise_exp_decay->GetParameter(0);
		fit_tau = linear_rise_exp_decay->GetParameter(1);
		fit_k = linear_rise_exp_decay->GetParameter(2);
		fit_t0 = linear_rise_exp_decay->GetParameter(3);

		// 4. add the cosine vibrate
		TF1 *linear_rise_exp_decay_cos_vibrate = new TF1(
			TString::Format("linear_rise_exp_decay_cos_vibrate%lld", entry),
			LinearRiseExpDecayCosineVibrate, 60, 900, 8
		);
		linear_rise_exp_decay_cos_vibrate->SetParName(0, "A");
		linear_rise_exp_decay_cos_vibrate->SetParName(1, "tau");
		linear_rise_exp_decay_cos_vibrate->SetParName(2, "k");
		linear_rise_exp_decay_cos_vibrate->SetParName(3, "t0");
		linear_rise_exp_decay_cos_vibrate->SetParName(4, "B");
		linear_rise_exp_decay_cos_vibrate->SetParName(5, "omega");
		linear_rise_exp_decay_cos_vibrate->SetParName(6, "phi");
		linear_rise_exp_decay_cos_vibrate->SetParName(7, "kappa");
		double lredcv_par[8] = {
			fit_a, fit_tau, fit_k, fit_t0, 100.0, 0.35, -2.0, 20.0
		};
		linear_rise_exp_decay_cos_vibrate->SetParameters(lredcv_par);
		linear_rise_exp_decay_cos_vibrate->SetNpx(5000);
		linear_rise_exp_decay_cos_vibrate->SetLineColor(kRed);
		g.Fit(linear_rise_exp_decay_cos_vibrate, "RQ");
		fit_a = linear_rise_exp_decay_cos_vibrate->GetParameter(0);
		fit_tau = linear_rise_exp_decay_cos_vibrate->GetParameter(1);
		fit_k = linear_rise_exp_decay_cos_vibrate->GetParameter(2);
		fit_t0 = linear_rise_exp_decay_cos_vibrate->GetParameter(3);
		fit_b = linear_rise_exp_decay_cos_vibrate->GetParameter(4);
		fit_omega = linear_rise_exp_decay_cos_vibrate->GetParameter(5);
		fit_phi = linear_rise_exp_decay_cos_vibrate->GetParameter(6);
		fit_kappa = linear_rise_exp_decay_cos_vibrate->GetParameter(7);

		// fill to tree
		opt.Fill();

		// fill first 20 trace to graph for checking
		if (graph_num < 20) {
			// add peak marker
			TMarker *marker_peak = new TMarker(
				peak_point, trace[peak_point], 20
			);
			marker_peak->SetMarkerColor(kRed);
			g.GetListOfFunctions()->Add(marker_peak);
			// add 90% peak point
			TMarker *marker_90 = new TMarker(
				point90, trace[point90], 20
			);
			marker_90->SetMarkerColor(kGreen);
			g.GetListOfFunctions()->Add(marker_90);
			// add 10% peak point
			TMarker *marker_10 = new TMarker(
				point10, trace[point10], 20
			);
			marker_10->SetMarkerColor(kGreen);
			g.GetListOfFunctions()->Add(marker_10);

			// add fitted functions
			g.GetListOfFunctions()->Add(exp_decay);
			g.GetListOfFunctions()->Add(linear_rise);
			// g.GetListOfFunctions()->Add(exp_rise_exp_decay);
			g.GetListOfFunctions()->Add(linear_rise_exp_decay);
			g.GetListOfFunctions()->Add(linear_rise_exp_decay_cos_vibrate);

			// save graph
			g.Write(TString::Format("g%ld", graph_num));
			++graph_num;
			// if (graph_num >= 20) break;
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// save graphs
	hist_trace.Write();
	hist_relative.Write();
	for (int i = 0; i < 20; ++i) {
		energy_trace[i].Write(TString::Format("ge%d", i));
	}
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}

}		// namespace ribll