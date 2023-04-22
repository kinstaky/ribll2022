#include "include/telescope/taf.h"

#include <fstream>

#include <TChain.h>
#include <TCutG.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <Math/Vector3D.h>

#include "include/event/dssd_event.h"
#include "include/event/csi_event.h"
#include "include/event/ta_event.h"

namespace ribll {

// thickness of tafd
const double tafd_thick[6] = {151.0, 152.0, 156.0, 150.0, 150.0, 150.0};
// center phi angle of tafd, to distinguish CsI index
const double tafd_center_phi[6] = {
	90.0*TMath::DegToRad(), 30.0*TMath::DegToRad(),
	-30.0*TMath::DegToRad(), -90.0*TMath::DegToRad(),
	-150.0*TMath::DegToRad(), 150.0*TMath::DegToRad()
};


Taf::Taf(unsigned int run, unsigned int index, const std::string &tag)
: Telescope(run, "taf"+std::to_string(index), tag)
, index_(index) {
}


int Taf::Track(double) {
	std::vector<TString> file_names;
	// add tafd
	file_names.push_back(
		TString::Format(
			"%s%s%s-merge-%s%04u.root",
			kGenerateDataPath,
			kMergeDir,
			(std::string("tafd")+name_[3]).c_str(),
			tag_.empty() ? "" : (tag_+"-").c_str(),
			run_
		)
	);
	// add tafcsi
	file_names.push_back(
		TString::Format(
			"%s%s%s-fundamental-%s%04u.root",
			kGenerateDataPath,
			kFundamentalDir,
			"tafcsi",
			tag_.empty() ? "" : (tag_+"-").c_str(),
			run_
		)
	);
	// first input file
	TFile *ipf = new TFile(file_names[0], "read");
	// first input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< file_names[0] << " failed.\n";
		return -1;
	}
	// add second tree
	if (!ipt->AddFriend("csi=tree", file_names[1])) {
		std::cerr << "Error: Get tree from "
			<< file_names[1] << " falied.\n";
		return -1;
	}

	// input tafd event
	AdssdMergeEvent tafd;
	// input tafcsi event
	CircularCsiFundamentalEvent tafcsi;
	// setup input branches
	tafd.SetupInput(ipt);
	tafcsi.SetupInput(ipt, "csi.");


	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "telescope");
	// output event
	TaEvent tele;
	// setup output branches
	tele.SetupOutput(&opt);

	long long total = 0;
	long long match = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Tracking %s telescope   0%%", name_.c_str());
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);

		// initialize telescope event
		tele.num = 0;
		for (size_t i = 0; i < 4; ++i) {
			tele.flag[i] = 0;
		}

		if (tafd.hit == 1 && tafd.energy[0] < 40.0) {
			++total;

			// record this tafd layer
			tele.flag[tele.num] |= 0x1;
			tele.energy[tele.num][0] = tafd.energy[0];
			tele.radius[tele.num][0] = tafd.radius[0];
			tele.theta[tele.num][0] = tafd.theta[0];
			tele.phi[tele.num][0] = tafd.phi[0];

			// track next layer
			// whether there two signal in CsI
			bool conflict = false;
			for (unsigned int i = 0; i < 2; ++i) {
				int csi_index = index_ * 2 + i;
				if (tafcsi.time[csi_index] < -9e4) continue;
				double max_phi = csi_index < 10
					? 120 - 30 * csi_index
					: 480 - 30 * csi_index;
				max_phi *= TMath::DegToRad();
				double min_phi = csi_index < 10
					? 90 - 30 * csi_index
					: 450 - 30 * csi_index;
				min_phi *= TMath::DegToRad();
				if (tafd.phi[0] > max_phi || tafd.phi[0] < min_phi) continue;
				// check if conflict
				if (
					(tele.flag[tele.num] & 0x2) != 0
					|| (tele.flag[tele.num] & 0x4) != 0
				) {
					conflict = true;
					break;
				}
				// tafcsi goes through the test
				tele.flag[tele.num] |= (0x2 << i);
				tele.energy[tele.num][1] = tafcsi.energy[csi_index];
			}
			if (!conflict) ++tele.num;
			else {
				tele.flag[tele.num] = 0;
			}
		}

		if (tele.num > 0 && tele.flag[0] > 0x1) ++match;

		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save and close files
	opt.Write();
	opf.Close();
	ipf->Close();

	std::cout << "Match " << name_ << "  "
		<< match << " / " << total
		<< "  " << double(match) / double(total) << "\n";

	return 0;
}


int Taf::ParticleIdentify() {
	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
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
	TaEvent tele;
	// setup output branches
	tele.SetupInput(ipt);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-particle-type-%s%04u.root",
		kGenerateDataPath,
		kParticleIdentifyDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "particle type");
	// output data
	ParticleTypeEvent type_event;
	// setup output branches
	type_event.SetupOutput(&opt);

	// read cut from files
	std::vector<ParticleCut> cuts[4];
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			std::string cut_tag;
			cut_tag = i == 0 ? "l8" : "g8";
			cut_tag += j == 0 ? "-a" : "-b";
			cuts[i*2+j].push_back({1, 1, ReadCut(cut_tag.c_str(), "1H")});
			cuts[i*2+j].push_back({1, 2, ReadCut(cut_tag.c_str(), "2H")});
			cuts[i*2+j].push_back({1, 3, ReadCut(cut_tag.c_str(), "3H")});
			cuts[i*2+j].push_back({2, 4, ReadCut(cut_tag.c_str(), "4He")});
		}
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Identifying %s particle   0%%", name_.c_str());
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		type_event.num = tele.num;
		// throw more than one particle or zero particle events
		if (type_event.num != 1) {
			type_event.num = 0;
			opt.Fill();
			continue;
		}
		// initialize
		type_event.mass[0] = 0;
		type_event.charge[0] = 0;
		type_event.layer[0] = 0;
		// only particle hits on TAFD and CsI
		// could be identified by dE-E method
		if (tele.flag[0] == 0x3 || tele.flag[0] == 0x5) {
			// delta energy lost in TAFD
			double &de = tele.energy[0][0];
			// energy lost in CsI
			double &e = tele.energy[0][1];
			// cut index to choose series of cuts
			size_t cut_index = 0;
			// add 2 if theta is over 0.8
			cut_index += tele.theta[0][0] > 0.67 ? 2 : 0;
			// add 1 if hit the CsI-B
			cut_index += tele.flag[0] == 0x5 ? 1 : 0;
			// loop the cuts to identify particle
			for (const auto &cut : cuts[cut_index]) {
				if (cut.cut->IsInside(e, de)) {
					type_event.charge[0] = cut.charge;
					type_event.mass[0] = cut.mass;
					type_event.layer[0] = 1;
					break;
				}
			}
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


int Taf::Calibrate(unsigned int length) {
	// if (AlphaCalibrate()) return -1;
	if (CsiCalibrate(length)) return -1;
	return 0;
}


// particle types
const size_t particle_types = 4;
// particle names
const char* const particle_names[particle_types] = {
	"1H", "2H", "3H", "4He"
};
// particle mass
const unsigned int particles_mass[particle_types] = {
	1, 2, 3, 4
};
// particle charge
const unsigned int particles_charge[particle_types] = {
	1, 1, 1, 2
};


// /// @brief L = pedo + a0 * E ^ ((a1+A) / (a2+A))
// /// @param[in] x energy (MeV)
// /// @param[in] par par parameters, par[0]=pedo, par[1]=a0, par[2]=a1,
// ///		par[3]=a2-a1>0, par[4]=A, par[5]=offset
// /// @returns channel
// double LinearExp(double *x, double *par) {
// 	double e = x[0] - par[5];
// 	return par[0] + par[1] * TMath::Power(
// 		e,
// 		((par[2]+par[4]) / (par[2]+par[3]+par[4]))
// 	);
// }

// /// @brief E = ((L - pedo) / a0) ^ ((a2+A) / (a1+A))
// /// @param[in] x channel
// /// @param[in] par parameters, par[0]=pedo, par[1]=a0, par[2]=a1,
// ///		par[3]=a2-a1>0, par[4]=A, par[5]=offset
// /// @returns energy (MeV)
// double LinearExpInverse(double *x, double *par) {
// 	return TMath::Power(
// 		(x[0] - par[5] - par[0]) / par[1],
// 		(par[2]+par[3]+par[4]) / (par[2]+par[4])
// 	);
// }

// /// @brief L = pedo + a0 * (E - a1 * A * Z^2 * ln(1 + E/(a1*A*Z^2)))
// /// @param[in] x Energy(MeV)
// /// @param[in] par parameters, par[0]=pedo, par[1]=a0, par[2]=a1,
// ///		par[3]=A, par[4]=offset
// /// @returns channel
// double LinearLog(double *x, double *par) {
// 	double t = par[2] * par[3] * 4.0;
// 	double e = x[0] - par[4];
// 	return par[0] + par[1] * (e - t * TMath::Log(e/t + 1.0));
// }


/// @brief fitting function of H iostopes in CsI, L = a0*E^a1+a2
/// @param[in] x energy (MeV)
/// @param[in] par parameters, par[0]-a0, par[1]-a1, par[2]-a2
/// @returns channel
///
double HFit(double *x, double *par) {
  return par[0] * pow(x[0], par[1]) + par[2];
}


/// @brief fitting function for He isotopes in CsI,
///		L = a0*(E-a1*A*Z^2*ln(1+E/(a1*A*Z^2)))+a2
/// @param[in] x energy (MeV)
/// @param[in] par parameters, par[0]=a0, par[1]=a1, par[2]=a2
/// @returns channel
///
double HeFit(double *x, double *par) {
	// a1 * A * Z^2, A = 4 and Z = 2
	double t = par[1] * 16.0;
	return par[0] * (x[0] - t * TMath::Log((x[0] + t)/ t)) + par[2];
}


int Taf::CsiCalibrate(unsigned int length) {
	// taf telescope event chain
	TChain ipt("taf", "chain of taf events");
	// particle type event chain
	TChain type_chain("type", "chain of particle type events");
	for (unsigned int i = run_; i < run_+length; ++i) {
		if (i == 628) continue;
		ipt.AddFile(TString::Format(
			"%s%s%s-telescope-%s%04u.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			name_.c_str(),
			tag_.empty() ? "" : (tag_+"-").c_str(),
			run_
		));
		type_chain.AddFile(TString::Format(
			"%s%s%s-particle-type-%s%04u.root/tree",
			kGenerateDataPath,
			kParticleIdentifyDir,
			name_.c_str(),
			tag_.empty() ? "" : (tag_+"-").c_str(),
			run_
		));
	}
	ipt.AddFriend(&type_chain);
	// input TA telescope event
	TaEvent ta_event;
	// input particle type event
	ParticleTypeEvent type_event;
	// setup output branches
	ta_event.SetupInput(&ipt);
	type_event.SetupInput(&ipt, "type.");

	// output root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-calibration-%04u-%04u.root",
		kGenerateDataPath,
		kCalibrationDir,
		(std::string("taf") + name_[3] + "csi").c_str(),
		run_,
		run_+length-1
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output calibration graphs, the first index is for CsI index
	// the second index is for phi index
	// the third index is for particle type
	TGraph gcali[2][2][particle_types];
	// pid graph with calculated CsI energy for checking
	TH2F pid("pid", "#DeltaE-E", 1000, 0, 200, 1000, 0, 15);

	// CsI energy calculator
	std::vector<elc::CsiEnergyCalculator> csi_calculators;
	for (size_t i = 0; i < particle_types; ++i) {
		csi_calculators.emplace_back(particle_names[i]);
	}

	// total number of entries
	long long entries = ipt.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling CsI calibration graph   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get events
		ipt.GetEntry(entry);
		if (ta_event.num != 1) continue;
		if (ta_event.flag[0] != 0x3 && ta_event.flag[0] != 0x5) continue;
		if (type_event.mass[0] < 1 || type_event.mass[0] > 4) continue;

		// calculate the first index in graph array case by phi angle
		size_t csi_index = ta_event.flag[0] == 0x3 ? 0 : 1;
		// calculate the second index in graph array case by theta angle
		size_t theta_index = ta_event.theta[0][0] < 0.67 ? 0 : 1;
		// calculate the third index with a trick using mass number
		size_t particle_index = type_event.mass[0] - 1;
		// CsI channel number
		double csi_channel = ta_event.energy[0][1];
		// calculated csi energy,
		// the trick using mass number - 1 as index is used
		double csi_energy = csi_calculators[particle_index].Energy(
			ta_event.theta[0][0], ta_event.energy[0][0], tafd_thick[index_]
		);
		gcali[csi_index][theta_index][particle_index].AddPoint(
			csi_energy, csi_channel
		);
		pid.Fill(csi_energy, ta_event.energy[0][0]);

		// // fill the appropriate graph, the first index is chosen by phi
		// // the second index is chosen ny theta and charge number
		// // 25000 offset in channel to disinguish different isotopes
		// // charge[0]-1 is the isotope, 0 for H, 1 for He
		// // i - isotope_index[charge[0]-1] is offset in isotope, this
		// // disinguishes between isotopes
		// gcali[csi_index][theta_index][charge[0]-1]->AddPoint(
		// 	csi_energy + 200 * (i - isotope_index[charge[0]-1]),
		// 	csi_channel
		// );
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fitting
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < 2; ++j) {
			// fit 1H (p)
			TF1 fitp (
				TString::Format("fp%c8%c", "lg"[j], "ab"[i]),
				HFit, 0, 600, 3
			);
			// set initial values
			fitp.SetParameter(0, 200.0);
			fitp.SetParameter(1, 0.9);
			fitp.FixParameter(2, 0.0);
			// fit
			gcali[i][j][0].Fit(&fitp, "QR+");
			// record parameters
			fitp.GetParameters(csi_calibrate_params_[i*2+j][0]);

			// fit 2H (d)
			TF1 fitd (
				TString::Format("fd%c8%c", "lg"[j], "ab"[i]),
				HFit, 0, 600, 3
			);
			// set initial values
			fitd.SetParameter(0, 200.0);
			fitd.SetParameter(1, 0.9);
			fitd.FixParameter(2, 0.0);
			// fit
			gcali[i][j][1].Fit(&fitd, "QR+");
			// record parameters
			fitd.GetParameters(csi_calibrate_params_[i*2+j][1]);

			// fit 3H (t)
			TF1 fitt (
				TString::Format("ft%c8%c", "lg"[j], "ab"[i]),
				HFit, 0, 600, 3
			);
			// set initial values
			fitt.SetParameter(0, 200.0);
			fitt.SetParameter(1, 0.1);
			fitt.FixParameter(2, 0.0);
			// fit
			gcali[i][j][2].Fit(&fitt, "QR+");
			// record parameters
			fitt.GetParameters(csi_calibrate_params_[i*2+j][2]);

			// fit 4He (alpha)
			TF1 fita (
				TString::Format("fa%c8%c",  "lg"[j], "ab"[i]),
				HeFit, 0, 600, 3
			);
			// set initial values
			fita.SetParameter(0, 200.0);
			fita.SetParameter(1, 0.9);
			fita.FixParameter(2, 0.0);
			// fit
			gcali[i][j][3].Fit(&fita, "QR+");
			// record parameters
			fita.GetParameters(csi_calibrate_params_[i*2+j][3]);

			// print result
			std::cout << "-----------" << "ab"[i] << "lg"[j] << "\n"
				<< "1H " << fitp.GetParameter(0) << ", " << fitp.GetParameter(1) << ", " << fitp.GetParameter(2) << "\n"
				<< "2H " << fitd.GetParameter(0) << ", " << fitd.GetParameter(1) << ", " << fitd.GetParameter(2) << "\n"
				<< "3H " << fitt.GetParameter(0) << ", " << fitt.GetParameter(1) << ", " << fitt.GetParameter(2) << "\n"
				<< "4He " << fita.GetParameter(0) << ", " << fita.GetParameter(1) << ", " << fita.GetParameter(2) << "\n";
		}
	}

	// save parameters
	if (WriteCsiCalibrateParameters()) {
		std::cerr << "Error: Write CsI calibrate parameters failed.\n";
	}

	opf.cd();
	// save graphs
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < 2; ++j) {
			for (size_t k = 0; k < particle_types; ++k) {
				gcali[i][j][k].Write(TString::Format(
					"g%c%c8%c", "pdta"[k], "lg"[j], "ab"[i]
				));
			}
		}
	}
	pid.Write("pid");
	// close files
	opf.Close();

	return 0;
}


int Taf::ReadCsiCalibrateParameters() {
	// calibrate parameters file name
	TString file_name;
	file_name.Form(
		"%s%staf%ucsi-cali-param.txt",
		kGenerateDataPath,
		kCalibrationDir,
		index_
	);
	// file input stream
	std::ifstream fin(file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open calibration parameters file "
			<< file_name << " failed.\n";
	}
	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			for (size_t k = 0; k < 3; ++k) {
				fin >> csi_calibrate_params_[i][j][k];
			}
		}
	}
	// close file
	fin.close();
	return 0;
}


int Taf::WriteCsiCalibrateParameters() {
	// calibrate parameters file name
	TString file_name;
	file_name.Form(
		"%s%staf%ucsi-cali-param.txt",
		kGenerateDataPath,
		kCalibrationDir,
		index_
	);
	// file output stream
	std::ofstream fout(file_name.Data());
	if (!fout.good()) {
		std::cerr << "Error: Open calibration parameters file "
			<< file_name << " failed.\n";
	}
	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			for (size_t k = 0; k < 3; ++k) {
				fout << csi_calibrate_params_[i][j][k] << " ";
			}
			fout << "\n";
		}
	}
	// close file
	fout.close();
	return 0;
}


/// @brief calculate He energy in Csi with calibration parameters
/// @param[in] channel CsI energy in channel
/// @param[in] par CsI energy calibration parameters
/// @returns energy of Csi if success, -1e5 for channel too small,
/// 	-2e5 for channel too big, -3e5 for reaching max loop
///
double CsiHeEnergy(double channel, double *par) {
	// the lower bound of CsI energy is 0 MeV
	double lower_bound = 0.0;
	// return an invalid energy if the channel is too small
	if (HeFit(&lower_bound, par) > channel) return -1e5;
	// assume the upper bound of CsI energy is 200 MeV
	double upper_bound = 200.0;
	// return an invalid energy if the channel is too large
	if (HeFit(&upper_bound, par) < channel) return -2e5;

	// binary search for the CsI energy
	// the tolerate accuracy
	double eps = 1e-3;
	// max loop avoiding infinite loops
	int max_loop = 1000;
	for (int loop = 0; loop < max_loop; ++loop) {
		double mid_energy = (lower_bound + upper_bound) / 2.0;
		double mid_channel = HeFit(&mid_energy, par);
		if (fabs(mid_channel - channel) < eps) return mid_energy;
		if (mid_channel < channel) {
			// calculated channel is smaller than request one,
			// increase the lower bound
			lower_bound = mid_energy;
		} else {
			// calculated channel is larger than request one,
			// decrease the upper bound
			upper_bound = mid_energy;
		}
	}
	// reach max loop, returns an invalid energy
	return -3e5;
}


/// @brief calculate CsI energy with calibration parameters
/// @param[in] channel CsI channel recored in DAQ 
/// @param[in] par CsI energy calibration parameters 
/// @param[in] index identfied particle types, 0-1H, 1-2H, 2-3H, 3-4He
/// @returns calibrated CsI energy
///
double CalibrateCsiEnergy(double channel, double *par, size_t index) {
	if (index < 3) {
		return TMath::Power((channel-par[2])/par[0], 1.0/par[1]);
	} else {
		return CsiHeEnergy(channel, par);
	}
}


int Taf::Rebuild() {
	// telescope file name
	TString telescope_file_name;
	telescope_file_name.Form(
		"%s%s%s-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// telescope file
	TFile telescope_file(telescope_file_name, "read");
	// input telescope tree
	TTree *ipt = (TTree*)telescope_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< telescope_file_name << " failed.\n";
		return -1;
	}
	// particle type file name
	TString particle_type_file_name;
	particle_type_file_name.Form(
		"%s%s%s-particle-type-%s%04u.root",
		kGenerateDataPath,
		kParticleIdentifyDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// add friend
	ipt->AddFriend("type=tree", particle_type_file_name);
	// input telescope event
	TaEvent ta_event;
	// input type event
	ParticleTypeEvent type_event;
	// setup input branches
	ta_event.SetupInput(ipt);
	type_event.SetupInput(ipt, "type.");

	// output file name
	TString particle_file_name;
	particle_file_name.Form(
		"%s%s%s-particle-%s%04u.root",
		kGenerateDataPath,
		kParticleDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile particle_file(particle_file_name, "recreate");
	// 2D histogram of TAF pid
	TH2F pid("pid", "#DeltaE-E", 1000, 0, 200, 1000, 0, 20);
	// particle tree
	TTree opt("tree", "taf particles");
	// output particle event
	ParticleEvent particle_event;
	// setup output branches
	particle_event.SetupOutput(&opt);

	// read CsI calibrate parameters
	if (ReadCsiCalibrateParameters()) {
		return -1;
	}

	// totla number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Rebuilding particle   0%%");
	fflush(stdout);
	// loop events
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		// initialize particle event
		particle_event.num = 0;
		if (
			type_event.num == 1
			&& type_event.layer[0] == 1
			&& type_event.charge[0] > 0
			&& type_event.mass[0] > 0
		) {
			particle_event.num = 1;
			// fill particle charge number
			particle_event.charge[0] = type_event.charge[0];
			// fill particle mass number
			particle_event.mass[0] = type_event.mass[0];

			// TAFD energy
			double si_energy = ta_event.energy[0][0];
			// calculate the first index of CsI calibration paramters
			size_t csi_index = ta_event.flag[0] == 0x3 ? 0 : 1;
			// calculate the second index case by theta angle
			size_t theta_index = ta_event.theta[0][0] < 0.67 ? 0 : 1;
			// calculate the third index with a trick using mass number
			// 0-1H, 1-2H, 2-3H, 3-4He
			size_t particle_index = type_event.mass[0] - 1;
			// calculate CsI energy from calibrate parameters
			double csi_energy = CalibrateCsiEnergy(
				ta_event.energy[0][1],
				csi_calibrate_params_[csi_index*2+theta_index][particle_index],
				particle_index
			);
			// calculate the total energy
			particle_event.energy[0] = si_energy + csi_energy;
			// fill to pid
			pid.Fill(csi_energy, si_energy);

			// calcuate the position
			ROOT::Math::Polar3DVector position(
				ta_event.radius[0][0], ta_event.theta[0][0], ta_event.phi[0][0]
			);
			particle_event.x[0] = position.X();
			particle_event.y[0] = position.Y();
			particle_event.z[0] = position.Z();

			// leave empty momentum since the reaction point is unknown
			particle_event.px[0] = 0.0;
			particle_event.py[0] = 0.0;
			particle_event.pz[0] = 0.0;
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histogram
	pid.Write();
	// save particle tree
	particle_file.cd();
	opt.Write();
	// close files
	particle_file.Close();
	telescope_file.Close();

	return 0;
}





}		// namespace ribll