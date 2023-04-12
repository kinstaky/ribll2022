#include "include/telescope/taf.h"

#include <fstream>

#include <TCutG.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>

#include "include/event/dssd_event.h"
#include "include/event/csi_event.h"
#include "include/event/telescope_event.h"

namespace ribll {

// thickness of tafd
const double tafd_thick[6] = {151.0, 152.0, 156.0, 150.0, 150.0, 150.0};
// center phi angle of tafd, to distinguish CsI index
const double tafd_center_phi[6] = {
	90.0*TMath::DegToRad(), 30.0*TMath::DegToRad(),
	-30.0*TMath::DegToRad(), -90.0*TMath::DegToRad(),
	-150.0*TMath::DegToRad(), 150.0*TMath::DegToRad()
};

// run used for alpha calibration
const unsigned int alpha_calibration_run[6] = {816, 0, 0, 0, 0, 0};
// energy of alpha source, in MeV
const double alpha_energy[3] = {5.157, 5.486, 5.805};
// histogram range of alpha source peaks
const double alpha_hist_range[6][2] = {
	{1700, 2200},
	{0, 60000},
	{0, 60000},
	{0, 60000},
	{0, 60000},
	{0, 60000}
};
// fit range of alpha source peaks
const double alpha_fit_range[6][6] = {
	{1850.0, 1880.0, 1970.0, 2000.0, 2090.0, 2110.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};

// particle types
const size_t particle_types = 5;
// particle names
const char* const particle_names[particle_types] = {
	"1H", "2H", "3H", "4He", "3He"
};
// particle mass
const unsigned int particles_mass[particle_types] = {
	1, 2, 3, 4, 3
};
// particle charge
const unsigned int particles_charge[particle_types] = {
	1, 1, 1, 2, 2
};
// isotope types
const size_t isotope_types = 2;
// isotope start index in particles array above
const size_t isotope_index[isotope_types] = {0, 3};
// phi angle types
const size_t theta_types = 2;
// phi angle names
const char* const theta_names[theta_types] = {
	"l08", "g08"
};


/// @brief L = pedo + a0 * E ^ ((a1+A) / (a2+A))
/// @param[in] x energy (MeV)
/// @param[in] par par parameters, par[0]=pedo, par[1]=a0, par[2]=a1,
///		par[3]=a2-a1>0, par[4]=A, par[5]=offset
/// @returns channel
double LinearExp(double *x, double *par) {
	double e = x[0] - par[5];
	return par[0] + par[1] * TMath::Power(
		e,
		((par[2]+par[4]) / (par[2]+par[3]+par[4]))
	);
}

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

/// @brief fitting function of H isotopes in CsI with offset
/// @param[in] x energy (MeV)
/// @param[in] par parameters, par[0]=pedo, par[1]=a0, par[2]=a1,
///		par[3]=a2-a1>0
/// @returns channel
double HFit(double *x, double *par) {
	double fpar[6];
	for(int i = 0; i < 4; ++i) fpar[i] = par[i];
	if (x[0] <= 200) {
		fpar[4] = 1.0;
		fpar[5] = 0.0;
	} else if (x[0] <= 400) {
		fpar[4] = 2.0;
		fpar[5] = 200.0;
	} else {
		fpar[4] = 3.0;
		fpar[5] = 400.0;
	}
  return LinearExp(x, fpar);
}


/// @brief L = pedo + a0 * (E - a1 * A * Z^2 * ln(1 + E/(a1*A*Z^2)))
/// @param[in] x Energy(MeV)
/// @param[in] par parameters, par[0]=pedo, par[1]=a0, par[2]=a1,
///		par[3]=A, par[4]=offset
/// @returns channel
double LinearLog(double *x, double *par) {
	double t = par[2] * par[3] * 4.0;
	double e = x[0] - par[4];
	return par[0] + par[1] * (e - t * TMath::Log(e/t + 1.0));
}

/// @brief fitting function for He isotope in CsI with offset
/// @param[in] x energy (MeV)
/// @param[in] par parameters, par[0]=pedo, par[1]=a0, par[2]=a1
/// @returns channel
double HeFit(double *x, double *par) {
	double fpar[6];
	for(int i = 0; i < 3; ++i) fpar[i] = par[i];
	if (x[0] <= 200.0) {
		fpar[3] = 4.0;
		fpar[4] = 0.0;
	} else {
		fpar[3] = 3.0;
		fpar[4] = 200.0;
	}
	return LinearLog(x, fpar);
}


Taf::Taf(unsigned int run, unsigned int index, const std::string &tag)
: Telescope(run, "taf"+std::to_string(index), tag)
, index_(index) {
}


int Taf::Track() {
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
	TFile *opf = new TFile(output_file_name, "recreate");
	// output tree
	TTree *opt = new TTree("tree", "telescope");
	// output event
	TelescopeEvent tele;
	// setup output branches
	tele.SetupOutput(opt);

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
		tele.particle = 0;
		for (size_t i = 0; i < 4; ++i) {
			tele.flag[i] = 0;
			tele.layer[i] = 0;
		}

		if (tafd.hit == 1) {
			++total;

			// record this tafd layer
			++tele.layer[tele.particle];
			tele.flag[tele.particle] |= 1;
			tele.energy[tele.particle][0] = tafd.energy[0];
			tele.radius[tele.particle][0] = tafd.radius[0];
			tele.theta[tele.particle][0] = tafd.theta[0];
			tele.phi[tele.particle][0] = tafd.phi[0];

			// track next layer
			bool conflict = false;
			for (unsigned int i = 0; i < 2; ++i) {
				unsigned int csi_index = index_ * 2 + i;
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
				if ((tele.flag[tele.particle] & 2) != 0) {
					conflict = true;
					break;
				}
				// tafcsi goes through the test
				++tele.layer[tele.particle];
				tele.flag[tele.particle] |= 2;
				tele.energy[tele.particle][1] = tafcsi.energy[csi_index];
			}
			if (!conflict) ++tele.particle;
			else {
				tele.flag[tele.particle] = 0;
				tele.layer[tele.particle] = 0;
			}
		}

		if (tele.particle > 0) ++match;

		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save and close files
	opt->Write();
	opf->Close();
	ipf->Close();

	std::cout << "Match " << name_ << "  "
		<< match << " / " << total
		<< "  " << double(match) / double(total) << "\n";

	return 0;
}

TCutG* ReadCut(const char *particle, const char *phi) {
	// cut file name
	TString cut_file_name;
	cut_file_name.Form(
		"%s%scut/taf0-%s-%s.txt",
		kGenerateDataPath,
		kParticleIdentifyDir,
		particle,
		phi
	);
	// open cut file to read points
	std::ifstream fin(cut_file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open file "
			<< cut_file_name << " failed.\n";
		return nullptr;
	}
	// result
	TCutG *result = new TCutG;
	// point index
	int point;
	// point positions
	double x, y;
	// loop to read points
	while (fin.good()) {
		fin >> point >> x >> y;
		result->SetPoint(point, x, y);
	}
	// close file
	fin.close();
	return result;
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
	TFile *ipf = new TFile(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input event
	TelescopeEvent tele;
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
	TFile *opf = new TFile(output_file_name, "recreate");
	// output tree
	TTree *opt = new TTree("tree", "particle type");
	// output data
	// particle number
	unsigned short particle_num;
	// A, mass number of particle
	unsigned int mass[4];
	// Z, change number of particle
	unsigned int charge[4];
	// setup output branches
	opt->Branch("particle", &particle_num, "p/s");
	opt->Branch("mass", mass, "a[p]/i");
	opt->Branch("charge", charge, "z[p]/i");

	// read cut from files
	TCutG *pid_cuts[particle_types*theta_types];
	for (size_t i = 0; i < theta_types; ++i) {
		for (size_t j = 0; j < particle_types; ++j) {
			TCutG *cut = ReadCut(particle_names[j], theta_names[i]);
			if (!cut) return -1;
			pid_cuts[i*particle_types+j] = cut;
		}
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Identifying particle   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		particle_num = tele.particle;
		if (tele.particle == 1 && tele.layer[0] == 2) {
			// initialize
			mass[0] = 0;
			// point position in pid graph
			double x = tele.energy[0][1];
			double y = tele.energy[0][0];
			// cut index offset in case theta angle difference
			size_t cut_offset = tele.theta[0][0] < 0.8 ? 0 : particle_types;
			// particle identification
			for (size_t i = 0; i < particle_types; ++i) {
				if (pid_cuts[cut_offset+i]->IsInside(x, y)) {
					mass[0] = particles_mass[i];
					charge[0] = particles_charge[i];
					break;
				}
			}
		}
		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save and close files
	opt->Write();
	opf->Close();
	ipf->Close();

	return 0;
}


int Taf::Calibrate() {
	// if (AlphaCalibrate()) return -1;
	if (CsiCalibrate()) return -1;
	return 0;
}


int Taf::AlphaCalibrate() {
	// calibrate tafd with alpha source
	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-merge-%04u.root",
		kGenerateDataPath,
		kMergeDir,
		(name_.substr(0, 3) + "d" + name_[3]).c_str(),
		alpha_calibration_run[index_]
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
	// input data
	DssdMergeEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-calibration-%04u.root",
		kGenerateDataPath,
		kCalibrationDir,
		(name_.substr(0, 3) + "d" + name_[3]).c_str(),
		alpha_calibration_run[index_]
	);
	// output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// output energy histogram
	TH1F *he = new TH1F(
		"he", "fit alpha source peaks",
		500, alpha_hist_range[index_][0], alpha_hist_range[index_][1]
	);

	// output parameters txt file
	TString param_file_name;
	param_file_name.Form(
		"%s%s%s-alpha-cali-param.txt",
		kGenerateDataPath,
		kCalibrationDir,
		(name_.substr(0, 3) + "d" + name_[3]).c_str()
	);
	std::ofstream fout(param_file_name.Data());
	if (!fout.good()) {
		std::cerr << "Error: Open file "
			<< param_file_name << " failed.\n";
		return -1;
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// fill energy to histogram
	for (long long entry = 0; entry < entries; ++entry) {
		ipt->GetEntry(entry);
		if (event.hit != 1) continue;
		he->Fill(event.energy[0]);
	}

	// fit alpha peaks
	TFitResultPtr fit_result[3];
	for (size_t i = 0; i < 3; ++i) {
		TF1 *fpeak = new TF1(
			TString::Format("fpeak%ld", i),
			"gaus",
			alpha_fit_range[index_][2*i],
			alpha_fit_range[index_][2*i+1]
		);
		fpeak->SetParameter(0, 1000);
		fpeak->SetParameter(1, alpha_fit_range[index_][2*i]+10);
		fit_result[i] = he->Fit(fpeak, "RQS+");
	}


	// calibrate
	TGraph *gcali = new TGraph;
	for (size_t i = 0; i < 3; ++i) {
		gcali->AddPoint(fit_result[i]->Parameter(1), alpha_energy[i]);
	}
	TF1 *fcali = new TF1(
		"fcali", "pol1",
		alpha_hist_range[index_][0], alpha_hist_range[index_][1]
	);
	TFitResultPtr cali_result = gcali->Fit(fcali, "QRS+");

	// show calibrate result
	std::cout << cali_result->Parameter(0) << " "
		<< cali_result->Parameter(1) << "\n";

	// output fit results
	fout << cali_result->Parameter(0) << " "
		<< cali_result->Parameter(1) << " "
		<< cali_result->Chi2() << " "
		<< cali_result->Ndf() << " "
		<< cali_result->Chi2() / cali_result->Ndf() << "\n";
	for (size_t i = 0; i < 3; ++i) {
		fout << fit_result[i]->Parameter(0) << " "
			<< fit_result[i]->Parameter(1) << " "
			<< fit_result[i]->Parameter(2) << " "
			<< fit_result[i]->Chi2() << " "
			<< fit_result[i]->Ndf() << " "
			<< fit_result[i]->Chi2() / fit_result[i]->Ndf() << "\n";
	}
	fout.close();

	// save histogram and close files
	he->Write();
	for (size_t i = 0; i < 3; ++i) {
		fit_result[i]->Write(TString::Format("r%ld", i));
	}
	gcali->Write("gcali");
	cali_result->Write("rcali");
	opf->Close();
	// close input file
	ipf->Close();

	return 0;
}


int Taf::CsiCalibrate() {
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
	TFile *ipf = new TFile(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input event
	TelescopeEvent tele;
	// setup output branches
	tele.SetupInput(ipt);

	// add pid file
	TString pid_file_name;
	pid_file_name.Form(
		"%s%s%s-particle-type-%s%04u.root",
		kGenerateDataPath,
		kParticleIdentifyDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	if (!ipt->AddFriend("pid=tree", pid_file_name)) {
		std::cerr << "Error: Add pid file "
			<< pid_file_name << " failed.\n";
		return -1;
	}
	// particle mass number
	unsigned int mass[4];
	// particle charge number
	unsigned int charge[4];
	// setup input branches
	ipt->SetBranchAddress("pid.mass", mass);
	ipt->SetBranchAddress("pid.charge", charge);

	// output root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-calibration-%04u.root",
		kGenerateDataPath,
		kCalibrationDir,
		(std::string("taf") + name_[3] + "csi").c_str(),
		run_
	);
	// output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// output calibration graphs, the first index is for CsI index
	// the second index is for phi index
	// the third index is for isotopes
	TGraph *gcali[2][theta_types][isotope_types];
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < theta_types; ++j) {
			for (size_t k = 0; k < isotope_types; ++k) {
				gcali[i][j][k] = new TGraph;
			}
		}
	}
	// pid graph with calculated CsI energy
	TGraph *pid = new TGraph;

	// get tafd calibration parameters
	TString param_file_name;
	param_file_name.Form(
		"%s%s%s-alpha-cali-param.txt",
		kGenerateDataPath,
		kCalibrationDir,
		(name_.substr(0, 3) + "d" + name_[3]).c_str()
	);
	std::ifstream fin(param_file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open file "
			<< param_file_name << " failed.\n";
		return -1;
	}
	// tafd calibration parameters
	double tafd_param[2];
	// read parameters from file
	fin >> tafd_param[0] >> tafd_param[1];
	// close file
	fin.close();

	// total number of entries
	long long entries = ipt->GetEntries();
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
		ipt->GetEntry(entry);
		if (tele.particle != 1 || tele.layer[0] != 2) continue;
		// calculate first index in graph array case by phi angle
		size_t csi_index = tele.phi[0][0] > tafd_center_phi[index_] ? 0 : 1;
		// calculate second index in graph array case by theta angle
		size_t theta_index = tele.theta[0][0] < 0.8 ? 0 : 1;
		// calculated CsI energy
		double csi_energy = -1e4;
		// CsI channel number
		double csi_channel = tele.energy[0][1];
		for (size_t i = 0; i < particle_types; ++i) {
			// jump if not the identified particle
			if (
				charge[0] != particles_charge[i]
				|| mass[0] != particles_mass[i]
			) continue;

			// calculate CsI energy
			csi_energy = CalculateCsiEnergy(
				particle_names[i],
				tele.theta[0][0],
				tafd_param[0] + tafd_param[1] * tele.energy[0][0],
				tafd_thick[index_]
			);
			// fill the appropriate graph, the first index is chosen by phi
			// the second index is chosen ny theta and charge number
			// 25000 offset in channel to disinguish different isotopes
			// charge[0]-1 is the isotope, 0 for H, 1 for He
			// i - isotope_index[charge[0]-1] is offset in isotope, this
			// disinguishes between isotopes
			gcali[csi_index][theta_index][charge[0]-1]->AddPoint(
				csi_energy + 200 * (i - isotope_index[charge[0]-1]),
				csi_channel
			);
			pid->AddPoint(csi_energy, tafd_param[0] + tafd_param[1] * tele.energy[0][0]);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fitting
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < theta_types; ++j) {
			// fit H isotopes
			TF1 *fitH = new TF1(
				TString::Format("fH_c%ld_%s", i, theta_names[j]),
				HFit, 0, 600, 4
			);
			fitH->FixParameter(0, 0.0);
			fitH->SetParLimits(2, 0, 1000);
			fitH->SetParLimits(3, 0, 100);
			gcali[i][j][0]->Fit(fitH, "QR+ rob=0.8");
std::cout << "----------------------------\n"
<< "csi " << i << " theta " << theta_names[j] << " H fitting\n"
<< fitH->GetParameter(0) << " " << fitH->GetParameter(1)
<< " " << fitH->GetParameter(2) << " " << fitH->GetParameter(3) << "\n";

			// fit He isotopes
			TF1 *fitHe = new TF1(
				TString::Format("fH_c%ld_%s", i, theta_names[j]),
				HeFit, 0, 400, 3
			);
			fitHe->FixParameter(0, 0.0);
			fitHe->SetParLimits(2, 0, 100);
			gcali[i][j][1]->Fit(fitHe, "QR+");
std::cout << "----------------------------\n"
<< "csi " << i << " theta " << theta_names[j] << " He fitting\n"
<< fitHe->GetParameter(0) << " " << fitHe->GetParameter(1)
<< " " << fitHe->GetParameter(2) << "\n";
		}
	}

	opf->cd();
	// save graphs
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < theta_types; ++j) {
			for (size_t k = 0; k < isotope_types; ++k) {
				gcali[i][j][k]->Write(TString::Format(
					"g_c%ld_%s_%ld", i, theta_names[j], k
				));
			}
		}
	}
	pid->Write("pid");
	// close files
	opf->Close();
	ipf->Close();

	return 0;
}



}		// namespace ribll