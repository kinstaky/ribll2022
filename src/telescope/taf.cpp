#include "include/telescope/taf.h"

#include <fstream>

#include <TCutG.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH1F.h>
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

// run used for alpha calibration
const unsigned int alpha_calibration_run[6] = {816, 816, 825, 825, 825, 825};
// energy of alpha source, in MeV
const double alpha_energy[3] = {5.157, 5.486, 5.805};
// histogram range of alpha source peaks
const double alpha_hist_range[6][3] = {
	{500, 1700, 2200},
	{700, 1800, 2500},
	{500, 11500, 16500},
	{500, 12000, 17000},
	{500, 11500, 16500},
	{500, 11500, 16500}
};
// fit range of alpha source peaks
const double alpha_fit_range[6][16][6] = {
	{
		{1835.0, 1860.0, 1955.0, 1980.0, 2070.0, 2100.0},
		{1865.0, 1890.0, 1990.0, 2015.0, 2100.0, 2130.0},
		{1830.0, 1860.0, 1950.0, 1975.0, 2060.0, 2095.0},
		{1840.0, 1865.0, 1960.0, 1985.0, 2070.0, 2095.0},
		{1810.0, 1835.0, 1925.0, 1955.0, 2040.0, 2065.0},
		{1820.0, 1845.0, 1940.0, 1965.0, 2050.0, 2075.0},
		{1830.0, 1860.0, 1950.0, 1980.0, 2070.0, 2095.0},
		{1800.0, 1830.0, 1920.0, 1945.0, 2030.0, 2060.0},
		{1780.0, 1805.0, 1895.0, 1920.0, 2005.0, 2035.0},
		{1820.0, 1850.0, 1940.0, 1970.0, 2050.0, 2080.0},
		{1790.0, 1820.0, 1910.0, 1940.0, 2020.0, 2060.0},
		{1770.0, 1800.0, 1890.0, 1915.0, 1995.0, 2030.0},
		{1765.0, 1800.0, 1880.0, 1915.0, 1990.0, 2030.0},
		{1750.0, 1790.0, 1870.0, 1900.0, 1980.0, 2015.0},
		{1740.0, 1780.0, 1855.0, 1890.0, 1865.0, 2005.0},
		{1770.0, 1810.0, 1890.0, 1930.0, 2005.0, 2045.0}
	},
	{
		{2050.0, 2080.0, 2180.0, 2220.0, 2305.0, 2345.0},
		{2045.0, 2085.0, 2180.0, 2230.0, 2310.0, 2355.0},
		{2005.0, 2040.0, 2130.0, 2170.0, 2260.0, 2300.0},
		{1995.0, 2025.0, 2125.0, 2155.0, 2250.0, 2280.0},
		{2025.0, 2060.0, 2145.0, 2195.0, 2270.0, 2320.0},
		{2025.0, 2075.0, 2160.0, 2205.0, 2285.0, 2330.0},
		{2010.0, 2050.0, 2140.0, 2185.0, 2265.0, 2300.0},
		{2040.0, 2070.0, 2175.0, 2210.0, 2300.0, 2330.0},
		{1970.0, 2010.0, 2100.0, 2135.0, 2220.0, 2260.0},
		{2035.0, 2075.0, 2165.0, 2200.0, 2270.0, 2330.0},
		{2045.0, 2085.0, 2180.0, 2220.0, 2300.0, 2345.0},
		{1945.0, 1995.0, 2075.0, 2115.0, 2195.0, 2245.0},
		{1940.0, 1970.0, 2060.0, 2105.0, 2180.0, 2225.0},
		{1975.0, 2040.0, 2105.0, 2165.0, 2225.0, 2295.0},
		{1920.0, 1965.0, 2050.0, 2090.0, 2165.0, 2210.0},
		{1885.0, 1935.0, 2010.0, 2050.0, 2130.0, 2170.0}
	},
	{
		{13880.0, 14010.0, 14810.0, 14950.0, 15700.0, 15820.0},
		{13390.0, 13560.0, 14280.0, 14430.0, 15140.0, 15270.0},
		{13240.0, 13380.0, 14120.0, 14250.0, 14980.0, 15100.0},
		{13130.0, 13320.0, 14010.0, 14220.0, 14830.0, 15050.0},
		{12840.0, 13030.0, 13690.0, 13900.0, 14510.0, 14710.0},
		{12670.0, 12850.0, 13510.0, 13670.0, 14310.0, 14490.0},
		{12770.0, 12970.0, 13630.0, 13810.0, 14450.0, 14620.0},
		{12830.0, 13040.0, 13680.0, 13900.0, 14500.0, 14690.0},
		{12370.0, 12570.0, 13200.0, 13420.0, 13980.0, 14190.0},
		{12480.0, 12660.0, 13320.0, 13520.0, 14100.0, 14310.0},
		{12540.0, 12730.0, 13350.0, 13570.0, 14140.0, 14360.0},
		{12610.0, 12800.0, 13460.0, 13690.0, 14260.0, 14480.0},
		{12580.0, 13050.0, 13540.0, 14000.0, 14340.0, 14810.0},
		{12240.0, 12420.0, 13070.0, 13250.0, 13840.0, 14070.0},
		{12510.0, 12690.0, 13310.0, 13520.0, 14120.0, 14330.0},
		{12570.0, 12790.0, 13400.0, 13640.0, 14210.0, 14440.0}
	},
	{
		{14300.0, 14600.0, 15300.0, 15500.0, 16200.0, 16400.0},
		{13500.0, 14200.0, 14500.0, 15100.0, 15500.0, 16000.0},
		{13400.0, 13700.0, 14400.0, 14550.0, 15250.0, 15400.0},
		{13350.0, 13550.0, 14250.0, 14450.0, 15150.0, 15400.0},
		{13750.0, 13950.0, 14650.0, 14850.0, 15550.0, 15700.0},
		{13250.0, 13500.0, 14200.0, 14500.0, 15100.0, 15300.0},
		{13200.0, 13600.0, 14200.0, 14500.0, 15000.0, 15400.0},
		{13400.0, 13600.0, 14200.0, 14500.0, 15000.0, 15400.0},
		{13600.0, 13800.0, 14500.0, 14700.0, 15300.0, 15600.0},
		{13100.0, 13500.0, 14100.0, 14400.0, 14900.0, 15200.0},
		{13300.0, 13700.0, 14300.0, 14600.0, 15100.0, 15500.0},
		{13200.0, 13500.0, 14100.0, 14400.0, 14900.0, 15300.0},
		{13300.0, 13600.0, 14200.0, 14500.0, 15000.0, 15400.0},
		{13200.0, 13500.0, 14100.0, 14500.0, 15000.0, 15400.0},
		{12500.0, 12900.0, 13500.0, 13580.0, 14300.0, 14600.0},
		{12800.0, 13100.0, 13700.0, 14000.0, 14600.0, 14900.0}
	},
	{
		{13800.0, 14100.0, 14700.0, 15000.0, 15600.0, 15900.0},
		{13300.0, 13500.0, 14200.0, 14400.0, 15000.0, 15300.0},
		{13600.0, 13800.0, 14500.0, 14700.0, 15400.0, 15600.0},
		{13400.0, 13700.0, 14300.0, 14600.0, 15200.0, 15500.0},
		{13400.0, 13700.0, 14300.0, 14600.0, 15200.0, 15600.0},
		{13100.0, 13300.0, 13900.0, 14200.0, 14700.0, 15100.0},
		{12800.0, 13100.0, 13600.0, 14000.0, 14500.0, 14800.0},
		{12900.0, 13200.0, 13800.0, 14100.0, 14600.0, 14900.0},
		{12800.0, 13100.0, 13700.0, 14000.0, 14500.0, 14900.0},
		{12700.0, 13000.0, 13500.0, 13800.0, 14300.0, 14700.0},
		{12900.0, 13200.0, 13900.0, 14200.0, 14700.0, 15000.0},
		{12900.0, 13200.0, 13700.0, 14100.0, 14600.0, 14900.0},
		{12600.0, 12900.0, 13500.0, 13800.0, 14300.0, 14600.0},
		{12700.0, 13000.0, 13600.0, 13900.0, 14400.0, 14700.0},
		{12400.0, 12700.0, 13200.0, 13600.0, 14100.0, 14400.0},
		{12500.0, 12800.0, 13300.0, 13700.0, 14100.0, 14500.0}
	},
	{
		{13600.0, 13800.0, 14500.0, 14700.0, 15400.0, 15550.0},
		{13300.0, 13500.0, 14200.0, 14400.0, 15100.0, 15300.0},
		{13400.0, 13600.0, 14300.0, 14500.0, 15100.0, 15300.0},
		{13000.0, 13200.0, 13900.0, 14100.0, 14700.0, 14900.0},
		{13300.0, 13500.0, 14200.0, 14500.0, 15100.0, 15400.0},
		{13000.0, 13200.0, 13800.0, 14100.0, 14700.0, 14900.0},
		{13200.0, 13400.0, 14100.0, 14400.0, 15000.0, 15200.0},
		{13100.0, 13400.0, 14000.0, 14200.0, 14900.0, 15100.0},
		{13000.0, 13300.0, 13900.0, 14200.0, 14700.0, 15000.0},
		{12700.0, 12900.0, 13500.0, 13800.0, 14400.0, 14600.0},
		{12700.0, 12900.0, 13500.0, 13800.0, 14400.0, 14600.0},
		{12900.0, 13400.0, 13800.0, 14300.0, 14700.0, 15200.0},
		{12900.0, 13200.0, 13800.0, 14100.0, 14600.0, 15000.0},
		{12400.0, 12700.0, 13200.0, 13500.0, 14100.0, 14300.0},
		{12200.0, 12500.0, 13100.0, 13300.0, 13900.0, 14200.0},
		{12200.0, 12500.0, 13100.0, 13300.0, 13900.0, 14200.0}
	}
};

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
// isotope types
const size_t isotope_types = 2;
// isotope start index in particles array above
const size_t isotope_index[isotope_types] = {0, 3};
// phi angle types
const size_t theta_types = 2;
// phi angle names
const char* const theta_names[theta_types] = {
	"l8", "g8"
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
			cut_index += tele.theta[0][0] > 0.8 ? 2 : 0;
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


int Taf::Calibrate() {
	// if (AlphaCalibrate()) return -1;
	if (CsiCalibrate()) return -1;
	return 0;
}


int Taf::AlphaCalibrate() {
	// calibrate tafd with alpha source
	// output root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-alpha-calibration-%04u.root",
		kGenerateDataPath,
		kCalibrationDir,
		(name_.substr(0, 3) + "d" + name_[3]).c_str(),
		alpha_calibration_run[index_]
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output energy histogram for each strip
	std::vector<TH1F> hist_energy;
	for (int i = 0; i < 16; ++i) {
		hist_energy.emplace_back(
			TString::Format("he%i", i), "fit alpha source peaks",
			alpha_hist_range[index_][0],
			alpha_hist_range[index_][1], alpha_hist_range[index_][2]
		);
	}
	if (index_ < 2) {
		// VME
		// input file name
		TString input_file_name;
		input_file_name.Form(
			"%s%s%04u.root",
			kCrate3Path,
			kCrate3FileName,
			alpha_calibration_run[index_]
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
		// energy values
		int madc[2][32];
		// setup input branches
		ipt->SetBranchAddress("madc", madc);

		// madc module of this tafd
		size_t mod = vtaf_front_module[index_];
		// madc channel of this tafd
		size_t ch = vtaf_front_channel[index_];

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filling histogram   0%%");
		fflush(stdout);
		// fill energy to histogram
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			ipt->GetEntry(entry);
			for (size_t i = 0; i < 16; ++i) {
				if (madc[mod][ch+i] > 1000 && madc[mod][ch+i] < 5000) {
					hist_energy[i].Fill(madc[mod][ch+i]);
				}
			}
		}
		// show finish
		printf("\b\b\b\b100%%\n");
		// close input file
		ipf.Close();

	} else {
		// XIA
		// input file name
		TString input_file_name;
		input_file_name.Form(
			"%s%s%s-map-%04u.root",
			kGenerateDataPath,
			kMappingDir,
			(name_.substr(0, 3) + "d" + name_[3]).c_str(),
			alpha_calibration_run[index_]
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
		// input map event
		DssdMapEvent event;
		// setup input branches
		event.SetupInput(ipt);

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filling histogram   0%%");
		fflush(stdout);
		// fill energy to histogram
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			ipt->GetEntry(entry);
			// jump if back side
			if (event.side == 1) continue;
			// jump if energy out of range
			if (event.energy < 10000 || event.energy > 20000) continue;
			hist_energy[event.strip].Fill(event.energy);

			// trick to prevent error in bad strip
			if (index_ == 5 && event.strip == 14) {
				hist_energy[15].Fill(event.energy);
			}
		}
		// show finish
		printf("\b\b\b\b100%%\n");
		// close input file
		ipf.Close();
	}


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

	// fit alpha peaks
	TFitResultPtr fit_result[16][3];
	for (size_t i = 0; i < 16; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			TF1 fpeak(
				TString::Format("s%ldp%ld", i, j),
				"gaus",
				alpha_fit_range[index_][i][2*j],
				alpha_fit_range[index_][i][2*j+1]
			);
			fpeak.SetParameter(0, 1000);
			fpeak.SetParameter(1, alpha_fit_range[index_][i][2*j]+10);
			fit_result[i][j] = hist_energy[i].Fit(&fpeak, "QRS+");
		}
	}

	// calibrate
	TGraph gcali[16];
	// calibrate results
	TFitResultPtr cali_result[16];
	// loop to calibrate
	for (size_t i = 0; i < 16; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			gcali[i].AddPoint(fit_result[i][j]->Parameter(1), alpha_energy[j]);
		}
		TF1 fcali(
			"fcali", "pol1",
			alpha_hist_range[index_][1], alpha_hist_range[index_][2]
		);
		cali_result[i] = gcali[i].Fit(&fcali, "QRS+");
	}
	// show calibrate result
	for (size_t i = 0; i < 16; ++i) {
		std::cout << cali_result[i]->Parameter(0) << " "
			<< cali_result[i]->Parameter(1) << "\n";
	}

	// output fit results
	for (size_t i = 0; i < 16; ++i) {
		fout << cali_result[i]->Parameter(0) << " "
			<< cali_result[i]->Parameter(1) << " "
			<< cali_result[i]->Chi2() << " "
			<< cali_result[i]->Ndf() << " "
			<< cali_result[i]->Chi2() / cali_result[i]->Ndf() << "\n";
	}
	for (size_t i = 0; i < 16; ++i) {
		fout << "-----------------\n";
		for (size_t j = 0; j < 3; ++j) {
			fout << fit_result[i][j]->Parameter(0) << " "
				<< fit_result[i][j]->Parameter(1) << " "
				<< fit_result[i][j]->Parameter(2) << " "
				<< fit_result[i][j]->Chi2() << " "
				<< fit_result[i][j]->Ndf() << " "
				<< fit_result[i][j]->Chi2() / fit_result[i][j]->Ndf() << "\n";
		}
	}
	// close file
	fout.close();

	// save histogram and close files
	opf.cd();
	for (size_t i = 0; i < 16; ++i)	hist_energy[i].Write();
	for (size_t i = 0; i < 16; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			fit_result[i][j]->Write(TString::Format("rs%ldp%ld", i, j));
		}
	}
	for (size_t i = 0; i < 16; ++i) {
		gcali[i].Write(TString::Format("gcali%ld", i));
	}
	for (size_t i = 0; i < 16; ++i) {
		cali_result[i]->Write(TString::Format("rcali%ld", i));
	}
	// close output file
	opf.Close();

	return 0;
}


int Taf::Particle() {
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
	// particle tree
	TTree opt("tree", "taf particles");
	// output particle event
	ParticleEvent particle_event;
	// setup output branches
	particle_event.SetupOutput(&opt);

	// CsI energy calculators
	elc::CsiEnergyCalculator csi_calculators[4]{
		elc::CsiEnergyCalculator("1H"),
		elc::CsiEnergyCalculator("2H"),
		elc::CsiEnergyCalculator("3H"),
		elc::CsiEnergyCalculator("4He")
	};

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
			// get calculator index, here is a trick.
			// In calculators array, 0-1H, 1-2H, 2-3H, 3-4He.
			// The mass number is index+1 in coincidence.
			size_t calculator_index = type_event.mass[0] - 1;
			// calculate CsI energy from TAFD energy
			double csi_energy = csi_calculators[calculator_index].Energy(
				ta_event.theta[0][0], si_energy, tafd_thick[index_]
			);
			// calculate the total energy
			particle_event.energy[0] = si_energy + csi_energy;

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

	// save particle tree
	particle_file.cd();
	opt.Write();
	// close files
	particle_file.Close();
	telescope_file.Close();

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
	TaEvent tele;
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

	// CsI energy calculator
	std::vector<elc::CsiEnergyCalculator> csi_calculators;
	for (size_t i = 0; i < particle_types; ++i) {
		csi_calculators.emplace_back(particle_names[i]);
	}

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
		if (tele.num != 1 || tele.flag[0] != 0x3) continue;
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
			csi_energy = csi_calculators[i].Energy(
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