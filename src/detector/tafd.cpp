#include "include/detector/tafd.h"

#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFitResult.h>

namespace ribll {

const ROOT::Math::Polar3DVector tafd_center(135.0, 0.0, 0.0);
const std::pair<double, double> tafd_radius_ranges[6] = {
	{68, 170.5},
	{68, 170.5},
	{68, 170.5},
	{68, 170.5},
	{68, 170.5},
	{68, 170.5}
};
const std::pair<double, double> tafd_phi_ranges[6] = {
	{117.6*TMath::DegToRad(), 62.4*TMath::DegToRad()},
	{57.6*TMath::DegToRad(), 2.4*TMath::DegToRad()},
	{-2.4*TMath::DegToRad(), -57.6*TMath::DegToRad()},
	{-62.4*TMath::DegToRad(), -117.6*TMath::DegToRad()},
	{-122.4*TMath::DegToRad(), -177.6*TMath::DegToRad()},
	{177.6*TMath::DegToRad(), 122.4*TMath::DegToRad()}
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


Tafd::Tafd(unsigned int run, unsigned int index, const std::string &tag)
: Adssd(run, "tafd"+std::to_string(index), tag)
, index_(index) {

	center_ = tafd_center;
	radius_range_ = tafd_radius_ranges[index];
	phi_range_ = tafd_phi_ranges[index];
}


int Tafd::MatchTrigger(double, double) {
	if (name_ == "tafd0" || name_ == "tafd1") {
		return Detector::VmeMatchTrigger<DssdFundamentalEvent>();
	}
	std::cerr << "Error: Use ExtractTrigger instead.\n";
	return -1;
}


int Tafd::ExtractTrigger(
	double window_left,
	double window_right
) {
	if (name_ == "tafd0" || name_ == "tafd1") {
		std::cerr << "Error: Use MatchTrigger instead.\n";
		return -1;
	}
	return Dssd::ExtractTrigger(window_left, window_right);
}


int Tafd::NormalizeSides(TChain *chain, bool iteration) {
	if (SideNormalize(chain, 0, 4, iteration)) {
		std::cerr << "Error: Normalize first side failed.\n";
		return -1;
	}
	if (SideNormalize(chain, 1, 1, iteration)) {
		std::cerr << "Error: Normalize second side failed.\n";
		return -1;
	}
	return 0;
}


bool Tafd::NormEnergyCheck(size_t, const DssdFundamentalEvent &event) const {
	if (index_ == 0) {
		if (event.front_energy[0] > 1e4 || event.back_energy[0] > 1e4) {
			return false;
		}
	} else if (index_ == 1) {
		if (event.front_energy[0] > 1e4 || event.back_energy[0] > 1e4) {
			return false;
		}
	} else {
		return false;
	}
	return true;
}



int Tafd::Calibrate() {
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


}