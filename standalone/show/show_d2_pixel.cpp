#include <iostream>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPaveText.h>
#include <TString.h>
#include <TTree.h>

#include "include/event/dssd_event.h"

using namespace ribll;

constexpr int fsmin = 9;
constexpr int fsmax = 19;
constexpr int bsmin = 13;
constexpr int bsmax = 22;


int ShowD2Pixel(int run) {
	// input file name
	TString input_file_name = TString::Format(
		"%s%st0d2-fundamental-%04d.root",
		kGenerateDataPath,
		kFundamentalDir,
		run
	);
	// input file
	TFile input_file(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)input_file.Get("tree");
	// input event
	DssdFundamentalEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sshow-t0d2-pixel-%04d.root",
		kGenerateDataPath,
		kShowDir,
		run
	);
	// output file
	TFile output_file(output_file_name, "recreate");
	// histogram of energy in one pixel
	TH1F *hist_energy[32][32];
	for (int fs = fsmin; fs <= fsmax; ++fs) {
		for (int bs = bsmin; bs <= bsmax; ++bs) {
			if (bs == 17) continue;
			hist_energy[fs][bs] = new TH1F(
				TString::Format("hef%db%d", fs, bs), "energy",
				1000, 0, 60000
			);
		}
	}
	// 2D graph of energy peak position
	TH2F hist_beam_mean("hbe", "beam peak position", 32, 0, 32, 32, 0, 32);
	// 2D graph of energy peak sigma
	TH2F hist_beam_sigma("hbs", "beam peak sigma", 32, 0, 32, 32, 0, 32);
	// 2D graph of beam energy resolution
	TH2F hist_beam_resolution("hbr", "beam energy resolution", 32, 0, 32, 32, 0, 32);
	// 2D graph of energy peak position and back strip in mirror
	TH2F hist_beam_mean_mirror(
		"hbem", "beam peak position", 32, 0, 32, 32, 0, 32
	);
	// 2D graph of energy peak sigma and back strip in mirror
	TH2F hist_beam_sigma_mirror(
		"hbsm", "beam peak sigma", 32, 0, 32, 32, 0, 32
	);

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100;
	// show start
	printf("Reading events   0%%");
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get entry
		ipt->GetEntry(entry);
		if (event.front_hit != 1 || event.back_hit != 1) continue;
		if (event.front_strip[0] < fsmin || event.front_strip[0] > fsmax) continue;
		if (event.back_strip[0] < bsmin || event.back_strip[0] > bsmax) continue;
		if (event.front_energy[0] < 1000) continue;
		hist_energy[event.front_strip[0]][event.back_strip[0]]->Fill(
			event.front_energy[0]
		);
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit beam energy
	for (int fs = fsmin; fs <= fsmax; ++fs) {
		for (int bs = bsmin; bs <= bsmax; ++bs) {
			if (bs == 17) continue;
			// find max value
			double max_position = 0.0;
			double max_value = 0.0;
			for (int i = 0; i < hist_energy[fs][bs]->GetNcells(); ++i) {
				if (hist_energy[fs][bs]->GetBinCenter(i) < 10000) continue;
				if (hist_energy[fs][bs]->GetBinContent(i) > max_value) {
					max_value = hist_energy[fs][bs]->GetBinContent(i);
					max_position = hist_energy[fs][bs]->GetBinCenter(i);
				}
			}

			TF1 *energy_fit = new TF1("f1", "gaus", max_position-10000, max_position+10000);
			energy_fit->SetParameter(0, max_value);
			energy_fit->SetParameter(1, max_position);
			energy_fit->SetParameter(2, 1000);
			hist_energy[fs][bs]->Fit(energy_fit, "RQ+");
			double mean = energy_fit->GetParameter(1);
			mean = mean < 0 ? 0.0 : mean;
			double sigma = energy_fit->GetParameter(2);
			sigma = mean < 0 ? 0.0 : sigma;
			hist_beam_mean.SetBinContent(fs+1, bs+1, mean);
			hist_beam_sigma.SetBinContent(fs+1, bs+1, sigma);
			hist_beam_resolution.SetBinContent(fs+1, bs+1, sigma / mean);
			hist_beam_mean_mirror.SetBinContent(fs+1, -(bs+1), mean);
			hist_beam_sigma_mirror.SetBinContent(fs+1, -(bs+1), sigma);
		}
	}

	// save histograms
	for (int fs = fsmin; fs <= fsmax; ++fs) {
		for (int bs = bsmin; bs <= bsmax; ++bs) {
			if (bs == 17) continue;
			hist_energy[fs][bs]->Write();
		}
	}
	// save histograms
	hist_beam_mean.Write();
	hist_beam_sigma.Write();
	hist_beam_resolution.GetXaxis()->SetRangeUser(fsmin-1, fsmax+1);
	hist_beam_resolution.GetYaxis()->SetRangeUser(bsmin-1, bsmax+1);
	hist_beam_resolution.Write();
	hist_beam_mean_mirror.Write();
	hist_beam_sigma_mirror.Write();
	// close files
	output_file.Close();
	input_file.Close();
	return 0;
}


int SummaryD2Pixel(int start_run, int end_run) {
	// output file name
	TString output_file_name = TString::Format(
		"%s%sshow-t0d2-pixel-%04d-%04d.root",
		kGenerateDataPath,
		kShowDir,
		start_run,
		end_run
	);
	// output file
	TFile output_file(output_file_name, "recreate");
	TCanvas resolution_canvas;
	TCanvas over8_canvas;
	TCanvas over5_canvas;

	for (int run = start_run; run <= end_run; ++run) {
		if (run == 628) continue;
		output_file.cd();
		// beam resolution over 8%
		TH2F hist_res_over8(
			TString::Format("h8r%d", run), "beam resolution over 8%%",
			32, 0, 32, 32, 0, 32
		);
		// beam resolution over 5%
		TH2F hist_res_over5(
			TString::Format("h5r%d", run), "beam resolution over 5%%",
			32, 0, 32, 32, 0, 32
		);

		// input file name
		TString input_file_name = TString::Format(
			"%s%sshow-t0d2-pixel-%04d.root",
			kGenerateDataPath,
			kShowDir,
			run
		);
		// input file
		TFile input_file(input_file_name, "read");
		// beam energy resolution histogram
		TH2F *hist_beam_resolution = (TH2F*)input_file.Get("hbr");
		if (!hist_beam_resolution) {
			std::cerr << "Error: Get histogram from "
				<< input_file_name << "\n";
			continue;
		}

		for (int fs = fsmin; fs < fsmax; ++fs) {
			for (int bs = bsmin; bs < bsmax; ++bs) {
				if (bs == 17) continue;
				double resolution = hist_beam_resolution->GetBinContent(fs+1, bs+1);
				hist_res_over8.SetBinContent(
					fs+1, bs+1, resolution > 0.08 ? 1.0 : 0.0
				);
				hist_res_over5.SetBinContent(
					fs+1, bs+1, resolution > 0.05 ? 1.0 : 0.0
				);
			}
		}

		hist_res_over8.GetXaxis()->SetRangeUser(fsmin-1, fsmax+1);
		hist_res_over8.GetYaxis()->SetRangeUser(bsmin-1, bsmax+1);
		hist_res_over5.GetXaxis()->SetRangeUser(fsmin-1, fsmax+1);
		hist_res_over5.GetYaxis()->SetRangeUser(bsmin-1, bsmax+1);
		// run number
		TPaveText text(fsmax-1, bsmin-1, fsmax+1, bsmin+1);
		text.AddText(TString::Format("%d", run));
		text.SetTextColor(kRed);
		text.SetBorderSize(0);
		text.SetFillColor(0);
		text.SetFillStyle(0);

		resolution_canvas.cd();
		hist_beam_resolution->Draw("colz");
		text.Draw();
		resolution_canvas.Print(TString::Format(
			"%s%st0d2-pixel-shift-%d-%d.gif%s",
			kGenerateDataPath,
			kShowDir,
			start_run,
			end_run,
			run == start_run ? "" :
				run == end_run ? "++" : "+50"
		));

		over8_canvas.cd();
		hist_res_over8.Draw("colz");
		text.Draw();
		over8_canvas.Print(TString::Format(
			"%s%st0d2-pixel-shift-over8-%d-%d.gif%s",
			kGenerateDataPath,
			kShowDir,
			start_run,
			end_run,
			run == start_run ? "" :
				run == end_run ? "++" : "+50"
		));

		over5_canvas.cd();
		hist_res_over5.Draw("colz");
		text.Draw();
		over5_canvas.Print(TString::Format(
			"%s%st0d2-pixel-shift-over5-%d-%d.gif%s",
			kGenerateDataPath,
			kShowDir,
			start_run,
			end_run,
			run == start_run ? "" :
				run == end_run ? "++" : "+50"
		));

		output_file.cd();
		hist_beam_resolution->Write();
	}
	// close file
	output_file.Close();
	return 0;
}

int main(int argc, char **argv) {
	if (argc < 2 || argc > 3) {
		std::cout << "Usage: " << argv[0] << " run [end_run]\n"
			<< "  run           run to show\n"
			<< "  end_run       summary runs\n";
		return -1;
	}

	// run number
	int run = atoi(argv[1]);
	if (argc == 2) {
		if (ShowD2Pixel(run)) {
			return -1;
		}
	} else {
		int end_run = atoi(argv[2]);
		if (SummaryD2Pixel(run, end_run)) {
			return -1;
		}
	}
	return 0;
}