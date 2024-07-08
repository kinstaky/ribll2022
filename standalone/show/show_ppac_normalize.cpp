#include <iostream>
#include <fstream>

#include <TString.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFile.h>

#include "include/defs.h"

using namespace ribll;

int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%sppac-normalize.root",
		kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// graph
	TGraph gx[3], gy[3];
	// set line color
	gx[0].SetLineColor(kBlack);
	gx[1].SetLineColor(kBlue);
	gx[2].SetLineColor(kRed);
	gy[0].SetLineColor(kBlack);
	gy[1].SetLineColor(kBlue);
	gy[2].SetLineColor(kRed);
	// set marker style
	gx[0].SetMarkerStyle(20);
	gx[1].SetMarkerStyle(21);
	gx[2].SetMarkerStyle(22);
	gy[0].SetMarkerStyle(20);
	gy[1].SetMarkerStyle(21);
	gy[2].SetMarkerStyle(22);
	// set marker color
	gx[0].SetMarkerColor(kBlack);
	gx[1].SetMarkerColor(kBlue);
	gx[2].SetMarkerColor(kRed);
	gy[0].SetMarkerColor(kBlack);
	gy[1].SetMarkerColor(kBlue);
	gy[2].SetMarkerColor(kRed);

	// multigraph
	TMultiGraph mgx, mgy;

	// output csv file
	TString csv_file_name = TString::Format(
		"%s%sppac-normalize.csv", kGenerateDataPath, kShowDir
	);
	std::ofstream csv_out(csv_file_name);
	csv_out << "run, x0, x1, x2, y0, y1, y2\n";

	for (int run = 618; run <= 716; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		
		// show file name
		TString normalize_file_name = TString::Format(
			"%s%sppac-normalize-%04d.txt",
			kGenerateDataPath, kNormalizeDir, run
		);
		// read file
		std::ifstream fin(normalize_file_name);
		if (!fin.good()) {
			std::cerr << "Error: Read file "
				<< normalize_file_name << " failed.\n";
		}
		double x[3], y[3];
		fin >> x[0] >> x[1] >> x[2] >> y[0] >> y[1] >> y[2];
		fin.close();

		csv_out << run << ", " << x[0] << ", " << x[1] << ", " << x[2]
			<< ", " << y[0] << ", " << y[1] << ", " << y[2] << "\n";

		for (int i = 0; i < 3; ++i) {
			gx[i].AddPoint(run, x[i]);
			gy[i].AddPoint(run, y[i]);
		}		
	}

	csv_out.close();

	mgx.Add(gx);
	mgx.Add(gx+1);
	mgx.Add(gx+2);
	mgy.Add(gy);
	mgy.Add(gy+1);
	mgy.Add(gy+2);

	gx[0].Write("gx0");
	gx[1].Write("gx1");
	gx[2].Write("gx2");
	gy[0].Write("gy0");
	gy[1].Write("gy1");
	gy[2].Write("gy2");
	mgx.Write("mgx");
	mgy.Write("mgy");
	opf.Close();
	return 0;
}