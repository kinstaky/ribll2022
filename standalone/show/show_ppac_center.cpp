#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <TF1.h>

#include "include/event/particle_event.h"

using namespace ribll;

constexpr double d1z = 100.0;
constexpr double ppac_xz[3] = {-695.2, -454.2, -275.2};
constexpr double ppac_yz[3] = {-689.2, -448.2, -269.2};
constexpr double xc[3] = {0.0, -2.29, -3.43};
constexpr double yc[3] = {0.0, 0.85, 1.82};

void FitAndFill(const double *x, const double *y, const double x0, const double y0, TH1F *h) {
	double a[3];
	for (size_t i = 0; i < 3; ++i) {
		a[i] = x[i] - x0;
	}
	double b[3];
	for (size_t i = 0; i < 3; ++i) {
		b[i] = y0 - y[i];
	}
	double k = -(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	for (size_t i = 0; i < 3; ++i) {
		h[i].Fill(a[i]*k + b[i]);
	}
}


double GausPol0(double *x, double *par) {
	return par[3] + par[0] * exp((x[0]-par[1])*(x[0]-par[1])/2.0/par[2]/par[2]);
}

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			<< "  run               Set the run number.\n";
		return -1;
	}
	unsigned int run = atoi(argv[1]);

	// t0 file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-particle-%04u.root",
		kGenerateDataPath, kParticleDir, run
	);
	// t0 file
	TFile t0_file(t0_file_name, "read");
	// t0 tree
	TTree *tree = (TTree*)t0_file.Get("tree");
	if (!tree) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// xppac file name
	TString xppac_file_name;
	xppac_file_name.Form(
		"%s%sxppac-particle-%04u.root",
		kGenerateDataPath, kParticleDir, run
	);
	tree->AddFriend("xppac=tree", xppac_file_name);
	// input t0 event
	ParticleEvent t0_event;
	// input particle type event
	ParticleEvent xppac_event;
	unsigned short xflag;
	unsigned short yflag;
	// setup branches
	t0_event.SetupInput(tree);
	xppac_event.SetupInput(tree, "xppac.");
	tree->SetBranchAddress("xppac.xflag", &xflag);
	tree->SetBranchAddress("xppac.yflag", &yflag);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sxppac-center-%04u.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// TH1F hd("hd", "distance distribution", 1000, -10, 10);
	TH1F hdx[3]{
		TH1F("hdx0", "#Deltax0", 1000, -10, 10),
		TH1F("hdx1", "#Deltax1", 1000, -10, 10),
		TH1F("hdx2", "#Deltax2", 1000, -10, 10)
	};
	TH1F hdy[3]{
		TH1F("hdy0", "#Deltay0", 1000, -10, 10),
		TH1F("hdy1", "#Deltay1", 1000, -10, 10),
		TH1F("hdy2", "#Deltay2", 1000, -10, 10)
	};


	for (long long entry = 0; entry < tree->GetEntriesFast(); ++entry) {
		tree->GetEntry(entry);
		if (
			t0_event.num != 1 || t0_event.mass[0] != 14
			|| xflag != 0x7 || yflag != 0x7
			// yflag != 0x7
		) continue;

		// double ratio = (ppac_xz[0] - ppac_xz[1]) / (ppac_xz[0] - d1z);
		// double x1 = xppac_event.x[0]+xc[0] - xppac_event.x[1]-xc[1];
		// double x2 = xppac_event.x[0]+xc[0] - t0_event.x[0];


		// double ratio = (ppac_yz[2] - ppac_yz[0]) / (ppac_yz[0] - d1z);
		// double y1 = xppac_event.y[2]+yc[2] - t0_event.y[0];
		// double y2 = xppac_event.y[0]+yc[0] - t0_event.y[0];

		// double ratio = (ppac_yz[1] - ppac_yz[0]) / (ppac_yz[0] - d1z);
		// double y1 = xppac_event.y[1]+yc[1] - t0_event.y[0];
		// double y2 = xppac_event.y[0]+yc[0] - t0_event.y[0];

		// double ratio = (ppac_yz[0] - ppac_yz[2]) / (ppac_yz[1] - ppac_yz[2]);
		// double y1 = xppac_event.y[0]+yc[0] - xppac_event.y[2]-yc[2];
		// double y2 = xppac_event.y[1]+yc[1] - xppac_event.y[2]-yc[2];
		// hd.Fill(y2*ratio - y1);

		FitAndFill(ppac_xz, xppac_event.x, d1z, t0_event.x[0], hdx);
		FitAndFill(ppac_yz, xppac_event.y, d1z, t0_event.y[0], hdy);
	}

	for (int i = 0; i < 3; ++i) {
		TF1 fx(TString::Format("fx%d", i), "gaus", -10, 10);
		fx.SetParameter(0, 20);
		fx.SetParameter(1, 0.0);
		fx.SetParameter(2, 1.0);
		// fx.SetParameter(3, 1.0);
		hdx[i].Fit(&fx, "QR+");
		std::cout << "x" << i << "  " << fx.GetParameter(1) << "\n";

		TF1 fy(TString::Format("fy%d", i), "gaus", -10, 10);
		fy.SetParameter(0, 20);
		fy.SetParameter(1, 0.0);
		fy.SetParameter(2, 1.0);
		// fy.SetParameter(3, 1.0);
		hdy[i].Fit(&fy, "QR+");
		std::cout << "y" << i << "  " << fy.GetParameter(1) << "\n";
	}

	// hd.Write();
	for (size_t i = 0; i < 3; ++i) {
		hdx[i].Write();
		hdy[i].Write();
	}
	// close files
	opf.Close();
	t0_file.Close();
	return 0;
}
