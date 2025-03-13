#define PPAC3
// #define PPAC2
// #define PPAC1AI
// #define PPAC1AF
// #define PPAC1FI
// #define PPAC1FF

void target_from_ppac() {
    // Open the ROOT file
	TString file_name = "/mnt/d/data/ribll2022/simulate/target-resolution.root";
    TFile *ipf = new TFile(file_name, "read");
#ifdef PPAC3
    TH1F *h3ptx = (TH1F*)ipf->Get("h3ptx");
	h3ptx->SetStats(0);
	h3ptx->Draw("hist");
#endif

#ifdef PPAC2
	TH1F *h2ptx0 = (TH1F*)ipf->Get("h2ptx0");
	TH1F *h2ptx1 = (TH1F*)ipf->Get("h2ptx1");
	TH1F *h2ptx2 = (TH1F*)ipf->Get("h2ptx2");
	// h2ptx0->SetStats(0);
	// h2ptx1->SetStats(0);
	// h2ptx2->SetStats(0);
	h2ptx0->SetLineColor(kRed);
	h2ptx1->SetLineColor(kBlue);
	h2ptx2->SetLineColor(kBlack);
	THStack *hs = new THStack("hs", "Stacked 1D histograms");
	hs->Add(h2ptx0);
	hs->Add(h2ptx1);
	hs->Add(h2ptx2);
	hs->Draw("nostack hist");
#endif

#ifdef PPAC1AI
	TH1F *h1patx0 = (TH1F*)ipf->Get("h1patx0");
	h1patx0->SetStats(0);
	h1patx0->Draw("hist");
#endif

#ifdef PPAC1AF
	TH1F *h1patx0 = (TH1F*)ipf->Get("h1patx0");
	TH1F *h1patx1 = (TH1F*)ipf->Get("h1patx1");
	TH1F *h1patx2 = (TH1F*)ipf->Get("h1patx2");
	h1patx0->SetLineColor(kRed);
	h1patx1->SetLineColor(kBlue);
	h1patx2->SetLineColor(kBlack);
	THStack *hs = new THStack("hs", "Stacked 1D histograms");
	hs->Add(h1patx0);
	hs->Add(h1patx1);
	hs->Add(h1patx2);
	hs->Draw("nostack hist");
#endif


#ifdef PPAC1FI
	TH1F *h1pftx00 = (TH1F*)ipf->Get("h1pftx00");
	h1pftx00->SetStats(0);
	h1pftx00->Draw("hist");
#endif


#ifdef PPAC1FF
	TH1F *h1pftx00 = (TH1F*)ipf->Get("h1pftx00");
	TH1F *h1pftx10 = (TH1F*)ipf->Get("h1pftx10");
	TH1F *h1pftx20 = (TH1F*)ipf->Get("h1pftx20");
	h1pftx00->SetLineColor(kRed);
	h1pftx10->SetLineColor(kBlue);
	h1pftx20->SetLineColor(kBlack);
	THStack *hs = new THStack("hs", "Stacked 1D histograms");
	hs->Add(h1pftx00);
	hs->Add(h1pftx10);
	hs->Add(h1pftx20);
	hs->Draw("nostack hist");

#endif

}
