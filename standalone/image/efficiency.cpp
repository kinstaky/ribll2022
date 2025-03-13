// #define T0
// #define TAF
#define ALL

void efficiency() {
	TFile *ipf = new TFile(
		"/mnt/d/data/ribll2022/simulate/efficiency-0002-10x-avg.root",
		"read"
	);
	TGraph *t0pg0 = (TGraph*)ipf->Get("t0pg0");
	TGraph *t0pg1 = (TGraph*)ipf->Get("t0pg1");
	TGraph *t0pg2 = (TGraph*)ipf->Get("t0pg2");
	TGraph *tafg0 = (TGraph*)ipf->Get("tafg0");
	TGraph *tafg1 = (TGraph*)ipf->Get("tafg1");
	TGraph *tafg2 = (TGraph*)ipf->Get("tafg2");
	TGraph *g0 = (TGraph*)ipf->Get("g0");
	TGraph *g1 = (TGraph*)ipf->Get("g1");
	TGraph *g2 = (TGraph*)ipf->Get("g2");


	t0pg0->SetLineColor(kRed);
	t0pg1->SetLineColor(kBlack);
	t0pg2->SetLineColor(kBlue);
	tafg0->SetLineColor(kRed);
	tafg1->SetLineColor(kBlack);
	tafg2->SetLineColor(kBlue);
	g0->SetLineColor(kRed);
	g1->SetLineColor(kBlack);
	g2->SetLineColor(kBlue);

#ifdef T0
	TMultiGraph *mg = new TMultiGraph;
	mg->Add(t0pg0);
	mg->Add(t0pg1);
	mg->Add(t0pg2);
	mg->Draw("APL");
#endif

#ifdef TAF
	TMultiGraph *mg = new TMultiGraph;
	mg->Add(tafg0);
	mg->Add(tafg1);
	mg->Add(tafg2);
	mg->Draw("APL");
#endif

#ifdef ALL
	TMultiGraph *mg = new TMultiGraph;
	mg->Add(g0);
	mg->Add(g1);
	mg->Add(g2);
	mg->Draw("APL");
#endif

}