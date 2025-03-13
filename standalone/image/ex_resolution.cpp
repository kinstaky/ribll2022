void ex_resolution() {
	TFile *ipf = new TFile(
		"/mnt/d/data/ribll2022/simulate/resolution-0002.root",
		"read"
	);
	TGraph *g0 = (TGraph*)ipf->Get("g0");

	int avg = 7;
	int n = g0->GetN();
	double *x = g0->GetX();
	double *y = g0->GetY();

	TGraph *g_avg = new TGraph;
	for (int i = 0; i < n-avg/2; ++i) {
		if (i < avg/2) {
			g_avg->AddPoint(x[i], y[i]);
		} else {
			double sum = 0.0;
			for (int j = -avg/2; j <= avg/2; ++j) {
				sum += y[i+j];
			}
			g_avg->AddPoint(x[i], sum / double(avg));
		}
	}
	g_avg->Draw("APL");

}