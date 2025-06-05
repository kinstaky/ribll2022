from scipy import stats

GENERATE_PATH = "/mnt/d/data/ribll2022/"
with open(GENERATE_PATH+'spectrum/significance-test-fixed.txt', 'r') as file:
	for line in file:
		peak, lamba = line.strip().split(',')
		print(peak, stats.norm.ppf(stats.chi2.cdf(float(lamba)*2, 1)))