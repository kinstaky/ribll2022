import ROOT
from scipy.optimize import fsolve
from array import array
from functools import partial
import numpy as np
from math import sqrt
from tqdm import tqdm

GENERATE_PATH = "/mnt/d/data/ribll2022/"

MASS_U = 931.494
MASS_2H = MASS_U * 2.0135531980
MASS_4He = MASS_U * 4.0015060943
MASS_10BE = MASS_U * 10.0113403769
MASS_14C = MASS_U * 13.9999505089

PPAC_XZ = [-695.2, -454.2, -275.2]
PPAC_YZ = [-689.2, -448.2, -269.2]

def momentum_from_kinetic(mass, kinetic):
	return sqrt((2.0 * mass + kinetic) * kinetic)

def equations(vars, be_x, be_y, he_x, he_y, dx, dy, px, py, bek, hek, dk, pxindex, pyindex):
	tx, ty, ck = vars
	# calculate 10Be
	be_direction = np.array([be_x-tx, be_y-ty, 100.0])
	be_direction = be_direction / np.linalg.norm(be_direction)
	be_momentum = momentum_from_kinetic(MASS_10BE, bek)
	bep = be_direction * be_momentum
	# calculate 4He
	he_direction = np.array([he_x-tx, he_y-ty, 100.0])
	he_direction = he_direction / np.linalg.norm(he_direction)
	he_momentum = momentum_from_kinetic(MASS_4He, hek)
	hep = he_direction * he_momentum
	# calculate 2H
	d_direction = np.array([dx-tx, dy-ty, 135.0])
	d_direction = d_direction / np.linalg.norm(d_direction)
	d_momentum = momentum_from_kinetic(MASS_2H, dk)
	dp = d_direction * d_momentum
	# calculate 14C
	c_direction = np.array([tx-px, (ty-py)*PPAC_XZ[pxindex]/PPAC_YZ[pyindex], -PPAC_XZ[pxindex]])
	c_direction = c_direction / np.linalg.norm(c_direction)
	c_momentum = momentum_from_kinetic(MASS_14C, ck)
	cp = c_direction * c_momentum
	# total momentum
	total_p = bep + hep + dp - cp
	return [total_p[0], total_p[1], total_p[2]]

def create_equations(be_x, be_y, he_x, he_y, dx, dy, px, py, bek, hek, dk, pxindex, pyindex):
	return partial(
		equations,
		be_x=be_x, be_y=be_y, he_x=he_x, he_y=he_y, dx=dx, dy=dy, px=px, py=py,
		bek=bek, hek=hek, dk=dk, pxindex=pxindex, pyindex=pyindex
	)

def main():
	input_file = GENERATE_PATH + "simulate/full-resolution-g.root"
	output_file = GENERATE_PATH + "simulate/single-ppac-generate-fsolve-sim.root"

	# Open the input ROOT file
	input_root_file = ROOT.TFile.Open(input_file, "READ")
	if not input_root_file:
		print(f"Error: Could not open input file {input_file}")
		return

	# Get the input TTree
	input_tree = input_root_file.Get("tree")
	if not input_tree:
		print(f"Error: Could not find tree 'tree' in file {input_file}")
		input_root_file.Close()
		return

	tx = array('d', [0.])
	ty = array('d', [0.])
	bex = array('d', [0.])
	bey = array('d', [0.])
	he_x = array('d', [0.])
	he_y = array('d', [0.])
	dx = array('d', [0.])
	dy = array('d', [0.])
	bek = array('d', [0.])
	hek = array('d', [0.])
	dk = array('d', [0.])
	ck = array('d', [0.])
	ppac_xflag = array('H', [0])
	ppac_yflag = array('H', [0])
	ppac_x = array('d', 3*[0.])
	ppac_y = array('d', 3*[0.])
	input_tree.SetBranchAddress('tx', tx)
	input_tree.SetBranchAddress('ty', ty)
	input_tree.SetBranchAddress('bex', bex)
	input_tree.SetBranchAddress('bey', bey)
	input_tree.SetBranchAddress('hex', he_x)
	input_tree.SetBranchAddress('hey', he_y)
	input_tree.SetBranchAddress('dx', dx)
	input_tree.SetBranchAddress('dy', dy)
	input_tree.SetBranchAddress('bek', bek)
	input_tree.SetBranchAddress('hek', hek)
	input_tree.SetBranchAddress('dk', dk)
	input_tree.SetBranchAddress('ck', ck)
	input_tree.SetBranchAddress('ppac_xflag', ppac_xflag)
	input_tree.SetBranchAddress('ppac_yflag', ppac_yflag)
	input_tree.SetBranchAddress('ppac_x', ppac_x)
	input_tree.SetBranchAddress('ppac_y', ppac_y)

	# Create the output ROOT file
	output_root_file = ROOT.TFile.Open(output_file, "RECREATE")
	if not output_root_file:
		print(f"Error: Could not create output file {output_file}")
		input_root_file.Close()
		return

	# Create the output TTree
	output_tree = ROOT.TTree("tree", "Solve single PPAC tracking equations with scipy")
	# output data
	ftx = array('d', 9*[0.])
	fty = array('d', 9*[0.])
	fck = array('d', 9*[0.])
	# Define branches
	output_tree.Branch("tx", tx, "tx/D")
	output_tree.Branch("ty", ty, "ty/D")
	output_tree.Branch("ftx", ftx, "ftx[9]/D")
	output_tree.Branch("fty", fty, "fty[9]/D")
	output_tree.Branch("ppac_xflag", ppac_xflag, "pxflag/s")
	output_tree.Branch("ppac_yflag", ppac_yflag, "pyflag/s")
	output_tree.Branch("bek", bek, "bek/D")
	output_tree.Branch("hek", hek, "hek/D")
	output_tree.Branch("dk", dk, "dk/D")
	output_tree.Branch("ck", fck, "ck[9]/D")


	# Fill the tree with some example data
	n_entries = input_tree.GetEntries()
	for i in tqdm(range(n_entries), desc="Processing"):
		input_tree.GetEntry(i)
		for i in range(9):
			fck[i] = ck[0]
			ftx[i] = tx[0]
			fty[i] = ty[0]

		# if (ppac_xflag[0] & ppac_yflag[0] & 1) == 1:
		# 	equations = create_equations(
		# 		bex[0], bey[0], he_x[0], he_y[0],
		# 		dx[0], dy[0], ppac_x[0], ppac_y[0],
		# 		bek[0], hek[0], dk[0], 0, 0
		# 	)
		# 	guess = [tx[0], ty[0], ck[0]]
		# 	solution = fsolve(equations, guess)
		# 	ftx[0] = solution[0]
		# 	fty[0] = solution[1]
		# 	fck[0] = solution[2]


		for ix in range(3):
			if (ppac_xflag[0] & (1 << ix)) == 0:
				continue
			for iy in range(3):
				if (ppac_yflag[0] & (1 << iy)) == 0:
					continue
				equations = create_equations(
					bex[0], bey[0], he_x[0], he_y[0],
					dx[0], dy[0], ppac_x[ix], ppac_y[iy],
					bek[0], hek[0], dk[0], ix, iy
				)
				guess = [tx[0], ty[0], ck[0]]
				solution = fsolve(equations, guess)
				ftx[ix*3+iy] = solution[0]
				fty[ix*3+iy] = solution[1]
				fck[ix*3+iy] = solution[2]

		output_tree.Fill()

	# Write the output TTree to the output file
	output_tree.Write()

	# Close the files
	input_root_file.Close()
	output_root_file.Close()

if __name__ == "__main__":
	main()
