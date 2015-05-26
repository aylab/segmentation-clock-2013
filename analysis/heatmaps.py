import numpy as np
from scipy import *
from pylab import *
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colorbar
import sys
import shared
import os

def ensure_dir(dir):
	if not os.path.exists(dir):
		os.makedirs(dir)

def main():
	filename = ""
	if len(sys.argv) == 8:
		directory = sys.argv[1] # directory where the mutation is
		first_run = shared.toInt(sys.argv[2])  
		last_run = shared.toInt(sys.argv[3])
		feature = sys.argv[4]
		picture = sys.argv[5]
		start = shared.toInt(sys.argv[6])
		end = shared.toInt(sys.argv[7])
	else:
		print("ofeatures.py requires 7 parameters --", len(sys.argv) - 1, "given.")
		exit(1)

	bins = 10
	miny = 999999
	maxy = 0
	x = []

	for r in range(first_run, last_run):
		filename = directory + "/run" + str(r) + "/" + feature

		f = file(filename, "r")

		index = 0
		for line in f:
			mylist = line.split()
			for i in range(len(mylist)):
				mylist[i] = float(mylist[i])
			num = min(mylist)
			if (num < miny):
				miny = num
			num = max(mylist)
			if (num > maxy):
				maxy = num
			if (r == first_run):
				x.append(mylist)
			else:
				x[index] = x[index] + mylist
			index += 1

	z = []
	minz = inf
	maxz = 0
	for line in x:
		newline = [0] * (bins + 1)
		for el in line:
			bin = ((el - start) * 100) / (end - start)
			newline[int(bin / float(100 / bins))] += 1
		if (max(newline) > maxz):
			maxz = max(newline)
		if (min(newline) < minz):
			minz = min(newline)
		z.append(newline)

	hold(True)

	cbar = plt.colorbar(imshow(z, aspect='auto', extent=[start, end , 0, 16], cmap="hot"))
	xlabel(feature[0:1].upper() + feature[1:-4])
	ylabel('Cells')
	cbar.ax.set_yticks([0, 1])
	cbar.ax.set_yticklabels(['Low', 'High'])

	savefig(picture)
	
main()
