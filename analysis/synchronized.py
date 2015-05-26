"""
Tests for synchronization between cells
Copyright (C) 2012 Ahmet Ay, Jack Holland, Adriana Sperlea

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import math
import shared

def main():
	# check the given arguments
	if len(sys.argv) < 3:
		usage()
	elif len(sys.argv) == 4:
		if sys.argv[1] == "-c" or sys.argv[1] == "--no-color":
			shared.terminalRed = ""
			shared.terminalReset = ""
			filename = sys.argv[2]
			startTime = sys.argv[3]
		else:
			usage()
	else:
		filename = sys.argv[1]
		startTime = sys.argv[2]
	
	# open the input file and ensure 'startTime' is a number
	f = shared.openFile(filename, "r")
	startTime = shared.toFlo(startTime)
	startIndex = 0
	
	# extract the width and height from the concentrations file
	width, height = shared.widthAndHeight(f.readline().split(), filename)
	
	cellCount = width * height
	sums = [0 for i in range(cellCount)]
	cells = []
	
	# for every line in the file, get the sum of the concentration levels as well as each individual value cast as a float, and calculate the starting index based off of the starting time
	for line in f:
		data = line.split()
		curdata = []
		
		if data[0] < startTime:
			startIndex += 1
		
		for cell in range(cellCount):
			sums[cell] += shared.toFlo(data[cell + 1])
			curdata.append(shared.toFlo(data[cell + 1]))
		
		cells.append(curdata)

	# ensure there was data to retrieve
	if len(cells) == 0:
		print shared.terminalRed + "Couldn't get any cell data! Make sure '" + filename + "' is properly formatted. Exit status 3.", shared.terminalReset
		exit(3)
	
	# calculate the mean of each cell
	means = []
	for cell in range(cellCount):
		means.append(sums[cell] / len(cells))
	
	# calculate the total average score
	avgscore = 0.0
	for cell in range(1, cellCount):
		numerator = 0
		sqr1 = 0
		sqr2 = 0
		for i in range(startIndex, len(cells)):
			tstep = cells[i]
			xi = tstep[0]
			yi = tstep[cell]
			numerator += ((xi - means[0]) * (yi - means[cell]))
			sqr1 += (xi - means[0]) ** 2
			sqr2 += (yi - means[cell]) ** 2
		
		sqr1 = math.sqrt(sqr1)
		sqr2 = math.sqrt(sqr2)
		if sqr1 == 0 or sqr2 == 0:
			r = 1
		else:
			r = numerator / (sqr1 * sqr2)
		avgscore += r
	
	# print the synchronization score
	print round(avgscore / (cellCount - 1), 10)
	f.close()

# print usage information
def usage():
	print "Usage: python synchronized.py [-c|--no-color] <file with concentration levels> <time point to start at>"
	exit(0)

main()
