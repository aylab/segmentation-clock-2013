"""
File comparator that prints how many unique and common lines there are between two files
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
import os
import shared

def main():
	# check the given arguments
	if len(sys.argv) < 3:
		usage()
	elif len(sys.argv) == 4:
		if sys.argv[1] == "-c" or sys.argv[1] == "--no-color":
			shared.terminalRed = ""
			shared.terminalReset = ""
			directory = sys.argv[2]
			filename = sys.argv[3]
		else:
			usage()
	else:
		directory = sys.argv[1]
		filename = sys.argv[2]
	
	# ensure the directory exists and open the output file
	directory = shared.ensureDir(directory)
	ofile = shared.openFile(filename, "w")
	
	fsize = 0
	filename = directory + "/run0.txt"
	f = shared.openFile(filename)
	width, height = shared.widthAndHeight(f.readline().split(), filename)
	
	for line in f:
		data = line.split()
		fsize = shared.toFlo(data[0])
	avg = [[-1 for i in range(w * h)] for j in range(int(fsize * 10) + 50)]
	end = 0
	f.close()

	for run in range(0, runs):
		filename = directory + "/run" + str(run)+ ".txt"
		f = open(filename)
		width, height = shared.widthAndHeight(f.readline().split(), filename)

		for line in f:
			lineList = line.split(",")
			temptime = shared.toFlo(lineList[0])
			place = int(temptime * 10)
			if float(place + 1) - (temptime * 10) < (temptime * 10 - float(place)):
				index = place + 1
			else:
				index = place
			if (index > end):
				end = index
			for cell in range(0, w * h):
				if (avg[index][cell] == -1):
					avg[index][cell] = shared.toFlo(lineList[cell + 1])
				else:
					avg[index][cell] += shared.toFlo(lineList[cell + 1])
	
	ofile.write(str(width) + " " + str(height) + "\n")
	
	for cell in range(0, width * height):
		i = 0
		minend = end
		while (i < end):
			if (avg[i + 1][cell] == -1):
				x1 = i
				y1 = avg[i][cell]
				y2 = 0
				x2 = -1
				for k in range(i + 2, end + 1):
					if (avg[k][cell] != -1):
						y2 = avg[k][cell]
						x2 = k
						break
				if (x2 != -1):
					m = (y2 - y1) / (x2 - x1)
					for k in range (i + 1, x2):
						avg[k][cell] = y1 + m * h * (k - i)
					i = x2
				else:
					end = i
					if (end < minend):
						minend = end
					break
			else:
				i += 1 	
	
	end = minend
	for i in range(1, end + 1):
		list = avg[i]
		ofile.write(str(float(i) / 10) + " ")
		
		for cell in range(0, w * h):
			ofile.write(str(list[cell] / runs) + " ")
		
		ofile.write("\n")
		index += 1
	ofile.close()

# print usage information
def usage():
	print "Usage: python average-runs.py [-c|--no-color] <input directory with run#.txt> <output filename>"
	exit(0)

main()

