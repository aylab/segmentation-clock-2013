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
import shared

def main():
	# check the given arguments
	if len(sys.argv) < 3:
		usage()
	elif len(sys.argv) == 4:
		if sys.argv[1] == "-c" or sys.argv[1] == "--no-color":
			shared.terminalRed = ""
			shared.terminalReset = ""
			filename1 = sys.argv[2]
			filename2 = sys.argv[3]
		else:
			usage()
	else:
		filename1 = sys.argv[1]
		filename2 = sys.argv[2]
	
	# open input files
	f1 = shared.openFile(filename1, "r")
	f2 = shared.openFile(filename2, "r")
	
	# make a list for each file containing each line of it as an element
	list1 = []
	list2 = []
	for line in f1:
		list1.append(line)
	for line in f2:
		list2.append(line)
	f1.close()
	f2.close()
	
	size = len(mylist1)
	unique1 = 0
	unique2 = 0

	# add lines that are unique and in common
	for line in list1:
		if line not in list2:
			unique1 += 1
			continue
	
	for line in list2:
		if line not in list1:
			unique2 += 1
			continue
	
	# print the results
	print filename1, "has", unique1, "unique lines."
	print filename2, "has", unique2, "unique lines."
	print "There are", size - unique1, "lines that appear in both files."

# print usage information
def usage():
	print "Usage: python compare-files.py [-c|--no-color] <file 1> <file 2>"
	exit(0)

main()

