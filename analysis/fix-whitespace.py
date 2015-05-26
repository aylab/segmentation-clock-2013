"""
Small utility to fix the way Xcode saves whitespace (converts spaces to tabs)
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
			filename = sys.argv[2]
			spaces = sys.argv[3]
		else:
			usage()
	else:
		filename = sys.argv[1]
		spaces = sys.argv[2]
	
	# open the input file and check to ensure 'spaces' is an integer
	f = shared.openFile(filename, "r")
	spaces = shared.toInt(spaces)
	
	# replace spaces with tabs
	ofile = ""
	for line in f:
		count = 0
		while line[:spaces] == " " * spaces:
			line = line[spaces:]
			count += 1
		ofile = ofile + "\t" * count + line
	
	f.close()
	newfile = shared.openFile(filename, "w")
	newfile.write(ofile)

# print usage information
def usage():
	print "Usage: python fix-whitespace.py [-c|--no-color] <file to fix> <spaces per tab>"
	exit(0)

main()

