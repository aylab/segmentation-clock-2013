"""
Creates snapshots of cell tissue colored using concentration levels
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

import Image, ImageDraw
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
			directory = sys.argv[3]
		else:
			usage()
	else:
		filename = sys.argv[1]
		directory = sys.argv[2]
	
	# ensure the given directory exists
	directory = shared.ensureDir(directory)
	if (directory[-1] != '/'):
		directory = directory + '/'
	
	minCon, maxCon = findMinMax(filename) # find the minimum and maximum values of the concentrations
	width, height, parsedFile = readData(filename) # find the width and height of the cell tissue and parse the concentrations file into an array
	edge, size = findSizes(width, height) # configure the hexagon edge and window size based on the grid size
	
	index = 0
	for line in parsedFile:
		if (index % 10 == 0):
			plotHexagons(directory, size, index, line[1:], minCon, maxCon, edge, width, height)
		index += 1

# class that contains the information required to draw a hexagon
class hexagon():
	def __init__(self, center, edge, color): 
		# creates a hexagon with the center at self.center and the edge of size self.edge
		self.center = center
		self.edge = edge
		self.color = color
	
	def draw(self):
		# draws the hexagon centered at self.center
		half = float(self.edge / 2.0)
		rad = float(self.edge * 1.73 / 2.0)
		i = self.center[0]
		j = self.center[1]
		lines = [(i, j - self.edge), (i + rad, j - half), (i + rad, j + half), (i, j + self.edge), (i - rad, j + half), (i - rad, j - half)]
		return lines

# finds the minimum and maximum concentrations given a properly formatted file
def findMinMax(filename):
	f = shared.openFile(filename, "r")
	minCon = float("inf")
	maxCon = 0
	first = True
	for line in f:
		if not first:
			lineNums = line.split()
			for i in range(len(lineNums)):
				lineNums[i] = shared.toFlo(lineNums[i])
			num = min(lineNums[1:])
			if (num < minCon):
				minCon = num
			num = max(lineNums[1:])
			if (num > maxCon):
				maxCon = num
		first = False
	f.close()
	return minCon, maxCon

# extracts the width and height from the concentrations file and parses each line into an array
def readData(filename):
	f = shared.openFile(filename, "r")
	parsedFile = []
	index = 0
	width = height = 0
	for line in f:
		line = line.strip()
		if (index == 0):
			width, height = shared.widthAndHeight(line.split(" "), filename)
			if width < 4 or height < 4 or (width % 2 != 0) or (height % 2 != 0):
				print shared.terminalRed + "The size of the tissue must be at least 4x4 and its width and height must be even numbers. Exit status 2.", shared.terminalReset
				exit(2)
		else:
			aux = line.split("\t")
			parsedFile.append(aux)
		index += 1
	f.close()
	return width, height, parsedFile

# finds the edge length and size of a hexagon, properly formatted for 300 pixels per inch (journal specifications)
def findSizes(width, height):
	edge = (960 / 1.73) / float(width)
	size = (960, int(2 * edge * (height / 2) + edge * (height / 2)))
	return edge, size

# plot the hexagons for the given line of data
def plotHexagons(output, size, index, protLevels, minCon, maxCon, edge, width, height):
	center = (float(edge * 1.73 / 2.0), 0) # place the center of the first hexagon
	im = Image.new("RGB", size, (255, 255, 255))
	draw = ImageDraw.Draw(im)
	drawGrid(draw, center, edge, width, height, protLevels, minCon, maxCon)
	fix = fixDigits(index / 10, 4)
	im.save(output + fix + ".tiff", "TIFF")
	del im

# draws a grid of hexagons with the first hexagon's center at 'begin'
def drawGrid(draw, begin, edge, width, height, protLevels, minCon, maxCon):
	center = begin	
	rad = float(edge * 1.73 / 2.0)
	for i in range(height):
		# draw a half of the last cell
		if (i % 2 == 1):
			drawCell(draw, protLevels, i * width + (width - 1), edge, (0, center[1]), minCon, maxCon)
		for j in range(width):
			cell = i * width + j
			drawCell(draw, protLevels, cell, edge, center, minCon, maxCon)			
			center = (center[0] + 2.0 * rad, center[1]) # move center to the right by 2 * edge * sqrt(3) / 2
			
		if (i % 2 == 0):
			center = (begin[0] + rad, center[1] + edge + (edge / 2.0))
		else:
			center = (begin[0], center[1] + edge + (edge / 2.0))
	
	# draw last row of half hexagons from the first row
	for j in range(width):
		cell = j;
		drawCell(draw, protLevels, cell, edge, center, minCon, maxCon)
		
		center = (center[0] + 2.0 * rad, center[1])

# finds the edge of a hexagon given its size, width, and height
def findEdge(size, width, height):
	OX = size[1] / (3 * (height / 2) + 0.5)
	OY = float(size[0] / (width * 1.73))
	return min(OX, OY)

# converts a triplet of decimal integers to a triplet of hexadecimal integers
hexDigits = "0123456789abcdef" # hexadecimal digits
def rgb(triplet):
    triplet = triplet.lower()
    return (hexDigits.index(triplet[0]) * 16 + hexDigits.index(triplet[1]),
            hexDigits.index(triplet[2]) * 16 + hexDigits.index(triplet[3]),
            hexDigits.index(triplet[4]) * 16 + hexDigits.index(triplet[5]))

# draws an individual hexagon
redShades = [rgb("FF8E99"), rgb("FF70B7"), rgb("F76F87"), rgb("FF5084"), rgb("D93657"), rgb("F73E5f"), rgb("B32D45"), rgb("8D2439"), rgb("9B001C"), rgb("750017")] # the gradients of red to use
def drawCell(draw, protLevels, cell, edge, center, minCon, maxCon):
	colorindex = ((float(protLevels[cell]) - minCon) * 100) / (maxCon - minCon)
	colorindex = int (colorindex/10)
	hex = hexagon(center, edge, rgb("fcc75e")) # border
	draw.polygon(hex.draw(), fill = hex.color)
	del hex
	
	hex = hexagon(center, edge - 2, redShades[colorindex - 1]) # actual hexagon
	draw.polygon(hex.draw(), fill = hex.color)
	del hex

# fix the digits of a number so that they are padded with 0's up to 'n'
def fixDigits(index, n):
	zero = ""
	fix = str(index)
	while (len(fix) + len(zero) < n):
		zero += "0"
	fix = zero + fix
	return fix

# print usage information
def usage():
	print "Usage: python tissue-snapshots.py [-c|--no-color] <file with concentration levels> <directory to store images>"
	exit(0)

main()

