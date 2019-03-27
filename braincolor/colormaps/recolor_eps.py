#! /usr/bin/env python
#                                                      22 October 2010
#		Andrew J. Worth
#		andy@neuromorphometrics.com
#
#		Neuromorphometrics, Inc.
#		22 Westminster Street
#		Somerville, MA  02144-1630  USA
#
#		http://neuromorphometrics.com
#
#*********************************************************************
#
#  (c) Copyright 2010 Neuromorphometrics, Inc.   All rights reserved
#
#*********************************************************************

import sys
import string
import os
import re
import array
import math
import glob
import datetime
from optparse import OptionParser
from xml.dom import minidom

usage = "%prog [options] inFile.eps outFile.eps\n\
\n\
This program identifies and recolors items in an eps file.\n\
\n\
Examples:\n\
  %prog -p --color 'RRR GGG BBB' inFile.eps preparedFile.eps \n\
  %prog --mapFile labelFileName.xml preparedFile.eps finalFile.eps"

parser = OptionParser(usage)

parser.add_option("-p", "--prepare",
		  action="store_true", dest="doPrepare", default=False,
                  help="prepare the file for coloring "
                  "[default: %default]")
parser.add_option("-c", "--color",
                  action="store", type="string", dest="findMeColor",
                  default='255 255 0',
                  help="color used to identify colored items "
                       "[default: %default]")
parser.add_option("-m", "--mapFile",
                  action="store", type="string", dest="mapFileName",
                  default='parcLabels.xml',
                  help="Name to color mapping file "
                       "[default: %default]")
parser.add_option("-o", "--overrideCount",
                  action="store", type="int", dest="overrideCount",
                  default='0',
                  help="Override number of colors found (for testing) "
                       "[default: %default]")

parser.add_option("-v", action="store_true", dest="verbose", default=True,
                  help="say what is going on "
                  "[default: %default]")
parser.add_option("-q", action="store_false", dest="verbose",
                  help="don't say what is going on")
(options, args) = parser.parse_args()

if len(args) != 2:
	parser.error("2 arguments are required")

# Get required arguments
inFileName  = args[0]
outFileName = args[1]

if options.verbose:
	print "inFileName  = %s"     % inFileName
	print "outFileName = %s"     % outFileName
	print "color       = \"%s\"" % options.findMeColor
	print "mapFile     = %s"     % options.mapFileName

# This program takes an eps file and can do one of two things:
# 1) Prepare an eps file to be recolored by identifying the items in
#    the eps file that are colored, or
# 2) Recolor a prepared eps file given a color mapping file.
# 
# To prepare an eps file, first "eps2eps" is run on it to remove the 
# cruft.  Then the file is written out for each colored item found with
# that colored item set to a given color (yellow by default), opened
# for viewing, and then asks for an identifier for that colored item.
# Once all desired colored items are identified, the eps file is written
# out with a comment indicating the identifier for each color.
# 
# To recolor the prepared eps file, a mapping file is read that gives
# the color for each of the identifiers that might be found in the
# eps file comments.  The mapping file is XML like this:
# 
# 	<?xml version="1.0"?>
# 	<LabelList>
#	<Label>
#	  <Name>Unlabeled</Name>
#	  <Number>1</Number>
#	  <RGBColor>0 0 0</RGBColor>
#	</Label>
#	...
#	<Label>
#	  <Name>Vitamin E Tablet</Name>
#	  <Number>74</Number>
#	  <RGBColor>236 217 151</RGBColor>
#	</Label>
#	</LabelList>
#
# The <Number> tag is not used here, but it is used by NVM!
#
# After reading the mapping file, the eps file is read and the colors
# in the mapping file are substituted for the original colors for each
# of the identifiers, and then the final recolored eps file is written
# out.

labelList = [ ]
labelNums = [ ]
labelColors = [ ]

if os.access(options.mapFileName,os.F_OK): # then file exists
	# Read color map file
	xmldoc = minidom.parse(options.mapFileName)
	xmlLabelList = xmldoc.getElementsByTagName("Label")
	for xmlLabel in xmlLabelList:
		name = xmlLabel.getElementsByTagName("Name")
		labelList.append(name[0].firstChild.data)
		uNumber = xmlLabel.getElementsByTagName("Number")
		number = int(uNumber[0].firstChild.data)
		labelNums.append(number)
		uColor = xmlLabel.getElementsByTagName("RGBColor")
		t = uColor[0].firstChild.data.split()
		color = t[0] + ' ' + t[1] + ' ' + t[2]
		labelColors.append(color)

	#for ii in range(len(labelList)):
	#	name = labelList[ii]
	#	color = labelColors[ii]
	#	print '\'' + color + '\' \'' + name + '\''


# --------------------------------------------------------------------
#
if options.doPrepare:
	# run eps2eps on the input file:
	cruftlessFile = "nocruft." + inFileName
	print "\nPreparing file \"%s\" ...\n" % inFileName
	s = 'eps2eps %s %s' % (inFileName, cruftlessFile)
	if options.verbose:
		print s
	os.system(s)

	# read the resulting cruftless input file:
	ff = open(cruftlessFile,'r')
	contents = ff.read()
	ff.close()
	count = contents.count(" rG\n")
	if options.verbose:
		print "Found %d colored items" % count
	contentLines = contents.split('\n')

	if options.overrideCount != 0:
		count = options.overrideCount
	regionNameList = [ ]
	for ii in range(count):
		# open temporary file for writing
		tempFileName = outFileName.replace('.eps','_'+str(ii)+'.eps')
		tf = open(tempFileName,'w')
		foundCount = 0
		for line in contentLines:
		        h = re.compile('.* rG$')
		        hS = h.search(line)
			# if this is a color line
			#print line
		        if hS:
				#print "ii= ", ii, ' found line: ', line
				# if this is the color to change change it
				if foundCount == ii:
					s = options.findMeColor + ' rG\n'
					#print 'changed color: ', s
					tf.write(s)
				else: # not THE color, print it out
					tf.write(line+'\n')
				foundCount = foundCount + 1
			else: # something other than color, print it out
				tf.write(line+'\n')
		tf.close()
		s = 'open ' + tempFileName
		os.system(s)

		go = True
		while go: # then a label file was read
			thisName = raw_input("Enter label for region: ")
			if options.verbose:
				print "You said, ", thisName
			pickList = [ ]
			print "   0) (re-enter text) "
			print "   1) %s" % thisName
			pickList.append(thisName)
			pickNum=2
			if len(labelList)>0: 
				for label in labelList:
	        			h = re.compile('.* '+thisName+' .*')
	        			hS = h.search(label)
	        			if hS: # if this is a color comment line
						print "   %d) %s" % \
							(pickNum, label)
						pickList.append(label)
						pickNum = pickNum + 1
			choice = raw_input("Which one do you want to use?  ")
			theChoice = int(choice)
			if theChoice == 0:
				go = True # stay in while loop
			else:
				go = False # drop out
				thisName = pickList[theChoice-1]
		print "Using: \"%s\"" % thisName
		regionNameList.append(thisName)

	# open output file for writing
	of = open(outFileName,'w')

	# Now print out the new file with the names as comments
	foundCount = 0
	for line in contentLines:
	        h = re.compile('.* rG$')
	        hS = h.search(line)
	        if hS: # if this is a color line
			if foundCount < count:
				s = '% recoloreps ' \
			   	+ regionNameList[foundCount]+'\n'
				of.write(s)
				foundCount = foundCount + 1
			of.write(line+'\n')
		else: # something other than color, print it out
			of.write(line+'\n')

	of.close()


# --------------------------------------------------------------------
#
else:
	print "\nColoring file \"%s\" ...\n" % inFileName

	# Read input file and write ouput, changing the colors after
	#  finding '% recoloreps LABELSTRING' comments
	ff = open(inFileName,'r')
	contents = ff.read()
	ff.close()
	contentLines = contents.split('\n')

	# open output file for writing
	of = open(outFileName,'w')

	skipNextLine = False
	for line in contentLines:
		if skipNextLine == False:
	        	h = re.compile('% recoloreps .*')
	        	hS = h.search(line)
	        	if hS: # if this is a color comment line
				print 'Found color comment: ', line
				of.write(line+'\n')
				toFind = line[13:]
				if toFind in labelList:
					print 'looking for color for ', toFind
					index = labelList.index(toFind)
					print 'index is ', index
					color = labelColors[index]
					print 'writing color: ', color
					of.write(color+' rG\n')
					skipNextLine = True
				else:
					print toFind + ' is not in labelList\n'
			else: # something other than color, print it out
				of.write(line+'\n')
		else:
			print 'skipped actual color line\n'
			skipNextLine = False

	of.close()

# --------------------------------------------------------------------
#
if options.verbose:
	print "All done, bye."

sys.exit()
