#!/usr/bin/env python3

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import os
import subprocess
from subprocess import call
import pysam
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import Bio
from Bio.Seq import Seq
import argparse
import sys

#############################
## DEFAULT ARGUMENT VALUES ##
#############################
SeqLogo=False
PingPong=False
Distribution=False
Unique=False
Unstranded=False
MinLength=17
MaxLength=35

###########################
## USER ARGUMENT PARSING ##
###########################
parser = argparse.ArgumentParser(description='Read arguments')
parser.add_argument('-s', '--sample', type=str, help='sample name')
parser.add_argument('-i', '--infile', type=str, help='input bam file')
parser.add_argument('-m', '--minlength', type=int, help='minimum sRNA length')
parser.add_argument('-M', '--maxlength', type=int, help='maximum sRNA length')
parser.add_argument('--seqlogo', help='plot seqlogo', action='store_true')
parser.add_argument('--pingpong', help='plot Ping-Pong overlap signature', action='store_true')
parser.add_argument('--distribution', help='plot sRNA length distribution', action='store_true')
parser.add_argument('--unique', help='plot unique reads only (i.e. collapse identical reads)', action='store_true')
parser.add_argument('--unstranded', help='ignore strand information when plotting distribution', action='store_true')
args = parser.parse_args()
# sample name parsing
Sample = args.sample
if Sample is not None:
	print('Sample name is ' + Sample)
else:
	print('ERROR: no sample name (-s) specified')
	sys.exit(0)
# input bam file parsing
InFile = args.infile
if InFile is not None:
	if InFile.endswith('.bam'):
		print('Input file is ' + InFile)
	else:
		print('Input file ' + InFile + ' does not have normal bam filename ending - please check file format')
		sys.exit(0)
else:
	print('ERROR: no input bam file (-i) specified')
	sys.exit(0)
# minimum sRNA length parsing
if args.minlength is None:
	print('Using default minimum sRNA length (' + str(MinLength) + ')')
else:
	MinLength = args.minlength
	if MinLength > 0:
		print('Minimum sRNA length set (' + str(MinLength) + ')')
	else:
		print('ERROR: minimum sRNA length (-m) specified but no length given')
		sys.exit(0)
# maximum sRNA length parsing
if args.maxlength is None:
	print('Using default maximum sRNA length (' + str(MaxLength) + ')')
else:
	MaxLength = args.maxlength
	if MaxLength > 0:
		if MaxLength >= MinLength:
			print('Maximum sRNA length set (' + str(MaxLength) + ')')
		else:
			print('ERROR: maximum sRNA length (-M) must be higher than minimum sRNA length (-m)')
			sys.exit(0)
	else:
		print('ERROR: maximum sRNA length (-M) specified but no length given')
		sys.exit(0)
# seqlogo option parsing
if args.seqlogo:
	SeqLogo=True
	print('Plotting seqlogo')
else:
	print('Not plotting seqlogo')
# pingpong option parsing
if args.pingpong:
	PingPong=True
	print('Plotting Ping-Pong signature')
else:
	print('Not plotting Ping-Pong signature')
# distribution option parsing
if args.distribution:
	Distribution=True
	print('Plotting sRNA length distribution')
else:
	print('Not plotting sRNA length distribution')
# read collapsing option parsing
if args.unique:
	Unique=True
	print('Collapsing reads to unique reads only')
else:
	print('Not collapsing reads')
# stranded option parsing
if args.unstranded:
	Unstranded=True
	print('Ignoring strand information')
else:
	print('Using strand information')

######################################################################################
# function to parse a bam file, returns a dict of 5' base counts for each seq length #
######################################################################################
def bam_5prime(file,minlength=0,maxlength=1000,unique=False):
	# empty dicts for read sequences and read counts
	sense_seqs = []
	sense_counts = []
	antisense_seqs = []
	antisense_counts = []
	# read in bamfile
	bamfile = pysam.AlignmentFile(file,"rb")

	for read in bamfile:
		# convert line into string to make it splittable
		line = str(read)
		# split line on tabs
		linesplit = line.split("\t")
		# this filters out reads with flag 4 (i.e. unmapped reads)
		cigar = linesplit[5]
		if cigar.endswith('M'):
			seq = linesplit[9]
			if unique==True:
				count = 1
			else:
				count = int(linesplit[0].split('-')[1])
			# this is where the length filtering happens
			if minlength <= len(seq) <= maxlength:
				# this splits up sense and antisense reads
				if linesplit[1] == '16':
					# convert seq string to Biopython Seq object
					tempseq=Seq(seq)
					# set reverse-complement of Seq object as string
					tempseqrevcomp = tempseq.reverse_complement()
					antisense_seqs.append(tempseqrevcomp)
					antisense_counts.append(count)
				else:
					sense_seqs.append(seq)
					sense_counts.append(count)

	# go through all of the sequences, creating a list of all read lengths
	lengths = []
	for i in sense_seqs:
		if len(i) not in lengths:
			lengths.append(len(i))
	for i in antisense_seqs:
		if len(i) not in lengths:
			lengths.append(len(i))
	# sort the lengths list in ascending order
	lengths.sort()
	
	# for each base, the total number of the shortest reads starting with that base will be appended, followed by the next length, and the next length...
	A = []
	C = []
	G = []
	T = []
	N = []
	# go through the lengths, looking at the first base of each sense read of that length, and totalling them up
	for readlength in lengths:
		Acount=0
		Ccount=0
		Gcount=0
		Tcount=0
		Ncount=0
		for i in range(len(sense_seqs)):
			read = sense_seqs[i]
			if len(read) == readlength:
				firstbase = read[0]
				if unique==True:
					count = 1
				elif unique==False:
					count = sense_counts[i]
				if firstbase == 'A':
					Acount = Acount + count
				elif firstbase == 'C':
					Ccount = Ccount + count
				elif firstbase == 'G':
					Gcount = Gcount + count
				elif firstbase == 'T':
					Tcount = Tcount + count
				elif firstbase == 'N':
					Ncount = Ncount + count
		# go through the lengths, looking at the first base of each antisense read of that length, and totalling them up
		for i in range(len(antisense_seqs)):
			read = antisense_seqs[i]
			if len(read) == readlength:
				firstbase = read[0]
				if unique==True:
					count = 1
				elif unique==False:
					count = antisense_counts[i]
				if firstbase == 'A':
					Acount = Acount + count
				elif firstbase == 'C':
					Ccount = Ccount + count
				elif firstbase == 'G':
					Gcount = Gcount + count
				elif firstbase == 'T':
					Tcount = Tcount + count
				elif firstbase == 'N':
					Ncount = Ncount + count
		A.append(Acount)
		C.append(Ccount)
		G.append(Gcount)
		T.append(Tcount)
		N.append(Ncount)
		
	print('Bases counted')

	# format the dataframe
	formatted = {}
	formatted['Length'] = lengths
	formatted['A'] = A
	formatted['C'] = C
	formatted['G'] = G
	formatted['T'] = T
	formatted['N'] = N	
	basecounts = pd.DataFrame(formatted,columns=['Length','A','C','G','T','N'])
	return(basecounts)
	
	
#########################################################
# function to make a stacked barplot, returns a barplot #
#########################################################
def stacked_barplot(dataframe,samplename='Plot',xlabel='xlabel',ylabel='ylabel'):
	# Create the general blog and the "subplots" i.e. the bars
	f, ax1 = plt.subplots(1, figsize=(6.597222,3.298611))
	# Set the bar width
	bar_width = 0.75

	# positions of the left bar-boundaries
	bar_l = [i+1 for i in range(len(dataframe['Length']))]

	# positions of the x-axis ticks (center of the bars as bar labels)
	tick_pos = [i+(bar_width/2) for i in bar_l]

	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the N data
			dataframe['N'],
			# set the width
			width=bar_width,
			# with the label N
			label='N',
			# with alpha 1
			alpha=1,
			# with color
			color='#000000')

	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the A data
			dataframe['A'],
			# set the width
			width=bar_width,
			# with N on the bottom
			bottom=dataframe['N'],
			# with the label A
			label='A',
			# with alpha 1
			alpha=1,
			# with color
			color='#32CD32')

	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the C data
			dataframe['C'],
			# set the width
			width=bar_width,
			# with N and A on the bottom
			bottom=[i+j for i,j in zip(dataframe['A'],dataframe['N'])],
			# with the label C
			label='C',
			# with alpha 1
			alpha=1,
			# with color
			color='#4169E1')
		
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the G data
			dataframe['G'],
			# set the width
			width=bar_width,
			# with N and A and C on the bottom
			# NB if the number of chunks in each stack changes, "i+j+k for i,j,k in" has to be altered to match the number of chunks e.g. to do four chunks, it would be "i+j+k+l for i,j,k,l in"
			bottom=[i+j+k for i,j,k in zip(dataframe['C'],dataframe['A'],dataframe['N'])],
			# with the label G
			label='G',
			# with alpha 1
			alpha=1,
			# with color
			color='#FFD700')

	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the T data
			dataframe['T'],
			# set the width
			width=bar_width,
			# with N and A and C and G on the bottom
			# NB if the number of chunks in each stack changes, "i+j+k for i,j,k in" has to be altered to match the number of chunks e.g. to do four chunks, it would be "i+j+k+l for i,j,k,l in"
			bottom=[i+j+k+l for i,j,k,l in zip(dataframe['G'],dataframe['C'],dataframe['A'],dataframe['N'])],
			# with the label U
			label='U',
			# with alpha 1
			alpha=1,
			# with color
			color='#8B0000')

	# set the x ticks with names
	plt.xticks(tick_pos, dataframe[xlabel])

	# Set the label, legends and title
	ax1.set_ylabel(ylabel)
	ax1.set_xlabel(xlabel)
	#plt.legend(loc='upper right')
	#plt.title(samplename)

	# Set a buffer around the edge
	plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
	print('Plot created')

	# remove top and right spines
	ax1.spines['right'].set_visible(False)
	ax1.yaxis.set_ticks_position('left')
	ax1.spines['top'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	# remove tick marks from x axis while keeping labels
	ax1.tick_params(axis='x',which='both',length=0)
	
	# output plot to file
	outname = samplename + '_unstranded.pdf'
	outfile = PdfPages(outname)
	plt.savefig(outfile, format='pdf')
	outfile.close()
	# clear plot window to prevent overplotting
	plt.clf()
	return('done')

######################################################################################
# function to parse a bam file, returns a dict of 5' base counts for each seq length #
######################################################################################
def bam_5prime_stranded(file,minlength=0,maxlength=1000,unique=False):
	# empty dicts for read sequences and read counts
	sense_seqs = []
	sense_counts = []
	antisense_seqs = []
	antisense_counts = []
	# read in bamfile
	bamfile = pysam.AlignmentFile(file,"rb")

	for read in bamfile:
		# convert line into string to make it splittable
		line = str(read)
		# split line on tabs
		linesplit = line.split("\t")
		# this filters out reads with flag 4 (i.e. unmapped reads)
		cigar = linesplit[5]
		if cigar.endswith('M'):
			seq = linesplit[9]
			if unique==True:
				count = 1
			else:
				count = int(linesplit[0].split('-')[1])
			# this is where the length filtering happens
			if minlength <= len(seq) <= maxlength:
				# this splits up sense and antisense reads
				if linesplit[1] == '16':
					# convert seq string to Biopython Seq object
					tempseq=Seq(seq)
					# set reverse-complement of Seq object as string
					tempseqrevcomp = tempseq.reverse_complement()
					antisense_seqs.append(tempseqrevcomp)
					antisense_counts.append(count)
				else:
					sense_seqs.append(seq)
					sense_counts.append(count)

	# go through all of the sequences, creating a list of all read lengths
	lengths = []
	for i in sense_seqs:
		if len(i) not in lengths:
			lengths.append(len(i))
	for i in antisense_seqs:
		if len(i) not in lengths:
			lengths.append(len(i))
	# sort the lengths list in ascending order
	lengths.sort()
	
	# for each base, the total number of the shortest reads starting with that base will be appended, followed by the next length, and the next length...
	senseA = []
	senseC = []
	senseG = []
	senseT = []
	senseN = []
	antisenseA = []
	antisenseC = []
	antisenseG = []
	antisenseT = []
	antisenseN = []
	# go through the lengths, looking at the first base of each sense read of that length, and totalling them up
	for readlength in lengths:
		Acount=0
		Ccount=0
		Gcount=0
		Tcount=0
		Ncount=0
		for i in range(len(sense_seqs)):
			read = sense_seqs[i]
			if len(read) == readlength:
				firstbase = read[0]
				if unique==True:
					count = 1
				elif unique==False:
					count = sense_counts[i]
				if firstbase == 'A':
					Acount = Acount + count
				elif firstbase == 'C':
					Ccount = Ccount + count
				elif firstbase == 'G':
					Gcount = Gcount + count
				elif firstbase == 'T':
					Tcount = Tcount + count
				elif firstbase == 'N':
					Ncount = Ncount + count
		senseA.append(Acount)
		senseC.append(Ccount)
		senseG.append(Gcount)
		senseT.append(Tcount)
		senseN.append(Ncount)
	# go through the lengths, looking at the first base of each antisense read of that length, and totalling them up
	for readlength in lengths:
		Acount=0
		Ccount=0
		Gcount=0
		Tcount=0
		Ncount=0
		for i in range(len(antisense_seqs)):
			read = antisense_seqs[i]
			if len(read) == readlength:
				firstbase = read[0]
				# NB the counts here are subtracted rather than added to make the plotting on the bottom of the plot work
				if unique==True:
					count = 1
				elif unique==False:
					count = antisense_counts[i]
				if firstbase == 'A':
					Acount = Acount - count
				elif firstbase == 'C':
					Ccount = Ccount - count
				elif firstbase == 'G':
					Gcount = Gcount - count
				elif firstbase == 'T':
					Tcount = Tcount - count
				elif firstbase == 'N':
					Ncount = Ncount - count
		antisenseA.append(Acount)
		antisenseC.append(Ccount)
		antisenseG.append(Gcount)
		antisenseT.append(Tcount)
		antisenseN.append(Ncount)
		
	print('Bases counted')
	
	# format the dataframe
	formatted = {}
	formatted['Length'] = lengths
	formatted['senseA'] = senseA
	formatted['senseC'] = senseC
	formatted['senseG'] = senseG
	formatted['senseT'] = senseT
	formatted['senseN'] = senseN	
	formatted['antisenseA'] = antisenseA
	formatted['antisenseC'] = antisenseC
	formatted['antisenseG'] = antisenseG
	formatted['antisenseT'] = antisenseT
	formatted['antisenseN'] = antisenseN	
	basecounts = pd.DataFrame(formatted,columns=['Length','senseA','senseC','senseG','senseT','senseN','antisenseA','antisenseC','antisenseG','antisenseT','antisenseN'])
	return(basecounts)
	
####################################################################################################################
# function to make a stacked barplot of stranded data, returns a barplot with sense on top and antisense on bottom #
####################################################################################################################
def stacked_barplot_stranded(dataframe,samplename='Plot',xlabel='xlabel',ylabel='ylabel'):
	# Create the general blog and the "subplots" i.e. the bars
	f, ax1 = plt.subplots(1, figsize=(10,5))
	# Set the bar width
	bar_width = 0.75
	# positions of the left bar-boundaries
	bar_l = [i+1 for i in range(len(dataframe['Length']))]
	# positions of the x-axis ticks (center of the bars as bar labels)
	tick_pos = [i+(bar_width/2) for i in bar_l]
	# plot sense reads (top of plot)
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the sense N data
			dataframe['senseN'],
			# set the width
			width=bar_width,
			# with the label N
			label='N',
			# with alpha 1
			alpha=1,
			# with color
			color='#000000')
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the C data
			dataframe['senseA'],
			# set the width
			width=bar_width,
			# with N on the bottom
			bottom=dataframe['senseN'],
			# with the label A
			label='A',
			# with alpha 1
			alpha=1,
			# with color
			color='#32CD32')
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the C data
			dataframe['senseC'],
			# set the width
			width=bar_width,
			# with N and A on the bottom
			bottom=[i+j for i,j in zip(dataframe['senseA'],dataframe['senseN'])],
			# with the label C
			label='C',
			# with alpha 1
			alpha=1,
			# with color
			color='#4169E1')	
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the G data
			dataframe['senseG'],
			# set the width
			width=bar_width,
			# with N and A and C on the bottom
			# NB if the number of chunks in each stack changes, "i+j+k for i,j,k in" has to be altered to match the number of chunks e.g. to do four chunks, it would be "i+j+k+l for i,j,k,l in"
			bottom=[i+j+k for i,j,k in zip(dataframe['senseC'],dataframe['senseA'],dataframe['senseN'])],
			# with the label G
			label='G',
			# with alpha 1
			alpha=1,
			# with color
			color='#FFD700')
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the T data
			dataframe['senseT'],
			# set the width
			width=bar_width,
			# with N and A and C and G on the bottom
			# NB if the number of chunks in each stack changes, "i+j+k for i,j,k in" has to be altered to match the number of chunks e.g. to do four chunks, it would be "i+j+k+l for i,j,k,l in"
			bottom=[i+j+k+l for i,j,k,l in zip(dataframe['senseG'],dataframe['senseC'],dataframe['senseA'],dataframe['senseN'])],
			# with the label U
			label='U',
			# with alpha 1
			alpha=1,
			# with color
			color='#8B0000')
	# plot antisense reads (bottom of plot)
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the antisense N data
			dataframe['antisenseN'],
			# set the width
			width=bar_width,
			# with the label N
			label='antisenseN',
			# with alpha 1
			alpha=1,
			# with color
			color='#000000')
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the antisense A data
			dataframe['antisenseA'],
			# set the width
			width=bar_width,
			# with N on the bottom
			bottom=dataframe['antisenseN'],
			# with the label A
			label='antisenseA',
			# with alpha 1
			alpha=1,
			# with color
			color='#32CD32')
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the antisense C data
			dataframe['antisenseC'],
			# set the width
			width=bar_width,
			# with N and A on the bottom
			bottom=[i+j for i,j in zip(dataframe['antisenseA'],dataframe['antisenseN'])],
			# with the label C
			label='antisenseC',
			# with alpha 1
			alpha=1,
			# with color
			color='#4169E1')	
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the antisense G data
			dataframe['antisenseG'],
			# set the width
			width=bar_width,
			# with N and A and C on the bottom
			# NB if the number of chunks in each stack changes, "i+j+k for i,j,k in" has to be altered to match the number of chunks e.g. to do four chunks, it would be "i+j+k+l for i,j,k,l in"
			bottom=[i+j+k for i,j,k in zip(dataframe['antisenseC'],dataframe['antisenseA'],dataframe['antisenseN'])],
			# with the label G
			label='antisenseG',
			# with alpha 1
			alpha=1,
			# with color
			color='#FFD700')
	# Create a bar plot, in position bar_1
	ax1.bar(bar_l,
			# using the antisense T data
			dataframe['antisenseT'],
			# set the width
			width=bar_width,
			# with N and A and C and G on the bottom
			# NB if the number of chunks in each stack changes, "i+j+k for i,j,k in" has to be altered to match the number of chunks e.g. to do four chunks, it would be "i+j+k+l for i,j,k,l in"
			bottom=[i+j+k+l for i,j,k,l in zip(dataframe['antisenseG'],dataframe['antisenseC'],dataframe['antisenseA'],dataframe['antisenseN'])],
			# with the label U
			label='U',
			# with alpha 1
			alpha=1,
			# with color
			color='#8B0000')		
	# set the x ticks with names
	plt.xticks(tick_pos, dataframe[xlabel])
	# Set the label, legends and title
	ax1.set_ylabel(ylabel)
	ax1.set_xlabel(xlabel)
 	# plt.legend(loc='upper right')
	plt.title(samplename)
	# Plot a black line along the x axis
	plt.axhline(0, color='black')
	# Set all y axis labels to absolute values
	ax1.set_yticklabels([str(abs(int(x))) for x in ax1.get_yticks()])
	# Set a buffer around the edge
	plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
	print('Plot created')
	# output plot to file
	outname = samplename + '.pdf'
	outfile = PdfPages(outname)
	plt.savefig(outfile, format='pdf')
	outfile.close()
	# clear plot window to prevent overplotting
	plt.clf()
	return('done')

###########################################################################
# function to parse a bam file and plot seqlogos for a given length range #
###########################################################################
def bam_seqlogo(file,outfilename,minlength=0,maxlength=1000,unique=False):
	# empty dicts for read names, sequences and read counts
	posnames = []
	posseqs = []
	poscounts = []
	negnames = []
	negseqs = []
	negcounts = []
	postotal = 0
	negtotal = 0
	mappedtotal = 0
	# read in bamfile
	bamfile = pysam.AlignmentFile(file,"rb")
	for read in bamfile:
		# convert line into string to make it splittable
		line = str(read)
		# split line on tabs
		linesplit = line.split("\t")
		# this filters out reads with flag 4 (i.e. unmapped reads)
		cigar = linesplit[5]
		if cigar.endswith('M'):
			seq = linesplit[9].replace('T','U')
			count = int(linesplit[0].split('-')[1])
			name = linesplit[0].split('-')[0]
			strand = linesplit[1]
			# this is where the length filtering happens
			if minlength <= len(seq) <= maxlength:
				mappedtotal += 1
				# this is where the strand filtering happens
				# negative sense
				if strand == '16':
					negtotal += 1
					negnames.append(name)
					# convert seq string to Biopython Seq object
					tempseq=Seq(seq)
					# set reverse-complement of Seq object as string
					tempseqrevcomp = tempseq.reverse_complement()
					negseqs.append(str(tempseqrevcomp))
					negcounts.append(count)
				# positive sense
				else:
					postotal += 1
					posnames.append(name)
					posseqs.append(seq)
					poscounts.append(count)
	
	print('Total mapped reads = ' + str(mappedtotal) + '\nSense reads = ' + str(postotal) + '\nAntisense reads = ' + str(negtotal) + '\nMissing mapped reads = ' + str(mappedtotal-(postotal+negtotal)))
	# because seqlogo needs all reads to be the same length, we now need to pad the shorter reads with Ns (which will be ignored)
	# positive sense
	# sort the sequence list by length, and set the max length as the length of the last (i.e. longest) seq
	posmaxlen = len(sorted(posseqs, key=len)[len(posseqs)-1])
	# go through each seq, add as many Ns onto end as needed to bring length up to maxlength
	for i in range(len(posseqs)):
		lengthdiff = posmaxlen-len(posseqs[i])
		# if the seq is shorter than maxlength
		if lengthdiff>0:
			# set the seq as the initial newseq
			newseq = posseqs[i]
			# add as many Ns as necessary to bring newseq up to maxlength
			for j in range(lengthdiff):
				newseq = newseq + 'N'
			# replace the initial seq (which was shorter than maxlength) with the new seq (which is now padded to maxlength)
			posseqs[i] = newseq
	# negative sense
	# sort the sequence list by length, and set the max length as the length of the last (i.e. longest) seq
	negmaxlen = len(sorted(negseqs, key=len)[len(negseqs)-1])
	# go through each seq, add as many Ns onto end as needed to bring length up to maxlength
	for i in range(len(negseqs)):
		lengthdiff = negmaxlen-len(negseqs[i])
		# if the seq is shorter than maxlength
		if lengthdiff>0:
			# set the seq as the initial newseq
			newseq = negseqs[i]
			# add as many Ns as necessary to bring newseq up to maxlength
			for j in range(lengthdiff):
				newseq = newseq + 'N'
			# replace the initial seq (which was shorter than maxlength) with the new seq (which is now padded to maxlength)
			negseqs[i] = newseq
	# output the positive sense reads and call Seqlogo
	posout = ''
	for i in range(len(posseqs)):
		if unique==True:
			posout = posout + '>' + posnames[i] + '\n' + posseqs[i] + '\n'
		elif unique==False:
			for j in range(poscounts[i]):
				posout = posout + '>' + posnames[i] + '_' + str(j) + '\n' + posseqs[i] + '\n'
	posoutfile = open('pos.fas','wt')
	posoutfile.write(posout)
	posoutfile.close()
	command = 'weblogo -f pos.fas -D fasta -o ' + outfilename + '_sense.pdf -F pdf -A RNA -a \'ACGU\' -c classic --yaxis 1 --errorbars NO'
	subprocess.run(command, shell=True)
	subprocess.run('rm pos.fas', shell=True)
	# output the negative sense reads and call Seqlogo
	negout = ''
	for i in range(len(negseqs)):
		if unique==True:
			negout = negout + '>' + negnames[i] + '\n' + negseqs[i] + '\n'
		elif unique==False:
			for j in range(negcounts[i]):
				negout = negout + '>' + negnames[i] + '_' + str(j) + '\n' + negseqs[i] + '\n'
	negoutfile = open('neg.fas','wt')
	negoutfile.write(negout)
	negoutfile.close()
	command = 'weblogo -f neg.fas -D fasta -o ' + outfilename + '_antisense.pdf -F pdf -A RNA -a \'ACGU\' -c classic --yaxis 1 --errorbars NO'
	subprocess.run(command, shell=True)
	subprocess.run('rm neg.fas', shell=True)
	return('SeqLogos constructed')

##################################################################
# function to plot z-score of overlap of sense & antisense reads #
##################################################################
def overlap_plot(outname,bamfile,minlength,maxlength,minoverlap,maxoverlap):
	# convert bam to sam (needs to be sam for signature calculations
	bamsamcmd = "samtools view -h " + bamfile + " > " + bamfile.replace(".bam",".sam")
	call(bamsamcmd, shell=True)
	# calculate overlap signatures
	sigcmd = "python3 ~/bin/signature.py " + bamfile.replace(".bam",".sam") + " " + str(minlength) + " " + str(maxlength) + " " + str(minoverlap) + " " + str(maxoverlap) + " " + "overlaptable.txt"
	call(sigcmd, shell=True)
	overlap = []
	zscore = []
	for i in open('overlaptable.txt'):
		overlap.append(i.split('\t')[0])
		zscore.append(i.split('\t')[2])
	# remove the intermediate output files
	call('rm ' + bamfile.replace(".bam",".sam"), shell=True)
	call('rm overlaptable.txt', shell=True)
	# construct plot
	plt.plot(overlap,zscore)
	plt.xlabel('Overlap')
	plt.ylabel('Z Score')
	# output plot to file
	outfile = PdfPages(outname)
	plt.savefig(outfile, format='pdf')
	# close outfile
	outfile.close()
	# clear plot window to prevent overplotting
	plt.clf()
	print('Overlap signature plotted')

##################
# plotting calls #
##################
if SeqLogo is True:
	print('Plotting seqlogo (sample='+Sample+', input file='+InFile+')')
	bam_seqlogo(file=InFile,outfilename=InFile.strip('.bam')+'_'+Sample,minlength=MinLength,maxlength=MaxLength,unique=Unique)
	print('Plotted seqlogo (sample=' + Sample + ', input file=' + InFile+')')

if PingPong is True:
	print('Plotting Ping-Pong signature (sample='+Sample+', input file='+InFile+')')
	overlap_plot(outname=InFile.strip('.bam')+'_'+Sample+'_PingPong.pdf',bamfile=InFile,minlength=MinLength,maxlength=MaxLength,minoverlap=1,maxoverlap=30)
	print('Plotted Ping-Pong signature (sample='+Sample+', input file='+InFile+')')

if Distribution is True:
	if Unstranded==False:
		print('Plotting sRNA length distribution (sample='+Sample+', input file='+InFile+')')
		fiveprime_host = bam_5prime_stranded(file=InFile,minlength=MinLength,maxlength=MaxLength,unique=Unique)
		stacked_barplot_stranded(dataframe=fiveprime_host,samplename=InFile.strip('.bam')+'_'+Sample+'_5prime',xlabel='Length',ylabel='Count')
		print('Plotted sRNA length distribution (sample='+Sample+', input file='+InFile+')')
	if Unstranded==True:
		print('Plotting sRNA length distribution (sample='+Sample+', input file='+InFile+')')
		fiveprime_host = bam_5prime(file=InFile,minlength=MinLength,maxlength=MaxLength,unique=Unique)
		stacked_barplot(dataframe=fiveprime_host,samplename=InFile.strip('.bam')+'_'+Sample+'_5prime',xlabel='Length',ylabel='Count')
		print('Plotted sRNA length distribution (sample='+Sample+', input file='+InFile+')')

print("Plotting complete")
