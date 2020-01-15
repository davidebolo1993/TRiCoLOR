#!/usr/bin/env python3

#usage:
#python entropy.py -g ref.fa -b known.reps.bed

import argparse
from argparse import HelpFormatter
import os
import math
import random
import csv
from shutil import which
import pysam
import numpy as np
import subprocess, shlex


def main():


	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Simulate long-reads BAM centered on a known STRs randomly picked from a STRs BED. Compute entropy for non-overlapping, sliding windows of a certain size. Define the optimal entropy treshold for ONT/PB''', epilog='''This script was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	parameters = parser.add_argument_group('Arguments')

	parameters.add_argument('-g', '--genome', help='Reference genome', required=True)
	parameters.add_argument('-b', '--bed', help='BED containing known tandem repetitions for the reference genome in -g/--genome', required=True)
	parameters.add_argument('-s', '--size', help='Size of the sliding window. Used to train the entropy treshold [20]', default=20, type=int)

	args = parser.parse_args()

	ONT_error="45:25:30"
	PB_error="15:50:35"

	with open(os.path.abspath(args.genome),'r') as file:

		assert(file.readline().startswith('>')) #stop on error

	assert(os.access(os.path.abspath(os.path.dirname(args.genome)), os.W_OK)) #stop on error

	assert(os.path.exists(os.path.abspath(args.bed))) #stop on error

	assert(which('VISOR') is not None) #stop on error

	if not os.path.exists('EntropySimulations'):

		os.makedirs('EntropySimulations')

	genomereps=[]

	with open(os.path.abspath(args.bed), 'r') as f:

		for line in csv.reader(f, delimiter='\t'):

			genomereps.append((line[0], int(line[1]), int(line[2]), line[3].split('x')[1]))

	print('Mean STRs length is ' + str(mean([c-b for (a,b,c,d) in genomereps])))


	#simulate

	number_of_simulations=500
	counter = 0

	all_ONT=[]
	all_PB=[]

	while counter < number_of_simulations:

		rep=random.choice(genomereps)

		if rep[1] - 10000 < 0: #discard STRs which start is too close to the beginning of the chromosome

			genomereps.remove(rep)
			continue

		genomereps.remove(rep)

		command1="samtools faidx " + os.path.abspath(args.genome) + " " + rep[0]
		
		if not os.path.exists(os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa')):

			with open(os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa'), 'w') as outfa:

				subprocess.call(shlex.split(command1), stdout=outfa)

		out=os.path.abspath('EntropySimulations/Sim' + str(counter+1))
		os.makedirs(out)
		os.makedirs(out + '/haplotype')

		with open(os.path.abspath(out + '/haplotype/' + rep[0] + '.fa'), 'w') as outfa:

			subprocess.call(shlex.split(command1), stdout=outfa)

		with open(os.path.abspath(out + '/sim.bed'), 'w') as outsim:

			outsim.write(rep[0] + '\t' + str(rep[1]-10000) + '\t' + str(rep[2]+10000) + '\t' + str(100.0) + '\t' + str(100.0) + '\n')

		command2="VISOR LASeR -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + ' -bed ' + os.path.abspath(out + '/sim.bed') + ' -s ' +  os.path.abspath(out + '/haplotype') + ' -c 10 --threads 6 --noaddtag -r ' +  ONT_error + ' --readstype ONT -o ' +  str(out + '/ONT_bam')
		command3="VISOR LASeR -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + ' -bed ' + os.path.abspath(out + '/sim.bed') + ' -s ' +  os.path.abspath(out + '/haplotype') + ' -c 10 --threads 6 --noaddtag -r ' +  PB_error + ' --readstype PB -o ' +  str(out + '/PB_bam')

		subprocess.call(shlex.split(command2))
		subprocess.call(shlex.split(command3))

		ONT_bam=os.path.abspath(out + '/ONT_bam/sim.srt.bam')
		PB_bam=os.path.abspath(out + '/PB_bam/sim.srt.bam')

		ONTen_=BamScanner(ONT_bam, rep[0], rep[1], rep[2], args.size)
		PBen_=BamScanner(PB_bam, rep[0], rep[1], rep[2], args.size)

		ONT_onlyen=[]
		PB_onlyen=[]

		for i in range(len(ONTen_)):

				for l in ONTen_[i]:

					ONT_onlyen.append(l[0])

		for i in range(len(PBen_)):

				for l in PBen_[i]:

					PB_onlyen.append(l[0])

		all_ONT.extend(ONT_onlyen)
		all_PB.extend(PB_onlyen)

		#one can insert here a condition to store a tsv file for plotting entropy signal for a specific STRs (we plot in R using ggplot2)

		counter +=1

	ONT_value=np.percentile(all_ONT,2)
	PB_value=np.percentile(all_PB,2)

	print('Treshold for ONT data is ' + str(ONT_value))
	print('Treshold for PB data is ' + str(PB_value))

	print(str(sum(i > ONT_value for i in all_ONT)/len(all_ONT)*100) + ' windows excluded in ONT')
	print(str(sum(i > PB_value for i in all_PB)/len(all_PB)*100) + ' windows excluded in PB')

	#strore tsv files for plotting in R

	#if len(all_ONT) <= 50000:

		#ONT_dist_downsampled=all_ONT

	#else:

		#ONT_dist_downsampled=random.sample(all_ONT, 50000)

	#if len(all_PB) <= 50000:

		#PB_dist_downsampled=all_PB

	#else:

		#PB_dist_downsampled=random.sample(all_PB, 50000)

	#with open(os.path.abspath('ONTdist.tsv'), 'w') as dfout:

		#for el in ONT_dist_downsampled:

			#dfout.write(str(el) + '\n')

	#with open(os.path.abspath('PBdist.tsv'), 'w') as dfout:

		#for el in PB_dist_downsampled:

			#dfout.write(str(el) + '\n')


class CustomFormat(HelpFormatter):


	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()

def mean(data):


	n = len(data)

	return sum(data)/n 

def _ss(data):


	c = mean(data)
	ss = sum((x-c)**2 for x in data)
	
	return ss

def stddev(data, ddof=0):


	n = len(data)
	ss = _ss(data)
	pvar = ss/(n-ddof)

	return pvar**0.5


def entropy(string):


	prob = [float(string.count(c)) / len(string) for c in dict.fromkeys(list(string))]
	entropy = - sum([p * math.log(p) / math.log(2.0) for p in prob])

	return entropy


def BamScanner(bamfile,chromosome,start,end,scansize):


	BamFile=pysam.AlignmentFile(bamfile,'rb')

	ent__=[]

	for read in BamFile.fetch(chromosome, start, end):

		ent_=[]

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			sequence=read.seq
			coords=read.get_reference_positions(full_length=True)

			start=0
			end=scansize
			terminal_ind=len(sequence)-1

			while terminal_ind > end:

				if terminal_ind-end >= scansize:

					ent_.append((entropy(sequence[start:end]), coords[start:end]))
					start+=scansize
					end+=scansize

				else:

					break

		ent__.append(ent_)

	return ent__



if __name__ == "__main__":

	main()
