#!/usr/bin/env python3

from __future__ import print_function
import subprocess, shlex
from shutil import which
import random
import csv
import os
import sys
import pandas as pd
import re 
import glob
import pysam

import argparse
from argparse import HelpFormatter



def main():

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Simulations for TRiCoLOR program''', epilog='''This script was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 
	
	required=parser.add_argument_group('Required I/O arguments')

	required.add_argument('-g', '--genome', help='reference genome', metavar='.fa',required=True)	
	required.add_argument('-b', '--bed', help='BED containing known tandem repetitions for the reference genome in -g/--genome', metavar='BED',required=True)	


	specific=parser.add_argument_group('Simulations parameters')

	specific.add_argument('-n', '--number', help='number of simulations [100]', type=int, default=1, metavar='')
	specific.add_argument('-c', '--coveragemin', help='minimum coverage for the simulation [5]', type=int, default=5, metavar='')
	specific.add_argument('-C', '--coveragemax', help='maximum coverage for the simulation [10]', type=int, default=10, metavar='')
	specific.add_argument('-a', '--accuracy', help='mean accuracy rate [0.9]', type=float, default=0.9, metavar='')
	specific.add_argument('-l', '--length', help='mean length for simulated reads [8000]', metavar='', default=8000, type=int)	
	specific.add_argument('-s', '--size', help='size of contraction/expansion [7]', type=int, metavar='', default=7)
	specific.add_argument('--threads', help='number of cores to use for alignments [1]', metavar='', type=int, default=1)
	specific.add_argument('--readstype', help='Type of long reads (ONT, PB) [ONT]', default='ONT', choices=['PB', 'ONT'], metavar='')

	args = parser.parse_args()

	with open(os.path.abspath(args.genome),'r') as file:

		assert(file.readline().startswith('>')) #stop on error

	assert(os.access(os.path.abspath(os.path.dirname(args.genome)), os.W_OK)) #stop on error

	assert(os.path.exists(os.path.abspath(args.bed))) #stop on error

	assert(which('VISOR') is not None) #stop on error

	if not os.path.exists('PrecisionRecallSimulations'):

		os.makedirs('PrecisionRecallSimulations')
		os.makedirs('PrecisionRecallSimulations/contractions')
		os.makedirs('PrecisionRecallSimulations/expansions')

	genomereps=[]

	with open(os.path.abspath(args.bed), 'r') as f:

		for line in csv.reader(f, delimiter='\t'):

			genomereps.append((line[0], int(line[1]), int(line[2]), line[3].split('x')[1], line[3].split('x')[0]))

	if args.readstype == "ONT":

		ratio="45:25:30"

	else:

		ratio="15:50:35"

	assert(args.size <= 10) #for contractions it can be a problem otherwise: minimum repeated period for STRs in GRCh38 is 15 

	counter=0
	acceptable=5

	progress = ProgressBar(args.number*2, fmt=ProgressBar.FULL)
	
	while counter < args.number:

		rep=random.choice(genomereps)

		if rep[0] == "chrY" or rep[1] - 50000 <= 0:

			genomereps.remove(rep)
			continue

		with open('h1.contraction.bed', 'w') as bed1out, open('h2.contraction.bed', 'w') as bed2out:

			bed1out.write(rep[0] + "\t" + str(rep[1]) + "\t" + str(rep[2]) + "\t" + "tandem repeat contraction" + "\t" + rep[3] + ":" + str(args.size) + "\t" + "0" + "\n") #contraction
			bed2out.write(rep[0] + "\t" + str(rep[1]) + "\t" + str(rep[2]) + "\t" + "tandem repeat contraction" + "\t" + rep[3] + ":0" + "\t" + "0" + "\n") #no contraction

		with open('h1.expansion.bed', 'w') as bed1out, open('h2.expansion.bed', 'w') as bed2out:

			bed1out.write(rep[0] + "\t" + str(rep[1]) + "\t" + str(rep[2]) + "\t" + "tandem repeat expansion" + "\t" + rep[3] + ":" + str(args.size) + "\t" + "0" + "\n") #expansion
			bed2out.write(rep[0] + "\t" + str(rep[1]) + "\t" + str(rep[2]) + "\t" + "tandem repeat expansion" + "\t" + rep[3] + ":0" + "\t" + "0" + "\n") #no expansion


		with open('sim.bed', 'w') as infoout:

			infoout.write(rep[0] + '\t' + str(rep[1]-50000) + '\t' + str(rep[2] + 50000) + '\t' + str(100.0) + '\t' + str(100.0) + '\n') #simulate around the repetition

		command1="samtools faidx " + os.path.abspath(args.genome) + " " + rep[0]
		
		if not os.path.exists(os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa')):

			with open(os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa'), 'w') as outfa:

				subprocess.call(shlex.split(command1), stdout=outfa)


		outhc="PrecisionRecallSimulations/contractions/simfasta" + str(counter+1)
		outhe="PrecisionRecallSimulations/expansions/simfasta" + str(counter+1)

		commandhackdel1="VISOR HACk -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed h1.contraction.bed -o " + outhc + "/h1"
		commandhackdel2="VISOR HACk -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed h2.contraction.bed -o " + outhc + "/h2"

		subprocess.call(shlex.split(commandhackdel1))
		subprocess.call(shlex.split(commandhackdel2))


		commandhackexp1="VISOR HACk -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed h1.expansion.bed -o " + outhe + "/h1"
		commandhackexp2="VISOR HACk -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed h2.expansion.bed -o " + outhe + "/h2"

		subprocess.call(shlex.split(commandhackexp1))
		subprocess.call(shlex.split(commandhackexp2))

		outlc="PrecisionRecallSimulations/contractions/simbam" + str(counter+1)
		outle="PrecisionRecallSimulations/expansions/simbam" + str(counter+1)

		with open('test.bed', 'w') as outtest:

			outtest.write(rep[0] + '\t' + str(rep[1]-500) + '\t' + str(rep[2] + 500) + '\n')

		donecont1=False
		donecont2=False

		while not donecont1:

			commandlaserdel1="VISOR LASeR -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed sim.bed -s " + outhc + "/h1" + " -c " + str(args.coveragemin*2) + " --threads " +  str(args.threads) + " -a " + str(args.accuracy) +" -r " + ratio + " -l " + str(args.length) + " --noaddtag -o " + outlc + "/h1"

			subprocess.call(shlex.split(commandlaserdel1))

			counter1,bam1=check_coverage(pysam.AlignmentFile(outlc + "/h1/sim.srt.bam", 'rb'), rep[0], rep[1]-500, rep[2]+500, args.coveragemin, args.coveragemax)

			if bam1:

				donecont1=True

			else:

				commandremovedel="rm -r " + outlc + "/h1"
				subprocess.call(shlex.split(commandremovedel))
				print('Coverage outside wanted range for current simulation. Repeat')

		while not donecont2:

			commandlaserdel2="VISOR LASeR -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed sim.bed -s " + outhc + "/h2" + " -c " + str(args.coveragemin*2) + " --threads " +  str(args.threads) + " -a " + str(args.accuracy) +" -r " + ratio + " -l " + str(args.length) + " --noaddtag -o " + outlc + "/h2"

			subprocess.call(shlex.split(commandlaserdel2))

			counter2,bam2=check_coverage(pysam.AlignmentFile(outlc + "/h2/sim.srt.bam", 'rb'), rep[0], rep[1]-500, rep[2]+500, args.coveragemin, args.coveragemax)

			if bam2:

				donecont2=True

			else:

				commandremovedel="rm -r " + outlc + "/h2"
				subprocess.call(shlex.split(commandremovedel))
				print('Coverage outside wanted range for current simulation. Repeat')


		doneexp1=False
		doneexp2=False

		while not doneexp1:

			commandlaserexp1="VISOR LASeR -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed sim.bed -s " + outhe + "/h1" + " -c " + str(args.coveragemin*2) + " --threads " +  str(args.threads) + " -a " + str(args.accuracy) +" -r " + ratio + " -l " + str(args.length) + " --noaddtag -o " + outle + "/h1"

			subprocess.call(shlex.split(commandlaserexp1))

			counter1,bam1=check_coverage(pysam.AlignmentFile(outle + "/h1/sim.srt.bam", 'rb'), rep[0], rep[1]-500, rep[2]+500, args.coveragemin, args.coveragemax)

			if bam1:

				doneexp1=True

			else:

				commandremoveexp="rm -r " + outle + "/h1"
				subprocess.call(shlex.split(commandremoveexp))
				print('Coverage outside wanted range for current simulation. Repeat')


		while not doneexp2:

			commandlaserexp2="VISOR LASeR -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed sim.bed -s " + outhe + "/h2" + " -c " + str(args.coveragemin*2) + " --threads " +  str(args.threads) + " -a " + str(args.accuracy) +" -r " + ratio + " -l " + str(args.length) + " --noaddtag -o " + outle + "/h2"

			subprocess.call(shlex.split(commandlaserexp2))

			counter2,bam2=check_coverage(pysam.AlignmentFile(outle + "/h2/sim.srt.bam", 'rb'), rep[0], rep[1]-500, rep[2]+500, args.coveragemin, args.coveragemax)

			if bam2:

				doneexp2=True

			else:

				commandremoveexp="rm -r " + outle + "/h2"
				subprocess.call(shlex.split(commandremoveexp))
				print('Coverage outside wanted range for current simulation. Repeat')


		minimumcontractionlen=len(rep[3])*(int(rep[4])-args.size-2) 
		minimumexpansionlen=len(rep[3])*int(rep[4])

		outtc="PrecisionRecallSimulations/contractions/test" + str(counter+1)
		outte="PrecisionRecallSimulations/expansions/test" + str(counter+1)

		contractiontricolor="TRiCoLOR REFER -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed test.bed -bam " + outlc + "/h1/sim.srt.bam " + outlc +  "/h2/sim.srt.bam -m " + str(len(rep[3])) + " --precisemotif -s "  + str(minimumcontractionlen) + " --readstype " + args.readstype + " -o " + outtc
		expansiontricolor="TRiCoLOR REFER -g " + os.path.abspath(os.path.dirname(args.genome) +'/' + rep[0] + '.fa') + " -bed test.bed -bam " + outle + "/h1/sim.srt.bam " + outle +  "/h2/sim.srt.bam -m " + str(len(rep[3])) + " --precisemotif -s "  + str(minimumexpansionlen) + " --readstype " + args.readstype + " -o " + outte

		subprocess.call(shlex.split(contractiontricolor))
		subprocess.call(shlex.split(expansiontricolor))


		#contractions

		try:

			tabh1=pd.read_csv(outtc + '/haplotype1/' + rep[0] +'.repetitions.bed', sep="\t")
			tabh2=pd.read_csv(outtc + '/haplotype2/' + rep[0] +'.repetitions.bed', sep="\t")

			m1=[]
			n1=[]

			for b1,b2,m,n in zip(tabh1['Start'], tabh1['End'], tabh1['Repeated Motif'], tabh1['Repetitions Number']):

				if abs(b1-rep[1]) <= acceptable or abs(b2-(rep[2])) <= acceptable:

					if m in possible_rotations(rep[3]):

						m1.append(m)

					else:

						m1.append('wrong')

					n1.append(n)

			m2=[]
			n2=[]

			for b1,b2,m,n in zip(tabh2['Start'], tabh2['End'], tabh2['Repeated Motif'], tabh2['Repetitions Number']):

				if abs(b1-rep[1]) <= acceptable or abs(b2-(rep[2])) <= acceptable:

					if m in possible_rotations(rep[3]):

						m2.append(m)

					else:

						m2.append('wrong')
						
					n2.append(n)


			outtab=pd.DataFrame({'H1N':n1,'H2N':n2, 'H1M': m1, 'H2M':m2, 'TN':[rep[4]], 'TM':[rep[3]]})
			outtab.to_csv(outtc + '/result.tsv',sep='\t',index=False)

		except:

			removefasta="rm -r " + outhc
			removebam="rm -r " + outlc
			removetricolor="rm -r " + outtc
			print('Wrong mapping or very few repetitions to detect. Retry')


		#expansions

		try:

			tabh1=pd.read_csv(outte + '/haplotype1/' + rep[0] +'.repetitions.bed', sep="\t")
			tabh2=pd.read_csv(outte + '/haplotype2/' + rep[0] +'.repetitions.bed', sep="\t")

			m1=[]
			n1=[]

			for b1,b2,m,n in zip(tabh1['Start'], tabh1['End'], tabh1['Repeated Motif'], tabh1['Repetitions Number']):

				if abs(b1-rep[1]) <= acceptable or abs(b2-(rep[2])) <= acceptable:

					if m in possible_rotations(rep[3]):

						m1.append(m)

					else:

						m1.append('wrong')

					n1.append(n)

			m2=[]
			n2=[]

			for b1,b2,m,n in zip(tabh2['Start'], tabh2['End'], tabh2['Repeated Motif'], tabh2['Repetitions Number']):

				if abs(b1-rep[1]) <= acceptable or abs(b2-(rep[2])) <= acceptable:

					if m in possible_rotations(rep[3]):

						m2.append(m)

					else:

						m2.append('wrong')
						
					n2.append(n)


			outtab=pd.DataFrame({'H1N':n1,'H2N':n2, 'H1M': m1, 'H2M':m2, 'TN':[rep[4]], 'TM':[rep[3]]})
			outtab.to_csv(outte + '/result.tsv',sep='\t',index=False)

		except:

			removefasta="rm -r " + outhe
			removebam="rm -r " + outle
			removetricolor="rm -r " + outte
			print('Wrong mapping or very few repetitions to detect. Retry')


		progress.current += 1
		progress()
		counter +=1


	os.remove(os.path.abspath('h1.contraction.bed'))
	os.remove(os.path.abspath('h2.contraction.bed'))
	os.remove(os.path.abspath('h1.expansion.bed'))
	os.remove(os.path.abspath('h2.expansion.bed'))
	os.remove(os.path.abspath('sim.bed'))
	os.remove(os.path.abspath('test.bed'))

	print('Done with simulations')

	C0SS,C0PR,C0C=Collect(args.size,'contractions', allowed=0,output="PrecisionRecallSimulations")

	with open('PrecisionRecallSimulations/contractions/SummaryC0.txt', 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(C0SS[0]) + '-' +str(C0SS[1]) + '\n' + 'Precision and Recall: ' + str(C0PR[0]) + '-' + str(C0PR[1]) + '\n' + 'Core_ESW: ' + str(C0C[0]) + '-' + str(C0C[1]) + '-' + str(C0C[2]) + '\n') 

	C1SS,C1PR,C1C=Collect(args.size,'contractions', allowed=1,output="PrecisionRecallSimulations")

	with open('PrecisionRecallSimulations/contractions/SummaryC1.txt', 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(C1SS[0]) + '-' +str(C1SS[1]) + '\n' + 'Precision and Recall: ' + str(C1PR[0]) + '-' + str(C1PR[1]) + '\n' + 'Core_ESW: ' + str(C1C[0]) + '-' + str(C1C[1]) + '-' + str(C1C[2]) + '\n') 

	C2SS,C2PR,C2C=Collect(args.size,'contractions', allowed=2,output="PrecisionRecallSimulations")

	with open('PrecisionRecallSimulations/contractions/SummaryC2.txt', 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(C2SS[0]) + '-' +str(C2SS[1]) + '\n' + 'Precision and Recall: ' + str(C2PR[0]) + '-' + str(C2PR[1]) + '\n' + 'Core_ESW: ' + str(C2C[0]) + '-' + str(C2C[1]) + '-' + str(C2C[2]) + '\n') 

	E0SS,E0PR,E0C=Collect(args.size,'expansions', allowed=0,output="PrecisionRecallSimulations")

	with open('PrecisionRecallSimulations/expansions/SummaryE0.txt', 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(E0SS[0]) + '-' +str(E0SS[1]) + '\n' + 'Precision and Recall: ' + str(E0PR[0]) + '-' + str(E0PR[1]) + '\n' + 'Core_ESW: ' + str(E0C[0]) + '-' + str(E0C[1]) + '-' + str(E0C[2]) + '\n') 

	E1SS,E1PR,E1C=Collect(args.size,'expansions', allowed=1,output="PrecisionRecallSimulations")

	with open('PrecisionRecallSimulations/expansions/SummaryE1.txt', 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(E1SS[0]) + '-' +str(E1SS[1]) + '\n' + 'Precision and Recall: ' + str(E1PR[0]) + '-' + str(E1PR[1]) + '\n' + 'Core_ESW: ' + str(E1C[0]) + '-' + str(E1C[1]) + '-' + str(E1C[2]) + '\n') 

	E2SS,E2PR,E2C=Collect(args.size,'expansions', allowed=2,output="PrecisionRecallSimulations")

	with open('PrecisionRecallSimulations/expansions/SummaryE2.txt', 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(E2SS[0]) + '-' +str(E2SS[1]) + '\n' + 'Precision and Recall: ' + str(E2PR[0]) + '-' + str(E2PR[1]) + '\n' + 'Core_ESW: ' + str(E2C[0]) + '-' + str(E2C[1]) + '-' + str(E2C[2]) + '\n') 

	print('Done')

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



class ProgressBar(object):

	DEFAULT = 'Progress: %(bar)s %(percent)3d%%'
	FULL = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'

	def __init__(self, total, width=40, fmt=DEFAULT, symbol='=',output=sys.stderr):
		
		assert len(symbol) == 1

		self.total = total
		self.width = width
		self.symbol = symbol
		self.output = output
		self.fmt = re.sub(r'(?P<name>%\(.+?\))d',r'\g<name>%dd' % len(str(total)), fmt)
		self.current = 0

	def __call__(self):

		percent = self.current / float(self.total)
		size = int(self.width * percent)
		remaining = self.total - self.current
		bar = '[' + self.symbol * size + ' ' * (self.width - size) + ']'
		args = {
			'total': self.total,
			'bar': bar,
			'current': self.current,
			'percent': percent * 100,
			'remaining': remaining
		}

		print('\r' + self.fmt % args, file=self.output, end='')

	def done(self):

		self.current = self.total
		self()
		print('', file=self.output)


def possible_rotations(word):

    p = []

    for i in range(len(word)):

        p.append(word[i:]+word[:i])
        
    return p


def check_coverage(pysam_AlignmentFile, chromosome, start, end, coveragemin, coveragemax): #check if we have at least the wanted coverage for the interval


	counter = 0

	for read in pysam_AlignmentFile.fetch(chromosome, start, end):
		
		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			if read.reference_start <= start and read.reference_end >= end: #make sure that the same read cover the entire region

				counter +=1

	print('Coverage is ' + str(counter))

	if counter >= coveragemin and counter <= coveragemax:

		return counter, True

	else:

		return counter, False


def Collect(size,mode,allowed, output):


	folders=glob.glob(os.path.abspath(output + '/' + mode) + '/test*')

	TP=0 #repetitions with expected contratction/expansion, calculated from h1
	FP=0 #repetitions with unexpected contraction/expansion, calculated from h2
	TN=0 #expected repetitions without contraction/expansion, calculated from h2
	FN=0 #unexpected reptitions without contraction/expansion, calculated from h1

	EC=0
	SC=0
	WC=0


	for f in folders:

		tab=pd.read_csv(os.path.abspath(f) + '/result.tsv', sep='\t')

		correctN=int(tab['TN'])
		correctM=str(tab['TM']).split('\n')[0].split('   ')[1].split(' ')[1]


		if mode == 'contractions':


			H1N=int(tab['H1N'])
			H1M=str(tab['H1M']).split('\n')[0].split('   ')[1].split(' ')[1]


			H2N=int(tab['H2N'])
			H2M=str(tab['H2M']).split('\n')[0].split('   ')[1].split(' ')[1]


			if abs(correctN-H1N) <= size + allowed:

				TP +=1

			else:

				FN +=1


			if abs(correctN - H2N) <= allowed:

				TN +=1

			else:

				FP +=1


			if H1M == 'wrong':

				WC +=1

			else:

				if H1M == correctM:

					EC +=1

				else:

					SC +=1


			if H2M == 'wrong':

				WC +=1

			else:

				if H2M == correctM:

					EC +=1

				else:

					SC +=1

		else:

			H1N=int(tab['H1N'])
			H1M=str(tab['H1M']).split('\n')[0].split('   ')[1].split(' ')[1]


			H2N=int(tab['H2N'])
			H2M=str(tab['H2M']).split('\n')[0].split('   ')[1].split(' ')[1]


			if abs(correctN-H1N) <= size + allowed:

				TP +=1

			else:

				FN +=1


			if abs(correctN - H2N) <= allowed:

				TN +=1

			else:

				FP +=1


			if H1M == 'wrong':

				WC +=1

			else:

				if H1M == correctM:

					EC +=1

				else:

					SC +=1


			if H2M == 'wrong':

				WC +=1

			else:

				if H2M == correctM:

					EC +=1

				else:

					SC +=1

	sensitivity=TP/(TP+FN) if (TP+FN) else 0
	specificity=TN/(TN+FP) if (TN+FP) else 0

	precision=TP/(TP+FP) if (TP+FP) else 0
	recall=TP/(TP+FN) if (TP+FN) else 0

	exact_core=EC/(EC+SC+WC) if (EC+SC+WC) else 0
	shifted_core=SC/(EC+SC+WC) if (EC+SC+WC) else 0
	wrong_core=WC/(EC+SC+WC) if (EC+SC+WC) else 0

	return (sensitivity,specificity), (precision,recall), (exact_core,shifted_core,wrong_core)



if __name__=='__main__':

	main()
