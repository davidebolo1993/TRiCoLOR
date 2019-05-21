from __future__ import print_function
import subprocess
import os
import sys
import pandas as pd
import re 
import glob
import logging
import pysam

import argparse
from argparse import HelpFormatter



def main():

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Simulations for TRiCoLOR program''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 
	
	required=parser.add_argument_group('Required I/O arguments')

	required.add_argument('-g', '--genome', help='reference genome', metavar='.fa',required=True)	

	specific=parser.add_argument_group('Simulations parameters')

	specific.add_argument('-n', '--number', help='number of simulations [100]', type=int, default=100, metavar='')
	specific.add_argument('-c', '--coveragemin', help='minimum coverage for the simulation [5]', type=int, default=5, metavar='')
	specific.add_argument('-C', '--coveragemax', help='maximum coverage for the simulation [10]', type=int, default=10, metavar='')
	specific.add_argument('-a', '--accuracy', help='mean accuracy rate [0.9]', type=float, default=0.9, metavar='')
	specific.add_argument('-s', '--size', help='size of contraction/expansion [7]', type=int, metavar='', default=7)

	args = parser.parse_args()


	logging.basicConfig(filename=os.path.abspath(os.getcwd() + '/logfile.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

	Simulate(args.genome,args.number, args.accuracy, args.coveragemin, args.coveragemax, args.size)

	C1SS,C1PR,C1C=Collect(args.size,'contractions', allowed=1)

	with open(os.path.abspath('contractions/SummaryC1.txt'), 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(C1SS[0]) + '-' +str(C1SS[1]) + '\n' + 'Precision and Recall: ' + str(C1PR[0]) + '-' + str(C1PR[1]) + '\n' + 'Core_ESW: ' + str(C1C[0]) + '-' + str(C1C[1]) + '-' + str(C1C[2]) + '\n') 

	C2SS,C2PR,C2C=Collect(args.size,'contractions', allowed=2)

	with open(os.path.abspath('contractions/SummaryC2.txt'), 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(C2SS[0]) + '-' +str(C2SS[1]) + '\n' + 'Precision and Recall: ' + str(C2PR[0]) + '-' + str(C2PR[1]) + '\n' + 'Core_ESW: ' + str(C2C[0]) + '-' + str(C2C[1]) + '-' + str(C2C[2]) + '\n') 


	E1SS,E1PR,E1C=Collect(args.size,'expansions', allowed=1)


	with open(os.path.abspath('expansions/SummaryE1.txt'), 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(E1SS[0]) + '-' +str(E1SS[1]) + '\n' + 'Precision and Recall: ' + str(E1PR[0]) + '-' + str(E1PR[1]) + '\n' + 'Core_ESW: ' + str(E1C[0]) + '-' + str(E1C[1]) + '-' + str(E1C[2]) + '\n') 

	E2SS,E2PR,E2C=Collect(args.size,'expansions', allowed=2)


	with open(os.path.abspath('expansions/SummaryE2.txt'), 'w') as textout:

		textout.write('Sensitivity and Specificity: ' + str(E2SS[0]) + '-' +str(E2SS[1]) + '\n' + 'Precision and Recall: ' + str(E2PR[0]) + '-' + str(E2PR[1]) + '\n' + 'Core_ESW: ' + str(E2C[0]) + '-' + str(E2C[1]) + '-' + str(E2C[2]) + '\n') 



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

	logging.info('Coverage is ' + str(counter))

	if counter >= coveragemin and counter <= coveragemax:

		return counter, True

	else:

		return counter, False



def Simulate(reference, number_of_simulations, accuracy, coveragemin, coveragemax, size):

	reference=os.path.abspath(reference)	
	acceptable=5 #look for repetitions close to true coordinates +- 5bp

	#contraction first

	out=os.path.abspath('contractions')

	if not os.path.exists(out):

		os.makedirs(out)

	progress = ProgressBar(number_of_simulations*2, fmt=ProgressBar.FULL)
	
	i=1

	while i < number_of_simulations+1:

		logging.info('Contractions. Iteration: ' + str(i))

		#contraction first

		proceed=False

		while not proceed:

			subprocess.call(['bash', 'contractknown.sh', str(size)], stderr=open(os.devnull, 'wb'))

			with open(os.path.abspath('random.region.bed'), 'r') as infoin:

				for line in infoin:

					chromosome=line.split('\t')[0]
					start=int(line.split('\t')[1])
					end=int(line.split('\t')[2])
					motif=line.split('\t')[3].split('x')[1].split('\n')[0]
					number=int(line.split('\t')[3].split('x')[0])

					if start - 50000 < 0:

						proceed=False

					else:

						proceed=True

		with open(os.path.abspath('sim.bed'), 'w') as infoout:

			infoout.write(chromosome + '\t' + str(start-50000) + '\t' + str(start + 50000) + '\t' + str(100.0) + '\t' + str(100.0) + '\n')

		#keep haplotype separated


		if not os.path.exists(os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa')):

			with open(os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), 'w') as fout:

				subprocess.call(['samtools', 'faidx', os.path.abspath(reference), chromosome], stdout=fout)


		subprocess.call(['VISOR', 'HACk', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('h1.bed'), '-o', os.path.abspath(out + '/simfasta' + str(i) + '/h1')])
		subprocess.call(['VISOR', 'HACk', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('h2.bed'), '-o', os.path.abspath(out + '/simfasta' + str(i) + '/h2')])


		done=False

		while not done:

		
			subprocess.call(['VISOR', 'LASeR', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('sim.bed'), '-s', os.path.abspath(out + '/simfasta' + str(i) + '/h1'), '-c', str(coveragemin*2), '-th', str(7), '--noaddtag','-o', os.path.abspath(out + '/simbam' + str(i) + '/h1')])
			subprocess.call(['VISOR', 'LASeR', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('sim.bed'), '-s', os.path.abspath(out + '/simfasta' + str(i) + '/h2'), '-c', str(coveragemin*2), '-th', str(7), '--noaddtag','-o', os.path.abspath(out + '/simbam' + str(i) + '/h2')])


			with open(os.path.abspath('test.bed'), 'w') as outtest:

				outtest.write(chromosome + '\t' + str(start-500) + '\t' + str(start + 500) + '\n')


			counter1,bam1=check_coverage(pysam.AlignmentFile(os.path.abspath(out + '/simbam' + str(i) + '/h1/sim.srt.bam')), chromosome, start-500, start+500, coveragemin, coveragemax)
			counter2,bam2=check_coverage(pysam.AlignmentFile(os.path.abspath(out + '/simbam' + str(i) + '/h2/sim.srt.bam')), chromosome, start-500, start+500, coveragemin, coveragemax)

			
			if counter1==0 or counter2 == 0:

				done=True
				continue


			if bam1 and bam2:

				done=True

			else:

				subprocess.call(['rm', '-r', os.path.abspath(out + '/simbam' + str(i))])
				logging.info('Coverage outside wanted range for current iteration. Repeat')


		if counter1==0 or counter2 == 0:

			subprocess.call(['rm', '-r', os.path.abspath(out + '/simfasta' + str(i))])
			subprocess.call(['rm', '-r', os.path.abspath(out + '/simbam' + str(i))])
			continue



		totlen=len(motif)*(number-size-2)

		if totlen > 15:

			totlen=15

		subprocess.call(['TRiCoLOR', 'REFER', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('test.bed'), '-bam', os.path.abspath(out + '/simbam' + str(i) + '/h1/sim.srt.bam'), os.path.abspath(out + '/simbam' + str(i) + '/h2/sim.srt.bam'), '-o', os.path.abspath(out + '/test' + str(i)), '-t', str(5), '-m', str(len(motif)), '--precisemotif', '-s', str(totlen)])



		try:

			tabh1=pd.read_csv(os.path.abspath(out + '/test' + str(i) + '/haplotype1/' + chromosome +'.repetitions.bed'), sep='\t')
			tabh2=pd.read_csv(os.path.abspath(out + '/test' + str(i) + '/haplotype2/' + chromosome +'.repetitions.bed'), sep='\t')

			m1=[]
			n1=[]


			for b1,b2,m,n in zip(tabh1['Start'], tabh1['End'], tabh1['Repeated Motif'], tabh1['Repetitions Number']):


				if abs(b1-start) <= acceptable or abs(b2-(end-1))<= acceptable:

					if m in possible_rotations(motif):

						m1.append(m)

					else:

						m1.append('wrong')

					n1.append(n)


			m2=[]
			n2=[]


			for b1,b2,m,n in zip(tabh2['Start'], tabh2['End'], tabh2['Repeated Motif'], tabh2['Repetitions Number']):

				if abs(b1-start) <= acceptable or abs(b2-(end-1)) <= acceptable:

					if m in possible_rotations(motif):

						m2.append(m)

					else:

						m2.append('wrong')
					
					n2.append(n)


			outtab=pd.DataFrame({'H1N':n1,'H2N':n2, 'H1M': m1, 'H2M':m2, 'TN':[number], 'TM':[motif]})
			outtab.to_csv(os.path.abspath(out+ '/test' + str(i) + '/result.tsv'),sep='\t',index=False)


			progress.current += 1

			progress()

			i+=1

		except:


			subprocess.call(['rm', '-r', os.path.abspath(out + '/simfasta' + str(i))])
			subprocess.call(['rm', '-r', os.path.abspath(out + '/simbam' + str(i))])
			subprocess.call(['rm', '-r', os.path.abspath(out + '/test' + str(i))])

			logging.info('Wrong minimap2 mapping or very few repetitions to detect. Repeat')

	
	out=os.path.abspath('expansions')


	if not os.path.exists(out):

		os.makedirs(out)

	i=1

	while i < number_of_simulations+1:

		logging.info('Expansions. Iteration: ' + str(i))

		#expansions then

		proceed=False

		while not proceed:

			subprocess.call(['bash', 'expandknown.sh', str(size)], stderr=open(os.devnull, 'wb'))

			with open(os.path.abspath('random.region.bed'), 'r') as infoin:

				for line in infoin:

					chromosome=line.split('\t')[0]
					start=int(line.split('\t')[1])
					end=int(line.split('\t')[2])
					motif=line.split('\t')[3].split('x')[1].split('\n')[0]
					number=int(line.split('\t')[3].split('x')[0])

					if start - 50000 < 0:

						proceed=False

					else:

						proceed=True

		with open(os.path.abspath('sim.bed'), 'w') as infoout:


			infoout.write(chromosome + '\t' + str(start-50000) + '\t' + str(start + 50000) + '\t' + str(100.0) + '\t' + str(100.0) + '\n')


		if not os.path.exists(os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa')):

			with open(os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), 'w') as fout:

				subprocess.call(['samtools', 'faidx', os.path.abspath(reference), chromosome], stdout=fout)


		subprocess.call(['VISOR', 'HACk', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('h1.bed'), '-o', os.path.abspath(out + '/simfasta' + str(i) + '/h1')])
		subprocess.call(['VISOR', 'HACk', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('h2.bed'), '-o', os.path.abspath(out + '/simfasta' + str(i) + '/h2')])


		done=False

		while not done:
		
			subprocess.call(['VISOR', 'LASeR', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('sim.bed'), '-s', os.path.abspath(out + '/simfasta' + str(i) + '/h1'), '-c', str(coveragemin*2), '-th', str(7), '--noaddtag','-o', os.path.abspath(out + '/simbam' + str(i) + '/h1')])
			subprocess.call(['VISOR', 'LASeR', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('sim.bed'), '-s', os.path.abspath(out + '/simfasta' + str(i) + '/h2'), '-c', str(coveragemin*2), '-th', str(7), '--noaddtag','-o', os.path.abspath(out + '/simbam' + str(i) + '/h2')])


			with open(os.path.abspath('test.bed'), 'w') as outtest:

				outtest.write(chromosome + '\t' + str(start-500) + '\t' + str(start + 500) + '\n')


			counter1,bam1=check_coverage(pysam.AlignmentFile(os.path.abspath(out + '/simbam' + str(i) + '/h1/sim.srt.bam')), chromosome, start-500, start+500, coveragemin, coveragemax)
			counter2,bam2=check_coverage(pysam.AlignmentFile(os.path.abspath(out + '/simbam' + str(i) + '/h2/sim.srt.bam')), chromosome, start-500, start+500, coveragemin, coveragemax)


			if counter1==0 or counter2 == 0:

				done=True
				continue


			if bam1 and bam2:

				done=True

			else:

				subprocess.call(['rm', '-r', os.path.abspath(out + '/simbam' + str(i))])
				logging.info('Coverage outside wanted range for current iteration. Repeat')

		if counter1==0 or counter2 == 0:

			subprocess.call(['rm', '-r', os.path.abspath(out + '/simfasta' + str(i))])
			subprocess.call(['rm', '-r', os.path.abspath(out + '/simbam' + str(i))])
			continue


		totlen=len(motif)*(number-2)

		if totlen > 15:

			totlen=15

		subprocess.call(['TRiCoLOR', 'REFER', '-g', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'), '-bed', os.path.abspath('test.bed'), '-bam', os.path.abspath(out + '/simbam' + str(i) + '/h1/sim.srt.bam'), os.path.abspath(out + '/simbam' + str(i) + '/h2/sim.srt.bam'), '-o', os.path.abspath(out + '/test' + str(i)), '-t', str(3), '-m', str(len(motif)), '--precisemotif', '-s', str(totlen)])

		try:

			tabh1=pd.read_csv(os.path.abspath(out + '/test' + str(i) + '/haplotype1/' + chromosome +'.repetitions.bed'), sep='\t')
			tabh2=pd.read_csv(os.path.abspath(out + '/test' + str(i) + '/haplotype2/' + chromosome +'.repetitions.bed'), sep='\t')

			m1=[]
			n1=[]


			for b1,b2,m,n in zip(tabh1['Start'], tabh1['End'], tabh1['Repeated Motif'], tabh1['Repetitions Number']):


				if abs(b1-start) <= acceptable or abs(b2-(end-1))<= acceptable:

					if m in possible_rotations(motif):

						m1.append(m)

					else:

						m1.append('wrong')

					n1.append(n)


			m2=[]
			n2=[]


			for b1,b2,m,n in zip(tabh2['Start'], tabh2['End'], tabh2['Repeated Motif'], tabh2['Repetitions Number']):

				if abs(b1-start) <= acceptable or abs(b2-(end-1)) <= acceptable:

					if m in possible_rotations(motif):

						m2.append(m)

					else:

						m2.append('wrong')
					
					n2.append(n)


			outtab=pd.DataFrame({'H1N':n1,'H2N':n2, 'H1M': m1, 'H2M':m2, 'TN':[number], 'TM':[motif]})
			outtab.to_csv(os.path.abspath(out+ '/test' + str(i) + '/result.tsv'),sep='\t',index=False)


			progress.current += 1

			progress()

			i+=1

		except:


			subprocess.call(['rm', '-r', os.path.abspath(out + '/simfasta' + str(i))])
			subprocess.call(['rm', '-r', os.path.abspath(out + '/simbam' + str(i))])
			subprocess.call(['rm', '-r', os.path.abspath(out + '/test' + str(i))])

			logging.info('Wrong minimap2 mapping or very few repetitions to detect. Repeat')

	
	os.remove(os.path.abspath('random.region.bed'))
	os.remove(os.path.abspath('sim.bed'))
	os.remove(os.path.abspath('test.bed'))
	os.remove(os.path.abspath('h1.bed'))
	os.remove(os.path.abspath('h2.bed'))




def Collect(size,mode,allowed):


	folders=glob.glob(os.path.abspath(mode) + '/test*')

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

	sensitivity=TP/(TP+FN)
	specificity=TN/(TN+FP)

	precision=TP/(TP+FP)
	recall=TP/(TP+FN)

	exact_core=EC/(EC+SC+WC)
	shifted_core=SC/(EC+SC+WC)
	wrong_core=WC/(EC+SC+WC)

	return (sensitivity,specificity), (precision,recall), (exact_core,shifted_core,wrong_core)



if __name__=='__main__':

	main()
