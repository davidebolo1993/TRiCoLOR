#!/usr/bin/python

import argparse
from argparse import HelpFormatter

import os
import glob
import math
import random
import pysam
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns


def main():


	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''Simulations for TRiCoLOR program''', epilog='''This program was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	parameters = parser.add_argument_group('Arguments')

	parameters.add_argument('-f', '--folder', help='One folder containing simulations (folder containing expansions and contractions subfolders generated with the simulator)', metavar='folder', required=True)
	parameters.add_argument('-s', '--size', help='Scansize to train on simulated data [20]', metavar='', default=20, type=int)

	args = parser.parse_args()

	path=os.path.abspath(args.folder)
	size=args.size

	Get_Entropy_From_Simulations(path,size)




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


def BamScanner(bamfile,scansize):

	BamFile=pysam.AlignmentFile(bamfile,'rb')

	ent__=[]

	for read in BamFile.fetch():

		if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:

			sequence=read.seq

			start=0
			end=scansize
			terminal_ind=len(sequence)-1

			while terminal_ind > end:

				if terminal_ind-end >= scansize:

					ent__.append(entropy(sequence[start:end]))
					start+=scansize
					end+=scansize

				else:

					break

	return ent__



def Get_Entropy_From_Simulations(path,scansize):


	expansions=glob.glob(os.path.abspath(path + '/expansions') + '/simbam*')
	ebams=[glob.glob(os.path.abspath(i+'/h1')+'/*.srt.bam') for i in expansions]
	
	contractions=glob.glob(os.path.abspath(path + '/contractions') + '/simbam*')
	cbams=[glob.glob(os.path.abspath(i+'/h1')+'/*.srt.bam') for i in contractions]

	normals=glob.glob(os.path.abspath(path + '/contractions') + '/simbam*')
	nbams=[glob.glob(os.path.abspath(i+'/h2')+'/*.srt.bam') for i in normals]

	
	#Expansions

	E=[]

	for bam in ebams:

		en_=BamScanner(bam[0],scansize)		
					
		E.append(en_)


	for a in E:

		x=[el+1 for el,x in enumerate(a)]

	E_dist=[]

	for el in E:

		E_dist.extend(el[0:-1])

	E_MS=mean(E_dist) - 3 *stddev(E_dist)

	E_dist_downsampled=random.sample(E_dist, 10000)

	ax=sns.distplot(E_dist_downsampled, rug=True, rug_kws={'color': 'green'}, kde_kws={'color': 'olive', 'lw': 3}, hist_kws={'histtype': 'step', 'linewidth': 3,'alpha': 1, 'color': 'darkgreen'})
	data_x, data_y = ax.lines[0].get_data()
	Y_E_MS = np.interp(E_MS,data_x, data_y)
	plt.xlabel('Entropy')
	plt.ylabel('Density')
	plt.plot([E_MS],[Y_E_MS],marker='o',color='black',markersize=7)
	plt.title('Entropy distplot for BAM with microsatellites expansions')
	plt.savefig(os.path.abspath(path + '/Distplot_BAM_expansions.pdf'))
	plt.savefig(os.path.abspath(path + '/Distplot_BAM_expansions.png'))
	plt.show()



	#Normals

	N=[]

	for bam in nbams:

		en_=BamScanner(bam[0],scansize)		
					
		N.append(en_)


	for a in N:

		x=[el+1 for el,x in enumerate(a)]


	N_dist=[]

	for el in N:

		N_dist.extend(el[0:-1])

	N_MS=mean(N_dist) - 3 *stddev(N_dist)

	N_dist_downsampled=random.sample(N_dist, 10000)


	ax=sns.distplot(N_dist_downsampled, rug=True, rug_kws={'color': 'blue'}, kde_kws={'color': 'teal', 'lw': 3}, hist_kws={'histtype': 'step', 'linewidth': 3,'alpha': 1, 'color': 'darkblue'})
	data_x, data_y = ax.lines[0].get_data()
	Y_N_MS = np.interp(N_MS,data_x, data_y)
	plt.xlabel('Entropy')
	plt.ylabel('Density')
	plt.plot([N_MS],[Y_N_MS],marker='o',color='black',markersize=7)
	plt.title('Entropy distplot for BAM without microsatellites modifications')
	plt.savefig(os.path.abspath(path + '/Distplot_BAM_normals.pdf'))
	plt.savefig(os.path.abspath(path + '/Distplot_BAM_normals.png'))
	plt.show()



	#Contractions

	C=[]

	for bam in cbams:

		en_=BamScanner(bam[0],scansize)		
					
		C.append(en_)


	for a in C:

		x=[el+1 for el,x in enumerate(a)]

	C_dist=[]

	for el in C:

		C_dist.extend(el[0:-1])

	C_dist_downsampled=random.sample(C_dist, 10000)


	C_MS=mean(C_dist) - 3 *stddev(C_dist)

	ax=sns.distplot(C_dist_downsampled, rug=True, rug_kws={'color': 'red'}, kde_kws={'color': 'indianred', 'lw': 3}, hist_kws={'histtype': 'step', 'linewidth': 3,'alpha': 1, 'color': 'darkred'})
	data_x, data_y = ax.lines[0].get_data()
	Y_C_MS = np.interp(C_MS,data_x, data_y)
	plt.xlabel('Entropy')
	plt.ylabel('Density')
	plt.plot([C_MS],[Y_C_MS],marker='o',color='black',markersize=7)
	plt.title('Entropy distplot for BAM with microsatellites contractions')
	plt.savefig(os.path.abspath(path + '/Distplot_BAM_contractions.pdf'))
	plt.savefig(os.path.abspath(path + '/Distplot_BAM_contractions.png'))
	plt.show()

	#Total

	T=[]

	T.extend(N_dist)
	T.extend(E_dist)
	T.extend(C_dist)

	T_MS=mean(T) - 3 *stddev(T)

	T_dist_downsampled=random.sample(T, 10000)

	ax=sns.distplot(T_dist_downsampled, rug=True, rug_kws={'color': 'lightgray'}, kde_kws={'color': 'darkgray', 'lw': 3}, hist_kws={'histtype': 'step', 'linewidth': 3,'alpha': 1, 'color': 'dimgray'})
	data_x, data_y = ax.lines[0].get_data()
	Y_T_MS = np.interp(T_MS,data_x, data_y)
	plt.xlabel('Entropy')
	plt.ylabel('Density')
	plt.plot([T_MS],[Y_T_MS],marker='o',color='black',markersize=7)
	plt.title('Entropy distplot for all BAM')
	plt.savefig(os.path.abspath(path + '/Distplot_BAM_all.pdf'))
	plt.savefig(os.path.abspath(path + '/Distplot_BAM_all.png'))
	plt.show()


	#write data


	with open(os.path.abspath(path + '/Summary.txt'), 'w') as fout:

		fout.write('Contractions, Mean-3SD: ' + str(C_MS) + '\n')
		fout.write('Expansions, Mean-3SD: ' + str(E_MS) + '\n')
		fout.write('Normals, Mean-3SD: ' + str(N_MS) + '\n')
		fout.write('Total, Mean-3SD: ' + str(T_MS) + '\n')



if __name__ == '__main__':

	main()


