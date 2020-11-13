#!/usr/bin/python env

#python 3 standard library

import sys
import argparse
from argparse import HelpFormatter
import os

def main():

	#TRiCoLOR v1.1. Added logo. Changes to each module are commented below.

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''TRiCoLOR: Tandem Repeats Caller fOr LOng Reads''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 
	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='SENSoR, REFER, SAGE, ApP')

	## SENSoR ##

	#changes v.1.0 --> v1.1
	#v1.1 accept haplotype-tagged BAM
	#v1.1 removed calls to subprocesses (bedtools still required in PATH)
	#v1.1 printing messages instead of creating .log file
	#v1.1 added .gz compression for final BED (TRiCoLOR.srt.bed.gz)

	parser_sensor = subparsers.add_parser('SENSoR', help='Shannon ENtropy ScanneR. Scan haplotype-resolved/tagged BAM, calculate Shannon entropy along chromosomes and identify putative repetitive regions')

	required = parser_sensor.add_argument_group('Required I/O arguments')

	required.add_argument('-bam', '--bamfile', help='a couple of splitted BAM haplotypes or a single HP-tagged BAM', metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-o', '--output', help='output folder', metavar='folder', required=True, type=str)

	algorithm = parser_sensor.add_argument_group('Parameters for BAM scanning')

	algorithm.add_argument('-s', '--scansize', type=int, help='scansize (#bp) for BAM scanning [20]', metavar='',default=20)
	algorithm.add_argument('-e', '--entropy', type=float, help='Shannon entropy treshold [1.23]', metavar='',default=1.23) # HAS BEEN TRAINED
	algorithm.add_argument('-c', '--call', type=int, help='minimum number of reads supporting the entropy drops [5]', metavar='',default=5)
	algorithm.add_argument('-l', '--length', type=int, help='minimum length of the entropy drops [30]', metavar='', default=30)
	algorithm.add_argument('--inner', type=int, help='maximum distance range (#bp) to merge nearby intervals on the same haplotype [100]', metavar='', default=100)
	algorithm.add_argument('--outer', type=int, help='maximum distance range (#bp) to merge nearby intervals on different haplotypes [100]', metavar='', default=100)

	additionals = parser_sensor.add_argument_group('Additional parameters')

	additionals.add_argument('--chromosomes', help='scan BAM only for chromosomes provided. If None, scan BAM using all chromosomes in BAM header [None]', metavar='',nargs='+', action='append', default=None)
	additionals.add_argument('--exclude', help='text file containing chromosomes (one per line) to exclude from analysis. Is highly recommended to exclude decoy chromosomes, if any [None]', metavar='', default=None)

	parser_sensor.set_defaults(func=run_subtool)

	## REFER ##

	#changes v.1.0 --> v1.1
	#v1.1 accept haplotype-tagged BAM
	#v1.1 re-alignment uses asm10 preset to align consensus sequences back to reference
	#v1.1 removed calls to all subprocesses but consensus.cpp (alignment is internally performed with mappy and bedtools still required in PATH)
	#v1.1 printing messages instead of creating .log file
	#v1.1 added .gz compression for final VCF/BED (TRiCoLOR.srt.vcf.gz/TRiCoLOR.srt.bed.gz)

	parser_refer = subparsers.add_parser('REFER', help='REpeats FindER. Search repetitions in regions from BED using a regular-expression approach modified to allow errors')

	required = parser_refer.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='FASTA',required=True, type=str)
	required.add_argument('-bed','--bedfile', help='BED generated with SENSoR or equivalent proprietary BED containing putative repetitive regions', metavar='BED',required=True, type=str)
	required.add_argument('-bam','--bamfile', help='a couple of splitted BAM haplotypes or a single HP-tagged BAM',metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-o','--output', help='output folder',metavar='folder',required=True, type=str)

	consensus = parser_refer.add_argument_group('Consensus-building parameters')

	consensus.add_argument('--match', type=int, help='reward for a matching base [5]', metavar='', default = 5)
	consensus.add_argument('--mismatch', type=int, help='penalty for a base not matching [-4]', metavar='', default = -4)
	consensus.add_argument('--gapopen', type=int, help='penalty for gap opening [-8]', metavar='', default = -8)
	consensus.add_argument('--gapextend', type=int, help='penalty for gap extending [-6]', metavar='', default = -6)

	algorithm = parser_refer.add_argument_group('Regex-search parameters')

	algorithm.add_argument('-m','--motif', type=int, help='minimum size of the repetition motif [1]',metavar='',default=1)
	algorithm.add_argument('-mm','--maxmotif', type=int, help='exclude motifs which size is greater than value [6]',metavar='', default=6)
	algorithm.add_argument('-t','--times', type=int, help='minimum number of consecutive times the motif must be repeated to be detected [5]',metavar='',default=5)
	algorithm.add_argument('-s','--size', type=int, help='minimum size the repeated region must have to be called [50]',metavar='',default=50)
	algorithm.add_argument('-e','--editdistance', type=int, help='allowed number of insertions, deletions or substitutions in repetitions [1]',metavar='',default=1)
	algorithm.add_argument('--overlapping', help='look for overlapping repeated motif', action='store_true')
	algorithm.add_argument('--precisemotif', help='coherce -m/--motif to find only repetitions with specified motif size', action='store_true')
	algorithm.add_argument('--precisetimes', help='coherce -t/--times to find only repetitions occuring specified number of times', action='store_true')	

	utilities = parser_refer.add_argument_group('Coverage and soft clipping tresholds')

	utilities.add_argument('-c', '--coverage', type=int, help='minimum number of reads to call a consensus sequence for region [5]', metavar='', default=5)
	utilities.add_argument('-sc', '--softclipping', type=float, help='maximum percentage of soft-clipped bases allowed in consensus sequences [40.0]', metavar='', default=40.0)
	
	additionals = parser_refer.add_argument_group('Additional parameters')

	additionals.add_argument('--samplename', help='sample name in VCF header [SAMPLE]', metavar='',default='SAMPLE', type=str)
	additionals.add_argument('--threads', help='number of cores [1]',metavar='', default=1, type=int)
	additionals.add_argument('--mmidir', help='if provided, store (and read) minimap2 indexes for the processed chromosomes into a different directory [None]', metavar='', default=None, type=str)
	additionals.add_argument('--exclude', help='text file containing chromosomes (one per line) to exclude from analysis. Is highly recommended to exclude decoy chromosomes, if any [None]', metavar='', default=None)

	parser_refer.set_defaults(func=run_subtool)
	
	## ApP ## Alignment Plotter ##
	#changes v.1.0 --> v1.1
	#v1.1 code slightly reworked
	#v1.1 use a single region instead of a BED for plotting
	#v1.1 printing messages instead of creating .log file

	parser_app = subparsers.add_parser('ApP', help='Alignment Plotter. Generate an interactive HTML highlighting alignments and repetitions')

	required = parser_app.add_argument_group('Required I/O arguments')

	required.add_argument('-g', '--genome', metavar='FASTA', help='reference genome', required=True)
	required.add_argument('-bam','--bamfile', help='consensus BAM generated with REFER for haplotype 1 and haplotype 2',metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-o', '--output', metavar='folder', help='output folder',required=True)
	required.add_argument('region', metavar='REGION', help='region to visualize in samtools format (CHROM:START-END)')

	tables = parser_app.add_argument_group('BED with repetitions to highlight')

	tables.add_argument('-gb', '--genomebed', metavar='', default=None, help='BED generated by REFER with repetitions in reference [None]')
	tables.add_argument('-h1b', '--haplotype1bed', metavar='', default=None, help='BED generated by REFER with repetitions in haplotype1 [None]')
	tables.add_argument('-h2b', '--haplotype2bed', metavar='', default=None, help='BED generated by REFER with repetitions in haplotype2 None]')

	parser_app.set_defaults(func=run_subtool)

	## SAGE ##
	#changes v.1.0 --> v1.1
	#v1.1 accept couples of haplotype-tagged BAM
	#v1.1 as this performs a rough comparison based on edit-distance, now this module skip the re-alignment step
	#v1.1 removed calls to all subprocesses but consensus.cpp
	#v1.1 printing messages instead of creating .log file

	parser_sage = subparsers.add_parser('SAGE', help='SAmple GEnotyper. Derive the rough genotype of repetitive regions from individuals related to the one genotyped with refer')

	required = parser_sage.add_argument_group('Required I/O arguments')

	required.add_argument('-vcf','--vcffile', help='VCF generated with REFER', metavar='BCF',required=True)
	required.add_argument('-bam', '--bamfile', help='comma-separated couples of haplotype-resolved BAM or a couple of HP-tagged BAM from individuals related to the one genotyped with REFER', dest='bamfile', metavar='BAM', type=BAM, nargs='+', required=True)
	required.add_argument('-o','--output', help='output folder',metavar='folder',required=True)

	consensus = parser_sage.add_argument_group('Consensus-building parameters')

	consensus.add_argument('--match', type=int, help='reward for a matching base [5]', metavar='', default = 5)
	consensus.add_argument('--mismatch', type=int, help='penalty for a base not matching [-4]', metavar='', default = -4)
	consensus.add_argument('--gapopen', type=int, help='penalty for gap opening [-8]', metavar='', default = -8)
	consensus.add_argument('--gapextend', type=int, help='penalty for gap extending [-6]', metavar='', default = -6)

	utilities = parser_sage.add_argument_group('Coverage treshold')

	utilities.add_argument('-c', '--coverage', type=int, help='minimum number of reads to call a consensus sequence in region [5]', metavar='', default=5)

	additionals = parser_sage.add_argument_group('Additional parameters')

	additionals.add_argument('--mendel', help='check mendelian consistency of genotyped repetitions', action='store_true')
	additionals.add_argument('--samplename', help='one name for each BAM sample [None]', metavar='', default=None, nargs='+', action='append')
	additionals.add_argument('--threads', help='number of cores [1]',metavar='', default=1, type=int)

	parser_sage.set_defaults(func=run_subtool)
	
	#print help if no subcommand nor --help provided

	print(r"""

	 ______ ______  __  ______  ______  __      ______  ______    
	/\__  _/\  == \/\ \/\  ___\/\  __ \/\ \    /\  __ \/\  == \   
	\/_/\ \\ \  __<\ \ \ \ \___\ \ \/\ \ \ \___\ \ \/\ \ \  __<   
	   \ \_\\ \_\ \_\ \_\ \_____\ \_____\ \_____\ \_____\ \_\ \_\ 
	    \/_/ \/_/ /_/\/_/\/_____/\/_____/\/_____/\/_____/\/_/ /_/ v1.1

	""")
	
	if len(sys.argv)==1:
    	
		parser.print_help(sys.stderr)
		sys.exit(1)

	
	#case-insensitive submodules
	
	if sys.argv[1].lower() == "sensor":

		sys.argv[1] = "SENSoR"

	elif sys.argv[1].lower() == "refer":

		sys.argv[1] = "REFER"

	elif sys.argv[1].lower() == "sage":

		sys.argv[1] = "SAGE"

	elif sys.argv[1].lower() == "app":

		sys.argv[1] = "ApP"

	args = parser.parse_args()
	args.func(parser, args)


class CustomFormat(HelpFormatter):

	'''
	Custom help format
	'''

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


def BAM(s):

	'''
	Deal with HP-tagged/splitted BAM for SAGE
	'''

	try:

		x,y= map(str, s.split(','))
		return os.path.abspath(x),os.path.abspath(y)
	
	except:

		return os.path.abspath(s)


def run_subtool(parser, args):

	'''
	Choose which command to run
	'''

	if args.command == 'SENSoR': #Shannon ENtropy ScanneR

		from .SENSoR import SENSoR as submodule
	
	elif args.command == 'REFER': #REpeats FindER

		from .REFER import REFER as submodule

	elif args.command == 'SAGE': #SAmple GEnotyper

		from .SAGE import SAGE as submodule

	elif args.command == 'ApP': #Alignment Plotter

		from .ApP import ApP as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':


	main()
