#!/usr/bin/python env

#python 3 standard library

import argparse
import sys
from argparse import HelpFormatter

def main():

	
	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''TRiCoLOR: Tandem Repeats Caller fOr LOng Reads''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='SENSoR, REFER, SAGE, ApP')

	## SENSoR ##

	parser_sensor = subparsers.add_parser('SENSoR', help='Shannon ENtropy ScanneR. Scan haplotype-resolved BAM, calculate Shannon entropy along chromosomes and identify putative repetitive regions')

	required = parser_sensor.add_argument_group('Required I/O arguments')

	required.add_argument('-bam', '--bamfile', help='one or two haplotype-resolved BAM', metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-o', '--output', metavar='folder', help='output folder',required=True)

	algorithm = parser_sensor.add_argument_group('Parameters for BAM scanning')

	algorithm.add_argument('-s', '--scansize', type=int, help='scansize (#bps) for BAM scanning [20]', metavar='',default=20)
	algorithm.add_argument('-e', '--entropy', type=float, help='Shannon entropy treshold [1.25]', metavar='',default=1.25) #? HAS BEEN TRAINED
	algorithm.add_argument('-c', '--call', type=int, help='minimum number of reads supporting the entropy drops [5]', metavar='',default=5)
	algorithm.add_argument('-l', '--length', type=int, help='minimum length of the entropy drops [30]', metavar='', default=30)
	algorithm.add_argument('-id', '--innerdistance', type=int, help='maximum distance range (#bps) to merge nearby intervals on the same haplotype [100]', metavar='', default=100)
	algorithm.add_argument('-od', '--outerdistance', type=int, help='maximum distance range (#bps) to merge nearby intervals on different haplotypes [100]', metavar='', default=100)


	additionals = parser_sensor.add_argument_group('Additional parameters')

	additionals.add_argument('--chromosomes', help='scan BAM only for chromosomes provided. If None, scan BAM using all chromosomes in BAM header [None]', metavar='',nargs='+', action='append', default=None)

	parser_sensor.set_defaults(func=run_subtool)

	## REFER ##

	parser_refer = subparsers.add_parser('REFER', help='REpeats FindER. Search repetitions in regions from BED using a regular-expression approach modified to allow errors')

	required = parser_refer.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='FASTA',required=True)
	required.add_argument('-bed','--bedfile', help='BED generated with SENSoR or equivalent proprietary BED containing putative repetitive regions', metavar='BED',required=True)
	required.add_argument('-bam','--bamfile', help='one or two haplotype-resolved BAM',metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-o','--output', help='output folder',metavar='folder',required=True)

	consensus = parser_refer.add_argument_group('Consensus-building parameters')

	consensus.add_argument('--match', type=int, help='reward for a matching base [5]', metavar='', default = 5)
	consensus.add_argument('--mismatch', type=int, help='penalty for a base not matching [-4]', metavar='', default = -4)
	consensus.add_argument('--gapopen', type=int, help='penalty for gap opening [-8]', metavar='', default = -8)
	consensus.add_argument('--gapextend', type=int, help='penalty for gap extending [-6]', metavar='', default = -6)

	algorithm = parser_refer.add_argument_group('Regex-search parameters')

	algorithm.add_argument('-m','--motif', type=int, help='minimum size of the repetition motif [1]',metavar='',default=1)
	algorithm.add_argument('-mm','--maxmotif', type=int, help='exclude motifs which size is greater than value [6]',metavar='', default=6)
	algorithm.add_argument('-t','--times', type=int, help='minimum number of consecutive times the motif must be repeated to be detected [3]',metavar='',default=3)
	algorithm.add_argument('-s','--size', type=int, help='minimum size the repeated region must have to be called [30]',metavar='',default=30)
	algorithm.add_argument('-e','--editdistance', type=int, help='allowed number of insertions, deletions or substitutions in repetitions [1]',metavar='',default=1)
	algorithm.add_argument('--overlapping', help='look for overlapping repeated motif', action='store_true')
	algorithm.add_argument('--precisemotif', help='coherce -m/--motif to find only repetitions with specified motif size', action='store_true')
	algorithm.add_argument('--precisetimes', help='coherce -t/--times to find only repetitions occuring specified number of times', action='store_true')	


	utilities = parser_refer.add_argument_group('Coverage treshold')

	utilities.add_argument('-c', '--coverage', type=int, help='minimum number of reads to call a consensus sequence for region [5]', metavar='', default=5)
	
	additionals = parser_refer.add_argument_group('Additional parameters')

	additionals.add_argument('--samplename', help='sample name in BCF header [sample]', metavar='',default='sample')
	additionals.add_argument('--readstype', help='long reads technology (ONT, PB) [ONT]', metavar='', default='ONT', choices=['ONT', 'PB'])
	additionals.add_argument('--threads', help='number of cores [1]',metavar='', default=1, type=int)

	parser_refer.set_defaults(func=run_subtool)

	## SAGE ##

	parser_sage = subparsers.add_parser('SAGE', help='SAmple GEnotyper. Derive the rough genotype of repetitive regions from individuals related to the one genotyped with REFER')

	required = parser_sage.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='FASTA',required=True)
	required.add_argument('-bcf','--bcffile', help='BCF generated with REFER', metavar='BCF',required=True)
	required.add_argument('-bam', '--bamfile', help='one or more comma-separated couples of haplotype-resolved BAM from individuals related to the one already genotyped with REFER', dest='bamfile', metavar='BAM', type=BAM, nargs='+', required=True)
	required.add_argument('-o','--output', help='output folder',metavar='folder',required=True)

	consensus = parser_sage.add_argument_group('Consensus-building parameters')

	consensus.add_argument('--match', type=int, help='reward for a matching base [5]', metavar='', default = 5)
	consensus.add_argument('--mismatch', type=int, help='penalty for a base not matching [-4]', metavar='', default = -4)
	consensus.add_argument('--gapopen', type=int, help='penalty for gap opening [-8]', metavar='', default = -8)
	consensus.add_argument('--gapextend', type=int, help='penalty for gap extending [-6]', metavar='', default = -6)

	utilities = parser_sage.add_argument_group('Coverage treshold')

	utilities.add_argument('-c', '--coverage', type=int, help='minimum number of reads to call a consensus sequence in region [5]', metavar='', default=5)

	additionals = parser_sage.add_argument_group('Additional parameters')

	additionals.add_argument('--samplename', help='one name for each couple of BAM specified in -bam/--bamfile [SAMPLE1, ...]', metavar='',default=None, nargs='+', action='append')
	additionals.add_argument('--readstype', help='long reads technology (ONT, PB) [ONT]', metavar='', default='ONT', choices=['ONT', 'PB'])
	additionals.add_argument('--threads', help='number of cores [1]',metavar='', default=1, type=int)
	additionals.add_argument('--store', help=argparse.SUPPRESS, action='store_true')

	parser_sage.set_defaults(func=run_subtool)

	## ApP ## Alignment Plotter ##

	parser_app = subparsers.add_parser('ApP', help='Alignment Plotter. Generate an interactive HTML highlighting alignments and repetitions')

	required = parser_app.add_argument_group('Required I/O arguments')

	required.add_argument('-g', '--genome', metavar='FASTA', help='reference genome', required=True)
	required.add_argument('-bam','--bamfile', help='one or two consensus BAM generated with REFER',metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-bed','--bedfile', metavar='BED', help='propietary BED (CHROM,START,END,LABEL) with regions to plot',required=True)	
	required.add_argument('-o', '--output', metavar='folder', help='output folder',required=True)

	tables = parser_app.add_argument_group('BED with repetitions to highlight')

	tables.add_argument('-gb', '--genomebed', metavar='', default=None, help='BED generated with REFER with repetitions in reference[None]')
	tables.add_argument('-hb', '--haplotypebed', metavar='', default=None, help='one or more ordered BED generated with REFER with repetitions in BAM to -bam/--bamfile [None]',nargs='+', action='append')
	
	parser_app.set_defaults(func=run_subtool)

	args = parser.parse_args()

	args.func(parser, args)


## CLASS


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


## FUNCTIONS


def BAM(s):


    try:

        x,y= map(str, s.split(','))
        
        return x,y
    
    except:

        raise argparse.ArgumentTypeError('BAM files to SAGE -bam/--bamfile must be given as couples of BAM files: BAM1h1,BAM1h2 BAM2h1,BAM2h2 ...')


def run_subtool(parser, args):


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
