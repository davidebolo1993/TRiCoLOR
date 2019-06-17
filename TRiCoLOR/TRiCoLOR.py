#!/usr/bin/python env

#python 3 standard library

import argparse
from argparse import HelpFormatter


def main():

	
	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''TRiCoLOR: Tandem Repeats Caller fOr LOng Reads''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='SENSoR, REFER, ApP')


	## SENSoR ##

	parser_sensor = subparsers.add_parser('SENSoR', help='Shannon ENtropy ScanneR. Scan haplotype-resolved BAM, calculate Shannon entropy along chromosomes and output putative repetitive regions in BED')

	required = parser_sensor.add_argument_group('Required I/O arguments')

	required.add_argument('-bam', '--bamfile', help='one or two haplotype-resolved BAM', metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-o', '--output', metavar='folder', help='output folder',required=True)

	algorithm = parser_sensor.add_argument_group('Parameters for BAM scanning')

	algorithm.add_argument('-s', '--scansize', type=int, help='scansize (#bps) for BAM scanning [20]', metavar='',default=20)
	algorithm.add_argument('-e', '--entropy', type=float, help='Shannon entropy treshold [1.25]', metavar='',default=1.25) #? HAS BEEN TRAINED
	algorithm.add_argument('-c', '--call', type=int, help='minimum number of reads supporting the entropy drops [5]', metavar='',default=5)
	algorithm.add_argument('-l', '--length', type=int, help='minimum length of the entropy drops [30]', metavar='', default=30)

	filters=parser_sensor.add_argument_group('BED for output filtering')

	filters.add_argument('-x', '--exclude', help='BED with propietary regions to exclude from SENSoR BED [None]', metavar='',default=None)  

	additionals = parser_sensor.add_argument_group('Additional parameters')

	additionals.add_argument('-chrs', '--chromosomes', help='scan BAM only for chromosomes provided. If None, scan BAM using all chromosomes in BAM header [None]', metavar='',nargs='+', action='append', default=None)

	parser_sensor.set_defaults(func=run_subtool)


	## REFER ##

	parser_finder = subparsers.add_parser('REFER', help='REpeats FindER. Search repetitions in regions from BED using a regular-expression approach modified to allow errors; output BED with repetitions for haplotype/s-reference and a standard VCF')

	required = parser_finder.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='FASTA',required=True)
	required.add_argument('-bed','--bedfile', help='BED generated with SENSoR or equivalent proprietary BED containing putative repetitive regions', metavar='BED',required=True)
	required.add_argument('-bam','--bamfile', help='one or two haplotype-resolved BAM',metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-o','--output', help='output folder',metavar='folder',required=True)

	algorithm = parser_finder.add_argument_group('Regex-search parameters')

	algorithm.add_argument('-m','--motif', type=int, help='minimum size of the repetition motif [1]',metavar='',default=1)
	algorithm.add_argument('-mm','--maxmotif', type=int, help='exclude motifs which size is greater than value [6]',metavar='', default=6)
	algorithm.add_argument('-t','--times', type=int, help='minimum number of consecutive times the motif must be repeated to be detected [3]',metavar='',default=3)
	algorithm.add_argument('--overlapping', help='look for overlapping repeated motif', action='store_true')
	algorithm.add_argument('--precisemotif', help='coherce -m/--motif to find only repetitions with specified motif size', action='store_true')
	algorithm.add_argument('--precisetimes', help='coherce -t/--times to find only repetitions occuring specified number of times', action='store_true')	
	algorithm.add_argument('-s','--size', type=int, help='minimum size the repeated region must have to be called [15]',metavar='',default=15)
	algorithm.add_argument('-edit','--editdistance', type=int, help='allowed number of insertions, deletions or substitutions in repetition [1]',metavar='',default=1)

	utilities = parser_finder.add_argument_group('Coverage treshold')

	utilities.add_argument('-c', '--coverage', type=int, help='minimum number of reads to call a consensus sequence in region [5]', metavar='', default=5)
	
	additionals = parser_finder.add_argument_group('Additional parameters')

	additionals.add_argument('--samplename', help='sample name in VCF header [sample]',metavar='',default='sample')
	additionals.add_argument('--readstype', help='long reads technology (ONT, PB) [ONT]',metavar='',default='ONT', choices=['ONT', 'PB'])
	additionals.add_argument('-th', '--threads', help='number of cores to use [1]',metavar='',default=1, type=int)

	parser_finder.set_defaults(func=run_subtool)


	## ApP ## Alignment Plotter ##


	parser_plotter = subparsers.add_parser('ApP', help='Alignment Plotter. Generate an interactive plot that highlights repetitions and alignments in region; output HTML')

	required = parser_plotter.add_argument_group('Required I/O arguments')

	required.add_argument('-g', '--genome', metavar='FASTA', help='reference genome', required=True)
	required.add_argument('-bam','--bamfile', help='one or two consensus BAM generated with REFER',metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-bed','--bedfile', metavar='BED', help='propietary BED (CHROM,START,END,LABEL) with regions to plot',required=True)	
	required.add_argument('-o', '--output', metavar='folder', help='output folder',required=True)

	tables = parser_plotter.add_argument_group('BED with repetitions to highlight')

	tables.add_argument('-gb', '--genomebed', metavar='', default=None, help='BED for repetitions in reference [None]')
	tables.add_argument('-hb', '--haplotypebed', metavar='', default=None, help='one or more ordered BED for repetitions in BAM to -bam/--bamfile [None]',nargs='+', action='append')
	
	parser_plotter.set_defaults(func=run_subtool)

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


## FUNCTION


def run_subtool(parser, args):


	if args.command == 'SENSoR': #Shannon ENtropy ScanneR

		from .SENSoR import SENSoR as submodule
	
	elif args.command == 'REFER': #REpeats FindER

		from .REFER import REFER as submodule

	elif args.command == 'ApP': #Alignment Plotter

		from .ApP import ApP as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':


	main()
