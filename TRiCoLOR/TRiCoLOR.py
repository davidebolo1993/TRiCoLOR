#!/usr/bin/python env

#python 3 standard library

import argparse
from argparse import HelpFormatter


def main():

	parser = argparse.ArgumentParser(prog='TRiCoLOR', description='''TRiCoLOR: Tandem Repeats Caller fOr LOng Reads''', epilog='''This program was developed by Davide Bolognini and Tobias Rausch at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='SENSoR, REFER, ApP') #three submodules

	
	## SENSoR ##

	parser_sensor = subparsers.add_parser('SENSoR', help='Shannon ENtropy ScanneR. Scan haplotype-resolved BAM, calculate Shannon entropy along chromosomes and output putative repetitive regions in BED')

	required = parser_sensor.add_argument_group('Required I/O arguments')

	required.add_argument('-b', '--bamfile', help='one or two haplotype-resolved BAM', metavar='BAM', nargs='+', action='append', required=True)
	required.add_argument('-o', '--output', metavar='folder', help='output folder',required=True)

	algorithm = parser_sensor.add_argument_group('Parameters for BAM scanning')

	algorithm.add_argument('-s', '--scansize', type=int, help='scansize (#bps) for BAM scanning [20]', metavar='',default=20)
	algorithm.add_argument('-e', '--entropy', type=float, help='Shannon entropy treshold [1.3]', metavar='',default=1.3)
	algorithm.add_argument('-c', '--call', type=int, help='minimum number of reads supporting each entropy drop [5]', metavar='',default=5)
	algorithm.add_argument('-l', '--length', type=int, help='minimum length of the entropy drops [20]', metavar='', default=20)

	filters=parser_sensor.add_argument_group('BED for output filtering')

	filters.add_argument('-x', '--exclude', help='BED with propietary regions to exclude from SENSoR BED [None]', metavar='',default=None)  

	additionals = parser_sensor.add_argument_group('Additional parameters')

	additionals.add_argument('--label', help='label to identify the output [sample]',metavar='',default='sample')
	additionals.add_argument('-chrs', '--chromosomes', help='scan BAM only for chromosomes provided by user. If None, scan BAM using chromosomes in BAM header.', metavar='',nargs='+', action='append', default=None)

	parser_sensor.set_defaults(func=run_subtool)

	
	## REFER ## #REpeats FindER

	parser_finder = subparsers.add_parser('REFER', help='REpeats FindER. Search repetitions in regions from .bed file using a regular-expression approach modified to allow errors; outputs are .bed files with repetitions for each haplotype and reference and a standard .vcf file')

	required = parser_finder.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='.fa',required=True)
	required.add_argument('-bed','--bedfile', help='.bed file generated with SENSoR or proprietary .bed file in the same format containing regions that identify putative tandem repetitions', metavar='.bed',required=True)
	required.add_argument('-bam1','--bamfile1', help='haplotype-resolved .bam file, first haplotype',metavar='.bam',required=True)
	required.add_argument('-bam2','--bamfile2', help='haplotype-resolved .bam file, second haplotype',metavar='.bam',required=True)
	required.add_argument('-O','--output', help='where the results will be saved',metavar='folder',required=True)

	algorithm = parser_finder.add_argument_group('Regex-search parameters')

	algorithm.add_argument('-m','--motif', type=int, help='minimum size of the motif of the repetition: use 0 or 1 for any [0]',metavar='',default=0)
	algorithm.add_argument('-t','--times', type=int, help='consencutive times a repetition must occur at least to be detected: use 0 for any [3]',metavar='',default=3)
	algorithm.add_argument('-s','--size', type=int, help='minimum size the repeated period must have at least to be called [15]',metavar='',default=15)
	algorithm.add_argument('-o','--overlapping', type=str2bool, help='check for overlapping repetitions [False]',metavar='', default='False')
	algorithm.add_argument('-mm','--maxmotif', type=int, help='exclude motifs which length is greater than value [6]',metavar='', default=6)
	algorithm.add_argument('-edit','--editdistance', type=int, help='allow specified number of insertions, deletions or substitutions in repetitive patterns [1]',metavar='', default=1)

	utilities = parser_finder.add_argument_group('Coverage treshold')

	utilities.add_argument('-c', '--coverage', type=int, help='minimum number of reads to call a valid consensus sequence [5]', metavar='', default=5)
	
	additionals = parser_finder.add_argument_group('Additional argument')

	additionals.add_argument('-sname','--samplename', help='sample name to use in .vcf header [Sample]',metavar='',default='Sample')

	parser_finder.set_defaults(func=run_subtool)



	## ApP ## Alignment Plotter ##


	parser_plotter = subparsers.add_parser('ApP', help='Alignment Plotter. Generate an interactive plot highlighting repetitions found in a user-defined region; output is an .html file that can be opened using the default browser')

	required = parser_plotter.add_argument_group('Required I/O arguments')

	required.add_argument('-g', '--genome', metavar='.fa', help='reference genome', required=True)
	required.add_argument('-mbam1', '--mergedbamfile1',metavar='.bam', help='haplotype 1 merged.srt.bam file for the wanted chromosome',required=True)
	required.add_argument('-mbam2','--mergedbamfile2', metavar='.bam', help='haplotype 2 merged.srt.bam file for the wanted chromosome',required=True)
	required.add_argument('-bed','--bedfile', metavar='.bed', help='propietary .bed file with regions with chromosome, start, end and label for the wanted repetition/s',required=True)	
	required.add_argument('-O', '--output', metavar='folder', help='where the .html file/s will be saved',required=True)


	tables = parser_plotter.add_argument_group('Regions with repetitions')

	tables.add_argument('-gb', '--genomebed', metavar='', default=None, help='reference genome repetitions.bed file for the wanted chromosome. If None plot alignment but do not highlight reference repetitions [None]')
	tables.add_argument('-h1b', '--hap1bed', metavar='', default=None, help='haplotype 1 repetitions.bed file for the wanted chromosome. If None plot alignment but do not highlight haplotype 1 repetitions [None]')
	tables.add_argument('-h2b', '--hap2bed', metavar='', default=None, help='haplotype 2 repetitions.bed file for the wanted chromosome. If None plot alignment but do not highlight haplotype 2 repetitions [None]')


	parser_plotter.set_defaults(func=run_subtool)


	## RIVAL ## Repeats Illumina VALidator ???

	args = parser.parse_args()

	args.func(parser, args)


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


def str2bool(v):

	if v.lower() == 'true':

		return True

	elif v.lower() == 'false':

		return False
	
	else:

		raise argparse.ArgumentTypeError('Boolean value expected for argument TRiCoLOR REFER -o/--overlapping. Use True or False.')


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
