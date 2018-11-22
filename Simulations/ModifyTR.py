#!/usr/bin/python

import os
import sys
import csv
import pyfaidx
import random
import argparse

#generate a modified reference fasta for the chosen chromosome with a repeat expansion/contraction

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("--reference", help="reference file")
	parser.add_argument("--bed", help="bed file containing all the microsatellites for classic chromosomes")
	parser.add_argument("--chromosome", help="chromosome to be modified")
  	parser.add_argument("--mode", help="expand or contract the repetition?")
	parser.add_argument("--position", type=int, help="get the nearest repetition to the specified position")
	parser.add_argument("--size", type=int, help="How many repetitions do you want to remove/add?")
	parser.add_argument("--mismatch",default=False, help="Change a nucleotide with one of the other three, randomly")
	parser.add_argument("--howmany", type=int, default=None, help="How many mismatches have to be generated?")
	parser.add_argument("--label", help="label for the fasta containing the repetition")
	parser.add_argument("--output", help="name of the directory where the result will be saved")
	args = parser.parse_args()


	if args.mode == "contraction":

		Contract_Repetitions(args.reference,args.bed,args.chromosome,args.position,args.size,args.mismatch,args.howmany,args.label,args.output)

	elif args.mode == "expansion":

		Expand_Repetitions(args.reference,args.bed,args.chromosome,args.position,args.size,args.mismatch,args.howmany,args.label,args.output)

	else:

		sys.exit("accepted mode parameters are contraction and expansion")



def TakeClosest(point,start_list):

	closest_start = min(start_list, key=lambda x:abs(x-point))
	ind_closest = [i for i,el in enumerate(start_list) if el==closest_start][0]

	return ind_closest


def ExtractSequence(genome,chromosome,start,end):

	ref=pyfaidx.Fasta(genome)
	chrom=ref[chromosome]
	seq=chrom[:len(chrom)].seq
	bef_=seq[:start]
	wanted=seq[start:end-1] #pyfaidx wau to get the true coordinates
	aft_=seq[end-1:]

	return bef_,wanted,aft_



def Random_Change_N_Char(word,howmany):

	length = len(word)
	word = list(word)
	k = random.sample(range(0,length),howmany)
	k.sort()
	nuc_list=["A","T","C","G"]
   
	for index in k:

		add_rem=word[index]
		nuc_list.remove(add_rem)
		word[index] = ''.join(random.sample(nuc_list,k=1))
		nuc_list.append(add_rem)
   
	return("" . join(word))



def FastaGenerator(header,seq,label,pathout):

	dirname=os.path.abspath(pathout)

	if not os.path.exists(dirname):

		os.makedirs(dirname)

	fname=os.path.abspath(dirname+"/"+label+".fa")
	f=open(fname,"w")
	f.write(header + "\n" + seq)
	f.close()



def Contract_Repetitions(genome,bed,chromosome,deletion_approximate_point,restriction_size, mismatch,howmany,label,pathout):

	chr_start=[]
	chr_end=[]
	num_rep=[]
	repetition=[]


	with open(bed) as tsv:

		for line in csv.reader(tsv, delimiter="\t"):

			if line[0] == chromosome:

				chr_start.append(int(line[1]))
				chr_end.append(int(line[2]))
				num_,rep_=line[3].rsplit("x")
				num_rep.append(int(num_))
				repetition.append(rep_)

	ind_closest=TakeClosest(deletion_approximate_point,chr_start)
	before,seq,after=ExtractSequence(genome,chromosome,chr_start[ind_closest],chr_end[ind_closest])
	occurences=num_rep[ind_closest]
	new_occurences=occurences-restriction_size
	new_seq=repetition[ind_closest]*new_occurences


	if not mismatch:

		seq_to_write=before+new_seq+after
		header=">" + chromosome + ":"+str(chr_start[ind_closest]+1)+"_rep:"+new_seq+"_core:"+repetition[ind_closest]+"_before:"+str(num_rep[ind_closest])+"_now:"+str(new_occurences)+"_mismatch="+str(mismatch)+"_number_of_mismatches=None"

	else:

		mod_seq=Random_Change_N_Char(new_seq,howmany)
		seq_to_write=before+mod_seq+after
		header=">" + chromosome + ":"+str(chr_start[ind_closest]+1)+"_rep:"+mod_seq+"_core:"+repetition[ind_closest]+"_before:"+str(num_rep[ind_closest])+"_now:"+str(new_occurences)+"_mismatch="+str(mismatch)+"_number_of_mismatches=" + str(howmany)

	FastaGenerator(header,seq_to_write,label,pathout)




def Expand_Repetitions(genome,bed,chromosome,deletion_approximate_point,restriction_size, mismatch,howmany,label,pathout):

	chr_start=[]
	chr_end=[]
	num_rep=[]
	repetition=[]


	with open(bed) as tsv:

		for line in csv.reader(tsv, delimiter="\t"):

			if line[0] == chromosome:

				chr_start.append(int(line[1]))
				chr_end.append(int(line[2]))
				num_,rep_=line[3].rsplit("x")
				num_rep.append(int(num_))
				repetition.append(rep_)

	ind_closest=TakeClosest(deletion_approximate_point,chr_start)
	before,seq,after=ExtractSequence(genome,chromosome,chr_start[ind_closest],chr_end[ind_closest])
	occurences=num_rep[ind_closest]
	new_occurences=occurences+restriction_size
	new_seq=repetition[ind_closest]*new_occurences


	if not mismatch:

		seq_to_write=before+new_seq+after
		header=">" + chromosome + ":"+str(chr_start[ind_closest]+1)+"_rep:"+new_seq+"_core:"+repetition[ind_closest]+"_before:"+str(num_rep[ind_closest])+"_now:"+str(new_occurences)+"_mismatch="+str(mismatch)+"_number_of_mismatches=None"

	else:

		mod_seq=Random_Change_N_Char(new_seq,howmany)
		seq_to_write=before+mod_seq+after
		header=">" + chromosome + ":"+str(chr_start[ind_closest]+1)+"_rep:"+mod_seq+"_core:"+repetition[ind_closest]+"_before:"+str(num_rep[ind_closest])+"_now:"+str(new_occurences)+"_mismatch="+str(mismatch)+"_number_of_mismatches=" + str(howmany)

	FastaGenerator(header,seq_to_write,label,pathout)






if __name__ == main():
   
	main()
