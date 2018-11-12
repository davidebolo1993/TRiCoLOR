#!/usr/bin/python

import sys
import os
import glob
import pysam
import subprocess #for now, just calling Alfred from within Python: in the future, a Python wrapper can be written
import gzip
import shutil


def Sub_None(list_of_coord):

	return [-999999 if v is None else v for v in list_of_coord] #return a very small number so that it won't be take into account


def Mean_Qual(list_of_quals):

	return float(sum(list_of_quals)) / len(list_of_quals)


def Get_Start_Ind(list_of_coord,start):

	if start in list_of_coord:

		return [i for i in range(len(list_of_coord)) if list_of_coord[i] == start]

	else:

		fake_start = min(Sub_None(list_of_coord), key=lambda x:abs(x-start))
		return [i for i in range(len(list_of_coord)) if list_of_coord[i] == fake_start]


def Get_End_Ind(list_of_coord,end):

	if end in list_of_coord:

		return [i for i in range(len(list_of_coord)) if list_of_coord[i] == end]

	elif any(i >= end for i in Sub_None(list_of_coord)): #check if any coord is higher in read than end, else return None. After, None will be filtered out

		fake_end = min(Sub_None(list_of_coord), key=lambda x:abs(x-end))
		return [i for i in range(len(list_of_coord)) if list_of_coord[i] == fake_end]

	else:

		return


def Assert_Size(start,end,size): #probably more assertations are needed

	try:

		assert end-start >= size

	except:

		raise ValueError('chosen size is greater than end-start difference')


def Extract_Info(pysam_AlignedSegment):


	identifier=pysam_AlignedSegment.query_name
	mapping_qual=pysam_AlignedSegment.mapping_quality
	read_qual=Mean_Qual(pysam_AlignedSegment.query_qualities)
	identifier_2=identifier+"_"+str(mapping_qual)+"_"+str(read_qual)

	if pysam_AlignedSegment.is_reverse:

		identifier_2=identifier+"_"+str(mapping_qual)+"_"+str(read_qual)+"_rev"

	else:

		identifier_2=identifier+"_"+str(mapping_qual)+"_"+str(read_qual)+"_for"

	seq=pysam_AlignedSegment.seq #give the correct result even for reverse strand, no need to reverse-complement

	return (identifier_2,seq)


def Extract_Coord(pysam_AlignedSegment,start,end):


	a_coord=pysam_AlignedSegment.get_reference_positions(full_length=True) #list that indicates the reference position on which each base alignes (None included for insertion/soft clipped bases)
	identifier=pysam_AlignedSegment.query_name
	mapping_qual=pysam_AlignedSegment.mapping_quality
	read_qual=Mean_Qual(pysam_AlignedSegment.query_qualities)
	identifier2=identifier+"_"+str(mapping_qual)+"_"+str(read_qual)

	if pysam_AlignedSegment.is_reverse:

		identifier_2=identifier+"_"+str(mapping_qual)+"_"+str(read_qual)+"_rev"

	else:

		identifier_2=identifier+"_"+str(mapping_qual)+"_"+str(read_qual)+"_for"

	ind_start=Get_Start_Ind(a_coord,start)
	ind_end=Get_End_Ind(a_coord,end)

	if ind_end is None:

		return

	else:

		return (identifier_2,ind_start[0],ind_end[0])


def Bamfile_Analyzer(bamfile,chromosome,start,end,size): #sliding window function to retrive informations


	Assert_Size(start,end,size)

	BamFile=pysam.AlignmentFile(bamfile,"rb")

	i=0 #counter
	pointer=start+1 #slicing window

	inf=[] #list for appending results
	coord=[] #list for appending coordinates

	while end-1 >= start :

		inf.append([])
		coord.append([])

		for read in BamFile.fetch(chromosome,start,pointer):

			inf[i].append(Extract_Info(read))

			if end > start+size:

				coord[i].append(Extract_Coord(read,start,start+size))

			else:

				coord[i].append(Extract_Coord(read,start,end))

		if end-start > size:

			start += size
			pointer += size

		elif end-start <= size and start < end-1:

			start = end-1
			pointer = end

		else: #break is not needed, but, just in case...

			break

		i+=1 #update counter in order to store results progressively


	return inf,coord[:len(coord)-1] #number of item in list is equal to number of iteration for inf, equal to the number of times size is included in end-start window for coord


def Reorder(l1,l2): #reorder and filtering

	l2_2=list(filter(None.__ne__,l2)) #filtering None for reads that do not have an end in sequence
	l2_el=[el[0] for el in l2_2]
	mapping = dict((x[0], x[1:]) for x in l1)
	l1[:] = [(x,) + mapping[x] for x in l2_el]

	return l1,l2_2


def Com_Seq(info_list,coord_list): # for each

	filtered_seq=[]
	filtered_coord=[]

	for i in range(len(info_list)-1):

		fip=info_list[i]
		sip=info_list[i+1]
		id_fip=[a for a,b in fip]
		id_sip=[a for a,b in sip]
		int_id=list(set(id_fip).intersection(id_sip))
		com_seq=[(a,b) for a,b in fip if a in int_id]
		l1,l2=Reorder(com_seq,coord_list[i]) # reordering seems not to be necessary but allows filtering too

		filtered_seq.append(l1)
		filtered_coord.append(l2)

	return filtered_seq,filtered_coord


def Fasta_Generator(filtered_seq,filtered_coord,out):

	for i in range(len(filtered_seq)):

		dirname=os.path.abspath(out)

		if not os.path.exists(dirname):

			os.makedirs(dirname)

		fname=os.path.abspath(dirname+"/Fasta"+str(i+1)+".unaligned.fa")

		f=open(fname,"a")

		for j in range(len(filtered_seq[i])):

			fasta_in_window=">"+filtered_seq[i][j][0]+"\n"+filtered_seq[i][j][1][filtered_coord[i][j][1]:filtered_coord[i][j][2]+1]+"\n" #get sequence from start to end (end included)
			f.write(fasta_in_window)

		f.close() #save many .fasta file in Directory



def MA(out,mmi_ref): #multiple alignment function


    fafile=glob.glob(os.path.abspath(out)+"/*.unaligned.fa")

    outfa=[]
    outall=[]

    for fa in fafile:

        outfa.append(fa.replace("unaligned.fa","consensus.fa.gz"))
        outall.append(fa.replace("unaligned.fa","alignment.txt.gz"))

        for fastain,consensusout,alignmentout in zip(fafile,outfa,outall):

            subprocess.call(['/home/bolognin/Downloads/alfred/src/alfred', 'consensus', '-a', alignmentout, '-c', consensusout, '-f', 'fasta', '-t', 'ont', fastain])

            with gzip.open(consensusout,'rb') as f_in,open(consensusout.replace(".fa.gz",".fa"), 'wb') as f_out:

                shutil.copyfileobj(f_in, f_out)

            with open(consensusout.replace(".fa.gz",".sam"), 'w') as f:

                subprocess.call(['minimap2', '-ax', 'map-ont', mmi_ref, consensusout.replace(".fa.gz",".fa")], stdout=f)

            with open(consensusout.replace(".fa.gz",".bam"), 'w') as f:

                    subprocess.call(['htsbox', 'samview', '-bS', consensusout.replace(".fa.gz",".sam")], stdout=f)

            with open(consensusout.replace(".fa.gz",".srt.bam"), 'w') as f:

                subprocess.call(['samtools', 'sort', consensusout.replace(".fa.gz",".bam")], stdout=f)

            subprocess.call(['samtools', 'index', consensusout.replace(".fa.gz",".srt.bam")])
