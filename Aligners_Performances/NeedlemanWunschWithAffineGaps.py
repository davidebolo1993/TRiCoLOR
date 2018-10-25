#!/usr/bin/python

#Pure Python, does not need any additional module

def QuickMatrixPrinter(list): #as this code works with lists, this is a function that allows to print list as they were tables, just to have a look at the disposition of the values inside them.

	s = [[str(e) for e in row] for row in list]
	lens = [max(map(len, col)) for col in zip(*s)]
	fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
	table = [fmt.format(*row) for row in s]
	print('\n'.join(table))



def MatchMismatch(list_of_nucleotides, matchscore=1, mismatchscore=-3): #given the penalty for match and mismathc, returns a dictionary with the score (nested dictionary)

	d = {} 

	for nuc in list_of_nucleotides:

		d[nuc] = dict()

		for nuc2 in list_of_nucleotides:

			d[nuc].update({nuc2:matchscore} if nuc2==nuc else {nuc2:mismatchscore}) #create nested dictionary

	return d



def IsDNA(inseq):

	return set(list(inseq)) <= set(["A","C","G","T"])


def IsDNAString(inseq):

	try:

		assert(isinstance(inseq,str))

	except:

		raise TypeError("input is not string")

	try:

		assert(IsDNA(inseq))

	except:

		raise ValueError("input is not DNA")


#initialize the three matrix needed for affine gap dynamic programming as it is explained here: http://pages.cs.wisc.edu/~bsettles/ibs08/lectures/02-alignment.pdf


#initialize matrix M

def Initialize_M(i, j): #M must have 0 in top left and "-inf" in top row and left column

	if j == 0 and i == 0:

		return 0

	else:

		if j == 0 or i == 0:

			return minimum_value

		else:

			return 0


#initialize matrix Ix

def Initialize_Ix(i,j):

	if i > 0 and j == 0:

		return minimum_value

	else:

		if j > 0:

			return h_penalty+g_penalty*j

		else:

			return h_penalty


#initialize matrix Iy

def Initialize_Iy(i,j):

	if j > 0 and i == 0:

		return minimum_value

	else:

		if i > 0:

			return h_penalty+g_penalty*i

		else:

			return h_penalty



def MatchMismatchValue(seq1, seq2, i, j): #return the value in dict for match and mismatch

	return substitution_matrix[seq2[i-1]][seq1[j-1]]
    


def Backtrace(seq1,seq2,Ix,Iy,M):


	a = ''
	b = ''
	i = len(seq2)
	j = len(seq1)
	while (i>0 or j>0):
		if (i>0 and j>0 and M[i][j] == M[i-1][j-1] + MatchMismatchValue(seq1, seq2, i, j)):
			a += seq1[j-1]
			b += seq2[i-1]
			i -= 1; j -= 1
		elif (i>0 and M[i][j] == Iy[i][j]):
			a += '-'
			b += seq2[i-1]
			i -= 1
		elif (j>0 and M[i][j] == Ix[i][j]):
			a += seq1[j-1]

			b += '-'
			j -= 1

	a = ' '.join([a[j] for j in range(-1, -(len(a)+1), -1)])
	b = ' '.join([b[j] for j in range(-1, -(len(b)+1), -1)])

	return [a, b]



def AffineGapNeedlemanWunsch(seq1, seq2, h_penalty, g_penalty, substitution_matrix): #h is penalty for open, g is penalty for extending

	IsDNAString(seq1)
	IsDNAString(seq2)

	row = len(seq2) + 1
	col = len(seq1) + 1


	#initialize matrixes
	
	M = [[Initialize_M(i, j) for j in range(0, col)] for i in range(0, row)] 
	#QuickMatrixPrinter(M)
	Iy=[[Initialize_Iy(i, j) for j in range(0, col)] for i in range(0, row)]
	#QuickMatrixPrinter(Iy)
	Ix=[[Initialize_Ix(i, j) for j in range(0, col)] for i in range(0, row)]
	#QuickMatrixPrinter(Ix)	

	for j in range(1, col):
		for i in range(1, row):
			Ix[i][j] = max((h_penalty + g_penalty + M[i][j-1]), (g_penalty + Ix[i][j-1]), (h_penalty + g_penalty + Iy[i][j-1]))
			Iy[i][j] = max((h_penalty + g_penalty + M[i-1][j]), (h_penalty + g_penalty + Ix[i-1][j]), (g_penalty + Iy[i-1][j]))
			M[i][j] = max(MatchMismatchValue(seq1, seq2, i, j) + M[i-1][j-1], Ix[i][j], Iy[i][j])

	#all the matrixes have been generated
	#QuickMatrixPrinter(M)
	#QuickMatrixPrinter(Iy)
	#QuickMatrixPrinter(Ix)

	#retrive optiamal alignment

	res=Backtrace(seq1,seq2,Ix,Iy,M)

	return res




##################################################################################################### EXAMPLE #####################################################################################################

seq1="ATCGTTTTTTAGCTGATCGTAGCTG"
seq2="ATGTGAC"

h_penalty=-10 #penalty for opening a gap
g_penalty=-1 #penalty for extending a gap
list_of_nucleotides=["A","C","T","G"]
substitution_matrix=MatchMismatch(list_of_nucleotides,matchscore=1,mismatchscore=-3)
minimum_value=float("-inf")

#changing the penalty for opening-extending of course change the alignment
#try with inverting the h_penalty and g_penalty

alignment=AffineGapNeedlemanWunsch(seq1, seq2, h_penalty, g_penalty, substitution_matrix)

