#!/usr/bin/python


import pandas as pd


File = "/home/davideb/HG00733.calls_plus_loci.bed" #bed file

Table = pd.read_table(File, delimiter='\t')

Chromosome = Table['#chrom']
Start = Table['tStart']
Hap =  Table['hap']
End = Table['tEnd']
SVtype = Table['svType']
SvAnn = Table['svAnn']

FakeChrom = []
FakeStart=[]
FakeEnd=[]
FakeHap = []
FakeSVtype = []
FakeSvAnn = []

## add fake regions in which no alterations were identified
## minimum length of these regions if 1000 bp (if Start-End differs of exactly 10000 bp)


for i in range(len(Start)-1):
	if (Start[i+1]-End[i])>100000:
		NewStart = End[i]+44500
		NewEnd = Start[i+1]-44500
		FakeChrom.append(Chromosome[i])
		FakeStart.append(NewStart)
		FakeEnd.append(NewEnd)
		FakeHap.append("Fake")
		FakeSVtype.append("Fake")
		FakeSvAnn.append("Fake")


DataFake = {'Chromosome': FakeChrom, 'Start': FakeStart, 'End':FakeEnd , 'Hap' : FakeHap, 'SVtype':FakeSVtype, 'SVann':FakeSvAnn}
FakeBed=pd.DataFrame(data=DataFake,columns=['Chromosome','Start','End','Hap','SVtype','SVann'])
FakeBedNOXY=FakeBed[FakeBed['Chromosome'].str.contains("^chr[0-9]+$")] #remove X/Y chromosomes



DataTrue = {'Chromosome': Chromosome, 'Start': Start, 'End':End , 'Hap' : Hap, 'SVtype': SVtype, 'SVann':SvAnn}
True_Bed=pd.DataFrame(data=DataTrue,columns=['Chromosome','Start','End','Hap','SVtype','SVann'])
TrueBedDeletions=True_Bed[True_Bed['SVtype'].str.contains("deletion") & True_Bed['SVann'].str.contains("TandemRepeat") & True_Bed['Chromosome'].str.contains("^chr[0-9]+$")] #given bed filtered for chr1-chr22, deletion and TandemRepeats type of alterations


chrID = ['chr{}'.format(x) for x in list(range(1,23))] #Use ID to re-order chromosome

Frames = [TrueBedDeletions,FakeBedNOXY]
ConcatBed = pd.concat(Frames)
ConcatBed['Chromosome'] = pd.Categorical(ConcatBed['Chromosome'], chrID)
TestBed=ConcatBed.sample(1000) #random sample of regions with deletions and regions without deletions

TestBed.sort_values(by=['Chromosome', 'Start'],inplace=True,ascending=True)

TestBed.to_csv("/home/davideb/TestBed.txt",sep="\t",index=False,header=False)
