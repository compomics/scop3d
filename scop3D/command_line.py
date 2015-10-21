#__author__ = 'elien'
#scop3D 2.0
#situation: protein + multiple sequences as input

from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from weblogolib import *
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
from Bio import SeqIO
from Bio.PDB import *

import sys
import shutil
import os
import subprocess
import math
from math import log

#multiple sequence alignment
#1. Retrieval of data input
print('!The analysis can be performed with maximum 500 sequences!')
myProtMultInput = raw_input('Enter a filename: ')
#   test to see if the file was loaded
with open(myProtMultInput, 'r') as f:
    first_line = f.readline()
    print('The first line of your input is printed as verification that the correct file is used in the analysis:')
    print(first_line)
myMuscleInput = myProtMultInput

myCount = 0
myInputanalyze = open(myProtMultInput, 'r')
myInputlines = myInputanalyze.readlines()
myInputanalyze.close()
for myInputline in myInputlines:
    if myInputline.startswith('>'):
        myCount = myCount + 1
print('number of sequences in analysis: '+str(myCount))

#2. Multiple sequence alignment with Muscle
#       Reference: Edgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.
#                   Nucleic Acids Res 32(5), 1792-1797
print('You have to specify your muscle location:')
myMuscleLoc = raw_input('Enter your muscle location: ')
print(myMuscleLoc)
myMuscleExe = myMuscleLoc+'\muscle.exe'
print(myMuscleExe)
myAlignmentOutput = 'output_multiple_sequence_alignment.txt'
muscle_cline = MuscleCommandline(myMuscleExe, input=myMuscleInput, out=myAlignmentOutput, clwstrict=True)
muscle_cline()
alignment = open(myAlignmentOutput, "r")
print_alignment = alignment.read()
print print_alignment
alignment.close()

#3. Determination of similarity and identity between the different sequences
os.mkdir('temp_sec')
os.mkdir('pairwise_alignments')
#for each sequence a fasta-file will be created with the aid of seqretsplit from emboss
#   Rice P, Longden I, Bleasby A (2000) EMBOSS: The European Molecular Biology Open Software Suite.
#   Trends in Genetics 16(6), 276-277
subprocess.call(['C:\mEMBOSS\seqretsplit.exe','-nofeature', '-sequence', myProtMultInput, '-nofirstonly', '-auto', '-osdirectory_outseq', 'temp_sec'])

myFastafiles = os.listdir('temp_sec')
myFastaIter = iter(myFastafiles)

for myFasta in myFastaIter:
    if myFasta.endswith('fasta'):
        print(myFasta)
        myFastab = 'temp_sec\\'+myFasta
        for myFasta2 in os.listdir('temp_sec'):
            myFasta2b = 'temp_sec\\'+myFasta2

#a pairwise alignment among all sequences is performed with the aid of needle
#   Needleman SB, Wunsch CD (1970) A general method applicable to the search for similarities in the amino acid sequence of two proteins
#   J Mol Biol 48(3), 443-453
            subprocess.call(['C:\mEMBOSS\\needle.exe', '-asequence', myFastab, '-bsequence', myFasta2b, '-noendweight', '-endopen', '10.0', '-endextend', '0.5', '-brief', '-aformat', 'srspair', '-auto', '-aname_outfile', myFasta+'_'+myFasta2, '-adirectory_outfile', 'pairwise_alignments'])

        os.remove(myFastab)
        print(myFastafiles)

os.rmdir('temp_sec')

#similarity and identity between the different sequences is retrieved
myIdent = open('identity','w')
myIdentrel = open('identity_percentage','w')
mySim = open('similarity','w')
mySimrel = open('similarity_percentage','w')

myIdentdict = {}
myIdentreldict = {}
mySimdict = {}
mySimreldict = {}

myAlignfiles = os.listdir('pairwise_alignments')
for myAlign in myAlignfiles:
    if myAlign.endswith('.needle'):
        print(myAlign)
        myAlignb = 'pairwise_alignments\\'+myAlign
        myAlignfile = open(myAlignb, 'r')
        mySeq1 = 'ND'
        mySeq2 = 'ND'
        for myLine in myAlignfile:
            if '# 1: ' in myLine:
                mySeq1 = myLine[5:-1]
            if '# 2: ' in myLine:
                mySeq2 = myLine[5:-1]
            if '# Identity' in myLine:
                mySplit = myLine.split()
                myIdentity = mySplit[2]
                myIdentityrel = mySplit[3][1:-1]
            if '# Similarity' in myLine:
                mySplit = myLine.split()
                mySimilarity = mySplit[2]
                mySimilarityrel = mySplit[3][1:-1]
        print('sequences: '+mySeq1+' and '+mySeq2+' identity: '+myIdentity+' identity percent: '+myIdentityrel+' similarity: '+mySimilarity+' similarity percent: '+mySimilarityrel)
        myIdentdict[mySeq1+'\t'+mySeq2] = myIdentity
        myIdentreldict[mySeq1+'\t'+mySeq2] = myIdentityrel
        mySimdict[mySeq1+'\t'+mySeq2] = mySimilarity
        mySimreldict[mySeq1+'\t'+mySeq2] = mySimilarityrel
        myAlignfile.close()

for myKey in myIdentdict:
    myIdent.write(myKey+'\t'+myIdentdict[myKey]+'\n')
for myKey in myIdentreldict:
    myIdentrel.write(myKey+'\t'+myIdentreldict[myKey]+'\n')
for myKey in mySimdict:
    mySim.write(myKey+'\t'+mySimdict[myKey]+'\n')
for myKey in mySimreldict:
    mySimrel.write(myKey+'\t'+mySimreldict[myKey]+'\n')

myIdent.close()
myIdentrel.close()
mySim.close()
mySimrel.close()

#4. Determination of the consensus sequence
#   Rice P, Longden I, Bleasby A (2000) EMBOSS: The European Molecular Biology Open Software Suite.
#   Trends in Genetics 16(6), 276-277
print("calculation of the consensus sequence")
subprocess.call(['C:\mEMBOSS\cons.exe', '-sequence', myAlignmentOutput, '-plurality', '1', '-identity', '0', '-auto', '-outseq', 'consensus', '-name', 'consensus'])
myConsensus = open('consensus', "r")
print_consensus = myConsensus.readlines()
print(print_consensus)
shutil.copy2('consensus', 'consensus.fasta')
myCons = ""
for myLine in print_consensus:
    if myLine.startswith('>'):
        myCons = ""
    else:
        myLine = myLine[:-1]
        myCons = "".join((myCons, myLine))
myConsensus.close()

#5. Calculation of the abundance matrices - biopython module
myMatrixinput = AlignIO.read(myAlignmentOutput, "clustal")
mySummary_align = AlignInfo.SummaryInfo(myMatrixinput)
myFrequency_matrix = mySummary_align.pos_specific_score_matrix(myCons)
myFrequency_matrix_str = str(myFrequency_matrix)

#Abundance matrix in absolute values
print(myFrequency_matrix_str)
myAbunout = open('abundance_matrix.txt', 'w')
myAbunout.write(myFrequency_matrix_str)
myAbunout.close()

#Abundance matrix in percentage
myAbunout2 = open('abundance_matrix.txt', 'r')
myLines = myAbunout2.readlines()
myAbunout2.close()

myAbunperout = open('abundance_matrix_percentage.txt', 'w')
for myLine in myLines:
    mySplit = myLine.split()
    if myLines.index(myLine) == 0:
        myLine = myLine[:-1]
        myAbunperout.write('\t'+"\t".join(mySplit))
    else:
        myLine = myLine[:-1]
        for mySplitel in mySplit:
            if mySplit.index(mySplitel) == 0:
                myAbunperout.write('\n'+mySplit[0]+'\t')
            else:
                myPercent = (float(mySplitel) * 100) / float(myCount)
                #print(myPercent)
                myAbunperout.write(str(myPercent)+'\t')
myAbunperout.close()

#6. Calculation of the entropy
#retrieval of the input needed (abundance_matrix.txt)
myAbinput = open('Abundance_matrix.txt', 'r')
myAbinputlines = myAbinput.readlines()
myAbinput.close()

myEntropyout = open('entropy.txt', 'w')

print('line with residues:'+myAbinputlines[0])
myAAlist = myAbinputlines[0].split()
print(myAAlist)
myN = len(myAAlist)
print(str(myN))

for myAbinputline in myAbinputlines:
    mySum = 0
    if myAbinputline == myAbinputlines[0]:
        print("already processed: "+ myAbinputline)
    else:
        print(myAbinputline)
        myEsplit = myAbinputline.split()
        for myE in myEsplit:
            if myE == myE[0]:
                print('residue: '+ myE)
                myEntropyout.write(myE)
            else:
                if myE != '0.0':
                    myPi = (float(myE) / float(myCount))
                    print('pi = '+str(myPi))
                    mySum = mySum + (myPi * math.log(myPi,2))
                    print("mySum = "+str(mySum))
        mySum = -float(mySum)
        myEntropy = (float(mySum) / math.log(float(myN),2)) * 100
        print('entropy: '+str(myEntropy))
        if myEntropy == -0.0:
            myEntropy = '0.00'
            print(myEntropy)
            myEntropyout.write('\t'+myEntropy+'\n')
        else:
            myRound = round(myEntropy, 2)
            print(myRound)
            myEntropyout.write('\t'+str(myRound)+'\n')

myEntropyout.close()

#retrieval of the percent value for the consensus residue
myAbperinput = open('abundance_matrix_percentage.txt', 'r')
myAbperinputlines = myAbperinput.readlines()
myAbperinput.close()

myConsAbperout = open('consensus_frequency_mostabundant.txt', 'w')
myConsVarperout = open('variation.txt', 'w')

for myAbperinputline in myAbperinputlines:
    if myAbperinputline == myAbperinputlines[0]:
        print(myAbperinputline)
    else:
        print(myAbperinputline)
        myApsplit = myAbperinputline.split()
        myConsAbperout.write(myApsplit[0])
        myConsVarperout.write(myApsplit[0])
        print(myApsplit[0])
        print(max(myApsplit[1:]))
        myVar = 100 - float(max(myApsplit[1:]))
        myConsAbperout.write('\t'+max(myApsplit[1:])+'\n')
        myConsVarperout.write('\t'+str(myVar)+'\n')

myConsAbperout.close()
myConsVarperout.close()

#6. Calculation of the Weblogo
#   Crooks GE, Hon G, Chandonia JM, Brenner SE (2004) WebLogo: a sequence logo generator.
#   Genome Research 14, 1188-1190

fin = open("output_multiple_sequence_alignment.txt")
seqs = read_seq_data(fin)
data = LogoData.from_seqs(seqs)
options = LogoOptions()
options.title = "sequence logo"
format = LogoFormat(data, options)
fout = open("logo.eps", "w")
eps = eps_formatter(data, format)
fout.write(eps)

myAAdict = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY', 'H':'HIS', 'I':'ILE',
			'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
			'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}

#7. Retrieval of the pdb-file
print('Do you have a pdb-file of the protein you would like to analyze? Yes/No')
myAnswer = raw_input('Your answer: ')
print(myAnswer)
myPDBfile = ()
if myAnswer == 'Yes':		#when the user has a pdb-file
	myPDBfile = raw_input('Enter your pdb-filename (without extensions): ')
	with open(myPDBfile+'.pdb', 'r') as f:
		myFirstline = f.readline()
		print(myFirstline)
elif myAnswer == 'No':		#when the user does not have a pdb-file, a blast against the pdb will be performed
	print('Blast (www) will be performed with the consensus sequence')
	mySeq = open('consensus').read()
	print(mySeq)
	print('The blast might take some time.')
	result_handle = NCBIWWW.qblast('blastp', 'pdb', mySeq)
	myBlastoutput = open('blast_results.xml', 'w')
	myBlastoutput.write(result_handle.read())
	myBlastoutput.close()
	result_handle.close()
	
	blastresult_handle = open('blast_results.xml')
	blastsummary = open('blast_summary.txt', 'w')
	blast_record = NCBIXML.read(blastresult_handle)
	global choose_pdbChoices
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			title = alignment.title
			mySplit = title.split('|')
			myPdb = mySplit[3]
			e_value = hsp.expect
			length = alignment.length
			identities = hsp.identities
			positives = hsp.positives
			gaps = hsp.gaps
			blastsummary.write("title: " + str(title) + "\n" + "score: " + str(e_value) + "\n" +
						"length: " + str(length) + "\n" +
						"identities: " + str(identities) + "\n" +
						"positives: " + str(positives) + "\n" +
						"gaps: " + str(gaps) + "\n" + myPdb + "\n\n" )
			print("title: " + str(title[1:100]) + " ..." +"\n" + "score: " + str(e_value) + "\n" +
						"length: " + str(length) + "\n" +
						"identities: " + str(identities) + "\n" +
						"positives: " + str(positives) + "\n" +
						"gaps: " + str(gaps) + "\n" + myPdb + "\n\n" )
	blastsummary.close()

	print('Which PDB-file would you like to work with?')
	print('As a rule of thumb, a structure with high sequence coverage, identity and similarity is preferred.')
	myPDBfile = raw_input('Your choice of PDB: ')

else:
	print('Your answer was neither Yes or No. Scop3D will end.')
	sys.exit()

#8. adjustment of the pdb-file
print('pdbfile which you will work with: '+myPDBfile)

#8a. selection of the chains to adjust
structure = ()

try:
	pdbl = PDB.PDBList()
	pdbl.retrieve_pdb_file(myPDBfile)
	parser = PDB.PDBParser()
	structure = parser.get_structure(myPDBfile, myPDBfile[1:3].lower()+'/pdb'+myPDBfile.lower()+'.ent')

	if myAnswer == 'No':
		io = PDBIO()
		io.set_structure(structure)
		io.save(myPDBfile+'.pdb')

except:
	print('you structure is not available in the PDB')
	parser = PDB.PDBParser()
	structure = PDB.PDBParser().get_structure(myPDBfile,myPDBfile+'.pdb')

chains = [chain for chain in structure.get_chains()]
print(chains)
ppb = PPBuilder()
for chain in chains:
		print(chain)
 		print(chain.get_id())
 		for chainseq in ppb.build_peptides(chain):
 			print(chainseq.get_sequence())
 			mySeq = chainseq.get_sequence()

 			myFastaA = '>'+myPDBfile+'_'+chain.get_id()+'\n'+chainseq.get_sequence()
 			print(myFastaA)
 			myFastaAfile = open(myPDBfile+'_'+chain.get_id()+'.fasta', 'w')
 			myFastaAfile.write('>'+myPDBfile+'_'+chain.get_id()+'\n'+str(mySeq))
 			myFastaAfile.close()

 			myFastafiles = os.listdir('.')
 			myFastaIter = iter(myFastafiles)

 			for myFasta in myFastaIter:
 				if myPDBfile+'_'+chain.get_id() in myFasta:
 					print(myFasta)

 					#a pairwise alignment among all sequences is performed with the aid of needle
 					#   Needleman SB, Wunsch CD (1970) A general method applicable to the search for similarities in the amino acid sequence of two proteins
 					#   J Mol Biol 48(3), 443-453
 					subprocess.call(['C:\mEMBOSS\\needle.exe', '-asequence', myFasta, '-bsequence', 'consensus', '-noendweight', '-endopen', '10.0', '-endextend', '0.5', '-brief', '-aformat', 'markx3', '-auto'])

myEntinput = open('entropy.txt', 'r')
myEntlines = myEntinput.readlines()
myEntinput.close()
myVarinput = open('variation.txt', 'r')
myVarlines = myVarinput.readlines()
myVarinput.close()

print('Do you know which chains have to be adjusted? Yes/No')
myAnswer = raw_input('Your answer: ')
print(myAnswer)
if myAnswer == 'Yes':
	print('Give the chains you would like to adjust. (Example: A,C,D)')
	myChainanswer = raw_input('Your answer: ')
	myChainlist = myChainanswer.split(',')
	print(myChainlist)

	for myEl in myChainlist:
		myReorganize = open(myPDBfile+'_'+myEl+'_sequences', 'w')
 		myIterator = SeqIO.parse(myPDBfile+'_'+myEl+'.needle', 'fasta')
 		mySeqlist = list(myIterator)
 		for mySeq in mySeqlist:
 			print(mySeq)
 			myReorganize.write('>'+str(mySeq.description)+'\n'+str(mySeq.seq)+'\n')
 		myReorganize.close()

 		myDatafile = open(myPDBfile+'_'+myEl+'_data.txt', 'w')
 		myReorganize = open(myPDBfile+'_'+myEl+'_sequences', 'r')
 		myReorganizelines = myReorganize.readlines()
 		myReorganize.close()
 		myDatafile.write(myReorganizelines[0][1:7]+'\t')
 		myDatafile.write(myReorganizelines[2][1:7]+'\n')
 		myListcons = list(myReorganizelines[3])
 		myListchain = list(myReorganizelines[1])
 		for i in range(len(myListchain)):
 			if myListchain[i] != '-':
				myDatafile.write(myListchain[i]+'\t')
				myDatafile.write(myListcons[i]+'\t')
				n = i
				print(n)
				if myListcons[i] == '-':
					myEntropy = '120.00'
					myDatafile.write(myEntropy+'\t')
					myVariation = '120.00'
					myDatafile.write(myVariation+'\n')
					n = i+1
				else:
					if n > len(myEntlines)-1:
						print('end of the data')
						break
					else:
						line = myEntlines[n][:-1]
						myList = line.split('\t')
						print(myList)
						print(myList[1])
						myEntropy = "%.2f" % float(myList[1])
						myDatafile.write(myEntropy+'\t')
						myLine = myVarlines[n][:-1]
						myVarlist = myLine.split('\t')
						print(myVarlist)
						myVariation = "%.2f" % float(myVarlist[1])
						myDatafile.write(myVariation+'\n')
		myDatafile.close()

else:
	myIdentdict = {}
	myFiles = os.listdir('.')
	myFilesiter = iter(myFiles)
	for myFile in myFilesiter:
		if myFile.endswith('.needle') and myFile.startswith(myPDBfile.lower()):
 			print('bestand: '+myFile)
			myChain = myFile[-8:-7].upper()
 			print('keten: '+myChain)
 			myAlignchain = open(myFile, 'r')
 			for myLine in myAlignchain:
 				if myLine.startswith('# Identity:'):
 					print('identiteitlijn: '+myLine)
 					mySplit = myLine.split('(')
 					myIdent = mySplit[1][:-3]
 					print('identiteit: '+myIdent)
 					myIdentdict[myChain] = myIdent
	print('chainsdictionary: ')
	print(myIdentdict)

	for key,val in myIdentdict.iteritems():
 		if val == max(myIdentdict.values()):
 			print(key)
 			myReorganize = open(myPDBfile+'_'+key+'_sequences', 'w')
 			myIterator = SeqIO.parse(myPDBfile+'_'+key+'.needle', 'fasta')
 			mySeqlist = list(myIterator)
 			for mySeq in mySeqlist:
 				print(mySeq)
 				myReorganize.write('>'+str(mySeq.description)+'\n'+str(mySeq.seq)+'\n')
 			myReorganize.close()

 			myDatafile = open(myPDBfile+'_'+key+'_data.txt', 'w')
 			myReorganize = open(myPDBfile+'_'+key+'_sequences', 'r')
 			myReorganizelines = myReorganize.readlines()
 			myReorganize.close()
 			myDatafile.write(myReorganizelines[0][1:7]+'\t')
 			myDatafile.write(myReorganizelines[2][1:7]+'\n')
 			myListcons = list(myReorganizelines[3])
 			myListchain = list(myReorganizelines[1])
 			for i in range(len(myListchain)):
 				if myListchain[i] != '-':
					myDatafile.write(myListchain[i]+'\t')
					myDatafile.write(myListcons[i]+'\t')
					n = i
					print(n)
					if myListcons[i] == '-':
						myEntropy = '120.00'
						myDatafile.write(myEntropy+'\t')
						myVariation = '120.00'
						myDatafile.write(myVariation+'\n')
						n = i+1
					else:
						if n > len(myEntlines)-1:
							print('end of the data')
							break
						else:
							line = myEntlines[n][:-1]
							myList = line.split('\t')
							print(myList)
							print(myList[1])
							myEntropy = "%.2f" % float(myList[1])
							myDatafile.write(myEntropy+'\t')
							myLine = myVarlines[n][:-1]
							myVarlist = myLine.split('\t')
							print(myVarlist)
							myVariation = "%.2f" % float(myVarlist[1])
							myDatafile.write(myVariation+'\n')
			myDatafile.close()

#8b. adjustment of the B-values (towards entropy and conservation)
myChainlist = []
myFiles = os.listdir('.')
myFilesiter = iter(myFiles)
for myFile in myFilesiter:
 	if myFile.endswith('_data.txt') and myFile.startswith(myPDBfile):
  		print('bestand: '+myFile)
 		myChain = myFile[-10:-9].upper()
  		print('keten: '+myChain)
 		myChainlist.extend(myChain)
print(myChainlist)

myPDBin = open(myPDBfile+'.pdb', 'r')
myLines = myPDBin.readlines()
myPDBin.close()

myPDBout = open(myPDBfile+'_entropy.pdb','w')
myPDBvarout = open(myPDBfile+'_variation.pdb','w')

for myLine in myLines:
 	if myLine.startswith('ATOM'):
  		print('line: '+myLine)
  		myChain = myLine[21:22]
  		myRespdb = myLine[17:20]
  		myBvalue = myLine[60:66]
  		myPos = myLine[22:26]
  		print('summary: '+myChain+'\t'+myRespdb+'\t'+myBvalue+'\t'+myPos)
  		myI = myLines.index(myLine)
  		print('line index: '+str(myI))

  		if myChain in myChainlist:
  			myBinfo = open(myPDBfile+'_'+myChain+'_data.txt', 'r')
  		 	myBlines = myBinfo.readlines()
  		 	myBinfo.close()

   			mySplit = myBlines[1].split('\t')
  		 	myRes = mySplit[0]
  		 	print('res from data: '+myRes)
			myEnt = mySplit[2]
  		 	myVar = mySplit[3].strip('\n').strip('\r')

  		 	if myRespdb == myAAdict[myRes.upper()]:
  				myPDBout.write((myLine[0:60]+' '+myEnt+'\n'))
				myPDBvarout.writelines((myLine[0:60]+' '+myVar+'\n'))

  				myIplus = myI + 1
  				myTest = myLines[myIplus]
 				print('control of next line: '+myTest)
  				if myTest[22:26] != myPos:
  		 	 		myBinfo = open(myPDBfile+'_'+myChain+'_data.txt', 'w')
  		 	 		counter = 0
  		 	 		for myBline in myBlines:
 						if counter == 1:
 						 	print('line to remove: '+myBline)
 							counter = counter + 1
 						else:
 							myBinfo.write(myBline)
 							counter = counter + 1
 					myBinfo.close()

  			else:
  				print('look at this!')
 				print(myRespdb+'-'+myRes+'-'+myAAdict[myRes])

  		else:
  			myPDBout.write(myLine[0:60]+'120.00'+'\n')
			myPDBvarout.write(myLine[0:60]+'120.00'+'\n')

  	else:
  		myPDBout.write(myLine)
		myPDBvarout.write(myLine)

myPDBout.close()
myPDBvarout.close()