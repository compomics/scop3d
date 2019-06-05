""" 
	Scop3D Worker - worker module for Scop3D 
	Copyright 2018 Lukasz Kreft <lukasz.kreft@vib.be>
"""
"""
 	This file is part of Scop3D.

    Scop3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Scop3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Scop3D. If not, see <http://www.gnu.org/licenses/>.
"""

import csv
import math
import os
import json
import subprocess
import tempfile
import operator
import urllib2
from lxml import etree

from Bio import AlignIO
from Bio import ExPASy
from Bio import SeqIO
from Bio import SeqUtils
from Bio import SwissProt
from Bio import Entrez
from Bio import Seq

from Bio.Align.AlignInfo import SummaryInfo
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.PDB import PDBList
from Bio.PDB import PDBIO
from Bio.PDB import Select
from Bio.PDB import Residue
from Bio.PDB import Atom
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqRecord import SeqRecord
from weblogolib import *

### Stage 1

myNoDataValue = 110
Entrez.email = "bits@vib.be"

# 0. Creating workdirs
def mkdirs(outputDir, verbose):
	try:
		os.mkdir(outputDir)
	except OSError:
		pass
	sequencesDir = tempfile.mkdtemp()
	alignmentsDir = tempfile.mkdtemp()
	chainSequencesDir = tempfile.mkdtemp()
	return (sequencesDir, alignmentsDir, chainSequencesDir)

# 1. Retrieval of data input
def downloadUniprotSequences(uniprotID, blastFile, sequencesFile, cutoff, verbose):
	print('Obtaining sequences from UniProt...')
	with open(blastFile, 'r') as f:
		records = NCBIXML.read(f)
	if verbose:
		print('Found ' + str(len(records.alignments)) + ' matches')
	with open(sequencesFile, 'w') as f:
		if uniprotID != None:
			sequence = urllib2.urlopen('http://www.uniprot.org/uniprot/' + uniprotID + '.fasta')
			f.write(sequence.read() + "\n");
		for idx, alignment in enumerate(records.alignments):
			for hsp in alignment.hsps:
				title = alignment.title
				words = title.split('|')
				if words[0] == 'gi':
			            seqID = words[3]
				elif words[0] == 'sp':
				    seqID = words[1]
				identityPercent = 100.0 * float(hsp.identities) / float(hsp.align_length)
				if (identityPercent >= float(cutoff)):
					try:
						sequence = urllib2.urlopen('http://www.uniprot.org/uniprot/' + seqID + '.fasta')
						f.write(sequence.read() + "\n")
						if verbose:
							print(seqID + " (identity " + str(identityPercent) + "% >= cutoff " + str(cutoff) + "%) - adding")
					except Exception as e:
						if verbose:
							print("WARNING: unable to download entry " + seqID + " from Uniprot: " + str(e))
							print("Trying NCBI protein...")
						handle = Entrez.efetch(db="protein", id=seqID, rettype="fasta", retmode="xml")
						erecords = Entrez.parse(handle)
						for erecord in erecords:
							r = SeqRecord( Seq.Seq(erecord['TSeq_sequence'], IUPAC.unambiguous_dna ), id=erecord['TSeq_accver'], description=erecord['TSeq_defline'])
							SeqIO.write(r, f, "fasta")
							f.write("\n")
							if verbose:
								print("OK")
						handle.close()
				else:
					if verbose:
						print(seqID + " (identity " + str(identityPercent) + "% < cutoff " + str(cutoff) + "%) - skipping")
	print('OK')

def downloadNCBISequences(blastFile, sequencesFile, cutoff, verbose):
	print('Obtaining sequences from NCBI...')
	with open(blastFile, 'r') as f:
		records = NCBIXML.read(f)
	if verbose:
		print('Found ' + str(len(records.alignments)) + ' matches')
	with open(sequencesFile, 'w') as f:
		sequences = []
		for idx, alignment in enumerate(records.alignments):
			for hsp in alignment.hsps:
				title = alignment.title
				words = title.split('|')
				seqID = words[3]
				identityPercent = 100.0 * float(hsp.identities) / float(hsp.align_length)
				if (identityPercent >= float(cutoff)):
					sequences.append(seqID);
					if verbose:
						print(seqID + " (identity " + str(identityPercent) + "% >= cutoff " + str(cutoff) + "%) - adding")
				else:
					if verbose:
						print(seqID + " (identity " + str(identityPercent) + "% < cutoff " + str(cutoff) + "%) - skipping")
		try:
			handle = Entrez.efetch(db="nuccore", id=",".join(sequences), rettype="fasta", retmode="xml")
			records = Entrez.parse(handle)
			DNAsequences = []
			for record in records:
				DNAsequences.append( SeqRecord( Seq.Seq(record['TSeq_sequence'], IUPAC.unambiguous_dna ), id=record['TSeq_accver'], description=record['TSeq_defline']) )
			SeqIO.write(DNAsequences, f, "fasta")
			handle.close()
		except Exception as e:
			print("WARNING: unable to download this entry: " + str(e))
	print('OK')

def loadSequences(sequencesFile, verbose):
	print('Loading sequences...')
	sequences = list(SeqIO.parse(sequencesFile, "fasta"))
	count = len(sequences)
	if count == 0:
		raise Exception('No sequences in the file')
	if verbose:
		print('Number of sequences in analysis: ' + str(count))
		for record in sequences:
			print("ID: " + record.id + " \t Length: " + str(len(record)))
	print('OK')
	return count

# 2. Multiple sequence alignment with Muscle
# Reference: Edgar RC (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.
# 			 Nucleic Acids Res 32(5), 1792-1797
def doAlignment(sequencesFile, alignmentFile, verbose):
	print('Performing sequence alignment...')
	params = ['muscle', '-in', sequencesFile, '-out', alignmentFile, '-fasta']
	if not verbose:
		params.append('-quiet')
	subprocess.check_call(params)
	print('OK')

# 3. Determination of similarity and identity between the different sequences
# Reference: Rice P, Longden I, Bleasby A (2000) EMBOSS: The European Molecular Biology Open Software Suite.
# 			 Trends in Genetics 16(6), 276-277
def doSplit(sequencesFile, sequencesDir, verbose):
	print('Splitting fasta file into individual sequences...')
	subprocess.check_call(['seqretsplit', '-feature', '0', '-sequence', sequencesFile, '-firstonly', '0', '-auto', '-osdirectory2', sequencesDir])
	print('OK');

# 4. A pairwise alignment among all sequences is performed with the aid of needle
# Reference: Needleman SB, Wunsch CD (1970) A general method applicable to the search for similarities in the amino acid sequence of two proteins
#			 J Mol Biol 48(3), 443-453
def doPairwiseAlignment(sequencesDir, alignmentsDir, verbose):
	print('Performing a pairwise alignment...')
	for sequenceFile1 in os.listdir(sequencesDir):
		if not sequenceFile1.endswith('fasta'):
			continue
		sequencePath1 = os.path.join(sequencesDir, sequenceFile1)
		for sequenceFile2 in os.listdir(sequencesDir):
			if not sequenceFile2.endswith('fasta'):
				continue
			if verbose:
				print("Alignment of " + sequenceFile1 + " with " + sequenceFile2)
			sequencePath2 = os.path.join(sequencesDir, sequenceFile2)
			subprocess.check_call(['needle', '-asequence', sequencePath1, '-bsequence', sequencePath2, '-noendweight', '-endopen', '10.0', '-endextend', '0.5', '-brief', '-aformat', 'srspair', '-auto', '-aname_outfile', sequenceFile1+'_'+sequenceFile2, '-adirectory_outfile', alignmentsDir])
		os.remove(sequencePath1)
	print('OK')

# 5. Similarity and identity between the different sequences is retrieved
def calcIdentitySimilarity(alignmentsDir, alignmentStatsFile, verbose):
	print('Getting identity and similarity...')
	with open(alignmentStatsFile, 'w') as stats:
		writer = csv.writer(stats, dialect='excel-tab')
		for alignmentFile in os.listdir(alignmentsDir):
			if not alignmentFile.endswith('.needle'):
				continue
			alignmentOutput = os.path.join(alignmentsDir, alignmentFile)
			alignment = AlignIO.read(alignmentOutput, "emboss")
			sequence1 = alignment[0].id
			sequence2 = alignment[1].id
			length = str(alignment.get_alignment_length())
			identity = str(alignment.annotations['identity'])
			identityPercent = str(round(100 * float(identity) / int(length), 1)) + "%"
			similarity = str(alignment.annotations['similarity'])
			similarityPercent = str(round(100 * float(similarity) / int(length), 1)) + "%"
			if verbose:
				print('Sequences: ' + sequence1 + ' and ' + sequence2 + ' - ' +
					'identity: ' + identity + '/' + length +' ( ' + identityPercent + ' ), ' +
					'similarity: ' + similarity + '/'+ length + ' ( ' + similarityPercent + ' )')
			writer.writerow([sequence1, sequence2, identity + '/' + length, identityPercent, similarity + '/' + length, similarityPercent])
	print('OK')

# 6. Determination of the consensus sequence
# Reference: Rice P, Longden I, Bleasby A (2000) EMBOSS: The European Molecular Biology Open Software Suite.
# 			 Trends in Genetics 16(6), 276-277
def calcConsensus(alignmentFile, consensusFile, verbose):
	print("Calculating of the consensus sequence...")
	subprocess.check_call(['cons', '-sequence', alignmentFile, '-plurality', '1', '-identity', '0', '-auto', '-outseq', consensusFile, '-name', 'consensus'])
	if verbose:
		consensusSeq = SeqIO.read(consensusFile, 'fasta')
		print('Consensus sequence: ')
		print(consensusSeq.seq)
	print('OK')

def translateSequences(consensusFile, translatedFile, verbose):
	print("Translating the sequence(s)...")
	sequences = SeqIO.parse(consensusFile, 'fasta')
	translatedSequences = []
	for seq in sequences:
		if verbose:
			print('> sequence: ')
			print(seq.seq)
		translation = Seq.translate(seq.seq)
		if verbose:
			print('> translation: ')
			print(translation)
		translatedSequences.append(SeqRecord( translation, id=seq.id, description=seq.description ))
	with open(translatedFile, "w") as f:
		SeqIO.write(translatedSequences, f, "fasta")
	print('OK')

# 7. Calculation of the abundance matrices - biopython module
def calcAbundance(alignmentFile, consensusFile, abundanceFile, abundancePercentFile, verbose):
	print('Calculating the abundance matrix...')
	alignment = AlignIO.read(alignmentFile, "fasta")
	summary = SummaryInfo(alignment)
	consensusSeq = SeqIO.read(consensusFile, 'fasta')
	if (len(consensusSeq) == alignment.get_alignment_length()):
		abundanceMatrix = summary.pos_specific_score_matrix(consensusSeq)
	else:
		with open(consensusFile, "w") as f:			
			SeqIO.write(SeqRecord( summary.dumb_consensus(), id='consensus'), f, "fasta")
		abundanceMatrix = summary.pos_specific_score_matrix()
	if verbose:
		print("Abundance matrix (absolute values):")
		print(str(abundanceMatrix))
	with open(abundanceFile, 'w') as f:
		f.write(str(abundanceMatrix))
	for pos, abundance in enumerate(abundanceMatrix):
		for res, value in abundance.items():
			abundanceMatrix[pos][res] = 100.0 * float(value) / float(len(alignment))
	if verbose:
		print("Abundance matrix (percentages):")
		print(str(abundanceMatrix))
	with open(abundancePercentFile, 'w') as f:
		f.write(str(abundanceMatrix))
	print('OK')

# 8. Calculation of the entropy
def calcEntropy(sequenceCount, abundanceFile, entropyFile, verbose):
	print('Calculating the entropy...')
	with open(abundanceFile, 'r') as f:
		abundanceLines = f.readlines()
	residuesCount = len(abundanceLines[0].split())
	if verbose:
		print('Residues: ' + abundanceLines[0][:-1])
		print('Total number of residues: ' + str(residuesCount))
	with open(entropyFile, 'w') as f:
		writer = csv.writer(f, dialect='excel-tab')
		writer.writerow(['residue', 'H', '%'])
		for line in abundanceLines[1:]:
			total = 0.0
			words = line.split()
			for (idx, value) in enumerate(words):
				if idx == 0:
					residue = value
				else:
					value = float(value)
					if value != 0:
						pi = float(value) / sequenceCount
						total += pi * math.log(pi, 2)
			if total != 0.0:
				total = -1 * total
			entropy = round(100 * total / math.log(residuesCount, 2), 1)
			total = str(round(total, 3))
			entropy = str(entropy)
			writer.writerow([residue, total, entropy])
			if verbose:
				print(residue + "\t entropy: " + total + " " + entropy + "%")
	print('OK')

# 9. Retrieval of the percent value for the consensus residue
def calcVariation(abundancePercentFile, consensusAbundanceVariationFile, verbose):
	print('Calculating maximal abundance and variation for consensus sequence...')
	with open(abundancePercentFile, 'r') as f:
		abundancePercentLines = f.readlines()

	with open(consensusAbundanceVariationFile, 'w') as f:
		writer = csv.writer(f, dialect='excel-tab')
		writer.writerow(['residue', 'cons %', 'var %'])
		for line in abundancePercentLines[1:]:
			words = line.split()
			residue = words[0]
			for i in range(1, len(words)): 
				words[i] = float(words[i])
			maximum = max(words[1:])
			variation = 100 - maximum
			writer.writerow([residue, round(maximum, 1), round(variation, 1)])
			if verbose:
				print(residue + "\t max abundance: " + str(round(maximum, 1)) + "\t variation: " + str(round(variation, 1)))
	print('OK')

# 10. Calculation of the Weblogo
# Reference: Crooks GE, Hon G, Chandonia JM, Brenner SE (2004) WebLogo: a sequence logo generator.
# 			 Genome Research 14, 1188-1190
def drawWeblogo(alignmentFile, logoFile, verbose):
	print('Drawing weblogo...')
	options = LogoOptions()
	options.title = "Sequence logo"
	with open(alignmentFile) as f:
		seqs = read_seq_data(f)
		data = LogoData.from_seqs(seqs)
	format = LogoFormat(data, options)
	with open(logoFile, "w") as f:
		png = png_formatter(data, format)
		f.write(png)
	print('OK')

### Stage 2

def doUniprotBlast(uniprotid, blastFile, verbose):
	print('Performing BLAST with the UniProt ID ' + uniprotid + '...')
	results = NCBIWWW.qblast('blastp', 'swissprot', uniprotid)
	with open(blastFile, 'w') as f:
		f.write(results.read())
	results.close()
	print('OK')

def doNCBIBlast(sequence, blastFile, verbose):
	print('Performing NCBI with the DNA sequence ' + sequence + '...')
	results = NCBIWWW.qblast('blastn', 'nr', sequence)
	with open(blastFile, 'w') as f:
		f.write(results.read())
	results.close()
	print('OK')


def doBlast(consensusFile, blastFile, verbose):
	print('Performing BLAST with the consensus sequence...')
	seq = open(consensusFile).read()
	results = NCBIWWW.qblast('blastp', 'pdb', seq)
	with open(blastFile, 'w') as f:
		f.write(results.read())
	results.close()
	print('OK')

def getBlastSummary(blastFile, blastSummaryFile, verbose):
	print('Analyzing BLAST results...')
	with open(blastFile, 'r') as f:
		records = NCBIXML.read(f)
	blastsummary = open(blastSummaryFile, 'w')
	print('Found ' + str(len(records.alignments)) + ' matches')
	for idx, alignment in enumerate(records.alignments):
		for hsp in alignment.hsps:
			title = alignment.title
			words = title.split('|')
			summary = (str(idx+1) + ": " + str(words[3]) + " - " + str(title[1:100]) + "...\n" +
				"score: " + str(hsp.expect) + "\t\t" + "length: " + str(alignment.length) + "\t\t" +
				"identities: " + str(hsp.identities) + "\t\t" + "positives: " + str(hsp.positives) + "\t\t" + "gaps: " + str(hsp.gaps) + "\n")
			blastsummary.write(summary)
			print(summary)
	blastsummary.close()
	return records.alignments

def getUniprotSummary(uniprotFile, blastSummaryFile, verbose):
	print('Analyzing UNIPROT entry...')
	with open(uniprotFile, 'r') as f:
		proteinJson = f.read()
	protein = json.loads(proteinJson)
	blastsummary = open(blastSummaryFile, 'w')
	if (verbose):
		print('Found ' + str(len(protein['dbReferences'])) + ' db references. Filtering on PDB entries.')
	idx = 0
	for ref in protein['dbReferences']:
		if ref['type'] == 'PDB':
			idx+=1
			summary = (str(idx) + ": " + str(ref['id']) + "\t\t" +
				"resolution: " + str(ref['properties']['resolution']) + "\t\t" + "chains: " + str(ref['properties']['chains']) + "\t\t" +
				"method: " + str(ref['properties']['method']) + "\n")
			blastsummary.write(summary)
			if (verbose):
				print(summary)
	blastsummary.close()
	print('OK')
	return idx


def choosePDB(blastSummary, outputDir):
	print('As a rule of thumb, a structure with high sequence coverage, identity and similarity is preferred.')
	pdbNr = raw_input('Your choice of PDB [ 1 - ' + str(len(blastSummary)) + ' ]: ')
	pdbTitle = blastSummary[int(pdbNr)-1].title
	pdbID = pdbTitle.split('|')[3]
	return downloadPDB(pdbID, outputDir)

def downloadPDB(pdbID, outputDir):
	print('Downloading PDB file...')
	pdb = PDBList(pdb=outputDir, obsolete_pdb=outputDir)
	pdbOutput = pdb.retrieve_pdb_file(pdbID, file_format='pdb', pdir=outputDir)
	print('OK')
	return pdbOutput

# A pairwise alignment among all sequences is performed with the aid of needle
# Reference: Needleman SB, Wunsch CD (1970) A general method applicable to the search for similarities in the amino acid sequence of two proteins
#   		 J Mol Biol 48(3), 443-453
def doChainAlignments(pdbID, structure, consensusFile, chainSequencesDir, verbose):
	print('Pairwise alignment of chain sequences with consensus...')
	chains = [chain for chain in structure.get_chains()]
	ppb = PPBuilder()
	for chain in chains:
		sequence = ""
		for pp in ppb.build_peptides(chain):
			sequence += pp.get_sequence()
		sequenceID = pdbID + '_' + chain.get_id()
		sequenceOutput = os.path.join(chainSequencesDir, sequenceID + '.fasta')
		if (len(sequence) == 0):
			print("ERROR: Unable to get chain sequence from PDB file.")
			exit("Selected PDB-file does not contain protein structure.")
		with open(sequenceOutput, 'w') as f:
			f.write('>' + sequenceID + '\n' + str(sequence))
		if verbose:
			print('Chain ' + chain.get_id())
			print(sequence)
		subprocess.check_call(['needle', '-asequence', sequenceOutput , '-bsequence', consensusFile, '-noendweight', '-endopen', '10.0', '-endextend', '0.5', '-brief', '-aformat', 'srspair', '-auto', '-aname_outfile', pdbID + '_' + chain.get_id(), '-adirectory_outfile', chainSequencesDir])
	print('OK')

def selectChains(pdbID, chainSequencesDir, verbose):
	print ('Automatically choosing chains to be adjusted...')
	chainIdentity = {}
	adjustChains = []
	files = os.listdir(chainSequencesDir)
	for file in files:
		if file.startswith(pdbID) and file.endswith('.needle'):
			chain = file[-8:-7]
			fileOutput = os.path.join(chainSequencesDir, file)
			with open(fileOutput, 'r') as f:
				for line in f.readlines():
					if line.startswith('# Identity:'):
						words = line.split()
						try:
							#for >=10
							chainIdentity[chain] = float(words[3][1:-2])
						except ValueError:
							#for <10
							chainIdentity[chain] = float(words[4][0:-2])
	if verbose:
		print('Chain identity values: ')
	for chain, value in chainIdentity.items():
		if verbose:
			print(chain + " - " + str(value) + "%")
		if (value == max(chainIdentity.values())):
			adjustChains.append(chain)
	if verbose:
		print('Chains to be adjusted: ' + ",".join(adjustChains))
	print('OK')
	return adjustChains

def getChainEntropyVariationData(pdbID, entropyFile, consensusVariationFile, chainSequencesDir, adjustChains, outputDir, verbose):
	print('Adjusting chains...')
	with open(entropyFile, 'r') as f:
		myEntrlines = f.readlines()
	with open(consensusVariationFile, 'r') as f:
		myConslines = f.readlines()
	entropyData = {}
	consevationData = {}
	for chain in adjustChains:
		if verbose:
			print('Adjusting chain: ' + chain)
		#chain = chain.upper()
		entropyData[chain] = []
		consevationData[chain] = []
		alignmentOutput = os.path.join(chainSequencesDir, pdbID + '_' + chain + '.needle')
		sequences = AlignIO.read(alignmentOutput, "emboss")
		dataOutput = os.path.join(outputDir, pdbID + '_' + chain + '_entropy_conservation.txt')
		with open(dataOutput, 'w') as f:
			writer = csv.writer(f, dialect='excel-tab')
			writer.writerow([sequences[0].id, sequences[1].id, 'H %', 'cons %'])
			sequenceResidues = sequences[0].seq
			consensusResidues = sequences[1].seq
			ci = 1;
			for i in range(sequences.get_alignment_length()):
				if sequenceResidues[i] != '-':
					if consensusResidues[i] == '-':
						entropy = 110.0
						conservation = 110.0
					else:
						line = myEntrlines[ci][:-1]
						line = myConslines[ci][:-1]
						myEntrList = line.split('\t')
						myConslist = line.split('\t')
						entropy = float(myEntrList[2])
						conservation = float(myConslist[1])
						ci+=1
					entropy = round(float(entropy), 2)
					conservation = round(float(conservation), 2)
					entropyData[chain].append([sequenceResidues[i], entropy])
					consevationData[chain].append([sequenceResidues[i], conservation])
					writer.writerow([sequenceResidues[i], consensusResidues[i], entropy, conservation])
				else:
					#TODO: what then ?
					if consensusResidues[i] != '-':
						ci+=1
	print('OK')
	return(entropyData, consevationData)

def getChainFrequencyTypeData(pdbID, frequencyTypeFile, chainSequencesDir, adjustChains, outputDir, verbose):
	print('Adjusting chains...')
	with open(frequencyTypeFile, 'r') as f:
		mySNPlines = f.readlines()
	frequencyData = {}
	typeData = {}
	typeNumber = {}
	typeList = ['-']
	for chain in adjustChains:
		if verbose:
			print('Adjusting chain: ' + chain)
		#chain = chain.upper()
		frequencyData[chain] = []
		typeData[chain] = []
		typeNumber[chain] = []
		alignmentOutput = os.path.join(chainSequencesDir, pdbID + '_' + chain + '.needle')
		sequences = AlignIO.read(alignmentOutput, "emboss")
		dataOutput = os.path.join(outputDir, pdbID + '_' + chain + '_frequency_type.txt')
		with open(dataOutput, 'w') as f:
			writer = csv.writer(f, dialect='excel-tab')
			writer.writerow([sequences[0].id, sequences[1].id, 'frequency', 'type'])
			sequenceResidues = sequences[0].seq
			consensusResidues = sequences[1].seq
			ci = 1;
			for i in range(sequences.get_alignment_length()):
				if sequenceResidues[i] != '-':
					typeScore = 0
					if consensusResidues[i] == '-':
						frequencyScore = 0
						types = '-'
					else:
						line = mySNPlines[ci][:-1]
						words = line.split('\t')
						frequencyScore = int(words[1])
						types = words[3]
						for t in types.split(','):
							p = 0
							try:
								p = typeList.index(t)
							except ValueError:
								typeList.append(t)
								p = typeList.index(t)
							typeScore += 2**p
						ci += 1
					frequencyData[chain].append([sequenceResidues[i], frequencyScore])
					typeData[chain].append([sequenceResidues[i], types])
					typeNumber[chain].append([sequenceResidues[i], typeScore])
					writer.writerow([sequenceResidues[i], consensusResidues[i], frequencyScore, types])
				else:
					#TODO: what then ?
					if consensusResidues[i] != '-':
						ci+=1
	print('OK')
	return(frequencyData, typeData, typeNumber, typeList)

def changeBfactors(structure, adjustChains, substitutionData, verbose, noDataValue):
	print('Changing b-factors of atoms in PDB...')
	for model in structure:
		for chain in model:
			chainID = chain.get_id()
			if chainID in adjustChains:
				residueID = 0
				for residue in chain:
					residueName = SeqUtils.seq1(residue.get_resname())
					if residueName != substitutionData[chainID][residueID][0]:
						if verbose:
							print("WARNING: Unexpected residue " + str(residueID) + " in chain " + str(chainID))
							print("Expected residue " + substitutionData[chainID][residueID][0] + " got " + residueName)
						continue
					(heteroFlag, sequenceID, insertionCode) = residue.get_id()
					if heteroFlag != ' ':
						continue
					for atom in residue:
						value = float(substitutionData[chainID][residueID][1])
						if atom.is_disordered():
							atom.__class__ = Atom.DisorderedAtom
							for altloc in atom.disordered_get_id_list():
								a = atom.disordered_get(altloc)
								if verbose:
									print("Chain: " + chainID + "\t residue: " + residueName + "\t atom: " + str(a.get_full_id()) + " \t b-factor: " + str(a.get_bfactor()) + " => " + str(value))
								a.set_bfactor(value)
						else:
							if verbose:
								print("Chain: " + chainID + "\t residue: " + residueName + "\t atom: " + str(atom.get_full_id()) + " \t b-factor: " + str(atom.get_bfactor()) + " => " + str(value))
							atom.set_bfactor(value)
					residueID += 1
					if (residueID >= len(substitutionData[chainID])):
						break
			else:
				for residue in chain:
					for atom in residue:
						atom.set_bfactor(float(noDataValue))
	print('OK')
	return structure

def extractPDBdata(structure, adjustChains, substitutionData, verbose):
	print('Extracting atoms details from PDB...')
	pdbData = {}
	for model in structure:
		for chain in model:
			chainID = chain.get_id()
			if chainID in adjustChains:
				pdbData[chainID] = {}
				residueID = 0
				for residue in chain:
					residueName = SeqUtils.seq1(residue.get_resname())
					if residueName != substitutionData[chainID][residueID][0]:
						continue
					(heteroFlag, sequenceID, insertionCode) = residue.get_id()
					if heteroFlag != ' ':
						continue
					value = substitutionData[chainID][residueID][1]
					if value != "-":
						pdbData[chainID][sequenceID] = value
					if verbose:
						print("Chain: " + chainID + "\t residue: " + residueName + " " + str(sequenceID) + "\t value: " + value)
					residueID += 1
					if (residueID >= len(substitutionData[chainID])):
						break
	print('OK')
	return pdbData

def getUniprotRecord(uniprotID, recordFile = None, sequenceFile = None):
	print("Downloading record for " + uniprotID + " from UniProt...")
	response = urllib2.urlopen('https://www.ebi.ac.uk/proteins/api/proteins/' + uniprotID)
	proteinJson = response.read()
	protein = json.loads(proteinJson)
	ensemblID = None
	organismName = ''
	if recordFile != None:
		with open(recordFile, 'w', 0) as f:
			proteinJson = json.dumps(protein, indent=4, sort_keys=True)
			f.write(proteinJson + "\n")
			f.flush()
	if sequenceFile != None:
		for ref in protein['dbReferences']:
			if ref['type'] == 'Ensembl':
				ensemblID = ref['id']
		if ensemblID == None:
			raise Exception('Protein record in Uniprot does not have reference to Ensembl')
		organism = protein['organism']['names']
		for o in organism:
			if o['type'] == 'scientific':
				organismName = o['value']
				break
		record = SeqRecord( Seq.Seq(protein['sequence']['sequence'], IUPAC.protein), id=protein['accession'], name=protein['id'], description=protein['id'])
		with open(sequenceFile, "w") as f:
			SeqIO.write([record], f, "fasta")
	print('OK')
	return ensemblID, organismName

def getUniprotVariants(uniprotID, variantsFile, verbose):
	print("Downloading record for " + uniprotID + " from UniProt...")
	response = urllib2.urlopen('https://www.ebi.ac.uk/proteins/api/variation/' + uniprotID)
	variantsJson = response.read()
	variants = json.loads(variantsJson)
	with open(variantsFile, 'w') as f:
		variantsJson = json.dumps(variants, indent=4, sort_keys=True)
		f.write(variantsJson + "\n");
	if verbose:
		print ("Protein " + variants['accession'] + " " + variants['entryName'] + " has " + str(len(variants['features'])) + " features.")
	features = {}
	xrefIndex = {}
	for feature in variants['features']:
		variant = {
			'id': "",
			'xrefs': {}
		}
		featureID = ""
		if 'ftId' in feature:
			variant['id'] = feature['ftId']
			featureID += feature['ftId']
		if 'genomicLocation' in feature:
			featureID += "(" + feature['genomicLocation'] + ")";
		if 'begin' not in feature:
			print("WARNING: feature " + featureID + " begin is not defined - skipping...")
			continue
		skip = False;
		if 'xrefs' in feature:
			for xref in feature['xrefs']:
				variant['xrefs'][xref['id']] = xref['id']
				xrefIndex[xref['id']] = xref['id']
		position = int(feature['begin'])
		variant['position'] = position
		variant['wildType'] = feature['wildType']
		variant['alternativeSequence'] = feature['alternativeSequence']
		variant['consequence']= feature['consequenceType']
		if (position not in features):
			features[position] = []
		features[position].append(variant)
		if verbose:
			print('Adding feature ' + featureID + " at position " + str(position) + " ids: " + repr(variant['xrefs']))
	print('OK')
	return (features, xrefIndex)

def getEnsemblVariants(organismName, ensemblID, variantsFile, verbose):
	print("Getting variants data from EnsemblID " + ensemblID)
	with open(variantsFile, 'w') as f:
		subprocess.check_call(['perl', os.path.dirname(os.path.abspath(__file__)) + '/../EnsemblWorker/worker.pl', organismName, ensemblID], stdout=f)
	variants = {}
	with open(variantsFile, 'r') as f:
		features = f.readlines()
	for feature in features:
		feature = feature[:-1].split("\t")
		if feature[1] == 'HGMD_MUTATION':
			continue
		try :
			position = int(feature[4].split(' ')[0])
		except ValueError:
			try :
				position = int(feature[4].split(' ')[1])
			except ValueError:
				continue
		print('Adding feature at position '+ str(position) + ' : ' + repr(feature))
		sequence = feature[5].split('/')
		if len(sequence) == 2:
			alternativeSequence = sequence[1]
		else:
			alternativeSequence = sequence[0]
		variant = {
			'id': feature[0],
			'dbSNP': feature[0],
			'position' : position,
			'wildType' : sequence[0],
			'alternativeSequence' : alternativeSequence,
			'consequence' : feature[3]
		}
		if (position not in variants):
			variants[position] = []
		variants[position].append(variant)
	print('OK')
	return variants

def getSnpSites(alignmentFile, fastaFile, vcfFile, verbose):
	dnaRecs = [f for f in SeqIO.parse(alignmentFile, 'fasta')]
	with open(fastaFile, "w") as f:
		SeqIO.write(dnaRecs, f, "fasta")
	params = ['snp-sites', '-rv', '-o', vcfFile, fastaFile]
	subprocess.check_call(params)
	

def callVariants(dnaSeqenceFile, protAlignmentFile, protConsSeqFile, outputFile, verbose):
	# load the proteins alignment
	protAlnRecs = [f for f in SeqIO.parse(protAlignmentFile, 'fasta')]
	strain_protAlnSeq = {rec.id:''.join([nt for nt in rec.seq]) for rec in protAlnRecs}
	
	# load raw dna sequences
	dnaRecs = [f for f in SeqIO.parse(dnaSeqenceFile, 'fasta')]
	strain_dnaSeq = {rec.id:''.join([nt for nt in rec.seq]) for rec in dnaRecs}

	# load protein consensus seq
	consensusProtSeq = SeqIO.read(protConsSeqFile, 'fasta')

	strainProtIdx = {rec.id:0 for rec in dnaRecs}
	strains = [rec.id for rec in dnaRecs]
	snps = {rec.id:{} for rec in dnaRecs}
	
	seq_len = 0
	for s in strain_protAlnSeq:
		seq_len = len(strain_protAlnSeq[s])
	
	with open(outputFile, 'w') as f:
		writer = csv.writer(f, dialect='excel-tab')
		writer.writerow(['residue','freq','variants', 'type', 'ids'])
		
		# superimpose  protein alignment with nucleotide sequences and calculate snps
		for i in range(0, seq_len):
			cdns = {}
			freq = 0
			variants = []
			types = []
			ids = []
			counts = [{}, {}, {}]
			consAA = consensusProtSeq.seq[i]
			for strain in strains:
				aa = strain_protAlnSeq[strain][i]
				if (aa != '-'):
					idx = strainProtIdx[strain]
					cdn = strain_dnaSeq[strain][idx*3 : (idx+1)*3]
					counts[0][cdn[0]] = counts[0].get(cdn[0], 0) + 1
					counts[1][cdn[1]] = counts[1].get(cdn[1], 0) + 1
					counts[2][cdn[2]] = counts[2].get(cdn[2], 0) + 1
					cdns[strain] = cdn
			consCdn = "".join([max(counts[0].iteritems(), key=operator.itemgetter(1))[0], max(counts[1].iteritems(), key=operator.itemgetter(1))[0], max(counts[2].iteritems(), key=operator.itemgetter(1))[0]])
			consCdnAA = Seq.Seq(consCdn, IUPAC.unambiguous_dna).translate()
			for strain in strains:
				aa = strain_protAlnSeq[strain][i].upper()
				if (aa != '-'):
					idx = strainProtIdx[strain]
					cdn = strain_dnaSeq[strain][idx*3 : (idx+1)*3]
					if aa != consAA.upper() or cdn[0] != consCdn[0] or cdn[1] != consCdn[1] or cdn[2] != consCdn[2]:
						freq+=1
						variants.append(aa)
						if aa == consAA.upper():
							types.append('synonymous_variant')
						else:
							types.append('missense_variant')
						ids.append(strain)
						#if verbose:
						#	print str(i) + " " + consAA + "->" + aa  + " : " + strain + " " + str(idx) + " " + consCdn + " -> " + cdn + "(" + consCdnAA + "->" + aa + ")"
					strainProtIdx[strain]+=1
			if freq > 0:
				variants = ",".join(variants)
				types = ",".join(types)
				ids = ",".join(ids)
			else:
				variants = "-"
				types = "-"
				ids = "-"
			if verbose:
				print("Position: " + str(i) + " residue: " + consAA + " frequency: " + str(freq) + " variants: " + variants + " type: " + types)
			writer.writerow([consAA, freq, variants, types, ids])

def calcFrequencyType(features, sequenceFile, outputFile, verbose):
	print("Calculating SNP frequency...")
	sequence = SeqIO.read(sequenceFile, "fasta")
	with open(outputFile, 'w') as f:
		writer = csv.writer(f, dialect='excel-tab')
		writer.writerow(['residue','freq','variants', 'type', 'ids'])
		for residueID, residue in enumerate(sequence.seq):
			position = residueID + 1;
			freq = 0
			change = []
			snptype = []
			ids = []
			if position in features:
				freq = len(features[position])
				for f in features[position]:
					#wildtype = "-";
					#if len(f['wildType']) > 0:
					#	wildtype = f['wildType'][0]
					#if residue != wildtype:
					#	print("WARNING: sequence does not match: residue " + residue + "<>" + wildtype + " feature " + repr(f))
					change.append(f['alternativeSequence'])
					snptype.append(f['consequence'])
					fid = f['id']
					if 'xrefs' in f:
						fid += '(';
						fid += ",".join(str(v) for v in f['xrefs'].values());
						fid += ')';
					ids.append(fid)
				if len(change) != freq:
					print "WARNING: frequency different than number of changes"
				if len(snptype) != freq:
					print "WARNING: frequency different than number of types"
				if len(ids) != freq:
					print "WARNING: frequency different than number of ids"
			else:
				ids.append('-')
				change.append('-')
				snptype.append('-')
			change = ",".join(change)
			snptype = ",".join(snptype)
			ids = ",".join(ids)
			if verbose:
				print("Position: " + str(position) + " residue: " + residue + " frequency: " + str(freq) + " variants: " + change + " type: " + snptype)
			writer.writerow([residue, freq, change, snptype, ids])
	print('OK')

def savePDB(structure, outputFile, description = "", header = []):
	print('Saving ' + description + ' PDB...')
	io = PDBIO()
	io.set_structure(structure)
	io.save(outputFile)
	with open(outputFile, 'r+') as f:
		content = f.read()
		f.seek(0, 0)
		f.write("".join(header).rstrip('\r\n') + '\n' + content)
	print('OK')
