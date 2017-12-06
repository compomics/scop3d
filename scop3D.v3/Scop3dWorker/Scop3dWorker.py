
""" Scop3D Worker - worker module for Scop3d """

import csv
import math
import os
import json
import subprocess
import tempfile
import urllib2

from Bio import AlignIO
from Bio import ExPASy
from Bio import SeqIO
from Bio import SeqUtils
from Bio import SwissProt

from Bio.Align.AlignInfo import SummaryInfo
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.PDB import PDBList
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from weblogolib import *

### Stage 1

myNoDataValue = 110

# 0. Creating workdirs
def mkdirs(outputDir, verbose):
	try:
		os.mkdir(outputDir)
	except OSError:
		pass
	sequencesDir = tempfile.mkdtemp()
	alignmentsDir = tempfile.mkdtemp()
	chainSequencesDir = tempfile.mkdtemp()
	if verbose:
		print("Using following directories:")
		print("output: " + outputDir)
		print("sequences: " + sequencesDir)
		print("alignments: " + alignmentsDir)
		print("chains: " + chainSequencesDir)
	return (sequencesDir, alignmentsDir, chainSequencesDir)

# 1. Retrieval of data input
def downloadUniprotSequences(uniprotID, blastFile, sequencesFile, cutoff, verbose):
	print('Obtaining sequences from UniProt...')
	with open(blastFile, 'r') as f:
		records = NCBIXML.read(f)
	if verbose:
		print('Found ' + str(len(records.alignments)) + ' matches')
	with open(sequencesFile, 'w') as f:
		sequence = urllib2.urlopen('http://www.uniprot.org/uniprot/' + uniprotID + '.fasta')
		f.write(sequence.read() + "\n");
		for idx, alignment in enumerate(records.alignments):
			for hsp in alignment.hsps:
				title = alignment.title
				words = title.split('|')
				seqID = words[3]
				identityPercent = 100.0 * float(hsp.identities) / float(alignment.length)
				if (identityPercent >= float(cutoff)):
					try:
						sequence = urllib2.urlopen('http://www.uniprot.org/uniprot/' + seqID + '.fasta')
						f.write(sequence.read() + "\n");
						if verbose:
							print(seqID + " (identity " + str(identityPercent) + "% >= cutoff " + str(cutoff) + "%) - adding")
					except Exception as e:
						print("WARNING: unable to download this entry: " + str(e))
				else:
					if verbose:
						print(seqID + " (identity " + str(identityPercent) + "% < cutoff " + str(cutoff) + "%) - skipping")
	print('OK')

def loadSequences(sequencesFile, verbose):
	print('Loading sequences from ' + sequencesFile + '...')
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
	params = ['muscle', '-in', sequencesFile, '-out', alignmentFile, '-clwstrict']
	if verbose:
		print('Writing alignment to: ' + alignmentFile)
	else:
		params.append('-quiet')
	subprocess.check_call(params)
	print('OK')

# 3. Determination of similarity and identity between the different sequences
# Reference: Rice P, Longden I, Bleasby A (2000) EMBOSS: The European Molecular Biology Open Software Suite.
# 			 Trends in Genetics 16(6), 276-277
def doSplit(sequencesFile, sequencesDir, verbose):
	print('Splitting fasta file into individual sequences...')
	if verbose:
		print('Output directory: ' + sequencesDir)
	subprocess.check_call(['seqretsplit', '-feature', '0', '-sequence', sequencesFile, '-firstonly', '0', '-auto', '-osdirectory2', sequencesDir])
	print('OK');

# 4. A pairwise alignment among all sequences is performed with the aid of needle
# Reference: Needleman SB, Wunsch CD (1970) A general method applicable to the search for similarities in the amino acid sequence of two proteins
#			 J Mol Biol 48(3), 443-453
def doPairwiseAlignment(sequencesDir, alignmentsDir, verbose):
	print('Performing a pairwise alignment...')
	if verbose:
		print('Output directory: ' + alignmentsDir)
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
	if verbose:
		print('Writing statistics to: ' + alignmentStatsFile)
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
	if verbose:
		print('Writing sequence to: ' + consensusFile)
	subprocess.check_call(['cons', '-sequence', alignmentFile, '-plurality', '1', '-identity', '0', '-auto', '-outseq', consensusFile, '-name', 'consensus'])
	if verbose:
		consensusSeq = SeqIO.read(consensusFile, 'fasta')
		print('Consensus sequence: ')
		print(consensusSeq.seq)
	print('OK')

# 7. Calculation of the abundance matrices - biopython module
def calcAbundance(alignmentFile, consensusFile, abundanceFile, abundancePercentFile, verbose):
	print('Calculating the abundance matrix...')
	alignment = AlignIO.read(alignmentFile, "clustal")
	summary = SummaryInfo(alignment)
	consensusSeq = SeqIO.read(consensusFile, 'fasta')
	abundanceMatrix = summary.pos_specific_score_matrix(consensusSeq)
	if verbose:
		print("Abundance matrix (absolute values):")
		print(str(abundanceMatrix))
		print("Saved as: " + abundanceFile)
	with open(abundanceFile, 'w') as f:
		f.write(str(abundanceMatrix))
	for pos, abundance in enumerate(abundanceMatrix):
		for res, value in abundance.items():
			abundanceMatrix[pos][res] = 100.0 * float(value) / float(len(alignment))
	if verbose:
		print("Abundance matrix (percentages):")
		print(str(abundanceMatrix))
		print("Saved as: " + abundancePercentFile)
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
			maximum = float(max(words[1:]))
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
	if verbose:
		print('Logo saved as: ' + logoFile)
	print('OK')

### Stage 2


def doUniprotBlast(uniprotid, blastFile, verbose):
	print('Performing BLAST with the UniProt ID ' + uniprotid + '...')
	results = NCBIWWW.qblast('blastp', 'swissprot', uniprotid)
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
						chainIdentity[chain] = words[3][1:-2]
	if verbose:
		print('Chain identity values: ')
	for chain, value in chainIdentity.items():
		if verbose:
			print(chain + " - " + value + "%")
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
		chain = chain.upper()
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
	for chain in adjustChains:
		if verbose:
			print('Adjusting chain: ' + chain)
		chain = chain.upper()
		frequencyData[chain] = []
		typeData[chain] = []
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
					if consensusResidues[i] == '-':
						frequencyScore = 0
						typeScore = 0
					else:
						line = mySNPlines[ci][:-1]
						words = line.split('\t')
						frequencyScore = int(words[1])
						types = words[3].split(',')
						typeScoreTable = {}
						typeScore = 0
						for t in types:
							typeScoreTable[t] = 1
						if 'missense' in typeScoreTable:
							typeScore += 1
						if 'stop_lost' in typeScoreTable:
							typeScore += 2
						if 'stop_gained' in typeScoreTable:
							typeScore += 4
						ci+=1
					frequencyData[chain].append([sequenceResidues[i], frequencyScore])
					typeData[chain].append([sequenceResidues[i], typeScore])
					writer.writerow([sequenceResidues[i], consensusResidues[i], frequencyScore, typeScore])
				else:
					#TODO: what then ?
					if consensusResidues[i] != '-':
						ci+=1
	print('OK')
	return(frequencyData, typeData)

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
						continue
					(heteroFlag, sequenceID, insertionCode) = residue.get_id()
					if heteroFlag != ' ':
						continue
					for atom in residue:
						value = float(substitutionData[chainID][residueID][1])
						if verbose:
							print("Chain: " + chainID + "\t residue: " + residueName + "\t atom: " + atom.get_id() + "\t b-factor: " + str(atom.get_bfactor()) + " => " + str(value))
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

def getUniprotVariants(uniprotID, variantsFile, sequenceFile, verbose):
	print("Downloading record for " + uniprotID + " from UniProt...")
	response = urllib2.urlopen('https://www.ebi.ac.uk/proteins/api/variation/' + uniprotID)
	variantsJson = response.read()
	variants = json.loads(variantsJson)
	with open(variantsFile, 'w') as f:
		variantsJson = json.dumps(variants, indent=4, sort_keys=True)
		f.write(variantsJson + "\n");
	record = SeqRecord( Seq(variants['sequence'], IUPAC.protein), id=variants['accession'], name=variants['entryName'], description=variants['entryName'])
	with open(sequenceFile, "w") as f:
  		SeqIO.write([record], f, "fasta")
	if verbose:
		print ("Protein " + variants['accession'] + " " + variants['entryName'] + " has " + str(len(variants['features'])) + " features.")
	print('OK')
	return variants

def calcFrequencyType(variants, outputFile, verbose):
	print("Calculating SNP frequency...")
	if verbose:
		print('Reading variants...')
	features = {};
	for feature in variants['features']:
		featureID = "";
		if 'ftId' in feature:
			featureID += feature['ftId'];
		if 'genomicLocation' in feature:
			featureID += "(" + feature['genomicLocation'] + ")";
		if ('begin' not in feature) or ('end' not in feature):
			print("WARNING: feature " + featureID + " begin and/or end are not defined - skipping...")
			continue
		if (feature['begin'] != feature['end']):
			print("WARNING: feature " + featureID + " begin and end are different - skipping...")
			continue
		position = int(feature['begin'])
		if (position not in features):
			features[position] = []
		features[position].append(feature)
		if verbose:
			print('Adding feature ' + featureID + " at position " + str(position))
	if verbose:
		print('Analyzing sequence...')
	with open(outputFile, 'w') as f:
		writer = csv.writer(f, dialect='excel-tab')
		writer.writerow(['residue','freq','variants', 'type'])
		for residueID, residue in enumerate(variants['sequence']):
			freq = 0
			change = []
			snptype = []
			if residueID in features:
				freq = len(features[residueID])
				for f in features[residueID]:
					change.append(f['alternativeSequence'])
					snptype.append(f['consequenceType'])
			else:
				change.append('-')
				snptype.append('-')
			change = ",".join(change)
			snptype = ",".join(snptype)
			if verbose:
				print("Position: " + str(residueID + 1) + " residue: " + residue + " frequency: " + str(freq) + " variants: " + change + " type: " + snptype)
			writer.writerow([residue, freq, change, snptype])
	print('OK')

def savePDB(structure, outputFile, description = ""):
	print('Saving ' + description + ' PDB...')
	io = PDBIO()
	io.set_structure(structure)
	io.save(outputFile)
	print('OK')
