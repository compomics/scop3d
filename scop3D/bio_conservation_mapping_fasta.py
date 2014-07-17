from Bio.Align.Applications import MuscleCommandline
import os
import math
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SeqIO
from bio_retrieve_fasta import spacer
from bio_retrieve_fasta import fasta_dir
from bio_retrieve_fasta import project_dir
from bio_retrieve_fasta import project_name
from bio_retrieve_fasta import cons_tresh
from bio_retrieve_fasta import muscle_loc

file_ali_in = fasta_dir
file_ali_out = project_dir + spacer + "alignments" + spacer + "alignment_" + project_name + ".clw"
frequency_all = []
consensus = ""
procent_complete = []
entropy_position_list=[]

def make_dir():
    # make two directories: "project/alignments", "project/matrices/"
    os.makedirs(project_dir + spacer + "alignments")
    os.makedirs(project_dir + spacer + "matrices")


def lenght_seq():
    # check the length of all sequences in fasta-file (same length is preferred).
    fastafile = fasta_dir
    fasta = open(fastafile)
    length_seq = []
    length_unique_seq = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        length_seq.append(len(seq_record))
    for i in length_seq:
        if i not in length_unique_seq:
            length_unique_seq.append(i)
    if len(length_unique_seq) != 1:
        info_file= open(project_dir + spacer +"info.txt", "a")
        info_file.write("Warning: \n \t Not all seqeunces are of the same length."
                        "\n \t differences in length: " + str(length_unique_seq))
        info_file.close()


def alignment():
    # align sequences with muscle, (http://www.drive5.com/muscle/)
    if muscle_loc:
        muscle_cline = MuscleCommandline(muscle_loc, input=file_ali_in, out=file_ali_out, clwstrict=True)
        muscle_cline()
    else:
        muscle_cline = MuscleCommandline(input=file_ali_in, out=file_ali_out, clwstrict=True)
        muscle_cline()
    alignment = open(file_ali_out, "r")
    print_alignment = alignment.read()
    print print_alignment
    alignment.close()


def frequency():
    global consensus
    global frequency_matrix
    global aadict
    #read fasta file and detmine number of sequences.
    fasta_file = AlignIO.read(file_ali_out, "clustal")
    summary_align = AlignInfo.SummaryInfo(fasta_file)
    amount_seq = len(fasta_file)

    for record in fasta_file:
        print record.seq, record.id

    # make simple consensus sequence with treshold value. gap is marked with "-"
    consensus = summary_align.dumb_consensus(cons_tresh, "-")
    print "\n", consensus
    info_file= open(project_dir + spacer +"info.txt", "a")
    info_file.write("\n consensus sequence: \n" + str(consensus))
    info_file.close()

    # establish the abundance of a certain residue at each position in consensus sequence.
    frequency_matrix = summary_align.pos_specific_score_matrix(consensus)
    frequency_matrix_str = str(frequency_matrix)
    print frequency_matrix_str

    #print tis matrix to a file
    print ("writing file of frequencies")
    frequencies_file = open(project_dir + spacer + "matrices" + spacer + "frequencies.txt", "w")
    frequencies_file.writelines(frequency_matrix_str)
    frequencies_file.close()

    #read first line of the file (determines the AA present in the alignment
    frequency_file = open(project_dir + spacer + "matrices" + spacer + "frequencies.txt", "r")
    aa = frequency_file.readline()

    # making list of the AA, deleting blank space.
    aadict = []
    for i in range (len(aa)):
        if aa[i] != " ":
            aadict.append(aa[i])
    frequency_file.close()

    #delete last item in list because it is not an AA (\n)
    del aadict[-1]
    info_file= open(project_dir + spacer +"info.txt", "a")
    info_file.write("\n Amino Acids found in sequences: \n \t" + str(aadict))
    info_file.close()
    print aadict

    #make quotient of all input by amount of entries
    lines = range(0, len(consensus), 1)
    number_aa = int(len(aadict))

    #make quotient of all input by amount of entries
    n = 0
    #make first line of csv file
    procent_line = []
    global procent_complete
    global entropy_position_list
    procent_file = open(project_dir + spacer + "matrices" + spacer + "procent.csv", "a")
    procent_file.write(str(aadict) + "\n")
    procent_file.close()
    entropy_file = open(project_dir + spacer + "matrices" + spacer + "entropy.csv", "a")
    entropy_file.write(str(aadict) + "\n")
    entropy_file.close()
    for number in lines:
        entropy_list=[]
        entropy_num = 0
        entropy_position = 0
        for AA in aadict:
            entropy_num += frequency_matrix[number][AA]

        for AA in aadict:
            #abundance
            procent = 1 - (frequency_matrix[number][AA] / amount_seq)
            procent_line.append(procent)

            #entropy
            # calculate p.
            quotient_entropy = frequency_matrix[number][AA]/entropy_num
            # due to error with log2(0) if is needed.
            if quotient_entropy != 0:
                log_entropy = math.log(quotient_entropy,2)
                entropy = quotient_entropy * log_entropy
                entropy = math.fabs(entropy)
                entropy_list.append(entropy)
            else:
                entropy_list.append(0.000)
        # set numbers to 3digits behind the comma.
        for entropy_item in entropy_list:
            entropy_position = float(entropy_position)
            entropy_position += float(entropy_item)
            entropy_position = float(entropy_position)
            entropy_position = "%.3f"%entropy_position
        # make list of all entropy values.
        entropy_position_list.append(entropy_position)

        #write line to csv file
        print(procent_line)
        procent_line_str = str(procent_line)
        procent_file = open(project_dir + spacer + "matrices" + spacer + "procent.csv", "a")
        procent_file.write(procent_line_str)
        procent_file.write("\n")
        procent_file.close()
        procent_complete.append(procent_line)
        procent_line = []
        #starting new line
        n = 0
    entropy_file = open(project_dir + spacer + "matrices" + spacer + "entropy.csv", "a")
    entropy_file.write(str(entropy_position_list))
    entropy_file.close()
    print entropy_position_list
    print procent_complete


def conservation_mapping():
    make_dir()
    lenght_seq()
    alignment()
    frequency()

conservation_mapping()




