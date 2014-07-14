from Bio.PDB import *
import wx
import wx.xrc
import os
import shutil
from Bio.PDB import *
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from bio_conservation_mapping_fasta import procent_complete
from bio_conservation_mapping_fasta import spacer
from bio_conservation_mapping_fasta import project_dir
from bio_conservation_mapping_fasta import consensus
from bio_conservation_mapping_fasta import entropy_position_list
from bio_blast_pdb_fasta import project_name
from bio_blast_pdb_fasta import pdbfile
from bio_retrieve_fasta import output2
from bio_retrieve_fasta import output3
from bio_retrieve_fasta import output4
from weblogolib import *
from bio_retrieve_fasta import muscle_loc


bigest_amount_complete = []
chains = ""
list_positions = []
print project_name
list_unique_chains = []
conversion_fail = ""

class MyFrame1 ( wx.Frame ):

	def __init__( self, parent ):
        # choose the chains that must be adjusted, all the other chains are filled with arbitrary values,
        # for better coloring.
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"Scop3D - chains", pos = wx.DefaultPosition, size = wx.Size( 307,250 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )

		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

		lbl_chain = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"choose the chains to change:" ), wx.VERTICAL )

		self.m_panel1 = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		gSizer1 = wx.GridSizer( 2, 1, 0, 0 )

		chainsChoices = list_unique_chains
		self.chains = wx.CheckListBox( self.m_panel1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, chainsChoices, 0 )
		gSizer1.Add( self.chains, 0, wx.ALL, 5 )

		self.button_select = wx.Button( self.m_panel1, wx.ID_ANY, u"Select", wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer1.Add( self.button_select, 0, wx.ALL|wx.ALIGN_RIGHT|wx.ALIGN_BOTTOM, 5 )


		self.m_panel1.SetSizer( gSizer1 )
		self.m_panel1.Layout()
		gSizer1.Fit( self.m_panel1 )
		lbl_chain.Add( self.m_panel1, 1, wx.EXPAND |wx.ALL, 5 )


		self.SetSizer( lbl_chain )
		self.Layout()

		self.Centre( wx.BOTH )

		# Connect Events
		self.button_select.Bind( wx.EVT_BUTTON, self.done )

	def __del__( self ):
		pass


	# Virtual event handlers, overide them in your derived class
	def done( self, event ):
		global chains
		chains = self.chains.GetCheckedStrings()
		print chains
		self.Close()


def bigest():
    # get the highest abundant values.
    for list in procent_complete:
        smallest = 1
        amount_aa = 0
        for item in list:
            if smallest > item:
                item = item
                smallest = item
        smallest = float(smallest)
        smallest = "%.3f"%smallest
        bigest_amount_complete.append([smallest])
    print (bigest_amount_complete)
    print len(bigest_amount_complete)
    print len(entropy_position_list)
    for x in range (0, len(entropy_position_list), 1):
        bigest_amount_complete[x].append(entropy_position_list[x])
    return bigest_amount_complete


def get_reference():
    # get a list of all the chains in the fasta-file
    list_chains = []
    global list_unique_chains
    pdb = open(str(pdbfile), "r")
    for lines in pdb:
        if lines[0:4] == "ATOM":
            if lines[13:15] == "CA":
                list_chains.append(lines[21])
    for i in list_chains:
        if i not in list_unique_chains:
            list_unique_chains.append(i)

    app = wx.App(0)
    MainApp = MyFrame1(None)
    MainApp.Show()
    app.MainLoop()


def get_alignment(chain):
    parser=PDBParser()
    structure=parser.get_structure("bla", pdbfile)
    ppb=CaPPBuilder()
    model = structure[0]
    chaint = model[chain]
    seq = ""
    for pp in ppb.build_peptides(chaint):
        seq_temp = pp.get_sequence()
        seq += seq_temp
    print str(seq)
    pairwise_aln = open(project_dir + spacer + "fasta_cons_" + str(chain) + ".fasta", "w")
    pairwise_aln.write(">consensus \n" + str(consensus) + "\n" + ">sequence " + str(chain) + "\n" + str(seq))
    pairwise_aln.close()
    if muscle_loc:
        muscle_cline = MuscleCommandline(muscle_loc, input=project_dir + spacer + "fasta_cons_" + str(chain) + ".fasta", out=project_dir + spacer + "pairwise_aln_cons_" + str(chain) + ".clw", clwstrict=True)
        muscle_cline()
    else:
        muscle_cline = MuscleCommandline(input=project_dir + spacer + "fasta_cons_" + str(chain) + ".fasta", out=project_dir + spacer + "pairwise_aln_cons_" + str(chain) + ".clw", clwstrict=True)
        muscle_cline()
    seq_alignment = []
    handle = open(project_dir + spacer + "pairwise_aln_cons_" + str(chain) + ".clw", "rU")
    for record in SeqIO.parse(handle, "clustal") :
        seq = str(record.seq)
        seq_alignment.append(seq)
    handle.close()
    print len(seq_alignment[0])
    print len(seq_alignment[1])
    print seq_alignment
    return seq_alignment


def get_positions(chain):
    bigest_amount_complete = bigest()
    alignment = get_alignment(chain)
    global list_positions
    list_positions = []
    pdb = open(str(pdbfile), "r")
    for lines in pdb:
        if lines[0:4] == "ATOM" and lines[21] == chain:
            if lines[13:15] == "CA":
                list_positions.append(lines[23:26])
    print list_positions
    global list_unique_positions
    list_unique_positions = []
    for i in list_positions:
        if i not in list_unique_positions:
            list_unique_positions.append(i)
    print list_unique_positions
    print len(bigest_amount_complete)
    bigest_alignment = []
    y = 0
    z = 0
    for i in range(0, len(consensus), 1):
        if consensus[i] == alignment[0][y]:
            bigest_alignment.append(bigest_amount_complete[z])
            y += 1
            z += 1
        else:
            if consensus[i] == "-":
                z += 1
            else:
                while consensus[i] != alignment[0][y]:
                    bigest_alignment.append([0.5, 1.0, "T"])
                    y += 1
                else:
                    bigest_alignment.append(bigest_amount_complete[z])
                    y += 1
                    z += 1
    if len(alignment[0]) != len(bigest_alignment):
        end_gaps = len(alignment[0])-len(bigest_alignment)
        print end_gaps
        for i in range (0,end_gaps,1):
            bigest_alignment.append([0.5, 1.0, "T"])
    print bigest_alignment
    position_align = []
    for position, item in enumerate(alignment[1]):
        if item == "-":
            position_align.append(position)
    for i in reversed(position_align):
        del (bigest_alignment[i])
    print bigest_alignment
    print len(bigest_alignment)
    print len(list_unique_positions)
    return bigest_alignment


def make_outputpdb(name, outputtype):
    shutil.copyfile(str(pdbfile), str(pdbfile + name))
    for chain in chains:
        bigest_amount_complete = get_positions(chain)
        n = 0
        x = 0
        pdb = open(str(pdbfile + name), "r")
        myOutfile = open(project_dir + spacer + "myOutfile_" + name + ".pdb", "w")
        for myPDBline in pdb:
            if len(myPDBline) >= 21:
                if myPDBline[0:4] == ("ATOM") and myPDBline[21] == chain:
                    position = myPDBline[23:26]
                    if list_unique_positions[x] == position:
                        myBfactor = myPDBline[60:66].strip()
                        print(myBfactor)
                        print(myPDBline[:60]+"  "+str(bigest_amount_complete[n][outputtype])+myPDBline[66:])
                        myOutfile.write(myPDBline[:60]+"  "+str(bigest_amount_complete[n][outputtype])+myPDBline[66:])
                    else:
                        n += 1
                        print(myPDBline[:60]+"  "+str(bigest_amount_complete[n][outputtype])+myPDBline[66:])
                        myOutfile.write(myPDBline[:60]+"  "+str(bigest_amount_complete[n][outputtype])+myPDBline[66:])
                        x += 1
                else:
                    myOutfile.write(myPDBline)
            else:
                myOutfile.write(myPDBline)
        myOutfile.close()
        pdb.close()
        os.remove(str(pdbfile + name))
        shutil.copyfile(project_dir + spacer + "myOutfile_" + name + ".pdb", str(pdbfile + name))
        os.remove(project_dir + spacer + "myOutfile_" + name + ".pdb")


def pdb_cleanup(name, output):
    pdb = open(str(pdbfile + name), "r")
    myOutfile = open(str(project_dir + spacer + "myOutfile_"+ name + ".pdb"), "w")
    for myPDBline in pdb:
        if myPDBline.startswith("ATOM") and myPDBline[21] not in chains:
            if output == "output3":
                myOutfile.write(myPDBline[:60]+ '  0.500' + myPDBline[66:])
            if output == "output4":
                myOutfile.write(myPDBline[:60]+ '  1.000' + myPDBline[66:])
        else:
            myOutfile.write(myPDBline)
    myOutfile.close()
    pdb.close()

def sequence_logo():
    fin = open(project_dir + spacer + "alignments" + spacer + "alignment_" + project_name + ".clw")
    seqs = read_seq_data(fin)
    data = LogoData.from_seqs(seqs)
    options = LogoOptions()
    options.title = "sequence logo"
    format = LogoFormat(data, options)
    eps=eps_formatter( data, format)
    fout = open(project_dir + spacer + "sequence_logo.eps", 'w')
    fout.write(eps)

def convert():
    input = project_dir + spacer + "sequence_logo.eps "
    output = project_dir + spacer + "sequence_logo.png "
    cmd = 'convert ' + input + output
    os.system(cmd)


def parse_pdb():
    get_reference()
    if output3 == 1:
        outputtype = 0
        name = "abundance"
        make_outputpdb(name, outputtype)
        output = "output3"
        pdb_cleanup(name, output)
    if output4 == 1:
        outputtype = 1
        name = "amount"
        make_outputpdb(name, outputtype)
        output = "output4"
        pdb_cleanup(name, output)
    if output2 == 1:
        sequence_logo()
        try:
            convert()
        except:
            global conversion_fail
            conversion_fail = "fail"

        # im = Image.open(project_dir + spacer + "sequence_logo.eps")
        # im.save(project_dir + spacer + "sequence_logo.png", quality=100, optimize=True)


parse_pdb()