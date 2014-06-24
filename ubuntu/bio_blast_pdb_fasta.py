import wx
import wx.xrc
import wx.richtext
from bio_conservation_mapping_fasta import consensus
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.PDB import PDBList
import shutil
from bio_conservation_mapping_fasta import spacer
from bio_conservation_mapping_fasta import project_dir
from bio_conservation_mapping_fasta import project_name
from bio_retrieve_fasta import define_pdb
from bio_retrieve_fasta import pdbfile
pdbfile = pdbfile
row = ""
choose_pdbChoices = []

class MyFrame1 ( wx.Frame ):

	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"Scop3D - BLAST", pos = wx.DefaultPosition, size = wx.Size( 567,328 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )

		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

		bSizer1 = wx.BoxSizer( wx.VERTICAL )

		gSizer3 = wx.GridSizer( 1, 1, 0, 0 )

		self.m_richText2 = wx.richtext.RichTextCtrl( self, wx.ID_ANY, "Press blast to blast this sequence \n\n" + str(consensus), wx.DefaultPosition, wx.DefaultSize, 0|wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER|wx.WANTS_CHARS )
		gSizer3.Add( self.m_richText2, 1, wx.EXPAND |wx.ALL, 5 )


		bSizer1.Add( gSizer3, 1, wx.EXPAND, 5 )

		self.blast = wx.ToggleButton( self, wx.ID_ANY, u"blast", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer1.Add( self.blast, 0, wx.ALL|wx.ALIGN_RIGHT, 5 )


		self.SetSizer( bSizer1 )
		self.Layout()

		self.Centre( wx.BOTH )

		# Connect Events
		self.blast.Bind( wx.EVT_TOGGLEBUTTON, self.done )

	def __del__( self ):
		pass


	# Virtual event handlers, overide them in your derived class
	def done( self, event ):
            result_handle = NCBIWWW.qblast('blastp', 'pdb', consensus)
            print "bla"
            save_file = open(project_dir + spacer + "my_blast.xml", "w")
            save_file.write(result_handle.read())
            save_file.close()
            result_handle.close()
            self.Close()


def blast():
    app = wx.App(0)
    MainApp =MyFrame1(None)
    MainApp.Show()
    app.MainLoop()


def download_pdb():
    result_handle = open(project_dir + spacer + "my_blast.xml")
    blast_record = NCBIXML.read(result_handle)
    global choose_pdbChoices
    i = 0
    list_blast_results = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
                title = alignment.title
                e_value = hsp.expect
                print e_value
                length = alignment.length
                identities = hsp.identities
                positives = hsp.positives
                gaps = hsp.gaps
                choice = ("title: " + str(title[0:100]) + " ..." +
                          "\n" + "score: " + str(e_value) + "  " +
                          "length: " + str(length) + "\n" +
                          "identities: " + str(identities) + "  " +
                          "positives: " + str(positives) + "  " +
                          "gaps: " + str(gaps) +  "\n" )
                if i <5:
                    choose_pdbChoices.append(choice)
                    list_blast_results.append(title)
                choice = ""
                i += 1
    print choose_pdbChoices

    app = wx.App(0)
    MainApp =MyFrame2(None)
    MainApp.Show()
    MainApp.Maximize()
    app.MainLoop()

    print "bla"
    print row
    title = list_blast_results[int(row)]
    info_file= open(project_dir + spacer +"info.txt", "a")
    info_file.write("\n you chose this blastresult: \n" + str(list_blast_results[int(row)]))
    info_file.close()

    print title[17:21]
    global pdbfile
    pdb_dir = project_dir + spacer + project_name +"_blast_pdb"
    pdbfile = project_dir + spacer + "tmppdb.pdb"
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(title[17:21], obsolete = False, pdir = pdb_dir)
    title_lower = title[17:21].lower()
    shutil.copyfile(pdb_dir + spacer + "pdb" + title_lower + ".ent", project_dir + spacer + "tmppdb.pdb")


class MyFrame2 ( wx.Frame ):
    def __init__( self, parent ):
        wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"Scop3D - BLAST results", pos = wx.DefaultPosition, size = wx.Size( 603,413 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )

        self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

        Results = wx.StaticBoxSizer( wx.StaticBox( self, wx.ID_ANY, u"Blast results" ), wx.VERTICAL )

        gSizer4 = wx.GridSizer( 5, 1, 0, 0 )

        self.choose_pdb = wx.RadioBox( self, wx.ID_ANY, u"choose your pdb:", wx.DefaultPosition, wx.DefaultSize, choose_pdbChoices, 1, wx.RA_SPECIFY_COLS )
        self.choose_pdb.SetSelection( 10 )
        gSizer4.Add( self.choose_pdb, 0, wx.ALL, 5 )


        Results.Add( gSizer4, 1, wx.EXPAND, 5 )

        self.select = wx.Button( self, wx.ID_ANY, u"select", wx.DefaultPosition, wx.DefaultSize, 0 )
        Results.Add( self.select, 0, wx.ALL|wx.ALIGN_RIGHT, 5 )


        self.SetSizer( Results )
        self.Layout()

        self.Centre( wx.BOTH )

        # Connect Events
        self.select.Bind( wx.EVT_BUTTON, self.done )

    def __del__( self ):
        pass

    def done( self, event ):
            global row
            row = self.choose_pdb.GetSelection()

            self.Close()




def to_blast():
    if define_pdb == 1:
        blast()
        download_pdb()

to_blast()