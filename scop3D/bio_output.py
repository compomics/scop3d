__author__ = 'stijn'
# make the output GUI

import wx
import wx.xrc
import wx.richtext
import os
from bio_retrieve_fasta import output2
from bio_retrieve_fasta import output3
from bio_retrieve_fasta import output4
from bio_blast_pdb_fasta import project_dir
from bio_retrieve_fasta import project_name
from bio_retrieve_fasta import ccp4mg_loc
from bio_parse_pdb_fasta import consensus
from bio_pictures_pdb_fasta import spacer
from sys import platform as _platform



consensus = consensus
file_info = open(project_dir + spacer + "info.txt", "r")
errorpictures = """
	Picture error:
	There are no pictures to be shown.
	2 things can have happend.
	If you're running Scop3D on a mac,
	we have not been able to make ccp4mg work from command line on this operating system.

	If not make sure ccp4mg is correctly installed.
	Try running it separately from Scop3D.
	If you cannot make it work we strongly recommend you try a visualization program like pymol.
	and visualize myOutfile_abundance.pdb and myOutfile_amount.pdb outside Scop3D.

	Please let us know all problems so we can improve Scop3D.
	We are sorry for the inconvenience.
	"""
# get all b-values of procentual abundance
if output3 == 1:
    file_abundance=open(project_dir + spacer + "myOutfile_abundance.pdb", "r")
    abundance = "Lines of pdb containing C-alpha coordinates and the abundance: \n"
    for line in file_abundance:
        if line.startswith("ATOM") and "CA" in line:
            abundance += line

# get all b-values of entropy
if output4 == 1:
    file_amount=open(project_dir + spacer + "myOutfile_amount.pdb", "r")
    amount = "Lines of pdb containing C-alpha coordinates and the entropy: \n"
    for line in file_amount:
        if line.startswith("ATOM") and "CA" in line:
            amount += line

# open info file
info = file_info.read()

# get consensus-sequence
from Bio import SeqIO
alignmentid = "record id:"
alignmentseq = "record sequence:"
instances = []
for record in SeqIO.parse(open(project_dir + spacer + "alignments" + spacer +"alignment_" + project_name + ".clw"), "clustal") :
    alignmentid += "\n" + record.id
    alignmentseq += "\n" + record.seq
alignmentid += "\n\n consensusseq"
alignmentseq += "\n\n" + consensus

class Scop3D ( wx.Frame ):

	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"Scop3D - Output", pos = wx.DefaultPosition, size = wx.Size( 569,423 ), style = wx.DEFAULT_FRAME_STYLE )

		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

		bSizer1 = wx.BoxSizer( wx.VERTICAL )

		self.notebook_1 = wx.Notebook( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )
		self.notebook_1_pane_1 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		gSizer2 = wx.GridSizer( 0, 1, 0, 0 )

		self.m_richText1 = wx.richtext.RichTextCtrl( self.notebook_1_pane_1, wx.ID_ANY, info, wx.DefaultPosition, wx.DefaultSize, 0|wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER|wx.WANTS_CHARS )
		gSizer2.Add( self.m_richText1, 1, wx.EXPAND |wx.ALL, 5 )


		self.notebook_1_pane_1.SetSizer( gSizer2 )
		self.notebook_1_pane_1.Layout()
		gSizer2.Fit( self.notebook_1_pane_1 )
		self.notebook_1.AddPage( self.notebook_1_pane_1, u"project", True )

		self.notebook_1_pane_2 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		gbSizer2 = wx.GridBagSizer( 0, 0 )
		gbSizer2.SetFlexibleDirection( wx.BOTH )
		gbSizer2.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

		self.m_richText6 = wx.TextCtrl( self.notebook_1_pane_2, wx.ID_ANY, str(alignmentid) , wx.DefaultPosition, wx.DefaultSize, wx.TE_MULTILINE|wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER|wx.WANTS_CHARS )
		self.m_richText6.SetFont( wx.Font( 11, 76, 90, 90, False, "Courier 10 Pitch" ) )
		self.m_richText6.SetMinSize( wx.Size( 200,-1 ) )
		self.m_richText6.SetMaxSize( wx.Size( 300,-1 ) )

		gbSizer2.Add( self.m_richText6, wx.GBPosition( 0, 0 ), wx.GBSpan( 1, 1 ), wx.EXPAND |wx.ALL, 5 )

		self.m_richText71 = wx.TextCtrl( self.notebook_1_pane_2, wx.ID_ANY, "   " + str(alignmentseq) + "   ", wx.DefaultPosition, wx.DefaultSize, wx.TE_MULTILINE|wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER|wx.WANTS_CHARS )
		# self.m_richText71.SetMinSize( wx.Size( 1,-1 ) )
		self.m_richText71.SetFont( wx.Font( 11, 76, 90, 90, False, "Courier 10 Pitch" ) )
		gbSizer2.Add( self.m_richText71, wx.GBPosition( 0, 1 ), wx.GBSpan( 1, 1 ), wx.EXPAND |wx.ALL, 5 )

		gbSizer2.AddGrowableRow( 0 )
		gbSizer2.AddGrowableCol(1)
		# gbSizer2.AddGrowableRow( 5 )
		self.notebook_1_pane_2.SetSizer( gbSizer2 )
		self.notebook_1_pane_2.Layout()
		gbSizer2.Fit( self.notebook_1_pane_2 )
		self.notebook_1.AddPage( self.notebook_1_pane_2, u"alignment", False )

		if output3 ==1:
		    self.notebook_1_pane_3 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		    gSizer4 = wx.GridSizer( 0, 1, 0, 0 )


		    self.m_richText3 = wx.richtext.RichTextCtrl( self.notebook_1_pane_3, wx.ID_ANY, str(abundance), wx.DefaultPosition, wx.DefaultSize, 0|wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER|wx.WANTS_CHARS )
		    self.m_richText3.SetFont( wx.Font( 11, 76, 90, 90, False, "Courier 10 Pitch" ) )
		    gSizer4.Add( self.m_richText3, 1, wx.EXPAND |wx.ALL, 5 )


		    self.notebook_1_pane_3.SetSizer( gSizer4 )
		    self.notebook_1_pane_3.Layout()
		    gSizer4.Fit( self.notebook_1_pane_3 )
		    self.notebook_1.AddPage( self.notebook_1_pane_3, u"pdb_abundance", False )
		if output4 == 1:
		    self.notebook_1_pane_4 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		    gSizer5 = wx.GridSizer( 0, 1, 0, 0 )

		    self.m_richText4 = wx.richtext.RichTextCtrl( self.notebook_1_pane_4, wx.ID_ANY, str(amount), wx.DefaultPosition, wx.DefaultSize, 0|wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER|wx.WANTS_CHARS )
		    self.m_richText4.SetFont( wx.Font( 11, 76, 90, 90, False, "Courier 10 Pitch" ) )
		    gSizer5.Add( self.m_richText4, 1, wx.EXPAND |wx.ALL, 5 )


		    self.notebook_1_pane_4.SetSizer( gSizer5 )
		    self.notebook_1_pane_4.Layout()
		    gSizer5.Fit( self.notebook_1_pane_4 )
		    self.notebook_1.AddPage( self.notebook_1_pane_4, u"pdb_entropy", False )
		if output2 == 1 and os.path.exists(project_dir + spacer + "sequence_logo.png")==True:
		    self.notebook_1_pane_5 = wx.ScrolledWindow( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.VSCROLL )
		    self.notebook_1_pane_5.SetScrollRate( 5, 5 )
		    gSizer6 = wx.GridSizer(0, 1, 0, 0)


		    self.png = wx.Image(project_dir + spacer + "sequence_logo.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_bitmap2 = wx.StaticBitmap(self.notebook_1_pane_5, -1, self.png, (10, 5), (self.png.GetWidth(), self.png.GetHeight()))
		    gSizer6.Add( self.m_bitmap2, 0, wx.ALL, 5 )
		    self.notebook_1_pane_5.SetSizer( gSizer6)
		    self.notebook_1_pane_5.Layout()
		    gSizer6.Fit( self.notebook_1_pane_5)
		    self.notebook_1.AddPage( self.notebook_1_pane_5, u"seq_logo", False )
		elif output2 == 1:
		    self.notebook_1_pane_5 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		    gSizer6 = wx.GridSizer( 0, 1, 0, 0 )
		    self.m_seqlogo_fail = wx.richtext.RichTextCtrl( self.notebook_1_pane_5, wx.ID_ANY, str("sequence logo can be found here: " + project_dir + spacer + "sequence_logo.eps "), wx.DefaultPosition, wx.DefaultSize, 0|wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER|wx.WANTS_CHARS )
		    gSizer6.Add( self.m_seqlogo_fail, 1, wx.EXPAND |wx.ALL, 5 )
		    self.notebook_1_pane_5.SetSizer( gSizer6 )
		    self.notebook_1_pane_5.Layout()
		    gSizer6.Fit( self.notebook_1_pane_5 )
		    self.notebook_1.AddPage( self.notebook_1_pane_5, u"seq_logo", False )



		if _platform != "darwin" and os.path.exists(project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + "code1.png"):
		    self.pictures_abundance = wx.ScrolledWindow( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		    self.pictures_abundance.SetScrollRate( 5, 5 )
		    gSizer7 = wx.GridSizer( 4, 2, 0, 0 )
		    self.code0_abu = wx.Image(project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + "code0.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code0_abu = wx.StaticBitmap(self.pictures_abundance, -1, self.code0_abu, (10, 5), (self.code0_abu.GetWidth(), self.code0_abu.GetHeight()))
		    gSizer7.Add( self.m_code0_abu, 0, wx.ALL, 5 )

		    self.code1_abu = wx.Image(project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + "code1.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code1_abu = wx.StaticBitmap(self.pictures_abundance, -1, self.code1_abu, (10, 5), (self.code1_abu.GetWidth(), self.code1_abu.GetHeight()))
		    gSizer7.Add( self.m_code1_abu, 0, wx.ALL, 5 )

		    self.code2_abu = wx.Image(project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + "code2.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code2_abu = wx.StaticBitmap(self.pictures_abundance, -1, self.code2_abu, (10, 5), (self.code2_abu.GetWidth(), self.code2_abu.GetHeight()))
		    gSizer7.Add( self.m_code2_abu, 0, wx.ALL, 5 )

		    self.code3_abu = wx.Image(project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + "code3.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code3_abu = wx.StaticBitmap(self.pictures_abundance, -1, self.code3_abu, (10, 5), (self.code3_abu.GetWidth(), self.code3_abu.GetHeight()))
		    gSizer7.Add( self.m_code3_abu, 0, wx.ALL, 5 )

		    self.code4_abu = wx.Image(project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + "code4.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code4_abu = wx.StaticBitmap(self.pictures_abundance, -1, self.code4_abu, (10, 5), (self.code4_abu.GetWidth(), self.code4_abu.GetHeight()))
		    gSizer7.Add( self.m_code4_abu, 0, wx.ALL, 5 )

		    self.code5_abu = wx.Image(project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + "code5.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code5_abu = wx.StaticBitmap(self.pictures_abundance, -1, self.code5_abu, (10, 5), (self.code5_abu.GetWidth(), self.code5_abu.GetHeight()))
		    gSizer7.Add( self.m_code5_abu, 0, wx.ALL, 5 )

		    gSizer7.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )

		    self.m_ccp4mg_abundance = wx.ToggleButton( self.pictures_abundance, wx.ID_ANY, u"see in ccp4mg", wx.DefaultPosition, wx.DefaultSize, 0 )
		    gSizer7.Add( self.m_ccp4mg_abundance, 0, wx.ALL|wx.ALIGN_CENTER|wx.ALIGN_CENTER, 5 )

		    self.pictures_abundance.SetSizer( gSizer7 )
		    self.pictures_abundance.Layout()
		    self.m_ccp4mg_abundance.Bind( wx.EVT_TOGGLEBUTTON, self.ccp4mg_abundance )

		    gSizer7.Fit( self.pictures_abundance )
		    self.notebook_1.AddPage( self.pictures_abundance, u"pictures abundance", False )

		    self.pictures_entropy = wx.ScrolledWindow( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		    self.pictures_entropy.SetScrollRate( 5, 5 )
		    gSizer8 = wx.GridSizer( 4, 2, 0, 0 )
		    self.code0_ent = wx.Image(project_dir + spacer + "Pictures" + spacer + "entropy" + spacer + "code0.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code0_ent = wx.StaticBitmap(self.pictures_entropy, -1, self.code0_ent, (10, 5), (self.code0_ent.GetWidth(), self.code0_ent.GetHeight()))
		    gSizer8.Add( self.m_code0_ent, 0, wx.ALL, 5 )

		    self.code1_ent = wx.Image(project_dir + spacer + "Pictures" + spacer + "entropy" + spacer + "code1.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code1_ent = wx.StaticBitmap(self.pictures_entropy, -1, self.code1_ent, (10, 5), (self.code1_ent.GetWidth(), self.code1_ent.GetHeight()))
		    gSizer8.Add( self.m_code1_ent, 0, wx.ALL, 5 )

		    self.code2_ent = wx.Image(project_dir + spacer + "Pictures" + spacer + "entropy" + spacer + "code2.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code2_ent = wx.StaticBitmap(self.pictures_entropy, -1, self.code2_ent, (10, 5), (self.code2_ent.GetWidth(), self.code2_ent.GetHeight()))
		    gSizer8.Add( self.m_code2_ent, 0, wx.ALL, 5 )

		    self.code3_ent = wx.Image(project_dir + spacer + "Pictures" + spacer + "entropy" + spacer + "code3.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code3_ent = wx.StaticBitmap(self.pictures_entropy, -1, self.code3_ent, (10, 5), (self.code3_ent.GetWidth(), self.code3_ent.GetHeight()))
		    gSizer8.Add( self.m_code3_ent, 0, wx.ALL, 5 )

		    self.code4_ent = wx.Image(project_dir + spacer + "Pictures" + spacer + "entropy" + spacer + "code4.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code4_ent = wx.StaticBitmap(self.pictures_entropy, -1, self.code4_ent, (10, 5), (self.code4_ent.GetWidth(), self.code4_ent.GetHeight()))
		    gSizer8.Add( self.m_code4_ent, 0, wx.ALL, 5 )

		    self.code5_ent = wx.Image(project_dir + spacer + "Pictures" + spacer + "entropy" + spacer + "code5.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
		    self.m_code5_ent = wx.StaticBitmap(self.pictures_entropy, -1, self.code5_ent, (10, 5), (self.code5_ent.GetWidth(), self.code5_ent.GetHeight()))
		    gSizer8.Add( self.m_code5_ent, 0, wx.ALL, 5 )

		    gSizer8.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )

		    self.m_ccp4mg_entropy = wx.ToggleButton( self.pictures_entropy, wx.ID_ANY, u"see in ccp4mg", wx.DefaultPosition, wx.DefaultSize, 0 )
		    gSizer8.Add( self.m_ccp4mg_entropy, 0, wx.ALL|wx.ALIGN_CENTER|wx.ALIGN_CENTER, 5 )

		    self.pictures_entropy.SetSizer( gSizer8 )
		    self.pictures_entropy.Layout()

		    self.m_ccp4mg_entropy.Bind( wx.EVT_TOGGLEBUTTON, self.ccp4mg_entropy )

		    gSizer8.Fit( self.pictures_entropy )
		    self.notebook_1.AddPage( self.pictures_entropy, u"pictures entropy", False )

		if _platform == "darwin" or os.path.exists(project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + "code1.png") == False:
		    self.pictures_error = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		    gSizer200 = wx.GridSizer( 0, 1, 0, 0 )

		    self.m_richText100 = wx.richtext.RichTextCtrl( self.pictures_error, wx.ID_ANY, errorpictures, wx.DefaultPosition, wx.DefaultSize, 0|wx.VSCROLL|wx.HSCROLL|wx.NO_BORDER|wx.WANTS_CHARS )
		    gSizer200.Add( self.m_richText100, 1, wx.EXPAND |wx.ALL, 5 )


		    self.pictures_error.SetSizer( gSizer200 )
		    self.pictures_error.Layout()

		    gSizer200.Fit( self.pictures_error )
		    self.notebook_1.AddPage( self.pictures_error, u"pictures", False )
		# self.notebook_1_pane_6 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		# self.notebook_1.AddPage( self.notebook_1_pane_6, u"6", False )
		# self.notebook_1_pane_7 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		# self.notebook_1.AddPage( self.notebook_1_pane_7, u"7", False )
		# self.notebook_1_pane_8 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		# self.notebook_1.AddPage( self.notebook_1_pane_8, u"8", False )
		# self.notebook_1_pane_9 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		# self.notebook_1.AddPage( self.notebook_1_pane_9, u"9", False )

		bSizer1.Add( self.notebook_1, 1, wx.EXPAND, 0 )


		self.SetSizer( bSizer1 )
		self.Layout()

	def __del__( self ):
		pass

	if _platform == "win32":
	    def ccp4mg_entropy (self, event):
                ccmp4g_loc = ccp4mg_loc + spacer + "winccp4mg.exe"
                code_in = project_dir + spacer + "Pictures" + spacer + "tmp-abundance" + spacer + "code0" + ".mgpic.py"
                pic_out = project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + "code0" + ".png"
                os.system("\"" + ccmp4g_loc + "\"" + " -norestore -picture " + code_in + " -R " + pic_out)
        def ccp4mg_abundance (self, event):
                ccmp4g_loc = ccp4mg_loc + spacer + "winccp4mg.exe"
                code_in = project_dir + spacer + "Pictures" + spacer + "tmp-entropy" + spacer + "code0" + ".mgpic.py"
                pic_out = project_dir + spacer + "Pictures" + spacer + "entropy" + spacer + "code0" + ".png"
                os.system("\"" + ccmp4g_loc + "\"" + " -norestore -picture " + code_in + " -R " + pic_out)
	if _platform == "linux" or _platform == "linux2":
	    def ccp4mg_entropy (self, event):
	        cmd = ccp4mg_loc + spacer + "bin" + spacer + "ccp4mg"+ " -norestore -picture " + project_dir + spacer + "Pictures" + spacer + "tmp-entropy" + spacer + 'code0' + ".mgpic.py" + " -R " + project_dir + spacer + "Pictures" + spacer + "entropy" + spacer + 'code0' + ".png"
	        os.system(cmd)

	    def ccp4mg_abundance (self, event):
	        cmd = ccp4mg_loc + spacer + "bin" + spacer + "ccp4mg"+ "  -norestore -picture " + project_dir + spacer + "Pictures" + spacer + "tmp-abundance" + spacer + 'code0' + ".mgpic.py" + " -R " + project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + 'code0' + ".png"
	        os.system(cmd)




def main():
    app = wx.App(0)
    MainApp = Scop3D(None)
    MainApp.Show()
    MainApp.Maximize()
    app.MainLoop()

main()
