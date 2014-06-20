#input: pdb-file
#output: sequence of the structure within the pdbfile

from Bio.PDB import PDBList
import os
import shutil
import wx
import wx.xrc
from sys import platform as _platform
#import bio_blast_pdb_fasta
#import bio_parse_pdb_fasta
#import bio_conservation_mapping_fasta



spacer = ""
fasta_dir = ""
project_dir = ""
project_name = ""
define_pdb = 0
pdbfile = ""
cons_tresh = ""
win_platform = ""
muscle_loc = ""
output1 = 0
output2 = 0
output3 = 0
output4 = 0


def operating_sys():
    #check for linux
    global spacer
    global win_platform
    if _platform == "linux" or _platform == "linux2":
        spacer = "/"
    elif _platform == "darwin":
        spacer = "/"
    #check if windows
    elif _platform == "win32":
        spacer = "\\"
        win_platform = "win"



class Scop3D ( wx.Frame ):

	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = u"Scop3D", pos = wx.DefaultPosition, size = wx.Size( 482,254 ), style = wx.DEFAULT_FRAME_STYLE )

		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )

		bSizer14 = wx.BoxSizer( wx.VERTICAL )

		self.notebook_1 = wx.Notebook( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )
		self.project = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		gSizer12 = wx.GridSizer( 4, 2, 0, 0 )

		self.lbl_project_name = wx.StaticText( self.project, wx.ID_ANY, u"project:", wx.DefaultPosition, wx.DefaultSize, 0)
		self.lbl_project_name.Wrap( -1 )
		gSizer12.Add( self.lbl_project_name, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.project_name = wx.TextCtrl( self.project, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer12.Add( self.project_name, 0, wx.ALL, 5 )

		self.lbl_fasta_file = wx.StaticText( self.project, wx.ID_ANY, u"Fasta file:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.lbl_fasta_file.Wrap( -1 )
		gSizer12.Add( self.lbl_fasta_file, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.fastafile = wx.FilePickerCtrl( self.project, wx.ID_ANY, wx.EmptyString, u"Select a file", u"*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		gSizer12.Add( self.fastafile, 0, wx.ALL, 5 )

		self.lbl_output_loc = wx.StaticText( self.project, wx.ID_ANY, u"Output location:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.lbl_output_loc.Wrap( -1 )
		gSizer12.Add( self.lbl_output_loc, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.output_loc = wx.DirPickerCtrl( self.project, wx.ID_ANY, wx.EmptyString, u"Select a folder", wx.DefaultPosition, wx.DefaultSize, wx.DIRP_DEFAULT_STYLE )
		gSizer12.Add( self.output_loc, 0, wx.ALL, 5 )


		gSizer12.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		gSizer12.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		self.project.SetSizer( gSizer12 )
		self.project.Layout()
		gSizer12.Fit( self.project )
		self.notebook_1.AddPage( self.project, u"project", True )
		self.pdb = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		gSizer2 = wx.GridSizer( 4, 2, 0, 0 )

		self.lbl_pdb = wx.StaticText( self.pdb, wx.ID_ANY, u"Define pdb entry:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.lbl_pdb.Wrap( -1 )
		gSizer2.Add( self.lbl_pdb, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

		define_pdbChoices = [ u"Use own pdb", u"BLAST" ]
		self.define_pdb = wx.Choice( self.pdb, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, define_pdbChoices, 0 )
		self.define_pdb.SetSelection( 0 )
		gSizer2.Add( self.define_pdb, 0, wx.ALL, 5 )

		self.lbl_pdb_loc = wx.StaticText( self.pdb, wx.ID_ANY, u"Location pdb", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.lbl_pdb_loc.Wrap( -1 )
		gSizer2.Add( self.lbl_pdb_loc, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5 )

		self.pdb_loc = wx.FilePickerCtrl( self.pdb, wx.ID_ANY, wx.EmptyString, u"Select a file", u"*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		gSizer2.Add( self.pdb_loc, 0, wx.ALL, 5 )


		gSizer2.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		gSizer2.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		gSizer2.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		gSizer2.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		self.pdb.SetSizer( gSizer2 )
		self.pdb.Layout()
		gSizer2.Fit( self.pdb )
		self.notebook_1.AddPage( self.pdb, u"pdb", False )
		self.settings = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		gSizer5 = wx.GridSizer( 4, 2, 0, 0 )

		self.lbl_cons_tresh = wx.StaticText( self.settings, wx.ID_ANY, u"treshold consensusseq:", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.lbl_cons_tresh.Wrap( -1 )
		gSizer5.Add( self.lbl_cons_tresh, 0, wx.ALL, 5 )

		self.cons_tresh = wx.TextCtrl( self.settings, wx.ID_ANY, u"0.7", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.cons_tresh.SetMaxLength( 4 )
		gSizer5.Add( self.cons_tresh, 0, wx.ALL, 5 )


		if win_platform == "win":
		    self.lbl_muscle_loc = wx.StaticText( self.settings, wx.ID_ANY, u"location muscle.exe:", wx.DefaultPosition, wx.DefaultSize, 0 )
		    self.lbl_muscle_loc.Wrap( -1 )
		    gSizer5.Add( self.lbl_muscle_loc, 0, wx.ALL, 5 )

		    self.muscle_loc = wx.FilePickerCtrl( self.settings, wx.ID_ANY, wx.EmptyString, u"Select muscle.exe", u"*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		    gSizer5.Add( self.muscle_loc, 0, wx.ALL, 5 )


		gSizer5.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		gSizer5.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		gSizer5.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		gSizer5.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )


		self.settings.SetSizer( gSizer5 )
		self.settings.Layout()
		gSizer5.Fit( self.settings )
		self.notebook_1.AddPage( self.settings, u"settings", False )
		self.output = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		gSizer3 = wx.GridSizer( 0, 2, 0, 0 )

		self.lbl_output = wx.StaticText( self.output, wx.ID_ANY, u"Output options", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.lbl_output.Wrap( -1 )
		gSizer3.Add( self.lbl_output, 0, wx.ALL, 5 )

		# self.output1 = wx.CheckBox( self.output, wx.ID_ANY, u"output1", wx.DefaultPosition, wx.DefaultSize, 0 )
		# gSizer3.Add( self.output1, 0, wx.ALL, 5 )
        #
        #
		# gSizer3.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )

		self.output2 = wx.CheckBox( self.output, wx.ID_ANY, u"Sequence logo", wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer3.Add( self.output2, 0, wx.ALL, 5 )


		gSizer3.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )

		self.output3 = wx.CheckBox( self.output, wx.ID_ANY, u"pdb abundance", wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer3.Add( self.output3, 0, wx.ALL, 5 )


		gSizer3.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )

		self.output4 = wx.CheckBox( self.output, wx.ID_ANY, u"pdb Entropy", wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer3.Add( self.output4, 0, wx.ALL, 5 )

		gSizer3.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )
		gSizer3.AddSpacer( ( 0, 0), 1, wx.EXPAND, 5 )

		self.output.SetSizer( gSizer3 )
		self.output.Layout()
		gSizer3.Fit( self.output )
		self.notebook_1.AddPage( self.output, u"output", False )
		# self.notebook_1_pane_5 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		# self.notebook_1.AddPage( self.notebook_1_pane_5, u"5", False )
		# self.notebook_1_pane_6 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		# self.notebook_1.AddPage( self.notebook_1_pane_6, u"6", False )
		# self.notebook_1_pane_7 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		# self.notebook_1.AddPage( self.notebook_1_pane_7, u"7", False )
		# self.notebook_1_pane_8 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		# self.notebook_1.AddPage( self.notebook_1_pane_8, u"8", False )
		# self.notebook_1_pane_9 = wx.Panel( self.notebook_1, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		# self.notebook_1.AddPage( self.notebook_1_pane_9, u"9", False )

		bSizer14.Add( self.notebook_1, 4, wx.EXPAND, 0 )

		gSizer4 = wx.GridSizer( 0, 1, 0, 0 )

		self.submit = wx.ToggleButton( self, wx.ID_ANY, u"run", wx.DefaultPosition, wx.DefaultSize, 0 )
		gSizer4.Add( self.submit, 0, wx.ALL|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL, 5 )


		bSizer14.Add( gSizer4, 1, wx.EXPAND, 5 )


		self.SetSizer( bSizer14 )
		self.Layout()

		# Connect Events
		self.submit.Bind( wx.EVT_TOGGLEBUTTON, self.done )

	def __del__( self ):
		pass


	# Virtual event handlers, overide them in your derived class
	def done( self, event ):
            #getting location for output
            global project_name
            global project_dir
            global fasta_dir
            global fasta_name
            global define_pdb
            global pdbfile
            global cons_tresh
            global muscle_loc
            # global output1
            global output2
            global output3
            global output4

            project_name = self.project_name.GetValue()
            output_loc = self.output_loc.GetPath()
            fastafile = self.fastafile.GetPath()
            if win_platform == "win":
                muscle_loc = self.muscle_loc.GetPath()
            project_dir = output_loc + spacer + project_name
            define_pdb = self.define_pdb.GetCurrentSelection()
            pdbfile = self.pdb_loc.GetPath()
            cons_tresh = float(self.cons_tresh.GetValue())
            # output1 = self.output1.Get3StateValue()
            output2 = self.output2.Get3StateValue()
            output3 = self.output3.Get3StateValue()
            output4 = self.output4.Get3StateValue()
            print project_name
            print define_pdb
            #making project directory
            print "making directory for", project_name
            print "at this location", project_dir
            os.makedirs(project_dir)

            #verplaatst fasta file naar project map onder project/fasta/
            fasta_dir = project_dir + spacer + "fasta"
            os.makedirs(fasta_dir)
            fasta_name = project_name + ".fasta"
            fasta_dir = fasta_dir + spacer + fasta_name
            shutil.copyfile(fastafile, fasta_dir)
            if pdbfile:
                shutil.copyfile(pdbfile, project_dir + spacer + "tmppdb.pdb")
                pdbfile = project_dir + spacer + "tmppdb.pdb"

            write_info()
            self.Close()


def write_info():
    info_file= open(project_dir + spacer +"info.txt", "a")
    info_file.write("project: " + project_name + "\n")
    info_file.write("project_location: " + project_dir + "\n")
    info_file.write("fasta_location: " + fasta_dir + "\n")
    info_file.write("fasta_filename: " + fasta_name + "\n")
    info_file.close()


def main():
    operating_sys()
    app = wx.App(0)
    MainApp = Scop3D(None)
    MainApp.Show()
    app.MainLoop()
    #bio_conservation_mapping_fasta.conservation_mapping()
    #bio_blast_pdb_fasta.to_blast()
    #bio_parse_pdb_fasta.parse_pdb()


main()
