import os
from bio_parse_pdb_fasta import project_dir
from bio_parse_pdb_fasta import spacer
from bio_retrieve_fasta import ccp4mg_loc
from bio_retrieve_fasta import output3
from bio_retrieve_fasta import output4

name = ""

code0 = ""
code1 = ""
code2 = ""
code3 = ""
code4 = ""
code5 = ""

def directories(name):
    if name == "myOutfile_abundance.pdb":
        os.makedirs(project_dir + spacer + "Pictures" + spacer + "abundance")
    if name == "myOutfile_amount.pdb":
        os.makedirs(project_dir + spacer + "Pictures" + spacer + "entropy")

def write_code(name):
    global code0
    global code1
    global code2
    global code3
    global code4
    global code5
    code0 = """MolData (
     filename = ['FULLPATH',
 '""" + name + """',
 '""" + project_dir + spacer + name + """'],
     name = 'filename',
     model_symmetry =  {
                    'MolDisp_visible' : 1  },
      )

MolDisp (
          colour_parameters =  {
                         'colour_mode' : 'bvalue'},
          style_parameters =  {
                         'style_mode' : 'BALLSTICK'  })

ColourSchemeManager(
          name = 'bvalue',
          ranges = [0.0, 0.0, 1.0, 0.0],
          colours = ['blue', 'blue', 'red', 'red'],
          colour_wheel_direction = 'clockwise',
          interpolate_mode = 'HSV'          ),"""

    code1= """import math
import numpy as np
theta = 0
phi = 0
psi = 0

MolData (
     filename = ['FULLPATH',
'"""+ name +"""',
'"""+ project_dir + spacer + name +"""'],
     model_symmetry =  {
                    'MolDisp_visible' : 1  },
      )

MolDisp (
          colour_parameters =  {
                         'colour_mode' : 'bvalue'},
          style_parameters =  {
                         'style_mode' : 'BALLSTICK'  })

ColourSchemeManager(
          name = 'bvalue',
          ranges = [0.0, 0.0, 1.0, 0.0],
          colours = ['blue', 'blue', 'red', 'red'],
          colour_wheel_direction = 'clockwise',
          interpolate_mode = 'HSV'          ),

matrix = np.array([[0,0,1],[0,1,0],[-1,0,0]])
if matrix.shape == (3, 3):
        try:
            theta = math.asin(- matrix[2][0])
            if math.cos(theta) != 0:
                try:
                    psi = math.atan2(matrix[2][1] / math.cos(theta), matrix[2][2] / math.cos(theta))
                except:
                    psi = math.atan(matrix[2][1])
                try:
                    phi = math.atan2(matrix[1][0] / math.cos(theta), matrix[0][0] / math.cos(theta))
                except:
                    phi = math.atan(matrix[1][0])
            else:
                phi = theta
                psi = math.arcsin(matrix[1][1])
        except:
            raise AlgebraicError, 'the matrix is not orthogonal'

# print theta
# print phi
# print psi

euler = ([int(phi), int(theta), int(psi)])
if len(euler) == 3:
    quat = np.array([math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.) - math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) - math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.)])
else:
     raise AlgebraicError, str(euler) + ' can not be converted to quaternion'

# print quat

view = View (
     orientation = quat
     )"""

    code2 = """import math
import numpy as np
theta = 0
phi = 0
psi = 0

MolData (
     filename = ['FULLPATH',
'"""+ name +"""',
'"""+ project_dir + spacer + name +"""'],
     model_symmetry =  {
                    'MolDisp_visible' : 1  },
      )

MolDisp (
          colour_parameters =  {
                         'colour_mode' : 'bvalue'},
          style_parameters =  {
                         'style_mode' : 'BALLSTICK'  })

ColourSchemeManager(
          name = 'bvalue',
          ranges = [0.0, 0.0, 1.0, 0.0],
          colours = ['blue', 'blue', 'red', 'red'],
          colour_wheel_direction = 'clockwise',
          interpolate_mode = 'HSV'          ),

matrix = np.array([[-1,0,0],[0,1,0],[0,0,-1]])
if matrix.shape == (3, 3):
        try:
            theta = math.asin(- matrix[2][0])
            if math.cos(theta) != 0:
                try:
                    psi = math.atan2(matrix[2][1] / math.cos(theta), matrix[2][2] / math.cos(theta))
                except:
                    psi = math.atan(matrix[2][1])
                try:
                    phi = math.atan2(matrix[1][0] / math.cos(theta), matrix[0][0] / math.cos(theta))
                except:
                    phi = math.atan(matrix[1][0])
            else:
                phi = theta
                psi = math.arcsin(matrix[1][1])
        except:
            raise AlgebraicError, 'the matrix is not orthogonal'

# print theta
# print phi
# print psi

euler = ([int(phi), int(theta), int(psi)])
if len(euler) == 3:
    quat = np.array([math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.) - math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) - math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.)])
else:
     raise AlgebraicError, str(euler) + ' can not be converted to quaternion'

# print quat

view = View (
     orientation = quat
     )
"""
    code3 = """import math
import numpy as np
theta = 0
phi = 0
psi = 0

MolData (
     filename = ['FULLPATH',
'"""+ name +"""',
'"""+ project_dir + spacer + name +"""'],
     model_symmetry =  {
                    'MolDisp_visible' : 1  },
      )

MolDisp (
          colour_parameters =  {
                         'colour_mode' : 'bvalue'},
          style_parameters =  {
                         'style_mode' : 'BALLSTICK'  })

ColourSchemeManager(
          name = 'bvalue',
          ranges = [0.0, 0.0, 1.0, 0.0],
          colours = ['blue', 'blue', 'red', 'red'],
          colour_wheel_direction = 'clockwise',
          interpolate_mode = 'HSV'          ),

matrix = np.array([[0,0,-1],[0,1,0],[1,0,0]])
if matrix.shape == (3, 3):
        try:
            theta = math.asin(- matrix[2][0])
            if math.cos(theta) != 0:
                try:
                    psi = math.atan2(matrix[2][1] / math.cos(theta), matrix[2][2] / math.cos(theta))
                except:
                    psi = math.atan(matrix[2][1])
                try:
                    phi = math.atan2(matrix[1][0] / math.cos(theta), matrix[0][0] / math.cos(theta))
                except:
                    phi = math.atan(matrix[1][0])
            else:
                phi = theta
                psi = math.arcsin(matrix[1][1])
        except:
            raise AlgebraicError, 'the matrix is not orthogonal'

# print theta
# print phi
# print psi

euler = ([int(phi), int(theta), int(psi)])
if len(euler) == 3:
    quat = np.array([math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.) - math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) - math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.)])
else:
     raise AlgebraicError, str(euler) + ' can not be converted to quaternion'

# print quat

view = View (
     orientation = quat
     )
"""
    code4 = """import math
import numpy as np
theta = 0
phi = 0
psi = 0

MolData (
     filename = ['FULLPATH',
'"""+ name +"""',
'"""+ project_dir + spacer + name +"""'],
     model_symmetry =  {
                    'MolDisp_visible' : 1  },
      )

MolDisp (
          colour_parameters =  {
                         'colour_mode' : 'bvalue'},
          style_parameters =  {
                         'style_mode' : 'BALLSTICK'  })

ColourSchemeManager(
          name = 'bvalue',
          ranges = [0.0, 0.0, 1.0, 0.0],
          colours = ['blue', 'blue', 'red', 'red'],
          colour_wheel_direction = 'clockwise',
          interpolate_mode = 'HSV'          ),

matrix = np.array([[1,0,0],[0,0,-1],[0,1,0]])
if matrix.shape == (3, 3):
        try:
            theta = math.asin(- matrix[2][0])
            if math.cos(theta) != 0:
                try:
                    psi = math.atan2(matrix[2][1] / math.cos(theta), matrix[2][2] / math.cos(theta))
                except:
                    psi = math.atan(matrix[2][1])
                try:
                    phi = math.atan2(matrix[1][0] / math.cos(theta), matrix[0][0] / math.cos(theta))
                except:
                    phi = math.atan(matrix[1][0])
            else:
                phi = theta
                psi = math.arcsin(matrix[1][1])
        except:
            raise AlgebraicError, 'the matrix is not orthogonal'

# print theta
# print phi
# print psi

euler = ([int(phi), int(theta), int(psi)])
if len(euler) == 3:
    quat = np.array([math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.) - math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) - math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.)])
else:
     raise AlgebraicError, str(euler) + ' can not be converted to quaternion'

# print quat

view = View (
     orientation = quat
     )
"""
    code5 = """import math
import numpy as np
theta = 0
phi = 0
psi = 0

MolData (
     filename = ['FULLPATH',
'"""+ name +"""',
'"""+ project_dir + spacer + name +"""'],
     model_symmetry =  {
                    'MolDisp_visible' : 1  },
      )

MolDisp (
          colour_parameters =  {
                         'colour_mode' : 'bvalue'},
          style_parameters =  {
                         'style_mode' : 'BALLSTICK'  })

ColourSchemeManager(
          name = 'bvalue',
          ranges = [0.0, 0.0, 1.0, 0.0],
          colours = ['blue', 'blue', 'red', 'red'],
          colour_wheel_direction = 'clockwise',
          interpolate_mode = 'HSV'          ),

matrix = np.array([[1,0,0],[0,-1,0],[0,0,-1]])
if matrix.shape == (3, 3):
        try:
            theta = math.asin(- matrix[2][0])
            if math.cos(theta) != 0:
                try:
                    psi = math.atan2(matrix[2][1] / math.cos(theta), matrix[2][2] / math.cos(theta))
                except:
                    psi = math.atan(matrix[2][1])
                try:
                    phi = math.atan2(matrix[1][0] / math.cos(theta), matrix[0][0] / math.cos(theta))
                except:
                    phi = math.atan(matrix[1][0])
            else:
                phi = theta
                psi = math.arcsin(matrix[1][1])
        except:
            raise AlgebraicError, 'the matrix is not orthogonal'

# print theta
# print phi
# print psi

euler = ([int(phi), int(theta), int(psi)])
if len(euler) == 3:
    quat = np.array([math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.) - math.sin(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.),
                                    math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.cos(euler[2] / 2.) + math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.sin(euler[2] / 2.),
                                    math.sin(euler[0] / 2.) * math.cos(euler[1] / 2.) * math.cos(euler[2] / 2.) - math.cos(euler[0] / 2.) * math.sin(euler[1] / 2.) * math.sin(euler[2] / 2.)])
else:
     raise AlgebraicError, str(euler) + ' can not be converted to quaternion'

# print quat

view = View (
     orientation = quat
     )
"""
    codes = [code0, code1, code2, code3, code4, code5]
    x = 0
    for code in codes:
        picture1 = open(project_dir + spacer + "Pictures" + spacer + "tmp" + spacer + "code" + str(x) + ".mgpic.py", "w")
        x = int(x)
        x += 1
        picture1.write(code)
        picture1.close()

def execute_code(code, name):
    if name == "myOutfile_abundance.pdb":
        cmd = ccp4mg_loc + spacer + "bin" + spacer + "ccp4mg"+ " -norestore -picture " + project_dir + spacer + "Pictures" + spacer + "tmp" + spacer + code + ".mgpic.py" + " -R " + project_dir + spacer + "Pictures" + spacer + "abundance" + spacer + code + ".png" + " -quit"
        os.system(cmd)
    if name == "myOutfile_amount.pdb":
        cmd = ccp4mg_loc + spacer + "bin" + spacer + "ccp4mg"+ " -norestore -picture " + project_dir + spacer + "Pictures" + spacer + "tmp" + spacer + code + ".mgpic.py" + " -R " + project_dir + spacer + "Pictures" + spacer + "entropy" + spacer + code + ".png" + " -quit"
        os.system(cmd)

def main():
    names = []
    global name
    os.makedirs(project_dir + spacer + "Pictures")
    os.makedirs(project_dir + spacer + "Pictures" + spacer + "tmp")
    if output3 == 1 and output4 == 1:
        names = ["myOutfile_abundance.pdb", "myOutfile_amount.pdb"]
    elif output3 == 1:
        names = ["myOutfile_abundance.pdb"]
    elif output4 == 1:
        names = ["myOutfile_amount.pdb"]

    for name_entrie in names:
        name = name_entrie
        directories(name)
        write_code(name)
        codes_str = ["code0", "code1", "code2", "code3", "code4", "code5"]
        for code in codes_str:

            execute_code(code, name)
main()