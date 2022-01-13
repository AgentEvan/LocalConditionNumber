'''
Description:
    Obtaining local condition number according to Sparse Matrix
    
Author: Xianwen Guo
Date: 2019/08/22

'''

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import inv as spInv
import os

def get_TimeGeneral(case_dir):
    file_name = case_dir+'/latestTime'
    time = np.loadtxt(file_name, dtype=int).item()
    return time

def get_FoamData_scalar(feat_str, data_dir, N=20193, skr=23):  
    filename = data_dir + feat_str
    dataframe = pd.read_table(filename, skiprows=skr, nrows=N, header=None, sep='\s+')
    dataframe = dataframe.apply(pd.to_numeric, errors='ignore')
    return dataframe

def get_FoamData_vector(feat_str, data_dir, N=20193, skr=23):
    filename = data_dir + feat_str
    dataframe = pd.read_table(filename, skiprows=skr, nrows=N, header=None, sep='\s+')
    dataframe[0] = dataframe[0].str.split('(', expand=True)[1]
    dataframe[2] = dataframe[2].str.split(')', expand=True)[0]
    dataframe = dataframe.apply(pd.to_numeric, errors='ignore')
    return dataframe

def get_FoamData_tensor(feat_str, data_dir, N=20193, skr=23):  
    filename = data_dir + feat_str
    dataframe = pd.read_table(filename, skiprows=skr, nrows=N, header=None, sep='\s+')
    dataframe[0] = dataframe[0].str.split('(', expand=True)[1]
    dataframe[8] = dataframe[8].str.split(')', expand=True)[0]
    dataframe = dataframe.apply(pd.to_numeric, errors='ignore')
    return dataframe

def get_FoamData_symmTensor(feat_str, data_dir, N=20193, skr=23):  
    filename = data_dir + feat_str
    dataframe = pd.read_table(filename, skiprows=skr, nrows=N, header=None, sep='\s+')
    dataframe[0] = dataframe[0].str.split('(', expand=True)[1]
    dataframe[5] = dataframe[5].str.split(')', expand=True)[0]
    dataframe = dataframe.apply(pd.to_numeric, errors='ignore')
    return dataframe

def get_SparseMatrix(fileName):
    A_sparse = np.loadtxt(fileName)
    row = np.array(A_sparse[:,0], dtype=int)
    col = np.array(A_sparse[:,1], dtype=int)
    data = A_sparse[:,2]
    A_spm = coo_matrix((data, (row, col)))
    return A_spm

def Header(var_name, var_type, time_name='0'):
    text_header = "/*--------------------------------*- C++ -*----------------------------------*\\"+'\n'
    text_header += "| =========                 |                                                 |\n"
    text_header += "| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
    text_header += "|  \\\    /   O peration     | Version:  v1806                                 |\n"
    text_header += "|   \\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n"
    text_header += "|    \\\/     M anipulation  |                                                 |\n"
    text_header += "\\*---------------------------------------------------------------------------*/\n"
    text_header += "FoamFile\n"
    text_header += "{\n"
    text_header += "    version     2.0;\n"
    text_header += "    format      ascii;\n"
    text_header += "    class       "
    text_header += var_type
    text_header += ";\n    location    "
    text_header += time_name
    text_header += ";\n    object      "
    text_header += var_name
    text_header += ";\n}\n"
    text_header += "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"
    return text_header


def Boundary_Scalar_Channel2D(value=0):
    bound_text = "\nboundaryField\n"
    bound_text += "{\n"
    bound_text += "INLET\n"
    bound_text += "{\n"
    bound_text += "    type            cyclic;\n"
    bound_text += "}\n"
    bound_text += "OUTLET\n"
    bound_text += "{\n"
    bound_text += "    type            cyclic;\n"
    bound_text += "}\n"
    bound_text += "WALL\n"
    bound_text += "{\n"
    bound_text += "    type            fixedValue;\n"
    bound_text += "    value           uniform "
    bound_text += str(value)
    bound_text += ";\n"
    bound_text += "}\n"
    bound_text += "CENTER\n"
    bound_text += "{\n"
    bound_text += "    type            symmetryPlane;\n"
    bound_text += "}\n"
    bound_text += "frontAndBackPlanes\n"
    bound_text += "{\n"
    bound_text += "    type            empty;\n"
    bound_text += "}\n"
    bound_text += "}\n\n\n"
    bound_text += "// ************************************************************************* //\n\n"
    return bound_text


def writeFoamScalar(varName, data):
    time = get_TimeGeneral('.')
    nCells = data.shape[0]
    f = open('./'+str(time)+'/'+varName, 'w')
    f.write(Header(varName, 'volScalarField'))
    f.write("dimensions      [0 0 0 0 0 0 0];\n\n\n")
    f.write("internalField   nonuniform List<scalar>\n")
    f.write(str(nCells)+'\n')
    f.write('(\n')
    for i in range(nCells):
        f.write(("%.8g" % data[i])+'\n')
    f.write(')\n')
    f.write(';\n')
    f.write(Boundary_Scalar_Channel2D())  
    f.close()


if (__name__ == '__main__'):
    
    print('Loading Sparse Matrix ...\n')
    A_Sparse_exp = get_SparseMatrix('./A_Sparse_fromFull_exp')
    A_Sparse_imp = get_SparseMatrix('./A_Sparse_fromFull_imp')
    
    print('Calculate Matirx Inverse ...\n')
    Ainv1 = spInv(A_Sparse_exp.tocsc()).todense()
    Ainv2 = spInv(A_Sparse_imp.tocsc()).todense()
    print('Calculate Matirx Inverse Done.\n')
    
    N = Ainv1.shape[0]
    time = str(get_TimeGeneral('.'))
    cellV = np.loadtxt('./cellVolume')
    divRS = get_FoamData_vector('divRS_DNS', './'+time+'/', N).values
    U = get_FoamData_vector('U_DNS', './'+time+'/', N).values
    Uinf = 0
    vol = 0
    for i in range(N):
        # divRS[i,:] *= cellV[i]
        Uinf += cellV[i]*U[i,0]
        vol += cellV[i]
    Uinf = Uinf/vol
    print('Uinf is ', str(Uinf), '\n')
    
    rj1 = np.zeros(N)
    rj2 = np.zeros(N)

    Kj1 = np.zeros(N)
    Kj2 = np.zeros(N)
    for i in range(N):
        rj1[i] = np.linalg.norm(Ainv1[i,:])
        rj2[i] = np.linalg.norm(Ainv2[i,:])

        Kj1[i] = rj1[i]*np.linalg.norm(divRS[:,0])/(Uinf+1e-12)
        Kj2[i] = rj2[i]*np.linalg.norm(divRS[:,0])/(Uinf+1e-12)
        
    print('Calculate Local Condition Number Done ...\n')    
    writeFoamScalar('rj_exp', rj1)
    writeFoamScalar('rj_imp', rj2)
    writeFoamScalar('Kj_exp', Kj1)
    writeFoamScalar('Kj_imp', Kj2)
    print('Write Local Condition Number Done ...')
    
    print('Calculate Volume-averaged Condition Number ...\n')
    Kjv1 = 0
    Kjv2 = 0
    cV = 0
    for i in range(N):
        Kjv1 += Kj1[i] * cellV[i]
        Kjv2 += Kj2[i] * cellV[i]
        cV += cellV[i]
    Kjv1 /= cV
    Kjv2 /= cV
    
    f = open('./Kj_vol_exp', 'w')
    f.write("%.8g" % Kjv1)    
    f.close()
    
    f = open('./Kj_vol_imp', 'w')
    f.write("%.8g" % Kjv2)    
    f.close()
    
    print('Finish.')





