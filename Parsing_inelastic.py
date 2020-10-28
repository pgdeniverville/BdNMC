import numpy as np
import math
from scipy.interpolate import interp1d
from scipy.interpolate import SmoothBivariateSpline

def gen_cross_section_kinetic(mv,eps,alpha_d):
    with open("data/xs_dis_Ap.txt") as infile:
        array1=infile.read()
    array2=array1.split("\n")[:-1]
    array3=[line.split() for line in array2]
    array4=np.array([[float(line[0].split('_')[1].split('p')[1]),float(line[0].split('_')[2].split('.')[0][1:]),float(line[1])] for line in array3])
    alpha_d_set = 1/4.0/math.pi
    epsilon_set = 1e-3
    DIS_Func=SmoothBivariateSpline(array4[:,0],array4[:,1],array4[:,2])
    with open("data/xs_dis_Ap_neutron.txt") as infile:
        array1=infile.read()
    array2=array1.split("\n")[:-1]
    array3=[line.split() for line in array2]
    array4=np.array([[float(line[0].split('_')[1].split('p')[1]),float(line[0].split('_')[2].split('.')[0][1:]),float(line[1])] for line in array3])
    alpha_d_set = 1/4.0/math.pi
    epsilon_set = 1e-3
    DIS_Func2=SmoothBivariateSpline(array4[:,0],array4[:,1],array4[:,2])
    return np.array([[En,DIS_Func(mv,En)*(eps/1e-3)**2*alpha_d/alpha_d_set,DIS_Func2(mv,En)*(eps/1e-3)**2*alpha_d/alpha_d_set] for En in range(3,301,1)])


def gen_cross_section_baryonic(mv,alpha_d):
    with open("data/xs_dis_B.txt") as infile:
        array1=infile.read()
    array2=array1.split("\n")[:-1]
    array3=[line.split() for line in array2]
    array4=np.array([[float(line[0].split('_')[1].split('p')[1]),float(line[0].split('_')[2].split('.')[0][1:]),float(line[1])] for line in array3])
    alpha_d_set = 1/4.0/math.pi
    DIS_Func=SmoothBivariateSpline(array4[:,0],array4[:,1],array4[:,2],kx=1,ky=1)
    return np.array([[En,DIS_Func(mv,En)*(alpha_d/alpha_d_set)**2,DIS_Func(mv,En)*(alpha_d/alpha_d_set)**2] for En in range(3,301,1)])


def gen_cross_section(mv,eps=1e-3,alpha_d=0.1,model="Inelastic_Nucleon_Scattering",output="data/DIS.dat"):
    if(model=="Inelastic_Nucleon_Scattering"):
        arr=gen_cross_section_kinetic(mv,eps,alpha_d)
    elif(model=="Inelastic_Nucleon_Scattering_Baryonic"):
        arr=gen_cross_section_baryonic(mv,alpha_d)
    arr=np.append([[eps,alpha_d,0]],arr,axis=0)
    np.savetxt(output,arr)
 
#gen_cross_section(0.1,eps=1e-3,alpha_d=0.1,model="Inelastic_Nucleon_Scattering",output="data/DIS.dat")
