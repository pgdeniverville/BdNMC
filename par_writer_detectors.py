from par_writer_constants import *

##################
#DETECTOR PRESETS#
##################
def Cylinder(f,xpos,ypos,zpos,radius,length,theta,phi):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))

def Cuboid(f,xpos,ypos,zpos,width,length,height,phi,theta,psi):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))

def miniboone_detector(f,xpos=0.0,ypos=-1.9,zpos=491.0,radius=5.0):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def miniboone_detector_numi(f,xpos=0.0,ypos=0.0,zpos=100.0,radius=5.0):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def miniboone_detector_full(f,xpos=0.0,ypos=-1.9,zpos=491.0,radius=6.106):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def dune_mpd_detector(f,xpos=0.0,ypos=0.0,zpos=581.0,radius=2.5, length=5, theta=0, phi=math.pi/2.0):
    f.write("\ndetector cylinder\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(High_pressure_argon_string)

def dune_detector(f,xpos=0.0, ypos=0.0, zpos=573, width=3, height =2, length =5, theta=0, phi=0,psi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    f.write(Argon_LAr_string)

def dune_detector_off_axis(f,xpos=30, ypos=0.0, zpos=573, width=3, height =2, length =5, theta=0, phi=0,psi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    f.write(Argon_LAr_string)

def test_sphere(f,xpos=0.0,ypos=0.0,zpos=0.0,radius=1.0):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def test_cylinder(f,xpos=0.0,ypos=0.0,zpos=0.0,radius=1.0, length=1.0, theta=0, phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def test_cuboid(f,xpos=0.0,ypos=0.0,zpos=0.0,height=1.0,length=1.0,width=1.0,theta=0.0,phi=0,psi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

def bebc_detector(f,xpos=0.0,ypos=0.0,zpos=405,height=3.49,radius=1.23,theta=math.pi/2.0,phi=0):
    Cylinder(f,xpos,ypos,zpos,height,radius,theta,phi)
    f.write('\n')
    f.write(BEBC_Hydrogen)
    f.write('\n')
    f.write(BEBC_Neon)

MINOS_absorber_z=270
MINOS_target_z=950
#Don't use this for actual event generation!
def MINOS_detector(f,xpos=0.0,ypos=0.0,zpos=MINOS_target_z,radius=2.2,length=1.7,theta=0,phi=0):
    print("This detector should not be used for event generation!")
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(MINOS_string)

def MINOS_absorber_detector(f):
    MINOS_detector(f,zpos=MINOS_absorber_z)

def ICARUS_detector(f,xpos=0,ypos=0,zpos=600,height=3.16,length=17.95,width=6,theta=0,phi=0,psi=0):
    f.write("\ndetector cuboid\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    f.write(Argon_LAr_string)

ICARUS_NuMI_angle=0.1
ICARUS_NuMI_distance=789
def ICARUS_detector_NuMI(f,xpos=ICARUS_NuMI_distance*math.sin(ICARUS_NuMI_angle),ypos=0,zpos=ICARUS_NuMI_distance*math.cos(ICARUS_NuMI_angle),theta=ICARUS_NuMI_angle):
    ICARUS_detector(f,xpos=xpos,ypos=ypos,zpos=zpos,theta=theta)

CHARM2_string="material Carbon\nnumber_density 7.00652e22\nproton_number 6\nneutron_number 6\nelectron_number 6\nmass 11.2593\n"

#Assuming this is made of carbon.
def CHARM2_detector(f,xpos=0.0,ypos=0,zpos=870.8,width=3.7,length=36,height=3.7,phi=0,theta=0,psi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5} \ndet-theta {6}\ndet-phi {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    f.write(CHARM2_string)

#Using the CHARM2_string, don't use for scattering
def CHARM_decay_detector(f,xpos=5,ypos=0,zpos=480,width=3,length=35,height=3,phi=0,theta=0,psi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5} \ndet-theta {6}\ndet-phi {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    f.write(CHARM2_string)

#NOvA_absorber_d=240
#NOvA_target_d=920
NOvA_absorber_d=447
NOvA_target_d=990
#NOvA_angle=0.0575959#3.3 degree
NOvA_Target_Angle=0.0122
NOvA_Absorber_Angle=0.024

def NOvA_detector(f,xpos=0.0,ypos=NOvA_target_d*math.sin(NOvA_Target_Angle),zpos=NOvA_target_d*math.cos(NOvA_Target_Angle),height=4.2,length=14.3,width=2.9,theta=-NOvA_Target_Angle,phi=0,psi=0):
    print("This NOvA detector is prelimary!")
    #print("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5}\ndet-phi {6}\ndet-theta {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    f.write(NOvA_string)

def NOvA_absorber_detector(f):
    NOvA_detector(f,xpos=0.0,ypos=NOvA_absorber_d*math.sin(NOvA_Absorber_Angle),zpos=NOvA_absorber_d*math.cos(NOvA_Absorber_Angle))

def SBND_detector(f,xpos=0.457,ypos=0,zpos=110,width=4,length=5,height=4,phi=0,theta=0,psi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5} \ndet-theta {6}\ndet-phi {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    #Need to figure out what it's made of!
    f.write(Argon_LAr_string)

def SBND_dump_detector(f,xpos=0.457,ypos=0,zpos=60,width=4,length=5,height=4,phi=0,theta=0,psi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5} \ndet-theta {6}\ndet-phi {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    #Need to figure out what it's made of!
    f.write(Argon_LAr_string)

SHIP_Side=1.92856

def ship_detector(f,xpos=0.0,ypos=0,zpos=56,width=1.87,length=2,height=0.69,phi=0,theta=0,psi=0):
    #def ship_detector(f,xpos=0.0,ypos=0,zpos=30.0,radius=0.655,length=2.645,theta=0,phi=0):
    f.write("\ndetector cuboid\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nwidth {3}\nlength {4}\nheight {5} \ndet-theta {6}\ndet-phi {7}\ndet-psi {8}".format(str(xpos),str(ypos),str(zpos),str(width),str(length),str(height),str(phi),str(theta),str(psi)))
    f.write('\n')
    #Need to figure out what it's made of!
    f.write(SHiP_string)

#temp detector until I implement proper geometry handling
#This is actually the P0D, with only the proper number of neutrons and protons. NO ATOMS IMPLEMENTED
#Double check the fiducial mass on the pod. Is it 3 tons of water? 13 tons of stuff?
#Need to update with cuboid shape!
def t2k_ND280(f):
    xpos=11;ypos=0;zpos=280;detphi=0;radius=0.9413;dettheta=0.0436332;length=1.7;
    f.write("\ndetector cylinder\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\ndet-theta {4}\ndet-phi {5}\nlength {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(dettheta),str(detphi),str(length),str(length)))
    f.write('\n')
    f.write(ND280_string)

#This is 2 tons of water in an FGD, really they only have 1.1 tons of material each.
def t2k_FGD(f):
    print("This detector doesn't have accurate geometry and is only 2 tons of material. Improve for physics running.")
    Cuboid(f,xpos=11,ypos=0,zpos=280,width=2.3,length=0.35,height=2.4,phi=0,theta=0.0436332,psi=0)
    f.write('\n')
    f.write(Water_string)

def t2k_superK1000(f):
    xpos=12867.7;ypos=0;zpos=294719;detphi=0;radius=190.5;dettheta=1.5708;length=410;
    f.write("\ndetector cylinder\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\ndet-theta {4}\ndet-phi {5}\nlength {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(dettheta),str(detphi),str(length),str(length)))
    print("WARNING: This detector is enlarged 1000 fold. Its efficiency must be suppressed in a corresponding manner.")
    f.write('\n')
    f.write(Water_string)

def lsnd_detector(f,xpos=0.0,ypos=-4.65,zpos=29.8,radius=5.7/2.0-0.35,length=8.3,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Carbon_string)
    f.write('\n')
    f.write(Hydrogen_string)

def coherent_detector_LAr(f,xpos=20.0,ypos=0.0,zpos=0.0,radius=0.48,length=0.96,theta=pi/2.0,phi=pi/2.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_STP_string)


def coherent_detector_LAr_29kg(f,xpos=19.3,ypos=0.0,zpos=-20.83,radius=0.124,length=0.425,theta=pi/2.0,phi=pi/2.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_LAr_string)

def lsnd_detector_test(f,xpos=0.0,ypos=-4.65,zpos=29.8,radius=5.7/2.0-0.35,length=8.3,theta=0,phi=0):
    f.write("\ndetector sphere\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(xpos),str(ypos),str(zpos),str(radius)))
    f.write('\n')
    f.write(Carbon_string)
    f.write('\n')
    f.write(Hydrogen_string)

#def coherent_detector_LAr(f,xpos=16.8,ypos=0.0,zpos=-23.6,radius=0.44,length=0.88,theta=pi/2.0,phi=pi/2.0):

def coherent_detector_NaI(f,xpos=20.0,ypos=0.0,zpos=0.0,radius=0.601,length=2*0.601,theta=pi/2.0,phi=pi/2.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Sodium_Iodide_string)
0.08
def coherent_detector_CsI(f,xpos=19.3,ypos=0.0,zpos=0.0,radius=0.08,length=2*0.08,theta=pi/2.0,phi=pi/2.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Cesium_Iodide_string)

def coherent_detector_CsI_1T(f,xpos=20,ypos=0.0,zpos=0.0,radius=0.328,length=2*0.328,theta=pi/2.0,phi=pi/2.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Cesium_Iodide_string)

def captain_detector(f,xpos=0.0,ypos=0.0,zpos=30.0,radius=1,length=1.15,theta=pi/2.0,phi=0.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

def captain_detector_off(f,xpos=30.0,ypos=0.0,zpos=0.0,radius=1,length=1.15,theta=pi/2.0,phi=0.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

#4500 kg -> 5 tons
def lanl_detector(f,xpos=9.50568,ypos=19.07841,zpos=-0.3589,radius=0.909967,length=1.227,theta=0.0,phi=0.0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_LAr_string)

def NA62_decay_vol(f,xpos=0.0,ypos=0.0,zpos=149.5,radius=1,length=135,theta=0,phi=0):
    Cylinder(f,xpos,ypos,zpos,radius,length,theta,phi)
    f.write('\n')
    f.write(Empty_string)

def nucal(f,xpos=0.0,ypos=0.0,zpos=75.5,radius=1.3,length=23,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
   #Need to replace this with Aluminum, Iron, Liquid Scintilator if using for scattering
    f.write('\n')
    f.write(Carbon_string)

pip2_distance=18

def pip2_coherent_det(f,pip2_angle=0):
    print("pip2_angle={}".format(pip2_angle))
    Cylinder(f,math.sin(pip2_angle)*pip2_distance,0,math.cos(pip2_angle)*pip2_distance,1.75,10,pip2_angle,0)
    f.write('\n')
    f.write(Argon_LAr_string)

def pip2_30(f,pip2_angle=30.0/180.0*math.pi):
    pip2_coherent_det(f,pip2_angle=pip2_angle)

def pip2_90(f,pip2_angle=90.0/180.0*math.pi):
    pip2_coherent_det(f,pip2_angle=pip2_angle)


def ship_detector_modular(f,radius=0.8268,length=3.34,theta=0,phi=0):
    #def ship_detector(f,xpos=0.0,ypos=0,zpos=30.0,radius=0.655,length=2.645,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(_DET_XPOS),str(_DET_YPOS),str(_DET_ZPOS),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

def OPERA_detector(f,xpos=0.0,ypos=0.0,zpos=7.3e6,radius=56.4,length=1,theta=0,phi=0):
    print("Don't use this for event generation!")
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(xpos),str(ypos),str(zpos),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(Argon_string)

def MINOS_detector_modular(f,xpos=0.0,radius=2,length=1.7,theta=0,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(_DET_XPOS),str(_DET_YPOS),str(_DET_ZPOS),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(MINOS_string)
'''
def NOvA_detector_modular(f,xpos=0.0,radius=2,length=14,theta=-NOvA_angle,phi=0):
    f.write("\ndetector cylinder\n");
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\nlength {4}\ndet-theta {5}\ndet-phi {6}\n".format(str(_DET_XPOS),str(_DET_YPOS),str(_DET_ZPOS),str(radius),str(length),str(theta),str(phi)))
    f.write('\n')
    f.write(MINOS_string)
'''
def miniboone_detector_modular(f,radius=5.0):
    f.write("\ndetector sphere\n")
    f.write("x-position {0}\ny-position {1}\nz-position {2}\nradius {3}\n".format(str(_DET_XPOS),str(_DET_YPOS),str(_DET_ZPOS),str(radius)))
    f.write('\n')
    f.write(Hydrogen_string)
    f.write('\n')
    f.write(Carbon_string)

argoneut_Angle = 0.0523599; #3 Degrees down.
argoneut_z = 1033;
argoneut_y = -61;
argoneut_x = 0;
def argoneut_detector(f):
    Cuboid(f,xpos=argoneut_x,ypos=argoneut_y,zpos=argoneut_z,width=0.47,height=0.4,length=0.9,phi=argoneut_Angle,theta=0,psi=0)
    f.write('\n')
    f.write(Argon_LAr_string)