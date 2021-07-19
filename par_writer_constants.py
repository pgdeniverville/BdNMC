import math

mb=1e-31
m_proton=0.938
mpi0=0.134
meta=0.547862
mrho=0.77549
momega=0.782
mphi=1.020
pi=math.pi

meson_per_pi0_miniboone = {'pi0_decay' : '1.0', 'eta_decay' : str(1.0/30.0), 'rho_decay' : str(1.0/20.0), 'omega_decay' : '0.046', 'phi_decay' : str(1.0/150.0), 'pi0_decay_baryonic' : '1.0', 'eta_decay_baryonic' : str(1.0/30.0), 'rho_decay_baryonic' : str(1.0/20.0), 'omega_decay_baryonic' : '0.046', 'phi_decay_baryonic' : str(1.0/150.0)}

meson_per_pi0_lanl = {'pi0_decay' : 1.0, 'piminus_capture' : '0.63'}

meson_per_pi0_lsnd = {'pi0_decay' : '1.0'}

meson_per_pi0_coherent = {'pi0_decay' : 1.0, 'piminus_capture' : '0.63'}

#meson_per_pi0_numi = {'pi0_decay' : '1.0', 'eta_decay' : str(0.078), 'rho_decay' : str(0.11), 'omega_decay' : '0.11', 'phi_decay' : str(0.02)}

meson_per_pi0_numi = {'pi0_decay' : '1.0', 'eta_decay' : '0.114', 'rho_decay' : str(0.11), 'omega_decay' : '0.11', 'phi_decay' : str(0.02)}

meson_per_pi0_ship = {'pi0_decay' : '1.0', 'eta_decay' : str(0.105), 'rho_decay' : str(0.11), 'omega_decay' : '0.11', 'phi_decay' : str(0.02)}

##################
#Material Presets#
##################

Empty_string = "material Vacuum\nnumber_density 0\nproton_number 0\nneutron_number 0\nelectron_number 0\nmass 0.945778\n"

Hydrogen_string = "material Hydrogen\nnumber_density 7.26942e22\nproton_number 1\nneutron_number 0\nelectron_number 1\nmass 0.945778\n"

Water_string = "material Oxygen\nnumber_density 3.34184e22\nproton_number 8\nneutron_number 8\nelectron_number 8\nmass 15.999\nmaterial Hydrogen\nnumber_density 6.68368e22\nproton_number 1\nneutron_number 0\nelectron_number 1\nmass 0.945778\n"

Carbon_string = "material Carbon\nnumber_density 3.63471e22\nproton_number 6\nneutron_number 6\nelectron_number 6\nmass 11.2593\n"

Argon_string = "material Argon\nnumber_density 2.11e22\nproton_number 18\nneutron_number 22\nelectron_number 18\nmass {0}\n".format(str(39.948*0.938))

#Density = 1400 kg/m^3
Argon_LAr_string = "material Argon\nnumber_density 2.23725e22\nproton_number 18\nneutron_number 22\nelectron_number 18\nmass {0}\n".format(str(39.948*0.938))

ND280_string = "material nd280stuff\nnumber_density 3.7e23\nproton_number 1\nneutron_number 1\nelectron_number 1\nmass 0.945778\n"
#1.78 g/cm^3
Argon_STP_string = "material Argon\nnumber_density 2.68e22\nproton_number 18\nneutron_number 22\nelectron_number 18\nmass {0}\n".format(str(39.948*0.938))

Argon_string_SBND = "material Argon\nnumber_density 2.24e22\nproton_number 18\nneutron_number 22\nelectron_number 18\nmass {0}\n".format(str(39.948*0.938))
#Temporary!
Steel_string = Carbon_string

#Satisfies 3.72 g/cm^3
SHiP_string = "material ship_stuff\nnumber_density 1.11e24\nproton_number 1\nneutron_number 1\nelectron_number 1\nmass 1.89\n"

ND280_string = "material nd280stuff\nnumber_density 3.7e23\nproton_number 1\nneutron_number 1\nelectron_number 1\nmass 1.89\n"

Sodium_Iodide_string = "material Sodium\nnumber_density 1.58e22\nproton_number 11\n neutron_number 23\nelectron_number 11\nmass 21.61\nmaterial Iodine\nnumber_density 1.58e22\nproton_number 53\nneutron_number 72\nelectron_number 53\nmass 119.03\n"

Cesium_Iodide_string = "material Cesium\nnumber_density 1.04e22\nproton_number 55\n neutron_number 78\nelectron_number 55\nmass 132.9\nmaterial Iodine\nnumber_density 1.04e22\nproton_number 53\nneutron_number 72\nelectron_number 53\nmass 119.03\n"

#Don't know what the atomic makeup of the MINOS detector is. Not using this for event generation, so it should be okay.
MINOS_string = "material Steel\nnumber_density 5e24\nproton_number 1\nneutron_number 1\nelectron_number 1\nmass 0\n"

NOvA_string = "material Liquid_Scintillator\nnumber_density 5.16e22\nproton_number 8\nneutron_number 6\nelectron_number 8\nmass 14.011"

High_pressure_argon_string="material Argon\nnumber_density 1.48e20\nproton_number 18\nproton_number 18\nneutron_number 22\nelectron_number 18\nmass {0}\n".format(str(39.948*0.938))

BEBC_Hydrogen = "material Hydrogen\nnumber_density 7.07e21\nproton_number 1\nneutron_number 1\nelectron_number 1\nmass 0.938"

BEBC_Neon = "material Neon\nnumber_density 2.01e22\nproton_number 10\nneutron_number 10\nelectron_number 10\nmass 18.9286"