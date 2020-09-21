import sys
import os
import numpy
import math
import random
import matplotlib.pyplot as plt
import decimal

# USER DEFINED #
print()
frac_beta = int(sys.argv[1])
grains = int(sys.argv[2])

print(f'--- Constructing {frac_beta} % beta,', 100-frac_beta, '% alpha mixed phase model ---')
frac_beta = frac_beta/100
frac_alpha = 1 - frac_beta

print(f' {grains} Total grains... ')
################

resolution = 32
voxels = resolution*resolution*resolution

phase = numpy.random.choice([1,2],p = [frac_alpha,frac_beta], size = voxels)

# ------- CRSS based on Chi-Toan et al. ------- #
CRSS_prismatic = 100.0e6
CRSS_basal = 100.0e6
CRSS_pyramidal = 700.0e6

CRSS_BCC = 33.3e6 
# --------------------------------------------- #

material_yaml = open("material.yaml", 'w')
material_yaml.write("# mixed phase %i percent beta %i percent alpha\n\
homogenization:\n\
  SX:\n\
    mech: {type: none}\n\
microstructure:\n"\
%(frac_beta*100, frac_alpha*100))

# ------------ Generate random quaternion orientations ----------- #

quaternions = numpy.ndarray((grains,4))
w = numpy.ndarray((grains))
x = numpy.ndarray((grains))
y = numpy.ndarray((grains))
z = numpy.ndarray((grains))
for n in range(0, grains): # (inclusive,exclusive)
    i = 0
    while i == 0:
        u1 = round(random.randint(1,999999999999999999)/1e18, 18)
        u2 = round(random.randint(1,999999999999999999)/1e18, 18)
        u3 = round(random.randint(1,999999999999999999)/1e18, 18)
    
        w[n] = round(decimal.Decimal(math.sqrt(1-u1)*math.sin(2*math.pi*u2)),18)
        x[n] = round(decimal.Decimal(math.sqrt(1-u1)*math.cos(2*math.pi*u2)),18)
        y[n] = round(decimal.Decimal(math.sqrt(u1)*math.sin(2*math.pi*u3)),18)
        z[n] = round(decimal.Decimal(math.sqrt(u1)*math.cos(2*math.pi*u3)),18)
        
        # found from : http://planning.cs.uiuc.edu/node198.html
        
        # normalize quaternion:
        if math.sqrt(w[n]**2 + x[n]**2 + y[n]**2 + z[n]**2) == float(1.000000000000000000000000000000000000) and -0.99<w[n]<0.99 and -0.99<x[n]<0.99 and -0.99<y[n]<0.99 and -0.99<z[n]<0.99:
            orientations = numpy.array([w[n], x[n], y[n], z[n]])
            quaternions[n] = orientations
            print("grain:",n+1, "sqrt(sum(quats)) = :", math.sqrt(w[n]**2 + x[n]**2 + y[n]**2 + z[n]**2))
            print("quat:", quaternions[n])
            i += 1
    
    if phase[n] == 1:
        material_yaml.write("# grain %d\n\
- constituents:\n\
  - fraction: 1.0\n\
    orientation: [%.18f, %.18f, %.18f, %.18f]\n\
    phase: Titanium_alpha\n\
  homogenization: SX\n"\
%(n+1,quaternions[n,0], quaternions[n,1], quaternions[n,2], quaternions[n,3]))

    if phase[n] == 2:
        material_yaml.write("# grain %d\n\
- constituents:\n\
  - fraction: 1.0\n\
    orientation: [%.18f, %.18f, %.18f, %.18f]\n\
    phase: Titanium_beta\n\
  homogenization: SX\n"\
%(n+1,quaternions[n,0], quaternions[n,1], quaternions[n,2], quaternions[n,3]))

# ---------------------------------------------------------------- #

material_yaml.write("phase:\n\
  Titanium_alpha:\n\
    elasticity: {C_11: 160.0e9, C_12: 90.0e9, C_13: 66.0e9, C_33: 181.7e9, C_44: 46.5e9, type: hooke}\n\
    generic:\n\
      output: [F, Fe, Fp, P, Lp, O]\n\
    lattice: hex\n\
    c/a: 1.587\n\
    plasticity:\n\
      N_sl: [3, 3, 0, 6]\n\
      a_sl: 2.0\n\
      atol_xi: 1.0\n\
      dot_gamma_0_sl: 0.001\n\
      h_0_sl_sl: 15e6\n\
      h_sl_sl: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]\n\
      n_sl: 10\n\
      output: [gamma_sl, xi_sl]\n\
      xi_0_sl: [%f, %f, 0.0, %f]\n\
      xi_inf_sl: [%f, %f, 0.0, %f]\n\
      type: phenopowerlaw\n\
  Titanium_beta:\n\
    elasticity: {C_11: 233.3e9, C_12: 135.5e9, C_44: 118.0e9, type: hooke}\n\
    generic:\n\
      output: [F, Fe, Fp, P, Lp, O]\n\
    lattice: bcc\n\
    plasticity:\n\
      N_sl: [12, 12]\n\
      a_sl: 2.0\n\
      atol_xi: 1.0\n\
      dot_gamma_0_sl: 0.001\n\
      h_0_sl_sl: 1000.0e6\n\
      h_sl_sl: [1, 1, 1.4, 1.4, 1.4, 1.4]\n\
      n_sl: 10\n\
      output: [gamma_sl, xi_sl]\n\
      xi_0_sl: [33300000.000000, 33300000.000000]\n\
      xi_inf_sl: [33300000.000000, 33300000.000000]\n\
      type: phenopowerlaw\n"\
%(CRSS_prismatic,CRSS_basal,CRSS_pyramidal, CRSS_prismatic,CRSS_basal,CRSS_pyramidal))
material_yaml.close()

print("normalized quaternions with a precsion of", 18 ,"d.p. per component.")
print("made material.yaml with       ", quaternions.shape[0], "grains")
      
# --------------------------------------- #
#print("quaternions array:\n",quaternions)
oris = open("quaternions.txt",'w')
oris.write("w x y z\n")
numpy.savetxt("quaternions.txt", quaternions, fmt='%.18f')
oris.close()
print("     quaternions.txt with     ", quaternions.shape[0], "quaternions")

# ---------------------------------------- #
