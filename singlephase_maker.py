import sys
import os
import numpy
import math
import random
import matplotlib.pyplot as plt
import decimal

print()
phase = sys.argv[1]
grains = int(sys.argv[2])
print(f'--- Constructing {phase} single phase model with {grains} grains ---')

resolution = 64
voxels = resolution*resolution*resolution


# ------- CRSS based on Chi-Toan et al. ------- #
CRSS_prismatic = 100.0e6
CRSS_basal = 100.0e6
CRSS_pyramidal = 700.0e6

CRSS_BCC = 33.3e6 
# --------------------------------------------- #

material_yaml = open("material.yaml", 'w')
material_yaml.write("homogenization:\n\
  SX:\n\
    mech: {type: none}\n\
microstructure:\n")

oris = open("quaternions.txt",'w')
oris.write("w x y z\n")
# ------------ Generate random quaternion orientations ----------- #
resetcount = 0
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
        
        # random quaternion generation equations from : http://planning.cs.uiuc.edu/node198.html
        
        # normalize quaternion:
        if math.sqrt(w[n]**2 + x[n]**2 + y[n]**2 + z[n]**2) == float(1.000000000000000000000000000000000000) and -0.99<w[n]<0.99 and -0.99<x[n]<0.99 and -0.99<y[n]<0.99 and -0.99<z[n]<0.99:
            orientations = numpy.array([w[n], x[n], y[n], z[n]])
            quaternions[n] = orientations
            print("grain:",n+1, "sqrt(sum(quats)) = :", math.sqrt(w[n]**2 + x[n]**2 + y[n]**2 + z[n]**2))
            print("quat:", quaternions[n])
            i += 1
            
            material_yaml.write("# grain %d \n\
- constituents:\n\
  - fraction: 1.0\n\
    orientation: [%.18f, %.18f, %.18f, %.18f]\n\
    phase: Titanium_%s\n\
  homogenization: SX\n"\
%(n+1,quaternions[n,0], quaternions[n,1], quaternions[n,2], quaternions[n,3], phase))

            oris.write("%.18f %.18f %.18f %.18f\n"\
%(quaternions[n,0], quaternions[n,1], quaternions[n,2], quaternions[n,3]))

        else:  
            i = 0 #restart the loop
            #print("grain:",n+1, "invalid quaternion. reset loop...")
            resetcount += 1


# ---------------------------------------------------------------- #

if phase == "alpha":
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
      xi_0_sl: [%d, %d, 0.0, %d]\n\
      xi_inf_sl: [%d, %d, 0.0, %d]\n\
      type: phenopowerlaw\n"\
%(CRSS_prismatic, CRSS_basal, CRSS_pyramidal, CRSS_prismatic, CRSS_basal, CRSS_pyramidal))

if phase == "beta":
    material_yaml.write("phase:\n\
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
      xi_0_sl: [%d, %d]\n\
      xi_inf_sl: [%d, %d]\n\
      type: phenopowerlaw\n"\
%(CRSS_BCC, CRSS_BCC, CRSS_BCC, CRSS_BCC))

material_yaml.close()
print("--- Created ", phase, "model. ---")
print("normalized quaternions with a precsion of 12 d.p. per component.\n\
re-normalized", resetcount, "quaternions until correct.")
print("made: material.yaml with       ", quaternions.shape[0], "grains.")

# --------------------------------------- #

#print("quaternions array:", quaternions)
oris.close()

print("      quaternions.txt with     ", quaternions.shape[0], "quaternions.")
print("output files placed into cwd.")
print("use seeds_fromRandom, geom_fromVoronoiTesselation to create corresponding .geom")
# ---------------------------------------- #
