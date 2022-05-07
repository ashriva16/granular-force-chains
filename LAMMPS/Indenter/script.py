import os
import shutil
import numpy as np


# def editpypot(i):
#     with open('py_pot2.py', 'r') as file:
#         # read a list of lines into data
#         line = file.readlines()

#     line[71] = '        disorder = '+str(i)+'\n'
#     with open('py_pot2.py', 'w') as file:
#         file.writelines(line)

# Running LAMMMPS Simulation ------------------------------------------------------------------
# fraction = .5
# disorder = 0
# rep_expnenet = hard  # M
# att_exponent = 6  # N
# var_cmd = "-var F "+str(fraction)
# var_cmd = var_cmd+" -var D "+str(disorder)
# var_cmd = var_cmd+" -var M "+str(rep_expnenet)
# var_cmd = var_cmd+" -var N "+str(att_exponent)

# file = "HARDNESS/10000"+"_"+str(disorder) + "_" + str(att_exponent) + "/"
# file = file + str(fraction)+"_"+str(disorder)+"_" + \
#     str(rep_expnenet)+"_"+str(att_exponent)+"/"
# if not os.path.exists(file):
#     os.makedirs(file)

# cmd = "mpirun -np 4 lmp < in.indent_displacement_controlled.lmp "+var_cmd
# print(cmd)
# returned_value = os.system(cmd)  # returns the exit code in unix

# shutil.move("dump.neigh", file+"dump.neigh")
# shutil.move("log.lammps", file+"log.lammps")
# shutil.move("visualize.du", file+"visualize.du")

#------------------------------------------------------------------------------------
with open('py_pot2.py', 'r') as file:
    # read a list of lines into data
    line = file.readlines()

disorder = 0.1

line[71] = '        disorder = '+str(disorder)+'\n'
with open('py_pot2.py', 'w') as file:
    file.writelines(line)

print(line[71])

var_cmd = " -var D "+str(disorder)
# cmd = "mpirun -np 4 lmp_mpi -in Hooke_deform_final_cohesion.lmp"+var_cmd
# cmd = "mpirun -np 4 lmp_mpi -in Hooke_deform_final.lmp"+var_cmd
cmd = "mpirun -np 4 lmp_mpi -in Hooke_H-wall_V-periodic_Particle-indent.lmp"+var_cmd
print(cmd)
returned_value = os.system(cmd)
