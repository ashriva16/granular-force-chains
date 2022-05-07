import os

with open('py_pot2.py', 'r') as file:
    # read a list of lines into data
    line = file.readlines()

disorder = 0

line[17] = '        disorder = '+str(disorder)+'\n'
with open('py_pot2.py', 'w') as file:
    file.writelines(line)

print(line[17])

var_cmd = " -var D "+str(disorder)
# cmd = "mpirun -np 4 lmp_mpi -in potential.lmp"+var_cmd
cmd = "mpirun -np 4 lmp_mpi -in gran_energy_check.lmp"+var_cmd
print(cmd)
returned_value = os.system(cmd)
