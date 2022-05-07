import h5py
import os
import shutil
import numpy as np
import lammps_logfile
import matplotlib.pyplot as plt
import pickle

def plot(Y,ylabel,label):
    plt.plot(Y, label=label)
    plt.xlabel('Step')
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(label+".png", format='png',dpi=1000,bbox_inches='tight')
    plt.close()


def shutil_copy(file, dst):
    # if os.path.exists(dst):
    #     shutil.rmtree(dst)
    try:
        shutil.copyfile(file, dst+file)
        print("File copied successfully.")
    except:
        print("Error occurred while copying "+file)

def shutil_move(file, dst):
    # if os.path.exists(dst):
    #     shutil.rmtree(dst)
    try:
        shutil.move(file, dst+file)
        print("File copied successfully.")
    except:
        print("Error occurred while moving "+file)


def editpypot(disorder, att_coeff, rep_coeff1, rep_coeff2):
    with open('py_pot2.py', 'r') as file:
        # read a list of lines into data
        line = file.readlines()

    line[71] = '        disorder = '+str(disorder)+'\n'
    line[72] = '        att_coeff = '+str(att_coeff)+'\n'
    line[73] = '        rep_coeff1 = '+str(rep_coeff1)+'\n'
    line[74] = '        rep_coeff2 = '+str(rep_coeff2)+'\n'
    with open('py_pot2.py', 'w') as file:
        file.writelines(line)

def store_data_in_hdffile_2(name_, data, hf):
    #     if (name_ not in hf):
    #         hf.create_dataset(name_, (np.append(1, data.shape)),
    #                           'float64')
    n, p = data.shape
    hf[name_] = data.reshape(1, n, p)

def read_pairwise(STEP):

    file1 = open('neigh.deform', 'r')
    lines = [line.rstrip() for line in file1] 

    for num, line in enumerate(lines, 0):
        if 'ITEM' in line:
            if 'ITEM: TIMESTEP' in line:
                curr_timestep = int(lines[num+1])
            if 'ITEM: NUMBER OF ENTRIES' in line:
                num_entries = int(lines[num+1])
                table = []
            if 'ITEM: ENTRIES' in line:
                if(curr_timestep != STEP):
                    continue
                x = lines[num+1:num+1+num_entries]
                for i in range(len(x)):
                    table.append(list(map(float, x[i].split(" "))))

                if(curr_timestep == STEP):
                    break

    print("read_pairwise:", curr_timestep, "STEP:", STEP, flush=True)
    pairwise = np.array(table.copy())

    return pairwise

def read_grain_data(STEP):
    file2 = open('visualize.deform', 'r')
    lines = [line.rstrip() for line in file2]

    for num, line in enumerate(lines, 0):
        if 'ITEM' in line:
            if 'ITEM: TIMESTEP' in line:
                curr_timestep = int(lines[num+1])
            if 'ITEM: NUMBER OF ATOMS' in line:
                num_entries = int(lines[num+1])
                table = []
            if 'ITEM: ATOMS' in line:
                if(curr_timestep != STEP):
                    continue
                x = lines[num+1:num+1+num_entries]
                for i in range(len(x)):
                    table.append(list(map(float, x[i].split(" "))))

                if(curr_timestep == STEP):
                    break

    print("read_grain_data:", curr_timestep, "STEP:", STEP, flush=True)

    grains = np.array(table.copy())
    grains = grains[grains[:, 1] != 3]
    grains = grains[grains[:, 0].argsort()]

    return grains

#------------------------------------------------------------------------------------
disorder = 0
att_coeff = 0
rep_coeff1 = 1
rep_coeff2 = 1

rep_diorder = rep_coeff1/rep_coeff2

# script = "Hooke_deform_final_cohesion"
script = "Hooke_deform_final"
# script = "Hooke_H-wall_V-periodic_Particle-indent"

with open('disorder.pckl', 'wb') as f:
    pickle.dump([script, 
                disorder,
                att_coeff,
                rep_diorder,
                rep_coeff1,
                rep_coeff2], f)

file = "../results/"
file += script+"/"
file += str(disorder) + "_"
file += str(att_coeff) + "_"
file += str(rep_diorder) + "/"
print(file)

if not os.path.exists(file):
    os.makedirs(file)

editpypot(disorder, att_coeff, rep_coeff1, rep_coeff2)

var_cmd = " -var D "+str(disorder)
cmd = "mpirun -np 4 lmp_mpi -in "+script+'.lmp'+var_cmd
print(cmd, flush=True)
returned_value = os.system(cmd)
# print(returned_value, flush=True)

log1 = lammps_logfile.File("log.lammps")
log2 = lammps_logfile.File("log.deform")
KE = np.append(log1.data_dict['KinEng'], log2.data_dict['KinEng'])
PE = np.append(log1.data_dict['PotEng'], log2.data_dict['PotEng'])
TE = np.append(log1.data_dict['TotEng'], log2.data_dict['TotEng'])
pressure = np.append(log1.data_dict['v_p2'], log2.data_dict['v_p2'])

plot(KE, 'Energy Units', 'KinEng')
plot(PE, 'Energy', 'PotEng')
plot(TE, 'Energy', 'TotEng')
plot(pressure, 'Pressure', 'Pressure')

STEP = int(log2.data_dict['Step'][-1])
with h5py.File('data.h5', 'a') as f:
    tmp = read_grain_data(0)
    store_data_in_hdffile_2('grains0', tmp, f)
    tmp = read_grain_data(STEP)
    store_data_in_hdffile_2('grains', tmp, f)
    tmp = read_pairwise(0)
    store_data_in_hdffile_2('pairwise0', tmp, f)
    tmp = read_pairwise(STEP)
    store_data_in_hdffile_2('pairwise', tmp, f)

shutil_copy(script+".lmp", file)
shutil_copy("py_pot2.py", file)
shutil_move("hooke11.table", file)
shutil_move("hooke22.table", file)
shutil_move("hooke33.table", file)
shutil_move("hooke12.table", file)
shutil_move("hooke13.table", file)
shutil_move("hooke23.table", file)
shutil_move("data.h5", file)
shutil_move("neigh.deform", file)
shutil_move("visualize.deform", file)
shutil_move("visualize.equil", file)
shutil_move("log.lammps", file)
shutil_move("log.deform", file)
shutil_move("KinEng.png", file)
shutil_move("PotEng.png", file)
shutil_move("TotEng.png", file)
shutil_move("Pressure.png", file)
shutil_move("disorder.pckl", file)
