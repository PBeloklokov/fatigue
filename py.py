import os
import numpy as np

# dir_name = "fatigue"
# file_names = os.listdir(dir_name)

# for file_name in file_names:
#     # file_number = int(file_name[7:-4])
#     # print(f"fatigue\output_{file_number:06d}")
#     # os.rename("fatigue/" + file_name, f"fatigue/output_{file_number:06d}.vtk")
#     with open(dir_name + "/" + file_name, "r"):

sigma_u = 340e6
sigma_v = 1160e6
gamma = 0.5
beta_L = 0.31

attributes_size = 1001
psi = np.linspace(0, 1, attributes_size)


def get_dN(file_name):
    file_lines = open(file_name, "r").read().split("\n")

    it = 1
    while file_lines[it] != "DATASET UNSTRUCTURED_GRID":
        it += 1
    it += 1
    line_split = file_lines[it].split(" ")
    points_num = int(line_split[1])

    it += 1
    while file_lines[it][:5] != "CELLS":
        it += 1
    line_split = file_lines[it].split(" ")
    cells_num = int(line_split[1])

    it += 1
    while file_lines[it].split(" ")[0] != "CELL_DATA":
        it += 1
    it += 3
    att_first_num = it    

    sigma_xx, sigma_xy, sigma_yy = np.zeros(cells_num), np.zeros(cells_num), np.zeros(cells_num)

    it += 1
    while file_lines[it] != "SCALARS sigma_xx double 1":
        it += 1
    it += 2
    for sigma_it in range(cells_num):
        sigma_xx[sigma_it] = (float(file_lines[it + 3*sigma_it]) + float(file_lines[it + 3*sigma_it + 1]) + float(file_lines[it + 3*sigma_it + 2])) / 3

    it += points_num + 2
    for sigma_it in range(cells_num):
        sigma_xy[sigma_it] = (float(file_lines[it + 3*sigma_it]) + float(file_lines[it + 3*sigma_it + 1]) + float(file_lines[it + 3*sigma_it + 2])) / 3

    it += points_num + 2
    for sigma_it in range(cells_num):
        sigma_yy[sigma_it] = (float(file_lines[it + 3*sigma_it]) + float(file_lines[it + 3*sigma_it + 1]) + float(file_lines[it + 3*sigma_it + 2])) / 3

    delta_N_n = 1e10
    delta_sigma_1 = np.sqrt((sigma_xx - sigma_yy)**2 + 4 * sigma_xy**2)
    sigma_1 = (sigma_xx + sigma_yy) / 2 + delta_sigma_1 / 2
    sigma_eq = np.sqrt(sigma_1 * delta_sigma_1 / 2)
    
    with open("sigma_py.vtk", "a") as file:
        for en in range(cells_num):
            file.write(f"en = {en} xx = {sigma_xx[en]} yy = {sigma_yy[en]} xy = {sigma_xy[en]} _1 = {sigma_1[en]} eq = {sigma_eq[en]}\n")

    B_n = np.zeros(cells_num)
    for en in range(cells_num):
        if (sigma_1[en] > 0 and sigma_eq[en] > sigma_u):
            B_n[en] = 1e-3 * pow((sigma_eq[en] - sigma_u) / (sigma_v - sigma_u), 1 / beta_L) / (2 * (1 - gamma))
            attribute_en = int(file_lines[att_first_num + en])
            psi_k_n = psi[attribute_en - 1]
            delta_N_n_en = 0.5 * ((1 - pow(psi_k_n, 1 - gamma)) / (1 - gamma) - (1 - pow(psi_k_n, 2 * (1 - gamma))) / (2 * (1 - gamma))) / B_n[en]

            if delta_N_n_en < delta_N_n: 
                delta_N_n = delta_N_n_en
            with open("deformation_py.vtk", "a") as file:
                file.write(f"en = {en} attr = {attribute_en} psi_k_n = {psi_k_n} B_n = {B_n[en]} dN_n = {delta_N_n}\n")

    return delta_N_n




file = open("deformation_py.vtk", "r+")
file.truncate(0)
file = open("sigma_py.vtk", "r+")
file.truncate(0)
dN = get_dN("fatigue/output_000000.vtk")
print(dN)