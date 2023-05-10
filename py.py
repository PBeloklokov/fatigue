import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# dir_name = "solid_circle_fixed"
# file_names = os.listdir(dir_name)

# for file_name in file_names:
#     file_number = int(file_name[7:-4])
#     if file_number % 20 != 0:
#         # print(file_number)
#         os.remove(dir_name + "/" + file_name)
    # print(f"fatigue\output_{file_number:06d}")
    # os.rename("fatigue/" + file_name, f"fatigue/output_{file_number:06d}.vtk")
    # with open(dir_name + "/" + file_name, "r"):

# sigma_u = 340e6
# sigma_v = 1160e6
# gamma = 0.5
# beta_L = 0.31

# attributes_size = 1001
# psi = np.linspace(0, 1, attributes_size)


# def get_dN(file_name):
#     file_lines = open(file_name, "r").read().split("\n")

#     it = 1
#     while file_lines[it] != "DATASET UNSTRUCTURED_GRID":
#         it += 1
#     it += 1
#     line_split = file_lines[it].split(" ")
#     points_num = int(line_split[1])

#     it += 1
#     while file_lines[it][:5] != "CELLS":
#         it += 1
#     line_split = file_lines[it].split(" ")
#     cells_num = int(line_split[1])

#     it += 1
#     while file_lines[it].split(" ")[0] != "CELL_DATA":
#         it += 1
#     it += 3
#     att_first_num = it    

#     sigma_xx, sigma_xy, sigma_yy = np.zeros(cells_num), np.zeros(cells_num), np.zeros(cells_num)

#     it += 1
#     while file_lines[it] != "SCALARS sigma_xx double 1":
#         it += 1
#     it += 2
#     for sigma_it in range(cells_num):
#         sigma_xx[sigma_it] = (float(file_lines[it + 3*sigma_it]) + float(file_lines[it + 3*sigma_it + 1]) + float(file_lines[it + 3*sigma_it + 2])) / 3

#     it += points_num + 2
#     for sigma_it in range(cells_num):
#         sigma_xy[sigma_it] = (float(file_lines[it + 3*sigma_it]) + float(file_lines[it + 3*sigma_it + 1]) + float(file_lines[it + 3*sigma_it + 2])) / 3

#     it += points_num + 2
#     for sigma_it in range(cells_num):
#         sigma_yy[sigma_it] = (float(file_lines[it + 3*sigma_it]) + float(file_lines[it + 3*sigma_it + 1]) + float(file_lines[it + 3*sigma_it + 2])) / 3

#     delta_N_n = 1e10
#     delta_sigma_1 = np.sqrt((sigma_xx - sigma_yy)**2 + 4 * sigma_xy**2)
#     sigma_1 = (sigma_xx + sigma_yy) / 2 + delta_sigma_1 / 2
#     sigma_eq = np.sqrt(sigma_1 * delta_sigma_1 / 2)
    
#     with open("sigma_py.vtk", "a") as file:
#         for en in range(cells_num):
#             file.write(f"en = {en} xx = {sigma_xx[en]} yy = {sigma_yy[en]} xy = {sigma_xy[en]} _1 = {sigma_1[en]} eq = {sigma_eq[en]}\n")

#     B_n = np.zeros(cells_num)
#     for en in range(cells_num):
#         if (sigma_1[en] > 0 and sigma_eq[en] > sigma_u):
#             B_n[en] = 1e-3 * pow((sigma_eq[en] - sigma_u) / (sigma_v - sigma_u), 1 / beta_L) / (2 * (1 - gamma))
#             attribute_en = int(file_lines[att_first_num + en])
#             psi_k_n = psi[attribute_en - 1]
#             delta_N_n_en = 0.5 * ((1 - pow(psi_k_n, 1 - gamma)) / (1 - gamma) - (1 - pow(psi_k_n, 2 * (1 - gamma))) / (2 * (1 - gamma))) / B_n[en]

#             if delta_N_n_en < delta_N_n: 
#                 delta_N_n = delta_N_n_en
#             with open("deformation_py.vtk", "a") as file:
#                 file.write(f"en = {en} attr = {attribute_en} psi_k_n = {psi_k_n} B_n = {B_n[en]} dN_n = {delta_N_n}\n")

#     return delta_N_n



# file = open("deformation_py.vtk", "r+")
# file.truncate(0)
# file = open("sigma_py.vtk", "r+")
# file.truncate(0)
# dN = get_dN("fatigue/output_000013.vtk")
# print(dN)
# """
attributes_size = 1001
file_lines_or =  open(f"delta_N/solid_circle_{attributes_size - 1}.txt", "r").read().split("\n")
# file_lines =  open(f"delta_N/solid_circle_fixed.txt", "r").read().split("\n")
# file_lines =  open(f"delta_N/ssf_sym.txt", "r").read().split("\n")
# file_lines = open(f"delta_N/vel_test_2.txt", "r").read().split("\n")
ntrh_name = "vel_test_2_smooth3e-5"
file_lines = open(f"delta_N/{ntrh_name}.txt", "r").read().split("\n")
dN_or = np.array([float(el) for el in file_lines_or[:-1]])
dN = np.array([float(el) for el in file_lines[:-1]])
len_dN_or = min(len(dN), len(dN_or))
dN_or = dN_or[:len_dN_or]
for en in range(1, len(dN)):
    dN[en] += dN[en - 1]
for en in range(1, len_dN_or):
    dN_or[en] += dN_or[en - 1]

plt.subplot()
x = np.arange(len(dN))
x_or = np.arange(len_dN_or)
plt.plot(x, dN, "-", c="b", label=f"{ntrh_name}")
# plt.plot(x, dN, "-", c="b", label="smooth")
plt.plot(x_or, dN_or, "-", c="k", label="original")
plt.xlabel("Итерации повреждений")
plt.ylabel(r"$\sum_{i \leq n} \Delta N^i$")
# plt.title(f"psi.size = {attributes_size}")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
# """
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Sigma_u = np.arange(350, 1160, 1)
# Psi = np.arange(0, 1, 0.01)

# dN_lim = 1000
# X, Y = [], []
# for sigma_u in Sigma_u:
#     for psi in Psi:
#         B_n = 1e-3 * np.power((sigma_u - 340) / (1160 - 340), 1 / 0.31)
#         dN = 0.5 * (1 - (0.5 * np.sqrt(psi) - psi)) / B_n
#         if dN <= dN_lim:
#             X.append(sigma_u)
#             Y.append(psi)

# plt.subplot()
# plt.xlim(350, 1160)
# plt.ylim(0, 1)
# plt.title(f"dN < {dN_lim}")
# plt.plot(X, Y, ".")

# sigma_u, psi = np.meshgrid(sigma_u, psi)
# B_n = 1e-3 * np.power((sigma_u - 340) / (1160 - 340), 1 / 0.31)
# dN = 0.5 * (1 - (0.5 * np.sqrt(psi) - psi)) / B_n

# surf = ax.plot_surface(sigma_u, psi, dN, cmap=mpl.cm.coolwarm, lw=0)
# ax.set_zlim(0, 1e2)
# ax.zaxis.set_major_locator(mpl.ticker.LinearLocator(10))
# plt.show()