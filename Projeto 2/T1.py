# UFF - UNIVERSIDADE FEDERAL FLUMINENSE
# Máquinas II - TEE00178, turma A1 2022/1
# Grupo 1, Projeto 2 - tarefa 1
 
from numpy import sqrt, sin, cos, arccos, imag, real
 
# Definindo parâmetros da máquina:
npol = 4; fnom = 60; pnom = 6 * 745.666
v_cc = 13.6; i_cc = 28
v_bl = 25; f_bl = 15; i_bl = 28; p_bl = 1060
v_vz = 220; i_vz = 8.17; p_vz = 420
 
# Calc. elementos do circ. eq. (Chapman, seção 6.11):
r1 = v_cc / (2 * i_cc)
 
fp      = p_bl / (sqrt(3) * v_bl * i_bl)
z_bl_   = v_bl / (sqrt(3) * i_bl)
z_bl    = (z_bl_ * cos((arccos(fp))) + 1j * (z_bl_ * sin((arccos(fp)))))
r_bl    = z_bl.real;
x_bl_   = z_bl.imag
r2      = r_bl - r1
x_bl    = (fnom / f_bl) * x_bl_
x1 = x2 = 0.5 * x_bl
 
p_ce = 3 * r1 * i_vz ** 2
p_nu = p_vz - p_ce
rc   = (3 * (v_vz / sqrt(3)) ** 2) / p_nu
z_eq = v_vz / (sqrt(3) * i_vz)
xm   = z_eq - x1
 
# Printando elementos:
print("Elementos do circuito equivalente:", "\nR1 =", round(r1, 5), "Ω", \
"\nR2 =", round(r2, 5), "Ω","\nX1 =", round(x1, 5), "Ω", "\nX2 =", \
round(x2, 5), "Ω", "\nXm =", round(xm, 5), "Ω", "\nRc =", round(rc, 5), "Ω",)