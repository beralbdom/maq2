# UFF - UNIVERSIDADE FEDERAL FLUMINENSE
# Máquinas II - TEE00178, turma A1 2022/1
# Grupo 1, Projeto 2 - tarefa 4
# Elaborado por Bernardo Albuquerque

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint

# Parâmetros da máquina:
vfase = 220 / np.sqrt(3); npol = 4; fnom = 60; pnom = 6 * 735
ns = (120 * fnom) / npol; ws = ns * 2 * np.pi / 60

# Parâmetros do circuito equivalente (tarefa 1):
r1 = 0.24286; r2 = 0.20782
x1 = 0.50047; x2 = 0.50047; xm = 15.04629
rc = 130.32871

# Definindo equivalentes de Thévenin:
vthv = vfase * (xm / np.sqrt((r1 ** 2) + ((x1 + xm) ** 2)))
zthv = (1j * xm * (r1 + 1j * x1)) / (r1 + 1j * (x1 + xm))
rthv = zthv.real; xthv = zthv.imag

# 2a Lei de Newton para movimentos rotacionais:
def dwdt_carga(x, t):
    wr = x[0]
    s = (ws - wr) / ws
    Jcarga = 0.1
    Jmotor = 0.018
    
    tcarga = 17.81
    tperdas = 6
    tind = (3 * (vthv ** 2) * r2 / s) / (ws * (((rthv + r2 / s) ** 2) + ((xthv + x2) ** 2)))
    dwdt = (tind - tcarga - tperdas) / (Jmotor + Jcarga)

    return dwdt

def dwdt_vazio(x, t):
    wr = x[0]
    s = (ws - wr) / ws
    Jmotor = 0.018
    
    tperdas = 6
    tind = (3 * (vthv ** 2) * r2 / s) / (ws * (((rthv + r2 / s) ** 2) + ((xthv + x2) ** 2)))
    dwdt = (tind - tperdas) / Jmotor

    return dwdt

# Condição inicial e amostragem temporal:
x = [0]; t = np.linspace(0, 7e-1, 200)

# Resolvendo EDO e obtendo wr durante partida direta à vazio:
wr_vazio = odeint(dwdt_vazio, x, t)
wr_carga = odeint(dwdt_carga, x, t)

def corrente(caso):
  if (caso == 'vazio'):
    s = (ws - wr_vazio) / ws
    r2_ = r2 * (1 - s) / s
    z_eq = 1 / ((1 / (1j * xm)) + (1 / (rc)) + (1 / (r2_ + 1j * x2)))
    z_eq = z_eq + r1 + 1j * x1
    i = vfase / z_eq
  elif (caso == 'carga'):
    s = (ws - wr_carga) / ws
    r2_ = r2 * (1 - s) / s
    z_eq = 1 / ((1 / (1j * xm)) + (1 / (rc)) + (1 / (r2_ + 1j * x2)))
    z_eq = z_eq + r1 + 1j * x1
    i = vfase / z_eq
  return i

def torque(caso):
  if (caso == 'vazio'):
    s = (ws - wr_vazio) / ws
    tind = (3 * (vthv ** 2) * r2 / s) / (ws * (((rthv + (r2 / s)) ** 2) + ((xthv + x2) ** 2)))
  elif (caso == 'carga'):
    s = (ws - wr_carga) / ws
    tind = (3 * (vthv ** 2) * r2 / s) / (ws * (((rthv + (r2 / s)) ** 2) + ((xthv + x2) ** 2)))
  return tind

# Velocidade durante partida direta à vazio:
plt.rc('figure', figsize = (12, 6))
plt.plot(t, wr_vazio); plt.plot(t, wr_carga)
plt.xlabel('Tempo (s)'); plt.ylabel('Velocidade (rad/s)')
plt.title('Velocidade durante partida direta')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.legend(["Velocidade com motor à vazio", "Velocidade com carga de 75% da nominal"], fancybox = True, shadow = True)
plt.show()

# Corrente durante partida direta à vazio:
plt.rc('figure', figsize = (12, 6))
plt.plot(t, corrente('vazio').real); plt.plot(t, corrente('carga').real)
plt.xlabel('Tempo (s)'); plt.ylabel('Velocidade (rad/s)')
plt.title('Corrente durante partida direta')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.legend(["Corrente com motor à vazio", "Corrente com carga de 75% da nominal"], loc="lower right", fancybox = True, shadow = True)
plt.show()

# Corrente durante partida direta à vazio:
plt.rc('figure', figsize = (12, 6))
plt.plot(t, torque('vazio')); plt.plot(t, torque('carga'))
plt.xlabel('Tempo (s)'); plt.ylabel('Torque (N∙m)')
plt.title('Torque durante partida direta')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.legend(["Torque com motor à vazio", "Torque com carga de 75% da nominal"], loc="lower right", fancybox = True, shadow = True)
plt.show()