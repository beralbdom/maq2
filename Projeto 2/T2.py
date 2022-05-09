# UFF - UNIVERSIDADE FEDERAL FLUMINENSE
# Máquinas II - TEE00178, turma A1 2022/1
# Grupo 1, Projeto 2 - tarefa 2

from numpy import sqrt, sin, cos, pi, imag, real, linspace, array
from matplotlib import pyplot as plt

# Parâmetros da máquina:
vfas = 220 / sqrt(3); npol = 4; fnom = 60; pnom = 6 * 745.666
nsin = (120 * fnom) / npol; wsin = nsin * 2 * pi / 60

# Parâmetros do circuito equivalente (tarefa 1):
r1 = 0.24286; r2 = 0.20782
x1 = 0.50047; x2 = 0.50047; xm = 15.04629
rc = 130.32871

# Definindo equivalentes de Thévenin:
vthv = vfas * (xm / sqrt((r1 ** 2) + ((x1 + xm) ** 2)))
zthv = (1j * xm * (r1 + 1j * x1)) / (r1 + 1j * (x1 + xm))
rthv = zthv.real; xthv = zthv.imag

# Definindo vetor de velocidades:
s = linspace(0.001, 1, 150)
wmec = []
for i in range(len(s)):
    wmec.append((1 - s[i]) * wsin)

# Definindo torques de partida, nominal e máximo:
tpar = (3 * (vthv ** 2) * r2) / (wsin * (((rthv + r2) ** 2) + ((xthv + x2) ** 2)))
tnom = pnom / wsin # !!!! Não sei se tá certo !!!!

smax = r2 / sqrt(rthv ** 2 + (xthv + x2) ** 2)
wmax = (1 - smax) * wsin
tmax = (3 * (vthv ** 2) * r2 / smax) / (wsin * (((rthv + r2 / smax) ** 2) + ((xthv + x2) ** 2)))

# Definindo vetor de torques (Chapman, seção 6.5):
tind = []
for i in range(len(s)):
    t = (3 * (vthv ** 2) * r2 / s[i]) / (wsin * (((rthv + r2 / s[i]) ** 2) + ((xthv + x2) ** 2)))
    tind.append(t)

# Definindo função para cálculo do torque com parâmetros variáveis:
def _tind(v = vfas, f = fnom, r1 = r1, r2 = r2, rc = rc, x1 = x1, x2 = x2, xm = xm):
    nsin = (120 * f) / npol; wsin = nsin * 2 * pi / 60
    vthv = v * (xm / sqrt((r1 ** 2) + ((x1 + xm) ** 2)))
    zthv = (1j * xm * (r1 + 1j * x1)) / (r1 + 1j * (x1 + xm))
    rthv = zthv.real; xthv = zthv.imag

    _tind = []
    for i in range(len(s)):
        t = (3 * (vthv ** 2) * r2 / s[i]) / (wsin * (((rthv + r2 / s[i]) ** 2) + ((xthv + x2) ** 2)))
        _tind.append(t)

    return _tind

# Plotando curva característica:
plt.rc('figure', figsize = (12, 6))
plt.plot(wmec, tind) 
plt.plot(s[0], tpar, 'go'); plt.plot(wmax, tmax, 'ro'); plt.plot(wsin - 4, tnom, 'bo');
plt.text(s[0], tpar - 5, 'Tpar'); plt.text(wmax - 7, tmax - 5, 'Tmax'); plt.text(wsin - 24, tnom + 5, 'Tnom')
plt.text(s[0], tpar - 10, '{} N*m'.format(round(tpar, 2)));
plt.text(wmax - 7, tmax - 15, '{} N*m\n{} rad/s'.format(round(tmax, 2), round(wmax, 2)));
plt.text(wsin - 24, tnom, '{} N*m'.format(round(tnom, 2))); 

plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Curva característica de conjugado vs. velocidade')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

# Plotando curva característica variando parâmetros: 
# Deve ter um jeito melhor de fazer isso... talvez adicionando o plot dentro da função? 🤔
plt.plot(wmec, tind); plt.plot(wmec, _tind(v = 200))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: Vfase = 250 V')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

plt.plot(wmec, tind); plt.plot(wmec, _tind(v = 300))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: Vfase = 300 V')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()
 
plt.plot(wmec, tind); plt.plot(wmec, _tind(f = 30))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: Freq = 30 Hz')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()
 
plt.plot(wmec, tind); plt.plot(wmec, _tind(f = 90))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: Freq = 90 Hz')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()
 
plt.plot(wmec, tind); plt.plot(wmec, _tind(r1 = 2 * r1))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: R1 dobrado')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()
 
plt.plot(wmec, tind); plt.plot(wmec, _tind(r1 = r1 / 2))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: R1 pela metade')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

plt.plot(wmec, tind); plt.plot(wmec, _tind(r2 = 2 * r2))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: R2 dobrado')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

plt.plot(wmec, tind); plt.plot(wmec, _tind(r2 = r2 / 2))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: R2 pela metade')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

plt.plot(wmec, tind); plt.plot(wmec, _tind(rc = 2 * rc))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: Rc dobrado')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

plt.plot(wmec, tind); plt.plot(wmec, _tind(rc = rc / 2))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: Rc pela metade')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

plt.plot(wmec, tind); plt.plot(wmec, _tind(x1 = 2 * x1))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: X1 e X2 dobrados')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

plt.plot(wmec, tind); plt.plot(wmec, _tind(x1 = x1 / 2))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: X1 e X2 pela metade')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

plt.plot(wmec, tind); plt.plot(wmec, _tind(xm = 2 * xm))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: Xm dobrado')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()

plt.plot(wmec, tind); plt.plot(wmec, _tind(xm = xm / 2))
plt.xlabel('Velocidade mecânica (rad/s)'); plt.ylabel('Torque induzido (N*m)')
plt.title('Conjugado vs. velocidade: Xm pela metade')
plt.grid(color = 'gray', linestyle = '--', linewidth = 1, alpha = 0.33)
plt.show()