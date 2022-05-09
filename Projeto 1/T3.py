# UFF - UNIVERSIDADE FEDERAL FLUMINENSE
# Máquinas II - TEE00178, turma A1 2022/1
# Grupo 1, Projeto 1: Máquina H

from matplotlib import pyplot
from numpy import array, arange, sin, cos, pi

# Definindo parâmetros
npol = 2                                # [DADO] N. de polos.
nbob = 14                               # [DADO] N. de ranhuras por agrupamento de fase.
fpas = (19 / 21)                        # [DADO] Passo da máquina.
nran = 42                               # [DADO] N. de ranhuras.
ppol = (360 * pi / 180) / npol          # Passo polar em radianos.
penc = (fpas * ppol)                    # Passo encurtado em radianos.
gamm = (360 * pi / 180) / nran          # Distância entre ranhuras em radianos.
freq = 60                               # Frequência da rede em Hz.
we   = 2 * pi * freq                    # Velocidade angular.
kd   = (sin(nbob * gamm / 2)) / \
       (nbob * sin(gamm / 2))           # Fator de distribuição.
har1 = ((2 * 1 * nran) / npol) - 1      # Índices das harmônicas mais problemáticas.
har2 = ((2 * 1 * nran) / npol) + 1

print("Fator de distribuição:", kd)
print("\nÍndices de harmônicas de freq. mais baixa:", har1, "e", har2)

indf = []                               # Vetor de tamanho nran com índices de Fourier não múltiplos de 3.
for i in range(1, nran * 3, 2):
    if i % 3 != 0: 
        indf.append(i)
indf = array(indf)

print("\nÍndices de Fourier não múlt. de 3:\n", indf)

teta = []                               # Vetor representando ângulos elétricos para uma volta completa com passo de 1 radiano.
for i in arange(0, 2 * pi, pi / 180):
    teta.append(i)
teta = array(teta)

kp = []                                 # Fator de passo para cada harmônico.
for i in range(len(indf)):
    kp.append(sin(penc * indf[i] / 2))
kp = array(kp)

print("\nFatores de passo para cada harmônico:\n", kp)

# Definindo funções para cálculo da FMM.
def fmm_c(we, teta):                    # FMM para enrolamentos concentrados.
    soma = 0
    for i in range(len(indf)):
        an = (1 / (npol * pi * indf[i])) * (3 * sin(pi * indf[i] / 2) - sin((3 * pi * indf[i]) / 2))
        fa = an * (cos(3 * indf[i] * teta) * cos(we))
        fb = an * (cos(3 * indf[i] * teta - (2 * pi / 3)) * cos(we - (2 * pi / 3)))
        fc = an * (cos(3 * indf[i] * teta + (2 * pi / 3)) * cos(we + (2 * pi / 3)))
        soma += (fa + fb + fc)
    return soma

def fmm_d(we, teta):                    # FMM para enrolamentos distribuídos.
    soma = 0
    for i in range(len(indf)):
        an = ((kd * kp[i]) / (npol * pi * indf[i])) * (3 * sin(pi * indf[i] / 2) - sin((3 * pi * indf[i]) / 2))
        fa = an * (cos(3 * indf[i] * teta) * cos(we))
        fb = an * (cos(3 * indf[i] * teta - (2 * pi / 3)) * cos(we - (2 * pi / 3)))
        fc = an * (cos(3 * indf[i] * teta + (2 * pi / 3)) * cos(we + (2 * pi / 3)))
        soma += (fa + fb + fc)
    return soma

def fmm_d_h(we, teta):                  # FMM para enrolamentos distribuídos sem harmônicas críticas.
    soma = 0
    for i in range(len(indf)):
        if indf[i] != int(har1) and indf[i] != int(har2):
            an = ((kd * kp[i]) / (npol * pi * indf[i])) * (3 * sin(pi * indf[i] / 2) - sin((3 * pi * indf[i]) / 2))
            fa = an * (cos(3 * indf[i] * teta) * cos(we))
            fb = an * (cos(3 * indf[i] * teta - (2 * pi / 3)) * cos(we - (2 * pi / 3)))
            fc = an * (cos(3 * indf[i] * teta + (2 * pi / 3)) * cos(we + (2 * pi / 3)))
            soma += (fa + fb + fc)
    return soma

# Plotando gráficos.
pyplot.plot(teta, fmm_c(we, teta))
pyplot.xlabel('Posição angular'); pyplot.ylabel('FMM'); pyplot.title('FMM: Enrolamentos concentrados')
pyplot.show()

pyplot.plot(teta, fmm_d(we, teta))
pyplot.xlabel('Posição angular'); pyplot.ylabel('FMM'); pyplot.title('FMM: Enrolamentos distribuídos')
pyplot.show()

pyplot.plot(teta, fmm_d_h(we, teta))
pyplot.xlabel('Posição angular'); pyplot.ylabel('FMM'); pyplot.title('FMM: Enrolamentos distribuídos s/ harmônicas 41 e 43')
pyplot.show()