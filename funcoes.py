#Equação cúbica
#(a)*x^3 + (b)*x^2 +(c)*x + d = 0

import cmath
import math
import Banco
import numpy as np
import csv
#from random import randint
#from scipy.interpolate import interp1d

def FRoot3(a1,b1,c1,d1):

    v=[]

    delta = (18 * a1 * b1 * c1 * d1) - (4 * (b1 ** 3) * d1) + ((b1 ** 2) * (c1 ** 2)) - (4 * a1 * (c1 ** 3)) - (
                27 * (a1 ** 2) * (d1 ** 2))

    delta0 = (b1 ** 2) - (3 * a1 * c1)
    delta1 = (2 * (b1 ** 3)) - (9 * a1 * b1 * c1) + (27 * (a1 ** 2) * d1)
    C1 = ((delta1 + ((delta1 ** 2) - (4 * (delta0 ** 3))) ** (1 / 2)) / 2) ** (1 / 3)
    C2 = ((delta1 - ((delta1 ** 2) - (4 * (delta0 ** 3))) ** (1 / 2)) / 2) ** (1 / 3)
    u = (-1 + (cmath.sqrt(-3))) / 2

    if C1.imag < (10 ** -10):
        for n in range(0, 3):
            raiz = ((b1 + (u ** n) * C1) + (delta0 / ((u ** n) * C1))) / (-3 * a1)
            v.append(raiz.real)
    else:
        for n in range(0, 3):
            raiz = ((b1 + (u ** n) * C2) + (delta0 / ((u ** n) * C2))) / (-3 * a1)
            v.append(raiz.real)

    v = sorted(v)
    return v

def Root3(a1,b1,c1,d1):

    v = np.roots([a1, b1, c1, d1])
    v = sorted(v.real)
    return v

def Van(P,T,R,Tc,Pc):

    a = (27 * ((R * Tc) ** 2)) / (64 * Pc)
    b = (R * Tc) / (8 * Pc)

    a1 = -P
    b1 = (R * T + P * b)
    c1 = -a
    d1 = a * b

    return a1,b1,c1,d1

def SRK(P,T,R,Tc,Pc,W):

    Tr = (T / Tc)
    Om = 0.48 + 1.574 * W - 0.176 * W ** 2
    a = 0.42748 * (((R ** 2) * (Tc ** 2)) / (Pc)) * ((1 + Om * (1 - (Tr ** (1 / 2)))) ** 2)
    b = 0.08664 * ((R * Tc) / Pc)

    a1 = -P
    b1 = (R * T)
    c1 = (R * T * b - a + P * (b ** 2))
    d1 = a * b

    return a1, b1, c1, d1

def Peng(P,T,R,Tc,Pc,W):

    Tr = (T / Tc)
    Om = 0.37464 + 1.54226 * W - 0.26992 * W ** 2
    a = 0.45724 * (((R ** 2) * (Tc ** 2)) / (Pc)) * ((1 + Om * (1 - (Tr ** (1 / 2)))) ** 2)
    b = 0.0778 * ((R * Tc) / Pc)

    a1 = -P
    b1 = (R * T) - (2 * b * P) + (b * P)
    c1 = R * T * 2 * b - a + P * (b ** 2)
    d1 = a * b - R * T * (b ** 2) + 2 * (b ** 2) * P - (b ** 3)

    return a1, b1, c1, d1

def Tsai(P,T,R,Tc,Pc,W):

    Tr = (T / Tc)
    Om = 0.20473 + 0.83548 * W - 0.18472 * W ** 2 + 0.16675 * W ** 3 - 0.09881 * W ** 4
    alfa = (1 + Om * (1 - Tr)) ** 2
    a = 0.45724 * (((R ** 2) * (Tc ** 2)) / Pc) * alfa
    b = 0.0778 * ((R * Tc) / Pc)

    TC1 = Banco.TC1['agua']
    TC2 = Banco.TC2['agua']
    TC3 = Banco.TC1['agua']
    t = ((R * Tc) / Pc) * (TC1 + TC2 * (1 - (Tr ** (2 / 3))) + TC3 * (1 - (Tr ** (2 / 3))) ** 2)

    a1 = -P
    b1 = (R * T) - (2 * b * P) + (b * P)
    c1 = R * T * 2 * b - a + P * (b ** 2)
    d1 = a * b - R * T * (b ** 2) + 2 * (b ** 2) * P - (b ** 3)

    return a1, b1, c1, d1, t

def Vera(P, T, R, Tc, Pc, W):
    Tr = (T / Tc)
    k0 = 0.378893 + 1.4897153 * W - 0.17131848 * W ** 2 + 0.0196554 * W ** 3
    Om = k0 + 0.0089 * (1 + (Tr ** .5)) * (.7 - Tr)
    a = 0.457235 * (((R ** 2) * (Tc ** 2)) / (Pc)) * ((1 + Om * (1 - (Tr ** (1 / 2)))) ** 2)
    b = 0.0778 * ((R * Tc) / Pc)

    a1 = 1
    b1 = b - ((R * T) / P)
    c1 = ((a - 2 * b * R * T) / P) - (3 * (b ** 2))
    d1 = (((b ** 2) * R * T - a * b) / P) + b ** 3

    return a1, b1, c1, d1

def PacoteCte(Subs):

    Subs = Subs
    R = 8.314462  # m³*kPa/K*mol
    Tc = Banco.Tc[Subs]  # K
    Pc = Banco.Pc[Subs] * (10 ** 3)  # kPa
    MM = Banco.MM[Subs]  # g/mol
    W = Banco.W[Subs]  # ω
    TC1 = Banco.TC1[Subs]
    TC2 = Banco.TC2[Subs]
    TC3 = Banco.TC1[Subs]

    return R,Tc,Pc,MM,W,TC1,TC2,TC3

def Goff(T):
    Psat = -7.90298 * (373.16 / T - 1) + 5.02808 * (math.log10(373.16 / T)) - (1.3816E-7) * ((10 ** (11.344 * (1 - T / 373.16))) - 1) + ((8.1328E-3) * (10 ** ((-3.49149) * (373.16 / T - 1)) - 1)) + math.log10(1013.246)
    Psat = (10 ** Psat) / 10
    return Psat

def PsatDIPPR(T):
    x = []
    y = []
    with open('PsatDIPPR.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=';', quoting=csv.QUOTE_NONE)
        for row in reader:
            x.append(float(row[0].replace(',', '.')))
            y.append(float(row[1].replace(',', '.')))


    #v = interp1d(x, y, kind='cubic')
    v = np.interp(T, x, y, period=360)

    return v

def VolLiqDIPPR(T):
    x = []
    y = []
    with open('VolLiqDIPPR.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=';', quoting=csv.QUOTE_NONE)
        for row in reader:
            x.append(float(row[0].replace(',', '.')))
            y.append(float(row[1].replace(',', '.')))

    #v = interp1d(x, y, kind='cubic')
    v=np.interp(T, x, y, period=360)

    return v

def VolLiqDIPPR_P(Tmin,Tmax):
    x = []
    y = []
    with open('VolLiqDIPPR.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=';', quoting=csv.QUOTE_NONE)
        for row in reader:
            aux=float(row[0].replace(',', '.'))
            if aux>Tmin and aux<Tmax:
                x.append(aux)
                y.append(float(row[1].replace(',', '.')))


    return [x, y]

def FuncCompa (P,T,esc):
    result=[]
    Subs = 'agua'
    R = 8.314462  # m³*kPa/K*mol
    Tc = Banco.Tc[Subs]  # K
    Pc = Banco.Pc[Subs] * (10 ** 3)  # kPa
    MM = Banco.MM[Subs]  # g/mol
    W = Banco.W[Subs]  # ω
    TC1 = Banco.TC1[Subs]
    TC2 = Banco.TC2[Subs]
    TC3 = Banco.TC1[Subs]

    if Subs =='agua':
        # Goff
        Psat = Goff(T)
        result.append(Psat)

        # DIPPR Psat
        if T>647.13:
            Psat='Temperatura Alta'
            result.append(Psat)
        elif T<273.15:
            Psat = 'Temperatura Baixa'
            result.append(Psat)
        else:
            Psat= PsatDIPPR(T)
            result.append(Psat)

    else:
        Psat = "Não é agua"
        result.append(Psat)
        result.append(Psat)
        #2x para manter o tamanho do vetor



    # Van
    a1, b1, c1, d1 = Van(P, T, R, Tc, Pc)
    v = Root3(a1, b1, c1, d1)
    v[0] = v[0] / MM
    v[2] = v[2] / MM
    result.append(v[esc])


    # SRK
    a1, b1, c1, d1 = SRK(P, T, R, Tc, Pc, W)
    v = Root3(a1, b1, c1, d1)
    v[0] = v[0] / MM
    v[2] = v[2] / MM
    result.append(v[esc])

    # Peng-Robinson

    a1, b1, c1, d1 = Peng(P, T, R, Tc, Pc, W)
    v = Root3(a1, b1, c1, d1)
    v[0] = v[0] / MM
    v[2] = v[2] / MM
    result.append(v[esc])

    # Tsai-Chen (Peng-Robinson com volume de translação)

    a1, b1, c1, d1, t = Tsai(P, T, R, Tc, Pc, W)
    v = Root3(a1, b1, c1, d1)
    #Veja se é menos mesmo:

    v[0] = (v[0]+t) / MM
    v[2] = (v[2]+t) / MM
    result.append(v[esc])

    # Vera
    a1, b1, c1, d1 = Vera(P, T, R, Tc, Pc, W)
    v = Root3(a1, b1, c1, d1)
    v[0] = v[0] / MM
    v[2] = v[2] / MM
    result.append(v[esc])

    # DIPPR

    if esc == 0 and T < 647.13 and T > 273.15:
        v = VolLiqDIPPR(T)
        result.append(v)




    return result

def Erro(vDIPPR,vCorr):

    Erro=abs((vDIPPR-vCorr)/vDIPPR)
    return Erro

def Grafico(Rmin,Rmax):

    Rmin = round(Rmin)
    Rmax = round(Rmax)

    # Qual substância?
    Subs = 'agua'
    R = 8.314462  # m³*kPa/K*mol
    Tc = Banco.Tc[Subs]  # K
    Pc = Banco.Pc[Subs] * (10 ** 3)  # kPa
    MM = Banco.MM[Subs]  # g/mol
    W = Banco.W[Subs]  # ω
    TC1 = Banco.TC1[Subs]
    TC2 = Banco.TC2[Subs]
    TC3 = Banco.TC1[Subs]

    Vanl = []
    Vanv = []
    SRKl = []
    SRKv = []
    PRl = []
    PRv = []
    TChel = []
    TChev = []
    DIPPRx = []
    DIPPRy = []

    for T in range(Rmin, Rmax):
        P = ac.Goff(T)
        # Van
        a1, b1, c1, d1 = ac.Van(P, T, R, Tc, Pc)
        v = ac.Root3(a1, b1, c1, d1)
        v[0] = v[0] / MM
        v[2] = v[2] / MM
        Vanl.append(v[0])
        Vanv.append(v[2])

        # SRK
        a1, b1, c1, d1 = ac.SRK(P, T, R, Tc, Pc, W)
        v = ac.Root3(a1, b1, c1, d1)
        v[0] = v[0] / MM
        v[2] = v[2] / MM
        SRKl.append(v[0])
        SRKv.append(v[2])

        # Peng-Robinson

        a1, b1, c1, d1 = ac.Peng(P, T, R, Tc, Pc, W)
        v = ac.Root3(a1, b1, c1, d1)
        v[0] = v[0] / MM
        v[2] = v[2] / MM
        PRl.append(v[0])
        PRv.append(v[2])

        # Tsai-Chen (Peng-Robinson com volume de translação)

        a1, b1, c1, d1, t = ac.Tsai(P, T, R, Tc, Pc, W)
        v = ac.Root3(a1, b1, c1, d1)
        v[0] = v[0] + t
        v[2] = v[2] + t
        v[0] = v[0] / MM
        v[2] = v[2] / MM
        TChel.append(v[0])
        TChev.append(v[2])

    plt.plot(PRl, range(Rmin, Rmax), label='Tsai-Chen ', color='r')
    plt.plot(PRv, range(Rmin, Rmax), color='r')
    plt.ylabel('Temperatura (K)')
    plt.xlabel('Volume específico (kg/m³)')
    plt.legend()
    plt.show()

    return (    Vanl,    Vanv, SRKl, SRKv , PRl,PRv ,TChel ,TChev)

def MudarCor(Contador_Cor):
    if Contador_Cor == 0:
        Cor = 'r'
    if Contador_Cor==1:
        Cor='b'
    if Contador_Cor==2:
        Cor='g'
    if Contador_Cor==3:
        Cor="m"
    if Contador_Cor==4:
        Cor="c"
    if Contador_Cor==5:
        Cor="y"
    if Contador_Cor==6:
        Cor="k"
    if Contador_Cor==7:
        Cor='xkcd:sky blue'
    if Contador_Cor==8:
        Cor='xkcd:bright green'
    if Contador_Cor==9:
        Cor='xkcd:mustard'
    if Contador_Cor==10:
        Cor='xkcd:baby blue'
    if Contador_Cor==11:
        Cor='xkcd:hunter green'
    if Contador_Cor==12:
        Cor='xkcd:crimson'
    if Contador_Cor==13:
        Cor='xkcd:baby blue'
    if Contador_Cor==14:
        Cor='xkcd:goldenrod'

    return Cor

def GoffDIPPR():
        x = []
        y = []
        with open('fotograf.csv') as csvfile:
            reader = csv.reader(csvfile, delimiter=';', quoting=csv.QUOTE_NONE)
            for row in reader:
                x.append(float(row[0].replace(',', '.')))
                y.append(float(row[1].replace(',', '.')))

        return x,y






