__author__ = 'LuizArthur'

import cmath as cm
from math import pow
import numpy as np
import matplotlib.pyplot as pp

class CALC:
    def __init__(self, Var_txt, count_plot, name):
        self.Var = dict()

        for var in Var_txt:
            if var[0] == "comp":
                self.Var[var[0]] = []
                for i in range(1,len(var)):
                    self.Var[var[0]].append(float(var[i])/100)
            else:
                self.Var[var[0]] = float(var[1])

        self.file_name = name
        self.r = self.Var["r"]
        self.L = self.Var["L"]
        self.C = self.Var["C"]
        self.g = self.Var["g"]
        self.l = self.Var["l"]
        self.V1 = self.Var["V1"]
        self.Vb = self.V1
        self.Sb = 100e6
        self.fp = self.Var["fp"]
        self.f = self.Var["f"]
        self.prec = self.Var["prec"]
        self.w = 2 * cm.pi * self.f
        self.count_plot = count_plot

        self.U1 = self.V1 / cm.sqrt(3).real
        self.Z = (self.r + 1j * self.w * self.L) * self.l
        self.Y = (self.g + 1j * self.w * self.C) * self.l
        self.Zc = cm.sqrt(self.Z / self.Y)
        self.gamal = cm.sqrt(self.Z * self.Y)

        self.A = cm.cosh(self.gamal)
        self.B = self.Zc * cm.sinh(self.gamal)
        self.C = 1 / self.Zc * cm.sinh(self.gamal)
        self.D = self.A

        self.result = self.calc_Vr(self.A, self.B)
        self.result2 = self.compSeM()
        self.result3 = self.compSeEx()
        #self.ang_pot = self.calc_ang_pot(self.A,self.B,np.array(self.result[0])/3,np.array(self.result[1])/np.sqrt(3))

        self.plot_all_Vr()

        """
        f = pp.figure(1)
        pp.plot(self.result[0],self.result[1],self.result[0],self.result[2])
        #print(self.result[4])
        pp.draw()
        """

    def calc_Vr(self, A, B):
        Pr = []
        Qr = []
        result = []
        result2 = []
        ang_pot = []

        absA = abs(A)
        alfa = cm.phase(A)
        absB = abs(B)
        beta = cm.phase(B)


        i = 0
        while i > -1:

            Praux = i * 10e5 / 3 * self.fp
            Qraux = cm.tan(cm.acos(self.fp)).real * Praux
            a = pow((absA /absB), 2)
            b = (2 * absA / absB * (Praux * cm.cos(beta - alfa) + Qraux * cm.sin(beta - alfa)
                                              - pow(self.U1, 2) / (2 * absA * absB))).real
            c = pow(Praux, 2) + pow(Qraux, 2)

            raiz = pow(b, 2) - 4 * a * c

            if raiz < 0:
                result2[len(result2) - 1] = result[len(result) - 1]
                break

            Pr.append(Praux)
            Qr.append(Qraux)

            coeff = [a, 0, b, 0, c]

            aux = np.roots(coeff)
            result.append(max(abs(aux * cm.sqrt(3).real)))
            result2.append(min(abs(aux * cm.sqrt(3).real)))
            ang_pot.append(self.calc_ang_pot(A,B,Praux,result[len(result) - 1]/np.sqrt(3).real))

            # print(result[len(result)-1])

            i = i + self.prec

        return [3 * np.array(Pr)/self.Sb, np.array(result)/self.Vb, np.array(result2)/self.Vb, 3*np.array(Qr)/self.Sb, ang_pot]

    def calc_ang_pot(self,A,B,Pr,Vr):
        absA = abs(A)
        alfa = cm.phase(A)
        absB = abs(B)
        beta = cm.phase(B)

        ang_pot = (beta - np.arccos((Pr + absA*np.power(Vr,2)/absB*np.cos(beta-alfa))/(self.U1*Vr/absB)))
        ang_pot = (np.degrees(ang_pot))
        #print(ang_pot)

        return ang_pot

    def compSeM(self):
        A1 = cm.sqrt((self.A+1)/2)
        B1 = self.B/(2*A1)
        C1 = (A1*A1-1)/B1

        Xcmax = (-self.B.imag/((-1j*(self.A+1)/2).imag))

        results = []

        for perc in self.Var["comp"]:
            Xc = perc*Xcmax
            Al = A1*A1 + C1*B1 - 1j*C1*A1*Xc
            Bl = self.B - 1j*(self.A+1)*Xc/2

            result = self.calc_Vr(Al, Bl)
            results.append(result)

        return results


    def compSeEx(self):
        A1 = cm.sqrt((self.A+1)/2)
        B1 = self.B/(2*A1)
        C1 = (A1*A1-1)/B1

        a = -self.C.imag
        b = (-1j*2*self.A).imag
        c = self.B.imag
        coeff = [a,b,c]
        Xcmax = max(np.roots(coeff))

        results = []

        for perc in self.Var["comp"]:
            Xc = perc*Xcmax
            Al = self.A - 1j*self.C*Xc
            Bl = self.B - 1j*2*self.A*Xc - self.C*np.power(Xc,2)
            result = self.calc_Vr(Al, Bl)
            results.append(result)

        return results


    def plot_all_Vr(self):
        i=0
        for perc in self.Var["comp"]:
            f = pp.figure(self.count_plot)
            for j in range(1,3):
                pp.plot(self.result2[i][0], self.result2[i][j], color="green")
                pp.plot(self.result3[i][0], self.result3[i][j], color="red")
                pp.plot(self.result[0], self.result[j], color="blue")
                pp.title(self.file_name+"\nPot x Vr")
            pp.text(self.result2[i][0][int(.65*len(self.result2[i][0])-1)],
                self.result2[i][1][0],"Curva azul: Sem comp.")
            pp.text(self.result2[i][0][int(.65*len(self.result2[i][0])-1)],
                self.result2[i][1][int(0.25*len(self.result2[i][1])-1)],
                "Curva verde: Comp. meio "+str(perc*100)+"%")
            pp.text(self.result2[i][0][int(.65*len(self.result2[i][0])-1)],
                self.result2[i][1][int(0.5*len(self.result2[i][1])-1)],
                "Curva vermelha: Comp. extr. "+str(perc*100)+"%")
            pp.xlim([0,self.result2[i][0][int(len(self.result2[i][0])-1)]+
                10+1])
            i = i+1
            self.count_plot = self.count_plot+1

        i = 0
        for perc in self.Var["comp"]:
            f = pp.figure(self.count_plot)
            pp.plot(self.result2[i][0], self.result2[i][4], color="green")
            pp.plot(self.result3[i][0], self.result3[i][4], color="red")
            pp.plot(self.result[0], self.result[4], color="blue")
            pp.title(self.file_name+"\nPot x Angulo de Pot")
            pp.text(self.result2[i][0][int(.1*len(self.result2[i][0])-1)],
                self.result[4][len(self.result[4])-1],
                "Curva azul: Sem comp.")
            pp.text(self.result2[i][0][int(.1*len(self.result2[i][0])-1)],
                self.result[4][int(0.995*len(self.result[4])-1)],
                "Curva verde: Comp. meio "+str(perc*100)+"%")
            pp.text(self.result2[i][0][int(.1*len(self.result2[i][0])-1)],
                self.result[4][int(0.98*len(self.result[4])-1)],
                "Curva vermelha: Comp. extr. "+str(perc*100)+"%")
            pp.xlim([0,self.result2[i][0][int(len(self.result2[i][0])-1)]+
                10+1])
            i = i+1
            self.count_plot = self.count_plot+1

        pp.draw()

    def show_count_plot(self):
        return self.count_plot
