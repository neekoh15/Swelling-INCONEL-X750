from math import pi, exp, sqrt, log
import tabulate
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


class CavitySwelling:

    # RR
    rR = [((2/200) * (773 - (i + 473)) + 5) if i > 200
          else
          ((8/200) * (773 - (i + 473)) + 2)
        for i in range(800)]

    N0 = 6e14

    def __init__(self,
                    he_fit,
                    dpa_fit,
                    z,
                    uf, # fmd rate  = 0.107 para H347 | = 0.108 para X750
                    omega = 1.14E-29, # fraccion atomica de vacacias
                    s = 1, # cavity surface energy ?
                    efv = 1.6, # energia de formacion de vacancia
                    rM = 5000, # aka R
                    r = 380,
                    e = 1.4,# energia de migracion de vacancias
                    f = 1, # eficiencia de cascada
                    teol = 56.25,
                    _N0 = 6e14,# densida de dislocacion ?
                    fi = 0.01): #fraccion inincial
        self.teol = teol
        self.Teol = teol * 365 * 24 * 3600
        self.he_fit = he_fit
        self.dpa_fit = dpa_fit
        self.z = z
        self.uf = uf #fraccion desconocida
        self.omega = omega
        self.s = s
        self.efv = efv
        self.rM = rM
        self.r = r
        self.e = e
        self.f = f
        self.fi = fi
        self.N0 = _N0


    def rateHe(self, t):
        x = self.he_fit
        tot = 0
        for i in range(1, len(x)):
            tot += i * x[i] * self.td(t)**(i - 1)
        return tot

    def he(self, t):
        x = self.he_fit
        tot = 0
        for i in range(0, len(x)):
            tot += x[i] * self.td(t)**i
        return tot

    def dpa(self, t):
        x = self.dpa_fit
        tot = 0
        for i in range(0, len(x)):
            tot += x[i] * self.td(t)**i
        return tot

    def rateDpa(self, t):
        x = self.dpa_fit
        tot = 0
        for i in range(1, len(x)):
            tot += i * x[i] * self.td(t)**(i - 1)
        return tot

    def polynomialFit(x, y, n):
        x = np.array(x)
        y = np.array(y)
        p = np.polyfit(x, y, n)
        return np.poly1d(p)

    def td(self, t):
        return self.teol * t

    #ck
    def GHe(self, f, t):
        # usar rate en casos no lineales
        return f * (self.rateHe(t)  / (24 * 365 * 3600 * 1E6))

    def integrate(f, a, b, n):
        h = (b - a) / n
        s = 0
        for i in range(n):
            s += f(a + i * h)
        return s * h

    #ck
    def heTot(self, f, t):
        # usar integracion en casos no lineales
        return (CavitySwelling.integrate(lambda tt: self.GHe(f,tt) , 0, t, 1000)) * self.Teol

    #ck
    def CI(self, f, z, rM, d, ro, t, r, e):
        rdv = self.Rd(z, t, e)
        rcv = self.Rc(d, ro, t)
        ssgbv = self.ssgb(rM, d, ro, z, t, e)
        return (self.G(f, t) / (((rdv * (1 + self.ba(z, t, e))) + rcv + ssgbv) * self.DI(z))) * self.Q(f, z, rM, d, ro, t, r, e) 

    #ck
    def CE(self, z, efv):
        return exp(-(efv) * 1.6E-19 / (1.38E-23 * z))

    #ck
    def C(self, f, z, rM, d, ro, t, r, e):
        return (
            self.G(f, t)/
            ((self.Rd(z, t, e) + self.Rc(d, ro, t) + self.ssgb(rM, d, ro, z, t, e)) * self.DV(z, e)) *
            self.Q(f, z, rM, d, ro, t, r, e)
        )

    #ck
    def G(self, f, t):
         return f * (self.rateDpa(t) / (24 * 365 * 3600))

    #ck - e
    def Rd(self, z, t, e):  # verificar e usada
        rRt = self.rR[int(z - 473)] * 1E14
        return (rRt  * 0.4) / (
            ((rRt - self.N0) / self.N0) * (e**(-self.rr(z) * t * 300)) + 1
        )

    #ck
    def Rc(self, d, ro, t):
        return self.ss(d, ro) * (self.dpa(t) / self.dpa(1)) #cambio 347H

    #ck
    def ssgb(self, rM, d, ro, z, t, e):
        # return (3 / (r * 1e-9)) * sqrt(ss(d, ro) * dpa(t))
        # redeclaracion
        return (3 / (rM * 1E-9)) * ((self.Rd(z, t, e) + self.Rc(d, ro, t))**0.5)

    #ck
    def DV(self, z, e):
        return 6 * 1E-6 * exp(-e * 1.6E-19 / (1.38E-23 * z))

    #ck
    def Q(self, f, z, rM, d, ro, t, r, e):
        nv = self.n(f, z, rM, d, ro, t, r, e)
        return (2/nv) * ((1 + nv)**0.5 - 1)

    #ck
    def rr(self, z):
        return 1 * (((z - 773)/z) + 1)

    #ck
    def ss(self, d, ro):
        return 4 * pi * ro * (d/2) * 1E-9 * 1E23

    #ck
    def n(self, f, z, rM, d, ro, t, r, e):
        rdv = self.Rd(z, t, e)
        rcv = self.Rc(d, ro, t)
        ssgbv = self.ssgb(rM, d, ro, z, t ,e)
        return (4 * self.a(z, r) * self.G(f, t)) / (
            ((rdv * (1 + self.ba(z, t,e)) + rcv + ssgbv) * (rdv + rcv + ssgbv) * (self.DV(z, e) * self.DI(z)))
        )

    #ck
    def a(self, z, r):
        return 1E19 * r * self.DI(z)

    #ck
    def ba(self, z, t, e):
        return (self.zIa(z, t, e) - self.zVa(z, t, e)) / self.zIa(z, t, e)

    #ck
    def DI(self, z):
        return 12E-6 * exp((-0.15 * 1.6E-19)/(1.38E-23 * z))

    #ck
    def zIa(self, z, t, e):
        return 2 * pi / log(2 * (((pi * self.Rd(z, t, e))**-0.5)/self.lIa(z)))

    #ck
    def zVa(self, z, t, e):
        return 2 * pi / log(2 * (((pi * self.Rd(z, t, e))**-0.5)/self.lVa(z)))

    #ck
    def lIa(self, z):
        return (573/z) * 3.6E-9

    #ck
    def lVa(self, z):
        return (573/z) * 6E-10


    def rho1(self ,z):
        return ((2.22662E31 * ((z - 273.0)**-11.5142))/1000.0) * 5.0

    #transpose matrix
    def transpose(matrix):
        return [list(i) for i in zip(*matrix)]

    # fin cjv
    # plot matplotlib

    def plot_graph(x, y, x_label, y_label, title, block=True, fn = 0):
        plt.figure(fn)
        plt.plot(x, y)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.show(block=block)

    #plot_graph([i + 573 for i in range(800)], rR, "i+573", "RRi", "J1")

    def run(self, silent=False):
        omega = self.omega # omega
        z = self.z  # temperatura en K
        rho1z = self.rho1(z)
        s = self.s  # sigma Surface Energy
        efv = self.efv
        rM = self.rM
        r = self.r
        e = self.e
        f = self.f #0.0448 * 1.6/2.53  # fraction
        uf = self.uf  # uv
        # initial radius cavity
        AGB = self.heTot(f, self.fi)
        AGBS = []
        AGBS.append(AGB)
        ybs = []
        ybs.append(0)
        pbs = []
        pbs.append(0)
        cgbs = []
        cgbs.append(0)
        cjvs = []
        cjvs.append(0)
        final = 100
        dpas = []
        for i in range(1, final + 2):
            hefi5 = self.heTot(f, i/100)
            vTerm = ((AGB * 6) / (pi * rho1z * 1E-4))**(1/3)
            YB = pi * ((2E-10)**3 / (6 * omega)) * (hefi5 / AGB)
            PB =  ((1 + YB + YB**2 - YB**3) / ((1 - YB) ** 3)) * ((hefi5 /(omega * 6.023E23 * AGB)) * 8.31 * z)#Resultado en Pascal #Ck
            if(PB < 0):
                PB = 0
            CGB = exp(-(efv * (1.6E-19) + ((PB - ((2*s) / ((vTerm/2) * 1E-9))) * omega)) / (1.38 * z * 1E-23))
            CJV = (self.Rc(vTerm, rho1z, i/100) * self.DV(z, e) *
            (self.C(uf, z, rM, vTerm, rho1z, i/100, r, e) + self.CE(z, efv) - CGB)) + ((-(self.Rc(vTerm, rho1z, i/100))) * self.DI(z) * self.CI(uf, z, rM, vTerm, rho1z, i/100, r, e))
            AGB = AGB + ((CJV * (self.Teol/100)) if CJV > 0 else 0)
            AGBS.append(AGB)
            ybs.append(YB)
            pbs.append(PB)
            cgbs.append(CGB)
            cjvs.append(CJV)
            dpas.append(self.dpa(i/100))

        rv = tabulate.tabulate(CavitySwelling.transpose([range(0,final + 2), AGBS, ybs, pbs, cgbs, cjvs, dpas ]), headers=["i", "Vol", "YB", "PB", "CGB", "CJV", "DPA"])
        if not silent:
            print(rv)
            CavitySwelling.plot_graph([i/100  for i in range(0, final + 2)], AGBS, "Time (s)", "Volume", "Cavity Swelling", True, 1)
        self.deol =(AGB*6/pi/rho1z/0.0001)**0.333
        return rv, AGBS, ybs, pbs, cgbs, cjvs, dpas

if __name__ == "__main__":
    a = CavitySwelling()
    print(a.rateHe(0.5))
    print(a.rateDpa(0.5))
    print(a.GHe(1, 0.5))
    print(a.heTot(1, 0.01))
    print(a.dpa(1))
    a.run()
    #input("Press Enter to continue...")

