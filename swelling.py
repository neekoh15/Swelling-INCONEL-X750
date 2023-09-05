""" 
Algoritmo para calculo de Cavity Swelling para aleaciones H347 o Inconel X750.

Autores: Griffiths, M. & Ramos Nervi, J.

Modificado por:
    Martinez, Nicolas Agustin. 

"""

from math import pi, exp, sqrt, log
import tabulate
import matplotlib.pyplot as plt
import numpy as np


class CavitySwelling:

    """
    Genera un objeto con parametros fisicos para el calculo de Swelling para una determinada temperatura.

        he_fit: coeficientes de un polinomio representativos de la curva de generacion de Helio a lo largo de su ciclo de vida.

        dpa_fit: coeficientes del polinomio representativo de la curva de dpa a lo largo de su ciclo de vida.

        z: temperatura en grados Kelvin.

        uf: fmd rate.

        omega: fraccion atomica de vacancias.

        s: cavity surface energy.

        efv: energia de formacion de vacancia.

        rM: R ?????

        r: parametro de recombinacion.

        f: eficiencia de cascada.

        teol: tiempo end-of-life

        _N0: densidad de dislocacion.

        fi: fraccion inicial (% de tiempo de vida inicial)
        """

    # RR ????
    rR = [((2/200) * (773 - (i + 473)) + 5) if i > 200  
            else ((8/200) * (773 - (i + 473)) + 2) 
            for i in range(800)]


    def __init__(self, he_fit, dpa_fit, z, uf, omega = 1.14E-29, s = 1, efv = 1.6, rM = 5000, r = 380, e = 1.4, f = 1, teol = 56.25, _N0 = 6e14, fi = 0.01): 
        
        self.teol = teol                        # 56.25                         Time End-Of-Life
        self.Teol = teol * 365 * 24 * 3600      # 1.773.900.000                 Time End-Of-Life (seconds)
        self.he_fit = he_fit                    # [[]]                          Helium Fit
        self.dpa_fit = dpa_fit                  # [[]]                          DPA Fit
        self.z = z                              # 200 + 273 inicial?            T° en Kelvin
        self.uf = uf                            # 0.107                         FMD Rate: 0.107 H347 Mode / 0.108 X750 Mode
        self.omega = omega                      # 1.14E-29                      Fraccion atomica de vacancias
        self.s = s                              # 1                             Cavity Surface Energy
        self.efv = efv                          # 1.6 eV? MeV?                  Energia de formacion de vacancia
        self.rM = rM                            # 5000                          R ?????????
        self.r = r                              # 380                           Parametro de recombinacion
        self.e = e                              # 1.4 eV? MeV?                  Energia de migracion de vacancia
        self.f = f                              # 1                             Eficiencia de Cascada
        self.fi = fi                            # 0.01                          Fraccion Inicial
        self.N0 = _N0                           # 6e14                          Densidad de dislocaciones

        self.rv = None                          # Tabulate Object with result values      
        self.AGBS = []                          # Swelling
        self.YB = []                            # EOS Constante
        self.PB = []                            # Presion
        self.CGB = []                           # Vacancy flux to cavities? / Vacancy emission term
        self.CJVS = []                          # Swelling rate
        self.DPA = []                           # Displacement per atom
        self.HELIO = []                         # Helio
        self.RADIO = []                         # Cavity radius


    def rateHe(self, t:float) -> float:                        # Tasa de crecimiento del Helio
        """ 
        Calcula la tasa de formacion de Helio para un determinado tiempo de vida derivando la curva de generacion de Helio.\n
        """
        x = self.he_fit
        tot = 0
        for i in range(1, len(x)):
            tot += i * x[i] * self.td(t)**(i - 1)           
        return tot


    def he(self, t:float) -> float:                            # Helio Total
        """
        Recrea la curva de generacion de Helio utilizando los coeficientes de un polinomio de interpolacion.\n
        Devuelve el Helio total generado en un determinado tiempo de vida.
        """

        x = self.he_fit
        tot = 0
        for i in range(0, len(x)):
            tot += x[i] * self.td(t)**i                     
        return tot


    def dpa(self, t: float) -> float:                           # DPA Total
        """
        Recrea la curva de generacion de DPA utilizando los coeficientes de un polinomio de interpolacion.\n
        Devuelve el DPA para un determinado tiempo de vida.
        """

        x = self.dpa_fit
        tot = 0
        for i in range(0, len(x)):
            tot += x[i] * self.td(t)**i
        return tot


    def rateDpa(self, t:float) -> float:                       # Tasa de crecimiento del DPA
        """ 
        Calcula la tasa de crecimiento de DPA para un determinado tiempo de vida derivando la curva de DPA.\n
        """

        x = self.dpa_fit    
        tot = 0
        for i in range(1, len(x)):
            tot += i * x[i] * self.td(t)**(i - 1)
        return tot


    def td(self, t:float) -> float:
        """
        Devuelve la cantidad de anios representativa de un determinado % de tiempo de vida, ejemplo:
            t= 1.00 -> TEOL (56.25 anios)
        """
        return self.teol * t


    def GHe(self, f: float, t: float) -> float:                # Generacion de Helio
        """
        PROTECTED-COG Page 10, Ecuation (22)
        He production
        """
        # usar rate en casos no lineales
        return f * (self.rateHe(t)  / (24 * 365 * 3600 * 1E6))


    def heTot(self, f: float, t: float) -> float:              # Helio total
        """
        PROTECTED-COG ecuacion 37 (segundo termino, no considera DENUDED ZONE)
        Helio total que migra a borde de grano
        """
        # usar integracion en casos no lineales
        return (CavitySwelling.integrate(lambda tt: self.GHe(f,tt) , 0, t, 1000)) * self.Teol


    def CI(self, f:float, z:int, rM:float, d:float, ro:float, t:float, r:float, e:float) -> float:
        """
        PROTECTED-PODG  Ecuacion 28, pagina 11
        Steady state point defect concentration (intersticiales)
        """
        rdv = self.Rd(z, t, e)
        rcv = self.Rc(d, ro, t)
        ssgbv = self.ssgb(rM, d, ro, z, t, e)
        return (self.G(f, t) / (((rdv * (1 + self.ba(z, t, e))) + rcv + ssgbv) * self.DI(z))) * self.Q(f, z, rM, d, ro, t, r, e) 


    def CE(self, z:int, efv:float) -> float:
        """
        A Rate Theory Model of Radiation-Induced Swelling in an Austenitic Stainless Steel.

        Pagina 502 - Ecuacion (15)

        el factor 1.6 E-19 es para pasar el electronvolts a Joules
        """
        return exp(-(efv) * 1.6E-19 / (1.38E-23 * z))


    def C(self, f:float, z:int, rM:float, d:float, ro:float, t:float, r:float, e:float) -> float:
        """
        PROTECTED-COG ecuacion 27, pagina 11
        Steady-state point defect concentration (vacancias)
        """
        return (
            self.G(f, t)/
            ((self.Rd(z, t, e) + self.Rc(d, ro, t) + self.ssgb(rM, d, ro, z, t, e)) * self.DV(z, e)) *
            self.Q(f, z, rM, d, ro, t, r, e)
        )


    def G(self, f:float, t:float) -> float:
        """
        - A Rate Theory Model of Radiation-Induced Swelling in an Austenitic Stainless Steel. Pagina 499.
        - PROTECTED-COG- Page 21

        Freely-migrating point defect generation rate in atom fraction per second.

        FMD rate -> con los valores de 0.107 y 0.108 da un valor medio de 5e-9

        G == phi
        """
        return f * (self.rateDpa(t) / (24 * 365 * 3600))


    def Rd(self, z:int, t:float, e:float) -> float:  # verificar e usada
        """
        PROTECTED-COG - Pagina 9 Ecuacion 15
        
        Dislocation density

        - Podria ser definida como constante de valor 4E14 segun el articulo.

        """
        rRt = self.rR[int(z - 473)] * 1E14
        return (rRt  * 0.4) / (
            ((rRt - self.N0) / self.N0) * (e**(-self.rr(z) * t * 300)) + 1
        )


    def Rc(self, d: float, ro:float, t:float) -> float:
        """
        PROTECTED-COG - pagina 9 Ecuacion 17
                
        Sink strength of the cavities within the matrix as a function of irradiation dose
        """
        return self.ss(d, ro) * (self.dpa(t) / self.dpa(1)) #cambio 347H


    def ssgb(self, rM:float, d:float, ro:float, z:int, t:float, e:float) -> float:
        """
        PROTECTED-COG - Pagina 9 Ecuacion 18

        Grain boundary sink strength as a function of time

        DV: 0nm (no considera la denuded zone)

        rM: tamaño medio de radio de grano (5000nm) 5 micrones
        """
        # return (3 / (r * 1e-9)) * sqrt(ss(d, ro) * dpa(t))
        # redeclaracion
        return (3 / (rM * 1E-9)) * ((self.Rd(z, t, e) + self.Rc(d, ro, t))**0.5)


    def DV(self, z:int, e:float) -> float:
        """
        Coeficiente de difusion para vacancias.
        """
        return 6 * 1E-6 * exp(-e * 1.6E-19 / (1.38E-23 * z))


    def Q(self, f:float, z:int, rM:float, d:float, ro:float, t:float, r:float, e:float) -> float:
        """
        A Rate Theory Model of Radiation-Induced Swelling in an Austenitic Stainless Steel.
        Pagina 499 - Ecuacion (4)

        F(n) = 2/n * {(1+ n)^0.5 - 1}

        F(n) == Q
        """
        n = self.n(f, z, rM, d, ro, t, r, e)

        return (2/n) * ((1 + n)**0.5 - 1)


    def rr(self, z):
        """
        ??????


        """
        return 1 * (((z - 773)/z) + 1)


    def ss(self, d:float, ro:float) -> float:
        """
        PROTECTED-COG - Pagina 9 ecuacion 16
        
        Sink strength of the cavities within the matrix
        """
        return 4 * pi * ro * (d/2) * 1E-9 * 1E23


    def n(self, f:float, z:int, rM:float, d:float, ro:float, t:float, r:float, e:float) -> float:
        """
        A Rate Theory Model of Radiation-Induced Swelling in an Austenitic Stainless Steel.

        Pagina 499 - Ecuacion (5)

        n = 4.a.(phi) / Dv.Di. Sum Kv^2 . Sum Ki^2

        a == self.a(z, r)

        phi == self.G(f,t)
            
        """
        rdv = self.Rd(z, t, e)
        rcv = self.Rc(d, ro, t)
        ssgbv = self.ssgb(rM, d, ro, z, t ,e)

        return (4 * self.a(z, r) * self.G(f, t)) / (
            ((rdv * (1 + self.ba(z, t,e)) + rcv + ssgbv) * (rdv + rcv + ssgbv) * (self.DV(z, e) * self.DI(z)))
        )


    def a(self, z:int, r:int) -> float:          # Parametro de recombinacion
        """
        PROTECTED-COG Pagina 10 Ecuacion 20

        Parametro de tasa de recombinacion.

            z: temperatura en grados kelvin.

            r: coeficiente de recombinacion.

            DI: Coeficiente de difusion para atomos intersticiales.

        Se utiliza un tamanio de grano de 31nm.

        1E19 == 1/a_0^2

        a_0: lattice parameter
        """
        
        # a = (1/a_0^2) * r * DI    donde a_0 = 36 nm = 36E-9
        return 1E19 * r * self.DI(z)


    def ba(self, z:int, t:float, e:float) -> float:
        """
        - A Rate Theory Model of Radiation-Induced Swelling in an Austenitic Stainless Steel - Pagina 500 - Ecuacion (9)
        - Heald and Speight Ecuacion 24

        ba == b
        """
        return (self.zIa(z, t, e) - self.zVa(z, t, e)) / self.zIa(z, t, e)


    def DI(self, z:int) -> float:
        """

        Coeficiente de difusion atomos intersticiales

        """
        return 12E-6 * exp((-0.15 * 1.6E-19)/(1.38E-23 * z))


    def zIa(self, z:int, t:float, e:float) -> float:
        """
        - A Rate Theory Model of Radiation-Induced Swelling in an Austenitic Stainless Steel. Pagina 500 - Ecuacion (10)
        - Heald and Speight - pagina 1392 ecuacion 23

        R == ((pi * self.Rd(z, t, e))**-0.5: mean distance between dislocations.

        li == self.lIa(z): radius of an effective trapping cilynder arround the dislocation core.
        """

        return 2 * pi / log(2 * (((pi * self.Rd(z, t, e))**-0.5)/self.lIa(z)))


    def zVa(self, z:int, t:float, e:float) -> float:
        """
        - A Rate Theory Model of Radiation-Induced Swelling in an Austenitic Stainless Steel. Pagina 500 - Ecuacion (10)
        - Heald and Speight - pagina 1392 ecuacion 23

        R == ((pi * self.Rd(z, t, e))**-0.5: mean distance between dislocations.

        lv == self.lVa(z): radius of an effective trapping cilynder arround the dislocation core.
        """

        return 2 * pi / log(2 * (((pi * self.Rd(z, t, e))**-0.5)/self.lVa(z)))


    def lIa(self, z):
        """
        - A Rate Theory Model of Radiation-Induced Swelling in an Austenitic Stainless Steel. Pagina 500
        - Heald and Speight - Pagina 1392 Definicion entre ecuaciones (21 - 22)

        Radius of an effective trapping cilynder arround the dislocation core. (intersticial)
        """

        return (573/z) * 3.6E-9


    def lVa(self, z):
        """
        - A Rate Theory Model of Radiation-Induced Swelling in an Austenitic Stainless Steel. Pagina 500
        - Heald and Speight - Pagina 1392 Definicion entre ecuaciones (21 - 22)
        
        Radius of an effective trapping cilynder arround the dislocation core. (vacancies)
        """

        return (573/z) * 6E-10


    def rho1(self ,z:int) -> float:
        """
        Cavity Number Density
        """
        return ((2.22662E31 * ((z - 273.0)**-11.5142))/1000.0) * 5.0  # Curvas de Bachatarya

    #----------------------- HERRAMIENTAS DE CALCULO Y GRAFICACION --------------------------


    def transpose(matrix: list) -> list:
        """ Devuelve una matriz traspuesta"""
        return [list(i) for i in zip(*matrix)]


    def integrate(f, a:float, b:float, n:int) -> float:
        """ Devuelve el valor de la integral de una funcion, dado los puntos y un paso determinado."""
        h = (b - a) / n
        s = 0
        for i in range(n):
            s += f(a + i * h)
        return s * h

    def plot_graph(x:list, y:list, x_label:str, y_label:str, title:str, t:int, block=True, fn:int = 0) -> None:
        """ Crea un grafico de Y vs X en la carpeta de origen. """
        plt.figure(fn)
        plt.plot(x, y)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title +" "+ str(t-273)+"°c")
        # plt.show(block=block)
        plt.savefig(title + '.png')

        return None

    #plot_graph([i + 573 for i in range(800)], rR, "i+573", "RRi", "J1")

    def run(self, silent:bool=False) -> None:

        """ Inicializa el algoritmo para el calculo de Swelling, todos los resultados se guardan como atributos dentro del objeto CavitySwelling """

        omega = self.omega              # Fraccion atomica de vacancias
        z = self.z                      # Temperatura en K
        rho1z = self.rho1(z)            # ??????
        s = self.s                      # Sigma Surface Energy == Cavity surface energy ??
        efv = self.efv                  # energia de formacion de vacancia
        rM = self.rM                    # Radio de grano 5um => 5000nm
        r = self.r                      # Parametro de recombinacion
        e = self.e                      # Energia de migracion de vacancias
        f = self.f                      # Eficiencia de Cascada
        uf = self.uf                    # FMD Rate: 0.107 H347 Mode / 0.108 X750 Mode
        PI = pi

        # Declaracion de variables:
        CGB = None                      # Vacancy emission term # chequearlo
        YB  = None                      # Parametro de ecuacion de estado (EOS parameter)
        AGB = None                      # Grain Boundary cavity fractional area coverage
        CJV = None                      # Incremental Swelling
        PB  = None                      # Pressure in Cavities (PASCAL)
        #inicializacion de parametros:
        AGB = self.heTot(f, self.fi)

        final = 100
        rango = final + 2

        #inicializacion de listas:
        self.AGBS = [0 for i in range(rango)]
        self.YB = [0 for i in range(rango)]                    
        self.PB = [0 for i in range(rango)]                
        self.CGB = [0 for i in range(rango)]                 
        self.CJVS = [0 for i in range(rango)]            
        self.DPA = [0 for i in range(rango)]                
        self.HELIO = [0 for i in range(rango)]               
        self.RADIO = [0 for i in range(rango)]
        self.SS = [0 for i in range(rango)]              

        self.AGBS[0] = AGB

        #Funciones y constantes que se pueden calcucar fuera del ciclo For (OPTIMIZACIONES)
        C1 = 6 / (PI * rho1z * 1E-4)
        C2 = (2E-10)**3 / (6 * omega)
        C3 = omega * 6.023E23
        C4 = (efv-0.1) * (1.6E-19)
        C5 = 1.38 * 1E-23
        C6 = self.Teol/100

        DV = self.DV(z, e)
        CE = self.CE(z, efv)
        DI = self.DI(z)

        for i in range(1, rango):

            hefi5 = self.heTot(f, i/100)

            vTerm = (AGB * C1)**(1/3)        # incremento de radio

            radio = (vTerm/2)*1E-9

            YB = PI * C2 * (hefi5 / AGB)

            PB =  ((1 + YB + YB**2 - YB**3) / ((1 - YB) ** 3)) * ((hefi5 /(C3 * AGB)) * 8.31 * z)

            if(PB < 0):
                PB = 0

            CGB = exp(-( C4 + ((PB - ((2*s) / ((vTerm/2) * 1E-9))) * omega)) / ( z * C5))
            """ CGB == Cv emit """

            #CJV = (self.Rc(vTerm, rho1z, i/100) * self.DV(z, e) *
            #(self.C(uf, z, rM, vTerm, rho1z, i/100, r, e) + self.CE(z, efv) - CGB)) + ((-(self.Rc(vTerm, rho1z, i/100))) * self.DI(z) * self.CI(uf, z, rM, vTerm, rho1z, i/100, r, e))
            
            CJV = (self.Rc(vTerm, rho1z, i/100) * DV *
            (self.C(uf, z, rM, vTerm, rho1z, i/100, r, e) + CE - CGB)) + ((-(self.Rc(vTerm, rho1z, i/100))) * DI * self.CI(uf, z, rM, vTerm, rho1z, i/100, r, e))
            
            AGB = AGB + ((CJV * C6) if CJV > 0 else 0)
            #print('AGB: ', AGB)
            
            self.AGBS[i] = AGB
            self.YB[i] = YB
            self.PB[i] = PB
            self.CGB[i] = CGB
            self.CJVS[i] = CJV
            self.DPA[i] = self.dpa(i/100)
            self.HELIO[i] = hefi5
            self.RADIO[i] = radio
            self.SS[i] = self.ss(vTerm, rho1z)

        SWELL = [i*100 for i in self.AGBS]                       #SWELL = [i*100 for i in self.AGBS]
        
        HE_SWELL = [i * 1E-4 for i in self.HELIO]                   #HE_SWELL = [i/10000 for i in self.HELIO] # HE * E-4 = [%]
        
        header = ["i", "Vol","AGBS [%]" , "YB", "PB", "CGB", "CJV", "DPA", "HELIO", "HE [%]", "RADIO", "SS"]
        
        self.rv = tabulate.tabulate(CavitySwelling.transpose([range(0,rango), self.AGBS , SWELL , self.YB , self.PB , self.CGB, self.CJVS , self.DPA , self.HELIO , HE_SWELL, self.RADIO, self.SS]), headers=header)
        
        if not silent:      #si silent ==False ejecuta este codigo.
            try:
                #CavitySwelling.plot_graph([i/100  for i in range(0, final + 2)], AGBS, "Time (s)", "Volume", "Cavity Swelling", z, True, 1)  
                pass

            except Exception as e:
                #print('error ({}), al tratar de añadir un nuevo grafico'.format(e))
                pass
                
        self.deol =(AGB*6/pi/rho1z/0.0001)**0.333
        
        return None


if __name__ == "__main__":
    import os
    os.system('python run_mtsf.py')
    #raise Exception('Este codigo no esta hecho para ser ejecutado. Se debe ejecutar "run_mtsf"')
