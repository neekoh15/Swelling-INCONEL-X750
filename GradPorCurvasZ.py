import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import random
import pandas as pd
from scipy.interpolate import interp1d
from math import cos, sin
from matplotlib import cm
import math
#from modulo_swelling import calculate_swelling

puntos = []

#GEOMETRIA DE LA PIEZA
r= 5
h= 10
paso_radial = 20
paso_altura = 20
paso_angular = 80
radios = [i *(r/paso_radial) for i in range(paso_radial+1)]

#CURVAS PARAMETRIZADAS DE TEMPERATURA ((r,z) VS T)
reference_curves = [
    (0, 0, interp1d([0,h], [400,650])),     # centro
    (0, 5, interp1d([0,h], [450,580])),     # 90 grados radio 5
    (5, 0, interp1d([0,h], [470,650])),     # 0 grados radio 5
    (-5, 0, interp1d([0,h], [280,620])),    # 180 grados radio 5
    (0, -5, interp1d([0,h], [450,620]))     # 270 grados radio 5
]

#SWELLING EN FUNCION DE LA TEMPERATURA:
S = []
T = []
with open('CS vs T (all fits)/Inc-1-raw.txt') as file:
    files = file.readlines()

    for i in range(1,len(files)):
        S.append(float(files[i].split()[1]))
        T.append(float(files[i].split()[0]))
        
    file.close()

swelling = interp1d(T,S)

#CALCULO DINAMICO DEL SWELLING (atenuacion de dpa en funcion de la distancia recorrida por el neutron)
def get_data():
    def readFloat(f:str) -> float:
        """ Lee una linea de texto proveniente del archivo de configuraciones y devuelve el valor asignado al parametro de esa linea"""
        
        try:
            return float(f.readline().split("#")[0].split("=")[1].strip())
        except:
            return None

    with open("datos.txt", 'r') as f:

            omega = readFloat(f)                # fraccion atomica de vacancias
            se = readFloat(f)                   # sigma surface energy
            efv = readFloat(f)                  # energia de formacion de vacancias
            rM = readFloat(f)                   # Radio del grano (um)
            r = readFloat(f)                    # parametro de recombinacion
            ee = readFloat(f)                   # energia de migracion de vacancias
            fr = readFloat(f)                   # eficiencia de cascada
            teol = readFloat(f)                 # tiempo de plena potencia en anios
            N0 = readFloat(f)                   # densidad de dislocacion

    f.close()

    data = pd.read_excel('MOD-TB-DPA&HEvsT.xlsx', engine='openpyxl')
    time = data['t'].tolist()
    data = data.set_index('t')
    material = '347-1'
    dpa = np.array(data[material + ' dpa'])
    he = np.array(data[material + ' He [appm]'])
    fmd_rate = 0.107 if material.split('-')[0] == '347' else 0.108

    return [omega, se, efv, rM, r, ee, fr, teol, N0, fmd_rate, dpa, he, time]

parametros = get_data()
app = 56.25


class Punto:
    global r,h,reference_curves, paso_altura, paso_radial, paso_angular, radios

    def __init__(self, x, y, z, r, phi) -> None:
        
        self.x = x
        self.y = y 
        self.z = z 
        self.r = r
        self.indice_en_radios = radios.index(self.r)

        self.rho = (x**2 + y**2 + z**2)**(1/2)
        self.phi = phi
        self.T = self.interpolate_temperature()

        self.area = math.radians(360/paso_angular)/2 * (radios[self.indice_en_radios]**2 - radios[self.indice_en_radios-1]**2)
        self.h = (h/paso_altura)
        self.vol = self.h*self.area

        #calculo de swelling para teol=56.25 y atenuacion = 0 (sin tener en cuenta el efecto de la atenuacion del danio en la penetracion)
        self.vol_after = self.swelling()

        #swelling teniendo en cuenta la penetracion neutronica
        #self.atenuacion = self.calculate_attenuation()
        #self.vol_after2 = self.vol * (1 + calculate_swelling(app, self.atenuacion, self.T, parametros))

    #Calcular la temperatura interpolando las 4 curvas mas cercanas, utilizando el metodo de interpolacion bilineal
    def interpolate_temperature(self):
            # Verificar si el punto coincide con alguno de los puntos de referencia
        for curva in reference_curves:
            if (self.x, self.y) == (curva[0], curva[1]):
                return curva[2](self.z)

            # Calcular la distancia desde el punto intermedio (x, y) a cada uno de los puntos de referencia
            distances = [math.sqrt((self.x - curva[0])**2 + (self.y - curva[1])**2) for curva in reference_curves]

            # Realizar la interpolación ponderada
            t_interpolated = sum([curva[2](self.z) / distance for curva, distance in zip(reference_curves, distances)]) / sum([1 / distance for distance in distances])

        return t_interpolated

        # Ejemplo de puntos y temperaturas de referencia (tuplas: (x, y, t))
    
    def swelling(self):
        
        return self.vol * (1 + swelling(self.T)/100)
        
    def calculate_attenuation(self):
        #calcula la atenuacion del dpa en base a la curva de penetracion

        fuente_x = 0
        fuente_y = 0
        fuente_z = 0
        d = ((self.x - fuente_x)**2 + (self.y - fuente_y)**2 + (self.z-fuente_z)**2)**0.5

        return (1 - self.r/r)
        
    def __repr__(self) -> str:
     
        return (f'(x: {self.x}, y: {self.y}, z: {self.z}, phi: {self.phi}), r: {self.r}, h: {self.h}, indice: {self.indice_en_radios}\nArea: {self.area}, volumen: {self.vol}\n')
    

class GenDistribucion:

    def __init__(self, radio, altura, pasos_angulares, pasos_radiales, pasos_altura) -> None:

        self.radio = radio
        self.altura = altura
        self.paso_angular = pasos_angulares
        self.paso_radial = pasos_radiales
        self.paso_h = pasos_altura

        self.construct_geometry()
        

    def construct_geometry(self):

        #por cada paso en altura ...
        for z in range(self.paso_h):
            #por cada paso en el radio ...
            for i in range(1,self.paso_radial+1):
                rad = (i+1)*(self.radio/self.paso_radial)
                #por cada paso de angulo ...
                for j in range(self.paso_angular):

                    delta_phi = 2*np.pi/self.paso_angular

                    angulo_actual = delta_phi*(j+1)*180/np.pi
                    

                    x= rad*cos(delta_phi*(j+1))
                    y= rad*sin(delta_phi*(j+1))
                    h= z *(self.altura/self.paso_h)

                    r = i * self.radio/self.paso_radial


                    p = Punto(x,y,h,r, angulo_actual)

                    puntos.append(p) #if 180>angulo_actual>0 else None


def plot_3d_points(puntos):
        # Extraer las coordenadas x, y, z del array
        x = [point.x for point in puntos]
        y = [point.y for point in puntos]
        z = [point.z for point in puntos]

        atenuacion = np.array([p.atenuacion for p in puntos])

        # Configuración de la figura 3D

        ax1 = plt.subplot(121, projection='3d')
        ax1.set_box_aspect((1, 1, 1))  # Igual escala en todos los ejes

        #cmap = cm.get_cmap('Blues')
        
        # Dibujo de los puntos en el gráfico 3D
        sc = ax1.scatter(x, y, z, c=atenuacion, cmap=matplotlib.colormaps['Greens_r'], marker='h')

        # Etiquetas y título

        cbar = plt.colorbar(sc)
        cbar.set_label('Atenuacion de la dosis (%)')

        plt.title('dosis vs distancia a la fuente')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_zlabel('Z')

        #plt.show()


def plot_3d_points_with_temperature(points):
    # Extraer las coordenadas x, y, z y temperatura del array de objetos
    x = np.array([point.x for point in points])
    y = np.array([point.y for point in points])
    z = np.array([point.z for point in points])
    temperatura = np.array([point.T for point in points])
    
    

    # Configuración de la figura 3D
    ax2 = plt.subplot(122, projection='3d')
    ax2.set_box_aspect((1, 1, 1))  # Igual escala en todos los ejes

    # Mapa de colores para la temperatura (de azul oscuro a rojo claro)
    cmap = matplotlib.colormaps['jet']

    # Normalizar los valores de la temperatura para que se ajusten al mapa de colores
    temperatura_normalizada =  (temperatura - min(temperatura)) / (max(temperatura) - min(temperatura))

    # Dibujo de los puntos en el gráfico 3D con colores representativos de la temperatura
    sc = ax2.scatter(x, y, z, c=temperatura, cmap=cmap, marker='h', alpha=1)

    # Barra de colores
    cbar = plt.colorbar(sc)
    cbar.set_label('Temperatura')

    # Etiquetas y título
    plt.title('Gradiente de temperatura')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')

    #plt.show()



distribucion = GenDistribucion(radio= r, altura= h, pasos_angulares=paso_angular, pasos_radiales=paso_radial, pasos_altura=paso_altura)

#plt.figure(figsize=(12,12))

#plot_3d_points(puntos)
#plot_3d_points_with_temperature(puntos)

#plt.show()

vol_cilindro = np.pi*(r**2)*h
vol_calculado = sum([p.vol for p in puntos])
error = abs((vol_cilindro - vol_calculado)/vol_cilindro)*100

vol_after_swelling = sum([p.vol_after for p in puntos])

print('volumen real del cilindro: ', vol_cilindro)
print('volumen calculado: ', vol_calculado)
print('error en el calculo del volumen (%): ', error)
print('volumen con swelling: ', vol_after_swelling)
print('aumento de volumen total (%): ', ((vol_after_swelling/vol_cilindro) -1)*100)

with open('puntos.txt', 'a') as p:
    for punto in puntos:
        p.write(f'x:{punto.x}, y:{punto.y}, z:{punto.z}, T:{punto.T}, V_before:{punto.vol}, V_after:{punto.vol_after}\n')
