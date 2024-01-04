# Sweeling-Atucha

Se modificaron las siguientes funciones en el codigo:
```python
# Antes:
def Rd(self, z, t, e):
    # e: energia de migracion de vacancias
    rRt = self.rR[int(z - 473)] * 1E14
    return (rRt  * 0.4) / (((rRt - self.N0) / self.N0) * (e**(-self.rr(z) * t * 300)) + 1)
 
# Despues:
def Rd(self, z, t, e):
    # e: numero neperiano
    rRt = self.rR[int(z - 473)] * 1E14
    e = 2.71
    return (rRt  * 0.4) / (((rRt - self.N0) / self.N0) * (e**(-self.rr(z) * t * 300)) + 1)
```
# Energia de formacion de vacancias
Se ajusto la energia de formacion de vacancias para que disminuyera a medida que aumenta la temperatura

![efv vs T](https://github.com/rmeza-gi/VS/assets/59147377/1af7ffae-a75d-465a-9039-6f2ee1acf283)

De esta manera se corrige la divergencia en el swelling a altas temperaturas:
# Curva Swelling vs Temperatura

Comparacion Swelling a efv cte             |  Swelling vs T- efv(T)
:-------------------------:|:-------------------------:
![SwellingVsTDivergente](https://github.com/rmeza-gi/VS/assets/59147377/54a155e6-8898-4259-ba91-4e74ba9efebf)  |  ![SwellingVsT](https://github.com/rmeza-gi/VS/assets/59147377/061a2c27-3524-4a85-984f-48cb60e6db17)


# Geometria de la pieza

Se simulo una pieza cilindrica de Inconel con radio=5u y altura = 10u 
(u: unidad de medida normalizada)

La pieza se particiono en volumenes discretos (32.000 partes)

El error (%) por la aproximacion en la particion del volumen fue de:  1.36e-12

Se utilizaron curvas parametricas para generar el siguiente gradiente de temperaturas. Las temperaturas de cada punto se calculo realizando una interpolacion bilineal entre las 4 curvas parametricas mas cercanas al punto.

Para el calculo de swelling de cada volumen, se utilizo la curva "Swelling vs T - efv(T)" a 56.25 app.

El codigo tambien admite la posibilidad de calcular swelling para un app distinto, y a su vez, admite que cada punto posea su propio dpa. Con esto se consigue simular la penetracion del daño por dpa en la pieza. Solo a modo de ejemplo al final de este archivo se mostrara la pieza simulando una penetracion arbitraria del dpa en la pieza. Para activar esta funcionalidad, activar la importacion de "modulo_swelling" y utilizar "self.vol_after2" para calcular el aumento de volumen de cada punto. (ESTE CALCULO PUEDE DEMORAR BASTANTE DEPENDIENDO DEL TAMAÑO DE LA PARTICION Y LA CAPACIDAD DE COMPUTO)

# Gradiente de temperatura

![GradienteT](https://github.com/rmeza-gi/VS/assets/59147377/45ddd6d5-e911-4430-a892-266d8c3601b9)

A partir de la curva de Swelling vs T, se calculo el aumento de volumen de cada particion y se obtuvieron los siguientes resultados:

volumen real del cilindro:  785.3981633974483

volumen calculado:  785.3981633974377

error en el calculo del volumen (%):  1.3606554285286765e-12


volumen con swelling:  870.962831194115

aumento de volumen total (%):  10.89443186708432

# Deformacion

En base al aumento de volumen de cada particion, se simulo la deformacion de la pieza
![Gradiente y Deformacion](https://github.com/rmeza-gi/VS/assets/59147377/75608ef9-775c-4017-81a6-63845133d6eb)


Bulon Inconel M12x60:

![Figure 2023-07-26 155007](https://github.com/rmeza-gi/VS/assets/59147377/d357078e-00f4-4011-a64b-9a1b70580d7b)


# Posibilidad de calcular swelling con dpa variable

Volumen atenuado en dpa             |  funcion arbitraria de atenuacion
:-------------------------:|:-------------------------:
![atenuacion](https://github.com/rmeza-gi/VS/assets/59147377/af295d96-9f24-4702-b487-27e81d80327f)  |  ![funcion atenuacion](https://github.com/rmeza-gi/VS/assets/59147377/de0ddce1-d7fa-4718-8d7d-8030a8f73f6d)


Con el DPA variable, se consiguen los siguientes valores de swelling:

volumen inicial:  785.3981633974531

volumen final:  852.3182627341536

aumento de volumen (%):  8.520531681309173
