from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count, freeze_support
from itertools import repeat
from Swelling_atucha_Voids import *

res = {}
def fun_(t, he, dpa, initialFraction, fmd_rate, omega, se, efv, rM, r, e, fr, teol, N0):
    cs = CavitySwelling(he, dpa, z = t + 273, fi = initialFraction, uf= fmd_rate,
                        omega=omega, s=se, efv=efv, rM=rM, r=r, e=e, f=fr, teol=teol, _N0=N0)
    b = cs.run(silent=False)
    return t, b[1][100]*100, cs.rho1(t + 273), cs.deol

def add(future):
    pass

def process(he, dpa, title, fi, fmd_rate, omega, se, efv, rM, r, e, fr, teol, N0):
    gbVal =  [i for i in range(200, 660, 10)]
    pool = Pool(processes=8)
    res = pool.starmap(fun_, zip(gbVal, repeat(he), repeat(dpa), repeat(fi), repeat(fmd_rate), repeat(omega), repeat(se), repeat(efv), repeat(rM), repeat(r), repeat(e),  repeat(fr),  repeat(teol),  repeat(N0)))
    pool.terminate()
    t = []; s = []; deol = []; rho1 = []
    for i in res:
        t.append(i[0])
        s.append(i[1])
        rho1.append(i[2])
        deol.append(i[3])

    import tabulate
    print(tabulate.tabulate(CavitySwelling.transpose([t,s,rho1,deol]), headers=["C", "%", "rho1","Deol"]))
    import matplotlib.pyplot as plt
    plt.clf() 
    plt.plot(t, s, marker='.', linestyle='None')
    plt.savefig(title + ".png")
    with open(title + ".txt", 'w') as ofile:
        ofile.write("Fraccion inincial = {:0.3f}\n".format(fi))
        ofile.write(tabulate.tabulate(CavitySwelling.transpose([t,s,rho1,deol]), headers=["C", "%", "rho1","Deol"]))
        ofile.flush()

def castAndFlip(strIn = "3 2 1 0"):
    return [float(i) for i in strIn.split()][::-1]

def readString(f):
    try:
        return f.readline().split("#")[0].split("=")[1].strip()
    except:
        return None

def readFloat(f):
    return float(readString(f))

def readFit(f):
    return castAndFlip(readString(f))

if __name__ == '__main__':
    freeze_support()
    outdir = "out_fiteos"
    index = 0
    import os
    try:
        os.mkdir(outdir)
    except:
        pass
    with open("datos.txt") as f:
        omega = readFloat(f)
        se = readFloat(f)
        efv = readFloat(f)
        rM = readFloat(f)
        r = readFloat(f)
        ee = readFloat(f)
        fr = readFloat(f)
        teol = readFloat(f)
        N0 = readFloat(f)
        modo = readString(f)
        fmd_rate =0.107 #H347 def
        if "X750" in modo:
            fmd_rate = 0.108
        while (line := readString(f)):
            index += 1
            start = timer()
            he = readFit(f)
            dpa = readFit(f)
            fi = 0.01
            runOk = False
            while (not runOk and fi < 0.1):
                try:
                    process(he, dpa, os.path.join(outdir, modo + "_" + line), fi, fmd_rate, omega, se, efv, rM, r, ee, fr, teol, N0)
                    runOk = True
                except Exception as e:
                    print("error ({}) en initialFraction {:0.3f} incrementando en 0.001".format(e, fi))
                    fi += 0.001
                    if fi >= 0.1:
                        print(line+ "_" + modo + ":" +  str(e))
                        print(he)
                        print(dpa)
            end = timer()
            time = end - start
            print(time)
