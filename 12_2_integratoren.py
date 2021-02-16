import matplotlib.pyplot as plt
import numpy as np

def f(e, z, b):
    return np.sqrt(e) * z*np.exp(-b*e) / (1-z*np.exp(-b*e))


def integral(z, b):
    e = 0
    F = 0
    de = 1e-4
    Fde = 1
    while ((e<5) or (Fde>1e-10)) and F<2:
        Fde = f(e, z, b)*de
        #print(Fde)
        F += Fde
        e += de
    print(e)
    return F



def get_new_z(b):
    dz = 0.25
    newz = 0.5
    one = 10
    while np.abs(one-1) > 1e-5:
        one = integral(newz, b)
        print(one, newz)
        if one > 1:
            newz -= dz
        if one < 1:
            newz += dz
        if one == 1:
            return newz, one
        if newz > 0.99999999:
            return newz, one
        dz *= 1/2
    return newz, one


def integral_2(z, b):
    e = 0
    F = 0
    de = 1e-4
    Fde = 1
    while ((e<5) or (Fde>1e-10)):
        Fde = e*f(e, z, b)*de
        print(Fde)
        F += Fde
        e += de
    print(e)
    return F


"""
Temperaturabhängigkeit der Fugazität z:
"""


z0 = 1e-4
zs = []
bs = []
ns = []
b = 2


while b < 5:
    zn = get_new_z(b)
    zs.append(zn[0])
    ns.append(zn[1])
    bs.append(b)
    print(bs[-1], zs[-1])
    b += 0.05

"""
plt.plot(bs, zs)
plt.xlabel(r"Inverse Temperatur $\beta$")
plt.ylabel(r"Fugazität $z$")
plt.show()
"""
np.savetxt("z_von_T_3.txt", np.array((bs, zs, ns)))


plt.plot(bs, ns)
plt.xlabel(r"$\beta$")
plt.ylabel(r"$n$")
plt.show()


a = np.loadtxt("z_von_T_3.txt")
bs = a[0]
Ts = 1/bs
zs = a[1]
ns = a[2]

"""

plt.plot(Ts, zs)
plt.xlabel(r"Temperatur $T=1/\beta$")
plt.ylabel(r"Fugazität $z$")
plt.show()

mus = np.log(zs)*Ts
print(mus)
plt.plot(Ts, mus)
plt.xlabel(r"$T$")
plt.ylabel(r"$\mu$")
plt.show()

"""




"""



i = 0
Es = []
while i < len(bs):
    Es.append(integral_2(zs[i], bs[i]))
    i += 1
print(Es)
np.savetxt("E_von_T.txt", np.array((Es)))

"""
"""
Es = np.loadtxt("E_von_T_3.txt")
print(Es)
"""



"""
#plt.plot(bs, Es, label=r"$E(\beta)$")
plt.plot(bs, ns, label = r"$n(\beta)$")
plt.plot(bs, zs, label=r"$z(\beta)$")



plt.xlabel(r"$\beta$")
plt.legend()
plt.show()



#plt.plot(Ts, Es, label=r"$E(T)$")
plt.plot(Ts, ns, label = r"$n(T)$")
plt.plot(Ts, zs, label=r"$z(T)$")



plt.xlabel(r"$T$")
plt.legend()
plt.show()


plt.plot(Ts, Es, label=r"$E(T)$")



plt.xlabel(r"$T$")
plt.legend()
plt.show()
"""



print(len(Es))
print(Es[38])
C = []
i = len(Ts)-1
while i > 0:
    C.append((Es[i]-Es[i-1])/(Ts[i]-Ts[i-1]))
    i -= 1
plt.plot(Ts[:-1], C)
plt.xlabel(r"$T$")
plt.ylabel(r"$C_V(T)$")
plt.show()










