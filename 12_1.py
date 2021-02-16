import numpy as np
import matplotlib.pyplot as plt
"""
mu(0)=E_F=(3/2)**(2/3)=1.31037
In den gewählten dimensionslosen Einheiten sollte T_F = E_F sein. (glaube ich)
"""

def f(mu, T, e):
    return np.sqrt(e)/(np.exp((e-mu)/T)+1)

def integral(mu, T):
    e = 0
    F = 0
    de = 1e-4
    Fde = 1
    while ((e<mu) or (Fde>1e-10)) and F<2:
        Fde = f(mu, T, e)*de
        F += Fde
        e += de
    print(e)
    return F



def get_new_mu(oldmu, T):
    dmu = oldmu/10
    newmu = oldmu
    one = 10
    while np.abs(one-1) > 1e-5:
        one = integral(newmu, T)
        #print(one, newmu)
        if one > 1:
            newmu -= dmu
        if one < 1:
            newmu += dmu
        if one == 1:
            return newmu
        dmu *= 1/2
    return newmu

"""
diese Grenze für T weil sich der Integrator bei 1.25 aufhängt.
"""

"""
mu0 = (3/2)**(2/3)
mus = [mu0]
Ts = [0]
T = 1e-2


while T < 1.24:
    mus.append(get_new_mu(mus[-1], T))
    
    Ts.append(T)
    print(Ts[-1], mus[-1])
    T += 1e-2
w = open("12_1.txt", "w")
w.write(str(Ts))
w.write(str(mus))
w.close()
plt.plot(Ts, mus)
plt.show()
"""
def integral_2(mu, T):
    e = 0
    F = 0
    de = 1e-4
    Fde = 1
    while ((e<mu) or (Fde>1e-10)):
        Fde = e*f(mu, T, e)*de
        F += Fde
        e += de
    print(e)
    return F


mus = [1.3103706971044482, 1.310338705632351, 1.3101307661404902, 1.3097789243820055, 1.309379211575492, 1.3087718335232474, 1.3081168085577386, 1.3072545245286602, 1.3063449358092003, 1.3052525941575528, 1.3040097999394826, 1.3026965576482836, 1.3012097152637605, 1.2996054405904007, 1.297955550870901, 1.2959354174146505, 1.2938867890387327, 1.2916992485079482, 1.2893577282637654, 1.286847321065791, 1.2841847158299429, 1.281394373063613, 1.277890560324767, 1.2752230885740894, 1.2718996006662457, 1.2684139871365876, 1.2647443812411854, 1.2609001264691024, 1.2568597665667813, 1.2526482430569839, 1.2484890594374585, 1.2437035911188135, 1.2389136926447166, 1.2339758835547694, 1.2288318068429294, 1.2235066690080993, 1.2178909645702694, 1.2123233355879508, 1.2064481846574282, 1.2003732349971352, 1.194109080236866, 1.187680831892329, 1.1810552198296487, 1.1742431032357825, 1.1671907642856842, 1.1600383013502422, 1.152688937600233, 1.1452032135425751, 1.1374026350128297, 1.1295718844495086, 1.1214227513475838, 1.1131818273203469, 1.1047161012654083, 1.096098992028243, 1.0873015666466688, 1.0783557505343653, 1.0692465930713553, 1.059888074916666, 1.0504044572541427, 1.0407492502679079, 1.0309287037535386, 1.0209239635254619, 1.0107919910847127, 1.0004644410781465, 0.9899737311796342, 0.9792425706170425, 0.967767071742624, 0.9577491704140384, 0.9464787822270372, 0.9351792674685161, 0.923740673401091, 0.9121375343272417, 0.9005576632859779, 0.8884981936006175, 0.8764809084625623, 0.8641981926066661, 0.8518871661079848, 0.8393979511849418, 0.8267434532636357, 0.8139567838875977, 0.8010300532255158, 0.787956601307028, 0.7748752905431417, 0.7613868608029741, 0.7478636932923746, 0.7342338044485575, 0.7204579578001514, 0.7065975849791916, 0.692469083463128, 0.6782849633209023, 0.6642423449396493, 0.6495498751028493, 0.6349365887281598, 0.6202335283177922, 0.6053636717550952, 0.5903773738263923, 0.5753152557685656, 0.5600755816350476, 0.5447336675092046, 0.5292401440070508, 0.5135800108396547, 0.4978942541609337, 0.4820919462896149, 0.46599082855220786, 0.4507005044903385, 0.43377722871137997, 0.41746822157721186, 0.40096721936203833, 0.3844234058678526, 0.3677362765799367, 0.35085775607284975, 0.333931610418554, 0.31684370379166704, 0.2997947896520959, 0.2825653722975785, 0.26520857355000255, 0.24775246236126208, 0.2301025188756624, 0.21244035287602658, 0.19455719035853297, 0.17666894868640806, 0.15862249162332379, 0.1404056898509577, 0.12226538441415917]
Ts = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.060000000000000005, 0.07, 0.08, 0.09, 0.09999999999999999, 0.10999999999999999, 0.11999999999999998, 0.12999999999999998, 0.13999999999999999, 0.15, 0.16, 0.17, 0.18000000000000002, 0.19000000000000003, 0.20000000000000004, 0.21000000000000005, 0.22000000000000006, 0.23000000000000007, 0.24000000000000007, 0.25000000000000006, 0.26000000000000006, 0.2700000000000001, 0.2800000000000001, 0.2900000000000001, 0.3000000000000001, 0.3100000000000001, 0.3200000000000001, 0.3300000000000001, 0.34000000000000014, 0.35000000000000014, 0.36000000000000015, 0.37000000000000016, 0.38000000000000017, 0.3900000000000002, 0.4000000000000002, 0.4100000000000002, 0.4200000000000002, 0.4300000000000002, 0.4400000000000002, 0.45000000000000023, 0.46000000000000024, 0.47000000000000025, 0.48000000000000026, 0.49000000000000027, 0.5000000000000002, 0.5100000000000002, 0.5200000000000002, 0.5300000000000002, 0.5400000000000003, 0.5500000000000003, 0.5600000000000003, 0.5700000000000003, 0.5800000000000003, 0.5900000000000003, 0.6000000000000003, 0.6100000000000003, 0.6200000000000003, 0.6300000000000003, 0.6400000000000003, 0.6500000000000004, 0.6600000000000004, 0.6700000000000004, 0.6800000000000004, 0.6900000000000004, 0.7000000000000004, 0.7100000000000004, 0.7200000000000004, 0.7300000000000004, 0.7400000000000004, 0.7500000000000004, 0.7600000000000005, 0.7700000000000005, 0.7800000000000005, 0.7900000000000005, 0.8000000000000005, 0.8100000000000005, 0.8200000000000005, 0.8300000000000005, 0.8400000000000005, 0.8500000000000005, 0.8600000000000005, 0.8700000000000006, 0.8800000000000006, 0.8900000000000006, 0.9000000000000006, 0.9100000000000006, 0.9200000000000006, 0.9300000000000006, 0.9400000000000006, 0.9500000000000006, 0.9600000000000006, 0.9700000000000006, 0.9800000000000006, 0.9900000000000007, 1.0000000000000007, 1.0100000000000007, 1.0200000000000007, 1.0300000000000007, 1.0400000000000007, 1.0500000000000007, 1.0600000000000007, 1.0700000000000007, 1.0800000000000007, 1.0900000000000007, 1.1000000000000008, 1.1100000000000008, 1.1200000000000008, 1.1300000000000008, 1.1400000000000008, 1.1500000000000008, 1.1600000000000008, 1.1700000000000008, 1.1800000000000008, 1.1900000000000008, 1.2000000000000008, 1.2100000000000009, 1.2200000000000009, 1.2300000000000009]
plt.plot(Ts, mus)
plt.xlabel(r"$T$")
plt.ylabel(r"$\mu(T)$")
plt.show()

Es = []
i = 1
while i < len(Ts):
    Es.append(integral_2(mus[i], Ts[i]))
    i += 1

w = open("12_2.txt", "w")
w.write(str(Ts))
w.write(str(Es))
w.close()



plt.plot(Ts[1:], Es)
plt.xlabel(r"$T$")
plt.ylabel(r"$E(T)$")
plt.show()








