import numpy as np
import matplotlib.pyplot as plt

fz = 14

omegas = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9])
analytical = omegas/2
energies = np.array([0,0.0500026480288108,0.0999988941902927,0.150000566401733,0.200007866265813,0.250000005735051,0.300031017716087,0.350001944342565,0.399992439128166,0.449974371858181,0.499999949866372,1.00000019270392,1.50002678004418,1.9999986910359,2.5000342404863,3.00007476053421,3.50003404450369,3.99993764481104,4.50008763297982])
rel_diff = np.abs((analytical[1:]-energies[1:])/analytical[1:])


plt.plot(omegas, analytical, label="Analytical solution")
plt.plot(omegas,energies, label="Numerical hiddenL = 2")
plt.title("$E_L$ vs. HO-frequency (1 particle, 1D)",fontsize = fz)
plt.xlabel("$\omega$",fontsize = fz)
plt.ylabel("$E_L$",fontsize = fz)
plt.legend()
plt.show()

plt.title("Relative error in numerical solution",fontsize = fz)
plt.plot(omegas[1:],rel_diff)
plt.xlabel("$\omega$", fontsize = fz)
plt.ylabel("Relative error", fontsize = fz)
plt.show()

hidden = np.array([1,2,3,4,5,6,7,8,9,10])
iters_r1 = np.array([5,12,42,18,18,2,8,4,9,7])
iters_r2 = np.array([9,30,42,8,3,2,18,23,1,10])
iters_r3 = np.array([8,10,13,0,36,1,7,1,35,23])
iters_r4 = np.array([88,3,4,6,22,17,50,5,5,5])

ener_r1 = np.array([3.0062527120414,3.00281746630548,3.01757693229415,3.01338919624095,3.01923730404959,3.00713848668016,3.00857283516593,3.00398101659322,2.99998934319632,2.99997382851243])
ener_r2 = np.array([2.9993931927932,3.00788751892436,3.00270027586509,3.01342620914303,3.02067792290327,3.01025931400546,3.00598446254538,2.99787513773223,3.0161075263796,2.99806959217254 ])
ener_r3 = np.array([3.01021273923847,3.0212124686853,3.00220894121847,3.00804653066758,3.01466829887061,3.02237800353325,3.01239784581274,3.00891824889974,3.01192817385259,3.00977545023747 ])
ener_r4 = np.array([3.00673388584726,3.00772557521105,3.00827164054594,3.00768326744341,3.01349283624169,3.00987095018036,3.00846847742672,2.9997411534417,3.0027009141199,3.01004645784099])

ana_ener = 3.0
joined_ener = (ener_r1 + ener_r2 + ener_r3 + ener_r4)/4
joined_its = (iters_r1+iters_r2+ iters_r3 + iters_r4)/4
rel_3d = np.abs((ana_ener - joined_ener)/ana_ener)

plt.title("Hidden L vs. Iterations")
plt.plot(hidden,joined_its)
plt.xlabel("H")
plt.ylabel("I")
plt.show()

plt.title("Hidden L vs. $E_L$")
plt.plot(hidden,joined_ener)
plt.xlabel("H")
plt.ylabel("$E_L$")
plt.show()

plt.title("Rel error")
plt.plot(hidden,rel_3d)
plt.xlabel("H")
plt.ylabel("$E_L$")
plt.show()
