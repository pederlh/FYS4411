import numpy as np
import matplotlib.pyplot as plt

fz = 14
#std vekter = 0.1
#conv = 10^-3

omegas = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0])
analytical = omegas/2

energies_anasig = np.array([0.0504311505850437,0.103450555429314,0.150695409594281,0.200191463232913,0.250290530897608,0.299980668207348,0.350138040650666,0.400214933818619,0.450103340910881,0.500051696974224,0.550019064943937,0.599922039711554,0.650243163068026,0.700136514823112,0.750190186636847,0.800008802552174,0.849995975994293,0.899926895997044,0.949897380777095,0.999950473537914,1.04996744442788,1.09994883775279,1.15000103267144,1.1999961782411,1.24996061163836 ,1.29927046293778 ,1.34998589553193,1.39999355111919,1.4498435158713,1.50004652786108])

energies_noanasig = np.array([0.254471399107738,0.26041768811502 ,0.272792008499217,0.291349672638832,0.31431546131787 ,0.340269601477183,0.373857118759387,0.409937408095608,0.452360964511717,0.499948534311143,0.552593863511945,0.610194963879811,0.6717607140595,0.738377700976957,0.811170895867028,0.888528315807943,0.971393396613974,1.05730944266038 ,1.15724857534042 ,1.25565378146653 ,1.35091158387586 ,1.45935459435315 ,1.57699272656817 ,1.68989166723952 ,1.81713265862859 ,1.93155282674913 ,2.05734911801872 ,2.22723383718099,2.35309094242096 ,2.49564146698901])


rel_diff_ana = np.abs((analytical-energies_anasig)/analytical)
rel_diff_noana = np.abs((analytical-energies_noanasig)/analytical)


plt.plot(omegas, analytical, label="Analytical solution")
plt.plot(omegas,energies_anasig,"*" ,label="Numerical $\sigma = 1/\sqrt{\omega}$")
plt.plot(omegas, energies_noanasig,"*" ,label="Numerical $\sigma = 1.0$")
#plt.title("$E_L$ vs. HO-frequency (1 particle, 1D)",fontsize = fz)
plt.xticks(fontsize = fz-2)
plt.yticks(fontsize = fz-2)
plt.xlabel("$\omega$",fontsize = fz+2)
plt.ylabel("$E_L$",fontsize = fz+2)
plt.legend(fontsize = fz)
plt.show()

plt.title("Relative error in numerical solution",fontsize = fz)
plt.plot(omegas,rel_diff_ana, label = "$\sigma = 1/\sqrt{\omega}$")
plt.plot(omegas, rel_diff_noana, label = "$\sigma = 1.0$")
plt.xlabel("$\omega$", fontsize = fz)
plt.ylabel("Relative error", fontsize = fz)
plt.legend()
plt.show()
