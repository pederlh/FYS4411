import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(font_scale=1.3, rc={'legend.facecolor': 'White', 'legend.framealpha': 0.5, 'lines.markersize':5})

# Gradient descent. eta = 0.1

# No interaction
h = np.array([1,2,3,4,5,6,7,8,9,10])
e_1 = np.array([2.00746737193835,2.0002131132685,2.01495378671377,2.03699160689562,2.0309648502258,2.03549592425661,2.00917283485631,2.00190383746058,2.02666419332231,2.02424991274717 ])
i_1 = np.array([15,15,12,10,169,219,326,447,411,418])

e_2 = np.array([2.00047693357054,2.01450017375426,2.0102719734061,2.0133530466665,2.0017033215601,2.02311872357203,1.99273193637732,2.01238340180815,2.01724670113913,2.05147439584335 ])
i_2 = np.array([12,11,12,177,9,16,256,369,314,413])

e_3 = np.array([2.00365469458715,2.01370187006535,2.00555595809565,2.00980672817809,2.00295057618727,2.0068859031555,1.99264989827676,2.0650441073482,2.00024427407747,2.03914022949817 ])
i_3 = np.array([13,13,8,132,173,168,236,278,445,436])

e_4 = np.array([2.00364869885244,2.01358982886694,2.00554343151456,2.00980657444476,2.00289361863933,2.00689604656289,1.99267154590314,2.06501725510505,2.0003629934268,2.03914876845424])
i_4 = np.array([15,14,12,9,238,8,240,394,342,480])

mean_it = (i_1 + i_2 + i_3 +i_4)/(4)
mean_e = np.mean(np.array([e_1, e_2, e_3, e_4]), axis=0 )


fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)


ax1.set_ylabel('Iterations')
ax1.plot(h, mean_it,"-o")

ax2.plot(h, mean_e, "-o", color="tab:orange")
ax2.set_ylim([1.98, 2.05])
ax2.set_ylabel(r'$\langle E_L \rangle$')
ax2.set_xlabel('Hidden nodes $N$')
ax2.axhline(2.0, ls="--", color="tab:red", alpha=0.8)

plt.tight_layout()
plt.savefig("hidden_layers.pdf")
plt.show()
