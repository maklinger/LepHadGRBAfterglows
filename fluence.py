import matplotlib.pyplot as plt
import numpy as np



plt.figure(figsize=(6,7.5))
plt.loglog()
GeV2erg=1.0/624

data=np.loadtxt("sensitivities/IceCube")
plt.plot(data[:,0],data[:,1]*GeV2erg,"--",color="purple",alpha=1,linewidth=1.5)

data=np.loadtxt("sensitivities/IceCubeGen2")
plt.plot(data[:,0],data[:,1]*GeV2erg,"--",color="black",alpha=1,linewidth=1.5)



data=np.loadtxt("sensitivities/grand3yr")
plt.plot(data[:,0],data[:,1]*1e8*4*3.14*GeV2erg,"-.",color="green",alpha=1,linewidth=1.5)

data=np.loadtxt("sensitivities/chant3yr")
#plt.plot(data[:,0],data[:,1]*1e8*4*3.14*GeV2erg,":",color="gray",alpha=1,linewidth=1.5)

data=np.loadtxt("sensitivities/ara3yr")
#plt.plot(data[:,0],data[:,1]*1e8*3.14*4*GeV2erg,"-",color="gray",alpha=1,linewidth=1.5)


data=np.loadtxt("sensitivities/POEMMA")
plt.fill_between(data[:6,0],data[:6,1]*GeV2erg, data[6:,1]*GeV2erg,color="pink")

plt.text(1e5,6e-4,"IceCube", color = "purple")
plt.text(1e9,8e-6,"IceCube-Gen2", color = "k")
plt.text(3e9,4.5e-4, "GRAND", color = "green")
plt.text(6e8,3e-2,"POEMMA", color = "#F1948A")

#plt.text(8e9,4e-10,"GRAND 3yr",color="gray",rotation=40,fontsize=8)
#plt.text(1e9,7e-10,"CHANT 3yr",color="gray",rotation=40,fontsize=8)
#plt.text(5e7,4e-9,"ARA/ARIANNA 3yr",color="gray",rotation=-10,fontsize=8)
#plt.text(5e8,3e-8,"IceCube(2018)",color="black",rotation=20,fontsize=8)
#plt.text(7e9,1.2e-8,"POEMMA 10 yr",color="gray",rotation=40,fontsize=8)

Neu1= np.loadtxt("pgammacascade_neutrinos.txt")
plt.plot(Neu1[:,0]/1e9*50,Neu1[:,1]*1e3, linewidth=2.5, color = "red", label = r"$p\gamma$-cas. scenario")

Neu2 = np.loadtxt("ppcascade_neutrinos.txt")
plt.plot(Neu2[:,0]/1e9*23, Neu2[:,1]*1e3, linewidth =2.5, color = "deepskyblue", label = r"$pp$-cas. scenario")

Neu3 = np.loadtxt("proton_syn_neutrinos.txt")
plt.plot(Neu3[:,0]/1e9*50, Neu3[:,1]*1e3, linewidth = 2.5, color = "magenta", label = r"$p$-syn scenario")


plt.xlabel(r"$E_\nu$ [GeV]", fontsize=12)
plt.ylabel(r"Neutrino fluence [$\rm erg~cm^{-2}$]", fontsize=12)

plt.xticks(fontsize=12)
plt.yticks(10.0**np.arange(-10,1,1),fontsize = 12 )

plt.xlim(1e4,1e11)
plt.ylim(1e-10,1e-1)
plt.legend(loc=2)

plt.gca().set_aspect('equal')
plt.grid(color="gray", alpha = 0.1, which="minor")
plt.grid(color="gray", alpha = 0.3, which="major")
#plt.title("Neutrino fluence (1 ks)")
plt.tight_layout()

plt.savefig("All_Neu_fluence.pdf")
#plt.show()
