import matplotlib.pyplot as plt
import numpy as np

def plot(com, z, rho_z):
    f, ax = plt.subplots(figsize=(7.85, 7.85))
    ax.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
    plt.plot(rho_z, z, c='blue', linewidth=4, label='Density of the Tag')
    plt.hlines(com, 0, 300, linestyles='dashed', colors='black',label='ELP Center of Mass', linewidth=4)

    ax.set_xlabel(r'$\rho(z)$ $Kg/m^3)$', fontsize=25, labelpad=15)

    plt.xticks([0,50,100,150,200,250,300], fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylabel(r'$z$ (nm)', fontsize=25)

    plt.tight_layout()

    plt.savefig('Car9_density.png',dpi=330)

def main():
    # Z-coordinate of the center-of-mass of the ELP domains at the end of the simulation. 
    com = np.loadtxt('com_ELP.xvg', comments='@', skiprows=13, unpack=True)[-1][-1]
    
    # Density profile of Car9 along the z axis.
    z, rho_z = np.loadtxt(f'./density.xvg', comments='@', skiprows=13, unpack=True)

    plot(com, z, rho_z)

if __name__ == "__main__":
    main()