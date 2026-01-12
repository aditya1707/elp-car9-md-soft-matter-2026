import numpy as np
import matplotlib.pyplot as plt

def plot_CDF():
    f, ax = plt.subplots(figsize=(6, 5))
    # c = [15,20,25]
    peptides = ['car9','p9ag10a']
    params = {
        'xlim':(0,400), 
        'xsmallticks':[0, 200, 400],
        'xbigticks':[0, 400, 800, 1200, 1600],
        }
    
    axins = ax.inset_axes([0.5,0.1,0.35,0.4],xlim=params['xlim'],ylim=(0.8,1.0), xticks=params['xsmallticks'],yticks=[0.8,0.85,0.90,0.95,1.0])

    for i in peptides:
        sz = np.loadtxt(f'{i}_cluster_sizes.xvg',comments=['#','@'])
        cluster_sizes = np.array(sorted(sz[:,1],reverse=True))
        c_cdf = np.cumsum(cluster_sizes) / np.sum(cluster_sizes)

        x = np.arange(1, len(cluster_sizes)+1)
        aa = ax.plot(x, c_cdf,linestyle='-', label=f'{i}',alpha=0.8)
        axins.plot(x, c_cdf,linestyle='-', label=f'{i}',alpha=0.8)

    ax.indicate_inset_zoom(axins, edgecolor="black")
    ax.set_xlabel("# of Clusters",fontsize=14)
    ax.set_xticks(params['xbigticks'])
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels([str(i) for i in params['xbigticks']],fontsize=14)
    ax.set_yticklabels([str(i) for i in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]],fontsize=14)
    ax.set_ylabel("Cumulative Fraction",fontsize=14)
    ax.set_title(f"CDF of Cluster Sizes (Cutoff = 0.15)",fontsize=14)
    ax.grid(True)

    ax.legend(loc='upper right',fontsize=14)#,bbox_to_anchor=[1.55, 0.8])
    plt.tight_layout()
    plt.savefig("cluster_size_CDF.png",dpi=330)

def main():
    plot_CDF()

if __name__ == "__main__":
    main()