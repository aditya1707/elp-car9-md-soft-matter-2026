import pandas as pd
import seaborn as sns
import plumed as pl
import matplotlib.pyplot as plt

def plot_kde(df):
    f, ax = plt.subplots(1,1,figsize=(4,4))

    sns.kdeplot(df[df.columns[1:]].values.flatten(),label="Car9", ax=ax, color='tab:orange')

    ax.legend(fontsize=14)
    ax.set_xticks([0,1,2,3,4,5,6])
    ax.set_xticklabels([str(l) for l in [0,1,2,3,4,5,6]], fontsize=14)
    ax.set_xlabel('End-to-end distance (nm)', fontsize=14)
    ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    ax.set_yticklabels([str(l) for l in [0,0.2,0.4,0.6,0.8,1.0]], fontsize=14)
    ax.set_ylabel('Probability density', fontsize=14)
    ax.set_title('Car9 Trial 1', fontsize=16)
    ax.grid(True)
    ax.set_xlim(0,6.5)
    plt.tight_layout()
    plt.savefig("figureS13_car9_kde.png", dpi=300)

def main():
    df = pl.read_as_pandas("car9_COLVAR")
    # Only do this for the last 100 ns of the simulation.
    df = df[-int(len(df[1:])/5):]
    plot_kde(df)

if __name__ == "__main__":
    main()