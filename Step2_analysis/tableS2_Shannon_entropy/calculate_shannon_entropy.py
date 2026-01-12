import numpy as np

def shannon_entropy(cluster_sizes, normalize=False):
    cluster_sizes = np.array(cluster_sizes)
    probabilities = cluster_sizes / np.sum(cluster_sizes)
    probabilities = probabilities[probabilities > 0]
    # print(len(probabilities))
    H = -np.sum(probabilities * np.log2(probabilities))
    # print(H)
    if normalize:
        if len(probabilities) <=1:
            return 0
        else:
            H /= np.log2(len(probabilities))
    return H

def main():

    sz = np.loadtxt("../figureS14_clustering/car9_cluster_sizes.xvg",comments=['#','@'])[:,1]
    h_s = shannon_entropy(sz,normalize=False)
    print(f"Shannon entropy: {h_s}")

if __name__ == "__main__":
    main()