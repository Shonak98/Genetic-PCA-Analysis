import os
import subprocess as sp
import sys
from matplotlib import pyplot as plt


def plot(eig_path):
    x_val = []
    y_val = []
    with open(eig_path) as eigen:
        for i in eigen:
            spl = i.split()
    #         if spl[0][:2]== 'HG':
            x_val.append(float(spl[2]))
            y_val.append(float(spl[3]))

    plt.scatter(x_val,y_val)
    plt.save_fig('pca.png')


if __name__ ==  '__main__':
    arguments = sys.argv
    if arguments[1] == 'data-test':
        os.system('chmod 700 process.sh')
        os.system('./process.sh')

        plot('data/eig.eigenvec')
