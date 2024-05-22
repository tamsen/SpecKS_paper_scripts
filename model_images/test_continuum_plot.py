import os
import unittest
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib import pyplot as plt
import numpy as np

class TestContinuumPlot(unittest.TestCase):

    #https://stackoverflow.com/questions/24976471/matplotlib-rectangle-with-color-gradient-fill
    #https://matplotlib.org/stable/gallery/images_contours_and_fields/image_transparency_blend.html
    #https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html

    def test_continuum_plot_2(self):

        out_folder= "/home/tamsen/Data/Specks_outout_from_mesx"

        xmin, xmax, ymin, ymax = (0, 100, 0, 100)

        Xred = [[[2*i-2*j, 0,0,3*i] for j in range(0,xmax)] for i in range(0,ymax)]
        Xblue = [[[0, 0,5*j,5*j] for j in range(0,xmax)] for i in range(0,ymax)]

        fig, ax = plt.subplots()
        ax.imshow(Xblue)
        ax.imshow(Xred)
        #ax.imshow(Xblue)
        ax.set_axis_off()
        plt.tight_layout()
        out_file_name=os.path.join(out_folder, "continuum_plot2.png")
        plt.savefig(out_file_name)
        plt.close()

def normal_pdf(x, mean, var):
    return np.exp(-(x - mean)**2 / (2*var))

def linear(x, m, b):
    return m*x-b
def exp(x, mean, var):
    return np.exp(-(x - mean)/ (2*var))

if __name__ == '__main__':
    unittest.main()
