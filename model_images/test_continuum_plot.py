import colorsys
import os
import unittest

import matplotlib
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib import pyplot as plt
import numpy as np

import two_d_colors


class TestContinuumPlot(unittest.TestCase):

    #https://stackoverflow.com/questions/24976471/matplotlib-rectangle-with-color-gradient-fill
    #https://matplotlib.org/stable/gallery/images_contours_and_fields/image_transparency_blend.html
    #https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html


    def test_continuum_plot_3(self):

        out_folder= "/home/tamsen/Data/Specks_outout_from_mesx"

        xmin, xmax, ymin, ymax = (0, 100, 0, 100)
        nice_blue = matplotlib.colors.to_rgb(two_d_colors.high_Ne)
        nice_orange = matplotlib.colors.to_rgb(two_d_colors.low_dT)
        #Xred = [[[*nice_orange,0*i] for j in range(0,xmax)] for i in range(0,ymax)]
        #Xblue = [[[*nice_blue ,sum(*nice_blue)*j] for j in range(0,xmax)] for i in range(0,ymax)]
        print(nice_blue)#(0.21568627450980393, 0.49411764705882355, 0.7215686274509804)
        print(nice_orange)#(1.0, 0.4980392156862745, 0.0)
        #Xblue = [[[0, 0, 5 * j, 5 * j] for j in range(0, xmax)] for i in range(0, ymax)]
        #Xblue = [[[0.21*50,.49*50,.72*50, 5*j] for j in range(0, xmax)] for i in range(0, ymax)]
        Xblue = [[[int(0.21*255.0),int(0.49*255.0),int(0.75*255.0), int(255*j/100)]
                  for j in range(0, xmax)] for i in range(0, ymax)]

        Xred = [[[int(1.0*255.0),int(0.49*255.0),int(0.0*255.0), int(255*i/100)]
                 for j in range(0,xmax)] for i in range(0,ymax)]

        Xblack = [[[int(.0*255.0),int(.0*255.0),int(.0*255.0), int(1.5*(i*i+j*j)/100)]
                 for j in range(0,xmax)] for i in range(0,ymax)]

        fig, ax = plt.subplots()
        ax.imshow(Xred)
        ax.imshow(Xblue)
        ax.imshow(Xblack)

        #ax.imshow(Xblue)
        ax.set_axis_off()
        plt.tight_layout()
        out_file_name=os.path.join(out_folder, "continuum_plot3.png")
        plt.savefig(out_file_name)
        plt.close()

    def test_print_nice_color_rgb(self):
        nice_blue = matplotlib.colors.to_rgb(two_d_colors.two_d_colors.high_Ne)
        nice_orange = matplotlib.colors.to_rgb(two_d_colors.two_d_colors.low_dT)
        print(nice_blue)#(0.21568627450980393, 0.49411764705882355, 0.7215686274509804)
        print(nice_orange)#(1.0, 0.4980392156862745, 0.0)

        light_blue= two_d_colors.lighten_color(nice_blue , 0.5)
        light_orange = two_d_colors.lighten_color(nice_orange,0.5)
        print(light_blue)#(0.5909098367380425, 0.7487652801706457, 0.8777176142423497)
        print(light_orange )#(1.0, 0.7490196078431373, 0.5)

        dark_blue= two_d_colors.lighten_color(nice_blue ,-0.5)
        dark_orange = two_d_colors.lighten_color(nice_orange, -0.5)
        print(dark_blue)#(1.4090901632619575, 1.2512347198293543, 1.1222823857576503)
        print(dark_orange)#(1.0, 1.2509803921568627, 1.5)
    def test_continuum_plot_2(self):

        out_folder= "/home/tamsen/Data/Specks_outout_from_mesx"

        xmin, xmax, ymin, ymax = (0, 100, 0, 100)
        #RGBA
        Xred = [[[2*i-2*j,0,0,3*i] for j in range(0,xmax)] for i in range(0,ymax)]
        Xred = [[[(2*i-2*j)*1,0,0,3*i] for j in range(0,xmax)] for i in range(0,ymax)]
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
