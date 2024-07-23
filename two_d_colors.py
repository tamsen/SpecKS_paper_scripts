# eventually this should be read in from an xml config or similar
import math
import xml.etree.ElementTree as ET


class two_d_colors:

    # A color blind/friendly color cycle for Matplotlib line plots
    # from https://gist.github.com/thriveth/8560036
    color_blind_friendly_color_cycle_analogs = {'blue': '#377eb8', 'orange': '#ff7f00', 'green': '#4daf4a',
                                                'pink': '#f781bf', 'brown': '#a65628', 'purple': '#984ea3',
                                                'gray': '#999999', 'red': '#e41a1c', 'yellow': '#dede00'}


    high_Ne=color_blind_friendly_color_cycle_analogs['blue'] #alpha
    low_dT = color_blind_friendly_color_cycle_analogs['orange'] #light vs dark

    high_Ne_high_dT=color_blind_friendly_color_cycle_analogs['blue']
    low_Ne_low_dT = color_blind_friendly_color_cycle_analogs['orange']
    high_Ne_low_dT = color_blind_friendly_color_cycle_analogs['brown']
    low_Ne_high_dT = color_blind_friendly_color_cycle_analogs['gray']


    nice_blue=(0.21568627450980393, 0.49411764705882355, 0.7215686274509804)
    nice_orange= (1.0, 0.4980392156862745, 0.0)

    light_blue = (0.5909098367380425, 0.7487652801706457, 0.8777176142423497)
    light_orange = (1.0, 0.7490196078431373, 0.5)

    dark_blue = (0.105, 0.25, 0.36)
    dark_orange = (1.0, 0.25, 0)

    high_dT_colors_by_Ne={
        "low":"white",
        "high":color_blind_friendly_color_cycle_analogs['blue']}

    medium_dT_colors_by_Ne={
        "low":light_orange,
        "high":dark_blue}

    low_dT_colors_by_Ne={
        "low":color_blind_friendly_color_cycle_analogs['orange'],
        "high":'darkgray'}

    colors_by_dT={
        "low": low_dT_colors_by_Ne,
        "high": high_dT_colors_by_Ne}
    #example sim names
    #['Allo6_S030W010.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo4_S050W040.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo2_S070W020.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo2_S070W060.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo1_S080W060.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo3_S060W040.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo7_S020W015.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo1_S080W075.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo1_S080W030.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo2_S070W050.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo4_S050W045.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo7_S020W010.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo6_S030W025.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo5_S040W030.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo4_S050W030.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo3_S060W050.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo5_S040W020.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo1_S080W070.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo3_S060W010.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo5_S040W035.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo3_S060W055.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo8_S010W005.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo6_S030W020.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo2_S070W065.hist_Ks_hist_fit2.0_sim37_N0p1', 'Auto1_S080W080.hist_Ks_hist_fit2.0_sim37_N0p1', 'Auto2_S070W070.hist_Ks_hist_fit2.0_sim37_N0p1', 'Auto3_S060W060.hist_Ks_hist_fit2.0_sim37_N0p1', 'Auto8_S010W010.hist_Ks_hist_fit2.0_sim37_N0p1', 'Auto6_S030W030.hist_Ks_hist_fit2.0_sim37_N0p1', 'Auto4_S050W050.hist_Ks_hist_fit2.0_sim37_N0p1', 'Auto5_S040W040.hist_Ks_hist_fit2.0_sim37_N0p1', 'Auto7_S020W020.hist_Ks_hist_fit2.0_sim37_N0p1', 'Allo6_S030W010.hist_Ks_hist_fit2.0_sim37_N1', 'Allo4_S050W040.hist_Ks_hist_fit2.0_sim37_N1', 'Allo2_S070W020.hist_Ks_hist_fit2.0_sim37_N1', 'Allo2_S070W060.hist_Ks_hist_fit2.0_sim37_N1', 'Allo1_S080W060.hist_Ks_hist_fit2.0_sim37_N1', 'Allo3_S060W040.hist_Ks_hist_fit2.0_sim37_N1', 'Allo7_S020W015.hist_Ks_hist_fit2.0_sim37_N1', 'Allo1_S080W075.hist_Ks_hist_fit2.0_sim37_N1', 'Allo1_S080W030.hist_Ks_hist_fit2.0_sim37_N1', 'Allo2_S070W050.hist_Ks_hist_fit2.0_sim37_N1', 'Allo4_S050W045.hist_Ks_hist_fit2.0_sim37_N1', 'Allo7_S020W010.hist_Ks_hist_fit2.0_sim37_N1', 'Allo6_S030W025.hist_Ks_hist_fit2.0_sim37_N1', 'Allo5_S040W030.hist_Ks_hist_fit2.0_sim37_N1', 'Allo4_S050W030.hist_Ks_hist_fit2.0_sim37_N1', 'Allo3_S060W050.hist_Ks_hist_fit2.0_sim37_N1', 'Allo5_S040W020.hist_Ks_hist_fit2.0_sim37_N1', 'Allo1_S080W070.hist_Ks_hist_fit2.0_sim37_N1', 'Allo3_S060W010.hist_Ks_hist_fit2.0_sim37_N1', 'Allo5_S040W035.hist_Ks_hist_fit2.0_sim37_N1', 'Allo3_S060W055.hist_Ks_hist_fit2.0_sim37_N1', 'Allo8_S010W005.hist_Ks_hist_fit2.0_sim37_N1', 'Allo6_S030W020.hist_Ks_hist_fit2.0_sim37_N1', 'Allo2_S070W065.hist_Ks_hist_fit2.0_sim37_N1', 'Auto1_S080W080.hist_Ks_hist_fit2.0_sim37_N1', 'Auto2_S070W070.hist_Ks_hist_fit2.0_sim37_N1', 'Auto3_S060W060.hist_Ks_hist_fit2.0_sim37_N1', 'Auto8_S010W010.hist_Ks_hist_fit2.0_sim37_N1', 'Auto6_S030W030.hist_Ks_hist_fit2.0_sim37_N1', 'Auto4_S050W050.hist_Ks_hist_fit2.0_sim37_N1', 'Auto5_S040W040.hist_Ks_hist_fit2.0_sim37_N1', 'Auto7_S020W020.hist_Ks_hist_fit2.0_sim37_N1', 'Allo6_S030W010.hist_Ks_hist_fit2.0_sim37_N5', 'Allo4_S050W040.hist_Ks_hist_fit2.0_sim37_N5', 'Allo2_S070W020.hist_Ks_hist_fit2.0_sim37_N5', 'Allo2_S070W060.hist_Ks_hist_fit2.0_sim37_N5', 'Allo1_S080W060.hist_Ks_hist_fit2.0_sim37_N5', 'Allo3_S060W040.hist_Ks_hist_fit2.0_sim37_N5', 'Allo7_S020W015.hist_Ks_hist_fit2.0_sim37_N5', 'Allo1_S080W075.hist_Ks_hist_fit2.0_sim37_N5', 'Allo1_S080W030.hist_Ks_hist_fit2.0_sim37_N5', 'Allo2_S070W050.hist_Ks_hist_fit2.0_sim37_N5', 'Allo4_S050W045.hist_Ks_hist_fit2.0_sim37_N5', 'Allo7_S020W010.hist_Ks_hist_fit2.0_sim37_N5', 'Allo6_S030W025.hist_Ks_hist_fit2.0_sim37_N5', 'Allo5_S040W030.hist_Ks_hist_fit2.0_sim37_N5', 'Allo4_S050W030.hist_Ks_hist_fit2.0_sim37_N5', 'Allo3_S060W050.hist_Ks_hist_fit2.0_sim37_N5', 'Allo5_S040W020.hist_Ks_hist_fit2.0_sim37_N5', 'Allo1_S080W070.hist_Ks_hist_fit2.0_sim37_N5', 'Allo3_S060W010.hist_Ks_hist_fit2.0_sim37_N5', 'Allo5_S040W035.hist_Ks_hist_fit2.0_sim37_N5', 'Allo3_S060W055.hist_Ks_hist_fit2.0_sim37_N5', 'Allo8_S010W005.hist_Ks_hist_fit2.0_sim37_N5', 'Allo6_S030W020.hist_Ks_hist_fit2.0_sim37_N5', 'Allo2_S070W065.hist_Ks_hist_fit2.0_sim37_N5', 'Auto1_S080W080.hist_Ks_hist_fit2.0_sim37_N5', 'Auto2_S070W070.hist_Ks_hist_fit2.0_sim37_N5', 'Auto3_S060W060.hist_Ks_hist_fit2.0_sim37_N5', 'Auto8_S010W010.hist_Ks_hist_fit2.0_sim37_N5', 'Auto6_S030W030.hist_Ks_hist_fit2.0_sim37_N5', 'Auto4_S050W050.hist_Ks_hist_fit2.0_sim37_N5', 'Auto5_S040W040.hist_Ks_hist_fit2.0_sim37_N5', 'Auto7_S020W020.hist_Ks_hist_fit2.0_sim37_N5', 'Allo6_S030W010.hist_Ks_hist_fit2.0_sim36_N10', 'Allo4_S050W040.hist_Ks_hist_fit2.0_sim36_N10', 'Allo2_S070W020.hist_Ks_hist_fit2.0_sim36_N10', 'Allo2_S070W060.hist_Ks_hist_fit2.0_sim36_N10', 'Allo1_S080W060.hist_Ks_hist_fit2.0_sim36_N10', 'Allo3_S060W040.hist_Ks_hist_fit2.0_sim36_N10', 'Allo7_S020W015.hist_Ks_hist_fit2.0_sim36_N10', 'Allo1_S080W075.hist_Ks_hist_fit2.0_sim36_N10', 'Allo1_S080W030.hist_Ks_hist_fit2.0_sim36_N10', 'Allo2_S070W050.hist_Ks_hist_fit2.0_sim36_N10', 'Allo4_S050W045.hist_Ks_hist_fit2.0_sim36_N10', 'Allo7_S020W010.hist_Ks_hist_fit2.0_sim36_N10', 'Allo6_S030W025.hist_Ks_hist_fit2.0_sim36_N10', 'Allo5_S040W030.hist_Ks_hist_fit2.0_sim36_N10', 'Allo4_S050W030.hist_Ks_hist_fit2.0_sim36_N10', 'Allo3_S060W050.hist_Ks_hist_fit2.0_sim36_N10', 'Allo5_S040W020.hist_Ks_hist_fit2.0_sim36_N10', 'Allo1_S080W070.hist_Ks_hist_fit2.0_sim36_N10', 'Allo3_S060W010.hist_Ks_hist_fit2.0_sim36_N10', 'Allo5_S040W035.hist_Ks_hist_fit2.0_sim36_N10', 'Allo3_S060W055.hist_Ks_hist_fit2.0_sim36_N10', 'Allo8_S010W005.hist_Ks_hist_fit2.0_sim36_N10', 'Allo6_S030W020.hist_Ks_hist_fit2.0_sim36_N10', 'Allo2_S070W065.hist_Ks_hist_fit2.0_sim36_N10', 'Auto1_S080W080.hist_Ks_hist_fit2.0_sim36_N10', 'Auto2_S070W070.hist_Ks_hist_fit2.0_sim36_N10', 'Auto3_S060W060.hist_Ks_hist_fit2.0_sim36_N10', 'Auto8_S010W010.hist_Ks_hist_fit2.0_sim36_N10', 'Auto6_S030W030.hist_Ks_hist_fit2.0_sim36_N10', 'Auto4_S050W050.hist_Ks_hist_fit2.0_sim36_N10', 'Auto5_S040W040.hist_Ks_hist_fit2.0_sim36_N10', 'Auto7_S020W020.hist_Ks_hist_fit2.0_sim36_N10', 'Allo6_S030W010.hist_Ks_hist_fit2.0_sim37_N20', 'Allo4_S050W040.hist_Ks_hist_fit2.0_sim37_N20', 'Allo2_S070W020.hist_Ks_hist_fit2.0_sim37_N20', 'Allo2_S070W060.hist_Ks_hist_fit2.0_sim37_N20', 'Allo1_S080W060.hist_Ks_hist_fit2.0_sim37_N20', 'Allo3_S060W040.hist_Ks_hist_fit2.0_sim37_N20', 'Allo7_S020W015.hist_Ks_hist_fit2.0_sim37_N20', 'Allo1_S080W075.hist_Ks_hist_fit2.0_sim37_N20', 'Allo1_S080W030.hist_Ks_hist_fit2.0_sim37_N20', 'Allo2_S070W050.hist_Ks_hist_fit2.0_sim37_N20', 'Allo4_S050W045.hist_Ks_hist_fit2.0_sim37_N20', 'Allo7_S020W010.hist_Ks_hist_fit2.0_sim37_N20', 'Allo6_S030W025.hist_Ks_hist_fit2.0_sim37_N20', 'Allo5_S040W030.hist_Ks_hist_fit2.0_sim37_N20', 'Allo4_S050W030.hist_Ks_hist_fit2.0_sim37_N20', 'Allo3_S060W050.hist_Ks_hist_fit2.0_sim37_N20', 'Allo5_S040W020.hist_Ks_hist_fit2.0_sim37_N20', 'Allo1_S080W070.hist_Ks_hist_fit2.0_sim37_N20', 'Allo3_S060W010.hist_Ks_hist_fit2.0_sim37_N20', 'Allo5_S040W035.hist_Ks_hist_fit2.0_sim37_N20', 'Allo3_S060W055.hist_Ks_hist_fit2.0_sim37_N20', 'Allo8_S010W005.hist_Ks_hist_fit2.0_sim37_N20', 'Allo6_S030W020.hist_Ks_hist_fit2.0_sim37_N20', 'Allo2_S070W065.hist_Ks_hist_fit2.0_sim37_N20', 'Auto1_S080W080.hist_Ks_hist_fit2.0_sim37_N20', 'Auto2_S070W070.hist_Ks_hist_fit2.0_sim37_N20', 'Auto3_S060W060.hist_Ks_hist_fit2.0_sim37_N20', 'Auto8_S010W010.hist_Ks_hist_fit2.0_sim37_N20', 'Auto6_S030W030.hist_Ks_hist_fit2.0_sim37_N20', 'Auto4_S050W050.hist_Ks_hist_fit2.0_sim37_N20', 'Auto5_S040W040.hist_Ks_hist_fit2.0_sim37_N20', 'Auto7_S020W020.hist_Ks_hist_fit2.0_sim37_N20']"


def get_color_from_name(name):

    dT=get_dT_from_name(name)
    Ne=get_Ne_from_name(name)

    dT_amount="high"
    if dT == 0:
        dT_amount = "low"

    Ne_amount="high"
    if Ne < 1:
        Ne_amount = "low"

    color= two_d_colors.colors_by_dT[dT_amount][Ne_amount]
    return color
def get_Ne_from_name(name):
    splat= name.split("N")
    important_part=splat[-1]
    print(important_part)
    as_decimal=float(important_part.replace("p","."))
    print(as_decimal)
    return as_decimal
def get_dT_from_name(name):
    splat= name.split("_")
    important_parts=splat[1].replace(".hist","").split("W")
    print(important_parts)
    SPC_part=important_parts[0].replace("S","")
    WGD_part=important_parts[1]
    SPC_decimal=float(SPC_part)
    WGD_decimal=float(WGD_part)
    print([SPC_decimal,WGD_decimal])
    return SPC_decimal - WGD_decimal

def get_color_from_dT(name):

    dT=get_dT_from_name(name)
    base_color=two_d_colors.low_dT

    amount=1
    if dT == 0:
        return base_color

    if dT > 30:
        return 'tan'

    return 'orange'

def get_color_from_Ne(name):

    Ne=get_Ne_from_name(name)
    final_color='darkblue'
    if Ne < 20:
        final_color=two_d_colors.high_Ne
    if Ne < 10:
        final_color='lightblue'

    return final_color
def get_alpha_from_Ne(name):

    #Ne goes from 0.1 to 20 in out sims
    Ne=get_Ne_from_name(name)
    alpha=0.5
    if Ne < 20:
        alpha = 0.4
    if Ne < 10:
        alpha = 0.3
    if Ne < 5:
        alpha = 0.2

    return alpha

def get_alpha_from_dT(name):

    #dT goes from 0 to 50
    dT=get_dT_from_name(name)
    alpha=0.5
    if dT < 40:
        alpha = 0.4
    if dT < 20:
        alpha = 0.3
    if dT < 10:
        alpha = 0.2

    return alpha



#from  https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
#from, https://gist.github.com/ihincks/6a420b599f43fcd7dbd79d56798c4e5a
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    result= colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
    positive_r=[max([0,r])for r in result]
    print(result)
    return positive_r


