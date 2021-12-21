from astropy.coordinates import Angle, get_moon, get_sun
from datetime import datetime
import pytz
from astropy.time import Time
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import numpy as np
import matplotlib.cbook as cbook
import matplotlib.transforms as transforms


#### Colors for grahas ####
#Sun - Surya - Red - circle 
#Moon - Chandra - White - semicircle
#Mars - Mangala - Red - triangle
#Mercury - Budha - Green - bow (dhanusha)
#Jupiter - Guru - Yellow - lotus
#Venus - Shukra - White - square
#Saturn - Shani - Black - snake
#Rahu - Blue - crocodile
#Ketu - Grey - sword



def plot_moon_phase(day,drawing_origin,radius,fig,ax):
    center = np.array([0,0]) 
    trans = (fig.dpi_scale_trans + transforms.ScaledTranslation(drawing_origin[0],drawing_origin[1], ax.transData))    
    if (0 < day <= 7.5):
        circle = mpatches.Circle(center,radius,ec="k",fc="w",transform=trans,clip_on=False)
        ellipse = mpatches.Ellipse(center,2*radius*np.sin(np.pi/2 - day*12*np.pi/180),2*radius,fc="k",ec=None,transform=trans,clip_on=False)
        wedge = mpatches.Wedge(center,radius,90,270,fc="k",ec=None,transform=trans,clip_on=False)
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if (7.5 < day <= 15):
        circle = mpatches.Circle(center,radius,ec="k",fc="w",transform=trans,clip_on=False)
        ellipse = mpatches.Ellipse(center,2*radius*np.sin(np.pi/2 - day*12*np.pi/180),2*radius,fc="w",ec=None,transform=trans,clip_on=False)
        wedge = mpatches.Wedge(center,radius,90,270,fc="k",ec=None,transform=trans,clip_on=False)
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if (15 < day <= 22.5):
        circle = mpatches.Circle(center,radius,ec="k",fc="w",transform=trans,clip_on=False)
        ellipse = mpatches.Ellipse(center,2*radius*np.sin(np.pi/2 - (day-15)*12*np.pi/180),2*radius,fc="w",ec=None,transform=trans,clip_on=False)
        wedge = mpatches.Wedge(center,radius,270,90,fc="k",ec=None,transform=trans,clip_on=False)
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if (22.5 < day <= 30):
        circle = mpatches.Circle(center,radius,ec="k",fc="w",transform=trans,clip_on=False)
        ellipse = mpatches.Ellipse(center,2*radius*np.sin(np.pi/2 - (day-15)*12*np.pi/180),2*radius,fc="k",ec=None,transform=trans,clip_on=False)
        wedge = mpatches.Wedge(center,radius,270,90,fc="k",ec=None,transform=trans,clip_on=False)
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)
    
    return 0


def plot_sun(drawing_origin,radius,fig,ax):
    center = np.array([0,0])
    trans = (fig.dpi_scale_trans + transforms.ScaledTranslation(drawing_origin[0],drawing_origin[1], ax.transData))     

    circle = mpatches.Circle(center,radius,ec="orange",fc="yellow",transform=trans,clip_on=False)
    ax.add_patch(circle)

    for theta in np.arange(0,2*np.pi,np.pi/6):
        triangle = mpatches.RegularPolygon([1.2*radius*np.cos(theta),1.2*radius*np.sin(theta)],3,0.2*radius,ec="yellow",fc="orange",orientation=theta+np.pi/6,transform=trans,clip_on=False)
        ax.add_patch(triangle)
    
    return 0

def plot_inner_graha_phase(graha,angle_to_sun,drawing_origin,radius,fig,ax):
    center = np.array([0,0])
    trans = (fig.dpi_scale_trans + transforms.ScaledTranslation(drawing_origin[0],drawing_origin[1], ax.transData))
    
    if (graha == "Budha"):
        bright_side_color = "green"
        dark_side_color = "darkgreen"
    if (graha == "Shukra"):
        bright_side_color = "white"
        dark_side_color = "grey"
    else:
        print("That is not an inner graha! Use the outer graha function to plot ",graha)
        return   

    if (0 < (angle_to_sun/12) <= 7.5):
        circle = mpatches.Circle(center,radius,ec=dark_side_color,fc=bright_side_color,transform=trans,clip_on=False)
        ellipse = mpatches.Ellipse(center,2*radius*np.sin(np.pi/2 - (angle_to_sun/12)*12*np.pi/180),2*radius,fc=dark_side_color,ec=None,transform=trans,clip_on=False)
        wedge = mpatches.Wedge(center,radius,90,270,fc=dark_side_color,ec=None,transform=trans,clip_on=False)
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if (7.5 < (angle_to_sun/12) <= 15):
        circle = mpatches.Circle(center,radius,ec=dark_side_color,fc=bright_side_color,transform=trans,clip_on=False)
        ellipse = mpatches.Ellipse(center,2*radius*np.sin(np.pi/2 - (angle_to_sun/12)*12*np.pi/180),2*radius,fc=bright_side_color,ec=None,transform=trans,clip_on=False)
        wedge = mpatches.Wedge(center,radius,90,270,fc=dark_side_color,ec=None,transform=trans,clip_on=False)
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if (15 < (angle_to_sun/12) <= 22.5):
        circle = mpatches.Circle(center,radius,ec=dark_side_color,fc=bright_side_color,transform=trans,clip_on=False)
        ellipse = mpatches.Ellipse(center,2*radius*np.sin(np.pi/2 - ((angle_to_sun/12)-15)*12*np.pi/180),2*radius,fc=bright_side_color,ec=None,transform=trans,clip_on=False)
        wedge = mpatches.Wedge(center,radius,270,90,fc=dark_side_color,ec=None,transform=trans,clip_on=False)
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if (22.5 < (angle_to_sun/12) <= 30):
        circle = mpatches.Circle(center,radius,ec=dark_side_color,fc=bright_side_color,transform=trans,clip_on=False)
        ellipse = mpatches.Ellipse(center,2*radius*np.sin(np.pi/2 - ((angle_to_sun/12)-15)*12*np.pi/180),2*radius,fc=dark_side_color,ec=None,transform=trans,clip_on=False)
        wedge = mpatches.Wedge(center,radius,270,90,fc=dark_side_color,ec=None,transform=trans,clip_on=False)
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    return 0

def plot_rahu(drawing_origin,scale,fig,ax):
    
    center = np.array([0,0])
    trans = (transforms.Affine2D().scale(scale) + transforms.ScaledTranslation(drawing_origin[0],drawing_origin[1], ax.transData))
    
    from matplotlib.path import Path

    verts = [
    (-1, 1),    # Vert1 
    (1., 1.),   # Vert2
    (1., -1.),  # Pivot1
    (0, -1.),   # Vert3
    (-1., -1.), # Pivot2 
    (-1., 1.),
    (-1., 1.),  
    ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    path = Path(verts, codes)

    patch = mpatches.PathPatch(path, facecolor='blue',edgecolor='lightblue',lw=scale/5,transform=trans,clip_on=False)
    ax.add_patch(patch)

    return 0

def plot_ketu(drawing_origin,scale,fig,ax):

    center = np.array([0,0])
    trans = (transforms.Affine2D().scale(scale) + transforms.ScaledTranslation(drawing_origin[0],drawing_origin[1], ax.transData))

    from matplotlib.path import Path

    verts = [
    (-1, 1),    # Vert1 
    (1., 1.),   # Vert2
    (1., -1.),  # Vert3
    (-1., -1.), # Vert4
    (0, 0),     # Vert5 
    (-1., 1.),  # Vert6
    (-1., 1.),
    ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    path = Path(verts, codes)

    patch = mpatches.PathPatch(path, facecolor='dimgrey',edgecolor='lightgrey',lw=scale/5,transform=trans,clip_on=False)
    ax.add_patch(patch)

    return 0

