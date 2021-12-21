########################################## ॐ  ##########################################

## Import required packages ##

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.coordinates import Angle, get_moon, get_sun
from datetime import datetime
from geopy.geocoders import Nominatim
import pytz
from astropy.time import Time, TimeDelta
from astropy.table import Table
import matplotlib.patches as mpatches

from plot_graha import *
from calc_rahu_ketu_pos import *

##########################################

## Sanskrit typesetting using XeLaTeX ##
mpl.use("pgf")

## TeX preamble ##
preamble = r"""\usepackage{fontspec}
               \usepackage{polyglossia}
               \usepackage[T1]{fontenc}
               \usepackage{tikz}
               \setmainlanguage{english}
               \setotherlanguages{sanskrit}
               \newfontfamily\devanagarifont[Script=Devanagari]{Shobhika-Regular}
               \newcommand{\sam}[1]{\begin{sanskrit}#1\end{sanskrit}}"""

params = {
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'pgf.texsystem': 'xelatex',
    'pgf.preamble': preamble,
}

mpl.rcParams.update(params)

##########################################

def get_lat_lon(ra,dec,inclination):
    """
    Convert (ra,dec) equatorial coordinates of an object to (lat,lon) ecliptic coordinates
    
    Inputs :
    ra = right ascension of object in deg
    dec = declination of object in deg
    inclination = inclination of ecliptic to the celectial equator in deg

    Outputs :
    lat = ecliptic latitude in deg
    lon = ecliptic longtude in deg
    """
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)
    inclination = np.deg2rad(inclination)

    beta_coord = np.arcsin(np.sin(dec)*np.cos(inclination) - np.cos(dec)*np.sin(inclination)*np.sin(ra))
    lambda_coord = np.arccos(np.cos(ra)*np.cos(dec)/np.cos(beta_coord))
    if (np.pi <= ra < 2*np.pi):
        lambda_coord = np.pi + np.arccos(-np.cos(ra)*np.cos(dec)/np.cos(beta_coord))
    if (np.pi/2 <= dec < 3*np.pi/2):
        beta_coord = (np.pi/2) + np.arccos(np.sin(dec)*np.cos(inclination) - np.cos(dec)*np.sin(inclination)*np.sin(ra))
    
    lat = np.rad2deg(beta_coord)
    lon = np.rad2deg(lambda_coord)
    
    return lat, lon

def calc_nakshatra_tithi(time,filename="nakshatra_at_test_time.pdf",tz_str="Asia/Calcutta",time_format="%Y-%m-%dT%H:%M:%S",inclination=23.4,plot_other_grahas=False):
    """
    Calculates nakshatra and tithi at input time and makes plot of grahas

    Inputs :
    time = time at which to calculate panchanga
    filename = name of file to write the plot of position of grahas in nakshatra and rashi
    tz_str = time zone of the location at which to calculate panchanga
    time_format = format of input time
    inclination = inclination of ecliptic to the celectial equator in deg

    Output :
    Plot saved at filename
    """
    naksh_names_tab = Table.read('nakshatram_names.tex',format="latex")
    nakshatram_names = naksh_names_tab['names'].data

    nakshatram_extent = 360/27 ## deg
    pAdam_extent = nakshatram_extent/4 ## deg

    start_coord = 23 + 46/60
    nakshatram_coords = np.linspace(start_coord,start_coord+nakshatram_extent*26,27)
    nakshatram_coord_labels = []

    for i in range(len(nakshatram_coords)):
        coord = nakshatram_coords[i]
        if (coord > 360):
            coord = nakshatram_coords[i] - 360
        deg_coord = coord - coord%1
        min_coord = 60*(coord%1) 
        nakshatram_coord_labels.append(r"%d$^{{\circ}}$\,%d$^{{\prime}}$"%(deg_coord,min_coord))

    nakshatram_coord_labels = np.array(nakshatram_coord_labels)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(7,5))
    ax.set_rmax(2)
    ax.set_yticks([1])
    ax.set_yticklabels([""])
    ax.grid(axis="x")

    fmt = time_format
    date_str = time

    tz = pytz.timezone(tz_str)
    test_date = datetime.strptime(date_str, fmt)
    day_of_the_week_at_t = test_date.strftime('%A')
    local_time = tz.localize(test_date,is_dst=None)
    test_date_utc = local_time.astimezone(pytz.utc)

    test_date_utc_time = Time(test_date_utc.strftime(fmt),format="isot",scale="utc")

    moon_ra = (get_moon(test_date_utc_time)).ra.deg
    sun_ra = (get_sun(test_date_utc_time)).ra.deg
    moon_dec =  (get_moon(test_date_utc_time)).dec.deg
    sun_dec =  (get_sun(test_date_utc_time)).dec.deg

    moon_beta, moon_lambda = get_lat_lon(moon_ra,moon_dec,inclination)
    sun_beta, sun_lambda = get_lat_lon(sun_ra,sun_dec,inclination)

    each_tithi = 360/30

    if (sun_lambda > moon_lambda):
        d_theta = 360 - (sun_lambda - moon_lambda)
    else:
        d_theta = (moon_lambda - sun_lambda)

    final_tithi = d_theta/each_tithi
    final_nakshatram = ""

    for coord_id, coord in enumerate(nakshatram_coords):
        ax.plot(np.deg2rad([coord,coord]),[1,2.2],"k",alpha=0.3)
        if (coord > 360) :
            nakshatram_coords[coord_id] = coord - 360
        text_theta = (2*nakshatram_coords[coord_id] + nakshatram_extent) * (np.pi/360)
        if (90 < coord < 270):
            if ( coord < moon_lambda < coord + nakshatram_extent):
                ax.text(text_theta,1.5,nakshatram_names[coord_id],rotation= text_theta*180/np.pi + 180, ha="center", va="center")
                final_nakshatram = nakshatram_names[coord_id]
            else:
                ax.text(text_theta,1.5,nakshatram_names[coord_id],rotation= text_theta*180/np.pi + 180, ha="center", va="center",alpha=0.3)
        else:
            if ( coord < moon_lambda < coord + nakshatram_extent):
                ax.text(text_theta,1.5,nakshatram_names[coord_id],rotation=text_theta*180/np.pi, ha="center", va="center")
                final_nakshatram = nakshatram_names[coord_id]
            else:
                ax.text(text_theta,1.5,nakshatram_names[coord_id],rotation=text_theta*180/np.pi, ha="center", va="center",alpha=0.3)
            
    ax.set_xticks(np.pi/180 * nakshatram_coords)
    ax.set_xticklabels(nakshatram_coord_labels,fontsize=5)

    plot_moon_phase(final_tithi,np.array([np.pi/180 * moon_lambda,2.0]),0.1,fig,ax)
    plot_sun(np.array([np.pi/180 * sun_lambda,2.0]),0.1,fig,ax)

    #####################################
    
    rashi_names_tab = Table.read('rashi_names.tex',format="latex")
    rashi_names = rashi_names_tab['names'].data

    rashi_start = 23 + 46/60
    rashi_extent = 360/12
    rashi_coords = np.linspace(rashi_start,rashi_start+rashi_extent*11,12)

    for coord_id,coord in enumerate(rashi_coords):
        ax.plot(np.deg2rad([coord,coord]), [0,1],"k",alpha=0.3)
        rashi_text_theta = (2*rashi_coords[coord_id] + rashi_extent) * (np.pi/360)
        if (77 < coord < 255):
            if ( coord < sun_lambda < coord + rashi_extent):
                ax.text(rashi_text_theta,0.625,rashi_names[coord_id],rotation=rashi_text_theta*180/np.pi + 180, ha="center", va="center")
            else:
                ax.text(rashi_text_theta,0.625,rashi_names[coord_id],rotation=rashi_text_theta*180/np.pi + 180, ha="center", va="center",alpha=0.3)
        else:
            if ( coord < sun_lambda < coord + nakshatram_extent):
                ax.text(rashi_text_theta,0.625,rashi_names[coord_id],rotation=rashi_text_theta*180/np.pi, ha="center", va="center")
            else:
                ax.text(rashi_text_theta,0.625,rashi_names[coord_id],rotation=rashi_text_theta*180/np.pi, ha="center", va="center",alpha=0.3)
    
    ################  rahu ketu pos  #####################       
        
    known_eclipse_time = Time("2009-07-22T02:36:25",format="isot",scale="utc") ### This eclipse happend when moon was in Ketu's postion
    moon_ra_at_known_eclipse = (get_moon(known_eclipse_time)).ra.deg
    moon_dec_at_known_eclipse = (get_moon(known_eclipse_time)).dec.deg

    moon_beta_known_eclipse, moon_lambda_known_eclipse = get_lat_lon(moon_ra_at_known_eclipse,moon_dec_at_known_eclipse,inclination)


    rahu_lambda, ketu_lambda = calc_rahu_ketu_pos(known_eclipse_time,test_date_utc_time,moon_lambda_known_eclipse)

    plot_rahu([np.deg2rad(rahu_lambda),1.8],5,fig,ax)
    plot_ketu([np.deg2rad(ketu_lambda),1.8],5,fig,ax)

    ################legend for grahas#####################
  
    inv = ax.transData.inverted()
    
    sun_legend_pos = inv.transform((750,400))
    sun_label_pos = inv.transform((790,395))
    moon_legend_pos = inv.transform((750,350))
    moon_label_pos = inv.transform((790,345))
    rahu_legend_pos = inv.transform((750,300))
    rahu_label_pos = inv.transform((790,295))
    ketu_legend_pos = inv.transform((750,250))
    ketu_label_pos = inv.transform((790,245))

    plot_sun(sun_legend_pos,0.1,fig,ax)
    plt.text(sun_label_pos[0],sun_label_pos[1],"\sam{सूर्यः } ")
    plot_moon_phase(final_tithi,moon_legend_pos,0.1,fig,ax)
    plt.text(moon_label_pos[0],moon_label_pos[1],"\sam{चन्द्रः } ")

    plot_rahu(rahu_legend_pos,10,fig,ax)
    plt.text(rahu_label_pos[0],rahu_label_pos[1],"\sam{राहुः} ")
    plot_ketu(ketu_legend_pos,10,fig,ax)
    plt.text(ketu_label_pos[0],ketu_label_pos[1], "\sam{केतुः} ")

    ################legend for panchanga#####################

    plt.text(inv.transform((-200,400))[0],inv.transform((-200,400))[1],tz_str)
    plt.text(inv.transform((-200,370))[0],inv.transform((-200,370))[1],date_str.split("T")[0])
    plt.text(inv.transform((-200,340))[0],inv.transform((-200,340))[1],date_str.split("T")[1])

    tithi_names_tab = Table.read('tithi_names.tex',format="latex")
    tithi_number = tithi_names_tab['tithi'].data
    tithi_names = tithi_names_tab['tithi_names'].data

    plt.text(inv.transform((-200,310))[0],inv.transform((-200,310))[1],tithi_names[int(final_tithi)])

    vaara_names_tab = Table.read('vaara_names.tex',format="latex")
    weekday = vaara_names_tab['weekday'].data
    vaara_names = vaara_names_tab['vaara'].data
    
    plt.text(inv.transform((-200,280))[0],inv.transform((-200,280))[1],vaara_names[np.where(weekday==day_of_the_week_at_t)][0])
    
    plt.text(inv.transform((-200,250))[0],inv.transform((-200,250))[1],final_nakshatram)

    plt.title(r'\sam{ॐ }', fontsize=20)
    plt.ylim(0,2.2)

    plt.savefig(filename)
    
    return 0

tz_str = "Asia/Calcutta"
date_str = "1995-08-06T10:05:00"

calc_nakshatra_tithi(date_str,tz_str=tz_str)
