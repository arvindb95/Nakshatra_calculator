########################################## ॐ  ##########################################

## Import required packages ##

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.coordinates import get_body, EarthLocation
from astropy.time import Time
from astropy.table import Table
from datetime import datetime
import pytz
from timezonefinder import TimezoneFinder

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
               \newfontfamily\devanagarifont[Script=Devanagari]{Sanskrit 2003}
               \newcommand{\sam}[1]{\begin{sanskrit}#1\end{sanskrit}}"""

params = {
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'pgf.texsystem': 'xelatex',
    'pgf.preamble': preamble,
}

mpl.rcParams.update(params)

plt.style.use('dark_background')

##########################################


def calc_nakshatra_tithi(location, time, time_format="%Y-%m-%d %H:%M:%S", filename="nakshatra_at_test_time.pdf"):
    """
    Calculates nakshatra and tithi at input time and makes plot of grahas

    Inputs :
    location = enter address of your location. Eg: "Mumbai, India"
    time = time at which to calculate panchanga
    time_format = format of input time
    filename = name of file to write the plot of position of grahas in nakshatra and rashi

    Output :
    Plot saved at filename
    """

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(7, 5))

    ax.set_rmax(2.5)
    ax.set_yticks([1])
    ax.set_yticklabels([""])
    ax.grid(axis="x")

    observing_location = EarthLocation.of_address(location)
    tz_obj = TimezoneFinder()
    timezone_location = tz_obj.timezone_at(
        lng=observing_location.lon.value, lat=observing_location.lat.value)

    fmt = time_format
    date_str = time

    tz = pytz.timezone(timezone_location)
    test_date = datetime.strptime(date_str, fmt)

    local_time = tz.localize(test_date, is_dst=None)
    test_date_utc = local_time.astimezone(pytz.utc)

    test_date_utc_time = Time(
        test_date_utc.strftime(fmt), format="iso", scale="utc")

    print("Time in UTC : ", test_date_utc_time)

    moon_coord = get_body("moon", test_date_utc_time)
    sun_coord = get_body("sun", test_date_utc_time)

    moon_lambda = moon_coord.geocentrictrueecliptic.lon.value
    sun_lambda = sun_coord.geocentrictrueecliptic.lon.value

    ### -------- Thithi -------- ###

    each_tithi = 360/30

    if (sun_lambda > moon_lambda):
        tithi_theta = 360 - (sun_lambda - moon_lambda)
    else:
        tithi_theta = (moon_lambda - sun_lambda)

    final_tithi = tithi_theta/each_tithi

    ### -------- Vaara -------- ###

    day_of_the_week_at_t = test_date.strftime('%A')

    ### -------- Nakshatram -------- ###

    naksh_names_tab = Table.read('nakshatram_names.tex', format="latex")
    nakshatram_names = naksh_names_tab['names'].data

    nakshatram_extent = 360/27  # deg
    pAdam_extent = nakshatram_extent/4  # deg

    # Ref Indian Astronomy : An Introduction - S Balachandra Rao - Pg 35
    start_coord = 23 + 47/60 + 14.1/3600
    nakshatram_coords = np.linspace(
        start_coord, start_coord+nakshatram_extent*26, 27)
    nakshatram_coord_labels = []

    for i in range(len(nakshatram_coords)):
        coord = nakshatram_coords[i]
        if (coord > 360):
            coord = nakshatram_coords[i] - 360
        deg_coord = coord - coord % 1
        min_coord = 60*(coord % 1)
        nakshatram_coord_labels.append(
            r"%d$^{{\circ}}$\,%d$^{{\prime}}$" % (deg_coord, min_coord))

    nakshatram_coord_labels = np.array(nakshatram_coord_labels)

    for i in range(len(nakshatram_coords)):
        if (nakshatram_coords[i] > 360):
            nakshatram_coords[i] = nakshatram_coords[i] - 360

    final_nakshatram = ""

    for coord_id, coord in enumerate(nakshatram_coords):
        ax.plot(np.deg2rad([coord, coord]), [1, 2.2], "grey", alpha=0.3)
        if (coord > 360):
            nakshatram_coords[coord_id] = coord - 360
        text_theta = (
            2*nakshatram_coords[coord_id] + nakshatram_extent) * (np.pi/360)
        if (90 < coord < 270):
            if (coord < moon_lambda < coord + nakshatram_extent):
                ax.text(text_theta, 1.5, nakshatram_names[coord_id],
                        rotation=text_theta*180/np.pi + 180, ha="center", va="center")
                final_nakshatram = nakshatram_names[coord_id]
            else:
                ax.text(text_theta, 1.5, nakshatram_names[coord_id], rotation=text_theta *
                        180/np.pi + 180, ha="center", va="center", alpha=0.3)
        else:
            if (coord < moon_lambda < coord + nakshatram_extent):
                ax.text(text_theta, 1.5, nakshatram_names[coord_id],
                        rotation=text_theta*180/np.pi, ha="center", va="center")
                final_nakshatram = nakshatram_names[coord_id]
            elif ((moon_lambda + 360 - coord) < nakshatram_extent):
                ax.text(text_theta, 1.5, nakshatram_names[coord_id],
                        rotation=text_theta*180/np.pi, ha="center", va="center")
                final_nakshatram = nakshatram_names[coord_id]
            else:
                ax.text(text_theta, 1.5, nakshatram_names[coord_id],
                        rotation=text_theta*180/np.pi, ha="center", va="center", alpha=0.3)
    ax.set_xticks(np.pi/180 * nakshatram_coords)
    ax.set_xticklabels(nakshatram_coord_labels, fontsize=5)

    ### -------- Yoga --------###
    yoga_names_tab = Table.read('yoga_names.tex', format="latex")
    yoga_names = yoga_names_tab['names'].data

    each_yoga = nakshatram_extent

    add_theta = sun_lambda + moon_lambda

    if (add_theta > 360):
        yoga_theta = add_theta - 360
    else:
        yoga_theta = add_theta

    final_yoga_id = ((yoga_theta - 2*start_coord)/each_yoga)

    final_yoga = yoga_names[int(final_yoga_id)]


# Plot sun and moon

    plot_moon_phase(final_tithi, np.array(
        [np.pi/180 * moon_lambda, 1]), 0.1, fig, ax)
    plot_sun(np.array([np.pi/180 * sun_lambda, 2.05]), 0.1, fig, ax)

    #####################################

    rashi_names_tab = Table.read('rashi_names.tex', format="latex")
    rashi_names = rashi_names_tab['names'].data

    rashi_start = 23 + 46/60
    rashi_extent = 360/12
    rashi_coords = np.linspace(rashi_start, rashi_start+rashi_extent*11, 12)

    for coord_id, coord in enumerate(rashi_coords):
        ax.plot(np.deg2rad([coord, coord]), [0, 1], "grey", alpha=0.3)
        rashi_text_theta = (
            2*rashi_coords[coord_id] + rashi_extent) * (np.pi/360)
        if (77 < coord < 255):
            if (coord < sun_lambda < coord + rashi_extent):
                ax.text(rashi_text_theta, 0.625,
                        rashi_names[coord_id], rotation=rashi_text_theta*180/np.pi + 180, ha="center", va="center")
            else:
                ax.text(rashi_text_theta, 0.625,
                        rashi_names[coord_id], rotation=rashi_text_theta*180/np.pi + 180, ha="center", va="center", alpha=0.3)
        else:
            if (coord < sun_lambda < coord + rashi_extent):
                ax.text(rashi_text_theta, 0.625,
                        rashi_names[coord_id], rotation=rashi_text_theta*180/np.pi, ha="center", va="center")
            elif ((sun_lambda + 360 - coord) < rashi_extent):
                ax.text(rashi_text_theta, 0.625,
                        rashi_names[coord_id], rotation=rashi_text_theta*180/np.pi, ha="center", va="center")
            else:
                ax.text(rashi_text_theta, 0.625,
                        rashi_names[coord_id], rotation=rashi_text_theta*180/np.pi, ha="center", va="center", alpha=0.3)

    ################ plot other grahas #####################

    ### Inner grahas ###

    mercury_lambda = get_body(
        "mercury", test_date_utc_time).geocentrictrueecliptic.lon.value
    venus_lambda = get_body(
        "venus", test_date_utc_time).geocentrictrueecliptic.lon.value

    if (sun_lambda > mercury_lambda):
        mercury_angle_to_sun = 360 - (sun_lambda - mercury_lambda)
    else:
        mercury_angle_to_sun = (mercury_lambda - sun_lambda)

    if (sun_lambda > venus_lambda):
        venus_angle_to_sun = 360 - (sun_lambda - venus_lambda)
    else:
        venus_angle_to_sun = (venus_lambda - sun_lambda)

    plot_inner_graha_phase("budha", mercury_angle_to_sun, [
                           np.deg2rad(mercury_lambda), 1.9], 0.05, fig, ax)
    plot_inner_graha_phase("shukra", venus_angle_to_sun, [
                           np.deg2rad(venus_lambda), 1.9], 0.05, fig, ax)

    ### Outer grahas ###

    mars_lambda = get_body(
        "mars", test_date_utc_time).geocentrictrueecliptic.lon.value
    jupiter_lambda = get_body(
        "jupiter", test_date_utc_time).geocentrictrueecliptic.lon.value
    saturn_lambda = get_body(
        "saturn", test_date_utc_time).geocentrictrueecliptic.lon.value

    plot_outer_graha("mangala", [np.deg2rad(mars_lambda), 1.9], 0.05, fig, ax)
    plot_outer_graha("guru", [np.deg2rad(jupiter_lambda), 1.9], 0.05, fig, ax)
    plot_outer_graha("shani", [np.deg2rad(saturn_lambda), 1.9], 0.05, fig, ax)

    ##  rahu and ketu  ##

    # This eclipse happens when moon was in Rahu's postion
    known_eclipse_time = Time("2023-04-20 04:16:49", format="iso", scale="utc")

    moon_lambda_known_eclipse = get_body(
        "moon", known_eclipse_time).geocentrictrueecliptic.lon.value
    rahu_lambda, ketu_lambda = calc_rahu_ketu_pos(
        known_eclipse_time, test_date_utc_time, moon_lambda_known_eclipse)

    rahu_lambda, ketu_lambda = ketu_lambda, rahu_lambda

    plot_rahu([np.deg2rad(rahu_lambda), 1.9], 5, fig, ax)
    plot_ketu([np.deg2rad(ketu_lambda), 1.9], 5, fig, ax)

    ################ legend for grahas #####################

    inv = ax.transData.inverted()

    sun_legend_pos = inv.transform((750, 450))
    sun_label_pos = inv.transform((790, 445))
    moon_legend_pos = inv.transform((750, 400))
    moon_label_pos = inv.transform((790, 395))
    mercury_legend_pos = inv.transform((750, 350))
    mercury_label_pos = inv.transform((790, 345))
    venus_legend_pos = inv.transform((750, 300))
    venus_label_pos = inv.transform((790, 295))
    mars_legend_pos = inv.transform((750, 250))
    mars_label_pos = inv.transform((790, 245))
    jupiter_legend_pos = inv.transform((750, 200))
    jupiter_label_pos = inv.transform((790, 195))
    saturn_legend_pos = inv.transform((750, 150))
    saturn_label_pos = inv.transform((790, 145))
    rahu_legend_pos = inv.transform((750, 100))
    rahu_label_pos = inv.transform((790, 95))
    ketu_legend_pos = inv.transform((750, 50))
    ketu_label_pos = inv.transform((790, 45))

    plot_sun(sun_legend_pos, 0.1, fig, ax)
    plt.text(sun_label_pos[0], sun_label_pos[1], "\sam{सूर्यः } ")
    plot_moon_phase(final_tithi, moon_legend_pos, 0.1, fig, ax)
    plt.text(moon_label_pos[0], moon_label_pos[1], "\sam{चन्द्रः } ")
    plot_inner_graha_phase("budha", mercury_angle_to_sun,
                           mercury_legend_pos, 0.05, fig, ax)
    plt.text(mercury_label_pos[0], mercury_label_pos[1], "\sam{बुधः } ")
    plot_inner_graha_phase("shukra", venus_angle_to_sun,
                           venus_legend_pos, 0.05, fig, ax)
    plt.text(venus_label_pos[0], venus_label_pos[1], "\sam{शुक्रः } ")
    plot_outer_graha("mangala", mars_legend_pos, 0.05, fig, ax)
    plt.text(mars_label_pos[0], mars_label_pos[1], "\sam{मङ्गलः } ")
    plot_outer_graha("guru", jupiter_legend_pos, 0.05, fig, ax)
    plt.text(jupiter_label_pos[0], jupiter_label_pos[1], "\sam{बृहस्पतिः } ")
    plot_outer_graha("shani", saturn_legend_pos, 0.05, fig, ax)
    plt.text(saturn_label_pos[0], saturn_label_pos[1], "\sam{शनैश्चरः } ")

    plot_rahu(rahu_legend_pos, 5, fig, ax)
    plt.text(rahu_label_pos[0], rahu_label_pos[1], "\sam{राहुः} ")
    plot_ketu(ketu_legend_pos, 5, fig, ax)
    plt.text(ketu_label_pos[0], ketu_label_pos[1], "\sam{केतुः} ")

    ################ legend for panchanga #####################

    plt.text(inv.transform((-200, 500))[0],
             inv.transform((-200, 400))[1], '\sam{'+location+'}')

    plt.text(inv.transform((-200, 470))
             [0], inv.transform((-200, 370))[1], '\sam{'+date_str.split(" ")[0]+'}')
    plt.text(inv.transform((-200, 440))
             [0], inv.transform((-200, 340))[1], '\sam{'+date_str.split(" ")[1]+'}')

    plt.text(inv.transform((-150, 340))[0], inv.transform((-150, 340))[
             1], "\sam{पञ्चाङ्ग }", bbox=dict(facecolor='none', edgecolor='white'))

    tithi_names_tab = Table.read('tithi_names.tex', format="latex")
    tithi_number = tithi_names_tab['tithi'].data
    tithi_names = tithi_names_tab['tithi_names'].data

    plt.text(inv.transform((-200, 310))
             [0], inv.transform((-200, 310))[1], tithi_names[int(final_tithi)])

    vaara_names_tab = Table.read('vaara_names.tex', format="latex")
    weekday = vaara_names_tab['weekday'].data
    vaara_names = vaara_names_tab['vaara'].data

    plt.text(inv.transform((-200, 280))[0], inv.transform((-200, 280))[
             1], vaara_names[np.where(weekday == day_of_the_week_at_t)][0])

    plt.text(inv.transform((-200, 250))
             [0], inv.transform((-200, 250))[1], final_nakshatram)

    plt.text(inv.transform((-200, 220))
             [0], inv.transform((-200, 220))[1], final_yoga)

    plt.title(r'\sam{ॐ }', fontsize=20)
    plt.ylim(0, 2.2)

    plt.savefig(filename)
    plt.close()

    return 0

# location = "Bangalore, India"
# date_str = "2023-06-20 16:30:00"

# location = "Hershey, Pennsylvania, USA"
# date_str = "2023-07-18 22:19:55"


location = "Bengaluru, India"
date_str = "2024-12-02 15:30:00"

calc_nakshatra_tithi(location, date_str, filename="nakshatra_at_test_time.pdf")
