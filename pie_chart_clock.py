#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 18:59:44 2025

@author: hannah

makes pie charts of clock angle distribution
"""

import matplotlib.pyplot as plt
import numpy as np
#import spiceypy as sp
import math
import pandas as pd
import datetime as dt
from datetime import datetime
import matplotlib.dates as mdates
from scipy.optimize import curve_fit
import matplotlib.ticker as plticker
from pyestimate import sin_param_estimate
import matplotlib.ticker
from matplotlib.patches import Patch


root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'
root_saves = '/Users/hannah/OneDrive - Lancaster University/aurora/'

# pie chart function

# data for each visit
visit_02_data = pd.read_csv(root_folder+'visit_02_data_aug.csv', delimiter=',')
visit_03_data = pd.read_csv(root_folder+'visit_03_data_aug.csv', delimiter=',')
visit_04_data = pd.read_csv(root_folder+'visit_04_data_aug.csv', delimiter=',')
visit_05_data = pd.read_csv(root_folder+'visit_05_data_aug.csv', delimiter=',')
visit_08_data = pd.read_csv(root_folder+'visit_08_data_aug.csv', delimiter=',')
visit_09_data = pd.read_csv(root_folder+'visit_09_data_aug.csv', delimiter=',')
visit_10_data = pd.read_csv(root_folder+'visit_10_data_aug.csv', delimiter=',')
visit_11_data = pd.read_csv(root_folder+'visit_11_data_aug.csv', delimiter=',')
visit_16_data = pd.read_csv(root_folder+'visit_16_data_aug.csv', delimiter=',')
visit_17_data = pd.read_csv(root_folder+'visit_17_data_aug.csv', delimiter=',')
visit_18_data = pd.read_csv(root_folder+'visit_18_data_aug.csv', delimiter=',')
visit_19_data = pd.read_csv(root_folder+'visit_19_data_aug.csv', delimiter=',')
visit_20_data = pd.read_csv(root_folder+'visit_20_data_aug.csv', delimiter=',')
visit_21_data = pd.read_csv(root_folder+'visit_21_data_aug.csv', delimiter=',')
visit_24_data = pd.read_csv(root_folder+'visit_24_data_aug.csv', delimiter=',')
visit_25_data = pd.read_csv(root_folder+'visit_25_data_aug.csv', delimiter=',')
visit_27_data = pd.read_csv(root_folder+'visit_27_data_aug.csv', delimiter=',')
visit_28_data = pd.read_csv(root_folder+'visit_28_data_aug.csv', delimiter=',')
visit_34_data = pd.read_csv(root_folder+'visit_34_data_aug.csv', delimiter=',')
visit_35_data = pd.read_csv(root_folder+'visit_35_data_aug.csv', delimiter=',')


# # high cml - 01, 12, 15, 26
visit_01_data = pd.read_csv(root_folder+'visit_01_data_aug.csv', delimiter=',')
visit_12_data = pd.read_csv(root_folder+'visit_12_data_aug.csv', delimiter=',')
visit_15_data = pd.read_csv(root_folder+'visit_15_data_aug.csv', delimiter=',')
visit_26_data = pd.read_csv(root_folder+'visit_26_data_aug.csv', delimiter=',')


# extract clock angle for each visit
clock_02 = visit_02_data['Clock_Angle'].to_numpy()
clock_03 = visit_03_data['Clock_Angle'].to_numpy()
clock_04 = visit_04_data['Clock_Angle'].to_numpy()
clock_05 = visit_05_data['Clock_Angle'].to_numpy()
clock_08 = visit_08_data['Clock_Angle'].to_numpy()
clock_09 = visit_09_data['Clock_Angle'].to_numpy()
clock_10 = visit_10_data['Clock_Angle'].to_numpy()
clock_11 = visit_11_data['Clock_Angle'].to_numpy()
clock_16 = visit_16_data['Clock_Angle'].to_numpy()
clock_17 = visit_17_data['Clock_Angle'].to_numpy()
clock_18 = visit_18_data['Clock_Angle'].to_numpy()
clock_19 = visit_19_data['Clock_Angle'].to_numpy()
clock_20 = visit_20_data['Clock_Angle'].to_numpy()
clock_21 = visit_21_data['Clock_Angle'].to_numpy()
clock_24 = visit_24_data['Clock_Angle'].to_numpy()
clock_25 = visit_25_data['Clock_Angle'].to_numpy()
clock_27 = visit_27_data['Clock_Angle'].to_numpy()
clock_28 = visit_28_data['Clock_Angle'].to_numpy()
clock_34 = visit_34_data['Clock_Angle'].to_numpy()
clock_35 = visit_35_data['Clock_Angle'].to_numpy()

# clock_13 = visit_13_data['Clock_Angle'].to_numpy()

# high cml
clock_01 = visit_01_data['Clock_Angle'].to_numpy()
clock_12 = visit_12_data['Clock_Angle'].to_numpy()
clock_15 = visit_15_data['Clock_Angle'].to_numpy()
clock_26 = visit_26_data['Clock_Angle'].to_numpy()

clocks = np.concatenate((clock_01,clock_02,clock_03,clock_04,clock_05,clock_08,clock_09,clock_10,clock_11,clock_12,clock_15,clock_16,clock_17,clock_18,clock_19,clock_20,clock_21,clock_24,clock_25,clock_26,clock_27,clock_28,clock_34,clock_35))

# --------- 

Bz_pos = []
Bz_neg = []
By_pos = []
By_neg = []

for angle in range(len(clocks)):
    if clocks[angle] > -45 and clocks[angle] < 45:
        Bz_pos.append(clocks[angle])
    
    elif clocks[angle] > 45 and clocks[angle] < 135:
        By_pos.append(clocks[angle])

    elif clocks[angle] < -45 and clocks[angle] > -135:
        By_neg.append(clocks[angle])
        
    else:
        Bz_neg.append(clocks[angle])
        
clock_dist = np.array([len(Bz_pos),len(By_pos),len(Bz_neg),len(By_neg)])
clock_labels = ['+Bz','+By  ','-Bz  ',' -By']

# plt.pie(clock_dist, labels = clock_labels)
# plt.show() 

# -------------


df = pd.read_csv(root_folder+'juno_data_big_df_ionotime_july.csv')
df_155 = pd.read_csv(root_folder+'juno_data_doy_155_ionotime_july.csv')
df_156_plus = pd.read_csv(root_folder+'juno_data_doy_156_plus_ionotime_july.csv')
df_175 = pd.read_csv(root_folder+'juno_data_doy_175_ionotime_july.csv')

clock_angle = df['Clock_Angle'].to_numpy()
clock_angle_155 = df_155['Clock_Angle'].to_numpy()
clock_angle_156_plus = df_156_plus['Clock_Angle'].to_numpy()
clock_angle_175 = df_175['Clock_Angle'].to_numpy()

clocks_whole = np.concatenate((clock_angle,clock_angle_155,clock_angle_156_plus,clock_angle_175))

Bz_pos_w = []
Bz_neg_w = []
By_pos_w = []
By_neg_w = []

for angle_w in range(len(clocks_whole)):
    if clocks_whole[angle_w] > -45 and clocks_whole[angle_w] < 45:
        Bz_pos_w.append(clocks_whole[angle])
    
    elif clocks_whole[angle_w] > 45 and clocks_whole[angle_w] < 135:
        By_pos_w.append(clocks_whole[angle_w])

    elif clocks_whole[angle_w] < -45 and clocks_whole[angle_w] > -135:
        By_neg_w.append(clocks_whole[angle_w])
        
    else:
        Bz_neg_w.append(clocks_whole[angle_w])
        
clock_dist_w = np.array([len(Bz_pos_w),len(By_pos_w),len(Bz_neg_w),len(By_neg_w)])
clock_labels_w = ['+Bz','+By  ','-Bz  ',' -By ']

# ---------

clock_02_median = np.median(clock_02)
clock_03_median = np.median(clock_03)
clock_04_median = np.median(clock_04)
clock_05_median = np.median(clock_05)
clock_08_median = np.median(clock_08)
clock_09_median = np.median(clock_09)
clock_10_median = np.median(clock_10)
clock_11_median = np.median(clock_11)
clock_16_median = np.median(clock_16)
clock_17_median = np.median(clock_17)
clock_18_median = np.median(clock_18)
clock_19_median = np.median(clock_19)
clock_20_median = np.median(clock_20)
clock_21_median = np.median(clock_21)
clock_24_median = np.median(clock_24)
clock_25_median = np.median(clock_25)
clock_27_median = np.median(clock_27)
clock_28_median = np.median(clock_28)
clock_34_median = np.median(clock_34)
clock_35_median = np.median(clock_35)

# hiigh cmls
clock_01_median = np.median(clock_01)
clock_12_median = np.median(clock_12)
clock_15_median = np.median(clock_15)
clock_26_median = np.median(clock_26)

#clock_medians = [clock_02_median,clock_03_median,clock_04_median, clock_05_median, clock_08_median, clock_09_median, clock_10_median,clock_11_median,clock_16_median, clock_17_median, clock_18_median, clock_19_median, clock_20_median, clock_21_median,clock_24_median,clock_25_median,clock_27_median,clock_28_median,clock_34_median,clock_35_median]
clock_medians = [clock_01_median,clock_02_median,clock_03_median,clock_04_median, clock_05_median, clock_08_median, clock_09_median,clock_10_median,clock_11_median,clock_12_median,clock_15_median,clock_16_median, clock_17_median, clock_18_median, clock_19_median, clock_20_median, clock_21_median,clock_24_median,clock_25_median,clock_26_median,clock_27_median,clock_28_median,clock_34_median,clock_35_median]
clock_medians = np.array(clock_medians)


'''
make sure to check which clock angles will need converting - any in which the median/mean causes large/small angles
'''

def clock_format_convert(clock_angle):
    new_clock_angle = []
    for i in range(len(clock_angle)):
        if clock_angle[i] < 0:
            clock = clock_angle[i] + 360
            new_clock_angle.append(clock)
        else:
            new_clock_angle.append(clock_angle[i])
    return new_clock_angle


new_clock_24 = clock_format_convert(clock_24)

new_clock_24_median = np.median(new_clock_24)
# # #new_clock_24_median = new_clock_24_median_1 - 360
clock_medians[17] = new_clock_24_median

Bz_pos_v = []
Bz_neg_v = []
By_pos_v = []
By_neg_v = []

for angle_v in range(len(clock_medians)):
    if clock_medians[angle_v] > -45 and clock_medians[angle_v] < 45:
        Bz_pos_v.append(clocks[angle_v])
    
    elif clock_medians[angle_v] > 45 and clock_medians[angle_v] < 135:
        By_pos_v.append(clocks[angle_v])

    elif clock_medians[angle_v] < -45 and clock_medians[angle_v] > -135:
        By_neg_v.append(clocks[angle_v])
        
    else:
        Bz_neg_v.append(clocks_whole[angle_v])
        
clock_dist_v = np.array([len(Bz_pos_v),len(By_pos_v),len(Bz_neg_v),len(By_neg_v)])
clock_labels_v = ['+Bz','+By  ','-Bz  ',' -By ']

# ---------



# fig = plt.figure(figsize=(10,30))
# ax1 = plt.subplot(1,3,1)
# patches, labels, percentages = ax1.pie(clock_dist_w, labels = clock_labels_w, textprops={'fontsize': 26},autopct='%1.1f%%',labeldistance=1.15,colors=['red', 'green', 'orange', 'blue'],wedgeprops = { 'linewidth' : 2, 'edgecolor' : 'white' },startangle=60)
# ax1.text(-1.5,0.98,'a',style='italic',fontsize=50)
# for label, percentage in zip(labels, percentages):
#     label.set_text(label.get_text() + '\n' + percentage.get_text())
#     percentage.remove()
# plt.subplots_adjust(bottom=0.1, right=2.5, top=0.3)
# ax2 = plt.subplot(1,3,2)
# patches, labels, percentages = ax2.pie(clock_dist, labels = clock_labels, textprops={'fontsize': 26},autopct='%1.1f%%', labeldistance=1.15,colors=['red', 'green', 'orange', 'blue'],wedgeprops = { 'linewidth' : 2, 'edgecolor' : 'white' },startangle=75)
# ax2.text(-1.5,0.96,'b',style='italic',fontsize=50)
# for label, percentage in zip(labels, percentages):
#     label.set_text(label.get_text() + '\n' + percentage.get_text())
#     percentage.remove()
# ax3 = plt.subplot(1,3,3)
# patches, labels, percentages = ax3.pie(clock_dist_v, labels = clock_labels_v, textprops={'fontsize': 26},autopct='%1.1f%%', labeldistance=1.15,colors=['red', 'green', 'orange', 'blue'],wedgeprops = { 'linewidth' : 2, 'edgecolor' : 'white' },startangle=75)
# ax3.text(-1.5,0.96,'c',style='italic',fontsize=50)
# for label, percentage in zip(labels, percentages):
#     label.set_text(label.get_text() + '\n' + percentage.get_text())
#     percentage.remove()


fig = plt.figure(figsize=(15, 10))

# Uniform pie chart size
pie_width = 0.25
pie_height = 0.4
y_top = 0.50
y_bottom = 0.07

# Carefully chosen positions for balance
ax1 = fig.add_axes([0.2, y_top, pie_width, pie_height])    # top left
ax2 = fig.add_axes([0.55, y_top, pie_width, pie_height])   # top right
ax3 = fig.add_axes([0.375, y_bottom, pie_width, pie_height])  # bottom center

# --- PIE 1 ---
patches, labels, percentages = ax1.pie(clock_dist_w, labels=clock_labels_w,
    textprops={'fontsize': 16}, autopct='%1.1f%%', labeldistance=1.15,
    colors=['red', 'green', 'orange', 'blue'],
    wedgeprops={'linewidth': 2, 'edgecolor': 'white'}, startangle=60)
ax1.text(-1.3, 1.2, 'a', style='italic', fontsize=30)
for label, percentage in zip(labels, percentages):
    label.set_text(label.get_text() + '\n' + percentage.get_text())
    percentage.remove()

# --- PIE 2 ---
patches, labels, percentages = ax2.pie(clock_dist, labels=clock_labels,
    textprops={'fontsize': 16}, autopct='%1.1f%%', labeldistance=1.15,
    colors=['red', 'green', 'orange', 'blue'],
    wedgeprops={'linewidth': 2, 'edgecolor': 'white'}, startangle=75)
ax2.text(-1.3, 1.2, 'b', style='italic', fontsize=30)
for label, percentage in zip(labels, percentages):
    label.set_text(label.get_text() + '\n' + percentage.get_text())
    percentage.remove()

# --- PIE 3 ---
patches, labels, percentages = ax3.pie(clock_dist_v, labels=clock_labels_v,
    textprops={'fontsize': 16}, autopct='%1.1f%%', labeldistance=1.15,
    colors=['red', 'green', 'orange', 'blue'],
    wedgeprops={'linewidth': 2, 'edgecolor': 'white'}, startangle=75)
ax3.text(-1.3, 1.2, 'c', style='italic', fontsize=30)
for label, percentage in zip(labels, percentages):
    label.set_text(label.get_text() + '\n' + percentage.get_text())
    percentage.remove()


'''
# save plot
saveloc = (f'{root_saves}clock_angle_pie_sept.jpg') 
plt.savefig(saveloc,bbox_inches='tight',dpi=400)
'''

centres_deg = [0, 90, 180, -90]


# def quadrant_circular_bar(ax, counts, labels, colors, show_percent=True):
#     centres_deg = [0, 90, 180, -90]
#     theta = np.deg2rad(centres_deg)
#     width = np.deg2rad(90)

#     # Draw bars
#     ax.bar(
#         theta,
#         counts,
#         width=width,
#         bottom=0,
#         color=colors,
#         edgecolor='white',
#         linewidth=2.5
#     )

#     # Clock-angle formatting
#     ax.set_theta_zero_location('N')
#     ax.set_theta_direction(-1)

#     # Percentages with labels
#     if show_percent:
#         total = sum(counts)
#         labels_with_percent = [
#             f"{lbl}\n{counts[i]/total*100:.1f}%" for i, lbl in enumerate(labels)
#         ]
#     else:
#         labels_with_percent = labels
        
#     for angle, label in zip(ax.get_xticks(), ax.get_xticklabels()):
#         angle_deg = np.rad2deg(angle)
#         if angle_deg in [0, 180]:  # left/right
#             label.set_y(label.get_position()[1] - 0.01)
#         elif angle_deg in [90,270]:
#             label.set_y(label.get_position()[1] - 0.5)
#     #     #else:                       # top/bottom
#     #        # label.set_y(label.get_position()[1] - 0.001)


#     # Set the tick labels with percentages
#     ax.set_thetagrids([0, 90, 180, 270], labels_with_percent, fontsize=25)
    
#     #ax.set_yticklabels([])   # hide numbers
#     ax.yaxis.grid(True, color='gray', linestyle='--', linewidth=0.8)  # keep grid lines
    
#     # Optional: adjust tick label padding to move them further out
#     ax.tick_params(axis='x', pad=25)

def quadrant_circular_bar(ax, counts, labels, colors, show_percent=True):
    centres_deg = [0, 90, 180, 270]  # top, right, bottom, left
    theta = np.deg2rad(centres_deg)
    width = np.deg2rad(90)

    # Draw bars
    ax.bar(
        theta,
        counts,
        width=width,
        bottom=0,
        color=colors,
        edgecolor='white',
        linewidth=2.5
    )

    # Clock-angle formatting
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)

    # Prepare labels with percentages
    if show_percent:
        total = sum(counts)
        labels_with_percent = [
            f"{lbl}\n{counts[i]/total*100:.1f}%" for i, lbl in enumerate(labels)
        ]
    else:
        labels_with_percent = labels

    # Set the tick labels first
    ax.set_thetagrids([0, 90, 180, 270], labels_with_percent, fontsize=25)

    # Hide radial ticks if you want
    ax.set_yticklabels([])

    # Keep the grid lines
    ax.yaxis.grid(True, color='gray', linestyle='--', linewidth=0.8)

    # Adjust padding individually for each label using tick_params
    # top (0°) and bottom (180°) smaller pad, left (270°) and right (90°) larger
    ticks = ax.get_xticklabels()
    for i, label in enumerate(ticks):
        angle = centres_deg[i]
        if angle in [0, 180]:          # top/bottom
            label.set_y(label.get_position()[1] - 0.12)  # small nudge out
        elif angle in [90, 270]:       # right/left
            label.set_y(label.get_position()[1] - 0.22)  # bigger nudge out





# fig = plt.figure(figsize=(12, 10))
# plt.tight_layout(pad=3.0)

# # Carefully chosen positions for balance
# ax1 = fig.add_axes([0.2, y_top, pie_width, pie_height])    # top left
# ax2 = fig.add_axes([0.55, y_top, pie_width, pie_height])   # top right
# ax3 = fig.add_axes([0.375, y_bottom, pie_width, pie_height])  # bottom center

# ax1 = fig.add_subplot(2, 2, 1, polar=True)
# ax2 = fig.add_subplot(2, 2, 2, polar=True)
# ax3 = fig.add_subplot(2, 2, 3, polar=True)


# colors = ['red', 'green', 'orange', 'blue']

# quadrant_circular_bar(ax1, clock_dist_w, clock_labels_w, colors)
# ax1.text(-1.3, 1.2, 'a', style='italic', fontsize=30)

# quadrant_circular_bar(ax2, clock_dist, clock_labels, colors)
# ax2.text(-1.3, 1.2, 'b', style='italic', fontsize=30)

# quadrant_circular_bar(ax3, clock_dist_v, clock_labels_v, colors)
# ax3.text(-1.3, 1.2, 'c', style='italic', fontsize=30)

# plt.tight_layout()
# plt.show()

fig = plt.figure(figsize=(20, 10))

# Keep pies the same size
pie_height = 0.38

# Top row (a & b)
y_top = 0.6
ax1 = fig.add_axes([0.18, y_top, pie_width, pie_height], polar=True)
ax2 = fig.add_axes([0.55, y_top, pie_width, pie_height], polar=True)

# Bottom row (c) — move lower to create more vertical gap
y_bottom = 0.0001  # was 0.0001
ax3 = fig.add_axes([0.365, y_bottom, pie_width, pie_height], polar=True)

# # Manually positioned polar axes
# ax1 = fig.add_axes([0.18, y_top, pie_width, pie_height], polar=True)
# ax2 = fig.add_axes([0.55, y_top, pie_width, pie_height], polar=True)
# ax3 = fig.add_axes([0.365, 0.05, pie_width, pie_height], polar=True)

colors = ['red', 'green', 'orange', 'blue']

quadrant_circular_bar(ax1, clock_dist_w, clock_labels_w, colors)
ax1.text(-0.15, 1.12, 'a', transform=ax1.transAxes,
         fontsize=40, fontstyle='italic')

quadrant_circular_bar(ax2, clock_dist, clock_labels, colors)
ax2.text(-0.15, 1.12, 'b', transform=ax2.transAxes,
         fontsize=40, fontstyle='italic')

quadrant_circular_bar(ax3, clock_dist_v, clock_labels_v, colors)
ax3.text(-0.15, 1.12, 'c', transform=ax3.transAxes,
         fontsize=40, fontstyle='italic')

#plt.show()

# save plot
saveloc = (f'{root_saves}clock_angle_hist.jpg') 
plt.savefig(saveloc,bbox_inches='tight',dpi=400)
