'''
This script plots failure envelopes for a given set of composite
material properties by calling the fail_initiation Fortran subroutine
which implements failure criteria developed by Catalanotti et al.
(see Fortran subroutine for further details).
Copyright (C) 2021 Rutger Wouter Kok

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA

'''
import numpy as np
from catalanotti import fail_initiation, rotate_stress
import matplotlib.pyplot as plt

# define the number of points
n = 250

# In-situ ply properties IM7-8552
# XT = 2.3235  # tensile strength fiber direction
# XC = 1.2001  # compressive strength fiber direction
# YT = 0.1602  # tensile strength transverse direction
# YC = 0.198  # compressive strength transverse direction
# SL = 0.1302  # shear strength
# g12 = 5.290  # shear modulus
# etaL = 0.5
# alpha0 = np.radians(53.0)

# Calculate parameters needed for failure criteria calculation
phiC = np.arctan((1.0 - np.sqrt(1.0 - 4.0 * (SL / XC) * ((SL / XC) + etaL))) /
                 (2.0 * ((SL / XC) + etaL)))  # Eq. 72
etaT = (etaL * ST) / SL  # Eq.10
kappa = (ST**2.0 - YT**2.0) / (ST * YT)  # Eq. 43
lmbda = ((2.0 * etaL * ST) / SL) - kappa  # Eq. 45

# Create jump tables to index stresses (Abaqus VUMAT stress vector format)
idx_jump_table = {11: 0, 22: 1, 33: 2, 12: 3, 23: 4, 13: 5}
# Create jump table to define stress ranges for different loadcases
strength_table = {11: [-XC, XT], 22: [-YC, YT], 33: [-YC, YT],
                  12: [0.0, 2 * SL], 13: [0.0, 2 * SL], 23: [0.0, ST]}


def create_envelope(dir1, dir2):
    '''
    Creates an array of coordinates defining the failure envelope of a
    composite material in the dir1-dir2 stress space.
    dir1 = stress to be plotted along x-axis (11,22,33,12,23,13)
    dir2 = stress to be plotted along y-axis (11,22,33,12,23,13)
    '''
    idx1 = idx_jump_table[dir1]  # index of dir1 stress in stress vector
    idx2 = idx_jump_table[dir2]  # index of dir2 stress in stress vector
    tol = 0.005  # tolerance to ensure the full stress range is plotted
    min_stress_1 = strength_table[dir1][0] - tol  # min. stress dir1
    max_stress_1 = strength_table[dir1][1] + tol  # max. stress dir1
    min_stress_2 = strength_table[dir2][0]  # min. stress dir2
    max_stress_2 = strength_table[dir2][1]  # max. stress dir2
    stress_inc_1 = (max_stress_1 - min_stress_1) / float(n)  # stress inc. dir1
    stress_inc_2 = (max_stress_2 - min_stress_2) / float(n)  # stress inc. dir2
    trial_stress = np.zeros((6, 1))  # initialize trial_stress vector
    failure_coords = np.zeros((n, 2))  # initialize failure coords array
    for i in range(n):
        for j in range(n):
            # assign trial stresses
            trial_stress[idx1] = min_stress_1 + i * stress_inc_1
            trial_stress[idx2] = min_stress_2 + j * stress_inc_2
            # initialize failure indices
            fi_lc, fi_m, fi_lt = 0.0, 0.0, 0.0
            if trial_stress[0] > 0.0:  # longitudinal tensile failure
                fi_lt = trial_stress[0] / XT
            elif trial_stress[0] < 0.0:  # longitudinal compressive failure
                # rotate stresses to misalignment frame for long. compression
                trial_stressP = rotate_stress(trial_stress, phiC, XC, g12)
                fi_lc, a = fail_initiation(trial_stressP, ST, SL, etaL,
                                           etaT, lmbda, kappa)
            # matrix failure
            fi_m, a = fail_initiation(trial_stress, ST, SL, etaL, etaT,
                                      lmbda, kappa)
            # if any failure index > 1 append to failure coords and break
            if (fi_m > 0.99) or (fi_lc > 0.99) or (fi_lt > 0.99):
                failure_coords[i] = trial_stress[idx1], trial_stress[idx2]
                break
            print('Add a new comment')
    return failure_coords


def plot_envelope(dir1, dir2, failure_coords):
    '''
    Plots failure envelopes.
    dir1 = stress to be plotted along x-axis (11,22,33,12,23,13)
    dir2 = stress to be plotted along y-axis (11,22,33,12,23,13)
    failure_coords = stress-stress array to be plotted
    '''
    fig, ax = plt.subplots()  # Create a figure and an axes.
    x = failure_coords[:, 0]
    y = failure_coords[:, 1]
    ax.plot(x, y)  # Plot data on the axes.
    ax.set_xlabel(str(dir1))  # Add an x-label to the axes.
    ax.set_ylabel(str(dir2))  # Add a y-label to the axes.
    plt.grid()
    plt.show()


if __name__ == '__main__':
    # 11-22 envelope not plotted since stress in 22 direction ranges from
    # -YC to +YT
    envelope1 = create_envelope(22, 12)
    plot_envelope(22, 12, envelope1)
    envelope2 = create_envelope(22, 13)
    plot_envelope(22, 13, envelope2)
    envelope3 = create_envelope(22, 23)
    plot_envelope(22, 23, envelope3)
    envelope4 = create_envelope(11, 12)
    plot_envelope(11, 12, envelope4)
    envelope5 = create_envelope(11, 23)
    plot_envelope(11, 23, envelope5)
    envelope6 = create_envelope(12, 13)
    plot_envelope(12, 13, envelope6)
    envelope7 = create_envelope(12, 23)
    plot_envelope(12, 23, envelope7)
