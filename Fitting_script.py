# -*- coding: utf-8 -*-
"""
Created on Mon May 22 08:11:10 2023

@author: manue
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from brainpy import isotopic_variants
import csv
import scipy.interpolate
from sklearn.metrics import mean_absolute_error as mae

def plot_experimental_interpolated(x_exp, y_exp,y_theo,theoretical_isotopic_pattern_x, theoretical_isotopic_pattern_y,molecule,mz_min, mz_max,inten_min=0, inten_max=1.01,lw=1,tick_spacing = 2):
    molecule_str = str(molecule)
    molecule_str = molecule_str.replace('{', '').replace('}', '').replace(':', '_')
    molecule_composition = ' '.join(f'{element}:{count}' for element, count in molecule.items())

    fig, ax = plt.subplots(figsize=(6,4))    
    
    ax.plot(x_exp, y_exp, 'o',color="green", label='Experimental data')
    
    ax.plot(x_exp, y_theo,lw=lw,label='Interpolated data')
        
    
    ax.plot(x_exp,y_exp,"green",lw=lw,zorder=1)
    markerline, stemline, baseline, =ax.stem(theoretical_isotopic_pattern_x,theoretical_isotopic_pattern_y,"black",label='Isotopic pattern', markerfmt='ko', basefmt=" ",use_line_collection=True)
   

    plt.setp(markerline, markersize = 8,color="red")
    plt.setp(stemline, linewidth = 2,color="black")
    
    
    ax.set_xlabel('m/z',fontsize=16)
    ax.set_ylabel('Relative intensity',fontsize=16)
    ax.set_xlim([mz_min, mz_max])
    ax.set_ylim([inten_min, inten_max])
    ax.set_title(f' {molecule_composition}', fontsize=16)

    ax.legend()
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    fig.tight_layout()
    
    fig.savefig(molecule_str.replace('_raw.csv', '_ciu50_iwat.png'), dpi=300)
    print('Plot saved to '+molecule_str.replace('_raw.csv', '_ciu50_iwat.png'))



    

def LeastSquaresFitter(observed, expected):
    exp_max = max(observed)
    theo_max = max(expected)

    sum_of_squared_errors = 0
    sum_of_squared_theoreticals = 0

    for e, t in zip(observed, expected):
        normed_expr = e / exp_max
        normed_theo = t / theo_max
        sum_of_squared_errors += (normed_theo - normed_expr) ** 2
        sum_of_squared_theoreticals += normed_theo ** 2
    
    return sum_of_squared_errors / sum_of_squared_theoreticals
    
    
def theo_iso(peptide,npeaks=15,charge=5,continous=False):

    # Generate theoretical isotopic pattern
    
    theoretical_isotopic_cluster = isotopic_variants(peptide, npeaks, charge)
    if continous:
    # produce a theoretical profile using a gaussian peak shape
    
        mz_grid = np.arange(theoretical_isotopic_cluster[0].mz - 1,
                    theoretical_isotopic_cluster[-1].mz + 1, 0.02)
        intensity = np.zeros_like(mz_grid)
        sigma = 0.00002
        for peak in theoretical_isotopic_cluster:
            # Add gaussian peak shape centered around each theoretical peak
            intensity += peak.intensity * np.exp(-(mz_grid - peak.mz) ** 2 / (2 * sigma)
            ) / (np.sqrt(2 * np.pi) * sigma)

        # Normalize profile to 0-100

        intensity = (intensity / intensity.max()) * 1
        mz=mz_grid
    else:
   
        mz = []
        intensity = []

        for peak in theoretical_isotopic_cluster:
            mz.append(peak.mz)
            intensity.append(peak.intensity)
        
        intensity = pd.DataFrame(intensity)
        intensity = (intensity / intensity.max()) 
        mz = pd.DataFrame(mz)

        #plt.stem(mz[0],intensity[0], bottom=0,basefmt=" ",use_line_collection=True)
        mz =mz[0]
        intensity= intensity[0]
    return mz,intensity


def detected_peaks(x_inset2, y_inset2,x_peaks, y_peaks,name):
        
    fig, ax = plt.subplots(figsize=(6,4))    
    ax.plot(x_inset2, y_inset2, label='Experimental Data')
    ax.plot(x_peaks, y_peaks, 'ro', label='Peaks')
    ax.set_xlabel('m/z',fontsize=18)
    ax.set_ylabel('Relative intensity',fontsize=18)
    ax.set_title('Experimental Data with Identified Peaks',fontsize=20)  
    ax.legend()
    ax.tick_params(axis='both', which='major', labelsize=16)
    fig.tight_layout()
    fig.savefig(name.replace('_raw.csv', '_ciu50_iwat.png'), dpi=300)
    print('Plot saved to '+name.replace('_raw.csv', '_ciu50_iwat.png'))

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def prepare_data(file):
    raw = np.genfromtxt(file,delimiter='')
    plt.plot(raw[0:,0],raw[0:,1])
    return raw[0:,0],raw[0:,1]  ## x,y



def plot_inset(x,y,ylim,name,mz_min,mz_max,peak1_min,peak1_max):
    
    result1 =np.where((x> peak1_min) & (x< peak1_max))
    y_norm=NormalizeData(y[result1])
    
    return x[result1],y_norm



def generate_molecules(h_range, zn_range,cu_range):
    fixed_elements = {
        'C': 259,
        'N': 77,
        'O': 102,
        'S': 21,
    }

    molecules = []
    for h in range(h_range[0], h_range[1] + 1):
        for zn in range(zn_range[0], zn_range[1] + 1):
            for cu in range(cu_range[0], cu_range[1] + 1):
                composition = fixed_elements.copy()
                composition['H'] = h
                composition['Zn'] = zn
                composition['Cu'] = cu
                molecules.append(composition)

    return molecules


def calculate_r2(x_exp, y_exp, theoretical_isotopic_pattern_x, theoretical_isotopic_pattern_y, molecule, mz_min, mz_max, inten_min, inten_max, intensity_threshold=0.2,tick_spacing=5,plot=False):
    
    # Filter theoretical isotopic pattern based on intensity threshold
    filtered_theoretical_indices = np.where(theoretical_isotopic_pattern_y > intensity_threshold)[0]
    filtered_theoretical_x = theoretical_isotopic_pattern_x[filtered_theoretical_indices]
    filtered_theoretical_y = theoretical_isotopic_pattern_y[filtered_theoretical_indices]

    # Interpolate the filtered theoretical isotopic pattern onto the x-values of the experimental data
    y_theo = np.interp(x_exp, filtered_theoretical_x, filtered_theoretical_y)

    # Calculate the R2 metric
    ss_residual = np.sum((y_exp - y_theo) ** 2)
    ss_total = np.sum((y_exp - np.mean(y_exp)) ** 2)
    r2 = 1 - (ss_residual / ss_total)
    if plot:
        plot_experimental_interpolated(x_exp, y_exp, y_theo, filtered_theoretical_x, filtered_theoretical_y, molecule, mz_min=mz_min, mz_max=mz_max, inten_min=inten_min, inten_max=inten_max, lw=2, tick_spacing=tick_spacing)

    return r2



def calculate_chi2(x_exp, y_exp, theoretical_isotopic_pattern_x , theoretical_isotopic_pattern_y, intensity_threshold=0.2):
    # Interpolate the theoretical spectrum onto the x-values of the experimental spectrum
       # Filter theoretical isotopic pattern based on intensity threshold
    filtered_theoretical_indices = np.where(theoretical_isotopic_pattern_y > intensity_threshold)[0]
    filtered_theoretical_x = theoretical_isotopic_pattern_x[filtered_theoretical_indices]
    filtered_theoretical_y = theoretical_isotopic_pattern_y[filtered_theoretical_indices]

    # Interpolate the filtered theoretical isotopic pattern onto the x-values of the experimental data
    y_theo_interp = np.interp(x_exp, filtered_theoretical_x, filtered_theoretical_y)


    # Calculate the chi-square value
    chi2 = np.sum((y_exp - y_theo_interp) ** 2 / (y_theo_interp+ 1e-9))

    return chi2

def match_peaks(peaks_sim, peaks_trace, mean_peak_dist):
    peaks_sim = list(peaks_sim)
    matched_peaks = pd.DataFrame({"trace" : peaks_trace, "sim" : np.NaN})
    for idx, row in matched_peaks.iterrows():
        if len(peaks_sim) > 1:
            idx_closest = np.argmin(np.abs(row["trace"] - peaks_sim))
            if np.abs(row["trace"] - peaks_sim[idx_closest]) < mean_peak_dist / 2:
                mz = peaks_sim.pop(idx_closest)
                matched_peaks.loc[idx,"sim"] = mz
    unmached_peak_dict = np.array(peaks_sim)
    return matched_peaks, unmached_peak_dict


def calculate_mae(x_exp, y_exp, theoretical_isotopic_pattern_x , theoretical_isotopic_pattern_y, intensity_threshold=0.2):
    
    filtered_theoretical_indices = np.where(theoretical_isotopic_pattern_y > intensity_threshold)[0]
    filtered_theoretical_x = theoretical_isotopic_pattern_x[filtered_theoretical_indices]
    filtered_theoretical_y = theoretical_isotopic_pattern_y[filtered_theoretical_indices]

    # Interpolate the filtered theoretical isotopic pattern onto the x-values of the experimental data
    y_theo_interp = np.interp(x_exp, filtered_theoretical_x, filtered_theoretical_y)

    matched_peaks, unmached_peak_dict = match_peaks(y_theo_interp, y_exp,0.1)
    matched_peaks
    matched_peaks=matched_peaks.dropna()
    
    error = np.mean(np.abs(matched_peaks["trace"] - matched_peaks["sim"]))

    return error
    
###############################################################################
###############################################################################
###############################################################################

### INPUT PARAMETERS



import os

## Replace by your directory##
os.chdir('G:/Figure_1')

os.listdir()

file3="Zn7MT3_qex.csv"
x3,y3=prepare_data(file3)

h_range = (412, 416)  # Range for H
zn_range = (7, 7)  # Range for Zn
cu_range = (0, 0)  # Range for Zn

min_mz=1470
max_mz=1478

tick_spacing=2
npeaks=30
charge=5
inten_min=-0.05
inten_max=1.05
name_output="inset_1474"
intensity_threshold=0.2



file2="Zn7MT3_4Cu.csv"
x2,y2=prepare_data(file2)

h_range = (413, 416)  # Range for H
zn_range = (0, 4)  # Range for Zn
cu_range = (4, 6)  # Range for Zn

## for region with 10 metals
min_mz=1509
max_mz=1515
name_output="inset_1512"

## for region with 9 metals
min_mz=1496
max_mz=1502
name_output="inset_1499"

## for region with 8 metals
min_mz=1484
max_mz=1489
name_output="inset_1486"


## for region with 7 metals
min_mz=1470
max_mz=1477
name_output="inset_1475"

## for region with 6 metals
min_mz=1458
max_mz=1463
name_output="inset_1459"

## for region with 5 metals
min_mz=1445
max_mz=1450
name_output="inset_1448"

## for region with 4 metals
min_mz=1431.5
max_mz=1437
name_output="inset_1435"



tick_spacing=2
npeaks=30
charge=5
inten_min=-0.05
inten_max=1.05


intensity_threshold=0.2


################ RUNNING CODE ################

################# 1) GENERATE CANDIDATES


molecules = generate_molecules(h_range, zn_range,cu_range)

print("Generated molecules:")
for molecule in molecules:
    print(molecule)


################# 2) LOAD EXPERIMENTAL DATA


x_inset2,y_inset2=plot_inset(x3,y3,1,"TEST",min_mz,max_mz,peak1_min=min_mz,peak1_max=max_mz)

x_inset2,y_inset2=plot_inset(x2,y2,1,"TEST",min_mz,max_mz,peak1_min=min_mz,peak1_max=max_mz)


#################  3) FIND PEAKS EXPERIMENTAL DATA

from scipy.signal import find_peaks

# Assuming x_exp and y_exp are the arrays representing the experimental data
# Sort the experimental data based on x-values
sorted_indices = np.argsort(x_inset2)
x_exp_sorted = x_inset2[sorted_indices]
y_exp_sorted = y_inset2[sorted_indices]

# Find the peaks in the y_exp_sorted array
peaks, _ = find_peaks(y_exp_sorted,height=intensity_threshold)

# Retrieve the corresponding x and maximum y values
x_peaks = x_exp_sorted[peaks]
y_peaks = y_exp_sorted[peaks]

# Print the x and corresponding maximum y values
for x, y in zip(x_peaks, y_peaks):
    print(f"x = {x}, max y = {y}")


# Plot detected peaks

detected_peaks(x_inset2, y_inset2,x_peaks, y_peaks,"detected_peaks")




################# 3) ITERATE OVER MULTIPLE CANDIDATES 


best_mae=float('-inf')
best_metric = float('-inf')
best_molecule_r2 = None
best_molecule_chi = None
best_chi2= 10000
plot=False


results = []

for molecule in molecules:
    
    
    mz_theo, inten_theo = theo_iso(molecule, npeaks=npeaks, charge=charge)
    metric = calculate_r2(x_peaks,y_peaks, mz_theo, inten_theo,molecule,mz_min=min_mz, mz_max=max_mz,inten_min=inten_min, inten_max=inten_max,intensity_threshold=intensity_threshold,tick_spacing=tick_spacing,plot=plot)
    chi2=calculate_chi2(x_peaks,y_peaks, mz_theo, inten_theo,intensity_threshold=intensity_threshold)

    mass_error=calculate_mae(x_peaks,y_peaks, mz_theo, inten_theo, intensity_threshold=0.2)

    
    results.append({'Formula': molecule, 'R2': metric, 'Chi-square': chi2,"MAE (Da)": mass_error})


    
    if metric > best_metric:
        best_metric = metric
        best_molecule_r2 = molecule
        best_mae=mass_error
    
    if chi2 < best_chi2:

        best_molecule_chi = molecule
        best_chi2 = chi2
        best_mae=mass_error
        
results = sorted(results, key=lambda x: x['R2'], reverse=True)
csv_file = 'results.csv'
fieldnames = ['Formula', 'R2', 'Chi-square', "MAE (Da)"]

with open(csv_file, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print("Results saved to", csv_file)
print("Best molecule_chi:", best_molecule_chi)        
print("Best molecule_r2:", best_molecule_r2)
print("Best r2:", best_metric)
print("Best chi2:", best_chi2)
print("Best mae:", best_mae)
