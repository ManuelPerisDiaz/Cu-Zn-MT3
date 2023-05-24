# Cu-Zn-MT3
To accurately determine the stoichiometry of metal-protein complexes to native mass spectrometry data, we prepared a python script that generates multiple theoretical protein/peptide isotopic distributions, and determines the molecule with the best R2 and chi-square.
The fitting_script.py
1.	GENERATE CANDIDATES
The generate_molecules function is called with the specified h_range, zn_range, and cu_range parameters to generate a list of molecule candidates. 
2.	LOAD EXPERIMENTAL DATA
The plot_inset function is called to plot and retrieve the x_inset2 and y_inset2 arrays representing the experimental data.
3.	FIND PEAKS IN EXPERIMENTAL DATA
The experimental data is sorted based on the x_inset2 values. Then, the find_peaks function from scipy.signal is used to find the peaks in the sorted y_inset2 array. The peak positions (x_peaks) and corresponding maximum peak heights (y_peaks) are extracted. The detected_peaks function is called to plot the detected peaks and visualize.
4.	ITERATE OVER MULTIPLE CANDIDATES
A loop is started to iterate over each molecule in the molecules list. Inside the loop the theo_iso function is called to generate the theoretical isotopic pattern (mz_theo and inten_theo) for the current molecule candidate. Then, the calculate_r2 and calculate_chi2 functions are called to calculate the R2 and chi-square metric for the current molecule candidate, using the experimental peaks (x_peaks and y_peaks) and the theoretical isotopic pattern. Also, the calculate_mae function is called to calculate the mean absolute mass error for the current molecule candidate, using the experimental peaks and the theoretical isotopic pattern. Here, the molecule candidate that provide the highest R2, lowest chi-square is printed separately.
5.	SORT AND SAVE RESULTS TO CSV
The results list is sorted based on the R2 metric in descending order and save into a CSV file named 'results.csv' with fieldnames defined as ['Formula', 'R2', 'Chi-square', "MAE (Da)"].
6.	PRINT FINAL RESULTS
The final results, including the best molecule based on chi-square, R2, and its values, along with the mean absolute mass error, are printed to the console. The program also plots and saves all of the fittings.
