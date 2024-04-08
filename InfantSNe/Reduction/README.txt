#### Reduction folder ####

plot_all.ipynb
- Plot SNID comparisons individually and on the same plot 
- Plot all the SNe spectra along with their types and some chemical elements
- Plot all the SNe spectra with respective IAU names
- Compute redshift from host galaxy emission lines present in the SNe spectra
- Check clean spectra after removal of host lines
- Merging R1000B and R1000R 
- Comparison plots --> might be useful for section 5 to compare all type Ia with 2011fe and 2020jgl with 2023bee, 2021aefx, 2019ein, etc. 

pypeit_tutorial.ipynb
- Tutorial on how to reduce GTC/OSIRIS data with PypeIt 
- For the fluxing it is needed to smooth the sensitivity function using the script smooth_sensfunc.ipynb

tell_corr.ipynb
- Script created to perform the telluric correction

"data" folder
- all raw and reduced data from different programmes 

"final_txt"
- fits_to_txt: all reduced spectra (telluric-corrected and not telluric-corrected) for this project and other projects, in txt files
- merging: all final spectra (reduced, telluric-corrected and merged) for this project. Sanity check plot for the merging process and the computation of redshifts 