A series of scripts forming a pipeline to subtract the host galaxy contribution in STELLA telescopic observations of the galaxy Mrk 1018, leaving only the nucleus.

# data_process.py
- converts from fitz files to fits, takes relevant info from fits header
- sections image to focus on galaxy in question (Mrk 1018 here) and 3 nearby stars
- subtracts and saves the local background
- does source extraction and finds Mrk 1018 and reference star 3 (used for calibration) from reference stars 1 & 2 (brighter)
- does photometry at the 4 positions
- writes all info to counts file

# combined_plot.py 
- takes count table, averages over obs per night, finds photon error, writes to obs file, converst counts to mags and adds zeropoint, saves mrk 1018 mags and mags error, plots and saves combined LC

*****************2017_07_30_VIMOS_mean/reproject***************

reproj.py - takes high def image,scales, rotates and conserves flux to match stella images

*****************2017_07_30_VIMOS_mean***************

stella_pipeline_star.py - contains function that does host subtraction
    inputs: count table, background images, high def host galaxy image, high def ref star 3 count, high def obs date
    first for loop: for each stella obs; takes stella seeing, convolves with host galaxy image, adds stella background, adds noise, subtracts bakckground from created image, find host centroid and do photometry, save host contribution for each obs, find ref star 3 counts from high def obs date
    after first for loop: find ratio between stella ref 3 star counts and high def star count, scale host according to ref star ratio, subtract host counts from combined counts
    second for loop: average over obs per night, return obs counts
    
multiruns.py - import function described above and runs it 1000 times, saves obs counts for each run

lc_fns.py - functions needed to create the light curve 
    group_dps: saves a table of dates for later, extracts all outputs from each run and arranges by date, sorts from highest to lowest and saves each range per by date
    agnonlystats: takes range of dps for each date and extracts expectation value and 1, 2, 3 sigma values, saves in table arranged by date
    plot_agnstats: convert to mags, separate errors from limits, plots and saves data points
    
call_fns_on_data.py
    
***************************2016_09_24_VIMOS_onesigmabright***********************

same procedure as in 2017_07_30_VIMOS_mean to find upper shaded region systematic error


****************************2017_09_14_GMOS_onesigmafaint**************************

same procedure as in 2017_07_30_VIMOS_mean to find upper shaded region systematic error

***************************final_plot****************************

plot_lcs.py - takes the three outputs for the three high def images, uses centre and 1 sigma error if non-zero, if zero adds lower limit, creates shaded region behind data points, adds grey block for sunblock period, saves figure

    
