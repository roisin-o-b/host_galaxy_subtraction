A series of scripts forming a pipeline to subtract the host galaxy contribution in STELLA telescopic observations of the galaxy Mrk 1018, leaving only the nucleus.

# data_process.py
- converts from fitz files to fits, takes relevant info from fits header
- sections image to focus on galaxy in question (Mrk 1018 here) and 3 nearby stars
- subtracts and saves the local background
- does source extraction and finds Mrk 1018 and reference star 3 (used for calibration) from reference stars 1 & 2 (brighter)
- does photometry at the 4 positions
- writes all info to counts file

# combined_plot.py 
- averages the counts over the observations per night, finds the error, and writes to file
- converts counts to magnitudes and adds the zeropoint value, finds these uncertainties, and writes to file
- plots and saves initial light curve

# reproj.py 
- takes an image and scales, rotates and conserves flux to match another

# stella_pipeline_star.py 
- contains the functions necessary for host subtraction
- the inputs needed are: a count table, background images, a high-definition host-galaxy image, a count value for a reference star in the high-def. image, and the high-def. image observation date
- the first for loop takes each STELLA observation and convolves with the reprojected host-galaxy image, adds the STELLA background, adds noise, subtracts the background from the created image, finds the host-galaxy centroid and does photometry, saves the host-galaxy contribution for each STELLA observation, and matches the count value for the reference star using the high-defintion image observation date
- after the first for loop, the script calculates the ratio between the reference star counts in the STELLA image and those in the high-def. image, scales the host-galaxy counts according to the reference star ratio, and subtracts the resultant host-galaxy counts from the combined count value
- the second for loop averages over the observations per night and returns the count values
    
# multiruns.py 
- imports function described above and runs it 1000 times, saving obsservation counts for each run

# lc_fns.py 
- functions needed to create the light curve 
- _group_dps_ saves a table of dates for later, extracts all outputs from each run and arranges by date, sorts from highest to lowest and saves each range per by date
- _agnonlystats_ takes the range of data points for each date, and extracts the expectation value and 1, 2, 3 sigma values, then saves these values in a table, arranged by date
- _plot_agnstats_ converts the counts to magnitudes, separates the errors from the limits, plots this, and saves data points
    
# call_fns_on_data.py
- calls the above functions on the data from the multiple runs
  
# plot_lcs.py 
- takes the three outputs for three high-def. images (these approximate the expectation value and 1 sigma errors in the choice of host galaxy image)
- uses the central value and, if non-zero, the 1-sigma error values
- if the 1-sigma error values are zero, the script adds the lower limit
- plots the data points, creates the shaded region behind the data points, adds the grey block for the sunblock period, and saves the figure
    
