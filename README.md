## **Nanothermometry analysis pipeline**
Welcome! In this repository you will find the code we developed to process the data presented in this [work](http://google.com)

The analysis pipeline was designed for an implementation of the Anti-Stokes nanothermometry. The scripts presented here process hyperspectral images acquired with a custom-built confocal microscope. Each hyperspectral image is composed of `NxN` photoluminescence (PL) spectra acquired at a specific spatial position of the sample, using the same wavelength range.

### The files
Here you will find 4 kind of files:
1. An **auxiliary** file which contains all the *small* functions needed to run the analysis. Filename: functions_for_photoluminescence.py
2. Stand-alone scripts that account for each **step** of the analysis pipeline. Filename: stepX_KEYWORD_OF_THE_STEP.py
3. The **main** script. You should execute this script to analyze the data. Filename: find_beta_or_Tzero.py
4. Definitions to **plot** some fancy graphics (matplotlib). Filename: for_confocal.mplstyle

### Organization
The analysis pipeline code is organized as follows:
+ The **main** script calls all steps. The order can't be altered. Each step needs information that the previous one has created. Step 0, must be executed before Step 1, 1 before 2, and so on...
+ Once you have run Step 0, several new files are created (processed data). So next time you run the main script won't be neccesary to run Step 0 again. For instance, you can continue only with Step 1 (or 1 and 2, or 1, 2 and 3, etc.).
+ The keywords in the filename of the step scripts are self-explanatory. Inside each script some comments may you find. Variables name are also self-explanatory.
+ Step 4 is written to gather all saved data and to do some statistics. This step is not critical.

### Input parameters
The **main** script has several variables that determine the behavior of the code. These variables are passed to each step script. These variable are arguments for specific functions. Along with its value, inside a subfunction called `run_analysis` in the **main** script, you will find that each variable has a brief comment that explains its meaning. Examples of these variables are: `totalbins`, `image_size_px`, `camera_px_length`, etc...

### Outputs
The analysis pipeline creates a new folder called "processed_data". This folder contains several subfolders which in turn contains `.dat`, `.txt` (both ASCII encoded) and plots:
1. One folder for each hyperspectral image you have acquired. Data inside:
  * Calibration data (calculated irradiance, counts per second, raw spectra, plots)
  * Raw PL spectra (optional, set input as `plot_flag = True` in **main** script)
  * Binned PL spectra. NUmber of bins can be set through `totalbins` variable in the **main** script. This subfolder contains:
    - All binned PL spectra in a single file (all_bins.dat)
    - Wavelength array used (londa.dat)
    - Irradiance for each bin (bin_irradiance.dat)
    - Distribution of spectra considering its irradiance (bin_distro.dat)
    - Relevant plots
  * Normalized Stokes PL spectra (to check Stokes proportionality with incident irradiance)
  * Plots of the Anti-Stokes PL spectra
  * Ratios between different Anti-Stokes PL spectra (antistokes_quotient folder)
  * Sum of all spectra (spr folder). This was used to extract information related to the Surface Plasmon Resonance)
  * A folder named *matrix* were most valuable information is saved:
    - A matrix of photothermal coefficients, a.k.a., beta (if beta was fitted)
    - A matrix of Tzero (if Tzero was fitted)
    - A matrix of goodess of the fit
    - A matrix for each variable error
    - A subfolder called temp_vs_irrad where NP's Temperature is calculated as a function of the incident Irradiance
2. A folder for all common plots (for fast comparison)
3. A folder that contains statistics about all hyperspectral images acquired (only created if multiple NPs are being analyized, i.e. `single_NP_flag = False`, and step4 is executed)

#### Disclaimer
Some plots were adjusted to fit our requirements (axis limits, scale, colors, dpi). Modify them if you need. Defaults can be edited at for_confocal.mplstyle and other fields should be edited inside each step script.

#### Contact
If you find a bug, have a question or just want to make a comment, you can open a new issue or reach me anywhere on the web.
