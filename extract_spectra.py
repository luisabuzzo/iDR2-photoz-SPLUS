#############################################################################
#
# INPUTS
# BC03 model *ised_ASCII files - both low resolution with 1221 wavelengths or high resolution with 6900 wavelengths
#
# OUTPUTS
# SEDs array - with dimensions (1222, 222) or (6901, 222) depending on lr/hr - the first row contains the ages which correspond to the column below it (the first entry in this row is 0 - the first column (seds[:,0]) contains the wavelength values and all the following columns contain the flux values for the corresponding wavelength at each age value in the ages_array. The array therefore looks something like this:
#                       0     a1      a2      a3      ...
#                       W1    f1,1    f2,1    f3,1    ...
#                       W2    f1,2    f2,2    f3,2   ...
#                       W3    f1,3    f2,3    f3,3   ...
#                       ...   ...     ...     ...     ...
#
#############################################################################

import numpy as np
import os.path
from astropy.table import Table,vstack,hstack
import pandas
import matplotlib.pyplot as plt

# used to find the index of the closest age available to a desired age (in years)
def find_closest_value(a,a_bin):

    da = abs(a_bin-a)
    idx = np.where(da == min(da))
    idx = idx[0]
    return idx

# read the extracted BC03 spectra file and return the spectra as arrays
def import_spectra(filename):
    
    df = pandas.read_table(filename,header=None,dtype=np.float32,delimiter=' ')
    df = np.array(df)
    specs_ages = df[0,1:]
    specs_fluxes = df[1:,1:]
    specs_wave = df[1:,0]
    nwaves = len(specs_wave)
    nages = len(specs_ages)
    return specs_ages, specs_fluxes, specs_wave, nwaves, nages

def create_arrays(data_array, num_ages, num_wave):

    #if data_array[0] > num_ages-1:
    new_array = data_array[1:int(data_array[0])+1]
    n_data =[]
    n_data = np.append(n_data, data_array[len(new_array)+1:])
    new_data = n_data.reshape(len(n_data),1)
    return new_array, new_data

def process_ised_ascii_file(path_bc03, sedfile, path_output, outputfile):

#Open the file loaded into the function and read each line
    file = open(path_bc03+sedfile, 'r')
    lines = file.readlines()
    file.close()
    print("File read... extracting the spectra... this takes a while...")

#Determine whether the file is high resolution or low resolution and therefore define the number of wavelengths which will be present
    if 'lr' in sedfile:
        num_wave = 1221
    else:
        num_wave = 6900

#The number of ages is constant across all the model files
    num_ages = 221

#Remove the 5 lines in the middle of each file which contain text defining the model parameters - these prevent the code from floating the data and making the correct arrays
    if 'chab' in sedfile:
        lines_step = lines[0:37]
        new_lines = np.append(lines_step, lines[42:])
    if 'salp' in sedfile:
        lines_step = lines[0:37]
        new_lines = np.append(lines_step, lines[43:])

#Split the lines so that the data contains one large array of all the data in a row
    initial_data = []
    nlines = len(new_lines)
    for n in range(0,nlines):
        numbers = new_lines[n].split()
        if np.int((n+1)/10000)*1.0 == ((n+1)/10000):
            print('processed', n+1,'of', nlines, 'lines of input...')
        initial_data = np.append(initial_data,numbers)

#Float the data and reshape the array into a single column of data.
    float_data = initial_data.astype(float)
    data = float_data.reshape(len(float_data),1)

#Extract the ages and wavelenth arrays using the above function
    age, working_data = create_arrays(data,num_ages,num_wave)
    wavelength, work_data = create_arrays(working_data,num_ages,num_wave)

# work_data now contains 221 lists of 1221 fluxes, each preceded by the number '1221', and followed by the number '52' and 52 numbers (dont know what they mean, not needed)
    final_output_array = np.zeros((num_wave+1, num_ages+1))
    final_output_array[1:,0] = wavelength.flatten() #.reshape(len(wavelength),1)
    final_output_array[0,1:num_ages+1] = age.flatten() 
    age_processed = 0
    for i in range(0,num_ages):
        get_vals = work_data[1:num_wave+1]
        print("Reading age",i+1,age[i],'(yr) of',num_ages,'ages with',len(get_vals),'of',work_data[0],' fluxes')
        final_output_array[1:,i+1] = get_vals.flatten() 
        # remove the number '1221', the 1221 fluxes that were just stored, the number '52' and the 52 values from the array
        if age_processed < num_ages:
            work_data = work_data[num_wave+1+1+52:-1]
        age_processed += 1

    np.savetxt(path_output + outputfile, final_output_array, fmt='%1.8E')
    print(path_output + outputfile, ' written.')

    return

def main():

    path_bc03 = './'        # location of the BC03 .ised_ASCII files
    path_output = './'      # location of the output

    # define the bc03 SED file
    res = 'lr'              # resolution
    Z = '62'                # metallicity
    sfh = 'tau5'            # star formation history
    imf = 'chab'            # imf
    dust = '_dust00'        # dust
    sedfile = 'bc2003_' + res + '_m' + Z + '_' + imf + '_' + sfh + dust + '.ised_ASCII'
    outputfile = 'extracted_' + sedfile

    # extract the spectra from the .ised_ASCII. It will skip this if it finds a previously made 
    # file extracted_<sedfile> in the foler specified by <path_output>
    if os.path.isfile(path_output+outputfile):
        print("Extracted file ", path_output+outputfile, " found. Good!")
    else:
        print("Extracted file ", path_output+outputfile," does not yet exist. Building it...")
        process_ised_ascii_file(path_bc03,sedfile,path_output,outputfile)

    # read the extracted file. It will return:
    # nages : the number of ages for which there is a spectrum   
    # specs_age : the list of ages of the spectra 
    # nwavelengths : the number of wavelengths in each spectrum   
    # specs_wavelengths : the list of wavelengths for each spectrum
    # specs_flux : an array where each column i corresponds to the fluxes of the ith age    
    specs_age, specs_flux, specs_wavelength, nwavelengths, nages = import_spectra(path_output+outputfile)
    print("Extracted file ", path_output+outputfile, " read.")
    print("Numer of age bins: ", len(specs_age))
    print("Numer of wavelength bins: ", len(specs_wavelength))

    # select a spectrum of a certain age in years
    age_wanted = 1.e8 # for 100 Myr
    age_bin = find_closest_value(specs_age,age_wanted)
    age_myr = specs_age[age_bin][0]/1.e6

    # this is your spectrum
    lam = specs_wavelength
    flux = (specs_flux[:,age_bin]).flatten()

    # plot the spectrum
    plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Wavelength ($\\AA$)')
    plt.ylabel('$F_\\lambda$ (units)')
    title = sedfile 
    age_info = "Age: " + np.str_(age_myr) + ' Myr'
    plt.plot(lam,flux)
    plt.text(1e4,0.1*np.max(flux),age_info,fontsize=20)
    plt.title(title)
    plt.show()

if __name__ == '__main__':
    main()

    
