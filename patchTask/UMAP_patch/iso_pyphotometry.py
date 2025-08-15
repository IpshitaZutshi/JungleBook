
# For analyzing pyphotometry data. Similar to pyphotometry_preprocessing.py, but plots fewer things 
# and uses the pyphotometry function preprocess_data. 


import os
import numpy as  np
import matplotlib.pyplot as plt
from scipy.signal import medfilt, butter, filtfilt
from scipy.stats import linregress
from scipy.optimize import curve_fit, minimize
import scipy


from data_import import import_ppd, preprocess_data

#data_folder = r'\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\test_random_pulses'
data_folder = r'\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\patchTask\N14\N14_250519_170548' 
file = 'N14_striatum-2025-05-19-170555' 
data_filename = file + '.ppd'
data = import_ppd(os.path.join(data_folder, data_filename))   #, low_pass=20, high_pass=0.001)
grabDA_raw = data['analog_1'] # green 
cherry_raw = data['analog_2'] # red
time_seconds = data['time']/1000
sampling_rate = data['sampling_rate']
time = data['time']

# for stim
pulse_times_1 = data['pulse_times_1']
stim_time = data['digital_1']
barcode = data['digital_2'] 




preprocessed_data = {}


processed_signal = preprocess_data(data_dict=data, 
                                   signal="analog_1", 
                                   control="analog_2", 
                                   low_pass=10,
                                   normalisation="dF/F",
                                   plot=True,
                                   fig_path = os.path.join(os.path.expanduser('~'), 'Downloads', 'N7_striatum.png'))  # change name you want to save as

# Save plot to the Downloads folder
#downloads_path = os.path.join(os.path.expanduser('~'), 'Downloads', 'N7_striatum_plot.png')
#plt.gcf().savefig(downloads_path)

# Plot raw signals
fig,ax1=plt.subplots()  # create a plot to allow for dual y-axes plotting
plot1=ax1.plot(time_seconds, grabDA_raw, 'g', label='grabDA') #plot grabDA on left y-axis
ax2=plt.twinx()# create a right y-axis, sharing x-axis on the same plot
plot2=ax2.plot(time_seconds, cherry_raw, 'r', label='mCherry') # plot mCherry on right y-axis

# Add labels for both axes
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('grabDA_raw', color='g')
ax2.set_ylabel('mCherry', color='r')

plt.show()


# Plot rewards times as ticks.
#reward_ticks = ax1.plot(reward_cue_times, np.full(np.size(reward_cue_times), 1.625), label='Reward Cue', color='w', marker="|", mec='k')


""" GET BARCODE TIMESTAMPS """

# barcodes on and off are giving indices of timestamp where the barcode goes on or off

timestampsOn = []
timestampsOff = []
for i,signal in enumerate(barcode):
    if signal == 1 and barcode[i-1]==0:
        timestamp_On = i
        timestampsOn.append(i)
    if signal == 0 and barcode[i-1]==1:
        timestamp_Off = i
        timestampsOff.append(i)

timestampsOnOff = []
for On,Off in zip(timestampsOn,timestampsOff):
    timestampsOnOff.append([On, Off])

timestampsOn = np.array(timestampsOn)
timestampsOff = np.array(timestampsOff)

preprocessed_data['sampling_rate'] = sampling_rate
preprocessed_data['barcodesOn'] = timestampsOn # formerly timestampsOn
preprocessed_data['barcodesOnOff'] = timestampsOnOff # formerly timestampsOnOff
preprocessed_data['highLow'] = barcode
preprocessed_data['timestamps'] = time_seconds
preprocessed_data['time_ms'] = time


#%% normalisation to save to matlab

"""DENOISING"""
# Lowpass filter - zero phase filtering (with filtfilt) is used to avoid distorting the signal.
b,a = butter(2, 10, btype='low', fs=sampling_rate)
grabDA_denoised = filtfilt(b,a, grabDA_raw)
cherry_denoised = filtfilt(b,a, cherry_raw)


"""PHOTOBLEACHING CORRECTION"""

#METHOD 1
# The double exponential curve we are going to fit.
def double_exponential(t, const, amp_fast, amp_slow, tau_slow, tau_multiplier):
    '''Compute a double exponential function with constant offset.
    Parameters:
    t       : Time vector in seconds.
    const   : Amplitude of the constant offset. 
    amp_fast: Amplitude of the fast component.  
    amp_slow: Amplitude of the slow component.  
    tau_slow: Time constant of slow component in seconds.
    tau_multiplier: Time constant of fast component relative to slow. 
    '''
    tau_fast = tau_slow*tau_multiplier
    return const+amp_slow*np.exp(-t/tau_slow)+amp_fast*np.exp(-t/tau_fast)

# Fit curve to grabDA signal.
max_sig = np.max(grabDA_denoised)
inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
bounds = ([0      , 0      , 0      , 600  , 0],
          [max_sig, max_sig, max_sig, 36000, 1])
grabDA_parms, parm_cov = curve_fit(double_exponential, time_seconds, grabDA_denoised, 
                                  p0=inital_params, bounds=bounds, maxfev=1000)
grabDA_expfit = double_exponential(time_seconds, *grabDA_parms)

# Fit curve to mCherry signal.
max_sig = np.max(cherry_denoised)
inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
bounds = ([0      , 0      , 0      , 600  , 0],
          [max_sig, max_sig, max_sig, 36000, 1])
cherry_parms, parm_cov = curve_fit(double_exponential, time_seconds, cherry_denoised, 
                                  p0=inital_params, bounds=bounds, maxfev=1000)
cherry_expfit = double_exponential(time_seconds, *cherry_parms)

grabDA_detrended = grabDA_denoised - grabDA_expfit
cherry_detrended = cherry_denoised - cherry_expfit


"""MOTION CORRECTION"""
slope, intercept, r_value, p_value, std_err = linregress(x=cherry_detrended, y=grabDA_detrended)

#print('Slope    : {:.3f}'.format(slope))
#print('R-squared: {:.3f}'.format(r_value**2))


grabDA_est_motion = intercept + slope * cherry_detrended
grabDA_corrected = grabDA_detrended - grabDA_est_motion


"""NORMALIZATION"""
 #METHOD 1: dF/F
grabDA_dF_F = 100*grabDA_corrected/grabDA_expfit

#METHOD 2: Z-SCORE
grabDA_zscored = (grabDA_corrected-np.mean(grabDA_corrected))/np.std(grabDA_corrected)

#%% SAVE DATA IN MATLAB FORMAT

preprocessed_data['grabDA_z'] = grabDA_zscored 
preprocessed_data['grabDA_df'] = grabDA_dF_F
preprocessed_data['grabDA_raw'] = grabDA_raw
preprocessed_data['pulse_times'] = pulse_times_1
preprocessed_data['stim_time'] = stim_time
scipy.io.savemat(os.path.join(data_folder, f'{file}_photometry.mat'), {'photometryData': preprocessed_data})




# photometry around stim


#%% Plot signal at pulse times

time_step = 1000 / sampling_rate  # Time step in ms

window_ms = 5000  # 5 seconds before and after (total window = 10 seconds)
window_samples = int(window_ms * sampling_rate / 1000)  # Convert ms to samples

avg_time = np.linspace(-window_ms, window_ms, 2 * window_samples)

averaged_signal = np.zeros(2 * window_samples)

# Loop through pulse times and extract segments
valid_segments = 0
for pulse_time in pulse_times_1:
    # Find the index of the current pulse time
    pulse_index = np.searchsorted(time, pulse_time)
    
    # Define start and end indices for the window
    start_index = pulse_index - window_samples
    end_index = pulse_index + window_samples
    
    # Ensure indices are within valid range
    if start_index >= 0 and end_index <= len(grabDA_raw):
        segment = grabDA_zscored[start_index:end_index] # grabDA_detrended, raw, zscore, or df grabDA_zscored
        averaged_signal += segment
        valid_segments += 1

# Compute the average of the segments
if valid_segments > 0:
    averaged_signal /= valid_segments

# Plot the averaged signal
plt.figure(figsize=(10, 6))
plt.plot(avg_time, averaged_signal, label='Averaged Signal')
plt.axvline(0, color='red', linestyle='--', label='Pulse Time')
plt.xlabel('time relative to pulse (ms)')
plt.ylabel('signal (z-score)')   ##### CHANGE TO CORRESPOND TO PROPER SIGNAL METRIC
plt.title('Averaged grabDA Signal 5s Before and After Stim')
plt.legend()
plt.grid(True)
plt.show()


