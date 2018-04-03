import numpy as np
from numpy import random
from scipy import fftpack, signal
import matplotlib.pyplot as plt

# Data parameters
N = 1000000             # Length of data, in bins
T = 374                 # Pulsar period, in bins
W = 5                   # Pulse FWHM, in bins
PHI = random.randint(T) # Pulsar phase, in bins
DT = 0.001              # Time bin width in seconds
AMP = 0.20              # Pulsar amplitude, arbitrary units but the noise is set to 1.


# Generate data
# -------------

data = random.randn(N).astype(float)    # Noise with unit variance
time = np.arange(N) * DT

# Gaussian pulse profile.
sigma = W / 2.355
pulse_len = int(6*sigma)
pulse_len = pulse_len + pulse_len % 2 + 1
pulse_len = max(pulse_len, 3)
pulse = AMP * signal.gaussian(pulse_len, sigma)

# Add in pulse periodically.
for ii in xrange(PHI, N, T):
    this_data = data[ii:ii + pulse_len]
    this_n = len(this_data)
    this_data[:] += pulse[:this_n]

# Search it
# ---------

# First step: FFT the data.

# Zero pad. This helps with alignment of harmonics.
data = np.concatenate((data, np.zeros_like(data)))
print len(data)
plt.plot(data)
plt.show()
exit()

nfreq = N - 1
data_w = fftpack.fft(data)[0:nfreq]
data_w /= np.sqrt(len(data) / 2)    # Keeps the noise unit variance
freq = fftpack.fftfreq(len(data), d=DT)[0:nfreq]
periods = np.empty_like(freq)
periods[0] = 100000. * DT * N    # Infinity
periods[1:] = 1/freq[1:]

#plt.figure()
#plt.plot(freq, data_w.real)
#plt.plot(freq, data_w.imag)

#plt.figure()
#plt.plot(freq, abs(data_w))

if False:

    plt.figure()
    plt.plot(abs(data_w))

    w, = np.where(abs(data_w) > 500)
    print w
    plt.figure()
    plt.plot((freq[w] * T), '.')

    plt.figure()
    plt.plot(np.round(freq[w] * T / 1000.),
            (freq[w] * T) / np.round(freq[w] * T / 1000.), '.')


# Second step: for each frequency stack harmonics.

NHARMONICS = 100

SEARCH_T_MIN = 0.2
SEARCH_T_MAX = 5.
#SEARCH_T_MIN = (T - 0.2/50) * DT
#SEARCH_T_MAX = (T + 0.2/50) * DT

# Need to do 1 stack for ever frequency in NHARMONICS * frequency range. This
# is because the high harmonics are very sensitive to alignment, i.e. getting
# the frequency right.
f_min_ind = np.argmin(abs(freq - NHARMONICS/SEARCH_T_MAX))
# Careful that this isn't higher than Nyquist
f_max_ind = np.argmin(abs(freq - NHARMONICS/SEARCH_T_MIN))

#print f_min_ind, f_max_ind
#print 1/(freq[f_min_ind] / NHARMONICS), 1/(freq[f_max_ind] / NHARMONICS)

frac_df = 1. / NHARMONICS
# The frequency of each harmonic relative to the highest harmonic.
harmonic_factors = np.arange(frac_df, 1 + frac_df, frac_df, dtype=float)
for ii in range(f_min_ind, f_max_ind):

    # For this frequency, pull out the harmonics
    harmonic_inds = ii * harmonic_factors
    harmonic_inds = np.round(harmonic_inds).astype(int)
    harmonics = data_w[harmonic_inds]

    # Now each harmonic has an unknown phase depending on the phase of the
    # pulsar. However, this phase depends linearly on harmonic number. So can
    # pull that out with an FFT.
    harmonics_w = fftpack.fft(harmonics)
    #harmonics_w = fftpack.fft(np.concatenate([harmonics, np.zeros_like(harmonics)]))
    # Luckily this is only one parameter, so once we know how the phase
    # increases with harmonic, we also know the intercept. Perform rotation
    # such that all the power is in the real part.
    harmonics_w *= np.exp(1j * np.arange(1, NHARMONICS + 1))
    #harmonics_w *= np.exp(1j * np.arange(1, 2 * NHARMONICS + 1))
    # Normalize to unit variance in real part.
    harmonics_w /= np.sqrt(len(harmonics) / 2)
    #print np.var(harmonics_w.real)

    #plt.figure()
    #plt.plot(harmonics.real)
    #plt.plot(harmonics.imag)

    #plt.figure()
    #plt.plot(np.angle(harmonics))

    #plt.figure()
    #plt.plot(freq, abs(data_w))
    #plt.plot(freq[harmonic_inds], abs(data_w[harmonic_inds]), 'x')

    #print (np.sum(np.abs(harmonics)**2) - NHARMONICS) / np.sqrt(NHARMONICS)


    if np.any(harmonics_w.real > 8):
        trig_p = 1. / freq[harmonic_inds[0]]
        trig_sn = np.max(harmonics_w.real)
        trig_sn_ic = (np.sum(np.abs(harmonics)**2) - NHARMONICS)
        trig_sn_ic /= np.sqrt(NHARMONICS)
        print ("Trigger: T = %7.5f, S/N = %3.1f, (incoherent S/N = %3.1f)"
               % (trig_p, trig_sn, trig_sn_ic))
        plt.figure()
        plt.plot(harmonics_w.real)
        plt.plot(harmonics_w.imag)


plt.show()



