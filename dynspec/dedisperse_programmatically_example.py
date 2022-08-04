import matplotlib.pyplot as plt
import dedisperse_dynspec as dd

# Load a "raw" dynamic spectrum from a CSV file
dynspec = dd.Dynspec('1342623792_dyn_dynamic_spectrum.csv', sample_time=0.5, freqlo=72.955, bw=1.28, time_offset=1342623792, transpose=True)

# Dedisperse it, using an arbitrary reference frequency
dynspec.dedisperse(285, freq_ref=72.955)
# ...the dedispersed spectrum is now in dynspec.dynspec (NumPy array)

# Plot it
fig, ax = plt.subplots(1,1)
dynspec.plot(ax)
plt.show()
