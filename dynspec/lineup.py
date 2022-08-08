import matplotlib.pyplot as plt
import dedisperse_dynspec as dd

# Incl. Parkes
#dlo = dd.Dynspec('1342623792_dyn_dynamic_spectrum.csv', sample_time=0.5, freqlo=72.955, bw=0.64, time_offset=1342623792, transpose=True)
#dhi = dd.Dynspec('1342623496_dyn_dynamic_spectrum.csv', sample_time=0.5, freqlo=200.955, bw=1.28, time_offset=1342623496, transpose=True)

# Two MWA ones
dlo = dd.Dynspec('1342096104_dyn_dynamic_spectrum.csv', sample_time=0.5, freqlo=200.475, bw=0.64, time_offset=1342096104, transpose=True)
dhi = dd.Dynspec('PKS_1342096266_dynamic_spectrum.csv', sample_time=0.1, freqlo=1216.0, bw=0.5, time_offset=1342096266.11, transpose=False)

DM = 288
dlo.dedisperse(DM, freq_ref='low')
dhi.dedisperse(DM, freq_ref='centre')

dlo.t = dlo.get_time_at_infinite_frequency()
dhi.t = dhi.get_time_at_infinite_frequency()

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
fig.suptitle("DM = {:.2f}".format(DM))
dlo.plot(ax[1])
dhi.plot(ax[0])
plt.xlim([dhi.t[0], dlo.t[-1]])
plt.savefig('lineup.png')

