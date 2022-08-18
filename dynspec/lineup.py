import matplotlib.pyplot as plt
import dedisperse_dynspec as dd

# Two MWA ones
yamllo = "1342623792.yaml"
yamlhi = "1342623496.yaml"
with open(yamllo, 'r') as yaml:
    params = dd.parse_yaml(yaml)
    dlo = dd.Dynspec(**params)
with open(yamlhi, 'r') as yaml:
    params = dd.parse_yaml(yaml)
    dhi = dd.Dynspec(**params)

dlo.set_freq_ref('low')
dhi.set_freq_ref('centre')

# Incl. Parkes
#dlo = dd.Dynspec(input='1342096104_dyn_dynamic_spectrum.csv', sample_time=0.5, freqlo=200.475, bw=0.64, t0=1342096104, transpose=True)
#dhi = dd.Dynspec(input='PKS_1342096266_dynamic_spectrum.csv', sample_time=0.1, freqlo=1216.0, bw=0.5, t0=1342096266.11, transpose=False)

DM = 288

dlo.add_dm_padding(DM, fill_value=0.0)
dhi.add_dm_padding(DM, fill_value=0.0)

dlo.dedisperse(DM)
dhi.dedisperse(DM)

dlo.t = dlo.get_time_at_infinite_frequency()
dhi.t = dhi.get_time_at_infinite_frequency()

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
fig.suptitle("DM = {:.2f}".format(DM))
dlo.plot(ax[1])
dhi.plot(ax[0])
plt.xlim([dhi.t[0], dlo.t[-1]])
plt.savefig('lineup.png')

