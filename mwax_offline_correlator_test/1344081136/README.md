# MWAX Offline Correlator Test - VCS observation

[MWA Observation Status page](http://ws.mwatelescope.org/observation/obs/?obs_id=1344081136)

#### Wed, 10 Aug 2022

- Started downloading from ASVO (JobID: 583749)

The surrounding correlator jobs selected an integration time of 0.5 s and a frequency resolution of 10 kHz.
Understandably, these settings have to be folded into the setting up of the offline correlator as well.
The instructions for doing so are found on the [MWAX Offline Correlator wiki page](https://wiki.mwatelescope.org/pages/viewpage.action?spaceKey=MP&title=MWAX+Offline+Correlator).

The first step is modifying the `.sub` file headers:

| Key | Existing Value | Set to | Notes |
| :-- | :------------- | :----- | :---- |
| MODE | MWAX_VCS | MWAX_CORRELATOR | |
| INT_TIME_MSEC |  | 500 | 500ms = 0.5s |
| FSCRUNCH_FACTOR |  | 50 | 50 x 200 Hz = 10 kHz |
| FINE_CHAN_WIDTH_HZ |  | 10000 | 10000 Hz = 10 kHz |
| NFINE_CHAN |  | 128 | 128 x 10 kHz per coarse channel |

To effect these changes on Pawsey (**don't try yet**):

```
module use /pawsey/mwa/software/python3/modulefiles
module load mwax_offline_correlator

mwax_update_subfile_header -s MODE=MWAX_CORRELATOR -s INT_TIME_MSEC=500 -s FSCRUNCH_FACTOR=50 -s FINE_CHAN_WIDTH_HZ=10000 -s NFINE_CHAN=128 *.sub
```
