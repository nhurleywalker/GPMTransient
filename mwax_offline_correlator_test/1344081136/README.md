# MWAX Offline Correlator Test - VCS observation

[MWA Observation Status page](http://ws.mwatelescope.org/observation/obs/?obs_id=1344081136)

### Wed, 10 Aug 2022

- Started downloading from ASVO (JobID: 583749)

The surrounding correlator jobs selected an integration time of 0.5 s and a frequency resolution of 10 kHz.
Understandably, these settings have to be folded into the setting up of the offline correlator as well.
The instructions for doing so are found on the [MWAX Offline Correlator wiki page](https://wiki.mwatelescope.org/pages/viewpage.action?spaceKey=MP&title=MWAX+Offline+Correlator).

#### Modifying subfile headers

The first step is modifying the `.sub` file headers:

| Key | Existing Value | Set to | Notes |
| :-- | :------------- | :----- | :---- |
| MODE               | MWAX_VCS   | MWAX_CORRELATOR | |
| INT_TIME_MSEC      | 500        | 500             | 500ms = 0.5s                    |
| FSCRUNCH_FACTOR    | 50         | 50              | 50 x 200 Hz = 10 kHz            |
| FINE_CHAN_WIDTH_HZ | 10000      | 10000           | 10000 Hz = 10 kHz               |
| NFINE_CHAN         | 128        | 128             | 128 x 10 kHz per coarse channel |
| EXPOSURE_SECS      | 64         | 64              | |
| OBS_ID             | 1344081136 | 1344081136      | |
| OBS_OFFSET         | 0 (etc.)   | 0 (etc.)        | |

In this case, the only change that needs to happen is to the `MODE`.
This can be effected using the [`mwax_update_subfile_header`](https://github.com/MWATelescope/mwax_user_tools) utility.

To effect these changes on Pawsey:

```
module use /pawsey/mwa/software/python3/modulefiles
module load mwax_offline_correlator/136T

mwax_update_subfile_header -s MODE=MWAX_CORRELATOR *.sub
```

#### Running the correlator

See [correlate.sbatch](correlate.sbatch) (also still a **work in progress**).

##### Notes (Log)

- `NINPUTS_XGPU` = 288, but `mwax_db2correlate2db` complains about it not equalling the xGPU configuration:
```
...
[2022-08-13-18:07:21] ERR: mwax_db2correlate2db_open: NINPUTS_XGPU [288] does not match the current xGPU configuration [272]
[2022-08-13-18:07:21] ERR: Error calling open function
[2022-08-13-18:07:21] ERR: mwax_db2correlate2db main: error during PSRDADA client read call
...
```
