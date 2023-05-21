

# PolConvert in EVN data

The European VLBI Network (EVN) provides data correlated with the SFXC software correlator in FITS-IDI format.
When one or some antennas in the array recorded linear polarization, this tool allows to run PolConvert on such data, converting the polarization basis to circular. Note that this process is done internally before the EVN data is stored in the EVN Archive.


## Components

The implementation of PolConvert for EVN data lies in two files located in this directory:

- _polconvert_evn.py_: script to be called in order to run PolConvert.
- _polconvert_inputs.toml_: template file containing all input parameters required to run the conversion.



## Instructions

1. Once you have PolConvert installed in your system, copy the template input file _polconvert_inputs.toml_ into your working directory.
2. Modify the file to adequate the inputs to your particular EVN data.
   The file also contains the detailed explanation of each field.
3. Run the full program as `polconvert_evn.py  [your-modified-template-file]`.

All the output files and logs will be stored inside a `polconvert_logs` folder by default. The converted FITS-IDI files will also exhibit the extension `.PCONVERT`.

One should pay attention to the `Cross-Gains_ff-*.png` and `FRINGE.PLOTS/` plots to confirm that the conversion was successful.






