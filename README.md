# pylephare

Python wrapper of LePhare SED fitting.

http://www.cfht.hawaii.edu/~arnouts/LEPHARE/lephare.html

[See LePhare Documentation](http://www.cfht.hawaii.edu/~arnouts/LEPHARE/DOWNLOAD/lephare_doc.pdf)

It includes Monte Carlo generated SED spectra.


# Running the Code: Simple LePhare run()

The library enables you to run LePhare directly from python. 

## Before we start...
### Input Data formating

The code expects input data in form of a dataframe with the following column format:
`filtername1, filtername1.err, filtername2, filtername2.err, ... filternameN, filternameN.err, CONTEXT, Z-SPEC, STRING`

where `filtername{}` is a known filter by LePhare with the format `instrument.band`, for instance `sdss.u` or `ps1.z`.
The file `pylephare/data/data_test.csv` is an example. 

### Config File
LePhare request that you provide a config file. We invite you to check LePhare config documentation. An example is provided here: `pylephare/config/default.config` which is copied to your pylephare environment, see below.

All config file parameters could be changed later on.

### pylephare environment

the first time you import a pylephare modul, this will create the pylephare working environment:
- the directory `$LEPHAREWORK/pylephare`
- the directory `$LEPHAREWORK/pylephare/config`
- copy the package default config file  and parameter out config file inside `$LEPHAREWORK/pylephare/config`

## Let's go

Let's load the example dataframe (`pylephare/data/data_test.csv`) and configfiles (`$LEPHAREWORK/pylephare/config/default.config`):
```python
import pandas
from pylephare import io
data = pandas.read_csv(io._PACKAGE_ROOT+"/data/data_test.csv", sep=" ")
configfile = None # default file stored in your $LEPHAREWORK/pylephare/config/default.config
```
and let's load the `LePhare` class:
```python
from sedkcorr import lephare
lp = lephare.LePhare(data, configfile)
```
Then, to run the SED fit, simply do:
```python
fileouts = lp.run()
```

`fileouts` is a dictionary containing the fullpath of the config file used for the fit, the input data, the out results and the resulting spectra.


### Changing the configuration 
if you want to check what is inside the config file, simply do:
```python
lp.config.get_config("configfile")
```
use the `lp.config.set_value(KEY, VALUE)` method to affect the configfile.
you can swith on and off a given key by using the `lp.config.switch_{on/off}_key(KEY)` methods.


# Monte Carlo Fit 

A montecarlo fitting proceedure has been implemented to better account for mesurement errors when deriving the SED fitted parameters by lephare. 

For that, let's use only one entry of the catalog used in the first example:

```python
import pandas
data = pandas.read_csv(io._PACKAGE_ROOT+"/data/data_test.csv", sep=" ")
target = data.iloc[0]

configfile = None # use the default one
```
The `MCLePhare` class is a wrapper that draw `ndraw` realisation of the data given their respective errors and runs lephare (as illustrated in the first example) on each. Here let's do 50 draws:

```python
from sedkcorr import montecarlo
mclp = montecarlo.MCLePhare(target, 50)
```
and let's load lephare with our config file:
```python
mclp.load_lephare(configfile)
```

For this example, let's say we will largely increate the magnitude errorfloor of the data. We had 5 filters, will therefore need 5 values:
```python
mclp.lephare.config.set_intrinsic_error([0.1,0.08,0.08,0.08,0.1])
```
and then, let's run lephare on these 50 draws:
```python
mclp.run()
```

and here is how it looks like:
```python
mclp.show()
```


The results are stored as:
```python
mclp.catout
```
```
	Z_BEST	CHI_BEST	MOD_BEST	EXTLAW_BEST	EBV_BEST	DIST_MOD_BEST	NBAND_USED	Z_SEC	CHI_SEC	MOD_SEC	...	MASS_MED	MASS_SUP	SFR_BEST	SFR_INF	SFR_MED	SFR_SUP	SSFR_BEST	SSFR_INF	SSFR_MED	SSFR_SUP
IDENT																					
0	0.0319	0.460957E+01	11	0	0.000	0.357280E+02	5	-99.0000	0.100000E+10	-999	...	0.711027E+01	0.727246E+01	-0.209793E+01	-0.209343E+01	-0.179313E+01	-0.148264E+01	-0.868508E+01	-0.871817E+01	-0.818073E+01	-0.763570E+01
1	0.0319	0.352857E+01	1	1	0.160	0.357280E+02	5	-99.0000	0.100000E+10	-999	...	0.710554E+01	0.726496E+01	-0.103552E+01	-0.195251E+01	-0.163073E+01	-0.134995E+01	-0.701936E+01	-0.844735E+01	-0.784873E+01	-0.738573E+01
2	0.0319	0.374409E+01	1	1	0.180	0.357280E+02	5	-99.0000	0.100000E+10	-999	...	0.710554E+01	0.726523E+01	-0.102300E+01	-0.198167E+01	-0.166232E+01	-0.137990E+01	-0.701936E+01	-0.850053E+01	-0.792148E+01	-0.743389E+01
```
### Get synthetic photometry
to get the synthetic photometry for a given filter, nothing's easier do:

```python
mclp.get_synthetic_photometry('sdss.u', restframe=True, influx=False)
```

