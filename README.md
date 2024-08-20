# NMFig

**An inverse-gamma source variance prior with factorized parameterization for audio source separation.**

MATLAB code for the Variational EM algorithm for Source Separation using unfactorised NMF via Inverse-Gamma priors (NMFIG), which is presented in the paper:

```
D. Kounades-Bastian, L. Girin, X. Alameda-Pineda, S. Gannot, R. Horaud, 
"An Inverse-Gamma Source Variane Prior with Factorized Parametrization 
for Audio Source Separation", ICASSP, 2016
```

### DEMO


Simply run `example.m` by typing in the matlab console

```python
>> example
```

`example.m` will first generate a stereo mix with 3 sources (by loading trueSrc1.wav, ..) 

Then will call `initNMF()` to provide initial values for the VEM parameters (through binary masking)

Lastly it will apply `nmfig()` to separate the sources.

Separated sources will be saved as .wav files (estimatedSrc1.wav, etc.) in the directory `./results/`



nmfig.m implements the NMFIG algorithm.

For documentation on NMFIG type in the matlab console:

```
>> help nmfig
```

The functions in the folder ./aux_tools are downloaded from [multi_nmf](http://www.irisa.fr/metiss/ozerov/Software/multi_nmf_toolbox.zip). 

### SLIDES

From [inria.fr](https://team.inria.fr/perception/research/nmfig/)








