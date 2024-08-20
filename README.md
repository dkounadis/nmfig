# NMFig

**An inverse-gamma source variance prior with factorized parameterization for audio source separation.**

MATLAB code for the Variational EM algorithm for Source Separation using unfactorised NMF via Inverse-Gamma priors (NMFIG), which is presented in the paper:

```
D. Kounades-Bastian, L. Girin, X. Alameda-Pineda, S. Gannot, R. Horaud, "An Inverse-Gamma Source Variane Prior with Factorized Parametrization for Audio Source Separation", ICASSP, 2016
```

### DEMO

```python
Simply run example.m by typing in the matlab console

   >> example

example.m will first generate a stereo mix with 3 sources (by loading trueSrc1.wav,..) 

Then call initNMF() to provide initial values for the VEM parameters (through binary masking)

Finally it will call nmfig() to separate that mix.

Separated sources will be written in .wav files (estimatedSrc1.wav, etc.) in the directory ./results/



nmfig.m implements the NMFIG algorithm.

For documentation on NMFIG type in the matlab console:
```matlab
>> help nmfig
```

The functions in the folder ./aux_tools are downloaded from [multi_nmf](http://www.irisa.fr/metiss/ozerov/Software/multi_nmf_toolbox.zip). 

### Slides

**[perception/inria.fr](https://team.inria.fr/perception/research/nmfig/)**








