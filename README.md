# NMFIG

##  <a src="https://inria.hal.science/INRIA/public/HAL_logotype_rvb_fond_clair_fr_1.png" href="https://inria.hal.science/hal-01253169v1/document"> <img src="https://inria.hal.science/INRIA/public/HAL_logotype_rvb_fond_clair_fr_1.png"></a>

MATLAB code for the Variational EM algorithm for Source Separation via `non-factorised` full-rank PSD (NMF applies in the Bayesian prior) via Inverse-Gamma priors from the paper:

```
D. Kounades-Bastian, L. Girin, X. Alameda-Pineda, S. Gannot and R. Horaud, 
"An inverse-gamma source variance prior with factorized parameterization for audio source separation,"
2016 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Shanghai, China, pp. 136-140, doi: 10.1109/ICASSP.2016.7471652.
```

### DEMO

```python
# In Matlab command line
>> example
```

`example.m` generates a stereo mix of 3 sources (by loading trueSrc1.wav, ..).
Then calls `initNMF()` to provide initialisation for the Variational EM params (via binary masking).
Then applies `nmfig()` to separate the sources. 

Separated sources will be saved as .wav files (estimatedSrc1.wav, etc.) in the directory `./results/`. File `nmfig.m` implements the NMFIG algorithm.

### SLIDES

Visit - [perception - inria.fr](https://team.inria.fr/perception/research/nmfig/)








