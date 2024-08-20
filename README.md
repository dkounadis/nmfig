# NMFIG

MATLAB code for the Variational EM algorithm for Source Separation using unfactorised NMF via Inverse-Gamma priors (NMFIG), which is presented in the paper:

```
D. Kounades-Bastian, L. Girin, X. Alameda-Pineda, S. Gannot and R. Horaud, "An inverse-gamma source variance prior with factorized parameterization for audio source separation," 2016 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Shanghai, China, 2016, pp. 136-140, doi: 10.1109/ICASSP.2016.7471652. keywords: {Source separation;Bayes methods;Estimation;Covariance matrices;Inference algorithms;Indexes;Shape;Audio modeling;local Gaussian model;PSD model;audio source separation},
```

### DEMO


Simply run `example.m` in MATLAB

```python
>> example
```

`example.m` will generate a stereo mix with 3 sources (by loading trueSrc1.wav, ..).
Then will call `initNMF()` to provide initialisation for the VEM params (via binary masking).
Then apply `nmfig()` to separate the sources. Separated sources will be saved as .wav files (estimatedSrc1.wav, etc.) in the directory `./results/`. File `nmfig.m` implements the NMFIG algorithm.
The functions in the folder ./aux_tools are downloaded from [multi_nmf](http://www.irisa.fr/metiss/ozerov/Software/multi_nmf_toolbox.zip). 

### SLIDES

Visit - [perception - inria.fr](https://team.inria.fr/perception/research/nmfig/)








