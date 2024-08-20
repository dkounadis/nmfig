function [W,H,Kj,Z] = initNMF(X,J,cPerSrc,fs)
%initNMF   NMF initialization via Binary Masking based on TDOA localization 
%        
% INPUTS
%
%  X       : [F x L x I]     F <freq> x L <frams> complex STFT of I sensors
%
%  J       : [1 x 1]         number of Sources to be separated (a scalar)
%               
%  cPerSrc : [1 x 1]         number of components to be given to each src.
%
%  fs      : [1 x 1]         sampling freq. of the mix-signal in samples/s
%
% OUTPUTS
%
%  W       : [F x K]         array of bases, K = J * cPerSrc 
%
%  H       : [K x L]         array of contributions
%
%  Kj      : {J x 1} x [cPerSrc x 1] indices on W,H corresp. to each source
%
%  Z       : [F x L x J]     binary masks (probabilistic soft masks)
%
% version 4 Aout 2016, 13:48 PM
fprintf('[nmf via tdoa localization] Aug. 4 2016\n');

% [1 x 1] number of possible TDOA, single-iteration EM
K = 1001;    micDist = 2;    F = size(X,1);

% [K x 1] grid of TDOA, via possible angles (azimuths) in -pi/2, pi/2
tau = micDist / 340.29 * sin( linspace(-pi/2,pi/2,K) );

% [F x L] assure denominator deflation
Phase = sign(X(:,:,1))./sign( X(:,:,2) + 1e-17 );

% [F x K] candidate GMM means, 2F-1 is the win_len
mu = exp( -2j * pi * fs / (2*F-1) * (0:F-1)' * tau );

% [F x L x K] (Phase-mu), allow a singleton for L
d = bsxfun(@minus, Phase, permute(mu,[1 3 2]) );

% [F x L x K] exponents of the responsibilities, sigma = 0.3 who knows; 
d = - d.*conj(d) / .3;

% [F x L x K] responsibilities
d = exp( bsxfun(@minus, d, log(sum(exp(d),3))) );

% [1 x 1 x K] update priors exponent, because of sum on FL p may be > 0
p = log( sum(sum(d,2),1) );

% [1 x 1 x K] assure negative values
p = p - max(p,[],3);

% [1 x 1 x K] normalize
p = exp( p - log(sum(exp(p),3)) );

% [? x 1] (indices in p)
[~,ind] = findpeaks( p(:) );

% [? x 1] supplement TDOA (if found peaks are less than J)
ind = [ ind ; randi([1 K], J-numel(ind), 1) ];

% [? x 1] sort ind based on proportions p
[~,srcPk] = sort( p(ind) , 'descend' );

% [F x L x J] responsibilities of the J peaks that have the largest p
Z = d(:,:, ind( srcPk(1:J) ) );

% [F x L x J] renormalise
Z = bsxfun(@rdivide,Z,sum(Z,3));

% NMF parameters via KL and mic-1
[W,H,Kj] = Init_KL_NMF_fr_sep_sources( bsxfun(@times,Z,X(:,:,1) ) ,cPerSrc );


