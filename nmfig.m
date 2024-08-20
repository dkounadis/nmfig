function [A,S,W,H,o] = nmfig(X,W,H,Kj,maxNumIter)
%NMFIG   Variational EM for audio source separation using NMF and IG priors
%        
% INPUTS
%
%  X   : [F x L x I]         F <freq> x L <frams> complex STFT of I sensors
%
%  W   : [F x K]             initial matrix of bases (of all K components)
%
%  H   : [K x L]             initial matrix of contributions (all K comp.)
%
%  Kj  : {J x 1} x [? x 1]   Kj is a cell-array with J elements (one for 
%                            each source). Let K = size(W,2) be the total
%                            number of components. Any element Kj{j}, j=1:J
%                            of Kj contains the indexes (of the columns of
%                            W and rows of H) that correspond to source j.
%                            For example Kj = { [1 2 3] , [4 5] , [6 7] }
%                            would set the number of sources to J = 3,
%                            and use W(:,[1 2 3]) * H([1 2 3],:) as the NMF
%                            for souce j=1, W(:,[4 5]) * H([4 5],:) as
%                            the NMF of source j=2, W(:,[6 7]) * H(:,[6 7])
%                            as the NMF for source j=3.
%
%  maxNumIter : [1 x 1]      number of EM iterations
%
% OUTPUTS
%
%  A  : [F x I x J]          estimated mixing matrices
%
%  S  : [F x L x J]          estimated monochannel sources
% 
%  W  : [F x K]              estimated matrix of bases
%
%  H  : [K x L]              estimated matrix of contributions
%
%  o  : [K x 1]              IG shape hyperparameters (component relevance)
%
%  References:
%    [1] D. Kounades-Bastian, L. Girin, X. Alameda-Pineda, S. Gannot, 
%        R. Horaud, An Inverse-Gamma Source Variane Prior with Factorized
%        Parametrization for Audio Source Separation, ICASSP, 2016.
%
% version 30 August 2016, 15:32 PM
fprintf('[nmfig] v1.0, Aug. 30 2016, 15:32 PM\n');



%    /\
%   /__\
%  /    \ R C H E T Y P E



%% A   Constants indexSets & functions
[F,L,I] = size(X);   J = numel(Kj);  K = size(W,2);

% [K x 1] inverse-search partition from comp. to sources
jk = zeros(K,1);   for j=1:J, jk( Kj{j} ) = j; end

% [1 x 1] sensor noise variance
v = X(:)' * X(:) * 1e3 / numel(X);

% [F x 1 x K] spectral patterns
W = permute(W,[1 3 2]);

% [1 x L x K] activation
H = permute(H,[3 2 1]);

% [F x L x K] prior variance of C, shape hyperparameter o
u = bsxfun(@times,W,H);   o = 1;

% {F} x [I x J] initialise filters
A = cell(F,1);   A(:) = { ones(I,J) };

% {F} x [J x J] filter square, M-X create it as cell, E-S cast it in array 
U = cellfun(@(A) A'*A, A, 'uniformoutput', false);

% [F x L] norms of X
normX = sum( sum( X .* conj(X) , 2 ) , 3 );

% f(x) multiply an [F x L x J x J] MATRIX A with a [F x L x J] VECTOR b
mtimes4D = @(A,b) sum( bsxfun(@times , A , permute(b,[1 2 4 3]) ) ,4);
%%






for iter = 1:maxNumIter
%   ____
%  |
%  |____
%  |
%  |____ - S    S T E P



%% E-S   Source inference

% UPDATE
%   p        : [F x L x J]     prior precision of S
%   B        : [F x L x J]     unscaled linear of S
%   S        : [F x L x J]     mean of S - source estimate
%   U        : [F x 1 x J x J] filter Gramian / v
%   Vs       : [F x L x J x J] posterior covariance of S

% [F x J x J] was {F} x [J x J] divide by v
U = bsxfun(@rdivide, permute( cat(3,U{:}) , [3 1 2] ) , v );

% [F x 1 x J x J] 
U = permute(U,[1 4 2 3]);

% S covariance

% {J x 1} x [F x L] prior psd of S via NMF
p = cellfun(@(Kj) sum(u(:,:,Kj),3), Kj, 'uniformoutput', false);

% [F x L x J] precision
p = 1./cat(3,p{:});

% [F x L x J x J] allocate S covariance
R = zeros(F,L,J,J);

% [F x L x J x J] set p on R's diagonal via indexingss
R(:,:,1:J+1:J*J) = p;

% [J x J x F x L ] R = diag(p) + U
R = permute( bsxfun(@plus,R,U) , [3 4 1 2] );

% {F x L} x [J x J] cellarise
Vs = cell(F,L);    for ind = 1:F*L, Vs{ind} = R(:,:,ind); end

% {F x L} x [J x J] S covariance
Vs = cellfun(@inv,Vs,'uniformoutput',false);

% [F x L x J x J] make array cat is FLJxJ with dim-2 the columns of a Vs
Vs = permute( reshape( cat(1,Vs{:}) , J,F,L,J) , [2 3 1 4] );

% calculate S

% [F x 1 x J x I] A hermitian, from {F} x [I x J]
B = conj( permute( cat(3,A{:}) , [3 4 2 1]) );

% [F x L x J] A^H * X/v
B = mtimes4D( B, bsxfun(@rdivide,X,v) );

% [F x L x J] Vs x unscaled part of mean
S = mtimes4D(Vs,B);
%%




%       ___
%      /__/\
%      \  \ \
%       \  \ \
%   ___  \  \ \
%  /__/\  \__\ \
%  \  \ \ /  / /
%   \  \ \  / /
%    \  \ \/ /
%     \  \  /
%      \__\/       S T E P



%% u   inference & IG hyperparameters update

% UPDATE
%   C   : [F x L x K] posterior mean of components
%   Qc  : [F x L x K] 2nd order psd of C
%   u   : [F x L x K] u estimate
%   W   : [F x 1 x K] spectral pattern
%   H   : [1 x L x K] activation
%   o   : [1 x 1 x K] IG prior shape hyperparameter

% [F x L x J] B - U/v * S, J dimensionak linear
r = B - mtimes4D(U,S);

% [F x L x K] posterior mean of C
C = u .* r(:,:,jk);

% [F x L x J]  real( diag(Vs*U) ) / p, J dimensional quadratic
r = real( mtimes4D(Vs,U) ) .* p;

% [F x L x K] 2nd psd
R = u .* (1 - u .* r(:,:,jk))   +  C .* conj(C);
        
% [F x L x K] u
u = bsxfun(@rdivide, R + bsxfun(@times,W,H) , o + 1);

% [F x 1 x K] W
W = bsxfun(@rdivide, L .* o ,  sum(bsxfun(@rdivide,H,u),2) );

% [1 x L x K] H
H = bsxfun(@rdivide, F .* o ,  sum(bsxfun(@rdivide,W,u),1) );

% [1 x 1 x K] o
o = F * L ./ sum(sum(log(1 + R./bsxfun(@times,W,H)) ,1) ,2);
%%




%  |\  /|
%  | \/ | - X    S T E P



%% M-X filter & noise update

% UPDATE
%   A   : {F} x [I x J] filter
%   v   : [F  x 1]      sensor noise

% [F x L x I x J] X*S^H, linear terms 
B = bsxfun(@times, X, permute(conj(S),[1 2 4 3]) );

% [I x J x F] sum L
B = permute( sum(B,2) , [3 4 1 2] );

% {F} x [I x J] linear
r = cell(F,1); for f = 1:F, r{f} = B(:,:,f); end

% [F x L x J x J] Vs + S*S^H
B = Vs + bsxfun(@times,S,permute(conj(S),[1 2 4 3]));

% [J x J x F] sum L, JxJ leading
B = permute( sum(B,2) , [3 4 1 2] );

% {F} x [J x J] quadratic
R = cell(F,1); for f = 1:F, R{f} = B(:,:,f); end

% {F} x [I x J] r * inv(R), filters
A = cellfun(@mrdivide, r, R, 'uniformoutput', false);

% {F} x [J x J] Gramian (ES recasts in array)
U = cellfun(@(A) A'*A, A,    'uniformoutput', false);

% [F x 1] residual sum on L
v = cellfun(@(A,r,U,R) -2*A(:)'*r(:) + U(:)'*R(:), A,r,U,R);

% [F x 1] sensor
v = real(normX + v) / (L*I) + 1e-7;

fprintf('pass %d \n',iter);
end




%% compactify user output
A = permute( cat(3,A{:}) , [3 1 2] );
W = permute( W , [1 3 2]);
H = permute( H , [3 2 1]);
o = o(:);
