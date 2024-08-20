% Example of using nmfig.m for the separation of (static) sound sources
clear
close all

% contains auxilliary functions (STFT, KL_NMF, etc.) and is downloaded from
% http://www.irisa.fr/metiss/ozerov/Software/multi_nmf_toolbox.zip
addpath aux_tools/;   mkdir results

J = 3;              % number of sources
cPerSrc = 20;       % number of components per source (NMF dimension)
stft_win_len = 512; % STFT analysis window
maxNumIter = 10;    % number of EM iterations

% {J} x [M x I] load J ground-truth source-images from ./data folder
[y,fs] = arrayfun(@(j) audioread(sprintf('trueSrc%d.wav',j)), 1:J, 'uniformoutput', false);

% [M x I x J] array of ground-truth source images
y = cat(3,y{:});

% [M x I] make the mixture signal by adding the true src-images
x = sum(y,3);     fs = fs{1};

% M is the number of time samples, I is the number of mikes
[M,I] = size(x);

% [F x L x I] calculate the STFT of x
X = stft_multi( transpose(x), stft_win_len);

% write mix in audio file
audiowrite('./results/mix.wav', x/max(abs(x(:))) ,fs);

% initialization of NMF parameters
[W,H,Kj] = initNMF(X,J,cPerSrc,fs);

fprintf('Applying separation ..\n\n');

% call nmfig() to do the separation
[A,S] = nmfig(X,W,H,Kj,maxNumIter);

% [F x L x I x J] estimated STFT of source images
Y = bsxfun(@times,  permute(A,[1 4 2 3])  ,  permute(S,[1 2 4 3]) );

% [I x M x J] reconstruct the (estimated) time-domain src-img
ye = zeros(M,I,J);

for j=1:J
    ye(:,:,j) = transpose( istft_multi( Y(:,:,:,j) , M ) );
end

% [M x I x J] normalise estimated src-img to avoid decliping when writing
ye = bsxfun(@rdivide,ye,max(abs(ye)));

arrayfun(@(j) audiowrite(sprintf('./results/estimatedSrc%d.wav',j),ye(:,:,j),fs), 1:J , 'uniformoutput',false)

fprintf('\nSeparation results are written in ./results\n');
