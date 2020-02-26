%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script implements a simplified version of the BM3D algorithm
% for attenuation of additive white Gaussian noise from images.
%
% It requires the mex file  bm3d_thr_double.dll,  which performs
% BM3D hard-thresholding with double-precision floats.
% The dll is compiled for Windows 32-bit
%
% 
%  The general version of the BM3D algorithm is available at
%  http://www.cs.tut.fi/~foi/GCF-BM3D/
%
%
% Reference:
%  K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian,
%  "Image denoising by sparse 3D transform-domain collaborative filtering,"
%  IEEE Trans. Image Process., vol. 16, no. 8, pp. 2080-2095, August 2007.
%  DOI http://dx.doi.org/10.1109/TIP.2007.901238
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2005-2011 Tampere University of Technology. All rights reserved.
% This work should only be used for nonprofit purposes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y_hat = function_BM3D_HT_double(z,sigma)


%%%% Select transforms ('dct', 'dst', 'hadamard', or anything that is listed by 'help wfilters'):
transform_2D_HT_name     =  'haar'; %% transform used for the HT filt. of size N1 x N1
transform_3rd_dim_name   =  'haar'; %% tranform used in the 3-rd dim, the same for HT and Wiener filt.

%%%% Hard-thresholding (HT) parameters:
N1                  = 8;   %% N1 x N1 is the block size used for the hard-thresholding (HT) filtering
Nstep               = 3;   %% sliding step to process every next refernece block
N2                  = 16;  %% maximum number of similar blocks (maximum size of the 3rd dimensiona of a 3D array)
Ns                  = 25;  %% length of the side of the search neighborhood for full-search block-matching (BM)
tau_match           = .35; %% threshold for the block distance (d-distance)
lambda_thr2D        = 0;   %% threshold parameter for the coarse initial denoising used in the d-distance measure
beta                = 1;   %% Kaiser

%%%% Block-matching parameters:
lambda_thr3D        = 1;  %% threshold parameter for the hard-thresholding in 3D DFT domain
thrToIncStep        = 3;
smallLN             = 3;
stepFS              = Nstep*3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Note: touch below this point only if you know what you are doing!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Make parameters compatible with the interface of the mex-functions
%%%%

decLevel                     = 1;    %% dec. levels of the dyadic wavelet 2D transform for blocks (0 means full decomposition, higher values decrease the dec. number)
decLevel3                    = 0;    %% dec. level for the wavelet transform in the 3rd dimension
coefLevelSmallestScale       = 1.07; %% scales the theshold for the smallest-scale-coefficients
coefLevelNextToSmallestScale = 1.01; %% scales the theshold for the next-coarser-scale-coefficients
coefCoarserLevels            = 1.0;  %% scales the theshold for the rest of the coefficients (i.e. the most coarse ones)

vMask = ones(N1,1)*coefCoarserLevels; vMask((end/4+1):end/2)= coefLevelNextToSmallestScale; vMask((end/2+1):end) = coefLevelSmallestScale;
tMask = vMask * vMask'; %% Finally, create the 2D mask of threshold coefficeints that will be used for HT.


[Tfor, Tinv]   = getTransfMatrix(N1, transform_2D_HT_name, decLevel); %% get (normalized) forward and inverse transform matrices

if (strcmp(transform_3rd_dim_name, 'haar') == 1 || strcmp(transform_3rd_dim_name, 'rbio1.1') == 1),
    %%% Fast internal transform is used, no need to generate transform matrices.
    hadper_trans_single_den         = {};
    inverse_hadper_trans_single_den = {};
else
    %%% Create transform matrices. The transforms are later computed by matrix multiplication with them
    for h=[1 2 4 8 16 32];
        [Tfor3rd, Tinv3rd]   = getTransfMatrix(h, transform_3rd_dim_name, decLevel3);
        hadper_trans_single_den{h}         = (Tfor3rd);
        inverse_hadper_trans_single_den{h} = (Tinv3rd');
    end
end

Wwin2D           = kaiser(N1, beta) * kaiser(N1, beta)'; % Kaiser window used in the hard-thresholding part



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Initial estimate by hard-thresholding filtering

y_hat = bm3d_thr_double(z, hadper_trans_single_den, Nstep, N1, N2, lambda_thr2D, lambda_thr3D, tau_match*N1, (Ns-1)/2, sigma, thrToIncStep, (Tfor), (Tinv)', inverse_hadper_trans_single_den, (tMask), Wwin2D, smallLN, stepFS );


return;


function [Tforward, Tinverse] = getTransfMatrix (N, transform_type, dec_levels)
%
% Create forward and inverse transform matrices, which allow for perfect
% reconstruction. The forward transform matrix is normalized so that the
% l2-norm of each basis element is 1.
%
% [Tforward, Tinverse] = getTransfMatrix (N, transform_type, dec_levels)
%
%  INPUTS:
%
%   N               --> Size of the transform (for wavelets, must be 2^K)
%
%   transform_type  --> 'dct', 'dst', 'hadamard', or anything that is
%                       listed by 'help wfilters' (bi-orthogonal wavelets)
%                       'DCrand' -- an orthonormal transform with a DC and all
%                       the other basis elements of random nature
%
%   dec_levels      --> If a wavelet transform is generated, this is the
%                       desired decomposition level. Must be in the
%                       range [0, log2(N)-1], where "0" implies
%                       full decomposition.
%
%  OUTPUTS:
%
%   Tforward        --> (N x N) Forward transform matrix
%
%   Tinverse        --> (N x N) Inverse transform matrix
%

if exist('dec_levels') ~= 1,
    dec_levels = 0;
end

if N == 1,
    Tforward = 1;
elseif strcmp(transform_type, 'hadamard') == 1,
    Tforward    = hadamard(N);
elseif (N == 8) & strcmp(transform_type, 'haar')==1  %  hardcoded transform so that the wavelet toolbox is not needed to generate it
    Tforward = [ 0.500000000000000   0.500000000000000   0.500000000000000   0.500000000000000                   0                   0                   0                   0
        0                   0                   0                   0   0.500000000000000   0.500000000000000   0.500000000000000   0.500000000000000
        0.500000000000000   0.500000000000000  -0.500000000000000  -0.500000000000000                   0                   0                   0                   0
        0                   0                   0                   0   0.500000000000000   0.500000000000000  -0.500000000000000  -0.500000000000000
        0.707106781186547  -0.707106781186547                   0                   0                   0                   0                   0                   0
        0                   0   0.707106781186547  -0.707106781186547                   0                   0                   0                   0
        0                   0                   0                   0   0.707106781186547  -0.707106781186547                   0                   0
        0                   0                   0                   0                   0                   0   0.707106781186547  -0.707106781186547  ];
elseif (N == 8) & strcmp(transform_type, 'bior1.5')==1 % hardcoded transform so that the wavelet toolbox is not needed to generate it
    Tforward =  [ 0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274;
        0.219417649252501   0.449283757993216   0.449283757993216   0.219417649252501  -0.219417649252501  -0.449283757993216  -0.449283757993216  -0.219417649252501;
        0.569359398342846   0.402347308162278  -0.402347308162278  -0.569359398342846  -0.083506045090284   0.083506045090284  -0.083506045090284   0.083506045090284;
        -0.083506045090284   0.083506045090284  -0.083506045090284   0.083506045090284   0.569359398342846   0.402347308162278  -0.402347308162278  -0.569359398342846;
        0.707106781186547  -0.707106781186547                   0                   0                   0                   0                   0                   0;
        0                   0   0.707106781186547  -0.707106781186547                   0                   0                   0                   0;
        0                   0                   0                   0   0.707106781186547  -0.707106781186547                   0                   0;
        0                   0                   0                   0                   0                   0   0.707106781186547  -0.707106781186547];
elseif (N == 8) & strcmp(transform_type, 'dct')==1 % hardcoded transform so that the signal processing toolbox is not needed to generate it
    Tforward = [ 0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274;
        0.490392640201615   0.415734806151273   0.277785116509801   0.097545161008064  -0.097545161008064  -0.277785116509801  -0.415734806151273  -0.490392640201615;
        0.461939766255643   0.191341716182545  -0.191341716182545  -0.461939766255643  -0.461939766255643  -0.191341716182545   0.191341716182545   0.461939766255643;
        0.415734806151273  -0.097545161008064  -0.490392640201615  -0.277785116509801   0.277785116509801   0.490392640201615   0.097545161008064  -0.415734806151273;
        0.353553390593274  -0.353553390593274  -0.353553390593274   0.353553390593274   0.353553390593274  -0.353553390593274  -0.353553390593274   0.353553390593274;
        0.277785116509801  -0.490392640201615   0.097545161008064   0.415734806151273  -0.415734806151273  -0.097545161008064   0.490392640201615  -0.277785116509801;
        0.191341716182545  -0.461939766255643   0.461939766255643  -0.191341716182545  -0.191341716182545   0.461939766255643  -0.461939766255643   0.191341716182545;
        0.097545161008064  -0.277785116509801   0.415734806151273  -0.490392640201615   0.490392640201615  -0.415734806151273   0.277785116509801  -0.097545161008064];
elseif (N == 8) & strcmp(transform_type, 'dst')==1 % hardcoded transform so that the PDE toolbox is not needed to generate it
    Tforward = [ 0.161229841765317   0.303012985114696   0.408248290463863   0.464242826880013   0.464242826880013   0.408248290463863   0.303012985114696   0.161229841765317;
        0.303012985114696   0.464242826880013   0.408248290463863   0.161229841765317  -0.161229841765317  -0.408248290463863  -0.464242826880013  -0.303012985114696;
        0.408248290463863   0.408248290463863                   0  -0.408248290463863  -0.408248290463863                   0   0.408248290463863   0.408248290463863;
        0.464242826880013   0.161229841765317  -0.408248290463863  -0.303012985114696   0.303012985114696   0.408248290463863  -0.161229841765317  -0.464242826880013;
        0.464242826880013  -0.161229841765317  -0.408248290463863   0.303012985114696   0.303012985114696  -0.408248290463863  -0.161229841765317   0.464242826880013;
        0.408248290463863  -0.408248290463863                   0   0.408248290463863  -0.408248290463863                   0   0.408248290463863  -0.408248290463863;
        0.303012985114696  -0.464242826880013   0.408248290463863  -0.161229841765317  -0.161229841765317   0.408248290463863  -0.464242826880013   0.303012985114696;
        0.161229841765317  -0.303012985114696   0.408248290463863  -0.464242826880013   0.464242826880013  -0.408248290463863   0.303012985114696  -0.161229841765317];
elseif strcmp(transform_type, 'dct') == 1,
    Tforward    = dct(eye(N));
elseif strcmp(transform_type, 'dst') == 1,
    Tforward    = dst(eye(N));
elseif strcmp(transform_type, 'DCrand') == 1,
    x = randn(N); x(1:end,1) = 1; [Q,R] = qr(x);
    if (Q(1) < 0),
        Q = -Q;
    end;
    Tforward = Q';
else %% a wavelet decomposition supported by 'wavedec'
    %%% Set periodic boundary conditions, to preserve bi-orthogonality
    dwtmode('per','nodisp');

    Tforward = zeros(N,N);
    for i = 1:N
        Tforward(:,i)=wavedec(circshift([1 zeros(1,N-1)],[dec_levels i-1]), log2(N), transform_type);  %% construct transform matrix
    end
end

%%% Normalize the basis elements
Tforward = (Tforward' * diag(sqrt(1./sum(Tforward.^2,2))))';

%%% Compute the inverse transform matrix
Tinverse = inv(Tforward);

return;
