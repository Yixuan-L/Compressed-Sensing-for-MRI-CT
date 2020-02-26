% ------------------------------------------------------------------------------------------
%
%     Demo software for Compressed Sensing Image Reconstruction
%               Public release ver. 1.0 (1 March 2011)
%
% ------------------------------------------------------------------------------------------
%
% The software reproduces the experiments published in the paper
%
%  K. Egiazarian, A. Foi, and V. Katkovnik,
%  "Compressed Sensing Image Reconstruction via Recursive Spatially Adaptive Filtering,"
%  Proc. IEEE ICIP 2007, San Antonio (TX), USA, pp. 549-552, Sept. 2007.
%  DOI http://dx.doi.org/10.1109/ICIP.2007.4379013
%
% ------------------------------------------------------------------------------------------
%
% authors:               Karen Egiazarian
%                        Alessandro Foi
%
% web page:              http://www.cs.tut.fi/~comsens
%
% contact:               firstname.lastname@tut.fi
%
% ------------------------------------------------------------------------------------------
% Copyright (c) 2006-2011 Tampere University of Technology.
% All rights reserved.
% This work should be used for nonprofit purposes only.
% ------------------------------------------------------------------------------------------
%
%
% Description
% -----------
%
% This demo software reproduces the following four experiments from the paper
%
%  K. Egiazarian, A. Foi, and V. Katkovnik,
%  "Compressed Sensing Image Reconstruction via Recursive Spatially Adaptive Filtering,"
%  Proc. IEEE ICIP 2007, San Antonio (TX), USA, pp. 549-552, Sept. 2007.
%  DOI http://dx.doi.org/10.1109/ICIP.2007.4379013
%
% Experiments 1-3
% Tomographic reconstruction from sparse projections
% (approximating Radon projections as radial lines in FFT domain)
%  We consider three illustrative inverse problems of compressed sensing for computerized
%  tomography. In particular, we show reconstruction examples of the Shepp-Logan phantom
%  from sparse Radon projections, with 22 and 11 radial lines in FFT-domain
%  (i.e., available Radon projections), and reconstruction from limited-angle projections,
%  with a reduced subset of 61 projections within a 90 degrees aperture. The available
%  portions of the spectrum and the initial back-projection estimates are shown below.
%
%   experimentCase=1  - Sparse projections: 22 radial lines
%   experimentCase=2  - Sparse projections: 11 radial lines
%   experimentCase=3  - Limited-angle (90 degrees aperture, 61 radial lines)
%
%
% Experiment 4
% Reconstruction from low-frequency portion of Fourier spectrum
%  Reconstruction of the Cameraman image (256x256 pixels) from the low-frequency portion
%  of its Fourier spectrum (a 128×128 square centered at the DC).
%
%   experimentCase=4
%
%
%
% For all the above experiments, the block-matching and 3D filtering algorithm (BM3D)
% serves as the spatially adaptive filter used in the recursive stochastic approximation
% algorithm. The separable 3D Haar wavelet decomposition is adopted as the transform
% utilized internally by the BM3D algorithm.
%
% The general version of the BM3D algorithm is available at
% http://www.cs.tut.fi/~foi/GCF-BM3D/
%
%  K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian,
%  "Image denoising by sparse 3D transform-domain collaborative filtering,"
%  IEEE Trans. Image Process., vol. 16, no. 8, pp. 2080-2095, August 2007.
%  DOI http://dx.doi.org/10.1109/TIP.2007.901238
%
%
% ------------------------------------------------------------------------------------------
%
% Disclaimer
% ----------
%
% Any unauthorized use of these routines for industrial or profit-oriented activities is
% expressively prohibited. By downloading and/or using any of these files, you implicitly
% agree to all the terms of the TUT limited license (included in the file Legal_Notice.txt).
% ------------------------------------------------------------------------------------------
%

%%

clear all
% close all
% warning off all



%% options for this demo script

experimentCase=2;       %  experimentCase=1  - Sparse projections: 22 radial lines
%                       %  experimentCase=2  - Sparse projections: 11 radial lines
%                       %  experimentCase=3  - Limited-angle (90 degrees aperture, 61 radial lines)
%                       %  experimentCase=4  - Reconstruction from low-frequency portion of Fourier spectrum

upperHalfOnly=0;        %  Keep only upper half of Fourier plane, like in the experiments from the L1Magic package: http://www.acm.caltech.edu/l1magic/

downsampleFactor=3;     %  Downsampling factor in order to use smaller image size for much faster examples (default =1, i.e. no downsampling).



figure_every_n_seconds=2;     % update figure every n seconds; set to negative to disable figures
text_every_n_seconds=1;       % printout text every n seconds; set to negative to disable text during iterations
number_of_digits_text=4;      % number of digits in printout
number_of_digits_fig=4;       % number of digits in figures
figTitleFontSize=18;          % font size in figures




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  END OF OPTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













%% experiment setup

if any(experimentCase==[1 2 3])
    n=round(256/downsampleFactor);
    theta=phantom(n);      % create Sheep-Logan phantom of size n x n
elseif experimentCase==4
    theta=im2double(imread('cameraman256.png'));  % image must have size n x n
    if downsampleFactor>1
        theta=imresize(theta,1/downsampleFactor);
    end
    n=size(theta,1);
end

if experimentCase==1        %% Sparse projections: 22 radial lines
    L=22;                   % number of radial lines in the Fourier domain
    aperture=(pi/180)*180;    % aperture encompassing all scanning angles (aperture<pi is limited angle)
    direction=(pi/180)*0;     % direction of the scanning beam (middle axis)
    k_final=10000;            % number of iterations
    alpha=(1+1/300)^2;        % parameter for exponential decay of excitation noise
    beta=500;                 % parameter for exponential decay of excitation noise

elseif experimentCase==2    %% Sparse projections: 11 radial lines
    L=11;                   % number of radial lines in the Fourier domain
    aperture=(pi/180)*180;    % aperture encompassing all scanning angles (aperture<pi is limited angle)
    direction=(pi/180)*0;     % direction of the scanning beam (middle axis)
    k_final=10000;            % number of iterations
    alpha=(1+1/300)^2;        % parameter for exponential decay of excitation noise
    beta=500;                 % parameter for exponential decay of excitation noise

elseif experimentCase==3    %% Limited-angle: aperture 90 degrees, 61 radial lines
    L=61;                   % number of radial lines in the Fourier domain
    aperture=(pi/180)*90;     % aperture encompassing all scanning angles (aperture<pi is limited angle)
    direction=(pi/180)*45;    % direction of the scanning beam (middle axis)
    k_final=100000;           % number of iterations
    alpha=(1+1/3000)^2;       % parameter for exponential decay of excitation noise
    beta=8000;              % parameter for exponential decay of excitation noise

elseif experimentCase==4    %% Reconstruction from low-frequency portion of Fourier spectrum
    n_small=n/2;              % size of central low-frequency portion
    k_final=62;               % number of iterations
    alpha=(1+1/300)^2;        % parameter for exponential decay of excitation noise
    beta=1120;                % parameter for exponential decay of excitation noise
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  NOTHING TO MODIFY HERE BELOW  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% build mask of spectrum coefficients
if any(experimentCase==[1 2 3])
    S=LineMaskLimitedAngle(L,n,aperture,direction);
elseif experimentCase==4
    S=zeros(n);
    S([1:ceil(n_small/2),n-floor(n_small/2)+1:n],[1:ceil(n_small/2),n-floor(n_small/2)+1:n])=1;
end

if upperHalfOnly
    S=fftshift(S)';
    S((n/2)*(n+1)+2:end)=0;
    S=ifftshift(S');
end


%% generate observations

y=fft2(theta);
Omega=find(S);
y_hat_0=zeros(n);
y_hat_0(Omega)=y(Omega);

if ~upperHalfOnly
    theta_hat_0=real(ifft2(y_hat_0));
else
    theta_hat_0=real(ifft2(y_hat_0*2))-y_hat_0(1)/(n*n);
end
clear y


%% theta is not needed below this point (can be kept for computing PSNR)
% clear theta

%% Printout experiment info
disp(' ');
disp('-------------------------------------------------------------------------------')
disp(['Image size ',num2str(size(theta_hat_0,2)),'x',num2str(size(theta_hat_0,1)),' pixel']);
disp(['Available ',num2str(numel(Omega)),' FFT coefficients out of ',num2str(numel(theta_hat_0)),' (',num2str(numel(Omega)/numel(theta_hat_0)*100,3),'%)'])

if exist('theta','var')
    err_theta_0=10*log10(1/(mean((theta_hat_0(:)-theta(:)).^2)));   %%% PSNR
    disp('-------------------------------------------------------------------------------')
    disp(['PSNR(dB) back-projection estimate theta_hat_0   : ',num2str(err_theta_0,number_of_digits_text)])
end
disp('-------------------------------------------------------------------------------')

if figure_every_n_seconds>=0
    figureHandle0=figure;  % initializes figures and handles
    figureHandle=figure;   % initializes figures and handles
    figure(figureHandle0);
    subplot(1,3,1)
    imshow(fftshift(S),[0,1]), title([{' available portion $\Omega$'};{'of the FFT spectrum'}],'interpreter','latex','fontsize',figTitleFontSize)
    subplot(1,3,2)
    imshow([theta_hat_0],[0 1]);
    if exist('theta','var')
        title(['$\hat{\theta}^{(0)}$, PSNR: ',num2str(err_theta_0,number_of_digits_fig)],'interpreter','latex','fontsize',figTitleFontSize)
    else
        title(['$\hat{\theta}^{(0)}$'],'interpreter','latex')
    end
end

%% excitation schedule
% std_excite is std of excitation noise to be inserted into the data
std_excite=sqrt(alpha.^(-(1:k_final)-beta));

maxPSNR=255;    % maximum value of PSNR that can be achieved with the current numerical precision (double) (default: maxPSNR=255)
min_std_excite=10^(-maxPSNR/20)/sqrt(1-numel(Omega)/numel(theta_hat_0));  % minimum standard-deviation value allowable with regard to precision errors (relevant when std_excite is approaching zero)
k_final2maxPSNR=ceil(-log(min_std_excite^2)/log(alpha)-beta);  % number of iterations necessary for std_excite to reach min_std_excite

disp([' Variance of excitation noise is ',num2str(alpha,6),'^(-k-',num2str(beta,6),')']);
disp([' ~',num2str(k_final2maxPSNR),' iterations necessary to reach minimal excitation']);
disp(['  (min. st.dev.  = ',num2str(min_std_excite),',  precision ~',num2str(maxPSNR),'dB)'])
disp('-------------------------------------------------------------------------------')
disp([' Total number of iterations is set to  k_final = ',num2str(k_final)]);
disp('-------------------------------------------------------------------------------')

%% initializes variables for error criteria
if exist('theta','var')
    err_theta_hat=zeros(k_final,1);
    err_filtered=zeros(k_final,1);
    err_excitation=zeros(k_final,1);
end

%% initialiazation of recursive system
randn('seed', 0); rand('seed', 0);  % fix seed of pseudo-random noise (optional)

k=0;                       % k=0
theta_hat_k=theta_hat_0;   % the initial estimate is the back-projection estimate
time_counter=now*[1 1 1];  % time counter used for controlling screen updates and execution time


%%
while k<k_final,

    %%%% ---------------------------------------------------------------------------------------------------------------------------
    %%%% ------------------  RECURSION BEGINS  -------------------------------------------------------------------------------------
    %%%% ---------------------------------------------------------------------------------------------------------------------------
    %%%%    variables used in the recursion:  Omega
    %%%%                                      y_hat_0
    %%%%                                      theta_hat_k
    %%%%                                      std_excite

    k=k+1;

    %% std_filt is std of noise to be removed by the filtering
    if  k==1   %% initialize std_filt using MAD
        det_coefs=conv2(conv2(theta_hat_k,[-1 1],'valid'),[-1;1],'valid')/(2*0.6745);
        std_filt=median(abs(det_coefs(:)));        clear det_coefs
    else   %%  use std of excitation noise
        std_filt=std_excite(k-1);
    end

    %% spatially adaptive filtering with BM3D
    filtered_invT_y_hat_k=function_BM3D_HT_double(theta_hat_k,std_filt);

    if k<k_final
        std_excite(k)=max(min_std_excite,std_excite(k));   %% to take care of precision errors when std_excite is approaching zero
    else
        std_excite(k)=0;  %% stop excitation noise at the last iteration (to avoid noise in final estimate)
    end

    %% Add excitation noise
    y_hat_k=fft2(filtered_invT_y_hat_k + std_excite(k)*randn(size(filtered_invT_y_hat_k)));
    %% Substitute FFT coefficients by those given
    y_hat_k(Omega)=y_hat_0(Omega);

    theta_hat_k=real(ifft2(y_hat_k));   % image estimate

    %%%% ------------------------------------------------------------------------------------------------------------------------------
    %%%% ------------------  RECURSION ENDS  ------------------------------------------------------------------------------------------
    %%%% ------------------------------------------------------------------------------------------------------------------------------


    %% secondary scripts for printout and figure  (removed from main script for the sake of clarity)
    subscript_text_iter    % printout printed very text_every_n_seconds
    subscript_figure_iter  % figure which is updated every figure_every_n_seconds

end


%%
theta_hat_k=min(1,max(0,theta_hat_k));  % constrain between 0 and 1 after last substitution

%% computes and displays errors
if exist('theta','var')
    err_theta_hat(k)=10*log10(1/(mean((theta_hat_k(:)-theta(:)).^2)));                   %%% PSNR
    disp('-------------------------------------------------------------------------------')
    if text_every_n_seconds>=0  % no need to print again this, unless there is regular printout during the iterations
        disp(['PSNR(dB) back-projection estimate theta_hat_0   : ',num2str(err_theta_0,number_of_digits_text)])
    end
    disp(['PSNR(dB) final image estimate theta_hat_k_final : ',num2str(err_theta_hat(k),number_of_digits_text)])
end
disp('-------------------------------------------------------------------------------')
disp(' ');

%% update first figure and create figure with log spectra
if figure_every_n_seconds>=0
    figure(figureHandle0);
    subplot(1,3,3)
    imshow([theta_hat_k],[0 1]);
    if exist('theta','var')
        title(['$\hat{\theta}^{(k_{\rm{final}})}$, PSNR: ',num2str(err_theta_hat(k),number_of_digits_fig)],'interpreter','latex','fontsize',figTitleFontSize)
    else
        title(['$\hat{\theta}^{(k_{\rm{final}})}$'],'interpreter','latex','fontsize',figTitleFontSize)
    end

    figure
    subplot(1,2,1)
    logAbsRange=[min(log(abs(y_hat_k(:)))) max(log(abs(y_hat_k(:))))];
    imshow(log(abs(fftshift(y_hat_0))),logAbsRange)
    title(['$\mathcal{T}(\hat{\theta}^{(0)})$'],'interpreter','latex','fontsize',figTitleFontSize)
    subplot(1,2,2)
    imshow(log(abs(fftshift(y_hat_k))),logAbsRange)
    title(['$\mathcal{T}(\hat{\theta}^{(k_{\rm{final}})})$'],'interpreter','latex','fontsize',figTitleFontSize)
end

%% end of demo software

