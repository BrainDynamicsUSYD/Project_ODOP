function [noise_norm,snr] = norm_noise(signal,Ymax,x,y,t)
% New noise implementation protocol
% Normalizing noise
% Ymax - maximum magnitude of BOLD response

% white noise according to Tryantafyllou et al. using a fix value of V in a
% fine grid.
% res = 0.125;
% V = res^3;
% [noise] = SNR(V);
% params.noise = noise;
% noise_term_BOLD = noise*randn(size(deconvResponses_Forward_2D.reconvBOLD(:,:,:)));
% noise_term_BOLD = noise*randn(size(signal(:,:,:)));
noise_term_BOLD = randn(size(signal(:,:,:)));
noise_term_BOLD_rms = rms(noise_term_BOLD(:));

% convert noise to frequency space via FT
% filter to kc - cutoff point
% convert back to image space using IFT
scale = 0.125; % scale to res 1 mm. Scale varies according to resolution
[filtered_noise,~] = Smoothing_sinc(noise_term_BOLD,x,y,t,scale);

% calculate rms value
noise_rms = rms(filtered_noise(:));

% calculate initial SNR
%[Ymax,~] = max(deconvResponses_Forward_2D.reconvBOLD(:));
snr = Ymax/noise_rms;

% rescale noise to have SNR of 200 at dx = 1 mm
noise_norm = (snr/200)*noise_term_BOLD;