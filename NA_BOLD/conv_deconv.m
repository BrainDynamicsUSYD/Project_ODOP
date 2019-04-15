function [deconvResponses_Wiener_2D,BOLD_noisy,deconvResponses_Wiener_2D_noisyBOLD] = conv_deconv(deconvResponses_Forward_2D,x,y,t,params,V) 

% function to do convolution and deconvolution
% Marilia M Oliveira, 2019


% choose type of noise to use
% range of noise, varying with voxel volume
%[noise] = SNR(V); % adding noise
% 1% of white noise
noise = 0.01; 
params.noise = noise;

deconvResponses_Wiener_2D = deconvolution_HybridWiener_2D(deconvResponses_Forward_2D.reconvBOLD, ...
                                             x, y, t, params); 

% noise: BOLD response
noise_term_BOLD = noise*randn(size(deconvResponses_Forward_2D.reconvBOLD(:,:,:)));

BOLD_noisy = deconvResponses_Forward_2D.reconvBOLD + noise_term_BOLD;

deconvResponses_Wiener_2D_noisyBOLD = deconvolution_HybridWiener_2D(BOLD_noisy, ...
                                             x, y, t, params);

end
