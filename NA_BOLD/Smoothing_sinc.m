function [filtered_coord,vWidth] = Smoothing_sinc(BOLD,x,y,t,scale)
% Smoothing of image.
% FFT image to transfer from image space to k-space. Apply rectangular
% filter in k-space (works as an sinc filter in image space). IFFT back to
% image space. Filter cuttoff point at +-kc.

Nkx = length(x);
Nky = length(y);
Nw = length(t);

dx = mean(diff(x));
dy = mean(diff(y));
dt = mean(diff(t));

% calculating sampling rate in frequency space 
kxsamp = (1/dx)*2*pi;  % in 1/m
kysamp = (1/dy)*2*pi;  % in 1/m
wsamp = (1/dt)*2*pi;   % in 1/s

% creating vectors of frequencies kx, ky, and w
[kx, ky, w] = generate_kw_2D(kxsamp, kysamp, wsamp, Nkx, Nky, Nw);

% new voxel width
vWidth = dx/scale;

% cutoff point
kc = 2*pi/(2*vWidth);
kc_xind = dsearchn(kx', [-kc, kc]');
kc_yind = dsearchn(ky', [-kc, kc]');

% converting the signal to frequency space via fourier transform
signal_freq = coord2freq_2D(BOLD, kx, ky, w);

% rectangular filter
filter = zeros(size(BOLD));
filter(kc_yind(1):kc_yind(2), kc_xind(1):kc_xind(2), :) = 1;

% calculting filtered signal in frequency space
filtered_freq = signal_freq.*filter;

% converting the filtered signal to coordinate space via inverse fourier transform
filtered_coord = real(freq2coord_2D(filtered_freq, kx, ky, w));

% % getting positive part of image
% [row,col] = size(filtered_coord(:,:,1));
% filtered_coord1 = filtered_coord((row/2):end-1,(col/2):end-1,:);
% 
% % applying hann window
% W = size(filtered_coord1,1);
% window = hann(W,'periodic');
% filtered_coordW = filtered_coord1.*window;

end

