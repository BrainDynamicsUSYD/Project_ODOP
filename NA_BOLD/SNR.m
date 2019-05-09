function [noise] = SNR(V)
% noise
% equation from Triantafyllou et al. (2005) - Comparison of physiological
% noise at 1.5T, 3T and 7T and optimization of fMRI acquisition parameters
% Values of lambda and kappa from: 
% D. Chaimow et al. (2011) - Modeling and analysis of mechanisms underlying
% fMRI-based deconding of information conveyed in cortical columns.
% D. Chaimow et al. (2018) - Optimization of functional MRI for detection,
% deconding and high-resolution imaging of the response patterns of
% cortical columns

% 1/sigma = tSNR(V) = (kappa*V)/(sqrt(1+lambda^2*kappa^2*V^2))

% where:
%    sigma: standard deviation of signal change
%    tSNR: time-course signal to noise ratio
%    V: voxel volume
%          V(0.125 mm) = 0.125^3 = 0.002 mm^3;
%          V(0.25 mm) = 0.25^3 = 0.0156 mm^3;
%          V(0.50 mm) = 0.50^3 = 0.125 mm^3;
%          V(0.75 mm) = 0.75^3 = 0.4219 mm^3;
%          V(1 mm) = 1^3 = 1 mm^3;
%    lambda: field and scanner independent constant governing the relation
%            between temporal SNR and image SNR (lambda = 0.01297)
%    kappa: proportionality constant between volume and image SNR that is
%           field strength and hardware dependent (kappa = 6.641)
% TR = 5.4s

% Chaimow et al. (2011). Probably these are for 3T
lambda = 0.01297;
kappa = 6.641;

% Chaimow et al. (2018) values for 7T
% lambda = 0.0113;
% kappa = 9.9632;

% considering a task/spontaneous_signal ratio of 10
noise = (sqrt(1+lambda^2*kappa^2*V^2))/(10*kappa*V);
end
