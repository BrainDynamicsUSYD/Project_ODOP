function [resampled_signal, x, y, t] = resampling_2D(signal, ...
                                    x_experiment, y_experiment, t_experiment, params) 
%% resampling_2D.m
%
% Resampling 2D signal. This function interpolates the signal to a desired 
% resolution based on params.Nkx, params.Nky, and params.Nw. It then pads 
% zeros to the interpolated signal to make it centered at t=0, x=0, and
% y=0.
%
% Inputs: signal            : array of 2D signal (x,y,t)
%                             size(signal) = [length(x_experiment),
%                                             length(y_experiment),
%                                             length(t_experiment)] 
%         x_experiment      : vector of distance to get signal 
%         t_experiment      : vector of time to get signal 
%         params            : instance of the class loadParameters of the 
%                             toolbox
%
% Output: resampled_signal  : array of interpolated and zero-padded signal.
%                             size(resampled_signal) = [params.Nky, params.Nkx, params.Nw]
%         x                 : vector of new distance x =[-x_end,...,0,...x_end]
%         y                 : vector of new distance y =[-y_end,...,0,...y_end]
%         t                 : vector of new time t =[-t_end,...,0,...t_end]
% 
% Example:
% >> params = loadParameters;
% >> load signal.mat   % assuming the data is stored in this mat file
% >> x_experiment = linspace(0,5,256)*1e-3;  % in mm
% >> y_experiment = linspace(0,5,256)*1e-3;  % in mm
% >> t_experiment = linspace(0,20,256);     % in s
% >> [resampled_signal, x, y, t] = resampling_2D(signal, x_experiment, y_experiment,
%                               t_experiment, params) % gives out the resampled
%                                                       signal and new time/distance vectors 
% 
% James Pang, University of Sydney, 2017

%%
% interpolate experimental response to increase resolution

% if the start of t_experiment has an offset from zero, t_interp becomes of
% size Nw/2 and then add a zero in the middle of t
% if the start of t_experiment is zero, t_interp becomes of size Nw/2 + 1 
% and there is no need to add a zero in the middle of t 
if t_experiment(1)~=0
    t_interp = linspace(t_experiment(1), t_experiment(end), params.Nw/2);
    t = [-t_interp(end:-1:1), 0, t_interp];
else
    t_interp = linspace(t_experiment(1), t_experiment(end), params.Nw/2 + 1);
    t = [-t_interp(end:-1:2), t_interp];
end

% case for x is similar as above
if x_experiment(1)~=0
    if x_experiment(1) == -x_experiment(end)
        x_interp = x_experiment(end)*(2*(0:params.Nkx)/params.Nkx - 1);
        x = x_interp;
    else
        x_interp = linspace(x_experiment(1), x_experiment(end), params.Nkx/2);
        x = [-x_interp(end:-1:1), 0, x_interp];
    end
else
    x_interp = linspace(x_experiment(1), x_experiment(end), params.Nkx/2 + 1);
    x = [-x_interp(end:-1:2), x_interp];
end

% case for y is similar as x
if y_experiment(1)~=0
    if y_experiment(1) == -y_experiment(end)
        y_interp = y_experiment(end)*(2*(0:params.Nky)/params.Nky - 1);
        y = y_interp;
    else
        y_interp = linspace(y_experiment(1), y_experiment(end), params.Nky/2);
        y = [-y_interp(end:-1:1), 0, y_interp];
    end
else
    y_interp = linspace(y_experiment(1), y_experiment(end), params.Nky/2 + 1);
    y = [-y_interp(end:-1:2), y_interp];
end

% creating matrices of time and distance
[x_interp_mat, y_interp_mat, t_interp_mat] = meshgrid(x_interp, y_interp, t_interp);
[x_experiment_mat, y_experiment_mat, t_experiment_mat] = meshgrid(x_experiment, y_experiment, t_experiment);

% 2D interpolation
signal_interp = interp3(x_experiment_mat, y_experiment_mat, t_experiment_mat, signal, ...
                        x_interp_mat, y_interp_mat, t_interp_mat);

% since we increased the resolutions, we need to pad 
% the interpolated signal with zeros                          
signal_padded = padarray(signal_interp, [length(y)-length(y_interp), ...
                         length(x)-length(x_interp), ...
                         length(t)-length(t_interp)], 'pre');

% reassigning the resulting interpolated and zero-padded signal to the 
% final variable resampled_signal                          
resampled_signal = signal_padded;

% Fourier transform requires that the x=0 and t=0 be in an offset to the
% right such that an N-sized vector F, where N is even, needs to have the 
% property F(1) = F(end-1) and F(N/2+1) corresponds to the center (either
% the x=0, y=0, or t=0)
x = x(1:end-1);
y = y(1:end-1);
t = t(1:end-1);
resampled_signal = resampled_signal(1:end-1, 1:end-1, 1:end-1);
