% Code used for the PLOS paper

% Marilia M Oliveira, 2018, 2019.

tic

params = loadParameters;

xmax = 16; ymax = 16; tmax = 23;
dx = 0.1; dy = 0.1; dt = 0.2;
x_experiment = (0:dx:xmax)*1e-3;       % need to multiply by 1e-3 to convert to mm
y_experiment = (0:dy:ymax)*1e-3;       % need to multiply by 1e-3 to convert to mm
t_experiment = 0:dt:tmax;
fontT = 20;
fontS = 18;
fontL = 2;
tick = 0.05;
line = 1;
cmin = -1;
cmax = 1;
y_axisy = 0.8;
y_axisx = 0.8;

cmap = jet;

[xx, yy, tt] = meshgrid(x_experiment, y_experiment, t_experiment);

% Pinwheel angle and positions
TH = [ 0,  45,  90,  135]; 
XC = [7.8, 7.5, 7.2, 7.5];
YC = [7.5, 7.8, 7.5, 7.2];

% stimulation for one stimulus only
orientation = [TH(1)];
positionx = [XC(1)]; 
positiony = [YC(1)];

ts = 3;
ton = 7;
toff = 10;


zz = zeros(size(xx)); % stimulus
I1z01 = zeros(1,size(positionx,2));
I2z01 = zeros(1,size(positionx,2));
tstart1 = zeros(1,size(positionx,2));
tend1 = zeros(1,size(positionx,2));

jj = 1;
figure('Name','neural activity simulated','NumberTitle','on'),
for ll = 1:size(orientation, 2)
    
    % Generation of patchy connections
    A1=orientation(ll);
    xp1=positionx(ll);
    yp1=positiony(ll);
    s_x1=8;
    s_y1=1;
    a=2;
    [z0]=pearl(xx*1e3,yy*1e3,A1,xp1,yp1,s_x1,s_y1,a);
    z01=z0(:,:,1);
    
    subplot(round(size(orientation,2)),round(size(orientation,2)),ll)
    imagesc(x_experiment*1e3, y_experiment*1e3, z01);
    set(gca, 'ydir', 'normal')
    pbaspect([1 1 1])
    title(['neural activity ',int2str(ll)]); 
    PlotNames    
    
    % stimulus over time
    tstart = ts + ((ton+toff)*(ll-1));
    tend = tstart + ton;
    tstart1(jj) = tstart(1);
    tend1(jj) = tend(1);
    ton_all = tstart1(jj):tend1(jj)-1;
    ton_all1(jj,:) = ton_all;
    
    tstart_ind = dsearchn(t_experiment',tstart);
    tend_ind = dsearchn(t_experiment',tend);
    
    zz(:,:,tstart_ind:tend_ind) = z0(:,:,1:tend_ind-tstart_ind+1);
    
    % maximum position of stimulus
    [max_z0, I1z0, I2z0, I3z0] = max_voxel(z0);
    I1z01(jj) = I1z0; 
    I2z01(jj) = I2z0;
    
    jj = jj+1;
end

%% Convolution and deconvolution of neural activity w/ and w/o noise

params.Gamma = 0.8; % wave damping rate
params.v_b = 0.5e-3; % propagation speed

params.Nkx = 2^8;
params.Nky = 2^8;
params.Nw = 2^6;
[resampled_zz, x, y, t] = resampling_2D(zz, ...
                                        x_experiment, y_experiment, t_experiment, params);
t0 = dsearchn(t', 0);
x0 = dsearchn(x', 0);
y0 = dsearchn(y', 0);

V = 0.002; % voxel volume (0.125^3 mm^3)

deconvResponses_Forward_2D = deconvolution_Forward_2D(resampled_zz, ...
                                             x, y, t, params);

% convolution, Wiener deconvolution and noise addition realized in function
% choose in the conv_deconv function the type of noise to use
% [deconvResponses_Wiener_2D,BOLD_noisy,deconvResponses_Wiener_2D_noisyBOLD] = conv_deconv(deconvResponses_Forward_2D,x,y,t,params,V); 

% calculating maximum time at BOLD Forw
BOLD_max = deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,:);
[Ymax,Imax] = max(BOLD_max(:));
[I1max,I2max,I3max] = ind2sub(size(BOLD_max),Imax);

time_spe = I3max;

% calculating maximum time at Neural activity
NA_max = deconvResponses_Forward_2D.neural(x0:end,y0:end,:);
[NAmax,INAmax] = max(NA_max(:));
[I1NAmax,I2NAmax,I3NAmax] = ind2sub(size(NA_max),INAmax);

time_NA = I3NAmax;

%fig = figure('Position', [160 278 640 320]);
fig = figure('Position', [160 278 750 380]);
subplot(1,2,1)
imagesc(x(x0:end)*1e3,y(y0:end)*1e3,deconvResponses_Forward_2D.neural(x0:end,y0:end,time_NA))
title({['NA','. t: ',int2str(t(time_NA)),'s'];['Res.: 0.125 mm.']},'fontsize',fontT);
set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
    'ylim', [y(y0)*1e3, y(end)*1e3]);
pbaspect([1 1 1])
xlabel('x (mm)', 'fontsize', fontS)
xticklabels({'0','5','10','15'})
ylabel('y (mm)', 'fontsize', fontS)
colormap(cmap)
colorbar
subplot(1,2,2)
imagesc(x(x0:end)*1e3,y(y0:end)*1e3,deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,time_spe))
title({['BOLD','. t: ',int2str(t(time_spe)),'s'];['Res.: 0.125 mm.']},'fontsize',fontT);
set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
    'ylim', [y(y0)*1e3, y(end)*1e3]);
pbaspect([1 1 1])
xlabel('x (mm)', 'fontsize', fontS)
xticklabels({'0','5','10','15'})
ylabel('y (mm)', 'fontsize', fontS)
colormap(cmap)
colorbar

%% Profile of max voxel along time 

[vectorf_na,~, ~,~,~] = max_voxel(deconvResponses_Forward_2D.neural);
[vectorf_bold,~, ~,~,~] = max_voxel(deconvResponses_Forward_2D.reconvBOLD);

figure('Name','Profile of maximum voxel along time - Res.: 0.125 mm','NumberTitle','on', ...
    'Position', [160 278 750 380])
subplot(1,2,1)
plot(t(t0:end),vectorf_na(t0:end,:),'k-','LineWidth',fontL);
set(gca, 'fontSize', fontS, 'xlim', [0 tmax],'TickLength',[tick tick],'LineWidth',line);
pbaspect([1 1 1])
xlabel('Time (s)', 'fontsize', fontS)
ylabel('Y(r,t)', 'fontsize', fontS)
title({['NA_{max} voxel along t'];['Res.: 0.125 mm']},'fontsize',fontT)
subplot(1,2,2)
plot(t(t0:end),vectorf_bold(t0:end,:),'k-','LineWidth',fontL);
set(gca, 'fontSize', fontS, 'xlim', [0 tmax],'TickLength',[tick tick],'LineWidth',line);
pbaspect([1 1 1])
xlabel('Time (s)', 'fontsize', fontS)
ylabel('Y(r,t)', 'fontsize', fontS)
title({['BOLD_{max} voxel along t'];['Res.: 0.125 mm']},'fontsize',fontT)

%% Plot BOLD profile x,y 0.125mm and modulation

if orientation == 0 || orientation == 90
    
    figure('Name','BOLD profile x,y direction - Res.: 0.125 mm','NumberTitle','on', ...
        'Position', [160 278 750 380]),
    
    [matrix_rowx,matrix_coly] = profile(deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,time_spe));
    [matrix_rowx_NA,matrix_coly_NA] = profile(deconvResponses_Forward_2D.neural(x0:end,y0:end,time_spe-2));
    
    [visibility_BOLD,Maximum] = modulation_na(matrix_rowx);
    [visibility_NA,Maximum_NA] = modulation_na(matrix_rowx_NA);
    
        subplot(1,3,1)
        plot(y(y0:end)*1e3,matrix_rowx, 'k-','LineWidth',fontL)
        title({['BOLD - Prof x'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximum+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])

        subplot(1,3,2)
        plot(x(x0:end)*1e3,matrix_coly, 'k-','LineWidth',fontL)
        title({['BOLD - Prof y'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximum+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        subplot(1,3,3)
        plot(y(y0:end)*1e3,matrix_rowx_NA, 'k-','LineWidth',fontL)
        title({['NA - Prof x'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA rec', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximum_NA+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
end

%% Plot BOLD profile diagonal res 0.125mm and modulation

if orientation == 45 || orientation == 135
    
    figure('Name','BOLD profile diagonal direction - Res.: 0.125 mm','NumberTitle','on', ...
        'Position', [160 278 750 380]),

    [~, matrix_dicr, matrix_didecr] = Diagonal(deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,time_spe));
    [~, matrix_dicr_NA, matrix_didecr_NA] = Diagonal(deconvResponses_Forward_2D.neural(x0:end,y0:end,time_spe-2));
    % matrix_dicr = diag(fliplr(matrix)); % diagonal crescent: left top to
    % right bottom  \  --> use for OP 135 degree
    % matrix_didecr = diag(matrix); % diagonal decrescent: left bottom
    % to right top / --> use for OP 45 degree
    [visibilityd_BOLD,Maximumd] = modulation_na(matrix_didecr);
    [visibilityd_NA,Maximumd_NA] = modulation_na(matrix_didecr_NA);

    subplot(1,3,1)
    plot(x(x0:end)*1e3,matrix_dicr,'k-','LineWidth',fontL)
    title({['BOLD - Prof diag \'];['Res: 0.125 mm']},'fontsize',fontT)
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('BOLD', 'fontsize', fontS)
    set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximumd+0.05)],'TickLength',[tick tick],'LineWidth',line);
    pbaspect([1 1 1])

    subplot(1,3,2)
    plot(x(x0:end)*1e3,matrix_didecr,'k-','LineWidth',fontL)
    title({['BOLD - Prof diag /'];['Res: 0.125 mm']},'fontsize',fontT)
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('BOLD', 'fontsize', fontS)
    set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximumd+0.05)],'TickLength',[tick tick],'LineWidth',line);
    pbaspect([1 1 1])
    
    subplot(1,3,3)
    plot(x(x0:end)*1e3,matrix_didecr_NA,'k-','LineWidth',fontL)
    title({['NA - Prof diag /'];['Res: 0.125 mm']},'fontsize',fontT)
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('NA rec', 'fontsize', fontS)
    set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximumd_NA+0.05)],'TickLength',[tick tick],'LineWidth',line);
    pbaspect([1 1 1])

end

%% Noise - normalizing

% Normalizing the noise
[noise_norm,snr] = norm_noise(deconvResponses_Forward_2D.reconvBOLD,Ymax,x,y,t);
noise_norm2 = noise_norm*(1/0.125); 
BOLD_noise =  deconvResponses_Forward_2D.reconvBOLD + noise_norm2;

fig = figure('Position', [160 278 750 380]);
subplot(1,2,1)
imagesc(x(x0:end)*1e3,y(y0:end)*1e3,noise_norm(x0:end,y0:end,time_spe))
title({['noise','. t: ',int2str(t(time_spe)),'s'];['Normalized noise']},'fontsize',fontT);
set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
    'ylim', [y(y0)*1e3, y(end)*1e3]);
pbaspect([1 1 1])
xlabel('x (mm)', 'fontsize', fontS)
xticklabels({'0','5','10','15'})
ylabel('y (mm)', 'fontsize', fontS)
colormap(cmap)
colorbar
subplot(1,2,2)
imagesc(x(x0:end)*1e3,y(y0:end)*1e3,BOLD_noise(x0:end,y0:end,time_spe))
title({['BOLD+noise','. t: ',int2str(t(time_spe)),'s'];['Norm noise - Res 0.125 mm']},'fontsize',fontT);
set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
    'ylim', [y(y0)*1e3, y(end)*1e3]);
pbaspect([1 1 1])
xlabel('x (mm)', 'fontsize', fontS)
xticklabels({'0','5','10','15'})
ylabel('y (mm)', 'fontsize', fontS)
colormap(cmap)
colorbar

%% Wiener - fine grid 0.125 mm resolution
params.noise = (snr/200)*(1/0.125); % for Wiener filter

% no noise
Wiener_2D_BOLD = deconvolution_HybridWiener_2D(deconvResponses_Forward_2D.reconvBOLD, ...
                                             x, y, t, params);

figure('Name','Wiener - no noise, noise - Res. 0.125 mm','NumberTitle','on', ...
    'Position', [160 160 1000 760]),
subplot(2,2,1)
imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_BOLD.reconvBOLD(x0:end,y0:end,time_spe))
title({['BOLD Wiener','. t: ',int2str(t(time_spe)),'s'];['Res.: 0.125 mm.']},'fontsize',fontT);
set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
    'ylim', [y(y0)*1e3, y(end)*1e3]);
pbaspect([1 1 1])
xlabel('x (mm)', 'fontsize', fontS)
xticklabels({'0','5','10','15'})
ylabel('y (mm)', 'fontsize', fontS)
colormap(cmap)
colorbar
subplot(2,2,2)
imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_BOLD.neural(x0:end,y0:end,time_NA))
title({['recovered NA','. t: ',int2str(t(time_NA)),'s'];['Res.: 0.125 mm.']},'fontsize',fontT);
set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
    'ylim', [y(y0)*1e3, y(end)*1e3]);
pbaspect([1 1 1])
xlabel('x (mm)', 'fontsize', fontS)
xticklabels({'0','5','10','15'})
ylabel('y (mm)', 'fontsize', fontS)
colormap(cmap)
colorbar

% added noise
Wiener_2D_noisyBOLD = deconvolution_HybridWiener_2D(BOLD_noise, x, y, t, params);

subplot(2,2,3)
imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_noisyBOLD.reconvBOLD(x0:end,y0:end,time_spe))
title({['noise BOLD Wiener','. t: ',int2str(t(time_spe)),'s'];['Res.: 0.125 mm.']},'fontsize',fontT);
set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
    'ylim', [y(y0)*1e3, y(end)*1e3]);
pbaspect([1 1 1])
xlabel('x (mm)', 'fontsize', fontS)
xticklabels({'0','5','10','15'})
ylabel('y (mm)', 'fontsize', fontS)
colormap(cmap)
colorbar
subplot(2,2,4)
imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_noisyBOLD.neural(x0:end,y0:end,time_NA))
title({['recovered NA','. t: ',int2str(t(time_NA)),'s'];['Res.: 0.125 mm.']},'fontsize',fontT);
set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
    'ylim', [y(y0)*1e3, y(end)*1e3]);
pbaspect([1 1 1])
xlabel('x (mm)', 'fontsize', fontS)
xticklabels({'0','5','10','15'})
ylabel('y (mm)', 'fontsize', fontS)
colormap(cmap)
colorbar

%% Modulation Wiener - without noise and with noise - 0.125 mm

% OP 0 or 90
if orientation == 0 || orientation == 90
    
    figure('Name','Wiener BOLD profile x,y direction - Res.: 0.125 mm','NumberTitle','on', ...
        'Position', [160 160 1000 760]),
    
    % Profile Wiener 0.125 mm no noise
    [Wiener_rowx,Wiener_coly] = profile(Wiener_2D_BOLD.reconvBOLD(x0:end,y0:end,time_spe));
    % Profile NA rec 0.125 mm with noise
    [Wiener_rowx_NA,Wiener_coly_NA] = profile(Wiener_2D_BOLD.neural(x0:end,y0:end,time_spe-2));
    % Profile Wiener 0.125 mm
    [Wiener_noise_rowx,Wiener_noise_coly] = profile(Wiener_2D_noisyBOLD.reconvBOLD(x0:end,y0:end,time_spe));
    % Profile NA rec 0.125 mm
    [Wiener_noise_rowx_NA,Wiener_noise_coly_NA] = profile(Wiener_2D_noisyBOLD.neural(x0:end,y0:end,time_spe-2));
    
    % Modulation Wiener 0.125 mm
    [visibility_Wiener,MaxWiener] = modulation_na(Wiener_rowx);
    % Modulation Wiener 0.125 mm
    [visibility_Wiener_noise,MaxWiener_noise] = modulation_na(Wiener_noise_rowx);
    
    % Modulation NA 0.125 mm
    [visibility_Wiener_NA,~] = modulation_na(Wiener_rowx_NA);
    % Modulation NA 0.125 mm
    [visibility_Wiener_noise_NA,MaxWi_NA] = modulation_na(Wiener_noise_rowx_NA);
    
        subplot(2,3,1)
        plot(y(y0:end)*1e3,Wiener_rowx, 'k-','LineWidth',fontL)
        title({['Wiener - Prof x'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximum+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])

        subplot(2,3,2)
        plot(x(x0:end)*1e3,Wiener_coly, 'k-','LineWidth',fontL)
        title({['Wiener - Prof y'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximum+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        subplot(2,3,3)
        plot(y(y0:end)*1e3,Wiener_noise_rowx, 'k-','LineWidth',fontL)
        title({['Wiener noise - Prof x'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximum+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])

        subplot(2,3,4)
        plot(x(x0:end)*1e3,Wiener_noise_coly, 'k-','LineWidth',fontL)
        title({['Wiener noise - Prof y'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximum+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        subplot(2,3,5)
        plot(y(y0:end)*1e3,Wiener_rowx_NA, 'k-','LineWidth',fontL)
        title({['NA rec - Prof x'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA rec', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxWi_NA+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        subplot(2,3,6)
        plot(y(y0:end)*1e3,Wiener_noise_rowx_NA, 'k-','LineWidth',fontL)
        title({['NA rec noise - Prof x'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA rec', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxWi_NA+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
    
end

% OP 45 or 135
if orientation == 45 || orientation == 135
   
    % Diagonal Wiener 0.125 mm
    [~, Wiener_dicr, Wiener_didecr] = Diagonal(Wiener_2D_BOLD.reconvBOLD(x0:end,y0:end,time_spe));
    % Diagonal NA 0.125 mm
    [~, Wiener_dicr_NA, Wiener_didecr_NA] = Diagonal(Wiener_2D_BOLD.neural(x0:end,y0:end,time_spe-2));

    % Diagonal noise Wiener 0.125 mm
    [~, Wiener_noise_dicr, Wiener_noise_didecr] = Diagonal(Wiener_2D_noisyBOLD.reconvBOLD(x0:end,y0:end,time_spe));
    % Diagonal noise NA 0.125 mm
    [~, Wiener_noise_dicr_NA, Wiener_noise_didecr_NA] = Diagonal(Wiener_2D_noisyBOLD.neural(x0:end,y0:end,time_spe-2));
    
    % Modulation Wiener 0.125 mm
    [visibilityd_Wiener,~] = modulation_na(Wiener_didecr);
    % Modulation noise Wiener 0.125 mm
    [visibilityd_Wiener_noise,~] = modulation_na(Wiener_noise_didecr);
    
    % Modulation NA 0.125 mm
    [visibilityd_Wiener_NA,~] = modulation_na(Wiener_didecr_NA);
    % Modulation noise NA 0.125 mm
    [visibilityd_Wiener_noise_NA,MaxWi_NA] = modulation_na(Wiener_noise_didecr_NA);
    

    figure('Name','Wiener BOLD profile diagonal direction - Res.: 0.125 mm','NumberTitle','on', ...
            'Position', [160 160 1000 760]),
    
        subplot(2,3,1)
        plot(x(x0:end)*1e3,Wiener_dicr,'k-','LineWidth',fontL)
        title({['Wiener - Prof diag \'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximumd+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])

        subplot(2,3,2)
        plot(x(x0:end)*1e3,Wiener_didecr,'k-','LineWidth',fontL)
        title({['Wiener - Prof diag /'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximumd+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])

        subplot(2,3,3)
        plot(x(x0:end)*1e3,Wiener_noise_dicr,'k-','LineWidth',fontL)
        title({['Wiener noise - Prof diag \'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximumd+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])

        subplot(2,3,4)
        plot(x(x0:end)*1e3,Wiener_noise_didecr,'k-','LineWidth',fontL)
        title({['Wiener noise - Prof diag /'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (Maximumd+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        subplot(2,3,5)
        plot(x(x0:end)*1e3,Wiener_didecr_NA,'k-','LineWidth',fontL)
        title({['NA rec - Prof diag /'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA rec', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxWi_NA+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        subplot(2,3,6)
        plot(x(x0:end)*1e3,Wiener_noise_didecr_NA,'k-','LineWidth',fontL)
        title({['NA rec noise - Prof diag /'];['Res: 0.125 mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA rec', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxWi_NA+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
    
end


%% Resize

% Actually, smoothing applying rectangular function at k-space, varying kc
% according to the required resolution. Quantity of pixels does not change.

% scale = 0.5;      % Resolution 0.25 mm
% scale = 0.25;     % Resolution 0.5 mm
% scale = 0.1667;   % Resolution 0.75 mm
% scale = 0.125;    % Resolution 1.0 mm
% scale = 0.0833;   % Resolution 1.5 mm

res = {'0.25', '0.5', '0.75', '1.0'};
res1 = {0.25, 0.5, 0.75, 1.0};

%% resize (Smoothing at kc) no noise - BOLD. 
% Quantity of pixels does not change. - Smoothing
% Quantity of pixels change - Smoothing + pixelation
% BOLD

visibility_BOLD_res_all = zeros(1,length(res));
visibility_BOLD_res_pi_all = zeros(1,length(res));
visibilityd_BOLD_res_all = zeros(1,length(res));
visibilityd_BOLD_res_pi_all = zeros(1,length(res));
f8 = figure('Name','Resize - smoothing at kc. BOLD','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
f9 = figure('Name','Resize BOLD profile','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
f10 = figure('Name','Resize - smoothing at kc + pixelation. BOLD','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
f7 = figure('Name','Resize BOLD pixelated profile','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
f2 = figure('Name','Resize NA pixelated profile','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
jj = 1;
for scale = [0.5, 0.25, 0.1667, 0.125]

BOLD_2D_res = deconvResponses_Forward_2D.reconvBOLD;
NA_2D_res = deconvResponses_Forward_2D.neural;

% Smoothing at k-space rect function
[BOLD_res] = Smoothing_sinc(BOLD_2D_res,x,y,t,scale);
[NA_res] = Smoothing_sinc(NA_2D_res,x,y,t,scale);

% Pixelation - bicubic
[BOLD_res_pi] = imresize(BOLD_res,scale);
[NA_res_pi] = imresize(NA_res,scale);

[row,col] = size(BOLD_res_pi(:,:,1));
BOLD_res_pi1 = BOLD_res_pi((row/2):end-1,(col/2):end-1,:);

[rowa,cola] = size(NA_res_pi(:,:,1));
NA_res_pi1 = NA_res_pi((rowa/2):end-1,(cola/2):end-1,:);


    figure(f8);
    subplot(2,4,jj*2-1)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,NA_res(x0:end,y0:end,time_NA))
    title({['NA','. t: ',int2str(t(time_NA)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    subplot(2,4,jj*2)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,BOLD_res(x0:end,y0:end,time_spe))
    title({['BOLD','. t: ',int2str(t(time_spe)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    
    figure(f10)
	subplot(2,4,jj*2-1)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,NA_res_pi1(:,:,time_NA))
    title({['NA','. t: ',int2str(t(time_NA)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    subplot(2,4,jj*2)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,BOLD_res_pi1(:,:,time_spe))
    title({['BOLD pi','. t: ',int2str(t(time_spe)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar

% Modulation BOLD res. no noise
% OP 0 or 90
if orientation == 0 || orientation == 90
    % Profile BOLD res
    [BOLD_res_rowx,BOLD_res_coly] = profile(BOLD_res(x0:end,y0:end,time_spe));
    % Profile BOLD res pixelated
    [BOLD_res_pi_rowx,BOLD_res_pi_coly] = profile(BOLD_res_pi1(:,:,time_spe));
    x_res = linspace(x(x0),x(end),length(BOLD_res_pi_rowx));
    y_res = linspace(y(y0),y(end),length(BOLD_res_pi_coly));
    
    % Profile NA res pixelated
    [NA_res_pi_rowx,NA_res_pi_coly] = profile(NA_res_pi1(:,:,time_spe-2));
    x_res_NA = linspace(x(x0),x(end),length(NA_res_pi_rowx));
    y_res_NA = linspace(y(y0),y(end),length(NA_res_pi_coly));
    
    % Modulation BOLD res
    [visibility_BOLD_res,MaxBOLD_res] = modulation_na(BOLD_res_rowx);
    visibility_BOLD_res_all(jj) = visibility_BOLD_res;
    MaxBOLD_res_all(jj) = MaxBOLD_res;
    % Modulation BOLD res pixelated
    [visibility_BOLD_res_pi,~] = modulation_na(BOLD_res_pi_rowx);
    visibility_BOLD_res_pi_all(jj) = visibility_BOLD_res_pi;
    
    % Modulation NA res pixelated
    [visibility_NA_res_pi,MaxNA_res] = modulation_na(NA_res_pi_rowx);
    visibility_NA_res_pi_all(jj) = visibility_NA_res_pi;
    
        figure(f9);
        subplot(2,4,jj*2-1)
        plot(y(y0:end)*1e3,BOLD_res_rowx, 'k-','LineWidth',fontL)
        title({['BOLD - Prof x'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxBOLD_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x(x0:end)*1e3,BOLD_res_coly, 'k-','LineWidth',fontL)
        title({['BOLD - Prof y'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxBOLD_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f7);
        subplot(2,4,jj*2-1)
        plot(y_res*1e3,BOLD_res_pi_rowx, 'k-','LineWidth',fontL)
        title({['BOLD pi - Prof x'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxBOLD_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res*1e3,BOLD_res_pi_coly, 'k-','LineWidth',fontL)
        title({['BOLD pi - Prof y'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxBOLD_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f2);
        subplot(2,4,jj*2-1)
        plot(y_res_NA*1e3,NA_res_pi_rowx, 'k-','LineWidth',fontL)
        title({['NA pi - Prof x'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res_NA*1e3,NA_res_pi_coly, 'k-','LineWidth',fontL)
        title({['NA pi - Prof y'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
   
end

if orientation == 45 || orientation == 135
    % Diagonal BOLD res
    [~, BOLD_res_dicr, BOLD_res_didecr] = Diagonal(BOLD_res(x0:end,y0:end,time_spe));
    % Diagonal BOLD res pixelated
    [~, BOLD_res_pi_dicr, BOLD_res_pi_didecr] = Diagonal(BOLD_res_pi1(:,:,time_spe));
    x_res = linspace(x(x0),x(end),length(BOLD_res_pi_dicr));
    y_res = linspace(y(y0),y(end),length(BOLD_res_pi_didecr));
    
    % Diagonal NA res pixelated
    [~, NA_res_pi_dicr, NA_res_pi_didecr] = Diagonal(NA_res_pi1(:,:,time_spe-2));
    x_res_NA = linspace(x(x0),x(end),length(NA_res_pi_dicr));
    y_res_NA = linspace(y(y0),y(end),length(NA_res_pi_didecr));
    
    % Modulation BOLD res
    [visibilityd_BOLD_res,MaxBOLD_resd] = modulation_na(BOLD_res_didecr);
    visibilityd_BOLD_res_all(jj) = visibilityd_BOLD_res;
    % Modulation BOLD res pixelated
    [visibilityd_BOLD_res_pi,~] = modulation_na(BOLD_res_pi_didecr);
    visibilityd_BOLD_res_pi_all(jj) = visibilityd_BOLD_res_pi;
    
    % Modulation NA res pixelated
    [visibilityd_NA_res_pi,MaxNA_res] = modulation_na(NA_res_pi_didecr);
    visibilityd_NA_res_pi_all(jj) = visibilityd_NA_res_pi;
    
        figure(f9)
        subplot(2,4,jj*2-1)
        plot(x(x0:end)*1e3,BOLD_res_dicr,'k-','LineWidth',fontL)
        title({['BOLD - Prof diag \'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxBOLD_resd+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x(x0:end)*1e3,BOLD_res_didecr,'k-','LineWidth',fontL)
        title({['BOLD - Prof diag /'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxBOLD_resd+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])

        figure(f7)
        subplot(2,4,jj*2-1)
        plot(y_res*1e3,BOLD_res_pi_dicr,'k-','LineWidth',fontL)
        title({['BOLD pi - Prof diag \'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxBOLD_resd+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res*1e3,BOLD_res_pi_didecr,'k-','LineWidth',fontL)
        title({['BOLD pi - Prof diag /'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxBOLD_resd+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f2)
        subplot(2,4,jj*2-1)
        plot(y_res_NA*1e3,NA_res_pi_dicr,'k-','LineWidth',fontL)
        title({['NA pi - Prof diag \'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res_NA*1e3,NA_res_pi_didecr,'k-','LineWidth',fontL)
        title({['NA pi - Prof diag /'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
end

jj = jj+1;
end


%% Resize - smoothing at kc - BOLD + noise. 
% Quantity of pixels does not change. 
% BOLD + noise
jj = 1;
figure('Name','Resize - smoothing at kc. BOLD + noise','NumberTitle','on', ...
    'Position', [160 160 1000 760]),
for scale = [0.5, 0.25, 0.1667, 0.125]
subplot(2,2,jj)

noise_norm2_res = noise_norm*(1/res1{jj}); 
BOLD_noise_res1 =  deconvResponses_Forward_2D.reconvBOLD + noise_norm2_res;

BOLD_noise_2D_res = BOLD_noise_res1;

[BOLD_noise_res] = Smoothing_sinc(BOLD_noise_2D_res,x,y,t,scale);
BOLD_noise_res_rms(jj) = rms(BOLD_noise_res(:));

imagesc(x(x0:end)*1e3,y(y0:end)*1e3,BOLD_noise_res(x0:end,y0:end,time_spe))
title({['BOLD+noise','. t: ',int2str(t(time_spe)),'s'];['Smoothing for Res.: ', res{jj},' mm']},'fontsize',fontT);
set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
    'ylim', [y(y0)*1e3, y(end)*1e3]);
pbaspect([1 1 1])
xlabel('x (mm)', 'fontsize', fontS)
xticklabels({'0','5','10','15'})
ylabel('y (mm)', 'fontsize', fontS)
colormap(cmap)
colorbar

jj = jj+1;
end


%% Resize - smoothing at kc - BOLD Wiener - NO noise

% resize (Smoothing at kc) NO noise added.
% Not changing quantity of pixels.
% BOLD Wiener

visibility_Wiener_res_all = zeros(1,length(res));
visibility_Wiener_res_pi_all = zeros(1,length(res));
visibilityd_Wiener_res_all = zeros(1,length(res));
visibilityd_Wiener_res_pi_all = zeros(1,length(res));
f11 = figure('Name','Resize - smoothing at kc. BOLD Wiener','NumberTitle','on', ...
    'Position', [160 160 1500 760]);
f12 = figure('Name','Resize Wiener profile','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
f5 = figure('Name','Resize - smoothing at kc + pixelation. BOLD Wiener','NumberTitle','on', ...
    'Position', [160 160 1500 760]);
f6 = figure('Name','Resize Wiener pixelated profile','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
f1 = figure('Name','Resize NA rec (Wiener no noise) pixelated profile','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
jj = 1;
for scale = [0.5, 0.25, 0.1667, 0.125]

params.noise = (snr/200)*(1/res1{jj}); % for the Wiener filter
Wiener_2D_BOLD = deconvolution_HybridWiener_2D(deconvResponses_Forward_2D.reconvBOLD, ...
                                             x, y, t, params);

Wiener_2D_BOLD_BOLD_res = Wiener_2D_BOLD.reconvBOLD;
Wiener_2D_BOLD_NA_res = Wiener_2D_BOLD.neural;

[Wiener_2D_BOLD_res,~] = Smoothing_sinc(Wiener_2D_BOLD_BOLD_res,x,y,t,scale);
[Wiener_2D_NA_res,~] = Smoothing_sinc(Wiener_2D_BOLD_NA_res,x,y,t,scale);

% Pixelation - bicubic
[Wiener_2D_BOLD_res_pi] = imresize(Wiener_2D_BOLD_res,scale);
[Wiener_2D_NA_res_pi] = imresize(Wiener_2D_NA_res,scale);

[row,col] = size(Wiener_2D_BOLD_res_pi(:,:,1));
Wiener_2D_BOLD_res_pi1 = Wiener_2D_BOLD_res_pi((row/2):end-1,(col/2):end-1,:);

[rowa,cola] = size(Wiener_2D_NA_res_pi(:,:,1));
Wiener_2D_NA_res_pi1 = Wiener_2D_NA_res_pi((rowa/2):end-1,(cola/2):end-1,:);

    figure(f11)
    subplot(2,4,jj*2-1)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_BOLD_res(x0:end,y0:end,time_spe))
    title({['BOLD Wiener','. t: ',int2str(t(time_spe)),'s'];['Sth at Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    subplot(2,4,jj*2)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_NA_res(x0:end,y0:end,time_NA))
    title({['recovered NA','. t: ',int2str(t(time_NA)),'s'];['Sth at Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    
    figure(f5)
	subplot(2,4,jj*2-1)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_BOLD_res_pi1(:,:,time_spe))
    title({['BOLD Wiener pi','. t: ',int2str(t(time_spe)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    subplot(2,4,jj*2)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_NA_res_pi1(:,:,time_NA))
    title({['recovered NA pi','. t: ',int2str(t(time_NA)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    
    
% Modulation Wiener BOLD res. no noise
% OP 0 or 90
if orientation == 0 || orientation == 90
    % Profile Wiener res
    [Wiener_res_rowx,Wiener_res_coly] = profile(Wiener_2D_BOLD_res(x0:end,y0:end,time_spe));
    % Profile Wiener res pixelated
    [Wiener_res_pi_rowx,Wiener_res_pi_coly] = profile(Wiener_2D_BOLD_res_pi1(:,:,time_spe));
    x_res = linspace(x(x0),x(end),length(Wiener_res_pi_rowx));
    y_res = linspace(y(y0),y(end),length(Wiener_res_pi_coly));
    % Profile NA recovered res pixelated
    [Wiener_NA_res_pi_rowx,Wiener_NA_res_pi_coly] = profile(Wiener_2D_NA_res_pi1(:,:,time_spe-2));
    x_res_NA = linspace(x(x0),x(end),length(Wiener_NA_res_pi_rowx));
    y_res_NA = linspace(y(y0),y(end),length(Wiener_NA_res_pi_coly));
    
    % Modulation Wiener res
    [visibility_Wiener_res,MaxWiener_res] = modulation_na(Wiener_res_rowx);
    visibility_Wiener_res_all(jj) = visibility_Wiener_res;
    MaxWiener_res_all(jj) = MaxWiener_res;
    % Modulation Wiener res
    [visibility_Wiener_res_pi,~] = modulation_na(Wiener_res_pi_rowx);
    visibility_Wiener_res_pi_all(jj) = visibility_Wiener_res_pi;
    % Modulation NA recovered res
    [visibility_Wiener_NA_res_pi,MaxNA_res] = modulation_na(Wiener_NA_res_pi_rowx);
    visibility_Wiener_NA_res_pi_all(jj) = visibility_Wiener_NA_res_pi;
    
        figure(f12);
        subplot(2,4,jj*2-1)
        plot(y(y0:end)*1e3,Wiener_res_rowx, 'k-','LineWidth',fontL)
        title({['Wiener - Prof x'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x(x0:end)*1e3,Wiener_res_coly, 'k-','LineWidth',fontL)
        title({['Wiener - Prof y'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f6);
        subplot(2,4,jj*2-1)
        plot(y_res*1e3,Wiener_res_pi_rowx, 'k-','LineWidth',fontL)
        title({['Wiener pi - Prof x'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res*1e3,Wiener_res_pi_coly, 'k-','LineWidth',fontL)
        title({['Wiener pi - Prof y'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f1);
        subplot(2,4,jj*2-1)
        plot(y_res_NA*1e3,Wiener_NA_res_pi_rowx, 'k-','LineWidth',fontL)
        title({['NA rec pi - Prof x'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA rec', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res_NA*1e3,Wiener_NA_res_pi_coly, 'k-','LineWidth',fontL)
        title({['NA rec pi - Prof y'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA rec', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
   
    
end

if orientation == 45 || orientation == 135
    % Diagonal Wiener res
    [~, Wiener_res_dicr, Wiener_res_didecr] = Diagonal(Wiener_2D_BOLD_res(x0:end,y0:end,time_spe));
    % Diagonal Wiener res pixelated
    [~, Wiener_res_pi_dicr, Wiener_res_pi_didecr] = Diagonal(Wiener_2D_BOLD_res_pi1(:,:,time_spe));
    x_res = linspace(x(x0),x(end),length(Wiener_res_pi_dicr));
    y_res = linspace(y(y0),y(end),length(Wiener_res_pi_didecr));
    % Diagonal NA rec res pixelated
    [~, Wiener_NA_res_pi_dicr, Wiener_NA_res_pi_didecr] = Diagonal(Wiener_2D_NA_res_pi1(:,:,time_spe-2));
    x_res_NA = linspace(x(x0),x(end),length(Wiener_NA_res_pi_dicr));
    y_res_NA = linspace(y(y0),y(end),length(Wiener_NA_res_pi_didecr));
    
    % Modulation Wiener 0.125 mm
    [visibilityd_Wiener_res,MaxWiener_res] = modulation_na(Wiener_res_didecr);
    visibilityd_Wiener_res_all(jj) = visibilityd_Wiener_res;
    % Modulation Wiener 0.125 mm pixelated
    [visibilityd_Wiener_res_pi,~] = modulation_na(Wiener_res_pi_didecr);
    visibilityd_Wiener_res_pi_all(jj) = visibilityd_Wiener_res_pi;
    % Modulation NA rec 0.125 mm pixelated
    [visibilityd_Wiener_NA_res_pi,MaxNA_res] = modulation_na(Wiener_NA_res_pi_didecr);
    visibilityd_Wiener_NA_res_pi_all(jj) = visibilityd_Wiener_NA_res_pi;
    
        figure(f12)
        subplot(2,4,jj*2-1)
        plot(x(x0:end)*1e3,Wiener_res_dicr,'k-','LineWidth',fontL)
        title({['Wiener - Prof diag \'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x(x0:end)*1e3,Wiener_res_didecr,'k-','LineWidth',fontL)
        title({['Wiener - Prof diag /'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f6)
        subplot(2,4,jj*2-1)
        plot(y_res*1e3,Wiener_res_pi_dicr,'k-','LineWidth',fontL)
        title({['Wiener pi - Prof diag \'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res*1e3,Wiener_res_pi_didecr,'k-','LineWidth',fontL)
        title({['Wiener pi - Prof diag /'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f1)
        subplot(2,4,jj*2-1)
        plot(y_res*1e3,Wiener_NA_res_pi_dicr,'k-','LineWidth',fontL)
        title({['NA rec pi - Prof diag \'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA rec', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res*1e3,Wiener_NA_res_pi_didecr,'k-','LineWidth',fontL)
        title({['NA rec pi - Prof diag /'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA rec', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])

end

jj = jj+1;
end


%% Resize - smoothing at kc - noise BOLD Wiener - Noise

% resize (Smoothing at kc) noise added.
% Not changing quantity of pixels.
% Noise BOLD Wiener

visibilityd_Wiener_noise_res_all = zeros(1,length(res));
visibilityd_Wiener_noise_res_pi_all = zeros(1,length(res));
visibility_Wiener_noise_res_all = zeros(1,length(res));
visibility_Wiener_noise_res_pi_all = zeros(1,length(res));
f13 = figure('Name','Resize - smoothing at kc. Noise BOLD Wiener','NumberTitle','on', ...
    'Position', [160 160 1500 760]);
f14 = figure('Name','Resize Noise BOLD Wiener profile','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
f3 = figure('Name','Resize - smoothing at kc + pixelation. Noise BOLD Wiener','NumberTitle','on', ...
    'Position', [160 160 1500 760]);
f4 = figure('Name','Resize Noise BOLD Wiener pixelated profile','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
f15 = figure('Name','Resize Noise NA rec Wiener pixelated profile','NumberTitle','on', ...
        'Position', [160 160 1500 760]);
jj = 1;
for scale = [0.5, 0.25, 0.1667, 0.125]

noise_norm2_res = noise_norm*(1/res1{jj}); 
BOLD_noise_res =  deconvResponses_Forward_2D.reconvBOLD + noise_norm2_res;
params.noise = (snr/200)*(1/res1{jj}); % for Wiener filter
Wiener_2D_noisyBOLD_res1 = deconvolution_HybridWiener_2D(BOLD_noise_res, x, y, t, params); 

Wiener_2D_noisyBOLD_BOLD_res = Wiener_2D_noisyBOLD_res1.reconvBOLD;
Wiener_2D_noisyBOLD_NA_res = Wiener_2D_noisyBOLD_res1.neural;

[Wiener_2D_noisyBOLD_res,vWidth] = Smoothing_sinc(Wiener_2D_noisyBOLD_BOLD_res,x,y,t,scale);
[Wiener_2D_noisyNA_res,vWidth1] = Smoothing_sinc(Wiener_2D_noisyBOLD_NA_res,x,y,t,scale);

% Pixelation - bicubic
[Wiener_2D_noisyBOLD_res_pi] = imresize(Wiener_2D_noisyBOLD_res,scale);
[Wiener_2D_noisyNA_res_pi] = imresize(Wiener_2D_noisyNA_res,scale);

[row,col] = size(Wiener_2D_noisyBOLD_res_pi(:,:,1));
Wiener_2D_noisyBOLD_res_pi1 = Wiener_2D_noisyBOLD_res_pi((row/2):end-1,(col/2):end-1,:);

[rowa,cola] = size(Wiener_2D_noisyNA_res_pi(:,:,1));
Wiener_2D_noisyNA_res_pi1 = Wiener_2D_noisyNA_res_pi((rowa/2):end-1,(cola/2):end-1,:);

    figure(f13)
    subplot(2,4,jj*2-1)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_noisyBOLD_res(x0:end,y0:end,time_spe))
    title({['noise BOLD Wiener','. t: ',int2str(t(time_spe)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    subplot(2,4,jj*2)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_noisyNA_res(x0:end,y0:end,time_NA))
    title({['recovered NA','. t: ',int2str(t(time_NA)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    
    figure(f3)
	subplot(2,4,jj*2-1)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_noisyBOLD_res_pi1(:,:,time_spe))
    title({['noise BOLD Wiener pi','. t: ',int2str(t(time_spe)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar
    subplot(2,4,jj*2)
    imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Wiener_2D_noisyNA_res_pi1(:,:,time_NA))
    title({['recovered NA pi','. t: ',int2str(t(time_NA)),'s'];['Res.: ', res{jj},' mm']},'fontsize',fontT);
    set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
    pbaspect([1 1 1])
    xlabel('x (mm)', 'fontsize', fontS)
    xticklabels({'0','5','10','15'})
    ylabel('y (mm)', 'fontsize', fontS)
    colormap(cmap)
    colorbar

% Modulation Noise Wiener BOLD res. noise
% OP 0 or 90
if orientation == 0 || orientation == 90
    % Profile Wiener res
    [Wiener_noise_res_rowx,Wiener_noise_res_coly] = profile(Wiener_2D_noisyBOLD_res(x0:end,y0:end,time_spe));
    % Profile Wiener res pixelated
    [Wiener_noise_res_pi_rowx,Wiener_noise_res_pi_coly] = profile(Wiener_2D_noisyBOLD_res_pi1(:,:,time_spe));
    x_res = linspace(x(x0),x(end),length(Wiener_noise_res_pi_rowx));
    y_res = linspace(y(y0),y(end),length(Wiener_noise_res_pi_coly));
    % Profile NA rec res pixelated
    [Wiener_NA_noise_res_pi_rowx,Wiener_NA_noise_res_pi_coly] = profile(Wiener_2D_noisyNA_res_pi1(:,:,time_spe-2));
    x_res_NA = linspace(x(x0),x(end),length(Wiener_NA_noise_res_pi_rowx));
    y_res_NA = linspace(y(y0),y(end),length(Wiener_NA_noise_res_pi_coly));
    
    % Modulation Wiener res
    [visibility_Wiener_noise_res,MaxWiener_noise_res] = modulation_na(Wiener_noise_res_rowx);
    visibility_Wiener_noise_res_all(jj) = visibility_Wiener_noise_res;
    MaxWiener_noise_res_all(jj) = MaxWiener_noise_res;
    % Modulation Wiener res
    [visibility_Wiener_noise_res_pi,~] = modulation_na(Wiener_noise_res_pi_rowx);
    visibility_Wiener_noise_res_pi_all(jj) = visibility_Wiener_noise_res_pi;
    % Modulation Wiener res
    [visibility_Wiener_NA_noise_res_pi,MaxNA_res] = modulation_na(Wiener_NA_noise_res_pi_rowx);
    visibility_Wiener_NA_noise_res_pi_all(jj) = visibility_Wiener_NA_noise_res_pi;
    
        figure(f14);
        subplot(2,4,jj*2-1)
        plot(y(y0:end)*1e3,Wiener_noise_res_rowx, 'k-','LineWidth',fontL)
        title({['N-Wiener - Prof x'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_noise_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x(x0:end)*1e3,Wiener_noise_res_coly, 'k-','LineWidth',fontL)
        title({['N-Wiener - Prof y'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_noise_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f4);
        subplot(2,4,jj*2-1)
        plot(y_res*1e3,Wiener_noise_res_pi_rowx, 'k-','LineWidth',fontL)
        title({['N-Wiener pi - Prof x'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_noise_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res*1e3,Wiener_noise_res_pi_coly, 'k-','LineWidth',fontL)
        title({['N-Wiener pi - Prof y'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_noise_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f15);
        subplot(2,4,jj*2-1)
        plot(y_res_NA*1e3,Wiener_NA_noise_res_pi_rowx, 'k-','LineWidth',fontL)
        title({['N-NA rec pi - Prof x'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res_NA*1e3,Wiener_NA_noise_res_pi_coly, 'k-','LineWidth',fontL)
        title({['N-NA rec pi - Prof y'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
       
end

if orientation == 45 || orientation == 135
    % Diagonal Wiener pixelated
    [~, Wiener_noise_res_dicr, Wiener_noise_res_didecr] = Diagonal(Wiener_2D_noisyBOLD_res(x0:end,y0:end,time_spe));
    % Diagonal Wiener res pixelated
    [~, Wiener_noise_res_pi_dicr, Wiener_noise_res_pi_didecr] = Diagonal(Wiener_2D_noisyBOLD_res_pi1(:,:,time_spe));
    x_res = linspace(x(x0),x(end),length(Wiener_noise_res_pi_dicr));
    y_res = linspace(y(y0),y(end),length(Wiener_noise_res_pi_didecr));
    % Diagonal NA rec res pixelated
    [~, Wiener_NA_noise_res_pi_dicr, Wiener_NA_noise_res_pi_didecr] = Diagonal(Wiener_2D_noisyNA_res_pi1(:,:,time_spe-2));
    x_res_NA = linspace(x(x0),x(end),length(Wiener_NA_noise_res_pi_dicr));
    y_res_NA = linspace(y(y0),y(end),length(Wiener_NA_noise_res_pi_didecr));
    
    % Modulation Wiener pixelated
    [visibilityd_Wiener_noise_res,MaxWiener_noise_res] = modulation_na(Wiener_noise_res_didecr);
    visibilityd_Wiener_noise_res_all(jj) = visibilityd_Wiener_noise_res;
    % Modulation Wiener pixelated
    [visibilityd_Wiener_noise_res_pi,~] = modulation_na(Wiener_noise_res_pi_didecr);
    visibilityd_Wiener_noise_res_pi_all(jj) = visibilityd_Wiener_noise_res_pi;
    % Modulation NA rec pixelated
    [visibilityd_Wiener_NA_noise_res_pi,MaxNA_res] = modulation_na(Wiener_NA_noise_res_pi_didecr);
    visibilityd_Wiener_NA_noise_res_pi_all(jj) = visibilityd_Wiener_NA_noise_res_pi;
    
        figure(f14)
        subplot(2,4,jj*2-1)
        plot(x(x0:end)*1e3,Wiener_noise_res_dicr,'k-','LineWidth',fontL)
        title({['N-Wiener - Prof diag \'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_noise_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x(x0:end)*1e3,Wiener_noise_res_didecr,'k-','LineWidth',fontL)
        title({['N-Wiener - Prof diag /'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_noise_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f4)
        subplot(2,4,jj*2-1)
        plot(y_res*1e3,Wiener_noise_res_pi_dicr,'k-','LineWidth',fontL)
        title({['N-Wiener pi - Prof diag \'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_noise_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res*1e3,Wiener_noise_res_pi_didecr,'k-','LineWidth',fontL)
        title({['N-Wiener pi - Prof diag /'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 (MaxWiener_noise_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        
        figure(f15)
        subplot(2,4,jj*2-1)
        plot(y_res_NA*1e3,Wiener_NA_noise_res_pi_dicr,'k-','LineWidth',fontL)
        title({['N-NA rec pi - Prof diag \'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])
        subplot(2,4,jj*2)
        plot(x_res_NA*1e3,Wiener_NA_noise_res_pi_didecr,'k-','LineWidth',fontL)
        title({['N-NA rec pi - Prof diag /'];['Res: ', res{jj},' mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        xticklabels({'0','5','10','15'})
        ylabel('NA', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [-0.1 (MaxNA_res+0.05)],'TickLength',[tick tick],'LineWidth',line);
        pbaspect([1 1 1])


end

jj = jj+1;
end

%% Visibility

% Values around 100 means no visibility

if orientation == 0 || orientation == 90
    % visibility: BOLD, Wiener, noise Wiener - at res 0.125 mm
    visibility = [visibility_BOLD, visibility_Wiener, visibility_Wiener_noise]*100;
    % visibility: res.: 0.25 0.50 0.75 1.0 (mm)
    %             BOLD
    %             Wiener
    %             noise Wiener
    % Smoothing only - no pixelation
    visibility_res_all = [visibility_BOLD_res_all, 
                          visibility_Wiener_res_all,
                          visibility_Wiener_noise_res_all]*100;
    % visibility: res.: 0.125 0.25 0.50 0.75 1.0 (mm)
    %             BOLD
    %             Wiener
    %             noise Wiener
    % Smoothing only - no pixelation
    visibility_all = cat(2,visibility',visibility_res_all);
    
    % Smoothing + bicubic - pixelated
    visibility_res_pi_all = [visibility_BOLD_res_pi_all, 
                          visibility_Wiener_res_pi_all,
                          visibility_Wiener_noise_res_pi_all]*100;
    visibility_pi_all = cat(2,visibility',visibility_res_pi_all);
    
    visibility_0125 = [visibility_NA,visibility_Wiener_NA, visibility_Wiener_noise_NA]*100;
    visibility_NA_res = [visibility_NA_res_pi_all,
                        visibility_Wiener_NA_res_pi_all,
                        visibility_Wiener_NA_noise_res_pi_all]*100;
    visibility_NA_all = cat(2,visibility_0125',visibility_NA_res)
                         
end

if orientation == 45 || orientation == 135
    % visibility: BOLD, Wiener, noise Wiener - at res 0.125 mm
    visibilityd = [visibilityd_BOLD, visibilityd_Wiener, visibilityd_Wiener_noise]*100;
    % visibility: res.: 0.25 0.50 0.75 1.0 (mm)
    %             BOLD
    %             Wiener
    %             noise Wiener
    % Smoothing only - no pixelation
    visibilityd_res_all = [visibilityd_BOLD_res_all, 
                           visibilityd_Wiener_res_all,
                           visibilityd_Wiener_noise_res_all]*100;
    % visibility: res.: 0.125 0.25 0.50 0.75 1.0 (mm)
    %             BOLD
    %             Wiener
    %             noise Wiener
    visibilityd_all = cat(2,visibilityd',visibilityd_res_all);
    
    % Smoothing + bicubic - pixelated
    visibilityd_res_pi_all = [visibilityd_BOLD_res_pi_all, 
                           visibilityd_Wiener_res_pi_all,
                           visibilityd_Wiener_noise_res_pi_all]*100;
    visibilityd_pi_all = cat(2,visibilityd',visibilityd_res_pi_all);
    
    visibilityd_0125 = [visibilityd_NA,visibilityd_Wiener_NA, visibilityd_Wiener_noise_NA]*100;
    visibilityd_NA_res = [visibilityd_NA_res_pi_all,
                        visibilityd_Wiener_NA_res_pi_all,
                        visibilityd_Wiener_NA_noise_res_pi_all]*100;
    visibilityd_NA = cat(2,visibilityd_0125',visibilityd_NA_res)
end
