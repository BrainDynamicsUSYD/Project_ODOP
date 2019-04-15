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
y_axisy = 0.22;
y_axisx = 0.22;

cmap = jet;

[xx, yy, tt] = meshgrid(x_experiment, y_experiment, t_experiment);

% Pinwheel angle and positions
TH = [ 0,  45,  90,  135]; 
XC = [7.8, 7.5, 7.2, 7.5];
YC = [7.5, 7.8, 7.5, 7.2];

% stimulation for one stimulus only
orientation = [TH(2)];
positionx = [XC(2)]; 
positiony = [YC(2)];

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
params.v_b = 2e-3; % propagation speed

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
[deconvResponses_Wiener_2D,BOLD_noisy,deconvResponses_Wiener_2D_noisyBOLD] = conv_deconv(deconvResponses_Forward_2D,x,y,t,params,V); 

% calculating maximum time at BOLD Forw
BOLD_max = deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,:);
[Cmax,Imax] = max(BOLD_max(:));
[I1max,I2max,I3max] = ind2sub(size(BOLD_max),Imax);

time_spe = I3max;
%% Plot neural activity and BOLD

figure;
imagesc(x(x0:end)*1e3, y(y0:end)*1e3, squeeze(deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,time_spe)));
set(gca, 'ydir', 'normal');
title('neural activity Forw')
PlotNames

% Struct neural activity and BOLD
titles_stim = {'neural activity','BOLD Forw','BOLD Wiener','ND Wiener', ...
          'nBOLD BOLD','nBOLD BOLD Wiener','nBOLD ND Wiener'};

g1 = 'ND'; w1 = deconvResponses_Forward_2D.neural(x0:end,y0:end,:);
g2 = 'BOLD'; w2 = deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,:);
g3 = 'Wiener'; w3 = deconvResponses_Wiener_2D.reconvBOLD(x0:end,y0:end,:);
g4 = 'ND_Wiener'; w4 = deconvResponses_Wiener_2D.neural(x0:end,y0:end,:);
    
g5 = 'BOLD_noisyB'; w5 = BOLD_noisy(x0:end,y0:end,:);
g6 = 'Wiener_noisyB'; w6 = deconvResponses_Wiener_2D_noisyBOLD.reconvBOLD(x0:end,y0:end,:);
g7 = 'ND_noisyB'; w7 = deconvResponses_Wiener_2D_noisyBOLD.neural(x0:end,y0:end,:);
Stim = struct(g1,w1,g2,w2,g3,w3,g4,w4,g5,w5,g6,w6,g7,w7);
fields_stim = fieldnames(Stim);
    
% Plot neural activity and BOLD
figure('Name','neural activity and BOLD, convoluted and deconvoluted - Res.: 0.125 mm','NumberTitle','on');
for ii=1:numel(fields_stim)
    subplot(2,4,ii)
    if ii == 1
        contourf(x(x0:end)*1e3,y(y0:end)*1e3,resampled_zz(x0:end,y0:end, time_spe))
        title([titles_stim{ii}],'fontsize',fontT);
        set(gca, 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
            'ylim', [y(y0)*1e3, y(end)*1e3]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
        colormap(cmap)
        colorbar
    else
        contourf(x(x0:end)*1e3,y(y0:end)*1e3,Stim.(fields_stim{ii})(:,:,time_spe))
        title([titles_stim{ii}],'fontsize',fontT);
        set(gca, 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
            'ylim', [y(y0)*1e3, y(end)*1e3]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
        colormap(cmap)
        colorbar
    end
end

% Plot neural activity and BOLD
figure('Name','neural activity and BOLD, convoluted and deconvoluted - Res.: 0.125 mm','NumberTitle','on');
for ii=1:numel(fields_stim)
    subplot(2,4,ii)
    if ii == 1
        imagesc(x(x0:end)*1e3,y(y0:end)*1e3,resampled_zz(x0:end,y0:end, time_spe))
        title({[titles_stim{ii},'. t: ',int2str(t(time_spe)),'s'];['Res.: 0.125 mm.']},'fontsize',fontT);
        set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
            'ylim', [y(y0)*1e3, y(end)*1e3]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
        colormap(cmap)
        colorbar
    else
        imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Stim.(fields_stim{ii})(:,:,time_spe))
        title({[titles_stim{ii},'. t: ',int2str(t(time_spe)),'s'];['Res.: 0.125 mm.']},'fontsize',fontT);
        set(gca, 'YDir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
            'ylim', [y(y0)*1e3, y(end)*1e3]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
        colormap(cmap)
        colorbar
    end
end


vectorf1 = zeros(numel(t),numel(fields_stim));
[vector, I1, I2,I3, C] = max_voxel(Stim.(fields_stim{1}));
% Plot neural activity and BOLD max
figure('Name','neural activity and BOLD, convoluted and deconvoluted - Res.: 0.125 mm','NumberTitle','on');
for ii=1:numel(fields_stim)
    
    if ii == numel(fields_stim)
        [vectorf, I1f, I2f,I3, C] = max_voxel(Stim.(fields_stim{ii}));
        vectorf1(:,ii) = vectorf;
        I31(ii) = I3;
        C1(ii) = C;
        subplot(2,4,ii)
        hold on
        imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Stim.(fields_stim{ii})(:,:,I3))
        set(gca, 'ydir', 'normal');
        plot(x(x0+I2-1)*1e3,y(y0+I1-1)*1e3,'k.','MarkerSize',15);
        title({[titles_stim{ii},'. t_{max}: ',int2str(t(I3)),'s'];['Res.: 0.125 mm. Ind.: [',int2str(I1),',',int2str(I2),']']},'fontsize',fontT);
        set(gca, 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
            'ylim', [y(y0)*1e3, y(end)*1e3]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
        colormap(cmap)
        colorbar
        hold off
    else
        [vectorf, I1f, I2f,I3, C] = max_voxel(Stim.(fields_stim{ii}));
        vectorf1(:,ii) = vectorf;
        I31(ii) = I3;
        C1(ii) = C;
        subplot(2,4,ii)
        hold on
        imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Stim.(fields_stim{ii})(:,:,I3))
        set(gca, 'ydir', 'normal');
        plot(x(x0+I2-1)*1e3,y(y0+I1-1)*1e3,'k.','MarkerSize',15);
        title({[titles_stim{ii},'. t_{max}: ',int2str(t(I3)),'s'];['Res.: 0.125 mm. Ind.: [',int2str(I1),',',int2str(I2),']']},'fontsize',fontT);
        set(gca, 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
            'ylim', [y(y0)*1e3, y(end)*1e3]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
        colormap(cmap)
        colorbar
        hold off
    end
    
end

%% Profile of max voxel along time 

figure('Name','Profile of maximum voxel along time - Res.: 0.125 mm','NumberTitle','on')
for ii = 1:numel(fields_stim)
    subplot(2,4,ii)
    plot(t(t0:end),vectorf1(t0:end,ii),'k-','LineWidth',fontL);
    set(gca, 'fontSize', fontS, 'xlim', [0 tmax],'TickLength',[tick tick],'LineWidth',line);
    xlabel('Time (s)', 'fontsize', fontS)
    ylabel([titles_stim{ii}], 'fontsize', fontS)
    title({['Max voxel along t'];['Res.: 0.125 mm']},'fontsize',fontT)
end

%% Mean

% % times that stimulus is on and calculation of average
tton_all1 = ton_all1';
rton_all1 = reshape(tton_all1,[],1);
ton_all_ind = dsearchn(t',rton_all1);

%% Plot BOLD profile x,y 0.125mm

titlesp = {'BOLD Forw','BOLD Wiener','nB BOLD','nB BOLD Wiener'};

wp2 = deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,time_spe);
wp3 = deconvResponses_Wiener_2D.reconvBOLD(x0:end,y0:end,time_spe);
wp5 = BOLD_noisy(x0:end,y0:end,time_spe);
wp6 = deconvResponses_Wiener_2D_noisyBOLD.reconvBOLD(x0:end,y0:end,time_spe);
Stimp = struct(g2,wp2,g3,wp3,g5,wp5,g6,wp6);
fields_stimp = fieldnames(Stimp);

jj=1;
figure('Name','BOLD profile x,y direction - Res.: 0.125 mm','NumberTitle','on'),
for ii = 1:numel(fields_stimp)
    
    if orientation == 0 || orientation == 90
        
        [matrix_rowx,matrix_coly] = profile(Stimp.(fields_stimp{ii}));

        % Modulation of profile
        jja = 1;
        for kka = 2:(length(matrix_rowx)-1)

            % find max's and store in new array
            store_pre_max = matrix_rowx(kka-1);
            store_pos_max = matrix_rowx(kka+1);
            if matrix_rowx(kka) > store_pre_max && matrix_rowx(kka) > store_pos_max
                store_max = matrix_rowx(kka);
                store_max1(ii,jja) = store_max;
                jja=jja+1;
            end
        end

            % find min's and store in new array
        jji = 1;
        for kki=2:(length(matrix_rowx)-1)
            store_pre_min = matrix_rowx(kki-1);
            store_pos_min = matrix_rowx(kki+1);
            if matrix_rowx(kki) < store_pre_min && matrix_rowx(kki) < store_pos_min
                store_min = matrix_rowx(kki);
                store_min1(ii,jji) = store_min;
                jji=jji+1;
            end
        end

        % Max value and indice of each row of stored max
        store_max1t = store_max1';
        [MaxValue,indMaxValue] = max(store_max1t);
        % Getting the max neighour values (at the right) close to the max
        for kkm=1:size(indMaxValue,2)
            store_max1_pos = store_max1(ii,indMaxValue(kkm)+1);
            store_max1_pos1(1,size(indMaxValue,2)) = store_max1_pos;
        end

        % Getting the minimum value close to the max
        kkmi = 1;
        for kkmi=1:size(indMaxValue,2)
            store_min1_min = store_min1(ii,indMaxValue(kkmi));
            store_min1_min1(1,size(indMaxValue,2)) = store_min1_min;
        end


        [Maximum,ind_max]=max(store_max1(ii,:));
        % modulation related to maximum value
        jjma = 1;
        for kkma=2:length(store_max1)
            if kkma <= ind_max && store_max1(ii,kkma) ~= 0 && store_max1(ii,kkma-1) ~= 0
               modul_max = 100-((store_max1(ii,kkma-1)*100)/store_max1(ii,kkma));
               modul1_max(ii,jjma) = modul_max;
               jjma=jjma+1;        
            end
        end

        subplot(2,4,((ii*2)-1))
        plot(y(y0:end)*1e3,matrix_rowx, 'k-','LineWidth',fontL)
        title({['',titlesp{ii},' - Prof x'];['Res: 0.125 mm, a = ',int2str(a),'mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 y_axisy],'TickLength',[tick tick],'LineWidth',line);

        subplot(2,4,(ii*2))
        plot(x(x0:end)*1e3,matrix_coly, 'k-','LineWidth',fontL)
        title({['',titlesp{ii},' - Prof y'];['Res: 0.125 mm, a = ',int2str(a),'mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 y_axisx],'TickLength',[tick tick],'LineWidth',line);

        jj=jj+1;
    end
end

% Visibility:
% ((0.5*(peak1+peak2))-valey1)/((0.5*(peak1+peak2))+valey1);
if orientation == 0 || orientation == 90
    visibility = zeros(1,size(indMaxValue,2));
    for vv=1:size(indMaxValue,2)
        visibility(vv) = ((0.5*(MaxValue(vv)+store_max1_pos1(vv)))-store_min1_min1(vv))/((0.5*(MaxValue(vv)+store_max1_pos1(vv)))+store_min1_min1(vv));
    end
end

%% Plot BOLD profile diagonal res 0.125mm
jj=1;
figure('Name','BOLD profile diagonal direction - Res.: 0.125 mm','NumberTitle','on'),
for ii = 1:numel(fields_stimp)
    
    if orientation == 45 || orientation == 135
    
        [matrix_new, matrix_dicr, matrix_didecr] = Diagonal(Stimp.(fields_stimp{ii}));
        % matrix_dicr = diag(fliplr(matrix)); % diagonal crescent: left top to
        % right bottom  \  --> use for OP 135 degree
        % matrix_didecr = diag(matrix); % diagonal decrescent: left bottom
        % to right top / --> use for OP 45 degree

        % Modulation of profile
        jja = 1;
        for kka = 2:(length(matrix_didecr)-1)

            % find max's and store in new array
            store_pre_maxd = matrix_didecr(kka-1);
            store_pos_maxd = matrix_didecr(kka+1);
            if matrix_didecr(kka) > store_pre_maxd && matrix_didecr(kka) > store_pos_maxd
                store_maxd = matrix_didecr(kka);
                store_max1d(ii,jja) = store_maxd;
                jja=jja+1;
            end
        end

            % find min's and store in new array
        jji = 1;
        for kki=2:(length(matrix_didecr)-1)
            store_pre_mind = matrix_didecr(kki-1);
            store_pos_mind = matrix_didecr(kki+1);
            if matrix_didecr(kki) < store_pre_mind && matrix_didecr(kki) < store_pos_mind
                store_mind = matrix_didecr(kki);
                store_min1d(ii,jji) = store_mind;
                jji=jji+1;
            end
        end

        % Max value and indice of each row of stored max
        store_max1dt = store_max1d';
        [MaxValued,indMaxValued] = max(store_max1dt);
        % Getting the max neighour values (at the right) close to the max
        kkm = 1;
        for kkm=1:size(indMaxValued,2)
            store_max1_posd = store_max1d(ii,indMaxValued(kkm)+1);
            store_max1_pos1d(1,size(indMaxValued,2)) = store_max1_posd;
        end
        % Getting the minimum value close to the max
        kkmi = 1;
        for kkmi=1:size(store_min1d,1)
            store_min1_mind = store_min1d(ii,indMaxValued(kkmi));
            store_min1_min1d(1,size(indMaxValued,2)) = store_min1_mind;
        end

        [Maximumd,ind_maxd]=max(store_max1d(ii,:));
        % modulation related to maximum value
        jjma = 1;
        for kkma=2:length(store_max1d)
            if kkma <= ind_maxd && store_max1d(ii,kkma) ~= 0 && store_max1d(ii,kkma-1) ~= 0
               modul_maxd = 100-((store_max1d(ii,kkma-1)*100)/store_max1d(ii,kkma));
               modul1_maxd(ii,jjma) = modul_maxd;
               jjma=jjma+1;        
            end
        end

        subplot(2,4,((ii*2)-1))
        plot(x(x0:end)*1e3,matrix_dicr,'k-','LineWidth',fontL)
        title({['',titlesp{ii},' - Prof diag \'];['Res: 0.125 mm, a = ',int2str(a), 'mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 y_axisy],'TickLength',[tick tick],'LineWidth',line);

        subplot(2,4,(ii*2))
        plot(x(x0:end)*1e3,matrix_didecr,'k-','LineWidth',fontL)
        title({['',titlesp{ii},' - Prof diag /'];['Res: 0.125 mm, a = ',int2str(a), 'mm']},'fontsize',fontT)
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('BOLD', 'fontsize', fontS)
        set(gca, 'fontSize', fontS, 'xlim', [0 xmax], 'ylim', [0 y_axisx],'TickLength',[tick tick],'LineWidth',line);

        jj=jj+1;
    end
end

% Visibility:
% ((0.5*(peak1+peak2))-valey1)/((0.5*(peak1+peak2))+valey1);
if orientation == 45 || orientation == 135
    visibilityd = zeros(1,size(indMaxValued,2));
    for vv=1:size(indMaxValued,2)
        visibilityd(vv) = ((0.5*(MaxValued(vv)+store_max1_pos1d(vv)))-store_min1_min1d(vv))/((0.5*(MaxValued(vv)+store_max1_pos1d(vv)))+store_min1_min1d(vv));
    end
end

%% resize

titles = {'neural activity resized','BOLD response Forw','BOLD response Wiener','neural activity Wiener', ...
          'noisyBOLD BOLD','noisyBOLD BOLD Wiener','noisyBOLD ND Wiener'};
      
titlesp_res = {'BOLD Forw resized','BOLD Wiener resized', ...
          'nB BOLD resized','nB BOLD Wiener resized'};  

Maxmax = zeros(5,4);
Maxmax_pos = zeros(5,4);
MaxMin_pos = zeros(5,4);
Maxmin = zeros(5,4);
Maxmaxd = zeros(5,4);
Maxmax_posd = zeros(5,4);
MaxMin_posd = zeros(5,4);
deconvResponses_Forward_2D_res = deconvResponses_Forward_2D;
res = {'0.25', '0.5', '0.75', '1.0', '1.5'};
res1 = {0.25, 0.5, 0.75, 1.0, 1.5};
jj = 1;      
for scale = [0.5, 0.25, 0.1667, 0.125, 0.0833]
    
    % calculate new voxel volume
    Vres = res1{jj}^3;
    % calculate BOLD and NA considering new V and noise
    [deconvResponses_Wiener_2D_res,BOLD_noisy_res,deconvResponses_Wiener_2D_noisyBOLD_res] = conv_deconv(deconvResponses_Forward_2D_res,x,y,t,params,Vres);
    
    % Change resolution
    [ND_r,ND_r_prof,rowND_r,ND_r_bary,ND_r_barx,ND_r_diagcr,ND_r_diagdecr,ND_r_mvec,ND_r_C,ND_r_I1,ND_r_I2,ND_r_I3] = ChangeResolutionNoise(deconvResponses_Forward_2D_res.neural,scale);
    [BOLD_res,BOLD_prof,rowB,BOLD_bary,BOLD_barx,BOLD_diagcr,BOLD_diagdecr,BOLD_mvec,BOLD_C,BOLD_I1,BOLD_I2,BOLD_I3] = ChangeResolutionNoise(deconvResponses_Forward_2D_res.reconvBOLD,scale);
    [Wiener_res,Wiener_prof,rowW,Wiener_bary,Wiener_barx,Wiener_diagcr,Wiener_diagdecr,Wiener_mvec,Wiener_C,Wiener_I1,Wiener_I2,Wiener_I3] = ChangeResolutionNoise(deconvResponses_Wiener_2D_res.reconvBOLD,scale);
    [ND_res,ND_prof,rowND,ND_bary,ND_barx,ND_diagcr,ND_diagdecr,ND_mvec,ND_C,ND_I1,ND_I2,ND_I3] = ChangeResolutionNoise(deconvResponses_Wiener_2D_res.neural,scale);

    [BOLD_noisyBOLD_res,BOLD_noisyBOLD_prof,rowB_noisyBOLD,BOLD_noisyBOLD_bary,BOLD_noisyBOLD_barx,BOLD_noisyBOLD_diagcr,BOLD_noisyBOLD_diagdecr,BOLD_noisyBOLD_mvec, ...
        BOLD_noisyBOLD_C,BOLD_noisyBOLD_I1,BOLD_noisyBOLD_I2,BOLD_noisyBOLD_I3] = ChangeResolutionNoise(BOLD_noisy_res,scale);
    [Wiener_noisyBOLD_res,Wiener_noisyBOLD_prof,rowW_noisyBOLD,Wiener_noisyBOLD_bary,Wiener_noisyBOLD_barx,Wiener_noisyBOLD_diagcr,Wiener_noisyBOLD_diagdecr,Wiener_noisyBOLD_mvec, ...
        Wiener_noisyBOLD_C,Wiener_noisyBOLD_I1,Wiener_noisyBOLD_I2,Wiener_noisyBOLD_I3] = ChangeResolutionNoise(deconvResponses_Wiener_2D_noisyBOLD_res.reconvBOLD,scale);
    [ND_noisyBOLD_res,ND_noisyBOLD_prof,rowND_noisyBOLD,ND_noisyBOLD_bary,ND_noisyBOLD_barx,ND_noisyBOLD_diagcr,ND_noisyBOLD_diagdecr,ND_noisyBOLD_mvec, ...
        ND_noisyBOLD_C,ND_noisyBOLD_I1,ND_noisyBOLD_I2,ND_noisyBOLD_I3] = ChangeResolutionNoise(deconvResponses_Wiener_2D_noisyBOLD_res.neural,scale);
    
    % Struct neural activity and BOLD                        % 'on' times (for mean 'on calculation')
    f7 = 'neural_resize';    v7 = ND_r_prof;                m7 = ND_r(:,:,ton_all_ind);                 
    f1 = 'BOLD';             v1 = BOLD_prof;                m1 = BOLD_res(:,:,ton_all_ind);
    f2 = 'Wiener';           v2 = Wiener_prof;              m2 = Wiener_res(:,:,ton_all_ind);
    f3 = 'neural';           v3 = ND_prof;                  m3 = ND_res(:,:,ton_all_ind);
    f4 = 'BOLD_noisyBOLD';   v4 = BOLD_noisyBOLD_prof;      m4 = BOLD_noisyBOLD_res(:,:,ton_all_ind);
    f5 = 'Wiener_noisyBOLD'; v5 = Wiener_noisyBOLD_prof;    m5 = Wiener_noisyBOLD_res(:,:,ton_all_ind);
    f6 = 'neural_noisyBOLD'; v6 = ND_noisyBOLD_prof;        m6 = ND_noisyBOLD_res(:,:,ton_all_ind);
    Resize = struct(f7,v7,f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6);
    Resize1(jj) = Resize;
    fields = fieldnames(Resize);
          
    Resize_mean = struct(f7,m7,f1,m1,f2,m2,f3,m3,f4,m4,f5,m5,f6,m6);
    Resize1_mean(jj) = Resize_mean;
    fields_mean = fieldnames(Resize_mean);
    
    % struct of max BOLD (value get from res 0.125mm)
    t7 = ND_r(:,:,I3max);
    t1 = BOLD_res(:,:,I3max);
    t2 = Wiener_res(:,:,I3max);
    t3 = ND_res(:,:,I3max);
    t4 = BOLD_noisyBOLD_res(:,:,I3max);
    t5 = Wiener_noisyBOLD_res(:,:,I3max);
    t6 = ND_noisyBOLD_res(:,:,I3max);
    Resizet = struct(f7,t7,f1,t1,f2,t2,f3,t3,f4,t4,f5,t5,f6,t6);
    fieldst = fieldnames(Resizet);
    
    % Struct coordinates max point and max time
    vi17 = ND_r_I1;             vi27 = ND_r_I2;             vi37 = ND_r_I3;             vc7 = ND_r_C;
    vi11 = BOLD_I1;             vi21 = BOLD_I2;             vi31 = BOLD_I3;             vc1 = BOLD_C;
    vi12 = Wiener_I1;           vi22 = Wiener_I2;           vi32 = Wiener_I3;           vc2 = Wiener_C;
    vi13 = ND_I1;               vi23 = ND_I2;               vi33 = ND_I3;               vc3 = ND_C;
    vi14 = BOLD_noisyBOLD_I1;   vi24 = BOLD_noisyBOLD_I2;   vi34 = BOLD_noisyBOLD_I3;   vc4 = BOLD_noisyBOLD_C;
    vi15 = Wiener_noisyBOLD_I1; vi25 = Wiener_noisyBOLD_I2; vi35 = Wiener_noisyBOLD_I3; vc5 = Wiener_noisyBOLD_C;
    vi16 = ND_noisyBOLD_I1;     vi26 = ND_noisyBOLD_I2;     vi36 = ND_noisyBOLD_I3;     vc6 = ND_noisyBOLD_C;
    Resi1 = struct(f7,vi17,f1,vi11,f2,vi12,f3,vi13,f4,vi14,f5,vi15,f6,vi16);            Resc = struct(f7,vc7,f1,vc1,f2,vc2,f3,vc3,f4,vc4,f5,vc5,f6,vc6);
    Resi11(jj) = Resi1;                                                                 Resc1(jj) = Resc;
    fieldsi1 = fieldnames(Resi1);                                                       fieldsc = fieldnames(Resc);
    Resi2 = struct(f7,vi27,f1,vi21,f2,vi22,f3,vi23,f4,vi24,f5,vi25,f6,vi26);
    Resi21(jj) = Resi2;
    fieldsi2 = fieldnames(Resi2);
    Resi3 = struct(f7,vi37,f1,vi31,f2,vi32,f3,vi33,f4,vi34,f5,vi35,f6,vi36);
    Resi31(jj) = Resi3;
    fieldsi3 = fieldnames(Resi3);
    
    
    % Plot neural activity and BOLD
    figure('Name',['neural activity and BOLD, convoluted and deconvoluted - Res.: ', res{jj},' mm'],'NumberTitle','on');
    for ii=1:numel(fields)
        subplot(2,4,ii)
        hold on
        imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Resize.(fields{ii}))
        set(gca, 'ydir', 'normal')
        title({[titles{ii},];['Indice: [', int2str((Resi2.(fieldsi2{ii}))),',',int2str((Resi1.(fieldsi1{ii}))),'] Res: ', res{jj},' mm'];['t: ',int2str(t((Resi3.(fieldsi3{ii})))),' s']},'fontsize',fontT);  
        set(gca, 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
        colormap(cmap)
        colorbar
        hold off
    end

    
    % Plot neural activity and BOLD
    figure('Name',['neural activity and BOLD, convoluted and deconvoluted - Res.: ', res{jj},' mm'],'NumberTitle','on');
    for ii=1:numel(fields)
        subplot(2,4,ii)
        hold on
        imagesc(x(x0:end)*1e3,y(y0:end)*1e3,Resizet.(fieldst{ii}))
        title({[titles{ii},];['t: ',int2str(t(I3max)),' s']},'fontsize',fontT);  
        set(gca, 'ydir', 'normal', 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
        'ylim', [y(y0)*1e3, y(end)*1e3]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
        colormap(cmap)
        colorbar
        hold off
    end


    % Struct profile x, y
    fy1 = 'BOLD'; vy1 = BOLD_bary; 
    fy2 = 'Wiener'; vy2 = Wiener_bary; 
    fy4 = 'BOLD_noisyBOLD'; vy4 = BOLD_noisyBOLD_bary;
    fy5 = 'Wiener_noisyBOLD'; vy5 = Wiener_noisyBOLD_bary;
    Resizebary = struct(fy1,vy1,fy2,vy2,fy4,vy4,fy5,vy5);
    Resizebary1(jj) = Resizebary;
    fields_bary = fieldnames(Resizebary);
    y_desired = linspace(y(y0),y(end),length(BOLD_bary));
    
    fx1 = 'BOLD'; vx1 = BOLD_barx; 
    fx2 = 'Wiener'; vx2 = Wiener_barx; 
    fx4 = 'BOLD_noisyBOLD'; vx4 = BOLD_noisyBOLD_barx;
    fx5 = 'Wiener_noisyBOLD'; vx5 = Wiener_noisyBOLD_barx;
    Resizebarx = struct(fx1,vx1,fx2,vx2,fx4,vx4,fx5,vx5);
    Resizebarx1(jj) = Resizebarx;
    fields_barx = fieldnames(Resizebarx);
    x_desired = linspace(x(x0),x(end),length(BOLD_barx));
    
    
    % Plot profile x, y
    if orientation == 0 || orientation == 90
        if jj <= numel(res)

            jja=1;
            jji=1;
            jjma=1;
            jjmil=1;
            jjmir=1;
            kk = 1;
            figure('Name',['BOLD profile x,y direction - Res.: ', res{jj},' mm'],'NumberTitle','on');
            for ii = 1:numel(fields_bary)

                subplot(2,4,((ii*2)-1))
                hold on
                plot(y_desired*1e3,Resizebary.(fields_bary{ii}),'k-','LineWidth',2)
                hold off
                title({['',titlesp_res{ii}];['Res: ', res{jj}, ' mm - Prof x']},'fontsize',fontT)
                xlabel('x (mm)', 'fontsize', fontS)
                ylabel('BOLD', 'fontsize', fontS)
                set(gca, 'fontSize', fontS, 'xlim', [0, (xmax)], 'ylim', [0, y_axisy],'TickLength',[tick tick],'LineWidth',line);
                box on

                subplot(2,4,(ii*2))
                hold on
                plot(x_desired*1e3,Resizebarx.(fields_barx{ii}),'k-','LineWidth',2)
                hold off
                title({['',titlesp_res{ii}];['Res: ', res{jj}, ' mm - Prof y']},'fontsize',fontT)
                xlabel('x (mm)', 'fontsize', fontS)
                ylabel('BOLD', 'fontsize', fontS)
                set(gca, 'fontSize', fontS, 'xlim', [0, (xmax)], 'ylim', [0, y_axisx],'TickLength',[tick tick],'LineWidth',line);
                box on

                Maximum_res(jj,ii) = max(Resizebary.(fields_bary{ii}));

                    % Modulation of profile
                    jja = 1;
                    for kka = 2:(length(Resizebary.(fields_bary{ii}))-1)

                    % find max's and store in new array
                        store_pre_max_res = Resizebary.(fields_bary{ii})(kka-1);
                        store_pos_max_res = Resizebary.(fields_bary{ii})(kka+1);
                        if Resizebary.(fields_bary{ii})(kka) > store_pre_max_res && Resizebary.(fields_bary{ii})(kka) > store_pos_max_res
                            store_max_res = Resizebary.(fields_bary{ii})(kka);
                            store_max1_res(jj,ii,jja) = store_max_res;
                            jja=jja+1;
                        end
                    end

                    % find min's and store in new array
                        jji=1;
                    for kki=2:(length(Resizebary.(fields_bary{ii}))-1)
                        store_pre_min_res = Resizebary.(fields_bary{ii})(kki-1);
                        store_pos_min_res = Resizebary.(fields_bary{ii})(kki+1);
                        if Resizebary.(fields_bary{ii})(kki) < store_pre_min_res && Resizebary.(fields_bary{ii})(kki) < store_pos_min_res
                            store_min_res = Resizebary.(fields_bary{ii})(kki);
                            store_min1_res(jj,ii,jji) = store_min_res;
                            jji=jji+1;
                        end
                    end

                    % finding maximum values of the max array
                    % finding max value neighbour to max value
                    for mm=1:size(store_max1_res,1)
                        for nn=1:size(store_max1_res,2)
                            [Maxmaxii,IndMax_res] = max(store_max1_res(mm,nn,:));
                            if (IndMax_res+1) <= size(store_max1_res,3) %%
                                Maxmaxii_pos = (store_max1_res(mm,nn,IndMax_res+1));
                                Maxmax(mm,nn) = Maxmaxii;
                                Maxmax_pos(mm,nn) = Maxmaxii_pos;
                            end %%
                            % finding minimum value between two max values
                            if mm <= size(store_min1_res,1) && nn <= size(store_min1_res,2) && IndMax_res <= size(store_min1_res,3) %%
                                MaxMinii = store_min1_res(mm,nn,IndMax_res);
                                MaxMin_pos(mm,nn) = MaxMinii;
                            end
                        end
                    end


                    [Maximum_res,ind_max_res]=max(store_max1_res(jj,ii,:));
                    % modulation related to maximum value
                    jjma = 1;
                    for kkma=2:length(store_max1_res)
                        if kkma <= ind_max_res && store_max1_res(jj,ii,kkma) ~= 0 && store_max1_res(jj,ii,kkma-1) ~= 0
                           modul_max_res = 100-((store_max1_res(jj,ii,kkma-1)*100)/store_max1_res(jj,ii,kkma));
                           modul1_max_res(jj,ii,jjma) = modul_max_res;
                           jjma=jjma+1;        
                        end
                    end

                kk=kk+1;
            end
    
        end
    end
    
    
    % Struct diagonal
    fc1 = 'BOLD'; vc1 = BOLD_diagcr;
    fc2 = 'Wiener'; vc2 = Wiener_diagcr; 
    fc4 = 'BOLD_noisyBOLD'; vc4 = BOLD_noisyBOLD_diagcr;
    fc5 = 'Wiener_noisyBOLD'; vc5 = Wiener_noisyBOLD_diagcr;
    Resizediagc = struct(fc1,vc1,fc2,vc2,fc4,vc4,fc5,vc5);
    Resizediagc1(jj) = Resizediagc;
    fields_diagc = fieldnames(Resizediagc);
    
    fd1 = 'BOLD'; vd1 = BOLD_diagdecr; 
    fd2 = 'Wiener'; vd2 = Wiener_diagdecr; 
    fd4 = 'BOLD_noisyBOLD'; vd4 = BOLD_noisyBOLD_diagdecr;
    fd5 = 'Wiener_noisyBOLD'; vd5 = Wiener_noisyBOLD_diagdecr;
    Resizediagd = struct(fd1,vd1,fd2,vd2,fd4,vd4,fd5,vd5);
    Resizediagd1(jj) = Resizediagd;
    fields_diagd = fieldnames(Resizediagd);
    
    % Plot profile diagonal
    kkd = 1;
    if orientation == 45 || orientation == 135
        if jj <= numel(res)
            figure('Name',['BOLD profile diagonal direction - Res.: ', res{jj},' mm'],'NumberTitle','on');
            for ii = 1:numel(fields_diagc)

                % Visibility and modulation
                y_maxd_res = dsearchn(y_desired',7.875*1e-3);
                x_maxd_res = Resizediagd.(fields_bary{ii})(y_maxd_res);
                y_minld_res = dsearchn(y_desired',6.125*1e-3);
                x_minld_res = Resizediagd.(fields_bary{ii})(y_minld_res);
                y_minrd_res = dsearchn(y_desired',9.625*1e-3);
                x_minrd_res = Resizediagd.(fields_bary{ii})(y_minrd_res);

                modulationld_res = 100-((x_minld_res*100)/x_maxd_res);
                modulationld1_res(jj,kkd) = modulationld_res;
                visibilityld_res = (x_maxd_res-x_minld_res)/(x_maxd_res+x_minld_res);
                visibilityld1_res(jj,kkd) = visibilityld_res;
                modulationrd_res = 100-((x_minrd_res*100)/x_maxd_res);
                modulationrd1_res(jj,kkd) = modulationrd_res;
                visibilityrd_res = (x_maxd_res-x_minrd_res)/(x_maxd_res+x_minrd_res);
                visibilityrd1_res(jj,kkd) = visibilityrd_res;

                jja = 1;
                    for kka = 2:(length(Resizediagd.(fields_diagd{ii}))-1)

                        % find max's and store in new array
                        store_pre_max_resd = Resizediagd.(fields_diagd{ii})(kka-1);
                        store_pos_max_resd = Resizediagd.(fields_diagd{ii})(kka+1);
                        if Resizediagd.(fields_diagd{ii})(kka) > store_pre_max_resd && Resizediagd.(fields_diagd{ii})(kka) > store_pos_max_resd
                            store_max_resd = Resizediagd.(fields_diagd{ii})(kka);
                            store_max1_resd(jj,ii,jja) = store_max_resd;
                            jja=jja+1;
                        end
                    end

                    % find min's and store in new array
                        jji=1;
                    for kki=2:(length(Resizediagd.(fields_diagd{ii}))-1)
                        store_pre_min_resd = Resizediagd.(fields_diagd{ii})(kki-1);
                        store_pos_min_resd = Resizediagd.(fields_diagd{ii})(kki+1);
                        if Resizediagd.(fields_diagd{ii})(kki) < store_pre_min_resd && Resizediagd.(fields_diagd{ii})(kki) < store_pos_min_resd
                            store_min_resd = Resizediagd.(fields_diagd{ii})(kki);
                            store_min1_resd(jj,ii,jji) = store_min_resd;
                            jji=jji+1;
                        end
                    end

                    % finding maximum values of the max array
                    % finding max value neighbour to max value
                    for mm=1:size(store_max1_resd,1)
                        for nn=1:size(store_max1_resd,2)
                            [Maxmaxiid,IndMax_resd] = max(store_max1_resd(mm,nn,:));
                            Maxmaxii_posd = (store_max1_resd(mm,nn,IndMax_resd+1));
                            Maxmaxd(mm,nn) = Maxmaxiid;
                            Maxmax_posd(mm,nn) = Maxmaxii_posd;
                            % finding minimum value between two max values
                            if mm <= size(store_min1_resd,1) && nn <= size(store_min1_resd,2)
                                MaxMiniid = store_min1_resd(mm,nn,IndMax_resd);
                                MaxMin_posd(mm,nn) = MaxMiniid;
                            end
                        end
                    end



                    [Maximum_resd,ind_max_resd]=max(store_max1_resd(jj,ii,:));
                    % modulation related to maximum value
                    jjma = 1;
                    for kkma=2:length(store_max1_resd)
                        if kkma <= ind_max_resd && store_max1_resd(jj,ii,kkma) ~= 0 && store_max1_resd(jj,ii,kkma-1) ~= 0
                           modul_max_resd = 100-((store_max1_resd(jj,ii,kkma-1)*100)/store_max1_resd(jj,ii,kkma));
                           modul1_max_resd(jj,ii,jjma) = modul_max_resd;
                           jjma=jjma+1;        
                        end
                    end

                subplot(2,4,((ii*2)-1))
                hold on
                plot(y_desired*1e3,Resizediagc.(fields_diagc{ii}),'k-','LineWidth',fontL)
                hold off
                title({['',titlesp_res{ii}];['Res: ', res{jj}, ' mm - Prof diag \']},'fontsize',fontT)
                xlabel('x (mm)', 'fontsize', fontS)
                ylabel('BOLD', 'fontsize', fontS)
                set(gca, 'fontSize', fontS, 'xlim', [-0.2, (xmax)], 'ylim', [0, y_axisy],'TickLength',[tick tick],'LineWidth',line);
                box on

                subplot(2,4,(ii*2))
                hold on
                plot(x_desired*1e3,Resizediagd.(fields_diagd{ii}),'k-','LineWidth',fontL)
                hold off
                title({['',titlesp_res{ii}];['Res: ', res{jj}, ' mm - Prof diag /']},'fontsize',fontT)
                xlabel('x (mm)', 'fontsize', fontS)
                ylabel('BOLD', 'fontsize', fontS)
                set(gca, 'fontSize', fontS, 'xlim', [-0.2, (xmax)], 'ylim', [0, y_axisx],'TickLength',[tick tick],'LineWidth',line);
                box on

                kkd = kkd + 1;

            end
    
        end
    end
    
    % Maximum voxel along time
    vmv7 = ND_r_mvec;
    vmv1 = BOLD_mvec;
    vmv2 = Wiener_mvec;
    vmv3 = ND_mvec;
    vmv4 = BOLD_noisyBOLD_mvec;
    vmv5 = Wiener_noisyBOLD_mvec;
    vmv6 = ND_noisyBOLD_mvec;
    mvec = struct(f7,vmv7,f1,vmv1,f2,vmv2,f3,vmv3,f4,vmv4,f5,vmv5,f6,vmv6);
    mvec1(jj) = mvec;
    fields_mvec = fieldnames(mvec);
    
    figure('Name',['Profile of maximum voxel along time - Res.: ', res{jj},' mm'],'NumberTitle','on');
    for ii = 1:numel(fields_mvec)
        subplot(2,4,ii)
        plot(t(t0:end),mvec.(fields_mvec{ii})(t0:end),'k-','LineWidth',fontL)
        xlabel('t (s)', 'fontsize', fontS)
        ylabel([titles_stim{ii}], 'fontsize', fontS)
        title({['Max voxel along t'];['Res: ', res{jj}, ' mm']},'fontsize',fontT)
        set(gca, 'fontSize', fontS, 'xlim', [0, tmax],'TickLength',[tick tick],'LineWidth',line);
    end

    
    % neural activity and BOLD response mean 'on'
    figure('Name',['neural activity and BOLD, convoluted and deconvoluted, mean on - Res.: ', res{jj},' mm'],'NumberTitle','on');
    for ii = 1:numel(fields_mean)
    
        res_average_on = squeeze(mean(Resize_mean.(fields_mean{ii}),3));
        subplot(2,4,ii)
        imagesc(x(x0:end)*1e3,y(y0:end)*1e3,res_average_on);
        set(gca, 'ydir', 'normal');
        title({[titles{ii}];['Mean on. Res.: ', res{jj},' mm']},'fontsize',fontT);
        set(gca, 'fontSize', fontS, 'xlim', [x(x0)*1e3, x(end)*1e3], ...
            'ylim', [y(y0)*1e3, y(end)*1e3]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
        colormap(cmap)
        colorbar

    end
    
    jj= jj+1;
end

% Visibility OP 0 or 90:
% ((0.5*(peak1+peak2))-valey1)/((0.5*(peak1+peak2))+valey1);
visibility_res = zeros(size(Maxmax));
for vv1 = 1:size(Maxmax,1)
    for vv2=1:size(Maxmax,2)
        visibility_res(vv1,vv2) = ((0.5*(Maxmax(vv1,vv2)+Maxmax_pos(vv1,vv2)))-MaxMin_pos(vv1,vv2))/((0.5*(Maxmax(vv1,vv2)+Maxmax_pos(vv1,vv2)))+MaxMin_pos(vv1,vv2));
    end
end

% Visibility OP 45 or 135:
% ((0.5*(peak1+peak2))-valey1)/((0.5*(peak1+peak2))+valey1);
visibility_resd = zeros(size(Maxmaxd));
for vv1 = 1:size(Maxmaxd,1)
    for vv2=1:size(Maxmaxd,2)
        visibility_resd(vv1,vv2) = ((0.5*(Maxmaxd(vv1,vv2)+Maxmax_posd(vv1,vv2)))-MaxMin_posd(vv1,vv2))/((0.5*(Maxmaxd(vv1,vv2)+Maxmax_posd(vv1,vv2)))+MaxMin_posd(vv1,vv2));
    end
end

%% memory usage

allvars = whos;
memused = sum([allvars.bytes])

toc
