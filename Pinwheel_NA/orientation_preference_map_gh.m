% Xiaochen 30/11/2017
% Marilia 18/01/2018

N=129; %number of points
% N=257; 

%{
pinwheel location
--------------------------
          y axis
            ^
 --  -   |  |  |   +   --
            |
------------|-------------> x axix
            |
 --  +   |  |  |   -   --   

---------------------------
'-' sign stands for clockwise rotation
'+' sign stands for counter clockwise rotation

%}
%==========================================================================
% STEP 1: Compute the right-top pinwheel matrix
% Size of the fundamental domain 
X_rt=(0:1/(N-1):1);
Y_rt=(0:1/(N-1):1);

% Pinwheel centre
x_pw = 0.5; % mm
y_pw = 0.5; % mm


% Shift coordinates
Xc = X_rt - x_pw;
Yc = Y_rt - y_pw;

%%

% Grid of the domain
[X, Y] = meshgrid(Xc, Yc);

% theta
A_theta        = atan2d(Y, X);

xy_idx = ceil(N/2);
% A_rt(xy_idx, xy_idx) = -inf;

cMap=colormap('hsv');
cMap=rgb2hsv(cMap);
cMap(:,2)=0.25;
colorMap=hsv2rgb(cMap);
colormap(colorMap);

figure(1);
imagesc(Xc, Yc, A_theta); hold on;
title('\theta(x, y) = atan2d((y-y_c) / (x-x_c))')
colormap(colorMap);
axis equal
xlim([-0.5, 0.5])
ylim([-0.5, 0.5])

xlabel('x-x_c [mm]')
ylabel('y-y_c [mm]')
set(gca, 'YDir', 'normal');
colorbar()
%%

% Shift range to 0 - 180 
% phi = 1/2 * (theta + 180)
A_rt = 0.5*(atan2d(Y, X) + 180);

%% Make flipped and rotated copies

% To get the unit cell in the right direction we need to flip the
% fundamental domain
%     y|
%      |
%      |
% (0,0) ---->x
A_rt_flipped = flipud(A_rt);

tr = A_rt_flipped;
tl = fliplr(A_rt_flipped);
bl = rot90(rot90(A_rt_flipped));
br = flipud(A_rt_flipped);

%% Concatenate

tr = tr(1:end-1, 2:end);
tl = tl(1:end-1, :); % This one keeps the origin
% bl does not change
br = br(:, 2:end);
PMM = [tl tr; bl br]; 

% Length of primitive cell (4 unit cells)
X_pc = linspace(-1,1,size(PMM,1)); 
Y_pc = linspace(-1,1,size(PMM,1)); 
[XX, YY] = meshgrid(X_pc, Y_pc);

%% Pinwheel
fontS = 18;

figure
imagesc(X_pc, Y_pc, PMM); hold on;
colormap(colorMap);
axis equal
xlim([-1.0, 1.0])
ylim([-1.0, 1.0])
plot([-1.0, 1.0], [0, 0], 'k')
plot([0.0, 0.0], [-1.0, 1.0], 'k')
plot([-1.0, 1.0], [1.0, -1.0], 'k--')
plot([-1.0, 1.0], [-1.0, 1.0], 'k--')
ylabel('y (mm)', 'fontsize', fontS)
xlabel('x (mm)', 'fontsize', fontS)
set(gca, 'YDir', 'normal','fontSize', fontS);
c = colorbar();
c.Label.String = 'OP (degrees)';
c.Ticks = [0,45,90,135,180];

        
%% coordinates and angles of segments in unit cell added in pinwheel from above

XC = [-0.8, -0.5, -0.5, -0.25, -0.8, -0.70,  -0.3,  -0.2, -0.5,   0];
YC = [-0.5, -0.8, -0.2, -0.3, -0.7,  -0.25,  -0.7, -0.5,     0, -0.7];
TH = [  90,  135,   45,   15,  105,     65,  155,     0,    45, 170];

for kk=1:length(XC)
    % Lines at 90
    [xcoo, ycoo] = oriented_segment(XC(kk), YC(kk), -TH(kk), 0.05);
    line(xcoo, ycoo, 'color', 'k', 'linewidth', 2)
    line(xcoo, ycoo+2*abs(YC(kk)),  'color', 'k', 'linewidth', 2)
    line(xcoo+2*abs(XC(kk)), ycoo,  'color', 'k', 'linewidth', 2)
    line(xcoo+2*abs(XC(kk)), ycoo+2*abs(YC(kk)),  'color', 'k', 'linewidth', 2)
end


%% Pinwheel rotated

% Fixing the angles of segments and pinwheel.
% changing the position of the pinwheel to match the coordinates and angles 
% of segments.

%% 4 unit cells

% Pinwheel rotated
% 4 unit cells that considering bilaterality, it is actually 1 unit cell
% for left and right eye.
PM = [br bl; tr tl];
XPM = linspace(-1,1,size(PM,1));
YPM = linspace(-1,1,size(PM,1));
[XXPM, YYPM] = meshgrid(XPM,YPM);

%% 16 unit cells

% Pinwheel rotated
% 16 unit cells ~= 4 unit cells
PM2 = [PM PM; PM PM];
XPM2 = linspace(-2,2,size(PM2,1));
YPM2 = linspace(-2,2,size(PM2,1));
[XXPM2, YYPM2] = meshgrid(XPM2,YPM2);

%% 144 (12 x 12) unit cells (3*PM2)

% Pinwheel rotated
% 144 (12 x 12) unit cells (3*PM2)

% Length of primitive cell (144 unit cells)
PM12 = [PM2 PM2 PM2; PM2 PM2 PM2; PM2 PM2 PM2]; % concatenate
XPM12 = linspace(-6,6,size(PM12,1)); 
YPM12 = linspace(-6,6,size(PM12,1)); 
[XXPM12, YYPM12] = meshgrid(XPM12, YPM12);

% Pinwheel
figure
imagesc(XPM12, YPM12, PM12); hold on;
colormap(colorMap);
axis equal
xlim([-6.0, 6.0])
ylim([-6.0, 6.0])
plot([-6.0, 6.0], [0, 0], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [-1, -1], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [1, 1], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [-2, -2], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [2, 2], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [-3, -3], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [3, 3], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [-4, -4], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [4, 4], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [-5, -5], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [5, 5], 'Color',[0 0 0]+0.6)
plot([0.0, 0.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([-1.0, -1.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([1.0, 1.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([-2.0, -2.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([2.0, 2.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([-3.0, -3.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([3.0, 3.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([-4.0, -4.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([4.0, 4.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([-5.0, -5.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([5.0, 5.0], [-6.0, 6.0], 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [-6.0, 6.0], '--', 'Color',[0 0 0]+0.6)
plot([-6.0, 6.0], [6.0, -6.0], '--', 'Color',[0 0 0]+0.6)
plot([-4.0, -6.0], [6.0, 4.0], '--', 'Color',[0 0 0]+0.6)
plot([-2.0, -6.0], [6.0, 2.0], '--', 'Color',[0 0 0]+0.6)
plot([0.0, -6.0], [6.0, 0.0], '--', 'Color',[0 0 0]+0.6)
plot([2.0, -6.0], [6.0, -2.0], '--', 'Color',[0 0 0]+0.6)
plot([4.0, -6.0], [6.0, -4.0], '--', 'Color',[0 0 0]+0.6)
plot([-4.0, 6.0], [-6.0, 4.0], '--', 'Color',[0 0 0]+0.6)
plot([-2.0, 6.0], [-6.0, 2.0], '--', 'Color',[0 0 0]+0.6)
plot([0.0, 6.0], [-6.0, 0.0], '--', 'Color',[0 0 0]+0.6)
plot([2.0, 6.0], [-6.0, -2.0], '--', 'Color',[0 0 0]+0.6)
plot([4.0, 6.0], [-6.0, -4.0], '--', 'Color',[0 0 0]+0.6)
plot([-4.0, -6.0], [-6.0, -4.0], '--', 'Color',[0 0 0]+0.6)
plot([-2.0, -6.0], [-6.0, -2.0], '--', 'Color',[0 0 0]+0.6)
plot([0.0, -6.0], [-6.0, 0.0], '--', 'Color',[0 0 0]+0.6)
plot([2.0, -6.0], [-6.0, 2.0], '--', 'Color',[0 0 0]+0.6)
plot([4.0, -6.0], [-6.0, 4.0], '--', 'Color',[0 0 0]+0.6)
plot([-4.0, 6.0], [6.0, -4.0], '--', 'Color',[0 0 0]+0.6)
plot([-2.0, 6.0], [6.0, -2.0], '--', 'Color',[0 0 0]+0.6)
plot([0.0, 6.0], [6.0, 0.0], '--', 'Color',[0 0 0]+0.6)
plot([2.0, 6.0], [6.0, 2.0], '--', 'Color',[0 0 0]+0.6)
plot([4.0, 6.0], [6.0, 4.0], '--', 'Color',[0 0 0]+0.6)
% ylabel('y (mm)','fontSize', fontS)
% xlabel('x (mm)','fontSize', fontS)
set(gca, 'YDir', 'normal','fontSize', fontS)
set(gca, 'fontSize', fontS, 'XTick', [-5 0 5], ...
            'YTick', [-5 0 5]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
c = colorbar();
c.Label.String = 'OP (degrees)';
c.Ticks = [0,45,90,135,180];

% coordinates and angles of segments in unit cell added in pinwheel from above
% fix position of angles
XC = [-0.8, -0.5, -0.5, -0.25, -0.8, -0.70,  -0.3,  -0.2, -0.5,   0];
YC = [-0.5, -0.8, -0.2, -0.3, -0.7,  -0.25,  -0.7, -0.5,     0, -0.7];
TH = [  90,  135,   45,   15,  105,     65,  155,     0,    45, 170];

for kk=1:length(XC)
    % Lines at 90
    [xcoo, ycoo] = oriented_segment(XC(kk), YC(kk), -TH(kk), 0.05);
    line(xcoo+8, ycoo+8, 'color', 'k', 'linewidth', 2)
    line(xcoo+8, ycoo+8+2*abs(YC(kk)),  'color', 'r', 'linewidth', 2)
    line(xcoo+8+2*abs(XC(kk)), ycoo+8,  'color', 'b', 'linewidth', 2)
    line(xcoo+8+2*abs(XC(kk)), ycoo+8+2*abs(YC(kk)),  'color', 'm', 'linewidth', 2)
end

% Pearl stimulation

% 0º
A=0;
xp=XC(8);
yp=YC(8);
s_x=(8);
s_y=1;
a=2;
[z0]=pearl(XXPM12,YYPM12,A,xp,yp,s_x,s_y,a);
contour(XXPM12,YYPM12,z0,'k')

% % 45º
% A=45;
% xp=XC(3);
% yp=YC(3);
% s_x=(8);
% s_y=(1);
% a=2;
% [z45]=pearl(XXPM12,YYPM12,A,xp,yp,s_x,s_y,a);
% contour(XXPM12,YYPM12,z45,'k')

%%
% Neural activity

fontS = 18;

A=0;
xp=XC(8);
yp=YC(8);
s_x=(8);
s_y=(1);
a=2;
[z0]=pearl(XXPM12,YYPM12,A,xp,yp,s_x,s_y,a);
figure,
%imagesc(XPM12,YPM12,z0)
contour(XXPM12,YYPM12,z0)
set(gca, 'fontSize', fontS, 'XTick', [-5 0 5], ...
            'YTick', [-5 0 5]);
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
colormap(jet)
%colormap(1-gray)
%colorbar('Ticks',[0.5])
colorbar

A=45;
xp=XC(3);
yp=YC(3);
s_x=(8);
s_y=(1);
a=2;
[z0]=pearl(XXPM12,YYPM12,A,xp,yp,s_x,s_y,a);
figure,
%imagesc(XPM12,YPM12,z0)
contour(XXPM12,YYPM12,z0)
set(gca, 'fontSize', fontS, 'XTick', [-5 0 5], ...
            'YTick', [-5 0 5],'ydir', 'normal');
        xlabel('x (mm)', 'fontsize', fontS)
        ylabel('y (mm)', 'fontsize', fontS)
colormap(jet)
%colormap(gray)
%colorbar('Ticks',[0.5])
colorbar

%%
% Pinwheel + neural activity

fontS1 = 14;
figure,
ax1 = axes;
imagesc(ax1,XPM12,YPM12,PM12);
view(2)
ax2 = axes;
% A=45;
% xp=XC(3);
% yp=YC(3);
A=0;
xp=XC(8);
yp=YC(8);
s_x=(8);
s_y=(1);
a=2;
[z0]=pearl(XXPM12,YYPM12,A,xp,yp,s_x,s_y,a);
contour(ax2,XPM12,YPM12,z0)
%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,colorMap)
colormap(ax2,'jet')
%%Then add colorbars and get everything lined up
set([ax1,ax2], 'fontSize', fontS1, 'XTick', [-5 0 5], ...
            'YTick', [-5 0 5],'YDir', 'normal');
xlabel(ax1,'x (mm)', 'fontsize', fontS1)
ylabel(ax1,'y (mm)', 'fontsize', fontS1)
set([ax1,ax2],'Position',[.20 .15 .63 .75]);
cb1 = colorbar(ax1,'Position',[.85 .15 .035 .75],'Ticks',[0,45,90,135,180]);
cb2 = colorbar(ax2,'Position',[.06 .15 .035 .75],'Ticks',[0.2,0.4,0.6,0.8]);
cb1.Label.String = 'OP (degrees)';
