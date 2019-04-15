% Marilia 20/12/2017

% function [z] = pearl(xx,yy,A,xp,yp,s_x,s_y,a)
% % Xiaochen equation. Mine has a slightly change and the results for
% convolution and deconvolution will vary.
% % A = angle in degree
% % B = angle in degree B=0
% % xp = x' displacement from the centre in x direction
% % yp = y' displacement from the centre in y direction
% % s_x = sqrt(8); 
% % s_y = 1; 
% % a = 2 width of the unit cell [mm]
% 
% x_g = ((xx-xp).*cosd(A)+(yy-yp).*sind(A));
% y_g = (-(xx-xp).*sind(A)+(yy-yp).*cosd(A));
% 
% ff=(1/(2*pi*s_x*s_y)).*(exp(-0.5*(((x_g.^2)./s_x.^2)+((y_g.^2)./s_y.^2))));
% 
% k = (2*pi)/a; % spatial period
% basic_fn=((0.5*((cos(k*(xx-xp)))+1))).*((0.5*((cos(k*(yy-yp)))+1)));
% 
% z=basic_fn.*ff;
% end

function [z] = pearl(xx,yy,A,xp,yp,s_x,s_y,a)

% A = angle in degree
% B = angle in degree B=0
% xp = x' displacement from the centre in x direction
% yp = y' displacement from the centre in y direction
% s_x = 16; % sigma_x^2; s_x = 8
% s_y = 4.5; % sigma_y^2; s_y = 1
% a = 2 width of the unit cell [mm]

% rotation
x_g = ((xx-xp).*cosd(A)+(yy-yp).*sind(A));
y_g = (-(xx-xp).*sind(A)+(yy-yp).*cosd(A));

% Gaussian envelope
ff=(exp(-0.5*(((x_g.^2)./s_x)+((y_g.^2)./s_y))));

k = (2*pi)/a; % spatial period
% modulation/repetition
basic_fn=((0.5*(cos(k*(xx-xp))+1))).*((0.5*(cos(k*(yy-yp))+1)));

z=basic_fn.*ff;
end
