function [res_func1,prof_func,row,bary_func,barx_func,diag_cr,diag_decr,max_vector,Cm,I1m,I2m,I3m] = ChangeResolutionNoise(func_name,scale)

% Function to change resolution of matrix - BOLD and neural activity.

% scale values for  params.Nkx = 2^8; 
%                   params.Nky = 2^8;
%                   params.Nw = 2^6;
% scale = 0.5;    (1/2)   % Resolution 0.25 mm (0.25/0.125 = 2)
% scale = 0.25;   (1/4)   % Resolution 0.5 mm  (0.5/0.125 = 4)
% scale = 0.1667; (1/6)   % Resolution 0.75 mm (0.75/0.125 = 6)
% scale = 0.125;  (1/8)   % Resolution 1.0 mm  (1/0.125 = 8)
% scale = 0.0833; (1/12)  % Resolution 1.5 mm  (1.5/0.125 = 12)

% Marilia M Oliveira, 2018

% resize
res_func = imresize(func_name,scale);
[row,col] = size(res_func(:,:,1));
res_func1 = res_func((row/2):end-1,(col/2):end-1,:);

% position of maximum magnitude
[C,I] = max(res_func(:));
[I1,I2,I3] = ind2sub(size(res_func),I);

% profile
prof_func = res_func((row/2):end-1,(col/2):end-1,I3);

% x,y profile
bary_func = res_func(I1,(col/2):end-1,I3);
barx_func = res_func((row/2):end-1,I2,I3);

% diagonal profile
% diag_cr = diag(fliplr(prof_func)); % diagonal crescent: left bottom to right top /
% diag_decr = diag(prof_func); % diagonal decrescent: left top to right bottom \

[matrix_new, diag_cr, diag_decr] = Diagonal(prof_func);

% maximum voxel of matrix
prof_func1 = res_func((row/2):end-1,(col/2):end-1,:);

[Cm,Im] = max(prof_func1(:));
[I1m, I2m,I3m] = ind2sub(size(prof_func1),Im);
vector = prof_func1(I1m,I2m,:);
max_vector = reshape(vector, size(vector,3), 1);

end
