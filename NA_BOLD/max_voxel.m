function [max_vector, I1, I2, I3, C] = max_voxel(matrix)

% Find max voxel of matrix and plot along z direction (time)
% matrix = deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,:);

[C,I] = max(matrix(:));
[I1,I2,I3] = ind2sub(size(matrix),I);
vector = matrix(I1,I2,:);
max_vector = reshape(vector, size(vector,3), 1);

% plot(t(t0:end),max_vector(t0:end));

end
