function [matrix_rowx,matrix_coly] = profile(matrix)

% d55_BF = deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,time_spe);
% [row, col] = find(ismember(d55_BF, max(d55_BF(:))));

%matrix = deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,time_spe);

[row, col] = find(ismember(matrix, max(matrix(:))));
matrix_rowx = matrix(row,:);
matrix_coly = matrix(:,col);

end