function [matrix_new, diagcr, diagdecr] = Diagonal(matrix)

% Calculate diagonal of matrix making sure that maximum element is at
% center position

% matrix = deconvResponses_Forward_2D.reconvBOLD(x0:end,y0:end,time_spe);

% Marilia M Oliveira, 2018

% find max row and col of matrix
[row, col] = find(ismember(matrix, max(matrix(:)))); 

% center of matrix
center_loc=num2cell(ceil(size(matrix)/2));
center = matrix(sub2ind(size(matrix),center_loc{:}));
[rowc,colc] = find(matrix==center);
%equivalent to matrix(center_loc{1},center_loc{2},...,center_loc{n})

% how many positions need to shift the matrix
rowm = (rowc-row);
colm = (colc-col);

% Shifted matrix with maximum value in the center; possible to calculate
% main diagonal
matrix_new = circshift(matrix,[rowm(1,:),colm(1,:)]);

% repeting same process because when we flip the matrix, not necessary the
% maximum value stay in the center anymore.
matrix_new1 = fliplr(matrix_new);
[row1, col1] = find(ismember(matrix_new1, max(matrix_new1(:))));
center_loc1=num2cell(ceil(size(matrix_new1)/2));
center1 = matrix_new1(sub2ind(size(matrix_new1),center_loc1{:}));
[rowc1,colc1] = find(matrix_new1==center1);
rowm1 = (rowc1-row1);
colm1 = (colc1-col1);
matrix_new2 = circshift(matrix_new1,[rowm1(1,:),colm1(1,:)]);


% diagonal crescent: left bottom to right top /
diagcr = diag(matrix_new2);
% diagonal decrescent: left top to right bottom \
diagdecr = diag(matrix_new);

end
