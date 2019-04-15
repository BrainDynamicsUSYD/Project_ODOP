% Xiaochen 30/11/2017

function [xcoords, ycoords] = oriented_segment(xc, yc, theta, len)
% Gives back the coordinates of an oriented line segment
    cosOrient = cosd(theta);
    sinOrient = sind(theta);
    xcoords = xc + len * [cosOrient -cosOrient];
    ycoords = yc + len * [-sinOrient sinOrient];
end