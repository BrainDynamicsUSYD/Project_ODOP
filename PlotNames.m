function PlotNames

% cmap_new = colormap_bluetored;
% cmap = abs(1-hot);
cmap = jet;
% cmap = 1-gray;

xlabel('mm')
ylabel('mm')
colormap(cmap)
colorbar
end