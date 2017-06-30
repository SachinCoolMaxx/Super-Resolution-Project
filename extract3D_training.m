function [features] = extract3D_training(conf, A1, A2, scale, filters)

% Compute one grid for all filters for input image X
grid = sampling_grid3D(A1, size(A1), conf.intensity_threshold,...
    conf.window, conf.overlap, conf.border, scale);
feature_size = prod(conf.window) * numel(conf.filters);
% Current image features extraction [feature x index]
if isempty(filters)
    X = A2 - A1;
    f = X(grid);
    features = reshape(f, [size(f, 1)*size(f, 2)*size(f,3) size(f, 4)]);
else
    X = A1;
    features = zeros([feature_size size(grid, 4)], 'single');
    for i = 1:numel(filters)
        disp('filter')
        i
        f = convn(X, filters{i}, 'same');
        f = f(grid);
        f = reshape(f, [size(f,1)*size(f,2)* size(f,3) size(f, 4)]);
        features((1:size(f, 1)) + (i-1)*size(f, 1), :) = f;
    end
end

