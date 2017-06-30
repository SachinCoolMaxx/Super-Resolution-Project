% Image construction from overlapping patches
function [result] = overlap_add3D(patches, img_size, grid)

result = zeros(img_size);
weight = zeros(img_size);
size(grid,4)
for i = 1:size(grid, 4)
    patch = reshape(patches(:, i), size(grid, 1), size(grid, 2), size(grid,3));
    result(grid(:, :, :, i)) = result(grid(:, :, :, i)) + patch;
    weight(grid(:, :,:, i)) = weight(grid(:, :,:, i)) + 1;
end

I = logical(weight);
result(I) = result(I) ./ weight(I);
