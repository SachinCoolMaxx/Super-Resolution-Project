% Returns a sampling grid for an image (using overlapping window)
%   GRID = SAMPLING_GRID(IMG_SIZE, WINDOW, OVERLAP, BORDER, SCALE)
function grid = sampling_grid3D(img, img_size,intensity_threshold, window, overlap, border, scale)

fg = zeros(size(img));

% set threshold for image intensity values to be considered fg. 
fg(img > intensity_threshold) = 1;

if nargin < 5
    scale = 1;
end

if nargin < 4
    border = [0 0];   
end

if nargin < 3
    overlap = [0 0];    
end

% Scale all grid parameters
% window = window.* scale;
% overlap = overlap.*scale;
% border = border * scale;

% Create sampling grid for overlapping window
index = reshape(1:prod(img_size), img_size);
grid = index(1:window(1), 1:window(2), 1:window(3)) - 1;

% Compute offsets for grid's displacement.
skip = window - overlap; % for small overlaps
skip(skip == 0) = 1;
offset = index(1+border(1):skip(1):img_size(1)-window(1)+1-border(1), ...
               1+border(2):skip(2):img_size(2)-window(2)+1-border(2),...
               1+border(3):skip(3):img_size(3)-window(3)+1-border(3));
fgoffset = fg(1+border(1):skip(1):img_size(1)-window(1)+1-border(1), ...
               1+border(2):skip(2):img_size(2)-window(2)+1-border(2),...
               1+border(3):skip(3):img_size(3)-window(3)+1-border(3));
           
offset = reshape(offset, [1 1 1 numel(offset)]);
fgoffset = reshape(fgoffset,[1 1 1 numel(fgoffset)]);
fg_idxs = find(fgoffset == 1);
offset_fg = offset(fg_idxs);

% Prepare 4D grid - should be used as: sampled_img = img(grid);
numel(offset_fg)
[window 1]
grid = repmat(grid, [1 1 1 numel(offset_fg)]) + repmat(offset_fg, [window 1]);
