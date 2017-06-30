function [conf] = funcLearnANRDict(A1, A2,upscaling, N, intensity_threshold, dict_size, border)



conf.scale = upscaling; % scale-up factor
conf.level = 1; % # of scale-ups to perform
conf.window = N; % low-res. window size
conf.border = border; % border of the image (to ignore)

% High-pass filters for feature extraction (defined for upsampled low-res.)
conf.upsample_factor = upscaling;
 
% 
% Gy = [1 Oy -1]; % Gradient
% Ly = [1 Oy -2 Oy 1]/2; % Laplacian
% 
% Gx = [1 Ox -1]; % Gradient
% Lx = [1 Ox -2 Ox 1]/2; % Laplacian
% 
Gy = [1 0 -1]; % Gradient
Ly = [1 0 -2 0 1]/2; % Laplacian

Gx = [1 0 -1]; % Gradient
Lx = [1 0 -2 0 1]/2; % Laplacian



Gz(:,:,1) = 1*ones(3,3);
Gz(:,:,2) = -0*ones(3,3);

% not used but can be added as features
Gz(:,:,3) = -1*ones(3,3);
Lz = fspecial3('laplacian');

conf.filters = {Gy, Gx.', Ly, Lx.'}; % 2D versions
conf.interpolate_kernel = 'bicubic';

conf.overlap = [1,1,1]; % partial overlap (for faster training)

conf.upsample_factor
conf.intensity_threshold = intensity_threshold; % ignore the values less than threshold


features = extract3D_training(conf,A1,A2, conf.upsample_factor, conf.filters);
patches = extract3D_training(conf,A1, A2,conf.scale,{});

size(patches)


if size(patches,2) > 10000
    
    r_idxs = randi(size(patches,2),[10000,1]);
    patches = patches(:,r_idxs);
    features = features(:,r_idxs);
end


% Set KSVD configuration
ksvd_conf.iternum = 20*2; % TBD
% ksvd_conf.iternum = 40; % TBD
ksvd_conf.memusage = 'normal'; % higher usage doesn't fit...
%ksvd_conf.dictsize = 5000; % TBD
ksvd_conf.dictsize = dict_size; % TBD
ksvd_conf.Tdata = 3; % maximal sparsity: TBD
ksvd_conf.samples = size(patches,2);


% PCA dimensionality reduction
C = double(features * features');
[V, D] = eig(C);
D = diag(D); % perform PCA on features matrix 
D = cumsum(D) / sum(D);
k = find(D >= 1e-3, 1); % ignore 0.1% energy
conf.V_pca = V(:, k:end); % choose the largest eigenvectors' projection
conf.ksvd_conf = ksvd_conf;
features_pca = conf.V_pca' * features;


% Combine into one large training set
clear C D V
ksvd_conf.data = double(features_pca);
clear features_pca
% Training process (will take a while)
tic;
fprintf('Training [%d x %d] dictionary on %d vectors using K-SVD\n', ...
    size(ksvd_conf.data, 1), ksvd_conf.dictsize, size(ksvd_conf.data, 2))
[conf.dict_lores, gamma] = ksvd(ksvd_conf); 
toc;



fprintf('Computing high-res. dictionary from low-res. dictionary\n');
% dict_hires = patches / full(gamma); % Takes too much memory...
patches = double(patches); % Since it is saved in single-precision.
dict_hires = (patches * gamma') * inv(full(gamma * gamma'));

conf.dict_hires = double(dict_hires); 

% 

if dict_size < 1024
    lambda = 0.01;
elseif dict_size < 2048
    lambda = 0.1;
elseif dict_size < 8192
    lambda = 1;
else
    lambda = 5;
end


if dict_size < 10000
    conf.ProjM = inv(conf.dict_lores'*conf.dict_lores+lambda*eye(size(conf.dict_lores,2)))*conf.dict_lores';
    conf.PP = (1+lambda)*conf.dict_hires*conf.ProjM;
else
    % here should be an approximation
    conf.PP = zeros(size(conf.dict_hires,1), size(conf.V_pca,2));
    conf.ProjM = [];
end
conf.points = [1:1:size(conf.dict_lores,2)];

conf.pointslo = conf.dict_lores(:,conf.points);
conf.pointsloPCA = conf.pointslo'*conf.V_pca';

conf.PPs = [];
if  size(conf.dict_lores,2) < 40
    clustersz = size(conf.dict_lores,2);
else
    clustersz = 10;
end
D = abs(conf.pointslo'*conf.dict_lores);

for i = 1:length(conf.points)
    [vals idx] = sort(D(i,:), 'descend');
    if (clustersz >= size(conf.dict_lores,2)/2)
        conf.PPs{i} = conf.PP;
    else
        Lo = conf.dict_lores(:, idx(1:clustersz));
        conf.PPs{i} = 1.01*conf.dict_hires(:,idx(1:clustersz))*inv(Lo'*Lo+0.01*eye(size(Lo,2)))*Lo';
    end
end

ANR_PPs = conf.PPs;


end