function S2hat = funcSelfSuperRes_ANR_gradients(conf, S1)

% count = 1;
% interpolated = cell(1,1);
% for z = 1: size(S1,3)
%     interpolated{z} = S1(:,:,z);
%     count = count + 1;
% end


threshold = conf.intensity_threshold;
features = extract3D(conf, S1, conf.scale, conf.filters);
features = double(features);

features = conf.V_pca'*features;
blocksize = 50000; %if not sufficient memory then you can reduce the blocksize


D = abs(conf.pointslo'*features);
[val idx] = max(D);
parfor l = 1:size(features,2)
    patches(:,l) = conf.PPs{idx(l)} * features(:,l);
end



img_size = size(S1);
%
grid = sampling_grid3D(S1,img_size,conf.intensity_threshold, ...
    conf.window, conf.overlap, conf.border, conf.scale);

result = overlap_add3D(patches, img_size, grid);

result = result + S1;
S2hat = result;
fprintf('.');



end
