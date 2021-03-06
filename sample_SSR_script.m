% this script will generate multiple rotations of the subject image and
% synthesize a high resolution version of each and then combine all
% together using sapiro's method. this experiment uses 3d patches.

%addpath('/home/amod/mimecs/code/code_basic_libs/MyTree/code/');
addpath('/home/amod/mimecs/code/code_basic_libs/');
%addpath('./make_icosahedron/');
addpath(genpath('/home/amod/mimecs/code/code_basic_libs/ScSR/'));
addpath(genpath('/home/roy/Programs/spams/'))
addpath(genpath('/home/amod/mimecs/code/code_basic_libs/AplusCodes_SR/'))
addpath(genpath('/home/amod/mimecs/ISBI2013_RegressionSynthesis/code/metrix_mux/'));

% load input images

addpath(genpath('/home/amod/mimecs/code/code_basic_libs/sepspyr/deps/matlabPyrTools-1.3/'))


data_dir = '/iacl/pg15/amod/SelfSuperRes_Experiments/data/MS_Data/processed/FLAIR2p2mm/'
data_list = dir([data_dir,'*wmp.nii']);

results_dir = '/bmap/pg15/amod/SelfSuperRes_Experiments/results/MS_Data/FLAIR2p2mm/'
temp_results_dir = '/bmap/pg15/amod/SelfSuperRes_Experiments/results/MS_Data/FLAIR2p2mm/temp/'
% truth_dir =  '/iacl/pg15/amod/SelfSuperRes_Experiments/data/NeuroMM/imgs/cropped/';
% true_list = dir([truth_dir,'*.nii']);
% interp_psnr = zeros(20,1);
% all_psnr = zeros(3,20,9)


% tmp = load_nii([data_dir,'3dt2_0p8mm_brain_axial.nii']);
% I = double(tmp.img);

% first, load the atlas images. A1 = downsampled in X + downsampled in Z
% image. A2 = just downsampled in Z. These are both axial images
for iter = 1:length(data_list)
    tmp = load_nii([data_dir, data_list(iter).name]);
    subjname = strtok(data_list(iter).name, '_');
    A2 = double(tmp.img);
    A2 = wmPeakNormalize(A2);
    A2(A2<0) = 0;
    dim_A2_orig = size(A2);
    
    
    voxel_size_z = 2.2;
    voxel_size_x = 0.828;
    voxel_size_y = 0.828;
    
    HRvoxel_size_z = 0.828;
    HRvoxel_size_x = 0.828;
    HRvoxel_size_y = 0.828;
    
    A2_X_spacing = voxel_size_x;
    A2_Y_spacing = voxel_size_y;
    A2_Z_spacing = voxel_size_z;
    [A2Ygrid,A2Xgrid, A2Zgrid] = meshgrid(0+0.5*A2_Y_spacing:A2_Y_spacing:(size(A2,1)-0.5)*A2_Y_spacing, 0+0.5*A2_X_spacing:A2_X_spacing:(size(A2,2)-0.5)*A2_X_spacing, 0+0.5*A2_Z_spacing:A2_Z_spacing:(size(A2,3)-0.5)*A2_Z_spacing) ;
    [newA2Ygrid,newA2Xgrid, newA2Zgrid] = meshgrid(0+0.5*HRvoxel_size_y:HRvoxel_size_y:(size(A2,1)-0.5)*A2_Y_spacing, 0+0.5*HRvoxel_size_x:HRvoxel_size_x:(size(A2,2)-0.5)*A2_X_spacing, 0+0.5*HRvoxel_size_z:HRvoxel_size_z:(size(A2,3))*A2_Z_spacing-0.5*HRvoxel_size_z) ;
    
    
    
    ZI = interp3(A2Ygrid, A2Xgrid, A2Zgrid, permute(A2,[2,1,3]), newA2Ygrid, newA2Xgrid, newA2Zgrid, 'spline');
    A2up = permute(ZI,[2,1,3]);
    A2up = wmPeakNormalize(A2up);
    % PAD ALL IMAGES TO MAKE THEM BIG.
    sz = size(A2up);
    maxdiff = max(sz) - min(sz);
    A2up = padarray(A2up,[maxdiff,maxdiff,maxdiff]);
    A2up(isnan(A2up)) = 0;
    
    
    Ny = size(A2up,1);
    Nx = size(A2up,2);
    Nz = size(A2up,3);
    
    lp_Nyq_freq_z = 1/(2*voxel_size_z);
    lp_Nyq_freq_y = 1/(2*voxel_size_y);
    lp_Nyq_freq_x = 1/(2*voxel_size_x);
    
    Nyq_u = 1/(2*HRvoxel_size_x); % Nyquist of data in first dimension
    Nyq_v = 1/(2*HRvoxel_size_y); % Nyquist of data in second dimension
    Nyq_w = 1/(2*HRvoxel_size_z); % Nyquist of data in second dimension
    
    
    UGridSize = Nx*HRvoxel_size_x;
    VGridSize = Ny*HRvoxel_size_y;
    WGridSize = Nz*HRvoxel_size_z;
    
    du = 1/UGridSize;   % Wavenumber increment
    dv = 1/VGridSize;   % Frequency increment
    dw = 1/WGridSize;
    
    U = -Nyq_u : du : Nyq_u-du;
    V = -Nyq_v : dv : Nyq_v-dv;
    W = -Nyq_w : dw : Nyq_w-dw;
    cutoff_u = (length(U)/2)*(lp_Nyq_freq_x / Nyq_u);
    cutoff_v = (length(V)/2)*(lp_Nyq_freq_y / Nyq_v);
    cutoff_w = (length(W)/2)*(lp_Nyq_freq_z / Nyq_w);
    
    H_A = lpfilter_3d('truncSincZ', length(U), length(V),length(W), 0,1,cutoff_u,cutoff_v, cutoff_w);
    H_A = fftshift(H_A);
    F_A2up = fftshift(fftn(A2up));
    
    all_axes = [1,0,0;
        0,1,0;
        1,1,0;
        -1,1,0;
        1,1,1;
        1,1,-1;
        -1,1,1;
        -1,1,-1
        0,0,1];
    %
    %     freq = 2;
    %     sphere = 0;
    %     make_plot = 1;
    %     faceopaque = 0;
    %     [ ax,ay, az, TRI] = make_icosahedron(freq, 1, sphere, make_plot, faceopaque);
    %
    %     all_axes  = calcFaceNormal(ax, ay, az, TRI);
    %     all_axes = [all_axes;1,0,0;0,1,0;0,0,1];
    
    angles = [pi/2];%, pi/3, pi/2];
    upsample_factor = [voxel_size_y/HRvoxel_size_y,voxel_size_x/HRvoxel_size_x,voxel_size_z/HRvoxel_size_z];
    windows = {[2,2,2],[3,3,3],[4,4,4]};% min size of patch
    orig_blur_dir = [0,0,1];
    for witer = 1%:3
        window = windows{witer}
        S2_r = cell(size(all_axes,1), 1);
        parfor ax_iter = 1:size(all_axes,1)
            curr_axes = all_axes(ax_iter,:);
            
                curr_angle = angles(1);
                
                
                ratM = rotationmat3D(curr_angle, curr_axes);
                
                
                
                % find x and y projections
                unit_axes = curr_axes./sqrt(sum((curr_axes.^2)));
                rot_blur_dir = abs(ratM*orig_blur_dir');
                rot_blur_magn = upsample_factor(3)*rot_blur_dir;
                nz_idxs  = find(rot_blur_magn > 0.01);
                onez = ones(size(rot_blur_magn));
                onez(nz_idxs) = 0;
                upscaling = rot_blur_magn + [1;1;1];
                upscaling = round(upscaling');
                
                upscaling
                
                
                
                A2_r = rotImg3(A2up, curr_angle, curr_axes, 'spline',0);
                A2_r(isnan(A2_r)) = 0;
                F_A2_r = fftshift(fftn(A2_r));
                
                % now apply this filter to A2_r to create S1
                F_A2_r_H_a = F_A2_r;%.*H_A;
                F_A2_r_H_a(isnan(F_A2_r_H_a)) = 0;
                A2_r_H_a = real((ifftn(ifftshift(F_A2_r_H_a))));
                S1_r = A2_r_H_a;
                S1_r = wmPeakNormalize(S1_r);
                
                % Now to create A1. This is done by applying H_r to A2up.
                %         H_r = rotImg3(H_A, curr_angle, curr_axes, 'linear',0);
                %         H_r(isnan(H_r)) = 0;
                h_A = fftshift(ifftn(ifftshift(H_A)));
                h_A_r = rotImg3(h_A, curr_angle, curr_axes, 'linear',0);
                h_A_r(isnan(h_A_r)) = 0;
                H_r = abs(fftshift(fftn(h_A_r)));
                % somehow make H_r to be the size of F_A2up...HERE
                
                
                F_A2_H_r = F_A2up.*H_r;
                A2_H_r = real((ifftn(ifftshift(F_A2_H_r))));
                A1 = A2_H_r;
                A1 = wmPeakNormalize(A1);
                
                % now we have A1, A2, S1. Pull off a synthesis.
                A1(isnan(A1)) = 0;
                
                S1_r(isnan(S1_r)) = 0;
                disp('synthesizing for axis')
                curr_axes
                
                lambda = 0.10
                dict_size = 128;
                overlap = [2,2,2];%[N1(1)-2,N1(2)-2,N1(3)-1];
                sampling = [1,1,1];
                N1 = upscaling.*window;
                border = [15, 15, 15];
                [conf] = funcLearnZeyDeDict(A1, A2up, upscaling, N1,N1, dict_size, border);
                
                conf.border = [1,1,1];
                conf.overlap = conf.window - [1 1 1];
                conf
                S2_r5 = funcSelfSuperRes_Zeyde_gradients(conf, S1_r);
                S2_r5 = wmPeakNormalize(S2_r5);
                
                
                
                S2_r{ax_iter} = (S2_r5);% + S2_r5)/2;
                tmpp = tmp;%tmparray{ax_iter};
                tmpp.img = S2_r5; %S2_r{ax_iter, ang_iter};
                tmpp.hdr.dime.dim = [3 size(S2_r5,1) size(S2_r5,2) size(S2_r5,3) 1 1 1 1]
                tmpp.hdr.dime.pixdim = [1 HRvoxel_size_y HRvoxel_size_x HRvoxel_size_z 0 0 0 0]
                outfile = [temp_results_dir,subjname,'_varPatchSize_Zeyde_Gradient_truncSinc_',num2str(ax_iter),'_window_',num2str(conf.window),'_',num2str(dict_size),'_unsharp0p0_noHA_9dir.nii'];
                %          outfile = ['test_5x5x1_Gradient_MyTree_spline_',num2str(ax_iter),'_overlap_',num2str(overlap),'_9directions','.nii'];
                tmpp.fileprefix = outfile
                if ax_iter == 9
                    
                    save_nii(tmpp, outfile);
                end
                    
                
 
        end
%         
        for p = [0.05, 0.1, 0.2, 0.8, 1.2, 2,4,8,12]
            S2_HR = funcCombineAllFourierOrientations_standard(S2_r, all_axes, p,sz,maxdiff);
            S2_HR = wmPeakNormalize(S2_HR);
            
  
            tmp.img = S2_HR;
            tmp.hdr.dime.dim = [3 size(S2_HR,1) size(S2_HR,2) size(S2_HR,3) 1 1 1 1]
            tmp.hdr.dime.pixdim = [1 HRvoxel_size_y HRvoxel_size_x HRvoxel_size_z 0 0 0 0]
            outfile = [results_dir,subjname, '_varPatchsize_Zeyde_truncSinc_FinalSuperRes_p',num2str(p),'_3dpatches_',num2str(window),'sharp_0p0_noHA_9dir.nii'];
            
            tmp.fileprefix = outfile
            save_nii(tmp, outfile)
        end
        
    end
    
    
end
% 




