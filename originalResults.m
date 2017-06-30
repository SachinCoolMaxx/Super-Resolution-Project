% this script will generate multiple rotations of the subject image and
% synthesize a high resolution version of each and then combine all
% together using sapiro's method. this experiment uses 3d patches.

addpath('/home/amod/mimecs/code/code_basic_libs/');
addpath('/home/roy/Programs/fileformats_toolbox/nii/');
addpath(genpath('/home/amod/mimecs/code/code_basic_libs/ScSR/'));
addpath(genpath('/home/roy/Programs/spams/'))
addpath(genpath('/iacl/pg17/sachin/for_sachin/code/AplusCodes_SR/'))
addpath(genpath('/home/amod/mimecs/ISBI2013_RegressionSynthesis/code/metrix_mux/'));
addpath(genpath('/home/amod/mimecs/code/code_basic_libs/sepspyr/deps/matlabPyrTools-1.3/'))
addpath(genpath('/iacl/pg17/sachin/for_sachin/code_backup/code'))

% input dirdata_dir = '/iacl/pg15/amod/SelfSuperRes_Experiments/data/NeuroMM/imgs/cropped/2mm/';
data_dir = '/iacl/pg17/sachin/for_sachin/data/NeuroMM/2mmNew/'
data_list = dir([data_dir,'*.nii']);

% store results
results_dir = '/iacl/pg17/sachin/for_sachin/results/NeuroMM/2mmNew/'
temp_results_dir = '/iacl/pg17/sachin/for_sachin/results/NeuroMM/2mmNew/temp/'

% comparison with ground truth
truth_dir =  '/iacl/pg17/sachin/for_sachin/data/NeuroMM/grndtruth/'
true_list = dir([truth_dir,'*.nii']);
interp_psnr = zeros(20,1);
%all_psnr = zeros(3,20,9)



% first, load the atlas images. A1 = downsampled in X + downsampled in Z
% image. A2 = just downsampled in Z. These are both axial images

for iter = 1: length(data_list)
    
    tmp = load_untouch_nii([data_dir, data_list(iter).name]);
    subjname = strtok(data_list(iter).name, '_');
    A2 = double(tmp.img);
    A2(A2<0) = 0;
    int_thresh = 10;
    A2 = wmPeakNormalize(A2, int_thresh); % wm peak normalize for PSNR calculations later
    
    dim_A2_orig = size(A2);
    
    
    voxel_size_z = tmp.hdr.dime.pixdim(4);
    voxel_size_x = tmp.hdr.dime.pixdim(3);
    voxel_size_y = tmp.hdr.dime.pixdim(2);
    
    HRvoxel_size_z = min([voxel_size_y, voxel_size_x, voxel_size_z]);
    HRvoxel_size_x = HRvoxel_size_z;
    HRvoxel_size_y = HRvoxel_size_z;
    
    % create the grid
    A2_X_spacing = voxel_size_x;
    A2_Y_spacing = voxel_size_y;
    A2_Z_spacing = voxel_size_z;
    [A2Ygrid,A2Xgrid, A2Zgrid] = meshgrid(0+0.5*A2_Y_spacing:A2_Y_spacing:(size(A2,1)-0.5)*A2_Y_spacing, 0+0.5*A2_X_spacing:A2_X_spacing:(size(A2,2)-0.5)*A2_X_spacing, 0+0.5*A2_Z_spacing:A2_Z_spacing:(size(A2,3)-0.5)*A2_Z_spacing) ;
    [newA2Ygrid,newA2Xgrid, newA2Zgrid] = meshgrid(0+0.5*HRvoxel_size_y:HRvoxel_size_y:(size(A2,1)-0.5)*A2_Y_spacing, 0+0.5*HRvoxel_size_x:HRvoxel_size_x:(size(A2,2)-0.5)*A2_X_spacing, 0+0.5*HRvoxel_size_z:HRvoxel_size_z:(size(A2,3))*A2_Z_spacing-0.5*HRvoxel_size_z) ;
    
    
    % recalculate A2 on this grid
    ZI = interp3(A2Ygrid, A2Xgrid, A2Zgrid, permute(A2,[2,1,3]), newA2Ygrid, newA2Xgrid, newA2Zgrid, 'spline');
    A2up = permute(ZI,[2,1,3]);
    A2up = wmPeakNormalize(A2up, int_thresh);
    
    % pad the image so that we don't run out of x,y,z extents after
    % rotation
    sz = size(A2up)
    maxdiff = max(sz) - min(sz);
%       md = floor([256-sz(1),256-sz(2),256-sz(3)]./2)
    A2up = padarray(A2up,[maxdiff,maxdiff,maxdiff]);
%     A2up = PadImage(A2up,256);
    A2up(isnan(A2up)) = 0;
    
    
    % model the initial blurring as a truncSinc in Fourier space
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
    
    % create the blurring filter in Fourier space
    H_A = lpfilter_3d('truncSincZ', length(U), length(V),length(W), 0,1,cutoff_u,cutoff_v, cutoff_w);
    H_A = fftshift(H_A);
    F_A2up = fftshift(fftn(A2up));
    
    % directions of different blurs that we will consider
    
    all_axes = [1,0,0;
        0,1,0;
        1,1,0;
        -1,1,0;
        1,1,1;
        1,1,-1;
        -1,1,1;
        -1,1,-1
        0,0,1];
        
    angles = [pi/2];
    
    % upsample_factor only used as a parameter for patch size while extracting features
    upsample_factor = [voxel_size_y/HRvoxel_size_y,voxel_size_x/HRvoxel_size_x,voxel_size_z/HRvoxel_size_z];
    windows = {[2,2,2],[3,3,3],[4,4,4]};% min size of patch -- just pick one. [2,2,2] works
    
    % original blur direction is z
    orig_blur_dir = [0,0,1];
    
    for witer = 1%:3
        window = windows{witer}
        
        % store the synthesized images for each direction in a cell array
        S2_r = cell(size(all_axes,1), 1);
        
        % parallelize synthesis for directions. remove par for debugging.
        for ax_iter = 1:size(all_axes,1)
            curr_axes = all_axes(ax_iter,:);
            
            curr_angle = angles(1);
            
            % create rotation matrix for rotation of 90 degrees about
            % curr_axes
            ratM = rotationmat3D(curr_angle, curr_axes);
            
            
            
            % find x and y projections: for deciding the x and y
            % dimensions of the patch feature we are going to use. the
            % less the blur the longer the patch dimension
            unit_axes = curr_axes./sqrt(sum((curr_axes.^2)));
            rot_blur_dir = abs(ratM*orig_blur_dir');
            rot_blur_magn = upsample_factor(3)*rot_blur_dir;
            nz_idxs  = find(rot_blur_magn > 0.01);
            onez = ones(size(rot_blur_magn));
            onez(nz_idxs) = 0;
            upscaling = rot_blur_magn + [1;1;1];
            upscaling = round(upscaling');
            
            upscaling
            
            
            % this can be very memory inefficient. if memory is falling
            % short, use 'linear' instead of 'spline'
            
            % A2_r = rotImg3o(A2up, curr_angle, curr_axes, 'spline',0);
            A2_r = rotImg3o(A2up, curr_angle, curr_axes, 'linear',0);
            
            A2_r(isnan(A2_r)) = 0;
            F_A2_r = fftshift(fftn(A2_r));
            
            % now apply this filter to A2_r to create S1.
            % Recent experiments have shown that we can get better results
            % PSNR-wise if H_A is NOT applied. THIS IS NOT INTUITIVE AND
            % NEEDS TO BE INVESTIGATED
            F_A2_r_H_a = F_A2_r;%.*H_A;
            
            F_A2_r_H_a(isnan(F_A2_r_H_a)) = 0;
            A2_r_H_a = real((ifftn(ifftshift(F_A2_r_H_a))));
            S1_r = A2_r_H_a;
            S1_r = wmPeakNormalize(S1_r, int_thresh);
            
            % Now to create A1. This is done by applying H_r to A2up.
            % h_A is the filter in image space, which is obtained by
            % inverse Fourier transforming of H_A. h_A_r is rotated filter
            % in the images space, H_r is the rotated filter in Fourier
            % space
            h_A = fftshift(ifftn(ifftshift(H_A)));
            h_A_r = rotImg3o(h_A, curr_angle, curr_axes, 'linear',0); % can be spline if memory permits
            h_A_r(isnan(h_A_r)) = 0;
            H_r = abs(fftshift(fftn(h_A_r)));
            
            % Create A1 by blurring A2up by H_r
            F_A2_H_r = F_A2up.*H_r;
            A2_H_r = real((ifftn(ifftshift(F_A2_H_r))));
            A1 = A2_H_r;
            A1 = wmPeakNormalize(A1, int_thresh);
            
            % now we have A1, A2, S1. Pull off a synthesis.
            A1(isnan(A1)) = 0;
            
            S1_r(isnan(S1_r)) = 0;
            disp('synthesizing for axis')
            curr_axes
            
            % parameters for Anchored neighborhood regression
            lambda = 0.10
            
            % dict_size: testing has showed that this does not matter much
            % for 128, 256, 512, 1024. Using 128 because it is fast
            
            dict_size = 128;
            
            % patches are synthesized
            overlap = [2,2,2];%[N1(1)-2,N1(2)-2,N1(3)-1];
            
            % sampling patches at this interval for training
            sampling = [1,1,1];
            
            % patch size
            N1 = upscaling.*window;
            border = [15, 15, 15];
            
            % can be increased when there are lots of interpolation
            % intensities
            intensity_threshold = 30;
            % learn the mapping using  funcLearnANRDict
            [conf] = funcLearnANRDict(A1, A2up, upscaling, N1, intensity_threshold, dict_size, border);
            
            % change the border and overlap parameters for synthesis
            conf.border = [1,1,1];
            conf.overlap = conf.window - [1 1 1];
            
            conf
            
            
            
            S2_r5 = funcSelfSuperRes_ANR_gradients(conf, S1_r);
            S2_r5 = wmPeakNormalize(S2_r5, int_thresh);
            
            
            S2_r{ax_iter} = (S2_r5);
            tmpp = tmp;%tmparray{ax_iter};
            tmpp.img = S2_r5; %S2_r{ax_iter, ang_iter};
            tmpp.hdr.dime.dim = [3 size(S2_r5,1) size(S2_r5,2) size(S2_r5,3) 1 1 1 1]
            tmpp.hdr.dime.pixdim = [1 HRvoxel_size_y HRvoxel_size_x HRvoxel_size_z 0 0 0 0]
            outfile = [temp_results_dir,subjname,'_varPatchSize_ANR_Gradient_truncSinc_',num2str(ax_iter),'_window_',num2str(conf.window),'_',num2str(dict_size),'_9dir_orderchang.nii'];
            %          outfile = ['test_5x5x1_Gradient_MyTree_spline_',num2str(ax_iter),'_overlap_',num2str(overlap),'_9directions','.nii'];
            tmpp.fileprefix = outfile
            save_untouch_nii(tmpp, outfile);
            
            
        end
    
        %%
      
      %  for p = [0.05, 0.1, 0.2, 0.8, 1.2, 2,4,8,12]
       for p = [35]
            S2_HR = funcCombineAllFourierOrientations_standardMaxo(S2_r, all_axes, p,sz,maxdiff);
            S2_HR = wmPeakNormalize(S2_HR, int_thresh);
            
            tmp.img = S2_HR;
            tmp.hdr.dime.dim = [3 size(S2_HR,1) size(S2_HR,2) size(S2_HR,3) 1 1 1 1];
            tmp.hdr.dime.pixdim = [1 HRvoxel_size_y HRvoxel_size_x HRvoxel_size_z 0 0 0 0];
            outfile = [results_dir,subjname, 'pMaximum(p=infi)',num2str(p),'_3dpatches_withHROldnewrotation999',num2str(window),'_9dirOrderChanged_.nii'];
            
            tmp.fileprefix = outfile;
            save_untouch_nii(tmp, outfile);
            
           
            tmptruth = load_untouch_nii([truth_dir,true_list(iter).name]);
            A_HR = double(tmptruth.img);
            A_HR = wmPeakNormalize(A_HR, int_thresh);
            fg_idxs = find(A_HR > 375); % removes backgr    ound and very dim voxels
           
            maxI = max(A_HR(:));
           % A2up = A2up(maxdiff+1:maxdiff+sz(1), maxdiff+1:maxdiff+sz(2), maxdiff+1:maxdiff+sz(3));
            mse = mean((A_HR(fg_idxs) - S2_HR(fg_idxs)).^2)
            psnr = 20*log10(maxI) - 10*log10(mse)
            iter
            interp_psnr(iter) = psnr
%             k=k+1;
%             pstr = {'0p05', '0p1', '0p2', '0p8', '1p2', '2', '4','8','12'};
            
       end
    end
end
%


