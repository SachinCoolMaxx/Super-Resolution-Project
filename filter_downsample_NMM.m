addpath('/home/amod/mimecs/code/code_basic_libs/')
data_dir = '/bmap/pg15/amod/SelfSuperRes_Experiments/data/NeuroMM/imgs/'
label_dir = '/bmap/pg15/amod/SelfSuperRes_Experiments/data/NeuroMM/labels/'
outlabel_dir = '/bmap/pg15/amod/SelfSuperRes_Experiments/data/NeuroMM/labels/cropped/';
outdata_dir = '/bmap/pg15/amod/SelfSuperRes_Experiments/data/NeuroMM/imgs/cropped/'
data_list = dir([data_dir,'*.nii']);
label_list = dir([label_dir,'*.nii']);

for iter = 1:length(data_list)
    % load the image, crop it, downsample it to 3 different slice
    % thicknesses 2mm, 3mm, 5mm. only crop the labels, nothing else.
    
    tmpimg = load_nii([data_dir,data_list(iter).name]);
    img = double(tmpimg.img);
    tmplabel = load_nii([label_dir,label_list(iter).name]);  
    labelfilename = strtok(label_list(iter).name, '.');
    subj_name = strtok(labelfilename, 'g');
    labelimg = double(tmplabel.img);
    
    [imgcrop, cropParams] = Crop3D(img,0,15,0); % make a padding of 15 pixels in all directions hopefully. pad with 0.
    labelcrop = labelimg(cropParams.boundingBox(1):cropParams.boundingBox(2),cropParams.boundingBox(3):cropParams.boundingBox(4), cropParams.boundingBox(5):cropParams.boundingBox(6)); 
    
    mod30 = mod(size(imgcrop,3), 30); % why 30? well our experiments are going to be for super-res of 2mm, 3mm, 5mm, slice thickness. So something that is their LCM.
    postPadLength = 30 - mod30;
    
    imgcrop = padarray(imgcrop,[0,0,postPadLength],0, 'post');
    labelcrop = padarray(labelcrop,[0,0,postPadLength],0, 'post');
    
    % save cropped labels
    tmplabel.img = labelcrop;
    tmplabel.hdr.dime.dim = [3 size(labelcrop,1) size(labelcrop,2) size(labelcrop,3) 1 1 1 1]
    tmplabel.hdr.dime.pixdim = [1 1 1 1 0 0 0 0]
    outfile = [outlabel_dir,labelfilename,'_cropped.nii'];
    tmplabel.fileprefix = outfile
    save_nii(tmplabel, outfile)
    
    % save cropped image
    tmpimg.img = imgcrop;
    tmpimg.hdr.dime.dim = [3 size(imgcrop,1) size(imgcrop,2) size(imgcrop,3) 1 1 1 1]
    tmpimg.hdr.dime.pixdim = [1 1 1 1 0 0 0 0]
    outfile = [outdata_dir,subj_name,'_NMM_cropped.nii'];
    tmpimg.fileprefix = outfile
    save_nii(tmpimg, outfile)
    
end

%%
% for the cropped images; create downsampled images of 2mm, 3mm, 5mm slice
% thicknesses
outdata_list = dir([outdata_dir,'*.nii']);
data_5mmdir = '/bmap/pg15/amod/SelfSuperRes_Experiments/data/NeuroMM/imgs/cropped/5mm/'
for iter  = 1:length(outdata_list)
    tmpimg = load_nii([outdata_dir,outdata_list(iter).name]);
    filename = strtok(outdata_list(iter).name, '.');
    
    A = double(tmpimg.img);
    
    voxel_size_z = 5;
    voxel_size_x = 1;
    voxel_size_y = 1;
    
    HRvoxel_size_z = 1;
    HRvoxel_size_x = 1;
    HRvoxel_size_y = 1;
    
    
    Ny = size(A,1);
    Nx = size(A,2);
    Nz = size(A,3);
    
    
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
    F_A2 = fftshift(fftn(A));
    
    
    F_A2h = F_A2.*H_A;
    F_A2h(isnan(F_A2)) = 0;
    A2h = real((ifftn(ifftshift(F_A2h))));
%     outfile = [data_2mmdir, ];
%     tmp.fileprefix = outfile;
%     tmp.img = A2h;
%     save_nii(tmp, outfile);

    % now downsample
    A2_X_spacing = 1;
    A2_Y_spacing = 1;
    A2_Z_spacing = 5;
    A2h_ds = zeros(size(A,1), size(A,2), size(A,3)/A2_Z_spacing);
    [A2Ygrid,A2Xgrid, A2Zgrid] = meshgrid(0+0.5*HRvoxel_size_y :HRvoxel_size_y:(size(A,1)-0.5)*HRvoxel_size_y, 0+0.5*HRvoxel_size_x:HRvoxel_size_x:(size(A,2)-0.5)*HRvoxel_size_x, 0+0.5*HRvoxel_size_z:HRvoxel_size_z:(size(A,3)-0.5)*HRvoxel_size_z) ;
    [LRA2Ygrid,LRA2Xgrid, LRA2Zgrid] = meshgrid(0+0.5*voxel_size_y:voxel_size_y:(size(A2h_ds,1)-0.5)*voxel_size_y, 0+0.5*voxel_size_x:voxel_size_x:(size(A2h_ds,2)-0.5)*voxel_size_x, 0+0.5*voxel_size_z:voxel_size_z:(size(A2h_ds,3))*voxel_size_z-0.5*voxel_size_z) ;

    ZI = interp3(A2Ygrid, A2Xgrid, A2Zgrid, permute(A2h,[2,1,3]), LRA2Ygrid, LRA2Xgrid, LRA2Zgrid, 'spline');
    A2h_ds = permute(ZI,[2,1,3]);
    
    
    
    output_dir = data_5mmdir; 
    outfile = [output_dir, filename,'_5mmZ.nii'];
    tmplr = tmpimg;
    tmplr.hdr.dime.dim = [3 size(A2h_ds,1) size(A2h_ds,2) size(A2h_ds,3) 1 1 1 1]
    tmplr.hdr.dime.pixdim = [1 voxel_size_y voxel_size_x voxel_size_z 0 0 0 0]
    tmplr.fileprefix = outfile;
    tmplr.img = A2h_ds;
    save_nii(tmplr, outfile);


end

