function S2_HR = funcCombineAllFourierOrientations_standard(S2_r, all_axes, p, orig_size, maxdiff)


sz = orig_size;
% 
% % 
% maxdiff = 58
% sz = [256 256 198]

%%%% - TEMP %%%%
axialS2 = cell(size(S2_r));
fftAxialS2 = cell(size(S2_r));
sumweights = zeros(sz);
for ax_iter = 1:size(all_axes, 1)
    curr_axes = all_axes(ax_iter,:);
    currS2_r = S2_r{ax_iter};
%     currS2_r = currS2_r(maxdiff+1:maxdiff+sz(1), maxdiff+1:maxdiff+sz(2), maxdiff+1:maxdiff+sz(3));
%     tmp.img = currS2_r;
%     tmp.hdr.dime.dim = [3 size(currS2_r,1) size(currS2_r,2) size(currS2_r,3) 1 1 1 1]
%     tmp.hdr.dime.pixdim = [1 0.8 0.8 0.8 0 0 0 0]
%     outfile = [results_dir,'img_3x3x1_Gradient_MyTree_spline_',num2str(curr_axes(1)),'_',num2str(curr_axes(2)),'_',num2str(curr_axes(3)),'.nii'];
%     tmp.fileprefix = outfile
%     save_nii(tmp, outfile)
    
    % rotate by -pi/2 to the curr_axes
    axS2 = rotImg3(currS2_r, -pi/2, curr_axes, 'linear',0);
    axS2(isnan(axS2)) = 0;
    axS2 = axS2(maxdiff+1:maxdiff+sz(1), maxdiff+1:maxdiff+sz(2), maxdiff+1:maxdiff+sz(3));
    
    axialS2{ax_iter} = axS2;
    fftAxialS2{ax_iter} = fftn(axialS2{ax_iter});
  
    
end
parfor iter = 1:length(axialS2)
  sumweights = sumweights + abs(fftAxialS2{iter}).^p; 
%   tmp.img =  axialS2{iter}
%   tmp.hdr.dime.dim = [3 size(axS2,1) size(axS2,2) size(axS2,3) 1 1 1 1]
%   tmp.hdr.dime.pixdim = [1 0.8 0.8 0.8 0 0 0 0]
%   outfile = [results_dir,'axial_img_varPatch_Gradient_MyTree_spline_',num2str(iter),'.nii'];
%   tmp.fileprefix = outfile
%   save_untouch_nii(tmp, outfile)
end

weights = cell(size(axialS2))
for iter = 1:length(axialS2)
    weights{iter} = (abs(fftAxialS2{iter}).^p)./sumweights;
%     tmp.img = log10(1 + abs(fftshift(fftAxialS2{iter})));
%     tmp.hdr.dime.dim = [3 size(fftAxialS2{ax_iter},1) size(fftAxialS2{ax_iter},2) size(fftAxialS2{ax_iter},3) 1 1 1 1]
%     tmp.hdr.dime.pixdim = [1 0.8 0.8 0.8 0 0 0 0]
%     outfile = [results_dir,'fft_varPatch_Gradient_MyTree_spline_',num2str(iter),'_overlap.nii'];
%     tmp.fileprefix = outfile
%     save_untouch_nii(tmp, outfile)
    
%     tmp.img = weights{iter};
%     tmp.hdr.dime.dim = [3 size(fftAxialS2{ax_iter},1) size(fftAxialS2{ax_iter},2) size(fftAxialS2{ax_iter},3) 1 1 1 1]
%     tmp.hdr.dime.pixdim = [1 0.8 0.8 0.8 0 0 0 0]
%     outfile = [results_dir,'weights_varPatch_Gradient_MyTree_spline_',num2str(iter),'_p_',num2str(p),'_iter2.nii'];
%     tmp.fileprefix = outfile
%     save_untouch_nii(tmp, outfile)
    
    
end
finalFFT = zeros(sz);
for iter = 1:length(axialS2)
   finalFFT = finalFFT + weights{iter}.*fftAxialS2{iter};
    
end

% tmp.img = log10(1 + abs(fftshift(finalFFT)));
% tmp.hdr.dime.dim = [3 size(finalFFT,1) size(finalFFT,2) size(finalFFT,3) 1 1 1 1]
% tmp.hdr.dime.pixdim = [1 0.8 0.8 0.8 0 0 0 0]
% outfile = [results_dir,'final_fft_varPatch_Gradient_Zeyde_spline','_p_',num2str(p),'_iter2.nii'];
% outfile = [results_dir,'final_fft_aa12_flair_5x5x1_Gradient_MyTree_spline','_p_',num2str(p),'_overlap_ScSR_noH_A.nii'];
% tmp.fileprefix = outfile
% save_nii(tmp, outfile)

S2_HR = ifftn(finalFFT);

end
