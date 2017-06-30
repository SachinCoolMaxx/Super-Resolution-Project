function S2_HR = funcCombineAllFourierOrientations_standardMaxo(S2_r, all_axes, p, orig_size, maxdiff)

sz = orig_size;
% % 
% 
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
    axS2 = rotImg3o(currS2_r, -pi/2, curr_axes, 'linear',0);
    axS2(isnan(axS2)) = 0;
  axS2 = axS2(maxdiff+1:maxdiff+sz(1), maxdiff+1:maxdiff+sz(2), maxdiff+1:maxdiff+sz(3));
    
    axialS2{ax_iter} = axS2;
    fftAxialS2{ax_iter} = fftn(axialS2{ax_iter});
  
    
end


finalFFT = zeros(sz);
indices=zeros(sz);
for iter = 1:length(axialS2)
    for i=1:sz(1)
        for j=1:sz(2)
            for k=1:sz(3)
                if abs(fftAxialS2{iter}(i,j,k))>=abs(finalFFT(i,j,k))
                    finalFFT(i,j,k)=fftAxialS2{iter}(i,j,k);
                    indices(i,j,k)=iter;               
                end
            end
        end
    end
end

% code for unpadding according to the image size, u need to add extra 1 to
% both sides, if image is odd sized, becz paddinf is like x+1 and x


S2_HR = ifftn(finalFFT);

end




%---------------------------------------------------------
%my random experimental stuff
% S2_HR = S2_HR(md(1)+2:md(1)+sz(1)+1, md(2)+2:md(2)+sz(2)+1, md(3)+1:md(3)+sz(3));


% indices;
% a=unique(indices);
% out = [a,histc(indices(:),a)]  ;  
% freq=out(:,2)  ;
% figure('Name','Barfor3direcns')
% bar(a,freq)
% sz



 
 
% 
%  clear M
%  colorRGB = [1,0,0;0,1,0;0,0,1];
% 
%  mid=indices(:,:,90);
%  figure('Name','3direcnsFrame90z');
%  RGB = ind2rgb(mid,colorRGB);
%  image(RGB)
%  
%   mid2=indices(:,:,135);
%  figure('Name','3direcnsFrame135z');
%  RGB = ind2rgb(mid2,colorRGB);
%  image(RGB)
%  
%  
%   mid3=indices(:,:,45);
%  figure('Name','3direcnsFrame45z');
%  RGB = ind2rgb(mid3,colorRGB);
%  image(RGB)
%  
%  
%   mid4=indices(:,:,180);
%  figure('Name','3direcnsFrame180z');
%  RGB = ind2rgb(mid4,colorRGB);
%  image(RGB)
%  
%  %
%  clear N
% figure
% clear frontMatrix
% for j=1:180
%     
%     front=indices(j,:,:)
%     for k=1:180
%        frontMatrix(k,:)=front(:,:,k)';
%     end
%     frontMatrix
%     N(j)= ind2rgb(frontMatrix,colorRGB);
%     
% end    
% 
%  movie(N,1,2)
%  v = VideoWriter('/home/sachin/Documents/Results/myavifile3directionsFront.avi');
%  open(v)
%  writeVideo(v,N)
%       
% middle=indices(:,:,90);
% figure('Name','middle')
% RGB = ind2rgb(middle,colorRGB);
% image(RGB);
% 
% top=indices(:,:,1)
% figure('Name','top')
% RGB = ind2rgb(top,colorRGB);
% image(RGB);
% 
% bottom=indices(:,:,180)
% figure('Name','bottom')
% RGB = ind2rgb(bottom,colorRGB);
% image(RGB);
% 
% 
% front=indices(1,:,:)
% for k=1:180
%     frontMat(k,:)=front(:,:,k)';
% end
% frontMat;
% figure('Name','front')
% RGB = ind2rgb(frontMat,colorRGB);
% image(RGB);
% 
% 
% side=indices(:,1,:)
% for k=1:180
%     sideMat(k,:)=side(:,:,k);
% end    
% figure('Name','side')
% RGB = ind2rgb(sideMat,colorRGB);
% image(RGB);

%N(k-2)=im2frame(RGB,colorRGB);


% [x,y,z] = meshgrid(sz);
% scatter3(x(:),y(:),z(:),5,indices(:))


% tmp.img = log10(1 + abs(fftshift(finalFFT)));
% tmp.hdr.dime.dim = [3 size(finalFFT,1) size(finalFFT,2) size(finalFFT,3) 1 1 1 1]
% tmp.hdr.dime.pixdim = [1 0.8 0.8 0.8 0 0 0 0]
% outfile = [results_dir,'final_fft_varPatch_Gradient_Zeyde_spline','_p_',num2str(p),'_iter2.nii'];
% outfile = [results_dir,'final_fft_aa12_flair_5x5x1_Gradient_MyTree_spline','_p_',num2str(p),'_overlap_ScSR_noH_A.nii'];
% tmp.fileprefix = outfile
% % save_nii(tmp, outfile)
% figure;
% imagesc(abs(finalFFT(:,:,90)));

% finalFFT = zeros(256,256,256);
% indices=zeros(256,256,256);
% for iter = 1:length(axialS2)
%     for i=1:256
%         for j=1:256
%             for k=1:256
%                 if abs(fftAxialS2{iter}(i,j,k))>=abs(finalFFT(i,j,k))
%                     finalFFT(i,j,k)=fftAxialS2{iter}(i,j,k);
%                     indices(i,j,k)=iter;               
%                 end
%             end
%         end
%     end
% end






