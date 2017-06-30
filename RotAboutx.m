function [ imx ] = RotAboutx( img ,teta )

% [nd1, nd2, nd3] = size(img);
% run='rotaboutx'
% % nx = nd1
% % ny = nd3
% % nz = nd2
% 
% 
% imx=zeros(size(img));
% if teta==pi/2
%     for k=1:nd3
%         for j=21:200
% % since img is 181 219 180, for y to z swapping, we will neglect 39 slices
% % of y direcn, 20 on one side 19 on the other, so from 21 to 200 ->1 to 180
%              imx(:,nd2+1-k-20,j-20) = (img(:,j,k));
%         end
%     end
%     
% elseif teta == -pi/2
%     for k=1:nd3
%         for j=21:200
%             imx(:,k+20,180+1-(j-20)) = (img(:,j,k));
%         end
%     end
% end
% 
% imx = reshape(imx,size(img));



%-------------------------------------------

sz = size(img)
%md = floor([256-sz(1),256-sz(2),256-sz(3)]./2)+[1,1,1]; 257 257 258
% md = floor([256-sz(1),256-sz(2),256-sz(3)]./2);
run='rotaboutx'
% imgx=padarray(img,md);
% imgx = PadImage(img,256); 

% figure('Name','beforerotatn');
% imagesc(imgx(:,:,128))

imgx=img;

if teta == pi/2
    for k=1:256
        for j=1:256
            imx(:,k,256+1-j) = imgx(:,j,k);
        end
    end                                      
    
elseif teta == -pi/2
    for k=1:256
        for j=1:256
            imx(:,256+1-k,j) = imgx(:,j,k);
        end
    end
end
    
%    figure('Name','afterrotatn');
%    imagesc(imx(:,:,256))
%    
  
% imx = imx(md(1)+1:md(1)+sz(1),md(2)+1:md(2)+sz(2),md(3)+1:md(3)+sz(3));
 sizeOftheImageAfterRec = size(imx)
  
% for k=1:sz(3)
%         test34(:,k) = imx(:,110,k)';
% end
     
%    figure('Name','afterrotatnSameSlice');
%    imagesc(test34);

end


%rec = zeros(nd3,nd2)
% %imx = zeros(nx,ny,nz);
% for k=1:nd1
%     for j=1:nd3
%         rec(j,:) = imagepad(k,:,j);
%     end
%     imx(k,:,:) = rot90(rec);
% end

    