function [ imz ] = RotAboutz( img ,teta )

% [nd1, nd2, nd3] = size(img);
run='rotaboutz'
% 
% % nx = nd1
% % ny = nd3
% % nz = nd2
% 
% 
% % for k=1:nx
% %     for j=1:nd3
% %         imx(k,ny+1-j,:) = imagepad(k,:,j)'; 
% %     end
% % end
% run='z'
% imz=zeros(nd1,nd2,nd3);
% if teta==pi/2
%     for k=1:nd3
%         %
%         for j=20:200              
%             for i=1:181
%                  imz(181+1-(j-19),i+19,k) = img(i,j,k);
%             end
%         end
%     end
%     
% elseif teta == -pi/2
%     for k=1:nd3
%         %
%         for j=20:200              
%             for i=1:181
%                  imz(j-19,220-(i+19),k) = img(i,j,k);
%             end
%         end
%     end
% end
% 
% size(imz)

sz = size(img)
% md = floor([256-sz(1),256-sz(2),256-sz(3)]./2);     
run='rotaboutz'
% imgz = PadImage(img,256);
imgz=img;
if teta == pi/2
    for k=1:256
        for j=1:256
            for i=1:256
                imz(256+1-j,i,k) = imgz(i,j,k);
            end
        end
    end
    
elseif teta == -pi/2
    for k=1:256
        for j=1:256
            for i=1:256
                imz(j,256+1-i,k) = imgz(i,j,k);
            end
        end
    end
   
end

% imz = imz(md(1)+1:md(1)+sz(1),md(2)+1:md(2)+sz(2),md(3)+1:md(3)+sz(3));
sizeOftheImageAfterRec = size(imz)
end





