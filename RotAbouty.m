function [ imy ] = RotAbouty( img ,teta )
% 
% [nd1, nd2, nd3] = size(img);
% run='roty'
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
% imy=zeros(nd1,nd2,nd3);
% if teta==pi/2
%     for k=1:nd3
%         %becz img is 181 219 180, u will have to neglect one row of x
%         for j=2:181              
%              imy(k,:,180+2-j) = (img(j,:,k));
%         end
%     end
%     
% elseif teta == -pi/2
%     for k=1:nd3
%         %becz img is 180 219 180, u will have to neglect one row of x
%         for j=1:180 
%         %since i earlier also neglected x=181 , hence here also x=181
%         %neglected instead of x=1
%              imy(181+0-k,:,j) = (img(j,:,k));
%         end
%     end
% end
% 
% size(imy);

sz = size(img)
% md = floor([256-sz(1),256-sz(2),256-sz(3)]./2);     
run='rotabouty'
% imgy = PadImage(img,256);
% sizeafterpaddn = size(imgy)


% md = floor([256-sz(1)+1,257-sz(2),256-sz(3)]./2)+ [1,1,1];

imgy=img;
if teta == pi/2
    for k=1:256
        for j=1:256
            imy(k,:,256+1-j) = imgy(j,:,k);
        end
    end
    
elseif teta == -pi/2
    for k=1:256
        for j=1:256
            imy(256+1-k,:,j) = imgy(j,:,k);
        end
    end
end
   
    
% imy = imy(md(1)+1:md(1)+sz(1),md(2)+1:md(2)+sz(2),md(3)+1:md(3)+sz(3));
sizeOftheImageAfterRec = size(imy)

end

