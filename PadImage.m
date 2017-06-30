function [ img ] = PadImage(img,pad)   


sz = size(img);
md = ([pad-sz(1),pad-sz(2),pad-sz(3)]./2);
  
    if(floor(md(1)) == md(1))
        img = padarray(img,[md(1),0,0]);
                
    else
        md(1) = floor(md(1));
        img = padarray(img,[md(1)+1,0,0],'pre');
        img = padarray(img,[md(1),0,0],'post');
    end
    
   
    if(floor(md(2)) == md(2))
        img = padarray(img,[0,md(2),0]);
                
    else
        md(2) = floor(md(2));
        img = padarray(img,[0,md(2)+1,0],'pre');
        img = padarray(img,[0,md(2),0],'post');
    end
       
  
    if(floor(md(3)) == md(3))
       img = padarray(img,[0,0,md(3)]);
                
    else
        md(3) = floor(md(3));
        img = padarray(img,[0,0,md(3)+1],'pre');
        img = padarray(img,[0,0,md(3)],'post');
    end                            
     
    sizeafterpadding = size(img)
    
end
