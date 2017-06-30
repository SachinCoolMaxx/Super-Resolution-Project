% addpath('/home/sachin/Documents/distmesh/');
% figure;
% unitsize = .; % the unit length of distmeshsurface must be 1/X
% fd = @(p) dsphere(p,0,0,0,1);
% [vertices,faces] = distmeshsurface(fd,@huniform,unitsize,[-1 -1 -1;1 1 1])
% 
% % display sphere of HistFaces
% displaySphere(faces,vertices,HistFaces,...
% % 'histogram of orientations ','savename.fig');


% Amp = 1;
% freqHz = 12000;
% fsHz = 65536;
% dt = 1/fsHz;
% index = 1;
% sine = [];
% sampleNumber = [];
% for t = 0:dt:1-dt
%    sine(index) = Amp * sin(2*pi*freqHz*t);
%    sampleNumber(index) = index;
%    index = index + 1;
% end
% 
% % alternative to the above loop
% % sine = Amp*sin(2*pi*freq*(0:dt:1-dt));
% 
% N = 65536;
% transform = fft(sine,N)/N;
% magTransform = abs(transform);
% imagesc((magTransform));
% 
% faxis = linspace(-fsHz/2,fsHz/2,N);
% plot(faxis/1000,fftshift(magTransform));
% axis([-40 40 0 0.6])
% xlabel('Frequency (KHz)')

% 
% 
% 
% B=imread('/home/sachin/Documents/Results/Week3/Maxp9direcns.png');
% C=fft2(B);
% A= [1,0,1;1,0,1;1,1,1];
% figure('Name','Rand')
% imagesc(abs(uint8(A)))






 
% 
% 
% 
% colorRGB = [0 0.5 1;de2bi([0:7])];   
% 
% 
% middle=[.1,.2,.3;.4,.5,.6;.7,.8,.9];
% figure('Name','middle')
% RGB = ind2rgb(middle,colorRGB);
% image(RGB);
% 
% middle=[.11,.2,.3;.4,.5,.1;.75,.8,1];
% figure('Name','contour')
% contour3(middle)     
% imagesc(middle)
% colormap




addpath('/iacl/pg17/sachin/for_sachin/code/');

I = imread('/home/sachin/Downloads/cameraman.jpg');
imshow(I);
size(I)
Ix = RotAboutx( I ,pi/2 );
figure
imshow(Ix);



