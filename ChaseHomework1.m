clear all; close all; clc;
load Testdata

%% Setup data
L=15; %spatial domain
n=64; %Fourier modes 2^6
x2=linspace(-L,L,n+1); %Space discretization
x=x2(1:n); %space I actually use 
%(first point is the same as last and FFT works on periodic boundaries)
y=x; %the marble is moving in 3 dimensions 
z=x;

k=(2*pi/(2*L))*[0:n/2-1 -n/2:-1]; %standard shifted frame
ks=fftshift(k);%shifting k to unshift it

[X,Y,Z]=meshgrid(x,y,z);
[Kxs,Kys,Kzs]=meshgrid(ks,ks,ks);
[Kx,Ky,Kz]=meshgrid(k,k,k);

%% averaging the signal

Utave=zeros(n,n,n);
for j=1:20 
   Un(:,:,:)=reshape(Undata(j,:),n,n,n); 
   %this essentially takes each row of data and chops it into a cube
    Unt=fftn(Un(:,:,:));
    Utave=Utave+Unt;    
end
Utave=fftshift(Utave)/j;
Utave=abs(Utave)/max(abs(Utave), [], 'all');
%scales everything as a ratio of the max value 

 %if the max point in the cube is the signal, 
 %after it is scaled to some portion of 1, the true signal 
 %should lie around the area at or near the max of 1
 close all, isosurface(Kxs,Kys,Kzs,Utave, .75:1) 
 %detection threshold is around .75
 axis([-7 7 -7 7 -7 7]), grid on, drawnow;
 xlabel ('wavenumber (kx)'), ylabel ('wavenumber (ky)'),zlabel ('wavenumber (kz)'), 
 pause(2)
 
 freqpos=find((Utave/max(Utave,[],'all'))==1);
 [freqrow,freqcol,freqpage]=ind2sub(size(Utave),freqpos);
 
 xfreqpos=ks(1)+(freqcol/n)*(2*ks(64))%x wavenumber
 yfreqpos=ks(1)+(freqrow/n)*(2*ks(64))%y wavenumber
 zfreqpos=ks(1)+(freqpage/n)*(2*ks(64)) %z wavenumber
 %% the filter
 
tau=0.15;
 filter=exp((-tau*(Kx-xfreqpos).^2)+(-tau*(Ky-yfreqpos).^2)+(-tau*(Kz-zfreqpos).^2)); 
 %gaussian filter centered on max value from averaging the noisy signal
 xslice=1.88;yslice=-1.05; zslice=0.05; 
 %the cut marks to show the filter distribution
 close all, slice(Kxs,Kys,Kzs,fftshift(filter), xslice,yslice,zslice)
 axis([-7 7 -7 7 -7 7]), grid on, drawnow 
 xlabel ('wavenumber (kx)'), ylabel ('wavenumber (ky)'),zlabel ('wavenumber (kz)')
 
 %% finding the marble
 
 for j=1:20
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   Unt=fftn(Un(:,:,:));
   Untf=filter.*Unt;
   Unf=abs(ifftn(Untf));
   path(j)=find((Unf/max(Unf,[],'all'))==1);
   close all, isosurface(X,Y,Z,Unf/max(Unf, [], 'all'), .5:1) 
   %detection threshold of .5
   axis([-20 20 -20 20 -20 20]), grid on, drawnow  
   xlabel ('distance(x)'), ylabel ('distance(y)'),zlabel ('distance(z)')
   pause(.2)
 end
 
 %% plotting the marble's path
for j=1:20
 [mrow(j),mcol(j),mpage(j)]=ind2sub(size(Unf),path(j));
end
%this for loop deconstructs path into vectors containing x,y, and z
%coordinates
 xpath=-L+(mcol/n)*(2*L);% x positions
 ypath=-L+(mrow/n)*(2*L);% y positions
 zpath=-L+((mpage/n)*(2*L));% z positions
  plot3(xpath,ypath,zpath, '-o')
  grid on 
 xlabel ('distance(x)'), ylabel ('distance(y)'),zlabel ('distance(z)')



 
 
 
 
 
 
