%%Compile Faces and Organize into a single  matrix
clear all; close all; clc

mainpath = 'C:\Users\gavin\Documents\MATLAB\yalefaces';
cropface=[];
cropaveface=[];

D=dir(fullfile(mainpath,'subject*'));

  figure(1)
    for j=1:8
         subplot(3,3,j)
         face = imread(fullfile(mainpath,D(j).name));
         imshow(face)
        blurface = imresize(face, [64,64]);
        aveface = double(blurface) - mean(double(blurface));
        skinnface = reshape(blurface,1,4096);
        skinnaveface = reshape(aveface,1,4096);
        cropaveface = [skinnaveface' cropaveface];
        cropface = [skinnface' cropface];
     end


%% Pull out first 4 face sets

face1 = cropaveface(:,1:11);
face2 = cropaveface(:,12:22);
face3 = cropaveface(:,23:33);
face4 = cropaveface(:,34:44);

% figure(1)
% for j=1:25
%     normface=reshape(face2(:,j),64,64);
%     subplot(5,5,j)
%     imshow(normface)
% end

%% SVD of whole set

face_wave = dc_wavelet(cropaveface);

[U,S,V] = svd(face_wave,'econ');

figure(2)
for j=1:12
  subplot(3,4,j) 
  ut1 = reshape(U(:,j),32,32); 
  ut2 = ut1(32:-1:1,:); 
  pcolor(ut2), colormap(jet)
  set(gca,'Xtick',[],'Ytick',[])
  colorbar
end
%% Plot and assess Singular Values

figure(3)
subplot(2,1,1) 
plot(diag(S),'ko','Linewidth',[2]) 
set(gca,'Fontsize',[14]) 
axis([0 165 0 2*10^4])
xlabel('Singular Values')
ylabel('\sigma')
subplot(2,1,2) 
semilogy(diag(S),'ko','Linewidth',[2]) 
set(gca,'Fontsize',[14])
axis([0 165 0 2*10^4])
xlabel('Singular Values')
ylabel('\sigma')

%Singular value assessment
sig=diag(S);
energy1=sig(1)/sum(sig)
energy2=sig(2)/sum(sig) 
eng1perc=sum(sig(1:2))/sum(sig)
eng10perc=sum(sig(1:17))/sum(sig)
eng25perc=sum(sig(1:41))/sum(sig)
eng50perc=sum(sig(1:83))/sum(sig)
eng=sum(sig(1:128))/sum(sig)

%% Analyze V
k=0;
figure(4)
for j=1:3
    k=k+1;
  subplot(3,4,k) 
  plot(1:11,V(1:11,j),'ko-') 
  k=k+1;
  subplot(3,4,k) 
  plot(12:22,V(12:22,j),'ko-')
  k=k+1;
  subplot(3,4,k) 
  plot(23:33,V(23:33,j),'ko-')
  k=k+1;
   subplot(3,4,k) 
  plot(34:44,V(34:44,j),'ko-')
end

subplot(3,4,1), set(gca,'Xlim',[1 11],'Fontsize',[14]), title('Subject 1') 
subplot(3,4,2), set(gca,'Xlim',[12 22],'Fontsize',[14]), title('Subject 2')
subplot(3,4,3), set(gca,'Xlim',[23 33],'Fontsize',[14]), title('Subject 3') 
subplot(3,4,4), set(gca,'Xlim',[34 44],'Fontsize',[14]), title('Subject 4')
subplot(3,4,5), set(gca,'Xlim',[1 11],'Fontsize',[14]) 
subplot(3,4,6), set(gca,'Xlim',[12 22],'Fontsize',[14])
subplot(3,4,7), set(gca,'Xlim',[23 33],'Fontsize',[14])
subplot(3,4,8), set(gca,'Xlim',[34 44],'Fontsize',[14])
subplot(3,4,9), set(gca,'Xlim',[1 11],'Fontsize',[14]) 
subplot(3,4,10), set(gca,'Xlim',[12 22],'Fontsize',[14])
subplot(3,4,11), set(gca,'Xlim',[23 33],'Fontsize',[14])
subplot(3,4,12), set(gca,'Xlim',[34 44],'Fontsize',[14])


%% Reproduction
% figure(6)
% ff1=zeros(32,32);
% ff2=ff1;
% for j=1:61
%      
%     ut1 = reshape(U(:,j),32,32); 
%     ffo = ut1*sig(j)*V(24,j);
%     ff1 = ff1 + ffo;
%     ff = ut1*sig(j)*V(70,j);
%     ff2 = ff2 + ff;
% %     imshow(uint8(ffo))
%     ff1a = ff1(32:-1:1,:); 
%     if j >= 50
%     subplot(3,4,(62-j)) 
%     pcolor(ff1a/j), colormap(gray)
%   set(gca,'Xtick',[],'Ytick',[])
%   colorbar
%     end
% end 

figure(5)
ff1=zeros(32,32);
ff2=ff1;
for j=1:12
    ut1 = reshape(U(:,j),32,32); 
    ffo = ut1*sig(j)*V(8,j);
    ff1 = ff1 + ffo;
    ffoa = ut1*sig(j)*V(46,j);
    ff2 = ff2 + ffoa;

    ff1a = ff1(32:-1:1,:); 
    ff2a = ff2(32:-1:1,:); 
     subplot(3,4,j) 
%     imshow(uint8(ffo))
     pcolor(ff2a/j), colormap(gray)
   set(gca,'Xtick',[],'Ytick',[])

end 

figure(7)
subplot(2,2,1)
pcolor(ff1a/j), colormap(gray)
set(gca,'Xtick',[],'Ytick',[])
normface=reshape(cropface(:,8),64,64);
subplot(2,2,2)
imshow(normface)

subplot(2,2,3)
pcolor(ff2a/j), colormap(gray)
set(gca,'Xtick',[],'Ytick',[])
normface=reshape(cropface(:,46),64,64);
subplot(2,2,4)
imshow(normface)


