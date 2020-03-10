%%Compile Faces and Organize into a single  matrix
clear all; close all; clc

fold=1:38;
mainpath = 'C:\Users\gavin\Documents\MATLAB\CroppedYale\yaleB';
cropface=[];
cropaveface=[];

%First for loop creates the path to the pictures
for k=1:length(fold)
    foldnum=string(k);
    newpath=strcat(mainpath,foldnum);
    D=dir(fullfile(newpath,'*.pgm'));

    %Second for loop pulls out the images from the list "reduced set", 
    %reduces the resolution, turns it into a row, and adds it to a 
    %collective matrix of all faces. Pixels are still uint8.  
    for j=1:numel(D)
        face = imread(fullfile(newpath,D(j).name));
        blurface = imresize(face, [64,64]);
        aveface = double(blurface) - mean(double(blurface));
        fatface = reshape(blurface,1,4096);
        fataveface = reshape(aveface,1,4096);
        cropaveface = [fataveface' cropaveface];
        cropface = [fatface' cropface];
    end
end

%% Pull out first 4 face sets

face1 = cropaveface(:,1:64);
face2 = cropaveface(:,65:128);
face3 = cropaveface(:,129:192);
face4 = cropaveface(:,193:256);

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
axis([0 1024 0 8*10^4])
xlabel('Singular Values')
ylabel('\sigma')
subplot(2,1,2) 
semilogy(diag(S),'ko','Linewidth',[2]) 
set(gca,'Fontsize',[14])
axis([0 1024 0 8*10^4])
xlabel('Singular Values')
ylabel('\sigma')

%Singular value assessment
sig=diag(S);
energy1=sig(1)/sum(sig)
energy2=sig(2)/sum(sig) 
eng1perc=sum(sig(1:3))/sum(sig)
eng10perc=sum(sig(1:102))/sum(sig)
eng25perc=sum(sig(1:256))/sum(sig)
eng50perc=sum(sig(1:512))/sum(sig)
eng=sum(sig(1:789))/sum(sig)

%% Analyze V
k=0;
figure(4)
for j=1:3
    k=k+1;
  subplot(3,4,k) 
  plot(1:64,V(1:64,j),'ko-') 
  k=k+1;
  subplot(3,4,k) 
  plot(65:128,V(65:128,j),'ko-')
  k=k+1;
  subplot(3,4,k) 
  plot(129:192,V(129:192,j),'ko-')
  k=k+1;
   subplot(3,4,k) 
  plot(193:256,V(193:256,j),'ko-')
end

subplot(3,4,1), set(gca,'Xlim',[0 64],'Fontsize',[14]), title('Subject 1') 
subplot(3,4,2), set(gca,'Xlim',[65 128],'Fontsize',[14]), title('Subject 2')
subplot(3,4,3), set(gca,'Xlim',[129 192],'Fontsize',[14]), title('Subject 3') 
subplot(3,4,4), set(gca,'Xlim',[193 256],'Fontsize',[14]), title('Subject 4')
subplot(3,4,5), set(gca,'Xlim',[0 64],'Fontsize',[14]) 
subplot(3,4,6), set(gca,'Xlim',[65 128],'Fontsize',[14])
subplot(3,4,7), set(gca,'Xlim',[129 192],'Fontsize',[14])
subplot(3,4,8), set(gca,'Xlim',[193 256],'Fontsize',[14])
subplot(3,4,9), set(gca,'Xlim',[0 64],'Fontsize',[14]) 
subplot(3,4,10), set(gca,'Xlim',[65 128],'Fontsize',[14])
subplot(3,4,11), set(gca,'Xlim',[129 192],'Fontsize',[14])
subplot(3,4,12), set(gca,'Xlim',[193 256],'Fontsize',[14])


%% Reproduction
figure(6)
ff1=zeros(32,32);
ff2=ff1;
for j=1:61
     
    ut1 = reshape(U(:,j),32,32); 
    ffo = ut1*sig(j)*V(24,j);
    ff1 = ff1 + ffo;
    ff = ut1*sig(j)*V(70,j);
    ff2 = ff2 + ff;
%     imshow(uint8(ffo))
    ff1a = ff1(32:-1:1,:); 
    if j >= 50
    subplot(3,4,(62-j)) 
    pcolor(ff1a/j), colormap(gray)
  set(gca,'Xtick',[],'Ytick',[])
  colorbar
    end
end 

figure(5)
ff1=zeros(32,32);
ff2=ff1;
for j=1:102
    ut1 = reshape(U(:,j),32,32); 
    ffo = ut1*sig(j)*V(24,j);
    ff1 = ff1 + ffo;
    ffoa = ut1*sig(j)*V(70,j);
    ff2 = ff2 + ffoa;
%     imshow(uint8(ffo))
    ff1a = ff1(32:-1:1,:); 
    ff2a = ff2(32:-1:1,:); 
%     subplot(3,4,j) 
%     pcolor(ff1a/j), colormap(gray)
%   set(gca,'Xtick',[],'Ytick',[])

end 

figure(7)
subplot(2,2,1)
pcolor(ff1a/j), colormap(gray)
set(gca,'Xtick',[],'Ytick',[])
normface=reshape(cropface(:,24),64,64);
subplot(2,2,2)
imshow(normface)

subplot(2,2,3)
pcolor(ff2a/j), colormap(gray)
set(gca,'Xtick',[],'Ytick',[])
normface=reshape(cropface(:,70),64,64);
subplot(2,2,4)
imshow(normface)


