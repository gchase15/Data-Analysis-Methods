%% Compiles elements of test data into individual matrices
clear all; close all; clc
     mainpath = 'C:\Users\gavin\Documents\MATLAB\HW4Music\Test1Samples';
D=dir(fullfile(mainpath,'ACDC*cropped.jpg'));  
Band1Samp=[];
Band1AveSamp=[];

    for j=1:numel(D)
           song = imread(fullfile(mainpath,D(j).name));
         blurspectro = double(imresize(song, [64,64]));
%          imshow(blurspectro)
         avespectro = double(blurspectro) - mean(double(blurspectro));
         fatspectro =reshape(blurspectro,1,4096);
         fatavespectro = reshape(avespectro,1,4096);
         Band1AveSamp = [fatavespectro' Band1AveSamp];
         Band1Samp = [fatspectro' Band1Samp];
    end

%     D=dir(fullfile(mainpath,'ThirdDay*cropped.jpg'));  
% Band2Samp=[];
% Band2AveSamp=[];
% 
%     for j=1:numel(D)
%            song = imread(fullfile(mainpath,D(j).name));
%          blurspectro = double(imresize(song, [64,64]));
% %          imshow(blurspectro)
%          avespectro = double(blurspectro) - mean(double(blurspectro));
%          fatspectro =reshape(blurspectro,1,4096);
%          fatavespectro = reshape(avespectro,1,4096);
%          Band2AveSamp = [fatavespectro' Band2AveSamp];
%          Band2Samp = [fatspectro' Band2Samp];
%     end
%     
    D=dir(fullfile(mainpath,'Williams*cropped.jpg'));  
Band3Samp=[];
Band3AveSamp=[];

    for j=1:numel(D)
           song = imread(fullfile(mainpath,D(j).name));
         blurspectro = double(imresize(song, [64,64]));
%          imshow(blurspectro)
         avespectro = double(blurspectro) - mean(double(blurspectro));
         fatspectro =reshape(blurspectro,1,4096);
         fatavespectro = reshape(avespectro,1,4096);
         Band3AveSamp = [fatavespectro' Band3AveSamp];
         Band3Samp = [fatspectro' Band3Samp];
    end   
    
        D=dir(fullfile(mainpath,'ArmySongs*cropped.jpg'));  
Band4Samp=[];
Band4AveSamp=[];

    for j=1:numel(D)
           song = imread(fullfile(mainpath,D(j).name));
         blurspectro = double(imresize(song, [64,64]));
%          imshow(blurspectro)
         avespectro = double(blurspectro) - mean(double(blurspectro));
         fatspectro =reshape(blurspectro,1,4096);
         fatavespectro = reshape(avespectro,1,4096);
         Band4AveSamp = [fatavespectro' Band4AveSamp];
         Band4Samp = [fatspectro' Band4Samp];
    end
    
    %% Wavelet Transform
    
    band1_wave = dc_wavelet(Band1AveSamp);
%     band2_wave = dc_wavelet(Band2AveSamp);
    band3_wave = dc_wavelet(Band3AveSamp);
    band4_wave = dc_wavelet(Band4AveSamp);
    feature = 36; %must be <= the number of samples 
    
   [U,S,V]=svd([band1_wave,band3_wave,band4_wave],0);
    
    %% U
    
    figure(2)
for j=1:12
  subplot(3,4,j) 
  ut1 = reshape(U(:,j),32,32); 
  ut2 = ut1(32:-1:1,:); 
  pcolor(ut2), colormap(jet)
  set(gca,'Xtick',[],'Ytick',[])
  colorbar
end
    

%% Analyze S

figure(3)
subplot(2,1,1) 
plot(diag(S),'ko','Linewidth',[2]) 
set(gca,'Fontsize',[14]) 
axis([0 36 0 2*10^4])
xlabel('Singular Values')
ylabel('\sigma')
subplot(2,1,2) 
semilogy(diag(S),'ko','Linewidth',[2]) 
set(gca,'Fontsize',[14])
axis([0 36 0 2*10^4])
xlabel('Singular Values')
ylabel('\sigma')

%Singular value assessment
sig=diag(S);
energy1=sig(1)/sum(sig)
energy2=sig(2)/sum(sig) 
%eng1perc=sum(sig(1:2))/sum(sig)
eng10perc=sum(sig(1:4))/sum(sig)
eng25perc=sum(sig(1:9))/sum(sig)
eng50perc=sum(sig(1:18))/sum(sig)
eng=sum(sig(1:35))/sum(sig)
    
%% Analyze V
k=0;
figure(4)
for j=1:3
    k=k+1;
  subplot(3,3,k) 
  plot(1:12,V(1:12,j),'ko-') 
  k=k+1;
  subplot(3,3,k) 
  plot(13:24,V(13:24,j),'ko-')
  k=k+1;
  subplot(3,3,k) 
  plot(25:36,V(25:36,j),'ko-')
end

subplot(3,3,1), set(gca,'Xlim',[1 12],'Fontsize',[14]), title('ACDC (Hard Rock)') 
subplot(3,3,2), set(gca,'Xlim',[13 24],'Fontsize',[14]), title('Third Day (CCM)')
subplot(3,3,3), set(gca,'Xlim',[25 36],'Fontsize',[14]), title('John Williams (Classical)') 
subplot(3,3,4), set(gca,'Xlim',[1 12],'Fontsize',[14]) 
subplot(3,3,5), set(gca,'Xlim',[13 24],'Fontsize',[14])
subplot(3,3,6), set(gca,'Xlim',[25 36],'Fontsize',[14])
subplot(3,3,7), set(gca,'Xlim',[1 12],'Fontsize',[14]) 
subplot(3,3,8), set(gca,'Xlim',[13 24],'Fontsize',[14])
subplot(3,3,9), set(gca,'Xlim',[25 36],'Fontsize',[14])

%% 3D version of V analyzation 

figure(5)
plot3(V(1:12,1),V(1:12,2),V(1:12,3),'ko','Linewidth',[2]) 
hold on
plot3(V(13:24,1),V(13:24,2),V(13:24,3),'ro','Linewidth',[2])
hold on
plot3(V(25:36,1),V(25:36,2),V(25:36,3),'go','Linewidth',[2])
% hold on
% plot3(V(37:48,1),V(37:48,2),V(37:48,3),'bo','Linewidth',[2])
grid on, xlabel('V1'), ylabel('V2'), zlabel('V3')

%% Classification Task

band1prelda = []; band2prelda = []; band3prelda = []; band4prelda = [];
band1prenb = []; band2prenb = []; band3prenb = []; band4prenb = [];
band1pretree = []; band2pretree = []; band3pretree = []; band4pretree = [];

n=50
for k=1:n
    %pick a new randomization for each test.
q1 = randperm(12); q2 = randperm(12); q3 = randperm(12); q4 = randperm(12);

xband1 = V(1:12,1:10); xband2 = V(13:24,1:10);
xband3 = V(25:36,1:10); 
% xband4 = V(37:48,1:16);

xtrain = [xband1(q1(1:7),:); xband2(q2(1:7),:); xband3(q2(1:7),:)];
xtest = [xband1(q1(8:end),:); xband2(q2(8:end),:); xband3(q3(8:end),:)];
ctrain = [ones(7,1); 2*ones(7,1); 3*ones(7,1)];

%the classifiers (LDA,Naive Bayes, CART)
[predictlda] = classify(xtest,xtrain,ctrain);
nb = fitcnb(xtrain,ctrain);
predictnb = nb.predict(xtest);
tree=fitctree(xtrain,ctrain);
predicttree = predict(tree,xtest);


% Collect the results
band1prelda = [band1prelda predictlda(1:5)];
band2prelda = [band2prelda predictlda(6:10)];
band3prelda = [band3prelda predictlda(11:15)];
% band4prelda = [band4prelda predictlda(16:20)];

band1prenb = [band1prenb predictnb(1:5)];
band2prenb = [band2prenb predictnb(6:10)];
band3prenb = [band3prenb predictnb(11:15)];
% band4prenb = [band4prenb predictnb(16:20)];

band1pretree = [band1pretree predicttree(1:5)];
band2pretree = [band2pretree predicttree(6:10)];
band3pretree = [band3pretree predicttree(11:15)];
% band4pretree = [band4pretree predicttree(16:20)];

end

success1lda = sum(band1prelda(:) == 1)/(n*5)
success2lda = sum(band2prelda(:) == 2)/(n*5)
success3lda = sum(band3prelda(:) == 3)/(n*5)
% success4lda = sum(band4prelda(:) == 4)/5000

success1nb = sum(band1prenb(:) == 1)/(n*5)
success2nb = sum(band2prenb(:) == 2)/(n*5)
success3nb = sum(band3prenb(:) == 3)/(n*5)
% success4nb = sum(band4prenb(:) == 4)/5000

success1tree = sum(band1pretree(:) == 1)/(n*5)
success2tree = sum(band2pretree(:) == 2)/(n*5)
success3tree = sum(band3pretree(:) == 3)/(n*5)
% success4tree = sum(band4pretree(:) == 4)/5000


%% Display classification Results

figure(6)
subplot(1,3,1)
histogram(band1prelda,3,'BinWidth',.8)
axis([1 3 0 5000]), xticks([1.2 2 2.8]),xticklabels({'AC/DC', 'John Williams', 'Army Songs'})
set(gca,'YTick',[1000 2000 3000 4000 5000], 'Fontsize',[12])
subplot(1,3,2)
histogram(band2prelda,3,'BinWidth',.8), axis([1 3 0 5000])
xticks([1.2 2 2.8]),xticklabels({'AC/DC', 'John Williams', 'Army Songs'})
set(gca,'YTick',[1000 2000 3000 4000 5000], 'Fontsize',[12])
subplot(1,3,3)
histogram(band3prelda,3,'BinWidth',.8), axis([1 3 0 5000])
xticks([1.2 2 2.8]),xticklabels({'AC/DC', 'John Williams', 'Army Songs'})
set(gca,'YTick',[1000 2000 3000 4000 5000], 'Fontsize',[12])
% subplot(1,4,4)
% histogram(band4prenb,4), axis([1 4 0 5000])
% xticks([1.2 2 3 3.8]),xticklabels({'AC/DC', 'Third Day', 'John Williams', 'Army Songs'})
% set(gca,'YTick',[1000 2000 3000 4000 5000], 'Fontsize',[12])







       