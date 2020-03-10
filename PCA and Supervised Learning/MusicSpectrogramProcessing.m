%%Compile Music 
clear all; close all; clc
mainpath = 'C:\Users\gavin\Documents\MATLAB\HW4Music\Bluegrass\newBlue';
D=dir(fullfile(mainpath,'*.m4a'));

%%  Construct Spectrogram images of music files
Lclip=5.4;
     for j=1:numel(D)
         musicfile = audioread(fullfile(mainpath,D(j).name));
         resampmusic = [];
         for i=1:(floor(length(musicfile)/2))
             resampmusic = [resampmusic musicfile(i*2)];
         end
         Fs = length(resampmusic)/Lclip;
         nclip = length(resampmusic);
         tp2=linspace(0,Lclip,nclip+1); 
         tp=tp2(1:nclip);
         if rem(length(resampmusic),2) == 0
             kp=(2*pi/Lclip)*[0:nclip/2-1 -nclip/2:-1];
             kps=fftshift(kp);
         else
             kp=(2*pi/Lclip)*[0:nclip/2 -nclip/2:-1];
             kps=fftshift(kp);
         end      
%          figure(13)
%          plot((1:length(resampmusic))/Fs,resampmusic); 
%          xlabel('Time [sec]'); 
%          ylabel('Amplitude');
%          title('AC/DC'); 
%          drawnow 
%          p8 = audioplayer(resampmusic,Fs);
%          playblocking(p8);           
     width=1214;
    step=zeros(1,length(tp));
    mask=ones(1,2*width+1);
    musicstept_spec=[]; 
    tslide=(width+1):800:(length(step)-width);
    for z=1:length(tslide)
         step=zeros(1,length(tp));
         step(tslide(z)-width:1:tslide(z)+width)=mask;
         musicstep=step.*resampmusic; 
         musicstept=fft(musicstep); 
         musicstept_spec=[musicstept_spec; 
         abs(fftshift(musicstept))];     
    end

tslidep=linspace(0,Lclip,length(tslide));
f=figure('visible', false);
pcolor(tslidep,kps/(2*pi),musicstept_spec.'), 
shading interp 
set(gca,'XTick',[],'YTick',[],'Ylim',[0 2000],'Fontsize',[14]) 
colormap(hsv)
newFileNameChar = 'Bluegrass' + string(j+13) + '.jpg';
print(newFileNameChar,'-djpeg'); 
close(f)            
filescomplete = j
    end
    
    
    %% Gets rid of white border and grayscales images
    
     mainpath = 'C:\Users\gavin\Documents\MATLAB\HW4Music\Test3Samples';
D=dir(fullfile(mainpath,'*.jpg'));   
%     figure(8)
    for j=1:numel(D)
        thisfile = fullfile(mainpath,D(j).name);
        SpectroImg = imread(fullfile(mainpath,D(j).name));
        SpectroGray = double(rgb2gray(SpectroImg));
        SpectroCrop = SpectroGray(75:584,120:785);
%         imshow(uint8(SpectroCrop))
%         pause(.2)
        newFileName = split(thisfile,".");
        filename = string(newFileName(1,1));
        newFileNameChar = filename + '_cropped.jpg';
        imwrite(uint8(SpectroCrop),newFileNameChar); 
    end
    
      