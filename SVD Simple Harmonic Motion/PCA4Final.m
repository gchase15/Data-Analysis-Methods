%% Load Data
clear all; close all; clc
addpath('HW3')
load cam1_1
load cam2_1 
load cam3_1
load cam1_2
load cam2_2 
load cam3_2
load cam1_3
load cam2_3 
load cam3_3
load cam1_4
load cam2_4 
load cam3_4
implay(vidFrames1_4)
implay(vidFrames2_4)
implay(vidFrames3_4)
%% Cutting the size of the video
%   Find the region that the paint can occupies during the entire video 
%   and cut each frame to include only that section. 
%   X positions move from left to right.
%   Y positions move from top to bottom. 

%Synchronize videos first by cutting initial frames on 2/3 videos
shortFrames1_4=vidFrames1_4(:,:,:,17:end);
shortFrames2_4=vidFrames2_4(:,:,:,20:end);
shortFrames3_4=vidFrames3_4(:,:,:,14:end);

%Find the minimum number of frames to make all videos the same length
F1=size(shortFrames1_4,4);
F2=size(shortFrames2_4,4);
F3=size(shortFrames3_4,4);
FrameLength=[F1 F2 F3];
FrameLength=min(FrameLength);

%Pick the "box" that you want to cut out of each frame so you only have the
%movement of the paint can.  
 imshow(shortFrames1_4(:,:,:,1));
 [x1,y1]=ginput(2);
 imshow(shortFrames2_4(:,:,:,1));
 [x2,y2]=ginput(2);
 imshow(shortFrames3_4(:,:,:,1));
 [x3,y3]=ginput(2);  close

 %Shrink each frame to fit the box picked above and cut any frames after
 %FrameLength ends.
 cutFrames1_4=shortFrames1_4(min(y1):max(y1),min(x1):max(x1),:,1:FrameLength);
 cutFrames2_4=shortFrames2_4(min(y2):max(y2),min(x2):max(x2),:,1:FrameLength);
 cutFrames3_4=shortFrames3_4(min(y3):max(y3),min(x3):max(x3),:,1:FrameLength);
 
 %Confirm you got it right and enjoy the show
 implay(cutFrames1_4)
 implay(cutFrames2_4)
 implay(cutFrames3_4)

 %% Track Video 1
 FrameLength=376;
 %First Part of For loop converts each Frame into grayscale. It preserves a
 %copy of this image as FrameGray1 for comparison later. 
 xpos14=[];
 ypos14=[];
 for j=1:FrameLength
 Frame=cutFrames1_4(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

%The for and if loops below filter out colors other than bright white. 
%The Detection threshold is set to 253 to guard against the rare case where 
%the light's color isn't shown at 255. 
 for j=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(j,i) > 244;
             FrameGray(j,i)=255;
         else FrameGray(j,i)=0;
         end
     end
 end

 %These figures serve to validate that what remains in the filtered frame 
 %are the viable candidates for tracking. They should be commented out once
 %validation is complete. 
 figure(1)
 imshow(FrameGray)
 figure(2)
 imshow(FrameGray1)
 
 %Pick out the location of the first white pixel. Congratulations, you've
%extracted out the x and y positions for a paint can. 
 [ylight, xlight]=ind2sub(size(FrameGray),find(FrameGray,1));
 xpos14=[xpos14;
     xlight];
 ypos14=[ypos14;
     ylight];
 end
 
  %Subtract the mean from each row
 xpos14=xpos14-mean(xpos14); 
 ypos14=ypos14-mean(ypos14);
 %% Track Video 2
 
 xpos24=[];
 ypos24=[];
  for j=1:FrameLength
 Frame=cutFrames2_4(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

 for j=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(j,i) > 250;
             FrameGray(j,i)=255;
         else FrameGray(j,i)=0;
         end
     end
 end
%  figure(1)
%  imshow(FrameGray)
%  figure(2)
%  imshow(FrameGray1)
white=find(FrameGray);

 [ylight, xlight]=ind2sub(size(FrameGray),white(end));
 xpos24=[xpos24;
     xlight];
  ypos24=[ypos24;
     ylight];
   end
 xpos24=xpos24-(sum(xpos24)/length(xpos24));
 ypos24=ypos24-(sum(ypos24)/length(ypos24));
  %% Track Video 3
 
 xpos34=[];
 ypos34=[];
  for j=1:FrameLength
 Frame=cutFrames3_4(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

 for j=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(j,i) > 230;
             FrameGray(j,i)=255;
         else FrameGray(j,i)=0;
         end
     end
 end
%  figure(1)
%  imshow(FrameGray)
%  figure(2)
%  imshow(FrameGray1)

white=find(FrameGray);

 [ylight, xlight]=ind2sub(size(FrameGray),white(end));
 xpos34=[xpos34;
     xlight];
  ypos34=[ypos34;
     ylight];
   end
 xpos34=xpos34-(sum(xpos34)/length(xpos34));
 ypos34=ypos34-(sum(ypos34)/length(ypos34));
 %% The SVD
 
 %create a matrix of the three video's x and y pos over time
 A4=[xpos14';
     ypos14';
     xpos24';
     ypos24';
     xpos34';
     ypos34'];
 
  [u4,s4,v4]=svd(A4','econ'); %U and V are rotations, s is the stretch
  A4Cov=cov(A4');
  
  variance4=(s4.^2)./(size(s4,1)-1);
  var4=diag(variance4);
  
  sig4=diag(s4); %Collects the sigma values from s into a vector
  energy1=sig4(1)/sum(sig4)%Proportion of first singular value amongst all of them
  energy2=sig4(2)/sum(sig4) 
  eng4=sum(sig4(1:3))/sum(sig4)
  %% Describe Motion in mean frame 
  
        meanFrames=double(rgb2gray(cutFrames2_4(:,:,:,1)));
    for j=2:size(cutFrames2_4,4)
        meanFrames=meanFrames + double(rgb2gray(cutFrames2_4(:,:,:,j)));
    end
     meanFrames=uint8(meanFrames/size(cutFrames2_4,4));
    figure(25)
    imshow(meanFrames)
    hold on
    [xcent,ycent]=ginput(1);
    plot([xcent,xcent+100*v4(3,1)],[ycent,ycent+100*v4(4,1)],'g','linewidth',2)
    plot([xcent,xcent+100*v4(3,3)],[ycent,ycent+100*v4(4,2)],'m','linewidth',2)
    plot([xcent,xcent+100*v4(3,3)],[ycent,ycent+100*v4(4,3)],'k','linewidth',2)
    legend('mode 1','mode 2','mode 3','Location','SouthEast') 
    %% Plotting SVD
  x4=linspace(1,376,376);
  
  figure(8)
 subplot(1,2,1)
 plot(sig4,'ko','Linewidth',[1.5])
 axis([0 7 0 1600])
 set(gca,'Fontsize',[13],'Xtick',[0 1 2 3 4 5 6])
 xlabel('singular values'), ylabel('\sigma')
 grid on
 text(6.1,1.25*10^3,'(a)','Fontsize',[13])
 subplot(1,2,2)
 plot(var4,'ko','Linewidth',[1.5])
 axis([0 7 0 4*10^5])
 set(gca,'Fontsize',[13],'Xtick',[0 1 2 3 4 5 6]) 
 xlabel('singular values'), ylabel('\sigma^2')
 grid on
 text(6.1,3.16*10^5,'(b)','Fontsize',[13])
 
figure(12)
 plot(x4,u4(:,3),'g-.',x4,u4(:,2),'r',x4,u4(:,1),'k','Linewidth',[2]) 
 set(gca,'Fontsize',[13])
 legend('mode 3','mode 2','mode 1','Location','SouthEast')
 axis([0 225 -.25 .25])
 xlabel('Frames'), ylabel('U')
  
   %% Reproduction
   
  t=1:length(x4);
  figure(21)
  subplot(1,2,1)
  plot(t,xpos14,t,ypos14,t,xpos24,t,ypos24,t,xpos34,t,ypos34); 
   axis([0 385 -150 150])
  legend('x1','y1','x2','y2','x3','y3','Location','SouthEast')
  
  for j=3:3
  ff4=u4(:,1:j)*s4(1:j,1:j)*v4(:,1:j)'; % modal projections 
    subplot(1,2,2)
    plot(x4,ff4)
      legend('x1','y1','x2','y2','x3','y3','Location','SouthEast')
    axis([0 385 -150 150])
  end
  
subplot(1,2,1), text(19,130,'(a)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean')
subplot(1,2,2), text(19,130,'(b)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean') 
% subplot(2,2,3), text(19,130,'(c)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean') 
% subplot(2,2,4), text(19,130,'(d)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean')
  
  motion14=(xpos14(12:end)+xpos24(12:end)+ypos34(12:end))/3;
  motion24=(ypos14(12:end)+ypos24(12:end)+xpos34(12:end))/3;
  PCAmot14=(ff4(12:end,1)+ff4(12:end,3)+ff4(12:end,6))/3;
  PCAmot24=(ff4(12:end,2)+ff4(12:end,4)+ff4(12:end,5))/3;
  tm4=1:length(motion14);
  
  figure(4)
  subplot(1,2,1)
  plot(tm4,motion14, 'k', tm4,motion24,'k-.', 'Linewidth',[2])
  axis([0 380 -80 80])
  xlabel('Frames'), ylabel('Displacement From Mean')
  legend('Primary Movement', 'Secondary Movement', 'Location','SouthEast')
  subplot(1,2,2)
  plot(tm4,PCAmot14, 'k', tm4,PCAmot24,'k-.', 'Linewidth',[2])
  axis([0 380 -80 80])
  xlabel('Frames'), ylabel('Displacement From Mean')
  legend('Primary Movement', 'Secondary Movement', 'Location','SouthEast')
  
  
 