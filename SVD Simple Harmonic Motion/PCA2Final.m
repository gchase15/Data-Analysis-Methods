%% Load Data
% clear all; close all; clc
addpath('HW3')
% load cam1_1
% load cam2_1 
% load cam3_1
load cam1_2
load cam2_2 
load cam3_2
% load cam1_3
% load cam2_3 
% load cam3_3
% load cam1_4
% load cam2_4 
% load cam3_4
 implay(vidFrames1_2)
 implay(vidFrames2_2)
 implay(vidFrames3_2)

%% Cutting the size of the video

%   Find the region that the paint can occupies during the entire video 
%   and cut each frame to include only that section. 
%   X positions move from left to right.
%   Y positions move from top to bottom. 

%Synchronize videos first by cutting initial frames videos
shortFrames1_2=vidFrames1_2(:,:,:,14:end);
shortFrames2_2=vidFrames2_2(:,:,:,:);
shortFrames3_2=vidFrames3_2(:,:,:,18:end);

%Find the minimum number of frames to make all videos the same length
F1=size(shortFrames1_2,4);
F2=size(shortFrames2_2,4);
F3=size(shortFrames3_2,4);
FrameLength=[F1 F2 F3];
FrameLength=min(FrameLength);

%Pick the "box" that you want to cut out of each frame so you only have the
%movement of the paint can.  
 imshow(shortFrames1_2(:,:,:,1));
 [x1,y1]=ginput(2);
 imshow(shortFrames2_2(:,:,:,1));
 [x2,y2]=ginput(2);
 imshow(shortFrames3_2(:,:,:,1));
 [x3,y3]=ginput(2);  close
 
 %Shrink each frame to fit the box picked above and cut any frames after
 %FrameLength ends.
 cutFrames1_2=shortFrames1_2(min(y1):max(y1),min(x1):max(x1),:,1:FrameLength);
 cutFrames2_2=shortFrames2_2(min(y2):max(y2),min(x2):max(x2),:,1:FrameLength);
 cutFrames3_2=shortFrames3_2(min(y3):max(y3),min(x3):max(x3),:,1:FrameLength);
 
 %Confirm you got it right and enjoy the show
 implay(cutFrames1_2)
 implay(cutFrames2_2)
 implay(cutFrames3_2)

 %% Track Video 1
 FrameLength=301;
 %First Part of For loop converts each Frame into grayscale. It preserves a
 %copy of this image as FrameGray1 for comparison later. 
 xpos12=[];
 ypos12=[];
 for j=1:FrameLength
 Frame=cutFrames1_2(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

%The for and if loops below filter out colors other than bright white. 
%The Detection threshold is set to 253 to guard against the rare case where 
%the light's color isn't shown at 255. 
 for j=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(j,i) > 250;
             FrameGray(j,i)=255;
         else FrameGray(j,i)=0;
         end
     end
 end

 %These figures serve to validate that what remains in the filtered frame 
 %are the viable candidates for tracking. They should be commented out once
 %validation is complete. 
%  figure(1)
%  imshow(FrameGray)
%  figure(2)
%  imshow(FrameGray1)
 
 %Pick out the location of the first white pixel. Congratulations, you've
%extracted out the x and y positions for a paint can. 
 [ylight, xlight]=ind2sub(size(FrameGray),find(FrameGray,1));
 xpos12=[xpos12;
     xlight];
  ypos12=[ypos12;
     ylight];
  end
  xpos12=xpos12-(sum(xpos12)/length(xpos12));
  ypos12=ypos12-(sum(ypos12)/length(ypos12));
 %% Track Video 2
 
 xpos22=[];
 ypos22=[];
  for j=1:FrameLength
 Frame=cutFrames2_2(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

 for j=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(j,i) > 247;
             FrameGray(j,i)=255;
         else FrameGray(j,i)=0;
         end
     end
 end
%  figure(1)
%  imshow(FrameGray)
%  figure(2)
%  imshow(FrameGray1)

 [ylight, xlight]=ind2sub(size(FrameGray),find(FrameGray,1));
 xpos22=[xpos22;
     xlight];
 ypos22=[ypos22;
     ylight];
 end
 xpos22=xpos22-(sum(xpos22)/length(xpos22));
 ypos22=ypos22-(sum(ypos22)/length(ypos22));
  %% Track Video 3
 
 xpos32=[];
 ypos32=[];
  for j=1:FrameLength
 Frame=cutFrames3_2(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

 for j=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(j,i) > 245;
             FrameGray(j,i)=255;
         else FrameGray(j,i)=0;
         end
     end
 end
%  figure(1)
%  imshow(FrameGray)
%  figure(2)
%  imshow(FrameGray1)

 [ylight, xlight]=ind2sub(size(FrameGray),find(FrameGray,1));
 xpos32=[xpos32;
     xlight];
  ypos32=[ypos32;
     ylight];
   end
 xpos32=xpos32-(sum(xpos32)/length(xpos32));
 ypos32=ypos32-(sum(ypos32)/length(ypos32));
 %% The SVD
 
 %create a matrix of the three video's x and y pos over time
 A2=[xpos12';
     ypos12';
     xpos22';
     ypos22';
     xpos32';
     ypos32'];
 
 [u2,s2,v2]=svd(A2','econ'); %U and V are rotations, s is the stretch
  A2Cov=cov(A2')
   variance2=(s2.^2)./(size(s2,1)-1);
 var2=diag(variance2);
  sig2=diag(s2); %Collects the sigma values from s into a vector
  energy1=sig2(1)/sum(sig2) %Proportion of first singular value amongst all of them
  energy2=sig2(2)/sum(sig2) 
  eng=sum(sig(1:3))/sum(sig)
  %%    Describe Motion in mean frame 
      meanFrames=double(rgb2gray(cutFrames1_2(:,:,:,1)));
    for j=2:size(cutFrames1_2,4)
        meanFrames=meanFrames + double(rgb2gray(cutFrames1_2(:,:,:,j)));
    end
     meanFrames=uint8(meanFrames/size(cutFrames1_2,4));
    figure(26)
    imshow(meanFrames)
    hold on
    [xcent,ycent]=ginput(1);
    plot([xcent,xcent+100*v2(1,1)],[ycent,ycent+100*v2(2,1)],'g','linewidth',2)
    plot([xcent,xcent+150*v2(1,2)],[ycent,ycent+150*v2(2,2)],'m','linewidth',2)
%      plot([xcent,xcent+100*v2(1,3)],[ycent,ycent+100*v2(2,3)],'k','linewidth',2)
    legend('mode 1','mode 2','mode 3','Location','SouthEast')
    %% Plotting SVD
  x2=linspace(0,301,301);
  figure(6)
  subplot(1,2,1)
  plot(sig2,'ko','Linewidth',[1.5])
  axis([0 7 0 1.6*10^3])
  set(gca,'Fontsize',[13],'Xtick',[0 1 2 3 4 5 6]) 
  xlabel('singular values'), ylabel('\sigma')
  grid on
  text(6.1,1.45*10^3,'(a)','Fontsize',[13])
  
 subplot(1,2,2)
  plot(var2,'ko','Linewidth',[1.5])
 axis([0 7 0 3.5*10^5])
 set(gca,'Fontsize',[13],'Xtick',[0 1 2 3 4 5 6]) 
 xlabel('singular values'), ylabel('\sigma^2')
 grid on
 text(6.1,3.16*10^5,'(b)','Fontsize',[13])

figure(10)
plot(x2,u2(:,2),'r',x2,u2(:,1),'k','Linewidth',[2]) 
set(gca,'Fontsize',[13])
legend('mode 3','mode 2','mode 1','Location','NorthEast')
axis([0 320 -.2 .15])
xlabel('Frames'), ylabel('U')

  %% Reproduction
  
  t=1:length(x2);
  figure(19)
  subplot(2,2,1)
  plot(t,xpos12,t,ypos12,t,xpos22,t,ypos22,t,xpos32,t,ypos32); 
   axis([0 310 -150 150])
  legend('x1','y1','x2','y2','x3','y3','Location','NorthEast')
  
  for j=1:3
  ff2=u2(:,1:j)*s2(1:j,1:j)*v2(:,1:j)'; % modal projections 
    subplot(2,2,j+1)
    plot(x2,ff2)
      legend('x1','y1','x2','y2','x3','y3','Location','NorthEast')
    axis([0 310 -150 150])
  end
  
subplot(2,2,1), text(19,130,'(a)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean')
subplot(2,2,2), text(19,130,'(b)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean') 
subplot(2,2,3), text(19,130,'(c)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean') 
subplot(2,2,4), text(19,130,'(d)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean')
  
  motion12=(xpos12(12:end)+xpos22(12:end)+ypos32(12:end))/3;
  motion22=(ypos12(12:end)+ypos22(12:end)+xpos32(12:end))/3;
  PCAmot12=(ff2(12:end,1)+ff2(12:end,3)+ff2(12:end,6))/3;
  PCAmot22=(ff2(12:end,2)+ff2(12:end,4)+ff2(12:end,5))/3;
  tm2=1:length(motion12);
  
  figure(2)
  subplot(1,2,1)
  plot(tm2,motion12, 'k', tm2,motion22,'k--', 'Linewidth',[2])
  axis([0 320 -80 80])
  xlabel('Frames'), ylabel('Displacement From Mean')
  legend('Primary Movement', 'Secondary Movement', 'Location','SouthEast')
  subplot(1,2,2)
  plot(tm2,PCAmot12, 'k', tm2,PCAmot22,'k--', 'Linewidth',[2])
  axis([0 320 -80 80])
  xlabel('Frames'), ylabel('Displacement From Mean')
  legend('Primary Movement', 'Secondary Movement', 'Location','SouthEast')
  
  
