%% Load Data
% clear all; close all; clc
addpath('HW3')
load cam1_1
load cam2_1 
load cam3_1
% load cam1_2
% load cam2_2 
% load cam3_2
% load cam1_3
% load cam2_3 
% load cam3_3
% load cam1_4
% load cam2_4 
% load cam3_4

%% Cutting the size of the video

%   Find the region that the paint can occupies during the entire video 
%   and cut each frame to include only that section. 
%   X positions move from left to right.
%   Y positions move from top to bottom. 

%Synchronize videos first by cutting initial frames on 2/3 videos
shortFrames1_1=vidFrames1_1(:,:,:,10:end);
shortFrames2_1=vidFrames2_1(:,:,:,19:end);
shortFrames3_1=vidFrames3_1(:,:,:,11:end);

%Find the minimum number of frames to make all videos the same length
F1=size(shortFrames1_1,4);
F2=size(shortFrames2_1,4);
F3=size(shortFrames3_1,4);
FrameLength=[F1 F2 F3];
FrameLength=min(FrameLength);

%Pick the "box" that you want to cut out of each frame so you only have the
%movement of the paint can.  
 imshow(shortFrames1_1(:,:,:,1));
 [x1,y1]=ginput(2);
 imshow(shortFrames2_1(:,:,:,1));
 [x2,y2]=ginput(2);
 imshow(shortFrames3_1(:,:,:,1));
 [x3,y3]=ginput(2);  close
 
 %Shrink each frame to fit the box picked above and cut any frames after
 %FrameLength ends.
 cutFrames1_1=shortFrames1_1(min(y1):max(y1),min(x1):max(x1),:,1:FrameLength);
 cutFrames2_1=shortFrames2_1(min(y2):max(y2),min(x2):max(x2),:,1:FrameLength);
 cutFrames3_1=shortFrames3_1(min(y3):max(y3),min(x3):max(x3),:,1:FrameLength);
 
 %Confirm you got it right and enjoy the show
 implay(cutFrames1_1)
 implay(cutFrames2_1)
 implay(cutFrames3_1)

 %% Track Video 1
 FrameLength=217;
 %First Part of For loop converts each Frame into grayscale. It preserves a
 %copy of this image as FrameGray1 for comparison later. 
 xpos1=[];
 ypos1=[];
 for j=1:FrameLength
 Frame=cutFrames1_1(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

%The for and if loops below filter out colors other than bright white. 
%The Detection threshold is set to 253 to guard against the rare case where 
%the light's color isn't shown at 255. 
 for k=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(k,i) > 248;
             FrameGray(k,i)=255;
         else FrameGray(k,i)=0;
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
 xpos1=[xpos1;
     xlight];
 ypos1=[ypos1;
     ylight];
 end
  xpos1=xpos1-mean(xpos1);
  ypos1=ypos1-mean(ypos1);
   %% Track Video 2
 
 xpos2=[];
 ypos2=[];
  for j=1:FrameLength
 Frame=cutFrames2_1(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

 for k=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(k,i) > 244;
             FrameGray(k,i)=255;
         else FrameGray(k,i)=0;
         end
     end
 end
%  figure(1)
%  imshow(FrameGray)
%  figure(2)
%  imshow(FrameGray1)

 [ylight, xlight]=ind2sub(size(FrameGray),find(FrameGray,1));
 xpos2=[xpos2;
     xlight];
  ypos2=[ypos2;
     ylight];
  end
  xpos2=xpos2-mean(xpos2);
   ypos2=ypos2-mean(ypos2);
  %% Track Video 3
 
 xpos3=[];
 ypos3=[];
  for j=1:FrameLength
 Frame=cutFrames3_1(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

 for i=1:size(FrameGray,1);
     for k=1:size(FrameGray,2);
         if FrameGray(i,k) > 239;
             FrameGray(i,k)=255;
         else FrameGray(i,k)=0;
         end
     end
 end
%  figure(1)
%  imshow(FrameGray)
%  figure(2)
%  imshow(FrameGray1)

 [ylight, xlight]=ind2sub(size(FrameGray),find(FrameGray,1));
 xpos3=[xpos3;
     xlight];
 ypos3=[ypos3;
     ylight];
  end
 xpos3=xpos3-mean(xpos3);
 ypos3=ypos3-mean(ypos3);
 %% The SVD
 A1=[xpos1';
     ypos1';
     xpos2';
     ypos2';
     xpos3';
     ypos3'];
 
 [u,s,v]=svd(A1','econ');
 
 A1Cov=cov(A1');
 
 variance=(s.^2)./(size(s,1)-1);
 var=diag(variance);
 
 sig=diag(s);
 energy1=sig(1)/sum(sig) 
%   energy2=sig(2)/sum(sig) 
%   energy3=sig(3)/sum(sig) 
 eng=sum(sig(1:3))/sum(sig)
    
  %% Describe Motion in mean frame
    
    meanFrames=double(rgb2gray(cutFrames1_1(:,:,:,1)));
    for j=2:size(cutFrames1_1,4)
        meanFrames=meanFrames + double(rgb2gray(cutFrames1_1(:,:,:,j)));
    end
    meanFrames=uint8(meanFrames/size(cutFrames1_1,4));
     
    figure(25)
    imshow(meanFrames)
    hold on
    [xcent,ycent]=ginput(1);
    plot([xcent,xcent+100*v(1,1)],[ycent,ycent+100*v(2,1)],'g','linewidth',2)
    plot([xcent,xcent+100*v(1,2)],[ycent,ycent+100*v(2,2)],'m','linewidth',2)
    plot([xcent,xcent+100*v(1,3)],[ycent,ycent+100*v(2,3)],'k','linewidth',2)
    legend('mode 1','mode 2','mode 3','Location','SouthEast')   
    
  %% Plotting SVD
 x=linspace(1,217,217);
  
 figure(5)
 subplot(1,2,1)
 plot(sig,'ko','Linewidth',[1.5])
 axis([0 7 0 1.4*10^3])
 set(gca,'Fontsize',[13],'Xtick',[0 1 2 3 4 5 6])
 xlabel('singular values'), ylabel('\sigma')
 grid on
 text(6.1,1.25*10^3,'(a)','Fontsize',[13])
 subplot(1,2,2)
 plot(var,'ko','Linewidth',[1.5])
 axis([0 7 0 3.5*10^5])
 set(gca,'Fontsize',[13],'Xtick',[0 1 2 3 4 5 6]) 
 xlabel('singular values'), ylabel('\sigma^2')
 grid on
 text(6.1,3.16*10^5,'(b)','Fontsize',[13])

 figure(9)
 plot(x,u(:,3),'g-.',x,u(:,2),'r',x,u(:,1),'k','Linewidth',[2]) 
 set(gca,'Fontsize',[13])
 legend('mode 3','mode 2','mode 1','Location','NorthEast')
 axis([0 220 -.2 .15])
 xlabel('Frames'), ylabel('U')

  %% Reproduction
  
  t=1:length(x);
  figure(18)
  subplot(2,2,1)
  plot(t,xpos1,t,ypos1,t,xpos2,t,ypos2,t,xpos3,t,ypos3); 
   axis([0 230 -150 150])
  legend('x1','y1','x2','y2','x3','y3','Location','SouthEast')
  
  for j=1:3
  ff=u(:,1:j)*s(1:j,1:j)*v(:,1:j)'; % modal projections 
    subplot(2,2,j+1)
    plot(x,ff)
    legend('x1','y1','x2','y2','x3','y3','Location','SouthEast')
    axis([0 230 -150 150])
  end
subplot(2,2,1), text(19,130,'(a)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean')
subplot(2,2,2), text(19,130,'(b)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean') 
subplot(2,2,3), text(19,130,'(c)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean') 
subplot(2,2,4), text(19,130,'(d)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean')

  motion1=(xpos1(12:end)+xpos2(12:end)+ypos3(12:end))/3;
  motion2=(ypos1(12:end)+ypos2(12:end)+xpos3(12:end))/3;
  PCAmot1=(ff(12:end,1)+ff(12:end,3)+ff(12:end,6))/3;
  PCAmot2=(ff(12:end,2)+ff(12:end,4)+ff(12:end,5))/3;
  tm=1:length(motion1);
  
  figure(1)
  subplot(1,2,1)
  plot(tm,motion2, 'k', tm,motion1,'k--', 'Linewidth',[2])
  axis([0 220 -80 80])
  xlabel('Frames'), ylabel('Displacement From Mean')
  legend('Primary Movement', 'Secondary Movement', 'Location','SouthEast')
  subplot(1,2,2)
  plot(tm,PCAmot2, 'k', tm,PCAmot1,'k--', 'Linewidth',[2])
  axis([0 220 -80 80])
  xlabel('Frames'), ylabel('Displacement From Mean')
  legend('Primary Movement', 'Secondary Movement', 'Location','SouthEast')
  
    
