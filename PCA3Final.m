%% Load Data
% clear all; close all; clc
addpath('HW3')
% load cam1_1
% load cam2_1 
% load cam3_1
% load cam1_2
% load cam2_2 
% load cam3_2
load cam1_3
load cam2_3 
load cam3_3
% load cam1_4
% load cam2_4 
% load cam3_4
 implay(vidFrames1_3)
 implay(vidFrames2_3)
 implay(vidFrames3_3)

%% Cutting the size of the video

%   Find the region that the paint can occupies during the entire video 
%   and cut each frame to include only that section. 
%   X positions move from left to right.
%   Y positions move from top to bottom. 

%Synchronize videos first by cutting initial frames videos
shortFrames1_3=vidFrames1_3(:,:,:,20:end);
shortFrames2_3=vidFrames2_3(:,:,:,7:end);
shortFrames3_3=vidFrames3_3(:,:,:,13:end);

%Find the minimum number of frames to make all videos the same length
F1=size(shortFrames1_3,4);
F2=size(shortFrames2_3,4);
F3=size(shortFrames3_3,4);
FrameLength=[F1 F2 F3];
FrameLength=min(FrameLength);

%Pick the "box" that you want to cut out of each frame so you only have the
%movement of the paint can.  
 imshow(shortFrames1_3(:,:,:,1));
 [x1,y1]=ginput(2);
 imshow(shortFrames2_3(:,:,:,1));
 [x2,y2]=ginput(2);
 imshow(shortFrames3_3(:,:,:,1));
 [x3,y3]=ginput(2);  close
 close
 %Shrink each frame to fit the box picked above and cut any frames after
 %FrameLength ends.
 cutFrames1_3=shortFrames1_3(min(y1):max(y1),min(x1):max(x1),:,1:FrameLength);
 cutFrames2_3=shortFrames2_3(min(y2):max(y2),min(x2):max(x2),:,1:FrameLength);
 cutFrames3_3=shortFrames3_3(min(y3):max(y3),min(x3):max(x3),:,1:FrameLength);
 
 %Confirm you got it right and enjoy the show
 implay(cutFrames1_3)
 implay(cutFrames2_3)
 implay(cutFrames3_3)

 %% Track Video 1
 FrameLength=220;
 %First Part of For loop converts each Frame into grayscale. It preserves a
 %copy of this image as FrameGray1 for comparison later. 
 xpos13=[];
 ypos13=[];
 for j=1:FrameLength
 Frame=cutFrames1_3(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

%The for and if loops below filter out colors other than bright white. 
%The Detection threshold is set to 253 to guard against the rare case where 
%the light's color isn't shown at 255. 
 for j=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(j,i) > 226;
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

white=find(FrameGray);

 [ylight, xlight]=ind2sub(size(FrameGray),white(end));
 xpos13=[xpos13;
     xlight];
  ypos13=[ypos13;
     ylight];
  end
 xpos13=xpos13-mean(xpos13);
 ypos13=ypos13-mean(ypos13);
  
 %% Track Video 2
 
 xpos23=[];
 ypos23=[];
  for j=1:FrameLength
 Frame=cutFrames2_3(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

 for j=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(j,i) > 251;
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
 xpos23=[xpos23;
     xlight];
 ypos23=[ypos23;
     ylight];
   end
 xpos23=xpos23-mean(xpos23);
 ypos23=ypos23-mean(ypos23);
  %% Track Video 3
 
 xpos33=[];
 ypos33=[];
  for j=1:FrameLength
 Frame=cutFrames3_3(:,:,:,j);
 FrameGray=rgb2gray(Frame);
 FrameGray1=FrameGray;

 for j=1:size(FrameGray,1);
     for i=1:size(FrameGray,2);
         if FrameGray(j,i) > 241;
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
 xpos33=[xpos33;
     xlight];
  ypos33=[ypos33;
     ylight];
   end
  xpos33=xpos33-mean(xpos33);
  ypos33=ypos33-mean(ypos33);
 %% The SVD
 
 %create a matrix of the three video's x and y pos over time
 A3=[xpos13';
     ypos13';
     xpos23';
     ypos23';
     xpos33';
     ypos33'];
 
  [u3,s3,v3]=svd(A3','econ'); %U and V are rotations, s is the stretch
  A3Cov=cov(A3');
  variance3=(s3.^2)./(size(s3,1)-1);
  var3=diag(variance3);
  
  sig3=diag(s3); %Collects the sigma values from s into a vector
  energy1=sig3(1)/sum(sig3) %Proportion of first singular value amongst all of them
  energy2=sig3(2)/sum(sig3) 
  eng3=sum(sig3(1:3))/sum(sig3)
  
%%  Describe Motion in mean frame 
      meanFrames=double(rgb2gray(cutFrames1_3(:,:,:,1)));
    for j=2:size(cutFrames1_3,4)
        meanFrames=meanFrames + double(rgb2gray(cutFrames1_3(:,:,:,j)));
    end
     meanFrames=uint8(meanFrames/size(cutFrames1_3,4));
    figure(25)
    imshow(meanFrames)
    hold on
    [xcent,ycent]=ginput(1);
    plot([xcent,xcent+100*v3(1,1)],[ycent,ycent+100*v3(2,1)],'g','linewidth',2)
    plot([xcent,xcent+150*v3(1,2)],[ycent,ycent+150*v3(2,2)],'m','linewidth',2)
    plot([xcent,xcent+100*v3(1,3)],[ycent,ycent+100*v3(2,3)],'k','linewidth',2)
    legend('mode 1','mode 2','mode 3','Location','SouthEast')  
    
    %% Plotting SVD
  x3=linspace(1,220,220);
  
 figure(7)
 subplot(1,2,1)
 plot(sig3,'ko','Linewidth',[1.5])
 axis([0 7 0 1*10^3])
 set(gca,'Fontsize',[13],'Xtick',[0 1 2 3 4 5 6])
 xlabel('singular values'), ylabel('\sigma')
 grid on
 text(6.1,1.25*10^3,'(a)','Fontsize',[13])
 subplot(1,2,2)
 plot(var3,'ko','Linewidth',[1.5])
 axis([0 7 0 1.4*10^5])
 set(gca,'Fontsize',[13],'Xtick',[0 1 2 3 4 5 6]) 
 xlabel('singular values'), ylabel('\sigma^2')
 grid on
 text(6.1,3.16*10^5,'(b)','Fontsize',[13])

figure(11)
 plot(x3,u3(:,3),'g-.',x3,u3(:,2),'r',x3,u3(:,1),'k','Linewidth',[2]) 
 set(gca,'Fontsize',[13])
 legend('mode 3','mode 2','mode 1','Location','SouthEast')
 axis([0 225 -.18 .15])
 xlabel('Frames'), ylabel('U')
  
   %% Reproduction
  
  t=1:length(x3);
  figure(20)
  subplot(2,2,1)
  plot(t,xpos13,t,ypos13,t,xpos23,t,ypos23,t,xpos33,t,ypos33); 
   axis([0 225 -100 100])
  legend('x1','y1','x2','y2','x3','y3','Location','SouthEast')
  
  for j=1:3
  ff3=u3(:,1:j)*s3(1:j,1:j)*v3(:,1:j)'; % modal projections 
    subplot(2,2,j+1)
    plot(x3,ff3)
    legend('x1','y1','x2','y2','x3','y3','Location','SouthEast')
    axis([0 225 -100 100])
  end
  
subplot(2,2,1), text(19,80,'(a)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean')
subplot(2,2,2), text(19,80,'(b)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean') 
subplot(2,2,3), text(19,80,'(c)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean') 
subplot(2,2,4), text(19,80,'(d)','Fontsize',[10]), xlabel('Frames'), ylabel('Displacement From Mean')
  
  motion13=(xpos13(12:end)+xpos23(12:end)+ypos33(12:end))/3;
  motion23=(ypos13(12:end)+ypos23(12:end)+xpos33(12:end))/3;
  PCAmot13=(ff3(12:end,1)+ff3(12:end,3)+ff3(12:end,6))/3;
  PCAmot23=(ff3(12:end,2)+ff3(12:end,4)+ff3(12:end,5))/3;
  tm3=1:length(motion13);
  
  figure(3)
  subplot(1,2,1)
  plot(tm3,motion13, 'k', tm3,motion23,'k-.', 'Linewidth',[2])
  axis([0 225 -80 80])
  xlabel('Frames'), ylabel('Displacement From Mean')
  legend('Primary Movement', 'Secondary Movement', 'Location','SouthEast')
  subplot(1,2,2)
  plot(tm3,PCAmot13, 'k', tm3,PCAmot23,'k-.', 'Linewidth',[2])
  axis([0 225 -80 80])
  xlabel('Frames'), ylabel('Displacement From Mean')
  legend('Primary Movement', 'Secondary Movement', 'Location','SouthEast')
  
  
  
  
  
  