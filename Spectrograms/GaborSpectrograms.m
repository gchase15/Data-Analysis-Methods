%% The Setup

load handel 
L=length(y)/Fs;
n=length(y);
t2=linspace(0,L,n+1); 
t=t2(1:n);
k=(2*pi/L)*[0:n/2 -n/2:-1];
ks=fftshift(k);

v2=y'/2; %y is the data that comes when loading handel

figure(1)
plot((1:length(v2))/Fs,v2); %Fs is the sampling rate that comes when loading handel 
xlabel('Time [sec]'); 
ylabel('Amplitude'); 
title('Signal of Interest, v(n)');

 p8 = audioplayer(v2,Fs); 
 playblocking(p8);
 
%% Gabor Transform with Gaussian Window
figure(4)
vgt_spec=[]; 
tslide=0:.15:9;
tau=200;

for j=1:length(tslide)
    g=exp(-tau*(t-tslide(j)).^2); % Filter
    vg=g.*v2; 
    vgt=fft(vg); 
    vgt_spec=[vgt_spec; 
    abs(fftshift(vgt))]; 
    subplot(3,1,1), plot(t,v2,'k',t,g,'r')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,2), plot(t,vg,'k')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,3), plot(ks/(2*pi),abs(fftshift(vgt))/max(abs(vgt))) 
    xlabel('Frequency [Hz]'); 
    ylabel('|FFT(v)|')
    axis([-2500 2500 0 1])
    drawnow
    pause(0.02)
end

figure(5)
pcolor(tslide,ks/(2*pi),vgt_spec.'), 
shading interp 
set(gca,'Ylim',[-0 2500],'Fontsize',[14]) 
colormap(hsv)
xlabel('Time [sec]'); 
ylabel('Frequency [Hz]'); 

%% STEP FILTERS

figure(8)
width=2500;
step=zeros(1,length(t));
mask=ones(1,2*width+1);
vstept_spec=[]; 
tslide=(width+1):250:(length(step)-width);
for j=1:length(tslide)
    step=zeros(1,length(t));
    step(tslide(j)-width:1:tslide(j)+width)=mask;
    vstep=step.*v2; 
    vstept=fft(vstep); 
    vstept_spec=[vstept_spec; 
    abs(fftshift(vstept))]; 
    subplot(3,1,1), plot(t,v2,'k',t,step,'r')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,2), plot(t,vstep,'k')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,3), plot(ks/(2*pi),abs(fftshift(vstept))/max(abs(vstept))) 
    xlabel('Frequency [Hz]'); 
    ylabel('|FFT(v)|')
    axis([-2500 2500 0 1])
    drawnow
end

tslide=linspace(0,L,length(tslide));
figure(9)
pcolor(tslide,ks/(2*pi),vstept_spec.'), 
shading interp 
set(gca,'Ylim',[0 2500],'Fontsize',[14]) 
colormap(hsv)
xlabel('Time [sec]'); 
ylabel('Frequency [Hz]'); 

%% Mexican Hat Filter

figure(11)
width=10;
vmhwt_spec=[]; 
tslide=0:0.05:9;
for j=1:length(tslide)
    mhw=(1-((t-tslide(j))*width).^2).*exp(-(width^2*(t-tslide(j)).^2)/2);
    vmhw=mhw.*v2; 
    vmhwt=fft(vmhw); 
    vmhwt_spec=[vmhwt_spec; 
    abs(fftshift(vmhwt))]; 
    subplot(3,1,1), plot(t,v2,'k',t,mhw,'r')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,2), plot(t,vmhw,'k')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,3), plot(ks/(2*pi),abs(fftshift(vmhwt))/max(abs(vmhwt))) 
    xlabel('Frequency [Hz]'); 
    ylabel('|FFT(v)|') 
    axis([-2500 2500 0 1])
    drawnow  
end

figure(12)
pcolor(tslide,ks/(2*pi),vmhwt_spec.'), 
shading interp 
set(gca,'Ylim',[0 2500],'Fontsize',[14]) 
colormap(hsv)
xlabel('Time [sec]'); 
ylabel('Frequency [Hz]'); 

%% Piano Setup

tr_piano=16; %record time in seconds 
yp=audioread('music1.wav');
yp=yp';
Fsp=length(yp)/tr_piano;
Lp=tr_piano;
np=length(yp);
tp2=linspace(0,Lp,np+1); 
tp=tp2(1:np);
kp=(2*pi/Lp)*[0:np/2-1 -np/2:-1];
kps=fftshift(kp);

 figure(13)
plot((1:length(yp))/Fsp,yp); 
xlabel('Time [sec]'); 
ylabel('Amplitude');
title('Mary had a little lamb (piano)'); 
drawnow 

   p8 = audioplayer(yp,Fsp);
   playblocking(p8);
%% Piano STEP FILTERS

figure(16)
width=7014;
step=zeros(1,length(tp));
mask=ones(1,2*width+1);
ypstept_spec=[]; 
tslide=(width+1):7000:(length(step)-width);

for j=1:length(tslide)
    step=zeros(1,length(tp));
    step(tslide(j)-width:1:tslide(j)+width)=mask;
    ypstep=step.*yp; 
    ypstept=fft(ypstep); 
    ypstept_spec=[ypstept_spec; 
    abs(fftshift(ypstept))]; 
    subplot(3,1,1), plot(tp,yp,'k',tp,step,'r')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,2), plot(tp,ypstep,'k')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,3), plot(kps/(2*pi),abs(fftshift(ypstept))/max(abs(ypstept))) 
    xlabel('Frequency [Hz]'); 
    ylabel('|FFT(v)|')
    axis([-1200 1200 0 1])
    drawnow
end

tslidep=linspace(0,Lp,length(tslide));
figure(17)
pcolor(tslidep,kps/(2*pi),ypstept_spec.'), 
shading interp 
set(gca,'Ylim',[0 1000],'Fontsize',[14]) 
colormap(hsv)
xlabel('Time [sec]'); 
ylabel('Frequency [Hz]'); 


%% Recorder

tr_rec=14; % record time in seconds
yrec=audioread('music2.wav');
yrec=yrec';
Fsrec=length(yrec)/tr_rec;
Lrec=tr_rec;
nrec=length(yrec);
trec2=linspace(0,Lrec,nrec+1); 
trec=trec2(1:nrec);
krec=(2*pi/Lrec)*[0:nrec/2-1 -nrec/2:-1];
ksrec=fftshift(krec);

figure(18) 
plot((1:length(yrec))/Fsrec,yrec); 
xlabel('Time [sec]'); 
ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
drawnow
 p8 = audioplayer(yrec,Fsrec); 
 playblocking(p8);
 
 %% recorder step function
 
figure(21)
width=7144;
step=zeros(1,length(trec));
mask=ones(1,2*width+1);
yrecstept_spec=[]; 
tslide=(width+1):8000:(length(step)-width);

for j=1:length(tslide)
    step=zeros(1,length(trec));
    step(tslide(j)-width:1:tslide(j)+width)=mask;
    yrecstep=step.*yrec; 
    yrecstept=fft(yrecstep); 
    yrecstept_spec=[yrecstept_spec; 
    abs(fftshift(yrecstept))]; 
    subplot(3,1,1), plot(trec,yrec,'k',trec,step,'r')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,2), plot(trec,yrecstep,'k')
    xlabel('Time [sec]'); 
    ylabel('Amplitude'); 
    subplot(3,1,3), plot(ksrec/(2*pi),abs(fftshift(yrecstept))/max(abs(yrecstept))) 
    xlabel('Frequency [Hz]'); 
    ylabel('|FFT(v)|')
    axis([-1400 1400 0 1])
    drawnow
end

tsliderec=linspace(0,Lrec,length(tslide));
figure(22)
pcolor(tsliderec,ksrec/(2*pi),yrecstept_spec.'), 
shading interp 
set(gca,'Ylim',[600 2400],'Fontsize',[14]) 
colormap(hsv)
xlabel('Time [sec]'); 
ylabel('Frequency [Hz]'); 
   