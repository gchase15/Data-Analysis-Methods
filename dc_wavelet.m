function dcData = dc_wavelet(dcfile) 
[m,n]=size(dcfile); % 4096 x 80 
nw=32*32; % wavelet resolution

nbcol = size(colormap(gray),1);
for i=1:n
  X=double(reshape(dcfile(:,i),64, 64)); 
  [cA,cH,cV,cD]=dwt2(X,'haar');
  cod_cH1 = wcodemat(cH,nbcol);
  cod_cV1 = wcodemat(cV,nbcol); 
  cod_edge=cod_cH1+cod_cV1; 
  dcData(:,i)=reshape(cod_edge,nw,1);
end