function [result,w,U,S,V,threshold,sortSound1,sortSound2,sortSound3,sortSound4]...
    =Test1_trainer(sound1,sound2,sound3,sound4,feature) 

n1=length(sound1(1,:));
n2=length(sound2(1,:));
n3=length(sound3(1,:));
n4=length(sound4(1,:));


[U,S,V]=svd([sound1,sound2,sound3,sound4],0);

sounds = S*V';
U = U(:,1:feature);
so1 = sounds(1:feature,1:n1);
so2 = sounds(1:feature,n1+1:n1+n2);
so3 = sounds(1:feature,n1+n2+1:n1+n2+n3);
so4 = sounds(1:feature,n1+n2+n3+1:n1+n2+n3+n4);
soTot = sounds(1:feature,1:n1+n2+n3+n4);

m1 = mean(so1,2);
m2 = mean(so2,2);
m3 = mean(so3,2);
m4 = mean(so4,2);
mtot = mean(soTot,2);

Sw=0;
for i=1:n1
    Sw = Sw + (so1(:,i)-m1)*(so1(:,i)-m1)';
end
for i=1:n2
    Sw = Sw + (so2(:,i)-m2)*(so2(:,i)-m2)';
end
for i=1:n3
    Sw = Sw + (so3(:,i)-m3)*(so3(:,i)-m3)';
end
for i=1:n4
    Sw = Sw + (so4(:,i)-m4)*(so4(:,i)-m4)';
end

Sb = (n1*(m1-mtot)*(m1-mtot)');
Sb = Sb + (n2*(m2-mtot)*(m2-mtot)');
Sb = Sb + (n3*(m3-mtot)*(m3-mtot)');
Sb = Sb + (n4*(m4-mtot)*(m4-mtot)');

[V2,D2] = eig(Sb,Sw);  % linear discriminant analysis
dd = diag(D2);
[lambda,ind] = max(abs(dd));
w = V2(:,ind);
w = w/norm(w,2);

vsound1 = w'*so1;
vsound2 = w'*so2;
vsound3 = w'*so3;
vsound4 = w'*so4;

result = [vsound1,vsound2,vsound3,vsound4];

if mean(vsound1)>mean(vsound2) 
    w = -w;
    vsound1 = -vsound1;
    vsound2 = -vsound2;
end
% dog < threshold < cat

sortSound1 = sort(vsound1);
sortSound2 = sort(vsound2);
sortSound3 = sort(vsound3);
sortSound4 = sort(vsound4);

% figure(4)
% subplot(2,2,1)
% hist(sortsound1,30); hold on, plot([18.22 18.22],[0 10],'r')
% set(gca,'Xlim',[-200 200],'Ylim',[0 10],'Fontsize',[14]), title('dog')
% subplot(2,2,2)
% hist(sortsound2,30,'r'); hold on, plot([18.22 18.22],[0 10],'r')
% set(gca,'Xlim',[-200 200],'Ylim',[0 10],'Fontsize',[14]), title('cat')
    
t1 = length(sortSound1);
t2 = 1;
while sortSound1(t1)>sortSound2(t2)
    t1 = t1-1;
    t2 = t2+1;
end
threshold = (sortSound1(t1)+sortSound2(t2))/2;
% t1/length(sortdog)

