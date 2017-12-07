(* ::Package:: *)

clear 
clc
close all

% Contour-based variable selection
% Er-wei Bai and Changming Cheng

N=1000; % data length

L=8;  % number of variabels

lambda=zeros(100,8); 
jj=0;  %  the number of correctly detecting right variables 
for ii=1: 100   % 100 Monte Carlo simulations

  x=randn(N,L);   % generate input data
  x(:,1)=x (:,3)/2;
  
  % generate y
  for k=1:N
    y(k)=x(k,1)+x(k,2)*x(k,2)+x(k,3)+0.1*randn;
  end

  c=(max(y)-min(y))/15;   % choose c

  Hc=zeros(L,L);
  H_mid=zeros(L,L);
  
  % calculate  Hc 
 iii=0;
   for i=2:N
     for j=1:i-1
        if abs(y(i)-y(j))<=c
            H_mid=(x(j,:)-x(i,:))'*(x(j,:)-x(i,:));
            Hc=Hc+H_mid;
            iii=iii+1;
        end
     end
   end

var_x=cov(x);  

Hc=Hc/iii;

EI=eye(L);
for k=1:L
    index(k)=EI(k,:)*(2*var_x-Hc)*EI(:,k);
end
index=abs(index);
for i=1:L
    aa=lasso ([1 0;0 0],[index(i);0],'Lambda',N^(-0.4))';
    in22(i)=aa(1);
end

ind=find(in22~=0);
if length(ind)==3 & ind==[1,2,3]
    jj=jj+1;
end
lambda(ii,:)=in22;
end
jj

max1=max(lambda)
min1=min(lambda)




