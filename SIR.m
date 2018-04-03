function index=SIR(x,y)

%  slice inverse regression based variable selection
% this algorithm can't detect the variable x_2, 
% because x_2^2  is a symmetric function, and the variable x_2 is also
% symmetric about 0.
% Er-wei Bai and Changming Cheng

[N,L]=size(x);  % N is data length, L is number of variables.

H=100;        % number of slices
var_x=cov(x);



[ysort,I]=sort(y);  % Reorder the data

for i=1:L
  xsort(:,i)=x(I,i);
end
for j=1:H
    x_bar(j,:)=mean(xsort((j-1)*N/H+1:j*N/H,:));
end
B=zeros(H,L,L);
for j=1:H
    A=zeros(L,L);
    for m=1:N/H
         A=A+(xsort((j-1)*N/H+m,:)-x_bar(j,:))'*(xsort((j-1)*N/H+m,:)-x_bar(j,:));
    end
    B(j,:,:)=A/(N/H-1);
end
C=sum(B,1)/H;

F=var_x-squeeze(C);  % Calculate F

EI=eye(L);
for k=1:L
    index(k)=EI(k,:)*F*EI(:,k);
end

index=abs(index)
for i=1:L
    aa=lasso([1 0;0 0],[index(i);0],'Lambda',N^(-0.4))';
    in22(i)=aa(1);
end

index=find(in22~=0);
disp(['The contributing variables  are: ', num2str(index)])


