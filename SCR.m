function index=SCR(x,y)


% SCR-based variable selection
% Er-wei Bai and Changming Cheng


[N,L]=size(x);  % N is data length, L is the number of variables.


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
    aa=lasso([1 0;0 0],[index(i);0],'Lambda',N^(-0.4))';
    if aa(1)~=0
       in22(i)=1;
    else 
        in22(i)=0;
    end
end
index=find(in22~=0);
disp(['The contributing variables  are: ', num2str(index)])




