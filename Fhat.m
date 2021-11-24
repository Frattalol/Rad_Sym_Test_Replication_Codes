function[U]=Fhat(x,X)
T=size(X,1);
l=size(x,1);
if l>0
for i=1:l
U(i,1)=sum(X<=x(i,1))/(T+1);
end
else
U=0;
end