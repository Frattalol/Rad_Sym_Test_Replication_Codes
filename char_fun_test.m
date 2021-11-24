function[Pval,R]=char_fun_test(n,d,nb,X,alphafun,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Characteristic function Radial Symmetry test 
% based on Bahraoui and Quessy [2017]
% By Lorenzo Frattarolo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n number of observations
% d dimension of the random vector
% nb number of bootstrap replicates
% X data
% alphfun function handle to the kernel function
% sigma Smoothing parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% Pval Pvalues
% R Test Statistic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xi=exprnd(1,nb,n);
Delta= xi./repmat(mean(xi,2),1,n)-ones(nb,n);
Delta0=1/sqrt(n)*ones(1,n);
U=tiedrank(X)/(n+1);
W=U-1/2;
for i=1:n 
a=repmat(W(i,:),n,1);
b=W;
[alpham,alpha1ratiom,alpha2ratiom]=alphafun(sigma*(a-b));
[alphap,alpha1ratiop,alpha2ratiop]=alphafun(sigma*(a+b));    
B(i,:)=Bomega_fast(alpham,alphap);
D0(i,:)=B(i,:);
D1(i,:,:)=Bprimeomega_fast(sigma,alpham,alphap,alpha1ratiom,alpha1ratiop);
II(i,:,:)= (repmat(W(i,:),n,1)<=W)-1/2;
for j=1:n 
D2(i,j,:,:)=Bsecondomega_fast(sigma,alpham(j,:),alpha1ratiom(j,:),alpha2ratiom(j,:),alphap(j,:),alpha1ratiop(j,:),alpha2ratiop(j,:));
end
end
Domega=D0;
for k=1:d
 Domega=Domega +  (1/n)*((II(:,:,k)*D1(:,:,k)+(II(:,:,k)*D1(:,:,k))'));
for kprime=1:d  
 Domega=Domega + (1/n^2)*(II(:,:,k)*D2(:,:,k,kprime)*II(:,:,kprime)');   
end
end
RB=diag(Delta*Domega*Delta')/n;
R=Delta0*D0*Delta0';
Pval=mean(RB>R);
end