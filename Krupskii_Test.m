function[Pval,G,Gb]=Krupskii_Test(N,D,NBoot,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Krupskii Radial Symmetry test 
% Based on Krupskii(2016)
% by Lorenzo Frattarolo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N number of observations
% D dimension of the random vector
% NBoot number of bootstrap replicates
% X data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% Pval Pvalues
% G Test Statistic
% Gb Test Statistic multiplier bootstrap replications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=4;
mkd= factorial(k)/factorial(k+D)*(-D)^D;
U=(tiedrank(X)-0.5)/N;
u=1/2- mean(U,2);
g=gkd(u,k,D);
G=sqrt(N)*mkd.*mean(g);
Gb(1,1)=G;
[bootstat,bootind ]= bootstrp(NBoot, [],X);
for b=1:NBoot
 Xb=X(bootind(:,b),:);   
for d=1:D 
U(:,d)=(Fhat(Xb(:,d),Xb(:,d))*(N+1)/N)-0.5/N;
end
u=1/2- mean(U,2);
g=gkd(u,k,D);
Gb(b,1)=sqrt(N)*mkd.*mean(g);
end
stg=std(Gb); 
stat= abs(G)/stg;
Pval= 2*(1-normcdf(stat));
end
