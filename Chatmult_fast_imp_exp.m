function[C,indout,Cd]=Chatmult_fast_imp_exp(u,U1,chi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Copula and Its multiplier bootstrap replicates
% Using matlab implicit expansion save memory
% By Lorenzo Frattarolo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u Emirical Copula argument
% U1 pseudo observations
% chi multiplier bootstrap random variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C Empirical copula values at u
% indout Indicators at u
% Cd Univariate Marginal at u_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=size(U1,1);
l=size(u,1);
D=size(U1,2);
NBoot=size(chi,2)-1;
indout=true(l,T);
Cd=zeros(l,NBoot+1,D);
for d=1:D 
ind=(U1(:,d)'<=u(:,d));
Cd(:,:,d)=ind*chi;
indout= indout.*ind;
end
C = indout*chi;
end


 