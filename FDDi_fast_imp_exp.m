function[dC]=FDDi_fast_imp_exp(U,u,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Copula Finite Difference Derivative
% Using matlab implicit expansion save memory
% By Lorenzo Frattarolo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u Emirical Copula argument
% U pseudo observations
% i component of the random vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dC Empirical copula Derivative with respect to u_i  at u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=size(U,1);
D=size(U,2);
n=size(u,1);
h=N^(-1/2);
inddim=[1:D];
DD=inddim(inddim~=i);
indout=true(n,N);    
for d=DD
ind=U(:,d)'<=u(:,d);
indout= indout.*ind;
end
inddC=((U(:,i)'<=(u(:,i)+h))-(U(:,i)'<=(u(:,i)-h)))/((2*h)*N); 
dC= sum(indout.*inddC,2);
end
