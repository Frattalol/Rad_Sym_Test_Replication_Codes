function[Pval,CVM,CVMFD]=Emp_Cop_Test_plain(N,D,NBoot,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Copula Radial Symmetry test plain version
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
% CVM Cramer Von Mises Test Statistic
% CVMFD Cramer Von Mises Test Statistic ultiplier bootstrap replications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xi=exprnd(1,NBoot,N);
Delta= 1/sqrt(N)*(xi./repmat(mean(xi,2),1,N)-ones(NBoot,N))';
Delta0=1/sqrt(N)*ones(1,N)';
chitilde=[Delta0,Delta];
U=tiedrank(X)/(N+1);
 gridpoints=U;
 U2=U;
 U2surv=1-U;
[alphaas,indout,Cad]=Chatmult_fast_imp_exp(gridpoints,U2,chitilde);
alphaas=alphaas';
[alphabs,indout,Cbd]=Chatmult_fast_imp_exp(gridpoints,U2surv,chitilde);
alphabs=alphabs';
corraFD=sparse(zeros(size(alphaas)));
corrbFD=sparse(zeros(size(alphabs)));
for d1=1:D
[ID(d1).dCdasFD]=FDDi_fast_imp_exp(U2,gridpoints,d1);
[ID(d1).dCdbsFD]=FDDi_fast_imp_exp(U2surv,gridpoints,d1);
corraFD=corraFD+[zeros(1,size(gridpoints,1));repmat(sparse(ID(d1).dCdasFD'),NBoot,1)].*squeeze(Cad(:,:,d1))' ;
corrbFD=corrbFD+[zeros(1,size(gridpoints,1));repmat(sparse(ID(d1).dCdbsFD'),NBoot,1)].*squeeze(Cbd(:,:,d1))' ;
end
ECasFD= alphaas-corraFD;
ECbsFD= alphabs-corrbFD;
CVMFD=mean(((ECasFD-ECbsFD)).^2,2);
Pval=mean(CVMFD(2:end,:)>CVMFD(1,:));
CVM=CVMFD(1,:);
end
