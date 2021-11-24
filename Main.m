%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Replication Script:
% Loop on the different folders contining simulations for the different
% copula families.
% Return a folder Pval with p-values for each different configuration
% of simulation in different .mat files
% By Lorenzo Frattarolo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
%Create directory where the pvalues in output are stored 
mkdir([pwd,'\Pval'])
%Number of Bootstap Relicates
nb=1000;
%Number of repetitions of the Experiment
nrep=1000;
%Reads file names and other charadteristics from the simulation folder
familyfolder={'simulationsarch','simulationsellip','simulationsskewt'};
fsimul= [1 2 2 3];
sqvec={'','','Squared',''};
for f= 1: size(fsimul,2) 
listing = dir(strcat(pwd,'/',familyfolder{fsimul(f)},'/*.mat'));
listcell=struct2cell(listing);
[date idate]=sort(cell2mat(listcell(4,:)));
ordlist=listing(idate);
% Smoothing Parameter in the Characteristic function test
sigmavec=[1];
% Loop on different model simulation files
for ifile=1:size(ordlist,1)
disp([familyfolder{fsimul(f)},' ',sqvec{f},'  ',num2str(f),' of ',num2str(size(fsimul,2))])    
disp([ordlist(ifile).name,'  ',num2str(ifile),' of ',num2str(size(ordlist,1))])
load(strcat(pwd,'/',familyfolder{fsimul(f)},'/',ordlist(ifile).name))
alphafunvec=[{'alphan'}];
% Loop on different number of observations
for n=[125,250,500]
%Get from simulation file an (nrep n) x d matrix    
if f==3
u=abs(1-2*u);
end
if f==4
u=U;
end
d=size(u,2);
uu=u(1:n*nrep,:);
%Reshape in a 3-D nrep X n x d array
X0=reshape(uu,nrep,n,d);
save X0.mat -v7.3;
clear X0 uu
mX0=matfile('X0.mat');    
% Loop on different kernels for characteristic function test (only one
% possibility is included)
for ialpha=1 
fun=eval(strcat('@',char(alphafunvec(ialpha)))); 
%Loop on different Smoothing Parameters  for characteristic function test
for isigma=1
sigma=sigmavec(1,isigma);     
%Loop on repetitions
for rep=1:nrep
%load simulation for repetition rep     
X=squeeze(mX0.X0(rep,:,:));
tic
% TESTS
[Pvalcf(rep,1),R(rep,1)]=char_fun_test(n,d,nb,X,fun,sigma);
if (ialpha==1)&&(isigma==1)
[PvalGN(rep,1),CVMGN(rep,1)]=Emp_Cop_Test_GN(n,d,nb,X);
[Pvalplain(rep,1),CVMplain(rep,1)]=Emp_Cop_Test_plain(n,d,nb,X);
[Pvalmid(rep,1),CVMmid(rep,1)]=Emp_Cop_Test_mid(n,d,nb,X);
[Pvalkrupskii(rep,1),Gkrupskii(rep,1)]=Krupskii_Test(n,d,nb,X);
end
t(rep)=toc;
end
%Output
outfile=strcat(pwd,'/Pval/Pvalstatcf_',sqvec{f},'_',char(alphafunvec(ialpha)),'_','_sigma_',num2str(sigma),'_n_',num2str(n),'_from_',ordlist(ifile).name);
save(outfile,'Pvalcf','R')
end
end
outfile=strcat(pwd,'/Pval/PvalstatEC_',sqvec{f},'_n_',num2str(n),'_from_',ordlist(ifile).name);
save(outfile,'PvalGN','CVMGN','Pvalplain','CVMplain','Pvalmid','CVMmid')
outfile=strcat(pwd,'/Pval/PvalstatKru_',sqvec{f},'_n_',num2str(n),'_from_',ordlist(ifile).name);
save(outfile,'Pvalkrupskii','Gkrupskii')
end
end
end
