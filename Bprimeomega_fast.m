function[Bprime]=Bprimeomega_fast(sigma,alpham,alphap,alpha1ratiom,alpha1ratiop)
%[alpham,alpha1ratiom,alpha2ratiom]=alphafun(sigma*(a-b));
%[alphap,alpha1ratiop,alpha2ratiop]=alphafun(sigma*(a+b));
Bprime=sigma/2*(alpha1ratiom.*prod(alpham,2)-alpha1ratiop.*prod(alphap,2));
end


