function[Bsecond]=Bsecondomega_fast(sigma,alpham,alpha1ratiom,alpha2ratiom,alphap,alpha1ratiop,alpha2ratiop)
%[alpham,alpha1ratiom,alpha2ratiom]=alphafun(sigma*(a-b));
%[alphap,alpha1ratiop,alpha2ratiop]=alphafun(sigma*(a+b));
%d=length(a);
% for k=1:d
% for k1=1:d
% if k==k1
% Bsecond(k,k)= -sigma^2/2*(alpha2ratiom(k)*prod(alpham)+alpha2ratiop(k)*prod(alphap));
% else
% Bsecond(k,k1)=-sigma^2/2*(alpha1ratiom(k)*alpha1ratiom(k1)*prod(alpham)+(alpha1ratiop(k)*alpha1ratiop(k1)*prod(alphap)));    
% end
% end
% end
Bsecond=-sigma^2/2*((alpha1ratiom'*alpha1ratiom)*prod(alpham)+(alpha1ratiop'*alpha1ratiop)*prod(alphap));
diagBsecond=-sigma^2/2*(alpha2ratiom*prod(alpham)+alpha2ratiop*prod(alphap));
Bsecond=Bsecond-diag(diag(Bsecond))+diag(diagBsecond);
end


