function[g]=gkd(u,k,d)
g=abs(u).^(k+d).*sign(u);
end