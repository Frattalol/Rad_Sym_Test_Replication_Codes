function[alpha,alpha1ratio,alpha2ratio]=alphan(a)
alpha=exp(-(a.^2)/2);
alpha1ratio=-a;
alpha2ratio=a.^2-1;
end