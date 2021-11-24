function[B]=Bomega_fast(alpham,alphap)
%alpham=alphafun(sigma*(a-b));
%alphap=alphafun(sigma*(a+b));
B=1/2*(prod(alpham,2)-prod(alphap,2));
end


