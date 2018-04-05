function x = invgamrnd(mean,alpha)
    x = (gamrnd(alpha,1/(mean*(alpha-1))))^(-1);
end
       