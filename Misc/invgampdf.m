 function p = invgampdf(x,mean,alpha)
    beta=mean*(alpha-1);
    p=beta^alpha/gamma(alpha)*x^(-alpha-1)*exp(-beta/x);
end