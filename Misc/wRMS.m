function wrms=wRMS(simulations,samples)
wrms = sqrt(sum(bsxfun(@rdivide,bsxfun(@minus,simulations,samples.y'),samples.sigma').^2,2)/numel(samples.y));
end