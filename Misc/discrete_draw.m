function draw=discrete_draw(p,N)
[~,draw]=histc(rand(1 , N), [0 cumsum(p)./sum(p)]);
end