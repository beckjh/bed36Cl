function [ proposal,MH_ratio ] = NormalProposal(a,b)
    proposal=normrnd(a,b);
    MH_ratio=1;
end
