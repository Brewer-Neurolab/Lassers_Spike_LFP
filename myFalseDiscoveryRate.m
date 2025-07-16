function [imq,h,largestSignificantP]=myFalseDiscoveryRate(rank,nTests,alpha,pval)

imq=(rank./nTests).*alpha;

largestSignificantP=max(pval(pval<imq));

if isempty(largestSignificantP)
    largestSignificantP=-1;
end

h=pval<largestSignificantP;