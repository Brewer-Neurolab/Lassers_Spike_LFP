function [imq,h,largestSignificantP]=myFalseDiscoveryRate(rank,nTests,alpha,pval)

imq=(rank./nTests).*alpha;

largestSignificantP=max(pval(pval<imq));

h=pval<largestSignificantP;