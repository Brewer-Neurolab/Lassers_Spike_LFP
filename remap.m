function mappedVals=remap(data,inMin,inMax,outMin,outMax)

mappedVals=(data-inMin)*(outMax-outMin)/(inMax-inMin)+outMin;

end