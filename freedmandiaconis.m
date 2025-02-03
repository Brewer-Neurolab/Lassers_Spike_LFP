function nbins=freedmandiaconis(x)

Qx=iqr(x);

nbins=ceil((max(x)-min(x))/(2*Qx*(numel(x)^(-1/3))));

end