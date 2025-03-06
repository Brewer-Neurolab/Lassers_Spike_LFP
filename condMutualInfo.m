function CMI=condMutualInfo(Pxyz,Pxz,Pyz,Px,Py,Pz)

if ~isnan(Pxyz)
    %MI Miller Madow adjustment for low data
    % dMI=(numel(Pxy(Pxy~=0))-1)^2./(2*N*log(2));

    % calculate shannon entropy px
    Hx=-sum(Px.*log2(Px+eps));
    %calculate shannon entropy py
    Hy=-sum(Py.*log2(Py+eps));
    %calculate shannon entropy pz
    Hz=-sum(Pz.*log2(Pz+eps));

    Hxz=-sum(sum(Pxz.*log2(Pxz+eps)));
    Hyz=-sum(sum(Pyz.*log2(Pyz+eps)));

    Hxyz=-sum(sum(sum(Pxyz.*log2(Pxyz+eps))));

    %calculate conditional mutual information
    CMI=Hxz+Hyz-Hxyz-Hz;
    
    % MI=MI-dMI;
else
    CMI=0;
end

end