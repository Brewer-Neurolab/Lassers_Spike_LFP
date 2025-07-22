function standardPlotParams_240422(myAxis, xscale, yscale, xmin, xmax, ymin, ymax)

box on
set(myAxis,'TickLength',[0.035,0.025])
myAxis.LineWidth=1.25;
set(myAxis,'FontSize',18)

if xscale=="log" && ~isempty(xmin) && ~isempty(xmax)
    xfloor=floor(xmin);
    xceil=ceil(xmax);
    nDecades=xceil-xfloor+1;
    myXticks=logspace(xfloor,xceil,nDecades);
    set(myAxis,'XTick',myXticks)
    %generate lables
    xticklabels=[];
    for i=1:nDecades
        % xticklabels=[xticklabels,{"10^"+{string(xticks(i))}}];
        xticklabels=[xticklabels,{string(myXticks(i))}];
    end
    set(myAxis,'XTickLabel',xticklabels)
    xlim([10^xmin,10^xmax])

    myAxis.XAxis.MinorTick="on";
    minorTicks=[];
    for i=1:length(myXticks)-1
        minorTicks=[minorTicks,linspace(myXticks(i),myXticks(i+1),10)];
    end
    myAxis.XAxis.MinorTickValues=unique(minorTicks,'stable');
    myAxis.XScale="log";
end

if yscale=="log" && ~isempty(ymin) && ~isempty(ymax)
    yfloor=floor(ymin);
    yceil=ceil(ymax);
    nDecades=yceil-yfloor+1;
    myYticks=logspace(floor(ymin),ceil(ymax),nDecades);
    set(myAxis,'YTick',myYticks)
    %generate lables
    yticklabels=[];
    for i=1:nDecades
        yticklabels=[yticklabels,{"10^"+string(myXticks(i))}];
    end
    set(myAxis,'YTickLabel',yticklabels)
    ylim([ymin,ymax])

    myAxis.YAxis.MinorTick="on";
    minorTicks=[];
    for i=1:length(myYticks)-1
        minorTicks=[minorTicks,linspace(myYticks(i),myYticks(i+1),10)];
    end
    myAxis.YAxis.MinorTickValues=unique(minorTicks,'stable');
    myAxis.YScale="log";
end

if xscale=="linear" && ~isempty(xmin) && ~isempty(xmax)
    xticks('auto')
    myAxis.XAxis.MinorTick="on";
    xlim([xmin,xmax])
end
if yscale=="linear" && ~isempty(ymin) && ~isempty(ymax)
    yticks('auto')
    myAxis.YAxis.MinorTick="on";
    ylim([ymin,ymax])
end

end