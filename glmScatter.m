function glmScatter(glmTblAll,glmx,glmy,subregion,alpha)

xlog10=[];
x=[];
ylog10=[];
y=[];
colorVar=string();
% varNames=string(glmTblAll.Properties.VariableNames);
glmTbl=glmTblAll(glmTblAll.source_reg==subregion,:);
coeff_rowNames=string(glmTblAll.mdl{1}.Coefficients.Properties.RowNames);
for nConnect=1:height(glmTbl)%find(table2array(glmTblAll(:,[1,2,3]))==table2array(unique_sources(nAxons,:)))
    if isstring(glmx)
        if glmx=="mdl"
            xlog10(nConnect)=-log10(glmTbl.mdlPVal(nConnect)+eps^20);
            x(nConnect)=glmTbl.mdlPVal(nConnect)+eps^20;
        end
    else
        xlog10(nConnect)=-log10(glmTbl.mdl{nConnect}.Coefficients.pValue(glmx)+eps^20);
        x(nConnect)=glmTbl.mdl{nConnect}.Coefficients.pValue(glmx)+eps^20;
    end
    if isstring(glmy)
        if glmy=="mdl"
            ylog10(nConnect)=-log10(glmTbl.mdlPVal(nConnect)+eps^20);
            y(nConnect)=glmTbl.mdlPVal(nConnect)+eps^20;
        end
    else
        ylog10(nConnect)=-log10(glmTbl.mdl{nConnect}.Coefficients.pValue(glmy)+eps^20);
        y(nConnect)=glmTbl.mdl{nConnect}.Coefficients.pValue(glmy)+eps^20;
    end
    colorVar(nConnect)=string(glmTbl.fi(nConnect))+glmTbl.source_elec(nConnect);
end

disp(subregion)
disp(xlog10'+" "+ylog10')

unique_sources=unique(colorVar);
unique_colors=distinguishable_colors(numel(unique_sources));
cVec=[];
for nCol=1:numel(unique_sources)
    cVec(find(colorVar==unique_sources(nCol)),:)=repmat(unique_colors(nCol,:),[numel(find(colorVar==unique_sources(nCol))),1]);
end

figure
hold on
for nAxons=1:length(unique_sources)
    % scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),20,cVec(colorVar==unique_sources(nAxons)),'filled')
    thisColor=cVec(colorVar==unique_sources(nAxons),:);
    scatter(xlog10(colorVar==unique_sources(nAxons)),ylog10(colorVar==unique_sources(nAxons)),40,"MarkerEdgeColor",thisColor(1,:),"LineWidth",3)
    % scatter(x(colorVar==unique_sources(nAxons)),y(colorVar==unique_sources(nAxons)),40,"MarkerEdgeColor",thisColor(1,:),"LineWidth",1.5)
end

hold off

if isstring(glmx)
    if glmx=="mdl"
        xlabel("-log10 p model")
    end
else
    xlabel("-log10 p "+coeff_rowNames(glmx))
end
if isstring(glmy)
    if glmy=="mdl"
        ylabel("-log10 p model")
    end
else
    ylabel("-log10 p "+coeff_rowNames(glmy))
end

ax=gca;
ax.XScale="log";
ax.YScale="log";

xticks(logspace(-10,10,21))
xticklabels("10^{"+[-10:10]+"}")

% ylim([0.001,50])
% xlim([0.001,500])

BonferroniP=alpha/height(glmTbl);

xline(-log10(BonferroniP),"LineWidth",1)
yline(-log10(BonferroniP),"LineWidth",1)

legend([unique_sources,'',''],'Location','eastoutside')
% ax.tick
% end
% hold off
nsig=sum(xlog10>-log10(BonferroniP) & ylog10>-log10(BonferroniP));
% nsig=sum(x<BonferroniP & y<BonferroniP);

title(subregion+" #Significant Connections: "+nsig+"/"+height(glmTbl))
set(gca,"FontSize",32)

ax.LineWidth=2;
ax.TickLength=[0.05,0.05];

axis square

% Debug routine
% figure
% hold on
% for nAxons=1:length(unique_sources)
%     % scatter(amp(colorVar==unique_sources(nAxons)),interaction(colorVar==unique_sources(nAxons)),20,cVec(colorVar==unique_sources(nAxons)),'filled')
%     thisColor=cVec(colorVar==unique_sources(nAxons),:);
%     % scatter(xlog10(colorVar==unique_sources(nAxons)),ylog10(colorVar==unique_sources(nAxons)),40,"MarkerEdgeColor",thisColor(1,:),"LineWidth",1.5)
%     scatter(x(colorVar==unique_sources(nAxons)),y(colorVar==unique_sources(nAxons)),40,"MarkerEdgeColor",thisColor(1,:),"LineWidth",1.5)
% end
% hold off
% if isstring(glmx)
%     if glmx=="mdl"
%         xlabel("p model")
%     end
% else
%     xlabel("p "+coeff_rowNames(glmx))
% end
% if isstring(glmy)
%     if glmy=="mdl"
%         ylabel("p model")
%     end
% else
%     ylabel("p "+coeff_rowNames(glmy))
% end
% 
% ax=gca;
% ax.XScale="log";
% ax.YScale="log";
% xline((BonferroniP))
% yline((BonferroniP))
% 
% legend([unique_sources,'',''],'Location','eastoutside')
% % ax.tick
% % end
% % hold off
% nsig=sum(xlog10>-log10(BonferroniP) & ylog10>-log10(BonferroniP));
% % nsig=sum(x<BonferroniP & y<BonferroniP);
% 
% title(subregion+" #Significant Connections: "+nsig+"/"+height(glmTbl))
% set(gca,"FontSize",20)
% 
% axis square

end