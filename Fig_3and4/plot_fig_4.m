figure()
[~, Pos]=max(d(LStruct(:,1)));
%EtoPlot(:,1:end)=EtoPlot(:,end:-1:1);
for PlotInd=1:length(ToPlot)
    subplot(2,2,PlotInd);
    Indices = x>t_cumsum(3)-400 & x<t_cumsum(Pos-1)+174;
    plot2=plot(x(Indices),EtoPlot(Indices,(PlotInd-1)*effSize(2)+1:PlotInd*effSize(2)));
    axislimit3=axis;
    for i=3:Pos-1
        line([sum(d(LStruct(1:i,1))) sum(d(LStruct(1:i,1)))],[0 axislimit3(4)],'color',[0.83, 0.82, 0.78]);
    end
    %title(['Normalized intensity variation for ', num2str(ToPlot(end-effStructInd+1)) ,'nm wavelength']) 
    xlabel('Position in Device (nm)');
    ylabel('Normalized Intensity |E|^2');
set(plot2(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot2(2),'LineWidth',1,...
    'Color','blue');%[0.313725501298904 0.313725501298904 0.313725501298904]);
set(plot2(3),'MarkerSize',1,'Marker','diamond','LineWidth',1,...
    'Color','red');%[0 0 0]);
%     set(plot3(1),...
%     'Color',[0.501960813999176 0.501960813999176 0.501960813999176]);
% set(plot3(3),'LineWidth',2);
end