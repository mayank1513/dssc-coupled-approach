figure()
subplot(2,2,1)
%LHE(:,1:effSize(2))=LHE(:,effSize(2):-1:1);
plot1 = plot(wavelength,LHE);
axislimit3=axis;
for a=1:length(ToPlot)
        line([ToPlot(a) ToPlot(a)],[0 axislimit3(4)],'color',[0.99, 0.99, 0],'LineWidth',2.5);
       % text(round((t_cumsum(i)+t_cumsum(i-1))/2),0,layers{i},'HorizontalAlignment','center','VerticalAlignment','bottom')
end
plot2 = line(wavelength,LHE(:,2:3));
%title('Light Harvesting Efficiency') 
    xlabel('Wavelength (nm)');
    ylabel('LHE');
set(plot1(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot2(1),'LineWidth',1,...
    'Color',[0.313725501298904 0.313725501298904 0.313725501298904]);
set(plot2(2),'MarkerSize',1,'Marker','diamond','LineWidth',1.5,...
    'Color',[0 0 0]);

subplot(2,2,2)
%R(:,1:effSize(2))=R(:,effSize(2):-1:1);
plot2=plot(wavelength,R);
%title('Reflectance') 
    xlabel('Wavelength (nm)');
    ylabel('Reflectance');

    set(plot2(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot2(2),'LineWidth',1,...
    'Color',[0.313725501298904 0.313725501298904 0.313725501298904]);
set(plot2(3),'MarkerSize',1,'Marker','diamond','LineWidth',1.5,...
    'Color',[0 0 0]);

%    [~, Pos]=max(d(LStruct(:,1)));
%EtoPlot(:,1:end)=EtoPlot(:,end:-1:1);
%{
for effStructInd=1:length(ToPlot)
    subplot(3,2,2+effStructInd);
    %EtoPlot(:,(i1-1)*3+2:i1*3)=EtoPlot(:,i1*3:-1:(i1-1)*3+2);
    %axes1 = axes('Parent',figure1);%'Position',[0.13 0.425692695214106 0.331199294532628 0.499307304785894]);
    %xlim(axes1,[0 3000]);%-------------------------------------------------------------Change this if the cell changes
    %ylim(axes1,[0 3]);
    %box(axes1,'on');
%hold(axes1,'all');
    plot3=plot(x(x<t_cumsum(Pos-1)+450),EtoPlot(x<t_cumsum(Pos-1)+450,(effStructInd-1)*effSize(2)+1:effStructInd*effSize(2)));
    axislimit3=axis;
    for i=1:Pos-1
        line([sum(d(LStruct(1:i,1))) sum(d(LStruct(1:i,1)))],[0 axislimit3(4)],'color',[0.83, 0.82, 0.78]);
    end
    %title(['Normalized intensity variation for ', num2str(ToPlot(end-effStructInd+1)) ,'nm wavelength']) 
    xlabel('Position in Device (nm)');
    ylabel('Normalized Intensity |E|^2');
    set(plot3(1),...
    'Color',[0.501960813999176 0.501960813999176 0.501960813999176]);
set(plot3(3),'LineWidth',2);
%{    
    % Create textbox
annotation('textbox',...
    [0.135003876578134 0.861828031846892 0.0303531353135313 0.0577889447236182],...
    'String',{'TCO'},...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.207349206349206 0.844246421021987 0.110111111111111 0.0753768844221106],...
    'String',{'Working Electrode'},...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.354615520282187 0.856821892839513 0.0298641975308642 0.0628140703517589],...
    'String',{'PC'},...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.396943562610229 0.861460957178841 0.0536737213403881 0.0579345088161209],...
    'String',{'Electrolyte'},...
    'FitBoxToText','off',...
    'EdgeColor','none');
%}
end
%}

subplot(2,2,3)
%R(:,1:effSize(2))=R(:,effSize(2):-1:1);
plot3=plot(wavelength,T);
%title('Reflectance') 
    xlabel('Wavelength (nm)');
    ylabel('Transmittance');

    set(plot3(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot3(2),'LineWidth',1,...
    'Color',[0.313725501298904 0.313725501298904 0.313725501298904]);
set(plot3(3),'MarkerSize',1,'Marker','diamond','LineWidth',1.5,...
    'Color',[0 0 0]);


subplot(2,2,4)
Indices=find(x>=t_cumsum(activeLayer-1)&x<=t_cumsum(activeLayer));
plot2=plot(x(Indices),PhotoGenerationRate(Indices,:));
axislimit3=axis;
for a=activeLayer-1:maxLayer-3
        line([sum(d(LStruct(1:a,1))) sum(d(LStruct(1:a,1)))],[0 axislimit3(4)],'color',[0.83, 0.82, 0.78]);
       % text(round((t_cumsum(i)+t_cumsum(i-1))/2),0,layers{i},'HorizontalAlignment','center','VerticalAlignment','bottom')
end

    set(plot2(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot2(2),'LineWidth',1,...
    'Color',[0.313725501298904 0.313725501298904 0.313725501298904]);
set(plot2(3),'MarkerSize',1,'Marker','diamond','LineWidth',1.5,...
    'Color',[0 0 0]);

%title('Generation Rate in Device') 
    xlabel('Position in Device (nm)');
    ylabel('Generation rate /(sec-m^3)');
    %}