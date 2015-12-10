figure()
subplot(2,2,1)
%LHE(:,1:effSize(2))=LHE(:,effSize(2):-1:1);
plot1 = plot(wavelength,LHE1);
axislimit3=axis;
%title('Light Harvesting Efficiency') 
    xlabel('Wavelength (nm)');
    ylabel('LHE');
set(plot1(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot1(3),'LineWidth',1,...
    'Color',[0 0 0]);
set(plot1(2),'MarkerSize',1,'Marker','diamond','LineWidth',1.5,...
    'Color',[0 0 0]);

subplot(2,2,2)
%R(:,1:effSize(2))=R(:,effSize(2):-1:1);
plot2=plot(wavelength,R1);
%title('Reflectance') 
    xlabel('Wavelength (nm)');
    ylabel('Reflectance');

set(plot2(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot2(3),'LineWidth',1,...
    'Color',[0 0 0]);
set(plot2(2),'MarkerSize',1,'Marker','diamond','LineWidth',1.5,...
    'Color',[0 0 0]);

subplot(2,2,3)
%R(:,1:effSize(2))=R(:,effSize(2):-1:1);
plot3=plot(wavelength,LHE2);
%title('Reflectance') 
    xlabel('Wavelength (nm)');
    ylabel('LHE');

    set(plot3(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot3(2),'LineWidth',1,...
    'Color',[0 0 0]);
set(plot3(3),'MarkerSize',1,'Marker','diamond','LineWidth',1.5,...
    'Color',[0 0 0]);

subplot(2,2,4)
%R(:,1:effSize(2))=R(:,effSize(2):-1:1);
plot3=plot(wavelength,R2);
%title('Reflectance') 
    xlabel('Wavelength (nm)');
    ylabel('Reflectance');

    set(plot3(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot3(2),'LineWidth',1,...
    'Color',[0 0 0]);
set(plot3(3),'MarkerSize',1,'Marker','diamond','LineWidth',1.5,...
    'Color',[0 0 0]);