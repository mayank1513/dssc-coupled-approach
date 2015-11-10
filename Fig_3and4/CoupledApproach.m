clear; clc;  
%% Getting inputs
wavelength=400:0.25:800;    ToPlot=[460, 500, 528.75, 680]; % Entry in ToPlot must be and element of wavelen %nm
ToPlot=sort(ToPlot);            x_Step=1; %in nm
%[n, d, LayerSets,n_repeat,VariThickness, ActiveLayers,EffStruct] = GetInput('InputFiles\Fig1_TemplateForLayerStruct.xlsx',wavelength);
getInputs_;
Thita_i=0;%pi/12;%:pi/(128):pi/2;
d=x_Step*floor(d/x_Step);
%cherence='8*mean(wavelength)';% Change only the multiplication factor

h=6.626e-34; % Js Planck's constant
c=2.998e8; %m/s speed of light
q=1.602e-19; %Coulomb electric charge
%% 
n_repeatingSets=size(LayerSets);
total_Layers=0; for Ind0=1:n_repeatingSets(1) ;total_Layers=total_Layers+length(LayerSets(Ind0,1):LayerSets(Ind0,2))*LayerSets(Ind0,3); end
LStruct(total_Layers+2,2)=0; LStruct(1)=1; LStructInd=2; 
for Ind0=1:n_repeatingSets
    for Ind1=1:LayerSets(Ind0,3)
        LySet=LayerSets(Ind0,1):LayerSets(Ind0,2);
        for Ind2=1:length(LySet)
            LStruct(LStructInd,1)=LySet(Ind2)+1; % k_1 is defined only for layers and not for incoming and outgoing media
            LStruct(LStructInd,2)=Ind0;
            LStructInd=LStructInd+1;
        end
    end
end
LStruct(end,1)=length(d);

a=size(LStruct); maxLayer=a(1);
%% Specifically for finding Electric field
%x_Step=1; %Don't delete--this is passed as argument to other functions
x=x_Step/2:x_Step:sum(d(LStruct(1:end-1,1)));
E(length(x))=0; E1=E*0;EtoPlot(length(x),3*length(ToPlot))=0;
activeLayer=ActiveLayers(1);%find(layer(:,3)==1);  --------------------------------------
t_cumsum=cumsum(d(LStruct(:,1)));
%%
    AM15_data=load('AM15.txt'); %mW/cm^2/nm
    AM15=interp1(AM15_data(:,1), AM15_data(:,2), wavelength, 'linear', 'extrap');
effSize=size(EffStruct);
LHE(length(wavelength),effSize(2))=0; R=LHE*0;T=LHE*0;
PhotoGenerationRate(length(x),effSize(2))=0;

for effStructInd=1:effSize(2)
    E(length(x))=0; E1=E*0;
    effStruct=EffStruct(LStruct(2:end-1,1),effStructInd); effStruct(2:end+1)=effStruct; effStruct(1)=0;
    Bmax=sum(effStruct)+1;
PlotInd=effStructInd;
for ThitaInd=1:length(Thita_i)
   
for a=1:length(wavelength)
    k=2*pi*n(LStruct(:,1),a)/wavelength(a);
    Thita_t=asin(real(k(1))*sin(Thita_i(ThitaInd))./real(k));%--------------
    rt_TE=zeros(2,2,Bmax); rt_TM=rt_TE*0;
    layerInd=2; bInd=1;
    while(layerInd<maxLayer)%calculating reflection and refraction coefficient
        if(effStruct(layerInd)>0)%RT=[rFromOneSide rFromOtherSide; tFromOneSide tFromOtherSide]---------------------------------------------
            rs=((real(k(layerInd-1))*cos(Thita_t(layerInd-1))-real(k(layerInd))*cos(Thita_t(layerInd)))/(real(k(layerInd-1))*cos(Thita_t(layerInd-1))+real(k(layerInd))*cos(Thita_t(layerInd))))^2;
            rt_TE(:,:,bInd)=[repmat(rs,1,2); repmat(1-rs,1,2)];
            rp=((real(k(layerInd-1))*cos(Thita_t(layerInd-1))-real(k(layerInd))*cos(Thita_t(layerInd)))/(real(k(layerInd-1))*cos(Thita_t(layerInd))+real(k(layerInd))*cos(Thita_t(layerInd-1))))^2;
            rt_TM(:,:,bInd)=[repmat(rp,1,2); repmat(1-rp,1,2)];
            bInd=bInd+1; layerInd=layerInd+1;
        else
            layerIndStart=layerInd;
            while(layerInd<maxLayer&&effStruct(layerInd)==0)
                layerInd=layerInd+1;
            end%--------------------------------------------------------------
            IndTrf=unique(LStruct(layerIndStart:layerInd-1,2));
            repeatingSetsTrf=cell(length(IndTrf),2); LyInd=1;
            for IndTrf1=1:length(IndTrf)
                repeatingSetsTrf{IndTrf1,1}=LyInd:LyInd+length(LayerSets(IndTrf(IndTrf1),1):LayerSets(IndTrf(IndTrf1),2))-1;  LyInd=LyInd+length(LayerSets(IndTrf(IndTrf1),1):LayerSets(IndTrf(IndTrf1),2));
                repeatingSetsTrf{IndTrf1,2}=n_repeat(IndTrf(IndTrf1),1);
            end
            [rt_TE(:,:,bInd),rt_TM(:,:,bInd),ETrf_TE,ETrf_TM]=Trf2(wavelength(a),Thita_i(ThitaInd),x_Step,n(unique(LStruct(layerIndStart-1:layerInd,1)),a),d(unique(LStruct(layerIndStart-1:layerInd,1))),repeatingSetsTrf);%first layer thickness is immeterial, thus layerIndStart should never be 1
            Indices=find(x>=t_cumsum(layerIndStart-1)&x<=t_cumsum(layerInd-1));
            E(Indices)=ETrf_TE(:,1);
            E1(Indices)=ETrf_TE(:,2);
            bInd=bInd+1; layerInd=layerInd+1;%It's correct Mayank at 2:12PM/13Oct204: 
        end
    end
    if(layerInd<=maxLayer)
        warning('I am here');
        pause;
    end
    RT=zeros(2,2,Bmax); RT1=zeros(2,2,Bmax); 
    RT(:,:,1)=rt_TE(:,:,1); RT1(:,:,end)=rt_TE(:,:,end); 
    k1=k(effStruct==1); d1=d(LStruct(effStruct==1)); Thita_t1=Thita_t(effStruct==1);
    for bInd=2:Bmax
        Abs=exp(-2*imag(k1(bInd-1)*d1(bInd-1)/cos(Thita_t1(bInd-1)))); %absCoeff=2*imag(k)
        RT(1,1,bInd)=RT(1,1,bInd-1)+RT(2,1,bInd-1)*rt_TE(1,1,bInd)*RT(2,2,bInd-1)*Abs^2/(1-RT(1,2,bInd)*rt_TE(1,1,bInd)*Abs^2);
        RT(1,2,bInd)=rt_TE(1,2,bInd)+rt_TE(2,2,bInd)*RT(1,2,bInd-1)*rt_TE(2,1,bInd)*Abs^2/(1-RT(1,2,bInd-1)*rt_TE(1,1,bInd)*Abs^2);
        RT(2,1,bInd)=RT(2,1,bInd-1)*rt_TE(2,1,bInd)*Abs/(1-RT(1,2,bInd-1)*rt_TE(1,1,bInd)*Abs^2);
        RT(2,2,bInd)=RT(2,2,bInd-1)*rt_TE(2,2,bInd)*Abs/(1-RT(1,2,bInd-1)*rt_TE(1,1,bInd)*Abs^2);
        %% Exclusively for calculating field intensity distribution
        bInd1=Bmax-bInd+1; Abs=exp(-2*imag(k1(bInd1)*d1(bInd1)/cos(Thita_t1(bInd1))));
        RT1(1,1,bInd1)=rt_TE(1,1,bInd1)+rt_TE(2,1,bInd1)*RT1(1,1,bInd1+1)*rt_TE(2,2,bInd1)*Abs^2/(1-rt_TE(1,2,bInd1)*RT1(1,1,bInd1+1)*Abs^2);
        RT1(1,2,bInd1)=RT1(1,2,bInd1+1)+RT1(2,2,bInd1+1)*rt_TE(1,2,bInd1)*RT1(2,1,bInd1+1)*Abs^2/(1-rt_TE(1,2,bInd1)*RT1(1,1,bInd1+1)*Abs^2);
        RT1(2,1,bInd1)=rt_TE(2,1,bInd1)*RT1(2,1,bInd1+1)*Abs/(1-rt_TE(1,2,bInd1)*RT1(1,1,bInd1+1)*Abs^2);
        RT1(2,2,bInd1)=rt_TE(2,2,bInd1)*RT1(2,2,bInd1+1)*Abs/(1-rt_TE(1,2,bInd1)*RT1(1,1,bInd1+1)*Abs^2);            
    end
    %{
    if(abs(RT(:,:,end)-RT1(:,:,1))>0.0001) % Testing the correctness of Algo
    %else
        RT
        RT1
        pause
    end
    %}
    k1(2:end+1)=k1; k1(1)=k(1); k1(end+1)=k(end); d1(2:end+1)=d1; d1(1)=0; d1(end+1)=0; Thita_t1(2:end+1)=Thita_t1; Thita_t1(1)=Thita_t(1); Thita_t1(end+1)=Thita_t(end);
    C=zeros(2,Bmax+1); C(:,1)=[1; RT(1,1,end)];
    for cInd=2:Bmax
        C(1,cInd)=RT(2,1,cInd-1)/(1-RT(1,2,cInd-1)*RT1(1,1,cInd)*exp(-4*imag(k1(cInd-1)*d1(cInd-1)/cos(Thita_t1(cInd-1)))));
        C(2,cInd)=RT(2,1,cInd-1)*RT1(1,1,cInd)*exp(-2*imag(k1(cInd-1)*d1(cInd-1)/cos(Thita_t1(cInd-1))))/(1-RT(1,2,cInd-1)*RT1(1,1,cInd)*exp(-4*imag(k1(cInd-1)*d1(cInd-1)/cos(Thita_t1(cInd-1)))));
    end
    C(:,end)=[RT(2,1,end); 0];
    
    layerInd=2; cInd=2;
    while(layerInd<maxLayer)%calculating the field intensity distribution
        if(effStruct(layerInd)>0)
            Indices=find(x>=t_cumsum(layerInd-1)&x<=t_cumsum(layerInd));
            E(Indices)=C(1,cInd)*exp(-2*imag(k(layerInd))*(x(Indices)-t_cumsum(layerInd-1)));
            E1(Indices)=C(2,cInd)*exp(-2*imag(k(layerInd))*(t_cumsum(layerInd)-x(Indices)));
            cInd=cInd+1; layerInd=layerInd+1;
        else
            layerIndStart=layerInd;
            while(layerInd<maxLayer&&effStruct(layerInd)==0)
                layerInd=layerInd+1;
            end
            %% Lord Damodar's causeless mercy
            Indices=find(x>=t_cumsum(layerIndStart-1)&x<=t_cumsum(layerInd-1));
            E(Indices)=C(1,cInd-1)*exp(-2*imag(k1(cInd-1)*d1(cInd-1)/cos(Thita_t1(cInd-1))))*E(Indices);
            E1(Indices)=C(2,cInd)*exp(-2*imag(k1(cInd)*d1(cInd)/cos(Thita_t1(cInd))))*E1(Indices);
        end
    end
    if(length(wavelength)==1)
        wavelen_step=1;
    else
        wavelen_step=(wavelength(end)-wavelength(1))/(length(wavelength)-1);
    end
    Indices=find(x>=t_cumsum(activeLayer-1)&x<=t_cumsum(activeLayer));
    PhotoGenerationRate(Indices,effStructInd)=PhotoGenerationRate(Indices,effStructInd) + 4*pi*imag(n(activeLayer,a))*real(n(activeLayer,a))*AM15(a)*(E(Indices)+E1(Indices))'*wavelen_step/(h*c); %/(sec-m^3)absorption per unit volume
    % PhotoGenerationRate=PhotoGenerationRate*1e-6; %/(sec-cm^3)
    % Unit conversions---convert all to SI units
    % AM15---W/m^2/nm--nm balanced by lambda_step
    % Result is in sec/m^3=1e-6 sec/cm^3----overall multiplication factor
    % of 1e-7
    R(a,effStructInd)=RT(1,1,end);
    T(a,effStructInd)=RT(2,1,end);
    LHE(a,effStructInd)=1-R(a,effStructInd)-T(a,effStructInd);
    if(PlotInd<=effSize(2)*length(ToPlot)&&wavelength(a)==ToPlot(floor((PlotInd-1)/effSize(2)+1)))
        EtoPlot(:,PlotInd)=(E+E1); PlotInd=PlotInd+effSize(2);
    end
    wavelength(a)
    %{
    if isempty(find(imag(E)))
    else
        pause
    end
    %}
end
Jsc(effStructInd)=q*sum(PhotoGenerationRate(:,effStructInd))*x_Step*1e-9*1e6; % Amp/m^2
end
end %ThintaInd ends here

%% Plots
%close all %EField optimised for fig2_2.xlsx
plot_fig_3;
plot_fig_4;
%%
figure()
subplot(3,2,1)
%LHE(:,1:effSize(2))=LHE(:,effSize(2):-1:1);
plot1 = plot(wavelength,LHE);
%title('Light Harvesting Efficiency') 
    xlabel('Wavelength (nm)');
    ylabel('LHE');
set(plot1(1),'LineWidth',0.25,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461]);
set(plot1(2),'LineWidth',1,...
    'Color',[0.313725501298904 0.313725501298904 0.313725501298904]);
set(plot1(3),'MarkerSize',1,'Marker','diamond','LineWidth',1.5,...
    'Color',[0 0 0]);

subplot(3,2,2)
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

    [~, Pos]=max(d(LStruct(:,1)));
%EtoPlot(:,1:end)=EtoPlot(:,end:-1:1);
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

subplot(3,2,6)
Indices=find(x>=t_cumsum(activeLayer-1)&x<=t_cumsum(activeLayer));
plot(x(Indices),PhotoGenerationRate(Indices,:))
axislimit3=axis;
for i=1:maxLayer-3
        line([sum(d(LStruct(1:i,1))) sum(d(LStruct(1:i,1)))],[0 axislimit3(4)],'color',[0.83, 0.82, 0.78]);
       % text(round((t_cumsum(i)+t_cumsum(i-1))/2),0,layers{i},'HorizontalAlignment','center','VerticalAlignment','bottom')
end

%title('Generation Rate in Device') 
    xlabel('Position in Device (nm)');
    ylabel('Generation rate /(sec-cm^3)');
    %}
%{
figure()    
T(:,1:3)=T(:,3:-1:1);
plot(wavelength,T);
title('Transmittance') 
    xlabel('Wavelength (nm)');
    ylabel('Fraction of light transmitted (Transmittance)');
    %}