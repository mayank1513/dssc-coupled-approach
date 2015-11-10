%[n, d, LayerSets,n_repeat,VariThickness, ActiveLayers,EffStruct] = GetInput('InputFiles\Fig1_TemplateForLayerStruct.xlsx',wavelength);
layer_names={'incoming_medium' 'TCO' 'WorkingElectrode' 'SiO2' 'TiO2' 'Electrolyte' 'TCO' 'outgoing_medium'};
%{
n=[1.6 1.8 1.95 1.43 1.92 1.433 1.8 1.6]; freq_dependent=5; ActiveLayers=5;
d=[0 400 1500 95 80 5000 400 0]'; %nm
EffStruct=[[0 0 0 0 0 0 0 0]' [0 0 0 0 0 1 0 0]' [0 0 1 0 0 1 0 0]'];
%}
% for fig 3
n=[1 1.4 1.6 1.8 1.95 1.43 1.92 1.433 1.8 1.6]; freq_dependent=5; ActiveLayers=5;
d=[0 400 30000 400 1500 95 80 50000 400 0]'; %nm
EffStruct=[[0 0 0 0 0 0 0 0 0 0]' [0 0 1 0 0 0 0 1 0 0]' [0 0 1 0 1 0 0 1 0 0]'];
%}
LayerSets=[1 1 1; 2 2 1; 3 3 1; 4 4 1; 5 6 3; 7 7 1; 8 8 1]; n_repeat=[1 1 1; 1 1 1; 1 1 1; 1 1 1; 3 3 1; 1 1 1; 1 1 1];
%repeatingSets={[1 2] [1 1 1]; [3 4] [3 3 1]; [5 6] [1 1 1]};

% computing n --- 
n=repmat(n',1,length(wavelength));
% Dye parameters
B0=0.004; lambda0=538; d_lambda=64.16;

z=(wavelength-lambda0)/d_lambda; 
n(freq_dependent,:)=n(freq_dependent,:)+1i*B0*exp(1-z-exp(-z));