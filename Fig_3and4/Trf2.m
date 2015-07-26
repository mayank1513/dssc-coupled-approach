function [rt_TE,rt_TM,ETrf_TE,ETrf_TM]=Trf2(wavelen,Thita_i,x_step,n,d,repeatingSets)
size_r_set=size(repeatingSets);
repeatingSets_rev=cell(size_r_set);
LayerInd=1;
for a=1:size_r_set(1)
    repeatingSets_rev{a,1}=LayerInd:(LayerInd+length(repeatingSets{size_r_set(1)+1-a,1})-1); LayerInd=LayerInd+length(repeatingSets{size_r_set(1)+1-a,1});
    repeatingSets_rev{a,2}=repeatingSets{size_r_set(1)+1-a,2};
end
[rt_TE(1,2),rt_TE(2,2),rt_TM(1,2),rt_TM(2,2),ETrf_TE(:,2),ETrf_TM(:,2)]=Trf(wavelen,Thita_i,x_step,n(end:-1:1),d(end:-1:1),repeatingSets_rev);
[rt_TE(1,1),rt_TE(2,1),rt_TM(1,1),rt_TM(2,1),ETrf_TE(:,1),ETrf_TM(:,1)]=Trf(wavelen,Thita_i,x_step,n,d,repeatingSets);

ETrf_TE(:,2)=ETrf_TE(end:-1:1,2);
ETrf_TM(:,2)=ETrf_TM(end:-1:1,2);

%figure();
%plot(ETrf_TE)
end

function [R_TE,T_TE,R_TM,T_TM,E_TE,E_TM]=Trf(wavelen,Thita_i,x_step,n,d,repeatingSets)
warning('off','MATLAB:nearlySingularMatrix');% suppressing unnecessary warnings remove it before...
d(1)=0;
% The numbering of layers excludes incoming and outgoing medium. Thus, d and k for layer1 are d(2) and k(2) respectively and the number of layers recorded in repeating sets info is 1.
n_repeatingSets=size(repeatingSets);
% Constants
%Impedence0=376.7; %omega
%% Computing transfer matrices for each layer and then computing R, T, LHE
total_Layers=0; for Ind0=1:n_repeatingSets(1) ;total_Layers=total_Layers+length(repeatingSets{Ind0,1})*repeatingSets{Ind0,2}(1); end
m_Layer1_TE(2,2,length(d)-2)=0; m_Layer2_TE=m_Layer1_TE*0; m_Layer1_TM=m_Layer1_TE*0; m_Layer2_TM=m_Layer1_TE*0;% Preallocating Memory for better performance
LHE_TE(length(wavelen),length(Thita_i))=0; R_TE=LHE_TE*0; T_TE=LHE_TE*0; LHE_TM=LHE_TE*0; R_TM=LHE_TE*0; T_TM=LHE_TE*0; C_TE(2,total_Layers)=0;
    
LStruct(total_Layers+2)=0; LStruct(1)=1; LStructInd=2; 
for Ind0=1:n_repeatingSets
    for Ind1=1:repeatingSets{Ind0,2}(1)
        for Ind2=1:length(repeatingSets{Ind0,1})
            LStruct(LStructInd)=repeatingSets{Ind0,1}(Ind2)+1; % k_1 is defined only for layers and not for incoming and outgoing media
            LStructInd=LStructInd+1;
        end
    end
end
LStruct(end)=length(d); d_cumsum=cumsum(d(LStruct));
x=x_step/2:x_step:d_cumsum(end-1); E_TE(length(x))=0; E_TM=E_TE*0; %Phase_TE=E_TE*0; AbsDepth=x*0; AbsByLayer(length(d_cumsum),length(wavelen))=0;
for wavelenInd=1:length(wavelen)
    k=2*pi*n(:,wavelenInd)/wavelen(wavelenInd);
    for ThitaInd=1:length(Thita_i)
        Thita_t=asin(k(1)*sin(Thita_i(ThitaInd))./k); %asin(real(k(1))*sin(Thita_i(ThitaInd))./real(k)); %---------------note that these angles are complex but the result obtained is exactly same as the result postulated
        m_incomingMedium=[1 1; cos(Thita_i(ThitaInd))*k(1) -cos(Thita_i(ThitaInd))*k(1)];
        m_incomingMedium_TM=[cos(Thita_i(ThitaInd)) cos(Thita_i(ThitaInd)); k(1) -k(1)];
        for LayerInd=1:length(k)-2
            m_Layer1_TE(:,:,LayerInd)=[1 1; cos(Thita_t(LayerInd+1))*k(LayerInd+1) -cos(Thita_t(LayerInd+1))*k(LayerInd+1)];
            m_Layer1_TM(:,:,LayerInd)=[cos(Thita_t(LayerInd+1)) cos(Thita_t(LayerInd+1)); k(LayerInd+1) -k(LayerInd+1)];
            
            k1=k(LayerInd+1)*cos(Thita_t(LayerInd+1));%+1i*imag(k(LayerInd+1))/cos(Thita_t(LayerInd+1));
            m_Layer2_TE(:,:,LayerInd)=[exp(1i*k1*d(LayerInd+1)) exp(-1i*k1*d(LayerInd+1)); cos(Thita_t(LayerInd+1))*k(LayerInd+1)*exp(1i*k1*d(LayerInd+1)) -cos(Thita_t(LayerInd+1))*k(LayerInd+1)*exp(-1i*k1*d(LayerInd+1))];
            m_Layer2_TM(:,:,LayerInd)=[cos(Thita_t(LayerInd+1))*exp(1i*k1*d(LayerInd+1)) cos(Thita_t(LayerInd+1))*exp(-1i*k1*d(LayerInd+1)); k(LayerInd+1)*exp(1i*k1*d(LayerInd+1)) -k(LayerInd+1)*exp(-1i*k1*d(LayerInd+1))];
        end
        m_OutgoingMedium=[1 1; cos(Thita_t(end))*k(end) -cos(Thita_t(end))*k(end)];
        m_OutgoingMedium_TM=[cos(Thita_t(end)) cos(Thita_t(end)); k(end) -k(end)];
        % Computing system matrix
        msys_TE=eye(2);  msys_TM=eye(2);
        for a=1:n_repeatingSets(1)
            msys1=eye(2); msys1_TM=eye(2);
            for Ind=1:length(repeatingSets{a,1})
                msys1=msys1*m_Layer1_TE(:,:,repeatingSets{a,1}(Ind))/m_Layer2_TE(:,:,repeatingSets{a,1}(Ind));
                msys1_TM=msys1_TM*m_Layer1_TM(:,:,repeatingSets{a,1}(Ind))/m_Layer2_TM(:,:,repeatingSets{a,1}(Ind));
            end
            msys_TE=msys_TE*(msys1^repeatingSets{a,2}(1));
            msys_TM=msys_TM*(msys1_TM^repeatingSets{a,2}(1));
            %LayerInd=LayerInd+length(repeatingSets{a,1});
        end
        msys_TE=m_incomingMedium\msys_TE*m_OutgoingMedium; %very imp
        msys_TM=m_incomingMedium_TM\msys_TM*m_OutgoingMedium_TM;
        %}
        % Computing R, T and LHE
        r_TE=msys_TE(2,1)/msys_TE(1,1); R_TE(wavelenInd,ThitaInd)=abs(r_TE)^2;
        t_TE=1/msys_TE(1,1); T_TE(wavelenInd,ThitaInd)=real(k(end)*cos(Thita_t(end)))*abs(t_TE)^2/(real(k(1)*cos(Thita_i(ThitaInd))));%k(end)*abs(t)^2/k(1);
        LHE_TE(wavelenInd,ThitaInd)=1-R_TE(wavelenInd,ThitaInd)-T_TE(wavelenInd,ThitaInd);
        
        r_TM=msys_TM(2,1)/msys_TM(1,1); R_TM(wavelenInd,ThitaInd)=abs(r_TM)^2;
        t_TM=1/msys_TM(1,1); T_TM(wavelenInd,ThitaInd)=k(end)*cos(Thita_t(end))*abs(t_TM)^2/(k(1)*cos(Thita_i(ThitaInd)));%k(end)*abs(t)^2/k(1);
        LHE_TM(wavelenInd,ThitaInd)=1-R_TM(wavelenInd,ThitaInd)-T_TM(wavelenInd,ThitaInd);
        
	%% computing E Field
        % Computing constants
        C_TE(2,total_Layers)=0; C_TM=C_TE*0;
        C_TE(:,end)=m_Layer2_TE(:,:,end)\m_OutgoingMedium*[t_TE;0];
        C_TM(:,end)=m_Layer2_TM(:,:,end)\m_OutgoingMedium_TM*[t_TM;0];
        cInd=total_Layers-1;
        for Ind0=n_repeatingSets:-1:1
            for Ind1=1:repeatingSets{Ind0,2}(1)
                for Ind2=length(repeatingSets{Ind0,1}):-1:2
                    C_TE(:,cInd)=m_Layer2_TE(:,:,repeatingSets{Ind0,1}(Ind2-1))\m_Layer1_TE(:,:,repeatingSets{Ind0,1}(Ind2))*C_TE(:,cInd+1);
                    C_TM(:,cInd)=m_Layer2_TM(:,:,repeatingSets{Ind0,1}(Ind2-1))\m_Layer1_TM(:,:,repeatingSets{Ind0,1}(Ind2))*C_TM(:,cInd+1);
                    cInd=cInd-1;
                end
                if(Ind1==repeatingSets{Ind0,2}(1))
                    if Ind0<=1
                        %r_varify=m_incomingMedium\m_Layer1_TE(:,:,repeatingSets{Ind0,1}(1))*C_TE(:,cInd+1);
                        %cInd=cInd-1;
                    else
                        C_TE(:,cInd)=m_Layer2_TE(:,:,repeatingSets{Ind0-1,1}(end))\m_Layer1_TE(:,:,repeatingSets{Ind0,1}(1))*C_TE(:,cInd+1);
                        C_TM(:,cInd)=m_Layer2_TM(:,:,repeatingSets{Ind0-1,1}(end))\m_Layer1_TM(:,:,repeatingSets{Ind0,1}(1))*C_TM(:,cInd+1);
                        cInd=cInd-1;
                    end
                else
                    C_TE(:,cInd)=m_Layer2_TE(:,:,repeatingSets{Ind0,1}(end))\m_Layer1_TE(:,:,repeatingSets{Ind0,1}(1))*C_TE(:,cInd+1);
                    C_TM(:,cInd)=m_Layer2_TM(:,:,repeatingSets{Ind0,1}(end))\m_Layer1_TM(:,:,repeatingSets{Ind0,1}(1))*C_TM(:,cInd+1);
                    cInd=cInd-1;
                end
            end
        end
        
        k_1=k(LStruct(2:end-1)); Thita_t1=Thita_t(LStruct(2:end-1));
        
        for Ind0=1:total_Layers %computing E
            k1=k_1(Ind0)*cos(Thita_t1(Ind0));%+1i*imag(k_1(LayerInd+1))/cos(Thita_t1(Ind0));
            %k1y=k_1(Ind0)*sin(Thita_t1(Ind0));
            xInd=find(x>=d_cumsum(Ind0)&x<d_cumsum(Ind0+1));
            E_TE(xInd)=abs(C_TE(1,Ind0)*exp(1i*(k1*(x(xInd)-d_cumsum(Ind0))))+C_TE(2,Ind0)*exp(1i*(-k1*(x(xInd)-d_cumsum(Ind0))))).^2;
            E_TM(xInd)=abs(C_TM(1,Ind0)*exp(1i*(k1*(x(xInd)-d_cumsum(Ind0))))+C_TM(2,Ind0)*exp(1i*(-k1*(x(xInd)-d_cumsum(Ind0))))).^2;
            %Absorption=imag(2*k1)*E_TE(xInd,1)*real(k_1(Ind0)*wavelen(wavelenInd)/(2*pi));%*Impedence0;
            %AbsDepth(xInd)=AbsDepth(xInd)+Absorption'*(wavelen(end)-wavelen(1))/(length(wavelen)-1)*AM15(wavelenInd);
            %AbsByLayer(Ind0,wavelenInd)=sum(Absorption)*x_step;
        end
    end % Thita ends
    %PercentComplete=wavelenInd/length(wavelen)*100
end % wavelen loop ends
warning('on','MATLAB:nearlySingularMatrix');% undo suppressing unnecessary warnings
end