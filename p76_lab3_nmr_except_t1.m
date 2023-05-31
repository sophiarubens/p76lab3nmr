len=2500;

%% error propag., T1 fit type, and avg. RF field formulas
del_t2s_horiz= @(delMxy,Mxy,delM0,M0) (delMxy./abs(Mxy))+(delM0./abs(M0));
del_t2_horiz= @(delMxy,Mxy,delM0,M0) (1/2)*((delMxy./abs(Mxy))+(delM0./abs(M0)));
BRF= @(gamma,pulselen) pi/(2*gamma*pulselen);
delBRF= @(delgamma,gamma,delpulselen,pulselen) (pi/(2*gamma*pulselen))*((delgamma/gamma)+(delpulselen/pulselen));
% del_t1_horiz_linfit= @(delM0,M0,delMz,Mz) ((delM0)./(abs(M0-Mz)))+((delMz.*abs(M0))./(abs(Mz.*(M0-Mz))));
% del_t1_horiz_zct= @(deltau) deltau/log(2);

%% mineral oil
% % T_1^*
samplex=linspace(0,1.6,50);
sampley=118*ones(50,1);
mineral_t1s=[1.5,118;1.4,118;1.3,118;1.2,118;1.1,118;1.0,118;0.9,118;0.8,118;0.7,118;0.6,118;0.5,116;...
             0.4,116;0.3,116;0.2,114;0.15,112;0.1,104;0.075,96;0.05,86;0.03,74;0.02,68;0.01,62];
len_min_t1s=length(mineral_t1s);
unc_min_t1s_p=[0.1,0.1,0.1,0.1,0.1,0.1,0.03,0.03,0.03,0.03,0.03,0.01,0.01,0.01,0.01,0.01,0.003,0.003,0.001,1e-3,5e-4]; %in seconds
unc_min_t1s_m=ones(len_min_t1s,1);
figure
plot(mineral_t1s(:,1),mineral_t1s(:,2),'MarkerSize',5,'LineStyle','none')
hold on
plot(samplex,sampley,'LineWidth',1)
errorbar(mineral_t1s(:,1),mineral_t1s(:,2),-unc_min_t1s_m,unc_min_t1s_m,-unc_min_t1s_p,unc_min_t1s_p,'.','MarkerSize',10)
title('heavy mineral oil T_1^*')
xlabel('period P (s)')
ylabel('FID amplitude (mV)')
legend('','saturation amplitude','data points','Location','southeast')
hold off

% T_2^*
mineral_t2star_prelim=readmatrix('F0020CH1.CSV'); % heavy mineral oil
mineral_t2s_t_prelim=mineral_t2star_prelim(:,4);
mineral_t2s_m_prelim=mineral_t2star_prelim(:,5);
unc_mineral_t2s_t=1e-7;
unc_mineral_t2s_m=1e-4;
m_m0=0.12;
figure
subplot(1,2,1)
errorbar(mineral_t2s_t_prelim,mineral_t2s_m_prelim,-unc_mineral_t2s_m*ones(len,1),unc_mineral_t2s_m*ones(len,1),...
         -unc_mineral_t2s_t*ones(len,1),unc_mineral_t2s_t*ones(len,1),'.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([-0.15e-3,2.5e-3,0,0.14])
title('heavy mineral oil T_2^* raw data')
mineral_t2s_t_lin=NaN;
mineral_t2s_m_lin=NaN;
j=1;
for i=1:1:len
    if ((mineral_t2s_t_prelim(i)>=(5.6e-5)) && (mineral_t2s_t_prelim(i)<=0.3e-3) && ((mineral_t2s_m_prelim(i)>0))) 
        % include only after the pulse where the magnetization is physically meaningful (positive; it would be a problem to have 
        % negative or equal-to-zero values because it is necessary to divide by another positive number and then take the natural
        % log, which is only defined with rational reals for positive arguments)
        mineral_t2s_t_lin(j)=mineral_t2s_t_prelim(i);
        mineral_t2s_m_lin(j)=mineral_t2s_m_prelim(i);
        j=j+1;
    end
end
err_min_t2s_t=unc_mineral_t2s_t*ones(j-1,1);
err_min_t2s_m_lin=unc_mineral_t2s_m*ones(j-1,1);
err_min_t2s_m_linearized=del_t2s_horiz(err_min_t2s_m_lin',mineral_t2s_m_lin,unc_mineral_t2s_m,m_m0);
linearized_min_m_t2s=-log(mineral_t2s_m_lin/m_m0);
mineral_t2s_mdl=fitlm(linearized_min_m_t2s,mineral_t2s_t_lin);
hmin_T2S_int=mineral_t2s_mdl.Coefficients{1,1};
hmin_T2S=mineral_t2s_mdl.Coefficients{2,1};
hmin_T2S_unc=mineral_t2s_mdl.Coefficients{2,2};
hmin_T2S_rsq=mineral_t2s_mdl.Rsquared.Ordinary;
min_t2s_eqn=['y=',num2str(hmin_T2S,4),'x+',num2str(hmin_T2S_int,4),'; R-squared=',num2str(hmin_T2S_rsq,4)];
subplot(1,2,2)
plot(mineral_t2s_mdl)
hold on
errorbar(linearized_min_m_t2s,mineral_t2s_t_lin,-err_min_t2s_t,err_min_t2s_t,-err_min_t2s_m_linearized,err_min_t2s_m_linearized,'.')
title('linearized heavy mineral oil T_2^* data')
legend('Data','Fit','Confidence bounds','')
xlabel('-ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
text(linearized_min_m_t2s(3),mineral_t2s_t_lin(3),min_t2s_eqn,'Position',[0.2,3e-4])
hold off

% T_2
m_m0=0.12;
mineral_t2_prelim=readmatrix('F0039CH1.CSV');
mineral_t2_t_prelim=mineral_t2_prelim(:,4);
mineral_t2_m_prelim=mineral_t2_prelim(:,5);
unc_mineral_t2_t=1e-7;
unc_mineral_t2_m=1e-4;
m_t2=[6e-05,0.14;0.00502,0.13;0.01002,0.104;0.015,0.088;0.02002,0.078;0.02502,0.064;0.03,0.058;0.03502,0.05;0.04,0.046;0.04502,0.04]; 
%these points were read from graph of the raw data...pick out the points delineating the envelope (wsyiwyg axes)
figure
subplot(1,2,1)
errorbar(mineral_t2_t_prelim,mineral_t2_m_prelim,-unc_mineral_t2_m*ones(len,1),unc_mineral_t2_m*ones(len,1),...
         unc_mineral_t2_t*ones(len,1),unc_mineral_t2_t*ones(len,1),'.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([0,0.05,-0.02,0.16])
title('heavy mineral oil T_2 raw data')
m_t2_m_lin=(-1/2)*log(m_t2(:,2)./m_m0);
m_t2_mdl=fitlm(m_t2_m_lin,m_t2(:,1));
hmin_T2_int=m_t2_mdl.Coefficients{1,1};
hmin_T2=m_t2_mdl.Coefficients{2,1};
hmin_T2_unc=m_t2_mdl.Coefficients{2,2};
hmin_T2_rsq=m_t2_mdl.Rsquared.Ordinary;
min_t2_eqn=['y=',num2str(hmin_T2,4),'x+',num2str(hmin_T2_int,4),'; R-squared=',num2str(hmin_T2_rsq,4)];
unc_m_t2_t=unc_mineral_t2_t*ones(length(m_t2),1);
unc_m_t2_m_lin=del_t2_horiz(unc_mineral_t2_m,m_t2(:,2),unc_mineral_t2_m,m_m0);
subplot(1,2,2)
plot(m_t2_mdl)
hold on
errorbar(m_t2_m_lin,m_t2(:,1),-unc_m_t2_t,unc_m_t2_t,-unc_m_t2_m_lin,unc_m_t2_m_lin,'.')
legend('Data','Fit','Confidence bounds')
xlabel('-0.5*ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
title('mineral oil T_2')
text(m_t2_m_lin(3),m_t2(3,1),min_t2_eqn,'Position',[0,0.04])
hold off

% average applied RF field
piover2len_mineral=3.60e-6;
delpiover2len_mineral=5e-9; %0.005e-6
gamma_proton=2.675221900e8; %MAKE SURE IT SHOULDN'T BE THE ONE FOR AN ELECTRON
delgamma_proton=0.018; %0.000 000 000 18 e8 = 0.018
BRF_mineral=BRF(gamma_proton,piover2len_mineral);
del_BRF_mineral=delBRF(delgamma_proton,gamma_proton,delpiover2len_mineral,piover2len_mineral);

%% water
% % T_1^*
w_t1s=[1.0,50.0;2.0,74.0;3.0,92.0;4.0,100.0;5.0,106;6.0,110;7.0,112;8.0,118;9.0,118;10.0,118];
samplex=linspace(-3,15,50);
sampley=linspace(118,118,50);
len_w_t1s=length(w_t1s);
unc_w_t1s_p=[0.1,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5];
unc_w_t1s_m=0.5*ones(len_w_t1s);
figure
errorbar(w_t1s(:,1),w_t1s(:,2),-unc_w_t1s_m,unc_w_t1s_m,-unc_w_t1s_p,unc_w_t1s_p,'.','MarkerSize',10)
hold on
plot(samplex,sampley,'LineWidth',1)
title('distilled water T_1^*')
xlabel('period P (s)')
ylabel('FID amplitude (mV)')
legend('','saturation amplitude','data points','Location','southeast')
axis([0,13,40,130])
hold off

% % T_2^*
water_t2star_prelim=readmatrix('F0024CH1.CSV');
water_t2s_t_prelim=water_t2star_prelim(:,4);
water_t2s_m_prelim=water_t2star_prelim(:,5);
unc_w_t2s_t=1e-7;
unc_w_t2s_m=1e-4;
figure
subplot(1,2,1)
errorbar(water_t2s_t_prelim,water_t2s_m_prelim,-unc_w_t2s_m*ones(len,1),unc_w_t2s_m*ones(len,1),-unc_w_t2s_t*ones(len,1),...
         unc_w_t2s_t*ones(len,1))
xlabel('time (s)')
ylabel('FID amplitude (mV)')
axis([-0.6e-3,0.2e-2,-0.01,0.12])
title('water T_2^* raw data')
w_m0=0.114;
water_t2s_t_lin=NaN;
water_t2s_m_lin=NaN;
j=1;
for i=1:1:len
    if( (water_t2s_t_prelim(i)>5.6e-5) && (water_t2s_t_prelim(i)<0.28e-3) && (water_t2s_m_prelim(i)>0) )
        water_t2s_t_lin(j)=water_t2s_t_prelim(i);
        water_t2s_m_lin(j)=water_t2s_m_prelim(i);
        j=j+1;
    end
end
w_t2s_m_linearized=-log(water_t2s_m_lin/w_m0);
w_t2s_mdl=fitlm(w_t2s_m_linearized,water_t2s_t_lin);
err_w_t2s_t=unc_w_t2s_t*ones(j-1,1);
err_w_t2s_m_lin=unc_w_t2s_m*ones(j-1,1);
err_w_t2s_m_linearized=del_t2s_horiz(err_w_t2s_m_lin',water_t2s_m_lin,unc_w_t2s_m,w_m0);
w_T2S_int=w_t2s_mdl.Coefficients{1,1};
w_T2S=w_t2s_mdl.Coefficients{2,1};
w_T2S_unc=w_t2s_mdl.Coefficients{2,2};
w_T2S_rsq=w_t2s_mdl.Rsquared.Ordinary;
w_t2s_eqn=['y=',num2str(w_T2S,4),'x+',num2str(w_T2S_int),'; R-squared=',num2str(w_T2S_rsq,4)];
subplot(1,2,2)
plot(w_t2s_mdl)
hold on
errorbar(w_t2s_m_linearized,water_t2s_t_lin,-err_w_t2s_t,err_w_t2s_t,-err_w_t2s_m_linearized,err_w_t2s_m_linearized,'.')
legend('Data','Fit','Confidence bounds','')
xlabel('-ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
text(w_t2s_m_linearized(3),water_t2s_t_lin(3),w_t2s_eqn,'Position',[0.05,2.5e-4])
title('linearized water T_2^* data')

% % T_2
water_t2_prelim=readmatrix('F0023CH1.CSV');
water_t2_t_prelim=water_t2_prelim(:,4);
water_t2_m_prelim=water_t2_prelim(:,5);
unc_w_t2_t=1e-6;
unc_w_t2_m=1e-4;
figure
subplot(1,2,1)
errorbar(water_t2_t_prelim,water_t2_m_prelim,-unc_w_t2_m*ones(len,1),unc_w_t2_m*ones(len,1),-unc_w_t2_t*ones(len,1),...
         unc_w_t2_t*ones(len,1))
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([0,0.5,-0.02,0.12])
title('water T_2 raw data')
w_t2=[0.015,0.108; 0.03,0.1; 0.045,0.098; 0.06,0.096; 0.075,0.09; 0.09,0.084; 0.105,0.076; 0.12,0.078; 0.15,0.08;...
      0.165,0.074; 0.18,0.072; 0.195,0.07; 0.21,0.07; 0.225,0.064; 0.24,0.062; 0.255,0.06; 0.27,0.056; 0.285,0.06; 0.3,0.056;...
      0.315,0.056; 0.33,0.054; 0.345,0.056; 0.36,0.05; 0.375,0.048; 0.39,0.046; 0.405,0.046; 0.42,0.044; 0.435,0.044 ;0.45,0.042;...
      0.465,0.04;0.48,0.044;0.495,0.04]; % exclude 1st data point bc it is at t=0 and that goes to infinity if you try to plug it 
                                         % into the linearization formula...there are still enough other data points that it doesn't 
                                         % really matter...the corresp. raw data plot also has wysiwyg axes
w_t2_m_lin=(-1/2)*log(w_t2(:,2)./w_m0);
w_t2_mdl=fitlm(w_t2_m_lin,w_t2(:,1));
w_T2_int=w_t2_mdl.Coefficients{1,1};
w_T2=w_t2_mdl.Coefficients{2,1};
w_T2_unc=w_t2_mdl.Coefficients{2,2};
w_T2_rsq=w_t2_mdl.Rsquared.Ordinary;
w_t2_eqn=['y=',num2str(w_T2,4),'x+',num2str(w_T2_int,4),'; R-squared=',num2str(w_T2_rsq,4)];
unc_w_t2_t=unc_w_t2_t*ones(length(w_t2),1);
unc_w_t2_m_lin=del_t2_horiz(unc_w_t2_m,w_t2(:,2),unc_w_t2_m,w_m0);
subplot(1,2,2)
plot(w_t2_mdl)
hold on
errorbar(w_t2_m_lin,w_t2(:,1),-unc_w_t2_t,unc_w_t2_t,-unc_w_t2_m_lin,unc_w_t2_m_lin,'.')
xlabel('-0.5*ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
title('linearized water T_2 data')
legend('Data','Fit','Confidence Bounds','')
text(w_t2_m_lin(3),w_t2(3,1),w_t2_eqn,'Position',[0.05,0.4])
hold off

%% FC-70
% T_1^*
fc70_t1s=[0.850,60;1.000,60;2.000,60;0.400,55.2;0.300,52.8;0.200,46.4;0.100,34.4;0.050,21.6];
unc_fc70_t1s_t=[0.03,0.1,0.2,0.01,0.01,0.01,0.005,0.003];
unc_fc70_t1s_m=0.1*ones(length(fc70_t1s));
samplex=linspace(0,2.5);
sampley=linspace(60,60);
figure
errorbar(fc70_t1s(:,1),fc70_t1s(:,2),-unc_fc70_t1s_m,unc_fc70_t1s_m,-unc_fc70_t1s_t,unc_fc70_t1s_t,'.','MarkerSize',10)
hold on
plot(samplex,sampley,'LineWidth',1)
xlabel('time (s)')
ylabel('FID amplitude (mV)')
axis([0,2.5,0,70])
title('FC70 T_1^* data')
legend('data points','saturation amplitude')
hold off

% T_2^*
fc70_t2star_prelim=readmatrix('F0025CH1.CSV'); % FC70
fc70_t2s_t_prelim=fc70_t2star_prelim(:,4);
fc70_t2s_m_prelim=fc70_t2star_prelim(:,5);
unc_fc70_t2s_t=1e-7;
unc_fc70_t2s_m=1e-5;
figure
subplot(1,2,1)
errorbar(fc70_t2s_t_prelim,fc70_t2s_m_prelim,-unc_fc70_t2s_m*ones(len,1),unc_fc70_t2s_m*ones(len,1),-unc_fc70_t2s_t*ones(len,1),...
         unc_fc70_t2s_t*ones(len,1),'.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([-0.8e-3,4.7e-3,-0.01,0.065])
title('FC70 T_2^* raw data')
fc70m0=0.0608;
fc70_t2s_t_lin=NaN;
fc70_t2s_m_lin=NaN;
j=1;
for i=1:1:len
    if ((fc70_t2s_t_prelim(i)>=6.0e-5) && (fc70_t2s_t_prelim(i)<0.4e-3) && (fc70_t2s_m_prelim(i)>0))
        fc70_t2s_t_lin(j)=fc70_t2s_t_prelim(i);
        fc70_t2s_m_lin(j)=fc70_t2s_m_prelim(i);
        j=j+1;
    end
end
fc70_t2s_m_linearized=-log(fc70_t2s_m_lin/fc70m0);
fc70_t2s_mdl=fitlm(fc70_t2s_m_linearized,fc70_t2s_t_lin);
err_fc70_t2s_t=unc_fc70_t2s_t*ones(j-1,1);
err_fc70_t2s_m_lin=unc_fc70_t2s_m*ones(j-1,1);
err_fc70_t2s_m_linearized=del_t2s_horiz(err_fc70_t2s_m_lin',fc70_t2s_m_lin,unc_fc70_t2s_m,fc70m0);
fc70_T2S_int=fc70_t2s_mdl.Coefficients{1,1};
fc70_T2S=fc70_t2s_mdl.Coefficients{2,1};
fc70_T2S_unc=fc70_t2s_mdl.Coefficients{2,2};
fc70_T2S_rsq=fc70_t2s_mdl.Rsquared.Ordinary;
fc70_t2s_eqn=['y=',num2str(fc70_T2S,4),'x+',num2str(fc70_T2S_int),'; R-squared=',num2str(fc70_T2S_rsq,4)];
subplot(1,2,2)
plot(fc70_t2s_mdl)
hold on
errorbar(fc70_t2s_m_linearized,fc70_t2s_t_lin,-err_fc70_t2s_t,err_fc70_t2s_t,-err_fc70_t2s_m_linearized,err_fc70_t2s_m_linearized,'.')
title('linearized FC70 T_2^* data')
legend('Data','Fit','Confidence bounds','')
xlabel('-ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
text(fc70_t2s_m_linearized(3),fc70_t2s_t_lin(3),fc70_t2s_eqn,'Position',[-0.2,3.5e-4])
ylabel('time (s)')
hold off

% T_2
fc70_t2_prelim=readmatrix('F0026CH1.CSV');
fc70_t2_t_prelim=fc70_t2_prelim(:,4);
fc70_t2_m_prelim=fc70_t2_prelim(:,5);
unc_fc70_t2_t=1e-6;
unc_fc70_t2_m=1e-5;
figure
subplot(1,2,1)
plot(fc70_t2_t_prelim,fc70_t2_m_prelim,'b-')
hold on
errorbar(fc70_t2_t_prelim,fc70_t2_m_prelim,-unc_fc70_t2_m*ones(len,1),unc_fc70_t2_m*ones(len,1),-unc_fc70_t2_t*ones(len,1),...
         unc_fc70_t2_t*ones(len,1),'.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([-0.5e-4,5e-2,-1e-3,6e-2])
title('FC70 T_2 raw data')
hold off
fc70_t2=[6e-05,0.0584;0.00222,0.0576;0.0044,0.0536;0.00664,0.0456;0.00884,0.0416;0.01104,0.036;0.01322,0.032;0.01542,0.028;...
         0.01984,0.0224;0.02204,0.0216;0.02418,0.0176;0.0264,0.0168;0.02856,0.0152;0.0308,0.0152;0.033,0.0144;0.03524,0.0136;...
         0.03742,0.0128;0.03962,0.0128;0.04182,0.012;0.04408,0.0112;0.04624,0.012]; %the corresp. raw data plot also has wysiwyg axes
fc70_t2_m_lin=(-1/2)*log(fc70_t2(:,2)./fc70m0);
fc70_t2_mdl=fitlm(fc70_t2_m_lin,fc70_t2(:,1));
fc70_T2_int=fc70_t2_mdl.Coefficients{1,1};
fc70_T2=fc70_t2_mdl.Coefficients{2,1};
fc70_T2_unc=fc70_t2_mdl.Coefficients{2,2};
fc70_T2_rsq=fc70_t2_mdl.Rsquared.Ordinary;
fc70_t2_eqn=['y=',num2str(fc70_T2,4),'x+',num2str(fc70_T2_int,4),'; R-squared=',num2str(fc70_T2_rsq,4)];
unc_fc70_t2_t=unc_fc70_t2_t*ones(length(fc70_t2),1);
unc_fc70_t2_m_lin=del_t2_horiz(unc_fc70_t2_m,fc70_t2(:,2),unc_fc70_t2_m,fc70m0);
subplot(1,2,2)
plot(fc70_t2_mdl)
hold on
errorbar(fc70_t2_m_lin,fc70_t2(:,1),-unc_fc70_t2_t,unc_fc70_t2_t,-unc_fc70_t2_m_lin,unc_fc70_t2_m_lin,'.')
xlabel('-0.5*ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
title('linearized FC70 T_2 data')
legend('Data','Fit','Confidence Bounds','')
text(fc70_t2_m_lin(3),fc70_t2(3,1),fc70_t2_eqn,'Position',[0.01,0.04])
hold off


%% FC-43
% T_2
fc43_t2_prelim=readmatrix('F0028CH1.CSV');
fc43_t2_t_prelim=fc43_t2_prelim(:,4);
fc43_t2_m_prelim=fc43_t2_prelim(:,5);
unc_fc43_t2_t=1e-6;
unc_fc43_t2_m=1e-5;
fc43m0=0.0584;
figure
subplot(1,2,1)
plot(fc43_t2_t_prelim,fc43_t2_m_prelim,'b-')
hold on
errorbar(fc43_t2_t_prelim,fc43_t2_m_prelim,-unc_fc43_t2_m*ones(len,1),unc_fc43_t2_m*ones(len,1),-unc_fc43_t2_t*ones(len,1),...
         unc_fc43_t2_t*ones(len,1),'r.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
title('FC43 T_2 raw data')
axis([-0.01,0.043,-0.0093,0.065])
hold off
fc43_t2=[4e-05,0.0552;0.00118,0.0592;0.0024,0.0592;0.00362,0.0552;0.0048,0.056;0.006,0.0512;0.0072,0.0512;0.0084,0.0472;...
        0.00962,0.0456;0.01082,0.0424;0.012,0.0408;0.0132,0.0376;0.0144,0.036;0.0144,0.036;0.0156,0.0336;0.01678,0.0312;...
        0.018,0.0304;0.0192,0.0288;0.0204,0.0272;0.02162,0.0248;0.02278,0.024;0.024,0.0224]; %corresp. raw data plot has wysiwyg axes
fc43_t2_m_lin=(-1/2)*log(fc43_t2(:,2)./fc43m0);
fc43_t2_mdl=fitlm(fc43_t2_m_lin,fc43_t2(:,1));
fc43_T2_int=fc43_t2_mdl.Coefficients{1,1};
fc43_T2=fc43_t2_mdl.Coefficients{2,1};
fc43_T2_unc=fc43_t2_mdl.Coefficients{2,2};
fc43_T2_rsq=fc43_t2_mdl.Rsquared.Ordinary;
fc43_t2_eqn=['y=',num2str(fc43_T2,4),'x+',num2str(fc43_T2_int,4),'; R-squared=',num2str(fc43_T2_rsq,4)];
unc_fc43_t2_t=unc_fc43_t2_t*ones(length(fc43_t2),1);
unc_fc43_t2_m_lin=del_t2_horiz(unc_fc43_t2_m,fc43_t2(:,2),unc_fc43_t2_m,fc43m0);
subplot(1,2,2)
plot(fc43_t2_mdl)
hold on
errorbar(fc43_t2_m_lin,fc43_t2(:,1),-unc_fc43_t2_t,unc_fc43_t2_t,-unc_fc43_t2_m_lin,unc_fc43_t2_m_lin,'.')
xlabel('-0.5*ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
title('linearized FC43 T_2 data')
legend('Data','Fit','Confidence Bounds','')
text(fc43_t2_m_lin(3),fc43_t2(3,1),fc43_t2_eqn,'Position',[0.01,0.025])
hold off

%% glycerin
% T_2
gl_t2_prelim=readmatrix('F0030CH1.CSV');
gl_t2_t_prelim=gl_t2_prelim(:,4);
gl_t2_m_prelim=gl_t2_prelim(:,5);
unc_gl_t2_t=1e-6;
unc_gl_t2_m=1e-5;
figure
subplot(1,2,1)
plot(gl_t2_t_prelim,gl_t2_m_prelim)
hold on
errorbar(gl_t2_t_prelim,gl_t2_m_prelim,-unc_gl_t2_m*ones(len,1),unc_gl_t2_m*ones(len,1),-unc_gl_t2_t*ones(len,1),...
         unc_gl_t2_t*ones(len,1),'r.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([-0.005,0.0975,-0.009,0.065])
title('glycerin T_2 raw data')
hold off
glm0=0.0616;
gl_t2=[4e-05,0.0616;0.0044,0.0592;0.0088,0.0552;0.0132,0.0488;0.0176,0.0456;0.022,0.04;0.0264,0.0352;0.03084,0.0328;...
       0.0352,0.0304;0.03964,0.0272;0.044,0.0248;0.0484,0.0232;0.0528,0.0208;0.05724,0.0184;0.06164,0.0168;0.066,0.0168;...
       0.0704,0.0152;0.0748,0.0136;0.0792,0.0136;0.0836,0.0112;0.088,0.0112;0.0924,0.0096]; % corresp. raw data plot has wysiwyg axes
gl_t2_m_lin=(-1/2)*log(gl_t2(:,2)./glm0);
gl_t2_mdl=fitlm(gl_t2_m_lin,gl_t2(:,1));
gl_T2_int=gl_t2_mdl.Coefficients{1,1};
gl_T2=gl_t2_mdl.Coefficients{2,1};
gl_T2_unc=gl_t2_mdl.Coefficients{2,2};
gl_T2_rsq=gl_t2_mdl.Rsquared.Ordinary;
gl_t2_eqn=['y=',num2str(gl_T2,4),'x+',num2str(gl_T2_int,4),'; R-squared=',num2str(gl_T2_rsq,4)];
unc_gl_t2_t=unc_gl_t2_t*ones(length(gl_t2),1);
unc_gl_t2_m_lin=del_t2_horiz(unc_gl_t2_m,gl_t2(:,2),unc_gl_t2_m,glm0);
subplot(1,2,2)
plot(gl_t2_mdl)
hold on
errorbar(gl_t2_m_lin,gl_t2(:,1),-unc_gl_t2_t,unc_gl_t2_t,-unc_gl_t2_m_lin,unc_gl_t2_m_lin,'.')
xlabel('-0.5*ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
title('linearized glycerin T_2 data')
legend('Data','Fit','Confidence Bounds','')
text(gl_t2_m_lin(3),gl_t2(3,1),gl_t2_eqn,'Position',[0.01,0.09])
hold off

%% silicone ultra strong
% T_2
sus_t2_prelim=readmatrix('F0031CH1.CSV');
sus_t2_t_prelim=sus_t2_prelim(:,4);
sus_t2_m_prelim=sus_t2_prelim(:,5);
unc_sus_t2_t=1e-7;
unc_sus_t2_m=1e-5;
figure
subplot(1,2,1)
plot(sus_t2_t_prelim,sus_t2_m_prelim)
hold on
errorbar(sus_t2_t_prelim,sus_t2_m_prelim,-unc_sus_t2_m*ones(len,1),unc_sus_t2_m*ones(len,1),-unc_sus_t2_t*ones(len,1),...
         unc_sus_t2_t*ones(len,1),'r.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([-0.0025,0.008,-0.0005,0.02])
title('silicone ultra strong T_2 raw data')
hold off
susm0=0.018;
sus_t2=[5.2e-05,0.018;0.000988,0.0146;0.002012,0.0124;0.003008,0.0094;0.00402,0.008;0.004996,0.0074;0.006028,0.0068;0.00702,0.007]; 
%!!time axis has values e-3 -- at least the data tips tool still tells you the right thing
sus_t2_m_lin=(-1/2)*log(sus_t2(:,2)./susm0);
sus_t2_mdl=fitlm(sus_t2_m_lin,sus_t2(:,1));
sus_T2_int=sus_t2_mdl.Coefficients{1,1};
sus_T2=sus_t2_mdl.Coefficients{2,1};
sus_T2_unc=sus_t2_mdl.Coefficients{2,2};
sus_T2_rsq=sus_t2_mdl.Rsquared.Ordinary;
sus_t2_eqn=['y=',num2str(sus_T2,4),'x+',num2str(sus_T2_int,4),'; R-squared=',num2str(sus_T2_rsq,4)];
unc_sus_t2_t=unc_sus_t2_t*ones(length(sus_t2),1);
unc_sus_t2_m_lin=del_t2_horiz(unc_sus_t2_m,sus_t2(:,2),unc_sus_t2_m,susm0);
subplot(1,2,2)
plot(sus_t2_mdl)
hold on
errorbar(sus_t2_m_lin,sus_t2(:,1),-unc_sus_t2_t,unc_sus_t2_t,-unc_sus_t2_m_lin,unc_sus_t2_m_lin,'.')
xlabel('-0.5*ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
title('linearized silicone ultra strong T_2 data')
legend('Data','Fit','Confidence Bounds','')
text(sus_t2_m_lin(3),sus_t2(3,1),sus_t2_eqn,'Position',[-0.05,6e-3])
hold off

%% beeswax with solvent
% T_2
bws_t2_prelim=readmatrix('F0034CH1.CSV');
bws_t2_t_prelim=bws_t2_prelim(:,4);
bws_t2_m_prelim=bws_t2_prelim(:,5);
unc_bws_t2_t=1e-7;
unc_bws_t2_m=1e-5;
figure
subplot(1,2,1)
plot(bws_t2_t_prelim,bws_t2_m_prelim)
hold on
errorbar(bws_t2_t_prelim,bws_t2_m_prelim,-unc_bws_t2_m*ones(len,1),unc_bws_t2_m*ones(len,1),-unc_bws_t2_t*ones(len,1),...
         unc_bws_t2_t*ones(len,1),'r.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([-0.0005,0.01,-0.002,0.05])
title('beeswax with solvent T_2 raw data')
hold off
bwsm0=0.0472;
bws_t2=[4.4e-05,0.0472;0.001028,0.0472;0.002008,0.0456;0.003028,0.0448;0.004016,0.048;...
        0.005016,0.0424;0.00602,0.0424;0.007016,0.0416;0.008012,0.0416;0.009012,0.04];
bws_t2_m_lin=(-1/2)*log(bws_t2(:,2)./bwsm0);
bws_t2_mdl=fitlm(bws_t2_m_lin,bws_t2(:,1));
bws_T2_int=bws_t2_mdl.Coefficients{1,1};
bws_T2=bws_t2_mdl.Coefficients{2,1};
bws_T2_unc=bws_t2_mdl.Coefficients{2,2};
bws_T2_rsq=bws_t2_mdl.Rsquared.Ordinary;
bws_t2_eqn=['y=',num2str(bws_T2,4),'x+',num2str(bws_T2_int,4),'; R-squared=',num2str(bws_T2_rsq,4)];
unc_bws_t2_t=unc_bws_t2_t*ones(length(bws_t2),1);
unc_bws_t2_m_lin=del_t2_horiz(unc_bws_t2_m,bws_t2(:,2),unc_bws_t2_m,bwsm0);
subplot(1,2,2)
plot(bws_t2_mdl)
hold on
errorbar(bws_t2_m_lin,bws_t2(:,1),-unc_bws_t2_t,unc_bws_t2_t,-unc_bws_t2_m_lin,unc_bws_t2_m_lin,'.')
xlabel('-0.5*ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
title('linearized beeswax with solvent T_2 data')
legend('Data','Fit','Confidence Bounds','')
text(bws_t2_m_lin(3),bws_t2(3,1),bws_t2_eqn,'Position',[-0.005,1e-2])
hold off

%% buna-N high strength
% T_2
bhs_t2_prelim=readmatrix('F0032CH1.CSV');
bhs_t2_t_prelim=bhs_t2_prelim(:,4);
bhs_t2_m_prelim=bhs_t2_prelim(:,5);
unc_bhs_t2_t=1e-7;
unc_bhs_t2_m=1e-5;
figure
subplot(1,2,1)
plot(bhs_t2_t_prelim,bhs_t2_m_prelim)
hold on
errorbar(bhs_t2_t_prelim,bhs_t2_m_prelim,-unc_bhs_t2_m*ones(len,1),unc_bhs_t2_m*ones(len,1),-unc_bhs_t2_t*ones(len,1),...
         unc_bhs_t2_t*ones(len,1),'r.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([-0.00012,0.01,-0.0015,0.04])
title('buna-N high strength T_2 raw data')
hold off
bhsm0=0.0372;
bhs_t2=[4.8e-05,0.0372;0.000616,0.022;0.001212,0.0172;0.001812,0.0132;0.002436,0.0116;0.003008,0.0104;0.003612,0.01;...
        0.004208,0.0092;0.0048,0.0088;0.005404,0.008;0.006008,0.008;0.006612,0.0076;0.007216,0.008;0.007812,0.0072;...
        0.0084,0.0064;0.009032,0.0056;0.009632,0.006];
bhs_t2_m_lin=(-1/2)*log(bhs_t2(:,2)./bhsm0);
bhs_t2_mdl=fitlm(bhs_t2_m_lin,bhs_t2(:,1));
bhs_T2_int=bhs_t2_mdl.Coefficients{1,1};
bhs_T2=bhs_t2_mdl.Coefficients{2,1};
bhs_T2_unc=bhs_t2_mdl.Coefficients{2,2};
bhs_T2_rsq=bhs_t2_mdl.Rsquared.Ordinary;
bhs_t2_eqn=['y=',num2str(bhs_T2,4),'x+',num2str(bhs_T2_int,4),'; R-squared=',num2str(bhs_T2_rsq,4)];
unc_bhs_t2_t=unc_bhs_t2_t*ones(length(bhs_t2),1);
unc_bhs_t2_m_lin=del_t2_horiz(unc_bhs_t2_m,bhs_t2(:,2),unc_bhs_t2_m,bhsm0);
subplot(1,2,2)
plot(bhs_t2_mdl)
hold on
errorbar(bhs_t2_m_lin,bhs_t2(:,1),-unc_bhs_t2_t,unc_bhs_t2_t,-unc_bhs_t2_m_lin,unc_bhs_t2_m_lin,'.')
xlabel('-0.5*ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
title('linearized buna-N high strength T_2 data')
legend('Data','Fit','Confidence Bounds','')
text(bhs_t2_m_lin(3),bhs_t2(3,1),bhs_t2_eqn,'Position',[0.01,0.7e-2])
hold off

%% vaseline
% T_2
v_t2_prelim=readmatrix('F0035CH1.CSV');
v_t2_t_prelim=v_t2_prelim(:,4);
v_t2_m_prelim=v_t2_prelim(:,5);
unc_v_t2_t=1e-6;
unc_v_t2_m=1e-5;
figure
subplot(1,2,1)
plot(v_t2_t_prelim,v_t2_m_prelim)
hold on
errorbar(v_t2_t_prelim,v_t2_m_prelim,-unc_v_t2_m*ones(len,1),unc_v_t2_m*ones(len,1),-unc_v_t2_t*ones(len,1),unc_v_t2_t*ones(len,1),'r.')
xlabel('time (s)')
ylabel('M_{xy} (V)')
axis([-2e-5,0.025,-0.005,0.08])
title('vaseline T_2 raw data')
hold off
vm0=0.0704;
v_t2=[5e-05,0.0704;0.00201,0.068;0.00302,0.0648;0.00502,0.0592;0.00603,0.056;0.00701,0.0544;0.00801,0.0512;0.00902,0.048;...
      0.01002,0.0488;0.01102,0.0472;0.01203,0.0448;0.01301,0.0448;0.01403,0.0424;0.01502,0.0416;0.01602,0.0416;0.01701,0.0392;...
      0.01801,0.0392;0.01902,0.0376;0.02003,0.036;0.02101,0.0376;0.02201,0.036;0.02301,0.0352;0.02402,0.0328];
v_t2_m_lin=(-1/2)*log(v_t2(:,2)./vm0);
v_t2_mdl=fitlm(v_t2_m_lin,v_t2(:,1));
v_T2_int=v_t2_mdl.Coefficients{1,1};
v_T2=v_t2_mdl.Coefficients{2,1};
v_T2_unc=v_t2_mdl.Coefficients{2,2};
v_T2_rsq=v_t2_mdl.Rsquared.Ordinary;
v_t2_eqn=['y=',num2str(v_T2,4),'x+',num2str(v_T2_int,4),'; R-squared=',num2str(v_T2_rsq,4)];
unc_v_t2_t=unc_v_t2_t*ones(length(v_t2),1);
unc_v_t2_m_lin=del_t2_horiz(unc_v_t2_m,v_t2(:,2),unc_v_t2_m,vm0);
subplot(1,2,2)
plot(v_t2_mdl)
hold on
errorbar(v_t2_m_lin,v_t2(:,1),-unc_v_t2_t,unc_v_t2_t,-unc_v_t2_m_lin,unc_v_t2_m_lin,'.')
xlabel('-0.5*ln(M_{xy}/M_0) (a dimensionless quantity related to magnetization)')
ylabel('time (s)')
title('linearized vaseline T_2 data')
legend('Data','Fit','Confidence Bounds','')
text(v_t2_m_lin(3),v_t2(3,1),v_t2_eqn,'Position',[0.01,0.02])
hold off