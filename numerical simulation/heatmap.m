%% ******* original version: Maximilian Gram, University of Wuerzburg *******
% 04.02.2024
% this code reproduces Fig. 2 of doi.org/10.1002/mrm.28585
% I added the effect of B1+ field deviations for the exciatioan and
% refocusing pulses. Both pulses are simulated as instantaneous operators.
% Relaxation effects during excitation and refocsuing are neglected. 
% T1rho dispersion effects from deviating B1+ fields are also neglected. 
%% ******* revised version: Cai Wan, University Medical Center Freiburg *******
% last revised 04.02.2024
% reference doi.org/10.1002/mrm.28585 and doi: 10.1002/nbm.4834 
clc;
clear all; 
heng=0;
shu=0;
for bf0 = 1200:-24:0
    shu =shu+1;
  for BB1 = 0:2:100
      heng = heng+1;
df0   = (bf0-600); % [Hz] off-resonance 50mT: 2.323*10^6*130/10^6
dB1   = 1+(BB1-50)/100;                 % [ ] deviation of B1+ field, 1 -> optimal field
T1rho = 100 *1e-3;	         % [s] T1rho relaxation time
T2rho = 100	*1e-3;           % [s] T2rho relaxation time
fsl   = 500;                 % [Hz] SL amplitude

dw0   = 2*pi* df0;             % [rad/s]
wSL   = 2*pi* fsl;             % [rad/s]
weff  = sqrt((wSL)^2+(dw0)^2); % [rad/s] effective field
theta = atan2(dw0,wSL);        % [rad]   field tilting

M0 = [0; 0; 1];

tSL = (1 : 1 : 200) *1e-3; % [s] SL preparation time

%% reference
Mz_ref = exp(-tSL / T1rho);

%% simulate: composite spin-lock
Mz_CSL = zeros(numel(tSL),1);

for j=1:numel(tSL)
    % Supporting Information Table S1 a)
    P = R90_1(dB1) * ...
        SL_2(tSL(j)/2, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_1(tSL(j)/2, T1rho, T2rho, wSL, theta) * ...
        R90_1(dB1);
    M = P * M0;
    Mz_CSL(j) = abs(M(end));
end

%% simulate: balanced spin-lock
Mz_BSL = zeros(numel(tSL),1);

for j=1:numel(tSL)
    % Supporting Information Table S1 a)
    P = R90_2(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_2(tSL(j)/2, T1rho, T2rho, wSL, theta) * ...
        R180_2(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R90_1(dB1);
    M = P * M0;
    Mz_BSL(j) = M(end);
end

%% simulate: Paired self-compensated spin-lock
Mz_PSCSL = zeros(numel(tSL),1);

for j=1:numel(tSL)
    % Supporting Information Table S1 a)
    P = R90_2(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R90_1(dB1);
    M = P * M0;
    Mz_PSCSL(j) = M(end);
end

%% simulate: Triple-refocused spin-lock
Mz_TRSL = zeros(numel(tSL),1);

for j=1:numel(tSL)
    % Supporting Information Table S1 a)
    P = R90_2(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_2(dB1) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R90_2(dB1);
    M = P * M0;
    Mz_TRSL(j) = abs(M(end));
end

%% simulate: Multi-refocused spin-lock
Mz_MRSL = zeros(numel(tSL),1);

for j=1:numel(tSL)

    %%
       P = R90_2(dB1) * ...
        SL_1(tSL(j)/8, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_2(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_2(dB1) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_1(tSL(j)/8, T1rho, T2rho, wSL, theta) * ...
        R90_1(dB1);
    M = P * M0;
    Mz_MRSL(j) = abs(M(end));
end
%% plot results
% figure()
% hold on
% plot(tSL, Mz_ref, 'k--')
% plot(tSL, Mz_CSL, 'r-')
% plot(tSL, Mz_BSL, 'b-')
% plot(tSL, Mz_TRSL, 'k-')
% legend('ref', 'CSL', 'BSL','TRCSL')


%% monoexponential SL
% [pks,locs] = findpeaks(Mz_TRSL);  % max peak
% y_max = pks;
% x_max = locs*1e-3;
% fun1 = @(x,x_max)((x(1)*exp(-x_max/x(2))));
% 
% nMz_BSL=1-Mz_TRSL;                % min peak
% [pksL,locsL] = findpeaks(nMz_BSL);
% y_min = pksL;
% x_min = locsL*1e-3;
% fun2 = @(x,x_min)((x(1)*exp(-x_min/x(2))));
% 
% tSL = (1 : 1 : 200) *1e-3; % [s] SL preparation time
% fun = @(x,tSL)((x(1)*exp(-tSL/x(2))));  % global fit 
x0 = [0.8,90];
% [T1r1_mxa] = lsqcurvefit(fun2,x0,x_max,y_max);
% [T1r1_min] = lsqcurvefit(fun1,x0,x_min,1-y_min);
% [T1r1] = lsqcurvefit(fun,x0,tSL,Mz_TRSL');

% figure(2)
% hold on
% box on;
% grid on;  
% set(gca, 'GridLineStyle', ':'); 
% plot(tSL,Mz_ref,'r-');
% fill([(tSL),fliplr(tSL)],[(fun2(T1r1_mxa,tSL)),fliplr(fun1(T1r1_min,tSL))],[255/255 204/255 255/255],'EdgeColor','none');
% plot(tSL,Mz_TRSL,'k-',tSL,fun(T1r1,tSL),'b--')
% legend('monoexponential reference','area of min/max T_1_\rho quantification','trajectory M_z(T_S_L)','global fit');
% xlabel('Spin lock time T_S_L [s]');
% ylabel('Prepared magnetization M_z');

%% Delta Q test
k=5;
len=length(tSL)/40;
reshp=reshape(tSL,40,len);
for time=1:20
    for k=1:5
        randCSL(k)=randsample(reshp(:,k),1);
    end
%     randCSL=randsample(tSL,k);
%     OrderCSL=sort(randCSL);
    YCSL = zeros(numel(k),1);
    for i=1:1:5
        j=randCSL(i)*1e3;
        YCSL(i)=Mz_CSL(round(j),1);
    end    
    fun3 = @(x,randCSL)((x(1)*exp(-randCSL/x(2))));
    [T1r1_fit] = lsqcurvefit(fun3,x0,randCSL,YCSL);
    M01(time)=T1r1_fit(1,1);
    lp1(time)=T1r1_fit(1,2);  %C-SL

    YBSL = zeros(numel(k),1);
    for i=1:1:5
        j=randCSL(i)*1e3;
        YBSL(i)=Mz_BSL(round(j),1);
    end    
    [T1r2_fit] = lsqcurvefit(fun3,x0,randCSL,YBSL);
    M02(time)=T1r2_fit(1,1);
    lp2(time)=T1r2_fit(1,2);  %B-SL

    YTRSL = zeros(numel(k),1);
    for i=1:1:5
        j=randCSL(i)*1e3;
        YTRSL(i)=Mz_TRSL(round(j),1);
    end    
    [T1r3_fit] = lsqcurvefit(fun3,x0,randCSL,YTRSL);
    M03(time)=T1r3_fit(1,1);
    lp3(time)=T1r3_fit(1,2);  %TR-SL
    
    YMRSL = zeros(numel(k),1);
    for i=1:1:5
        j=randCSL(i)*1e3;
        YMRSL(i)=Mz_MRSL(round(j),1);
    end    
    [T1r4_fit] = lsqcurvefit(fun3,x0,randCSL,YMRSL);
    M04(time)=T1r4_fit(1,1);
    lp4(time)=T1r4_fit(1,2);  %MR-SL
end
t1rT=0.1*ones(1,20);
deltaq1=((lp1./t1rT)-1)*100;
deltaq2=((lp2./t1rT)-1)*100;
deltaq3=((lp3./t1rT)-1)*100;
deltaq4=((lp4./t1rT)-1)*100;
mea1(heng,shu)=mean(abs(deltaq1));
mea2(heng,shu)=mean(abs(deltaq2));
mea3(heng,shu)=mean(abs(deltaq3));
mea4(heng,shu)=mean(abs(deltaq4));

t1rM=ones(1,20);
deltaM1=((M01./t1rM)-1)*100;
deltaM2=((M02./t1rM)-1)*100;
deltaM3=((M03./t1rM)-1)*100;
deltaM4=((M04./t1rM)-1)*100;
meaM1(heng,shu)=mean(abs(deltaM1));
meaM2(heng,shu)=mean(abs(deltaM2));
meaM3(heng,shu)=mean(abs(deltaM3));
meaM4(heng,shu)=mean(abs(deltaM4));
  end
heng = 0;
end

    figure(1)
    subplot(2,4,1);
    imshow((mea1),[0 10]);
    colormap('jet');
%     colorbar;
    axis on;
%     hold on
%     plot(x0,y0,'r.','markersize',30)

    xlabel('B_1 inhomogeneity [%]');
    ylabel('B_0 inhomogeneity [ppm]'); 

    subplot(2,4,2);  
    imshow((mea2),[0 10]);
    colormap('jet');
%     colorbar;
    axis on;
    xlabel('B_1 inhomogeneity [%]');
    ylabel('B_0 inhomogeneity [ppm]'); 
%     set(gca,'ColorScale','log')
    subplot(2,4,3);
    imshow((mea3),[0 25]);
    colormap('jet');
%     colorbar;
    axis on;
    xlabel('B_1 inhomogeneity [%]');
    ylabel('B_0 inhomogeneity [ppm]'); 
    subplot(2,4,4);
    imshow((mea4),[0 25]);
    colormap('jet');
    colorbar;
    axis on;
    xlabel('B_1 inhomogeneity [%]');
    ylabel('B_0 inhomogeneity [ppm]'); 
    subplot(2,4,5);
    imshow((meaM1),[0 25]);
    colormap('jet');
    axis on;
        xlabel('B_1 inhomogeneity [%]');
    ylabel('B_0 inhomogeneity [ppm]'); 
    subplot(2,4,6);  
    imshow((meaM2),[0 25]);
    colormap('jet');
    axis on;
        xlabel('B_1 inhomogeneity [%]');
    ylabel('B_0 inhomogeneity [ppm]'); 
    subplot(2,4,7);
    imshow((meaM3),[0 25]);
    colormap('jet');
    axis on;
        xlabel('B_1 inhomogeneity [%]');
    ylabel('B_0 inhomogeneity [ppm]'); 
    subplot(2,4,8);
    imshow((meaM4),[0 25]);
    colormap('jet');
    colorbar;
    axis on;
        xlabel('B_1 inhomogeneity [%]');
    ylabel('B_0 inhomogeneity [ppm]'); 

%% function for matrix propagators P

% Supporting Information Eq. 8
function P = Rx(phi)
    P = [ 1    0        0;
           0 cos(phi) -sin(phi);
           0 sin(phi) cos(phi) ];
end

% Supporting Information Eq. 9
function P = Ry(phi)
    P = [ cos(phi) 0 sin(phi);
               0    1     0;
          -sin(phi) 0 cos(phi)];
end

% Supporting Information Eq. 10
function P = Rsl(tau, T1rho, T2rho)
    P = [ exp(-tau/(T1rho))        0                    0;
                  0            exp(-tau/(T2rho))          0;
                  0                  0            exp(-tau/(T2rho))];
end

% Supporting Information Eq. 11;
function P = R90_1(dB1)
    P = Ry(pi/2 *dB1); % tip down pulse with dB1 deviation
end
% Supporting Information Eq. 11;
function P = R90_2(dB1)
    P = Ry(-pi/2 *dB1); % tip up pulse with dB1 deviation
end

% Supporting Information Eq. 12;
function P = R180_1(dB1)
    P = Rx(pi *dB1); % 1st refocusing pulse
end
% Supporting Information Eq. 12;
function P = R180_2(dB1)
    P = Rx(-pi *dB1); % 2nd refocusing pulse
end

% Supporting Information Eq. 13
function P = SL_1(tau, T1rho, T2rho, wSL, theta)
    P = Ry(-theta) * Rsl(tau, T1rho, T2rho) * Rx(wSL*tau) * Ry(theta); % 1st SL pulse
end
% Supporting Information Eq. 14
function P = SL_2(tau, T1rho, T2rho, wSL, theta)
    P = Ry(theta-pi) * Rsl(tau, T1rho, T2rho) * Rx(wSL*tau) * Ry(pi-theta); % 2nd SL pulse
end