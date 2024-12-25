%% ******* original version: Maximilian Gram, University of Wuerzburg *******
% 04.02.2024
% this code reproduces Fig. 2 of doi.org/10.1002/mrm.28585
% I added the effect of B1+ field deviations for the exciatioan and
% refocusing pulses. Both pulses are simulated as instantaneous operators.
% Relaxation effects during excitation and refocsuing are neglected. 
% T1rho dispersion effects from deviating B1+ fields are also neglected. 
%% ******* revised version: Cai Wan, University Medical Center Freiburg *******
% last revised 25.12.2024
% reference doi.org/10.1002/mrm.28585 and doi: 10.1002/nbm.4834 

clc;
clear all; 

df0   = 350; % [Hz] off-resonance 50mT: 2.323*10^6*130/10^6
dB1   = 0.75;                 % [ ] deviation of B1+ field, 1 -> optimal field
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

%% simulate: simple spin-lock
Mz_SSL = zeros(numel(tSL),1);

for j=1:numel(tSL)
    % Supporting Information Table S1 a)
    P = R90_2(dB1) * ...
        SL_1(tSL(j), T1rho, T2rho, wSL, theta) * ...
        R90_1(dB1);
    M = P * M0;
    Mz_SSL(j) = abs(M(end));
end

%% simulate: RE spin-lock
Mz_RESL = zeros(numel(tSL),1);

for j=1:numel(tSL)
    % Supporting Information Table S1 a)
    P = R90_2(dB1) * ...
        SL_2(tSL(j)/2, T1rho, T2rho, wSL, theta) * ...
        SL_1(tSL(j)/2, T1rho, T2rho, wSL, theta) * ...
        R90_1(dB1);
    M = P * M0;
    Mz_RESL(j) = abs(M(end));
end

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
        R180_2(dB1) * ...
        SL_2(tSL(j)/2, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R90_1(dB1);
    M = P * M0;
    Mz_BSL(j) = abs(M(end));
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
    P = R90_1(dB1) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_2(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R90_1(dB1);
    M = P * M0;
    Mz_TRSL(j) = abs(M(end));
end

%% simulate: Multi-refocused spin-lock
Mz_MRSL = zeros(numel(tSL),1);

for j=1:numel(tSL)

    %%
       P = R90_2(dB1) * ... 
        SL_1(tSL(j)/8, T1rho, T2rho, wSL, theta) * ...
        R180_2(dB1) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_2(dB1) * ...
        SL_1(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_2(tSL(j)/4, T1rho, T2rho, wSL, theta) * ...
        R180_1(dB1) * ...
        SL_1(tSL(j)/8, T1rho, T2rho, wSL, theta) * ...
        R90_1(dB1);   
    M = P * M0;
    Mz_MRSL(j) = abs(M(end));
end


%% plot results

% figure(1)
% plot(tSL, Mz_ref, 'k-')
% hold on
% plot(tSL, Mz_SSL, 'k--')
% legend('ref', 'SSL')
% figure(2)
% plot(tSL, Mz_ref, 'k-')
% hold on
% plot(tSL, Mz_RESL, 'r-')
% legend('ref', 'RESL')
figure(3)
plot(tSL, Mz_ref, 'k-')
hold on
plot(tSL, Mz_CSL, 'r-')
legend('ref', 'CSL')
figure(4)
plot(tSL, Mz_ref, 'k-')
hold on
plot(tSL, Mz_BSL, 'b--')
legend('ref', 'BSL')
figure(5)
plot(tSL, Mz_ref, 'k-')
hold on
plot(tSL, Mz_TRSL, 'k-')
legend('ref', 'TRSL')
figure(6)
plot(tSL, Mz_ref, 'k-')
hold on
plot(tSL, Mz_MRSL, 'r--')
legend('ref', 'MRSL')

%% monoexponential SL
% [pks,locs] = findpeaks(Mz_CSL);  % max peak
% y_max = pks;
% x_max = locs*1e-3;
% fun1 = @(x,x_max)((x(1)*exp(-x_max/x(2))));
% 
% nMz_BSL=1-Mz_CSL;                % min peak
% [pksL,locsL] = findpeaks(nMz_BSL);
% y_min = pksL;
% x_min = locsL*1e-3;
% fun2 = @(x,x_min)((x(1)*exp(-x_min/x(2))));
% 
% tSL = (1 : 1 : 200) *1e-3; % [s] SL preparation time
fun = @(x,tSL)((x(1)*exp(-tSL/x(2))));  % global fit 
x0 = [0.9,90];
% [T1r1_mxa] = lsqcurvefit(fun2,x0,x_max,y_max);
% [T1r1_min] = lsqcurvefit(fun1,x0,x_min,1-y_min);
% [T1r1] = lsqcurvefit(fun,x0,tSL,Mz_CSL');
% 
% figure(2)
% hold on
% box on;
% grid on;  
% set(gca, 'GridLineStyle', ':'); 
% plot(tSL,Mz_ref,'r-');
% fill([(tSL),fliplr(tSL)],[(fun2(T1r1_mxa,tSL))-0.003,fliplr(fun1(T1r1_min,tSL))],[255/255 220/255 255/255],'EdgeColor','none');
% plot(tSL,Mz_CSL,'k-',tSL,fun(T1r1,tSL),'b--')
% legend('monoexponential reference','area of min/max T_1_\rho quantification','trajectory M_z(T_S_L)','global fit');
% xlabel('Spin lock time T_S_L [s]');
% ylabel('Prepared magnetization M_z');

%% Delta Q test
k=8;
tSL1 = (1 : 1 : 25) *1e-3; % [s] SL preparation time
for time=1:100
    for i=1:1:k
     randCSL(i)=randsample(tSL1+(i-1)*(25*1e-3),1);
    end
    YCSL = zeros(numel(k),1);
    for i=1:1:k
        j=randCSL(i)*1e3;
        YCSL(i)=Mz_CSL(round(j),1);
        YBSL(i)=Mz_BSL(round(j),1);
        YTRSL(i)=Mz_TRSL(round(j),1);
        YMRSL(i)=Mz_MRSL(round(j),1);
    end
    [T1r1_fit] = lsqcurvefit(fun,x0,randCSL,YCSL);
    M01(time)=T1r1_fit(1,1);
    lp1(time)=T1r1_fit(1,2);  %C-SL
    
    [T1r2_fit] = lsqcurvefit(fun,x0,randCSL,YBSL);
    M02(time)=T1r2_fit(1,1);
    lp2(time)=T1r2_fit(1,2);  %B-SL
    
    [T1r3_fit] = lsqcurvefit(fun,x0,randCSL,YTRSL);
    M03(time)=T1r3_fit(1,1);
    lp3(time)=T1r3_fit(1,2);  %TR-SL
    
    [T1r4_fit] = lsqcurvefit(fun,x0,randCSL,YMRSL);
    M04(time)=T1r4_fit(1,1);
    lp4(time)=T1r4_fit(1,2);  %QR-SL
end
t1rT=0.1*ones(1,100);
deltaq1=((lp1./t1rT)-1)*100;
mea1=mean(mean(abs(deltaq1)));

deltaq2=((lp2./t1rT)-1)*100;
mea2=mean(mean(abs(deltaq2)));

deltaq3=((lp3./t1rT)-1)*100;
mea3=mean(mean(abs(deltaq3)));

deltaq4=((lp4./t1rT)-1)*100;
mea4=mean(mean(abs(deltaq4)));

mx=max(deltaq1);
mn=min(deltaq1);
a00=sum( deltaq1 >= 30);
a01=sum(deltaq1 < 30 & deltaq1 >= 20);
a10=sum(deltaq1 < 20 & deltaq1 >= 10);
a11=sum(deltaq1 < 10 & deltaq1 >= 0);
a12=sum(deltaq1 < 0 & deltaq1 >= -10);
a13=sum(deltaq1 < -10 & deltaq1 >= -20);
a14=sum(deltaq1 < -20 & deltaq1 >= -30);
a15=sum(deltaq1 < -30);

t1rM=ones(1,100);
deltaM1=((M01./t1rM)-1)*100;
meaM1=mean(mean(abs(deltaM1)));

deltaM2=((M02./t1rM)-1)*100;
meaM2=mean(mean(abs(deltaM2)));

deltaM3=((M03./t1rM)-1)*100;
meaM3=mean(mean(abs(deltaM3)));

deltaM4=((M04./t1rM)-1)*100;
meaM4=mean(mean(abs(deltaM4)));

mxM=max(deltaM1);
mnM=min(deltaM1);
b00=sum( deltaM1 >= 10);
b000=sum(deltaM1 < 10 & deltaM1 >= 0);
b01=sum(deltaM1 < 0 & deltaM1 >= -10);
b10=sum(deltaM1 < -10 & deltaM1 >= -20);
b11=sum(deltaM1 < -20 & deltaM1 >= -30);
b12=sum(deltaM1 < -30 & deltaM1 >= -40);
b13=sum(deltaM1 < -40 & deltaM1 >= -50);
b14=sum(deltaM1 < -50 );

x1 = [35,25,15,5,-5,-15,-25,-35];
y1 = [0,0,2,59,39,0,0,0];
figure
bar(x1,y1,1);
hold on 
legend('errors of C-SL');
title('The deviation of fitting - T_1_\rho');
xlabel('T_1_\rho quantification errors [%]');
ylabel('Relative frequency [%]');

x = [15,5,-5,-15,-25,-35,-45,-55];
y = [0,0,32,68,0,0,0,0];
figure
bar(x,y,1);
hold on 
legend('errors of C-SL');
title('The deviation of fitting - M_z');
xlabel('T_1_\rho quantification errors [%]');
ylabel('Relative frequency [%]');

%% function for matrix propagators P

% Supporting Information Eq. 8
function P = Rx(phi)
    P = [ 1    0        0;
           0 cos(phi) sin(phi);
           0 -sin(phi) cos(phi) ];
end

% Supporting Information Eq. 9
function P = Ry(phi)
    P = [ cos(phi) 0 -sin(phi);
               0    1     0;
          sin(phi) 0 cos(phi)];
end

% Supporting Information Eq. 10
function P = Rsl(tau, T1rho, T2rho)
    P = [ exp(-tau/(T2rho))        0                    0;
                  0            exp(-tau/(T1rho))          0;
                  0                  0            exp(-tau/(T2rho))];
end

% Supporting Information Eq. 11;
function P = R90_1(dB1)
    P = Rx(pi/2 *dB1); % tip down pulse with dB1 deviation
end
% Supporting Information Eq. 11;
function P = R90_2(dB1)
    P = Rx(-pi/2 *dB1); % tip up pulse with dB1 deviation
end

% Supporting Information Eq. 12;
function P = R180_1(dB1)
    P = Ry(pi *dB1); % 1st refocusing pulse
end
% Supporting Information Eq. 12;
function P = R180_2(dB1)
    P = Ry(-pi *dB1); % 2nd refocusing pulse
end

% Supporting Information Eq. 13
function P = SL_1(tau, T1rho, T2rho, wSL, theta)
    P = Rx(-theta) * Rsl(tau, T1rho, T2rho) * Ry(wSL*tau) * Rx(theta); % 1st SL pulse
end
% Supporting Information Eq. 14
function P = SL_2(tau, T1rho, T2rho, wSL, theta)
    P = Rx(theta-pi) * Rsl(tau, T1rho, T2rho) * Ry(wSL*tau) * Rx(pi-theta); % 2nd SL pulse
end