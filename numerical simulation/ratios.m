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
y = logspace(1.3,2.7,25);
for ff1 = 1:1:25
df0   = 420; % [Hz] off-resonance 50mT: 2.323*10^6*130/10^6  360 425
dB1   = 0.85;                 % [ ] deviation of B1+ field, 1 -> optimal field 0.85
T1rho = 100*1e-3;	         % [s] T1rho relaxation time
T2rho = y(1,ff1)*1e-3;           % [s] T2rho relaxation time
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
% figure()
% hold on
% plot(tSL, Mz_ref, 'k--')
% plot(tSL, Mz_CSL, 'r-')
% plot(tSL, Mz_BSL, 'b-')
% plot(tSL, Mz_TRSL, 'k-')
% legend('ref', 'CSL', 'BSL','TRCSL')


%% monoexponential SL
[pks,locs] = findpeaks(Mz_TRSL);  % max peak
y_max = pks;
x_max = locs*1e-3;
fun1 = @(x,x_max)((x(1)*exp(-x_max/x(2))));

nMz_BSL=1-Mz_TRSL;                % min peak
[pksL,locsL] = findpeaks(nMz_BSL);
y_min = pksL;
x_min = locsL*1e-3;
fun2 = @(x,x_min)((x(1)*exp(-x_min/x(2))));

tSL = (1 : 1 : 200) *1e-3; % [s] SL preparation time
fun = @(x,tSL)((x(1)*exp(-tSL/x(2))));  % global fit 
x0 = [1,90];
[T1r1_mxa] = lsqcurvefit(fun2,x0,x_max,y_max);
[T1r1_min] = lsqcurvefit(fun1,x0,x_min,1-y_min);
[T1r1] = lsqcurvefit(fun,x0,tSL,Mz_TRSL');


%% Delta Q test
k=8;
tSL1 = (1 : 1 : 25) *1e-3; % [s] SL preparation time
for time=1:100
    for i=1:1:k
     randCSL(i)=randsample(tSL1+(i-1)*(25*1e-3),1);
    end
    OrderCSL=sort(randCSL);
    YCSL = zeros(numel(k),1);
    for i=1:1:8
        j=OrderCSL(i)*1e3;
        YCSL(i)=Mz_CSL(round(j),1);
    end    
    fun3 = @(x,OrderCSL)((x(1)*exp(-OrderCSL/x(2))));
    [T1r1_fit] = lsqcurvefit(fun3,x0,OrderCSL,YCSL);
    M01(time)=T1r1_fit(1,1);
    lp1(time)=T1r1_fit(1,2);  %C-SL

    YBSL = zeros(numel(k),1);
    for i=1:1:8
        j=OrderCSL(i)*1e3;
        YBSL(i)=Mz_BSL(round(j),1);
    end    
    [T1r2_fit] = lsqcurvefit(fun3,x0,OrderCSL,YBSL);
    M02(time)=T1r2_fit(1,1);
    lp2(time)=T1r2_fit(1,2);  %B-SL

    YTRSL = zeros(numel(k),1);
    for i=1:1:8
        j=OrderCSL(i)*1e3;
        YTRSL(i)=Mz_TRSL(round(j),1);
    end    
    [T1r3_fit] = lsqcurvefit(fun3,x0,OrderCSL,YTRSL);
    M03(time)=T1r3_fit(1,1);
    lp3(time)=T1r3_fit(1,2);  %TR-SL
    
    YMRSL = zeros(numel(k),1);
    for i=1:1:8
        j=OrderCSL(i)*1e3;
        YMRSL(i)=Mz_MRSL(round(j),1);
    end    
    [T1r4_fit] = lsqcurvefit(fun3,x0,OrderCSL,YMRSL);
    M04(time)=T1r4_fit(1,1);
    lp4(time)=T1r4_fit(1,2);  %MR-SL
end
t1rT=0.1*ones(1,100);
deltaq1=((lp1./t1rT)-1)*100;
deltaq2=((lp2./t1rT)-1)*100;
deltaq3=((lp3./t1rT)-1)*100;
deltaq4=((lp4./t1rT)-1)*100;
a10(ff1)=sum(abs(deltaq1) < 1 & abs(deltaq1) >= 0);
a11(ff1)=sum(abs(deltaq2) < 1 & abs(deltaq2) >= 0);
a12(ff1)=sum(abs(deltaq3) < 1 & abs(deltaq3) >= 0);
a13(ff1)=sum(abs(deltaq4) < 1 & abs(deltaq4) >= 0);

t1rM=ones(1,100);
deltaM1=((M01./t1rM)-1)*100;
deltaM2=((M02./t1rM)-1)*100;
deltaM3=((M03./t1rM)-1)*100;
deltaM4=((M04./t1rM)-1)*100;
b10(ff1)=sum(abs(deltaM1) < 1 & abs(deltaM1) >= 0);
b11(ff1)=sum(abs(deltaM2) < 1 & abs(deltaM2) >= 0);
b12(ff1)=sum(abs(deltaM3) < 1 & abs(deltaM3) >= 0);
b13(ff1)=sum(abs(deltaM4) < 1 & abs(deltaM4) >= 0);

end

sl=logspace(1.3,2.7,25);
figure(1)
hold on
plot(sl,a10,'k.-');
hold on
plot(sl,a11,'r.-');
hold on
plot(sl,a12,'b.-');
hold on
plot(sl,a13,'g.-');

legend('C-SL','B-SL','TR-SL','MR-SL');
xlabel('Relaxation times ratio T_2_\rho: T_1_\rho');
ylabel('Proportions meeting criterion of T_1_\rho \DeltaQ < 1%');  
title('Spin lock amplitude {\itf}_S_L =500Hz');
figure(2)
 hold on
plot(sl,b10,'k.-');
hold on
plot(sl,b11,'r.-');
hold on
plot(sl,b12,'b.-');
hold on
plot(sl,b13,'g.-');
legend('C-SL','B-SL','TR-SL','MR-SL');
xlabel('Relaxation times ratio T_2_\rho: T_1_\rho');
ylabel('Proportions meeting criterion of T_1_\rho-\DeltaQ < 0.5%');   
title('Spin lock amplitude {\itf}_S_L =500Hz');
box on;
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