close all
%% measurements from LHO - April 11, 2018 single bounce ITMY (?)
%[distance from edge of ISI [in], WX [um], WY [um]
% data = [10.25, 2619.92, 2133.15;
%     10 + 5/8, 2639.72, 2149.09;
%     12, 2568.43, 2083.96;
%     43.75, 1800.01, 1640.22;
%     47.125, 1788.59, 1452.71;
%     62.125 + 31.75, 994.95, 760.50; % dont use
%     62.125 + 34 + 3/8, 978.71, 485.59; % dont use
%     62.125 + 35.5, 936.45, 1117.25]; % dont use
% 

data = [-0.75, 2879.07, 2414.86;
    0, 2859.86, 2378.52;
    10.25, 2619.92, 2133.15;
    10 + 5/8, 2639.72, 2149.09;
    12, 2568.43, 2083.96;
%     43.75, 1800.01, 1640.22;
%     47.125, 1788.59, 1452.71;
    80.875, 1112.01, 1130;
    81.5, 1156.72, 1109.64;
    82.625, 1138.34, 1094.66;
    83.625, 1118.87, 1096.32;
    84.75, 1097.05, 1085.24;
    85.8125, 1124.53, 1102.31;
    87, 1061.07, 1086.64;
    88.125, 1045.31, 1099.34;
    100.562, 533.9*2, 573.5*2; % fitting backwards from PO2
    89.94 ,560.5*2, 548.79*2; % fitting backwards from PO2
    148.5, 1272.11, 1175.81;
    150.125, 1250.18, 1164.03;
    151.375, 1297.03, 1177.37;
    152.625, 1341.32, 1183.41;
    153.75, 1361.92, 1194.04;
    154.5, 1390.59, 1210.05;
    198.142, 2710.2, 2246.11;
    199.267, 2757.12, 2267.24;
    200.267, 2813.74, 2306.42;
    201.267, 2866.35, 2330.37;
    202.267, 2929.56, 2356.73;
    203.017, 2836.75, 2370.95;
    203.892, 2892.64, 2387.72;
    206.892, 2972.45, 2495.21];

edgeToOM1 = 1578E-3;
%SRMtoOM1 = 3.46;
SRMtoOM1 = 3.45;
SRMtoEdge = SRMtoOM1 - edgeToOM1;

data2 = data;
data2(:, 1) = data(:, 1)*25.4E-3 + SRMtoEdge; % distances in m from SRM
data2(:, 2:3) = 0.5*data(:, 2:3)*1E-6; % beam radius in um

figure(100)
plot(data2(:, 1), data2(:, 2), 'o');
hold all
plot(data2(:, 1), data2(:, 3), 'o');
grid on



% % beam measurements from jan 9 2018
% za  = [4.48 4.43 4.38 4.33 2.96 2.91 2.86 1.82 2.52];
% xwa = [952  926  912  891  1860 1920 1946 2881 2267]*1e-6/2;
% ywa = [1034 973  943  915  1745 1813 1871 2979 2252]*1e-6/2;
% % beam measurements from jan 15 2018
% zb  = [3.84 4.00 4.10 4.20 4.30 5.08 5.23 5.34 5.39 5.49 5.54 6.57 7.20 7.25];
% xwb = [1191 1016 964  920  935  1146 1167 1184 1328 1300 1329 1704 2523 2636]*1e-6/2;
% ywb = [949  875  877  880  876  1410 1394 1402 1423 1403 1411 1730 2196 2231]*1e-6/2;
% xdata = ([xwa xwb] + [ywa ywb])/2;
% zdata = [za zb];

% nominal guess of input beam
aNom = [-5, 2E-3]; % [-4, 2.1E-3];
aNom_design = [-3.88, 2.1E-3] % T1200410
% fit the input beam parameters
zdata = data2(:, 1);
xdata = (data2(:, 2) + data2(:, 3))/2;
xdata = data2(:, 2);
ydata = data2(:, 3);
    for kk = 2:2
    %xdata = (data2(:, kk));

    [betax,resnorm,resid,exitflag,output,lambda0,J] = lsqcurvefit(@SRMtoOMCBeamSizeRW, aNom, zdata, xdata);
    betax
    [betay,resnorm,resid,exitflag,output,lambda0,J] = lsqcurvefit(@SRMtoOMCBeamSizeRW, aNom, zdata, ydata);
    betay
    cix = nlparci(betax,resid,'jacobian',J)
    ciy = nlparci(betay,resid,'jacobian',J)
    figure(99)
    z = 0:0.01:8;
    %plot(z, SRMtoOMCBeamSizeRW_design(aNom,z),'k','LineWidth',3); hold all
    plot(z, SRMtoOMCBeamSizeRW_design(aNom_design,z),'k','LineWidth',3); hold all
    plot(z, SRMtoOMCBeamSizeRW(betax,z),'LineWidth',3); hold all
    plot(z, SRMtoOMCBeamSizeRW(betay,z),'LineWidth',3); hold all

    %plot(z, SRMtoOMCBeamSizeRW(ci(:,1),z));
    %plot(z, SRMtoOMCBeamSizeRW(ci(:,2),z));
    %plot(z, SRMtoOMCBeamSizeRW([ci(1,1), ci(2,2)],z));
    %plot(z, SRMtoOMCBeamSizeRW([ci(1,2), ci(2,1)],z));

    grid on
    hold all
    plot(zdata, xdata, 'bo')
    hold all
    plot(zdata, ydata, 'ro')
    xlabel('Distance from SRM [m]', 'FontSize',16)
    ylabel('Beam size [m]', 'FontSize',16)
    legend('design (T1200410)', 'x best fit','y best fit', 'x measurements','y measurements')
    % plot(zdata, [xwa xwb], 'o')
    % plot(zdata, [ywa ywb], 'o')
    end