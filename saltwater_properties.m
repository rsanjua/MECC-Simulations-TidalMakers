% % salt water properties
% clc;clear;close all;
% % properties are for salt water:
% % the purpose is to obtain correlations based on salinity [g/kg]
% % the correlations for rho, mu, Pi are from: https://doi.org/10.1016/j.desal.2023.116407
% % the correlation for D is from https://doi.org/10.1016/j.memsci.2018.11.067
% 
% MW_s = 58.44; % [g/mol] molar mass of salt
% v = 2;% [-] dissociation coefficient of sodium chloride
% R = 8.314e-2;% [L bar/mol/K] universal gas constant
% T = 298.15;% [K] temperature
% m3toLiter = 0.001;% convert from m3 to L
% 
% %% polynomials with independent variable S [g/kg]
% 
% syms X S m C rho D phi mu Nu Pi;%
% % Description:
% % X: salt mass fraction [kg salt/kg total]
% % S: mass based salinity [g salt/kg total]
% % C: volume based salinity [kg/m^3]
% % m: molality [mol/kg total]
% % rho: Density [kg/m^3]
% % phi: osmotic coefficient [-]
% % mu: dynamic viscosity [kg/m/s]
% % Nu: kinematic viscosity [m^2/s]
% % Pi: osmotic pressure [bar]
% rho(X) = 765*X+995;
% rho_0 = eval(rho(0));% [kg/m^3] density of pure water.
% 
% D(X) = 153*X^4-122*X^3+30.1*X^2-2*X+1.51;
% phi(X) = 4.92*X^2+(8.89e-2)*X+0.918;
% mu(X) = (2.15e-3)*X+9.8e-4;
% Nu(X) = mu(X)/rho(X);
% Pi(X,m) = v*rho_0*R*T*phi(X)*m*m3toLiter;
% 
% m(S) = S/MW_s;
% X(S) = S/1000;
% Nu(S) = Nu(X(S));
% Pi(S) = Pi(X(S),m(S));
% rho(S) = rho(X(S));
% Sal(S) = rho(X(S))*X(S);
% D(S) = D(X(S));
% S_arr = linspace(0,260,500);% maximum salinity is 260 g/L at the saturation limit
% 
% rho_arr = zeros(size(S_arr));
% D_arr = zeros(size(S_arr));
% Pi_arr = zeros(size(S_arr));
% Nu_arr = zeros(size(S_arr));
% C_arr = zeros(size(S_arr));
% for ii = 1:length(S_arr)
%     rho_arr(ii) = eval(rho(S_arr(ii)));
%     C_arr(ii) = eval(Sal(S_arr(ii)));
%     D_arr(ii) = eval(D(S_arr(ii)));
%     Pi_arr(ii) = eval(Pi(S_arr(ii)));
%     Nu_arr(ii) = eval(Nu(S_arr(ii)));
% 
% end
% polyStoC=polyfit(S_arr,C_arr,2);
% polyRho_S = polyfit(S_arr,rho_arr,1);
% polyD_S = (1e-9)*polyfit(S_arr,D_arr,4);
% polyPi_S = polyfit(S_arr,Pi_arr,3);
% polyNu_S = polyfit(S_arr,Nu_arr,2);
% 
% polyCtoS=polyfit(C_arr,S_arr,5);
% polyRho_C = polyfit(C_arr,rho_arr,2);
% polyD_C = (1e-9)*polyfit(C_arr,D_arr,5);
% polyPi_C = polyfit(C_arr,Pi_arr,4);
% polyNu_C = polyfit(C_arr,Nu_arr,3);

% % results:
prop.polyRho_S = [0.765000000000001	995.000000000000];
prop.polyD_S=[1.53000000000003e-19	-1.22000000000001e-16	3.01000000000001e-14	-2.00000000000001e-12	1.51000000000000e-09];
prop.polyPi_S = [4.15291273858314e-06	7.50394191992022e-05	0.774872742686858	3.57206643910082e-15];
prop.polyNu_S = [-8.16665437845683e-13	1.37820254928011e-09	9.85451522449316e-07];
prop.polyRho_C  = [-0.000342243249588213	0.743467265313842	995.573913014884];
prop.polyD_C  = [-2.22116080836760e-22	2.84179455408789e-19	-1.40056985615032e-16	2.97363528206224e-14	-1.94820728160552e-12	1.50958714648752e-09];
prop.polyPi_C  = [-2.85490752893339e-09	3.40027146155346e-06	-0.000369607629368109	0.772652431155934	0.0568215674979449];
prop.polyNu_C  = [1.49997482652253e-15	-1.77816990194392e-12	1.38705250179577e-09	9.85248016699939e-07];
prop.polyStoC = [7.64999999999208e-07	0.995000000000000	-1.28508606274885e-14];
prop.polyCtoS = [9.19432125777173e-24	-2.32090173468475e-18	1.20015588572223e-12	-7.76590713528739e-07	1.00502512562814	-2.88810488163488e-14];
prop.polyDerRho_S = polyder(prop.polyRho_S);

% % % figure;box on; hold on;set(gca,'FontSize',12,'LineWidth',1);
% % % xlabel('Salinity [g/L]');ylabel('Osmotic pressure [bar]');
% % % C_arr= linspace(0,250,250);
% % % Pi_arr = polyval(prop.polyPi_C,C_arr);
% % % plot(C_arr,Pi_arr,'k-','LineWidth',2);
% % % % yline(70,'k--');text(200,70,'P=70 bar','VerticalAlignment','bottom');
% % % yline(41,'k--','LineWidth',1);text(200,41,'P=41 bar','VerticalAlignment','bottom','FontSize',12);
% % % yline(83,'k--','LineWidth',1);text(200,83,'P=83 bar','VerticalAlignment','bottom','FontSize',12);
% % % yline(120,'k--','LineWidth',2);text(200,120,'P=120 bar','VerticalAlignment','bottom','FontSize',12);
% % % ylim([0 220]);
% % % % yline(120)