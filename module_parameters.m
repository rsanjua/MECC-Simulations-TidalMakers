% module properties
% properties from:
% https://wateranywhere.com/2-5-x40-700-gpd-99-4-salt-rejection-watermaker-filmtec-seawater-membrane/
% module.H_ch = (28/1000)*(25.4/1000);% [m] channel height. It uses 28 mil spacer.
% module.epsil = 0.90; % [-] spacer porosity
% module.memArea = 2.8;% [m^2] membrane area.
% module.memLength = 40*25.4/1000;% [m] membrane lengtg. (40 inch)
% module.W = module.memArea/module.memLength;% [m] membrane width
% module.crossArea = 0.5*module.W*module.H_ch*module.epsil;%[m^2] cross section area of the channel in membrane
% module.Vol_mem = module.crossArea*module.memLength;% [m^3] channel volume
% module.dHydraulic=(2*module.epsil*module.H_ch)/(1+4*(1-module.epsil));% Hydraulic diameter, % unit=[m]

% properties for SW30-2514
module=struct;
module.memArea= 0.60387;% membrane area from manuals
module.memLength=(14-2*1.19)*0.3048;% membrane length from https://doi.org/10.1016/0011-9164(87)90107-X
module.H_ch=28*25.4/1000/1000;% channel height
module.epsil=0.9;% channel porosity
module.W=module.memArea/(module.memLength*2);
module.crossArea=module.W*module.H_ch*module.epsil;% cross section area
module.dHydraulic=4*module.epsil/(2/module.H_ch+(1-module.epsil)*8/module.H_ch);% hydraulic diameter
module.Vol_mem = module.crossArea*module.memLength;%module volume
module.L_p=1.4746/3600/1000;% m/s/bar. water permeability coefficient
module.B_s=0.13159/3600/1000;% m/s. salt permeability coefficient
module.P_max= 69;% maximum allowable pressure for module

% properties for SW30-2514
module1=struct;
module1.memArea= 0.60387;% membrane area from manuals
module1.memLength=(14-2*1.19)*0.3048;% membrane length from https://doi.org/10.1016/0011-9164(87)90107-X
module1.H_ch=28*25.4/1000/1000;% channel height
module1.epsil=0.9;% channel porosity
module1.W=module1.memArea/(module1.memLength*2);
module1.crossArea=module1.W*module1.H_ch*module1.epsil;% cross section area
module1.dHydraulic=4*module1.epsil/(2/module1.H_ch+(1-module1.epsil)*8/module1.H_ch);% hydraulic diameter
module1.Vol_mem = module1.crossArea*module1.memLength;%module volume
module1.L_p=1.4746/3600/1000;% m/s/bar. water permeability coefficient
module1.B_s=0.13159/3600/1000;% m/s. salt permeability coefficient
module1.P_max= 69;% maximum allowable pressure for module

% properties for SW30-2514
module2=struct;
module2.memArea= 0.60387;% membrane area from manuals
module2.memLength=(14-2*1.19)*0.3048;% membrane length from https://doi.org/10.1016/0011-9164(87)90107-X
module2.H_ch=28*25.4/1000/1000;% channel height
module2.epsil=0.9;% channel porosity
module2.W=module2.memArea/(module2.memLength*2);
module2.crossArea=module2.W*module2.H_ch*module2.epsil;% cross section area
module2.dHydraulic=4*module2.epsil/(2/module2.H_ch+(1-module2.epsil)*8/module2.H_ch);% hydraulic diameter
module2.Vol_mem = module2.crossArea*module2.memLength;%module volume
module2.L_p=1.4746/3600/1000;% m/s/bar. water permeability coefficient
module2.B_s=0.13159/3600/1000;% m/s. salt permeability coefficient
module2.P_max= 69;% maximum allowable pressure for module
