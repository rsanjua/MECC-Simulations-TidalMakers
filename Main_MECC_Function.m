%% Main iteration function
% This is the main function that handles the constant configuration, time step iteration, and algorithm logic for the system.


%% default algorithm function header
% function [y, ydot, ydotdot, yp2, yp2dot, yp2dotdot, Fb, F1, F2, Pp1, Pm1, Pp2, Pm2, ...
%      Cf1, Cf2, theta, thetadot, thetadotdot, J1, J2, J1s, J2s, Hres] ...
%     = uncap_upward_pp(t, wh, ... known variables
%     yn1, ydotdotn1, yp2n1, Fbn1, F1n1, F2n1, Pp1n1, Pm1n1, Pp2n1, Pm2n1, ...
%     Cf1n1, Cf2n1, thetan1, J1n1, J2n1, J1sn1, J2sn1, Hresn1, ... Variables at previous time step
%     ...
%     yn2, yn3, yp2n2, yp2n3, thetan2, thetan3, ... Extra old info past original time step
%     ...
%     delt, mb, mp2, mp1, mlev, rho, mu, g, R, T, phib, req, r, Dpipe, Apipe,...
%     Ap1, A1R, Ap2, A2R, Am1, Am2, Ares, Vfeed1, Vfeed2, Wdepth, Lo, Lo2, ...
%     Lc1, Lc2, H1, H2, L1, L2, Ls1, Ls2, Ls3, Ls4, z1, z2, z2mem1, z1mem1, iv, ...
%     z1mem2, yeqp1, yeqp2, I, w1, w2, Pw1, Pw2, Ps1, Ps2, Cseawater, Pb1, Pres, ...
%     Pb2, ff1, ff2, Lm1, Lm2, Dm1, Dm2, alpha, Pep1, Pep2, Ce1, Ce2) % constants

% For now, will treat the system as a prototype being tested in a wave tank

%% Wave parameters
wave_period = 8; % Wave period [s]
wave_height = 0.1905; % Significant wave height [m] (measured from peak to trough) {guess based on max piston chamber} (0.1905)
water_depth = 1; % Ocean depth [m] (1)
initial_salinity = 35; % Ocean water salinity [g/L]
temperature = 17; % Ocean temperature [°C]
rho = 1025; % density of seawater [kg/m^3]
mu = 0.001159; % Viscosity of seawater [Pa*m]
molmass = 58.44; % Molar mass of salt [g/mol]
saltwater_properties;


%% Constants
g = 9.8;
R = 0.0831446261815324; % ideal gas constant [barL/Kmol]
vant_hoff_constant = 2; % vant hoff constant for salt water
R = R/molmass; % convert to [barL/K*g_salt]

%% Scale factor
piston_scale = 1;
lever_scale = 1;
number_of_membranes = 1; % number of membranes in parallel per pass (only implemented in ryan's code)

%% Customizable Parameters (These are all initial guesses that we change to figure things out)
piston1_bore_diameter = 0.1016 * piston_scale; % bore diameter of piston 1 [m] {manufacturing drawings}
piston1_rod_diameter = 0.0508 * piston_scale; % rod diameter of piston 1 [m] {manufacturing drawings}
piston1_area = piston1_bore_diameter^2 *pi/4; % Area of the piston face under the buoy [m^2]
piston1_mass = 2 * piston_scale^3; % Mass of piston 1 [kg] {guess}
piston1_rodside_area = piston1_area - piston1_rod_diameter^2*pi/4; % Area of the rod side of the piston [m^2]
piston1_equilibrium_elevation = 0.2413 * piston_scale; % height of middle of piston 1 range of motion from the seafloor [m] {diya notebook +}
    % equilibrium elevation = yeq + piston1 elevation + piston1 unused
    % length
piston1_elevation = 0.0508 * piston_scale; % height of piston 1 chamber bottom from seafloor [m]
piston1_thickness = 0.1016 * piston_scale; % Thickness of piston 1 disk

    
% PLEASE LETS ASSUME BOTH PISTONS ARE STARTING AT THE SAME HEIGHT PLEASE

piston2_bore_diameter = 0.0508 * piston_scale; % bore diameter of piston 2 [m] {diya notebook}
piston2_rod_diameter = 0.03302 * piston_scale; % rod diameter of piston 2 [m] {diya notebook}
piston2_area = piston2_bore_diameter^2 *pi/4; % Area of the piston face under the lever arm [m^2]
piston2_mass = 1 * piston_scale^3; % Mass of piston 2 [kg]
piston2_rodside_area = piston1_area - piston2_rod_diameter^2*pi/4; % Area of the rod side of the piston [m^2]
piston2_equilibrium_elevation = piston1_equilibrium_elevation; % please
piston2_elevation = 0.0508 * piston_scale; % height of piston 2 chamber bottom from seafloor [m]
piston2_thickness = 0.0508 * piston_scale; % Thickness of piston 2 disk

buoy_mass = 100; % Mass of buoy and piston 1 [kg] (100)
buoy_radius_small = 0.1524; % Radius of bottom of buoy [m] (0.1524)
buoy_cone_angle = 53.1; % Angle of taper [°]

lever_length = 0.8636 * lever_scale; % Length of lever arm [m]
pistdist = lever_length; % distance between the center of the two pistons [m]
lever_mass = 1 * lever_scale^3; % Mass of lever [kg]
lever_width = 0.1524 * lever_scale; % Width of lever plate [m]
lever_thickness = 0.003175 * lever_scale; % Thickness of one plate in the lever [m]
pivot_point_height = 0.635*piston_scale; % Height of pivot point from piston chamber midpoint [m]
pivot_piston1_horizontal_distance = 0.381 * lever_scale; % horizontal distance between center of piston 1 and lever pivot point [m]
pivot_piston2_horizontal_distance = pistdist - pivot_piston1_horizontal_distance; % horizontal distance between center of piston 2 and lever pivot point [m]
lever_length_piston1_side = 0.381 * lever_scale; % Length of lever from pivot to tip (facing piston 1) [m]
lever_length_piston2_side = lever_length - lever_length_piston1_side; % Length of lever from pivot to tip (facing piston 2) [m]
connecting_rod_1_length = pivot_point_height-piston1_equilibrium_elevation; %0.254; % Length of rod connecting piston 1 and lever arm [m]
connecting_rod_2_length = pivot_point_height-piston2_equilibrium_elevation; %0.254; % Length of rod connecting piston 2 and lever arm [m] 

pipe_diameter = 0.0508; % Pipe diameter [m] (prototype has 0.0508)
pipe_length_segment_1 = 0.1524 ; % Characteristic pipe length between inlet and piston 2 chamber (segment 1) [m]
pipe_length_segment_2 = 0.1524 * lever_scale; % Characteristic pipe length between piston 2 chamber and membrane 1 inlet (segment 2) [m]
pipe_length_segment_3 = 0.2032 * lever_scale; % Characteristic pipe length between membrane 1 outlet and piston 1 chamber (segment 3) [m]
pipe_length_segment_4 = 0.2032; % Characteristic pipe length between piston 1 chamber and membrane 2 inlet (segment 4) [m]
segment1_loss = 5.3; % Minor head loss coefficient for segment 1
segment2_loss = 5; % Minor head loss coefficient for segment 2 (assume laminar flow)
segment3_loss = 4.5; % Minor head loss coefficient for segment 3
segment4_loss = 5; % Minor head loss coefficient for segment 4
segment1_elevation = 0.0508*piston_scale; % elevation difference across segment 1 [m]
segment2_elevation = 0.0508*piston_scale; % elevation difference across segment 2 [m]
segment3_elevation = 0.0508 * piston_scale; % elevation difference across segment 3 [m]
segment4_elevation = 0.0508 * piston_scale; % elevation difference across segment 4 [m]

reservoir_area = 2 * piston_scale; % Area of base of reservoir between membrane 1 and piston 1 [m^2]
reservoir_start_height = 4 * piston_scale; % inital height of reservoir [m]
brine_recycle_percent = 0.5; % percent of second pass brine recycled back to the reservoir [%]

relaxation_constant = 0.5; % SOR omega value
residual_cutoff = 0.2; % residual cutoff value

membrane_cell_count = 100; % number of cells in membrane CFD

%% Membrane properties
module_parameters;
mem1_feed_qmax = 0.0003333*number_of_membranes; % Max feed volumetric flow rate [m^3/s]
mem1_permeate_qmax = 0.6/24/60/60*number_of_membranes; % Max permeate volumetric flow rate [m^3/s]
mem1_length = 0.3556; % Length of membrane module [m]
mem1_diameter = 0.0635; % Diameter of membrane module [m]
mem1_volume = module1.Vol_mem*number_of_membranes; % Volume on feed side [m^3]
mem1_area = module1.memArea*number_of_membranes; % Surface area of membrane [m^2]
mem1_perm_salt = module1.B_s; % Salt permeability of membrane 1 [m/s] STOLEN FROM SW30 25-40 MEMBRANE (UNITS ARE WEIRD)
mem1_perm_water = module1.L_p; % Water permeability of membrane 2 [m/s/bar] STOLEN FROM SW30 25-40 MEMBRANE
mem1_friction_factor = 0.1; % Membrane 1 friction factor
mem1_exit_pressure = 1; % Initial membrane 1 permeate exit pressure [bar]


% its 2 am and i need results im pretending theyre the same
mem2_feed_qmax = 0.0003333*number_of_membranes; % Max feed volumetric flow rate [m^3/s]
mem2_permeate_qmax = 0.6/24/60/60*number_of_membranes; % Max permeate volumetric flow rate [m^3/s]
mem2_length = 0.3556; % Length of membrane module [m]
mem2_diameter = 0.0635; % Diameter of membrane module [m]
mem2_volume = module2.Vol_mem*number_of_membranes; % Volume on feed side [m^3]
mem2_area = module2.memArea*number_of_membranes; % Surface area of membrane [m^2]
mem2_perm_salt = module2.B_s; % Salt permeability of membrane [m/s] STOLEN FROM SW30 25-40 MEMBRANE (UNITS ARE WEIRD)
mem2_perm_water = module2.L_p; % Water permeability of membrane [m/s/bar] STOLEN FROM SW30 25-40 MEMBRANE
mem2_friction_factor = 0.1; % Membrane 2 friction factor
mem2_exit_pressure = 1; % Initial membrane 2 permeate exit pressure [bar]


%% Calculatable values
buoy_eq_height = (((buoy_mass)/rho + (pi*buoy_radius_small^3)/(3*tand(buoy_cone_angle))) *3/(pi*tand(buoy_cone_angle)*tand(buoy_cone_angle)))^(1/3) - buoy_radius_small/tand(buoy_cone_angle)+0.462; % % Height from bottom of buoy to water level at equilibrium position [m]
buoy_eq_radius = buoy_radius_small + buoy_eq_height*tand(buoy_cone_angle); % Radius of buoy cross section at equilibrium height [m]
hydrostatic_pressure = water_depth * g * rho*10^(-5); % Hydrostatic pressure at inlet [bar]
piston1_amplitude = wave_height;
piston1_unused_length = piston1_equilibrium_elevation-piston1_amplitude-piston1_elevation; % height of unused portion at the bottom of piston 1 [m]
maxtheta = thetacalcy(piston1_amplitude,lever_length_piston1_side,connecting_rod_1_length,pivot_piston1_horizontal_distance); 
piston2_amplitude = -lever_length_piston2_side*sind(maxtheta) + connecting_rod_2_length*(1-sqrt(1-((pivot_piston2_horizontal_distance-lever_length_piston2_side*cosd(maxtheta))/connecting_rod_2_length)^2));
piston2_unused_length = piston2_equilibrium_elevation-piston2_amplitude-piston2_elevation; % height of unused portion at the bottom of piston 1 [m]
reservoir_pressure = hydrostatic_pressure; % Pressure at the top of the reservoir [bar]

piston1_equilibrium_height = piston1_amplitude; % Height of piston 1 in chamber when lever arm is horizontal (from top of unused portion) [m]
piston2_equilibrium_height = piston2_amplitude; % Height of piston 2 in chamber when lever arm is horizontal (from top of unused portion) [m]
mem1_exit_salt_concentration = 2.5; % Salt concentration at permeate exit in membrane 1 [g/L] {guesses}
mem2_exit_salt_concentration = 0.5; % Salt concentration at permeate exit in membrane 2 [g/L] {guesses}
mem1_brine_pressure = 1; % pressure of brine at exit of membrane 1 [bar]
mem2_brine_pressure = 1; % pressure of brine at exit of membrane 2 [bar]
segment1_alpha = 2; % alpha constant in segment 1 --> laminar flow = 2, turbulent = 1
segment2_alpha = 2; % alpha constant in segment 2
segment3_alpha = 2; % alpha constant in segment 3
segment4_alpha = 2; % alpha constant in segment 4
pipe_area = pipe_diameter^2/4*pi; % Cross sectional area of pipe
mem1_inlet_area = pipe_area*number_of_membranes; % Cross sectional area of membrane inlet hole
mem2_inlet_area = pipe_area*number_of_membranes; % Cross sectional area of membrane inlet hole
moment_of_inertia = 1/12*lever_mass*(lever_length^2+lever_width^2) + (lever_width*lever_length)*((lever_length_piston1_side+lever_length_piston2_side)/2 - lever_length_piston1_side);
piston1_velocity_max_pp = mem2_permeate_qmax/piston1_area; % restriction on the velocity  of piston 1 during permeate production
piston1_velocity_max_fl = mem2_feed_qmax/piston1_area; % restriction on the velocity of piston 1 during flushing
piston2_velocity_max_pp = mem1_permeate_qmax/piston2_area; % restriction on the velocity of piston 2 during permeate production
piston2_velocity_max_fl = mem1_feed_qmax/piston2_area; % restriction on the velocity of piston 2 during flushing
piston1_negative_velocity_max_pp = piston2_velocity_max_pp * lever_length_piston2_side/lever_length_piston1_side; % restriction on the velocity of piston 1 based on the velocity of piston 2 during permeate produciton
piston1_negative_velocity_max_fl = piston2_velocity_max_fl * lever_length_piston2_side/lever_length_piston1_side; % restriction on the velocity of piston 1 based on the velocity of piston 2 during flushing
piston1_disc_eq_height_seafloor = piston1_equilibrium_height + piston1_unused_length + piston1_elevation;
piston2_disc_eq_height_seafloor = piston2_equilibrium_height + piston2_unused_length + piston2_elevation;
mem1_delta_x = mem1_length/membrane_cell_count;
mem2_delta_x = mem2_length/membrane_cell_count;




%% Property structure
c=struct; %all constants

c.mb = buoy_mass;
c.mp2 = piston2_mass;
c.mp1 = piston1_mass;
c.mlev = lever_mass;
c.rho = rho;
c.mu = mu;
c.g = g;
c.R = R;
c.T = temperature;
c.phib = buoy_cone_angle;
c.yeq = buoy_eq_height;
c.req = buoy_eq_radius;
c.r = buoy_radius_small;
c.Dpipe = pipe_diameter;
c.Apipe = pipe_area;
c.Ap1 = piston1_area;
c.A1R = piston1_rodside_area;
c.Ap2 = piston2_area;
c.A2R = piston2_rodside_area;
c.Am1 = mem1_area;
c.Am2 = mem2_area;
c.Am1in = mem1_inlet_area;
c.Am2in = mem2_inlet_area;
c.Ares = reservoir_area;
c.Vfeed1 = mem1_volume;
c.Vfeed2 = mem2_volume;
c.Wdepth = water_depth;
c.Lo = piston1_unused_length;
c.Lo2 = piston2_unused_length;
c.Lc1 = connecting_rod_1_length;
c.Lc2 = connecting_rod_2_length;
c.H1 = connecting_rod_1_length;
c.H2 = connecting_rod_2_length;
c.L1 = pivot_piston1_horizontal_distance;
c.L2 = pivot_piston2_horizontal_distance;
c.Ls1 = pipe_length_segment_1;
c.Ls2 = pipe_length_segment_2;
c.Ls3 = pipe_length_segment_3;
c.Ls4 = pipe_length_segment_4;
c.tp1 = piston1_thickness;
c.tp2 = piston2_thickness;
c.z1 = piston1_elevation;
c.z2 = piston2_elevation;
c.z2mem1 = segment2_elevation;
c.z1mem1 = segment3_elevation;
c.iv = vant_hoff_constant;
c.z1mem2 = segment4_elevation;
c.yeqp1 = piston1_equilibrium_height;
c.yeqp2 = piston2_equilibrium_height;
c.I = moment_of_inertia;
c.w1 = lever_length_piston1_side;
c.w2 = lever_length_piston2_side;
c.Pw1 = mem1_perm_water;
c.Pw2 = mem2_perm_water;
c.Ps1 = mem1_perm_salt;
c.Ps2 = mem2_perm_salt;
c.Pinlet = hydrostatic_pressure;
c.Cseawater = initial_salinity;
c.Pb1 = mem1_brine_pressure;
c.Pres = reservoir_pressure;
c.Pb2 = mem2_brine_pressure;
c.ff1 = mem1_friction_factor;
c.ff2 = mem2_friction_factor;
c.Lm1 = mem1_length;
c.Lm2 = mem2_length;
c.Dm1 = mem1_diameter;
c.Dm2 = mem2_diameter;
c.alpha1 = segment1_alpha;
c.alpha2 = segment2_alpha;
c.alpha3 = segment3_alpha;
c.alpha4 = segment4_alpha;
c.omega = relaxation_constant;
c.residual_cutoff = residual_cutoff;
c.qmaxf1 = mem1_feed_qmax;
c.qmaxp1 = mem1_permeate_qmax;
c.qmaxf2 = mem2_feed_qmax;
c.qmaxp2 = mem2_permeate_qmax;
c.vmaxp1 = piston1_velocity_max_pp;
c.vmaxf1 = piston1_velocity_max_fl;
c.vmaxp2 = piston2_velocity_max_pp;
c.vmaxf2 = piston2_velocity_max_fl;
c.vmaxnegp1 = piston1_negative_velocity_max_pp;
c.vmaxnegf1 = piston1_negative_velocity_max_fl;
c.heqp1 = piston1_disc_eq_height_seafloor;
c.heqp2 = piston2_disc_eq_height_seafloor;
c.Yrec = brine_recycle_percent;
c.N = membrane_cell_count;
c.delx1 = mem1_delta_x;
c.delx2 = mem2_delta_x;



%% Time step discretization
numsteps = 4000; % Number of time grid points
total_time = 64; % Total simulation time [s]
dt = total_time / (numsteps - 1); % Time steps
c.delt = dt;


%% Arrays
t = zeros(1, numsteps); % Time grid [s]
wh = zeros(1, numsteps); % Wave height grid [m] (0 = sea level)
y = zeros(1, numsteps); % Buoy height grid [m] (0 = sea level)
ydot = zeros(1, numsteps); % Buoy velocity grid [m/s]
ydotdot = zeros(1, numsteps); % Buoy acceleration grid [m/s^2]
yp2 = zeros(1, numsteps); % Piston 2 height [m]
yp2dot = zeros(1, numsteps); % Piston 2 velocity [m/s]
yp2dotdot = zeros(1, numsteps); % Piston 2 acceleration
Fb = zeros(1, numsteps); % Upward force of buoy [N] THIS IS THE MOORING FORCE
F1 = zeros(1, numsteps); % Tension force in rod connecting piston 1 and lever arm [N]
F2 = zeros(1, numsteps); % Tension force in rod connecting piston 2 and lever arm [N]
Pp1 = zeros(1, numsteps); % Pressure at face of piston1 [bar]
Pm1 = zeros(1, numsteps); % Pressure at the entrance of membrane 1 [bar]
Pep1 = zeros(1, numsteps); % Pressure of permeate at membrane 1 exit [bar]
Pp2 = zeros(1, numsteps); % Pressure at face of piston 2 [bar]
Pm2 = zeros(1, numsteps); % Pressure at the entrance of membrane 2 [bar]
Pep2 = zeros(1, numsteps); % Pressure of permeate at membrane 2 exit [bar]
Cf1 = zeros(1, numsteps); % Salt concentration on feed side at membrane 1 [g/L]
Ce1 = zeros(1, numsteps); % Salt concentration at exit of membrane 1 [g/L]
Cf2 = zeros(1, numsteps); % Salt concentration on feed side at membrane 2 [g/L]
Ce2 = zeros(1, numsteps); % Salt concentration at exit of membrane 2 [g/L]
Cres = zeros(1, numsteps); % Salt concentration in reservoir [g/L]
theta = zeros(1, numsteps); % Angle of lever arm (0 is horizontal, clockwise rotation is positive) [degrees]
thetadot = zeros(1, numsteps); % Angular velocity of lever arm [degrees/s]
thetadotdot = zeros(1, numsteps); % Angular acceleration of lever arm [degrees/s2]
J1 = zeros(1, numsteps); % Permeate flux through membrane 1
J2 = zeros(1, numsteps); % Permeate flux through membrane 2
J1s = zeros(1, numsteps); % Salt flux through membrane 1
J2s = zeros(1, numsteps); % Salt flux through membrane 2
Hres = zeros(1, numsteps); % Height of reservoir
Chigh1 = zeros(1, numsteps); % Concentration at beginning of flushing stage in membrane 1
Chigh2 = zeros(1, numsteps); % Concentration at beginning of flushing stage in membrane 2
Cf1arr = ones(1, membrane_cell_count); % Feedside concentration array in membrane 1
Cf2arr = ones(1, membrane_cell_count); % Feedside concentration array in membrane 2

%% Totals
feed_water_drawn = 0; % volume of feed water used [m^3]
fresh_water_produced = 0; % volume of fresh water produced [m^3]
%brine_produced_mem1 = 0; % volume of brine produced from membrane 1 [m^3]
%brine_produced_mem2 = 0; % volume of brine produced from membrane 2 [m^3]
brine_produced = 0; % total brine produced
flushcount1 = 0; % amount of water flushed from membrane 1 since start of flushing cycle
flushcount2 = 0; % amount of water flushed from membrane 2 since start of flushing cycle

%% Misc programming variables
b = 1;
n = 1; % timepoint
switching = false; % switch states during runtime?


%% System states
states = struct;
states.buoy_up = true; % direction buoy is moving
states.mem1_permeate_step = true;  % State of membrane 1 (permeate production is true, flushing is false)
states.mem2_permeate_step = true; % State of membrane 2 (permeate production is true, flushing is false)
states.restricted_flow = true; % Buoy restriction

%% Bounds
bounds1 = struct;% this structure stores the bounds for the parameters
bounds1.S_min = 0;% [g/kg] minimum salinity
bounds1.S_max = 260;% [g/kg] maximum salinity
bounds1.P_min = 0;% [bar g] minimum pressure
bounds1.P_max = 120;% [bar g] maximum pressure
bounds1.Q_min = 0;% [m/s] minimum flow rate
bounds1.Q_max = 1.4/3600;% [m^3/s] maximum flow rate
%bounds1.Q_max = mem1_feed_qmax;
bounds1.u_min = bounds1.Q_min/module1.crossArea;% [m/s] minimum velocity
bounds1.u_max = bounds1.Q_max/module1.crossArea;% [m/s] maximum velocity
bounds1.Js_min = 0; %[g/m^2/s] minimum salt flux
bounds1.Js_max = module1.B_s*bounds1.S_max*polyval(prop.polyRho_S,bounds1.S_max); %[g/m^2/s] maximum salt flux
bounds1.Jw_min = 0; %[m/s] minimum water flux
bounds1.Jw_max = module1.L_p*bounds1.P_max; %[m/s] maximum water flux

bounds2 = struct;% this structure stores the bounds for the parameters
bounds2.S_min = 0;% [g/kg] minimum salinity
bounds2.S_max = 260;% [g/kg] maximum salinity
bounds2.P_min = 0;% [bar g] minimum pressure
bounds2.P_max = 120;% [bar g] maximum pressure
bounds2.Q_min = 0;% [m/s] minimum flow rate
bounds2.Q_max = 1.4/3600;% [m^3/s] maximum flow rate
%bounds2.Q_max = mem2_feed_qmax;
bounds2.u_min = bounds2.Q_min/module2.crossArea;% [m/s] minimum velocity
bounds2.u_max = bounds2.Q_max/module2.crossArea;% [m/s] maximum velocity
bounds2.Js_min = 0; %[g/m^2/s] minimum salt flux
bounds2.Js_max = module2.B_s*bounds2.S_max*polyval(prop.polyRho_S,bounds2.S_max); %[g/m^2/s] maximum salt flux
bounds2.Jw_min = 0; %[m/s] minimum water flux
bounds2.Jw_max = module2.L_p*bounds2.P_max; %[m/s] maximum water flux


%% Initialize time grid
for b = 1:numsteps
    t(b) = (b - 1) * (total_time / (numsteps - 1));
end



%% Initialize constants (may become variables in the future)
for b = 1:numsteps
    %Ce1(b) = mem1_exit_salt_concentration;
    %Ce2(b) = mem2_exit_salt_concentration;
    Pep1(b) = mem1_exit_pressure;
    Pep2(b) = mem2_exit_pressure;
end



%% INITIAL CONDITIONS at t = 0 (buoy in equilibrium) Start at restricted pp upward
y(1) = 0;
ydot(1) = 0;
% ydot(1) = piston1_velocity_max_pp;
% yn1start = y(1)-c.delt*ydot(1);
% yn2start = y(1)-2*c.delt*ydot(1);
% yn3start = y(1)-3*c.delt*ydot(1);
yn1start = 0;
yn2start = 0;
yn3start = 0;

Fb(1) = 0;
ydotdot(1) = 0;
F1(1) = 0;
F2(1) = 0;
yp2(1) = 0;
yp2dot(1) = 0;
% yp2dot(1) = -piston1_velocity_max_pp*lever_length_piston2_side/lever_length_piston1_side;
yp2dotdot(1) = 0;
% yp2n1start = yp2(1)-c.delt*yp2dot(1);
% yp2n2start = yp2(1)-2*c.delt*yp2dot(1);
% yp2n3start =yp2(1)-3*c.delt*yp2dot(1);
yp2n1start = 0;
yp2n2start = 0;
yp2n3start = 0;
Hres(1) = reservoir_start_height;
Pp1(1) = hydrostatic_pressure - (rho*g*piston1_equilibrium_height)*10^(-5);
Pm1(1) = hydrostatic_pressure - (rho*g*(segment2_elevation - segment1_elevation))*10^(-5);
Pp2(1) = Pm1(1) + c.rho*c.g*(-c.z2mem1-c.Lo2-c.yeqp2)*10^(-5);
Pm2(1) = Pp2(1) - rho*g*(segment4_elevation+c.Lo2+c.yeqp2)*10^-5;
Cf1(1) = initial_salinity;
%Cf1(1) = 70;
Ce1(1) = 0.001;
Cres(1) = Ce1(1);
Cf2(1) = Cres(1);
%Cf2(1) = 0.04;
Ce2(1) = 0.001;
theta(1) = 0;
J1(1) = 0;
J2(1) = 0;
J1s(1) = 0;
J2s(1) = 0;
Chigh1(1) = Cf1(1);
Chigh2(1) = Cres(1);
Cf1arr = Cf1arr * Cf1(1);
Cf2arr = Cf2arr * Cf2(1);


%% Calculate equilibrium height
yeqstop = 15
y(yeqstop) = 20;
count = 0;
cutoffa = 0.00001
while abs(y(yeqstop)) > cutoffa && count < 120
    count = count+1
    n = 2;
    [y(n), ydot(n), ydotdot(n), yp2(n), yp2dot(n), yp2dotdot(n), Fb(n), F1(n), F2(n), Pp1(n), Pm1(n), Pp2(n), Pm2(n), ...
         Cf1(n), Cf2(n), Ce1(n), Ce2(n), Cres(n), theta(n), thetadot(n), thetadotdot(n), J1(n), J2(n), J1s(n), J2s(n), Hres(n),Chigh1(n),Chigh2(n),states.buoy_up,states.restricted_flow,Cf1arr,Cf2arr] ...
        = butter_churner(n, wh(n), ... known variables
        0, 0, 0, Fb(n-1), F1(n-1), F2(n-1), Pp1(n-1), Pm1(n-1), Pp2(n-1), Pm2(n-1), ...
        Cf1(n-1), Cf2(n-1), Ce1(n-1), Ce2(n-1), Cres(n-1), 0, J1(n-1), J2(n-1), J1s(n-1), J2s(n-1), Hres(n-1), Chigh1(n-1), Chigh2(n-1), ... Variables at previous time step
        Cf1arr,Cf2arr, ...
        0, 0, 0, 0, 0, 0, ... Extra old info past original time step
        c, states, module1, module2, bounds1, bounds2, prop, ...
        Pep1(n), Pep2(n));
    
    n = 3;
    [y(n), ydot(n), ydotdot(n), yp2(n), yp2dot(n), yp2dotdot(n), Fb(n), F1(n), F2(n), Pp1(n), Pm1(n), Pp2(n), Pm2(n), ...
         Cf1(n), Cf2(n), Ce1(n), Ce2(n), Cres(n), theta(n), thetadot(n), thetadotdot(n), J1(n), J2(n), J1s(n), J2s(n), Hres(n),Chigh1(n),Chigh2(n),states.buoy_up,states.restricted_flow,Cf1arr,Cf2arr] ...
        = butter_churner(n, wh(n), ... known variables
        y(n-1), ydotdot(n-1), yp2(n-1), Fb(n-1), F1(n-1), F2(n-1), Pp1(n-1), Pm1(n-1), Pp2(n-1), Pm2(n-1), ...
        Cf1(n-1), Cf2(n-1), Ce1(n-1), Ce2(n-1), Cres(n-1), theta(n-1), J1(n-1), J2(n-1), J1s(n-1), J2s(n-1), Hres(n-1), Chigh1(n-1), Chigh2(n-1), ... Variables at previous time step
        Cf1arr,Cf2arr,...
        0, 0, 0, 0, 0, 0, ... Extra old info past original time step
        c, states, module1, module2, bounds1, bounds2, prop, ...
        Pep1(n), Pep2(n));

    n = 4;
    [y(n), ydot(n), ydotdot(n), yp2(n), yp2dot(n), yp2dotdot(n), Fb(n), F1(n), F2(n), Pp1(n), Pm1(n), Pp2(n), Pm2(n), ...
         Cf1(n), Cf2(n), Ce1(n), Ce2(n), Cres(n), theta(n), thetadot(n), thetadotdot(n), J1(n), J2(n), J1s(n), J2s(n), Hres(n),Chigh1(n),Chigh2(n),states.buoy_up,states.restricted_flow,Cf1arr,Cf2arr] ...
        = butter_churner(n, wh(n), ... known variables
        y(n-1), ydotdot(n-1), yp2(n-1), Fb(n-1), F1(n-1), F2(n-1), Pp1(n-1), Pm1(n-1), Pp2(n-1), Pm2(n-1), ...
        Cf1(n-1), Cf2(n-1), Ce1(n-1), Ce2(n-1), Cres(n-1), theta(n-1), J1(n-1), J2(n-1), J1s(n-1), J2s(n-1), Hres(n-1), Chigh1(n-1), Chigh2(n-1), ... Variables at previous time step
        Cf1arr,Cf2arr,...
        y(n-2), 0, yp2(n-2), 0, theta(n-2), 0, ... Extra old info past original time step
        c, states, module1, module2, bounds1, bounds2, prop, ...
        Pep1(n), Pep2(n));

    for n = 5:yeqstop
    [y(n), ydot(n), ydotdot(n), yp2(n), yp2dot(n), yp2dotdot(n), Fb(n), F1(n), F2(n), Pp1(n), Pm1(n), Pp2(n), Pm2(n), ...
         Cf1(n), Cf2(n), Ce1(n), Ce2(n), Cres(n), theta(n), thetadot(n), thetadotdot(n), J1(n), J2(n), J1s(n), J2s(n), Hres(n),Chigh1(n),Chigh2(n),states.buoy_up,states.restricted_flow,Cf1arr,Cf2arr] ...
        = butter_churner(n, wh(n), ... known variables
        y(n-1), ydotdot(n-1), yp2(n-1), Fb(n-1), F1(n-1), F2(n-1), Pp1(n-1), Pm1(n-1), Pp2(n-1), Pm2(n-1), ...
        Cf1(n-1), Cf2(n-1), Ce1(n-1), Ce2(n-1), Cres(n-1), theta(n-1), J1(n-1), J2(n-1), J1s(n-1), J2s(n-1), Hres(n-1), Chigh1(n-1), Chigh2(n-1), ... Variables at previous time step
        Cf1arr,Cf2arr,...
        y(n-2), y(n-3), yp2(n-2), yp2(n-3), theta(n-2), theta(n-3), ... Extra old info past original time step
        c, states, module1, module2, bounds1, bounds2, prop, ...
        Pep1(n), Pep2(n));
    end
    c.yeq = c.yeq-y(n);
    % if count > 300
    %     c.yeq = c.omega*(c.yeq-y(n))+(1-c.omega)*(c.yeq);
    %     c.omega = 0.8-(0.001*count);
    % end
    
    c.yeq
end


% % test flushing
% states.mem1_permeate_step = false;
% states.mem2_permeate_step = false;

%% Initialize wave position
for b = 1:numsteps
    wh(b) = (wave_height / 2) * sin(2 * pi / wave_period * t(b)); % WAVE PROFILE (Change if needed)
end

%% startup iterations

n = 2;
[y(n), ydot(n), ydotdot(n), yp2(n), yp2dot(n), yp2dotdot(n), Fb(n), F1(n), F2(n), Pp1(n), Pm1(n), Pp2(n), Pm2(n), ...
     Cf1(n), Cf2(n), Ce1(n), Ce2(n), Cres(n), theta(n), thetadot(n), thetadotdot(n), J1(n), J2(n), J1s(n), J2s(n), Hres(n),Chigh1(n),Chigh2(n),states.buoy_up,states.restricted_flow,Cf1arr,Cf2arr] ...
    = butter_churner(n, wh(n), ... known variables
    y(n-1), ydotdot(n-1), yp2(n-1), Fb(n-1), F1(n-1), F2(n-1), Pp1(n-1), Pm1(n-1), Pp2(n-1), Pm2(n-1), ...
    Cf1(n-1), Cf2(n-1), Ce1(n-1), Ce2(n-1), Cres(n-1), theta(n-1), J1(n-1), J2(n-1), J1s(n-1), J2s(n-1), Hres(n-1), Chigh1(n-1), Chigh2(n-1), ... Variables at previous time step
    Cf1arr,Cf2arr, ...
    yn1start, yn2start, yp2n1start, yp2n2start, theta(n-1), theta(n-1), ... Extra old info past original time step
    c, states, module1, module2, bounds1, bounds2, prop, ...
    Pep1(n), Pep2(n));
n;
states;
feed_water_drawn = feed_water_drawn + max(0,yp2dot(n)*c.Ap2/c.delt);
fresh_water_produced = fresh_water_produced + states.mem2_permeate_step * J2(n)*c.Am2/c.delt;
brine_produced = brine_produced + (1-states.mem2_permeate_step) * max(0,(1-c.Yrec)*-ydot(n)*c.Ap1/c.delt) + (1-states.mem1_permeate_step) * max(0,-yp2dot(n)*c.Ap2/c.delt);    
flushcount1 = flushcount1 + (1-states.mem1_permeate_step) * max(0,-yp2dot(n)*c.Ap2/c.delt);
flushcount2 = flushcount2 + (1-states.mem2_permeate_step) * max(0,(1-c.Yrec)*-ydot(n)*c.Ap1/c.delt);


n = 3;
[y(n), ydot(n), ydotdot(n), yp2(n), yp2dot(n), yp2dotdot(n), Fb(n), F1(n), F2(n), Pp1(n), Pm1(n), Pp2(n), Pm2(n), ...
     Cf1(n), Cf2(n), Ce1(n), Ce2(n), Cres(n), theta(n), thetadot(n), thetadotdot(n), J1(n), J2(n), J1s(n), J2s(n), Hres(n),Chigh1(n),Chigh2(n),states.buoy_up,states.restricted_flow,Cf1arr,Cf2arr] ...
    = butter_churner(n, wh(n), ... known variables
    y(n-1), ydotdot(n-1), yp2(n-1), Fb(n-1), F1(n-1), F2(n-1), Pp1(n-1), Pm1(n-1), Pp2(n-1), Pm2(n-1), ...
    Cf1(n-1), Cf2(n-1), Ce1(n-1), Ce2(n-1), Cres(n-1), theta(n-1), J1(n-1), J2(n-1), J1s(n-1), J2s(n-1), Hres(n-1), Chigh1(n-1), Chigh2(n-1), ... Variables at previous time step
    Cf1arr,Cf2arr,...
    y(n-2), yn1start, yp2(n-2), yp2n1start, theta(n-2), theta(n-2), ... Extra old info past original time step
    c, states, module1, module2, bounds1, bounds2, prop, ...
    Pep1(n), Pep2(n));
n;
states;
feed_water_drawn = feed_water_drawn + max(0,yp2dot(n)*c.Ap2/c.delt);
fresh_water_produced = fresh_water_produced + states.mem2_permeate_step * J2(n)*c.Am2/c.delt;
brine_produced = brine_produced + (1-states.mem2_permeate_step) * max(0,(1-c.Yrec)*-ydot(n)*c.Ap1/c.delt) + (1-states.mem1_permeate_step) * max(0,-yp2dot(n)*c.Ap2/c.delt);   
flushcount1 = flushcount1 + (1-states.mem1_permeate_step) * max(0,-yp2dot(n)*c.Ap2/c.delt);
flushcount2 = flushcount2 + (1-states.mem2_permeate_step) * max(0,(1-c.Yrec)*-ydot(n)*c.Ap1/c.delt);


n = 4;
[y(n), ydot(n), ydotdot(n), yp2(n), yp2dot(n), yp2dotdot(n), Fb(n), F1(n), F2(n), Pp1(n), Pm1(n), Pp2(n), Pm2(n), ...
     Cf1(n), Cf2(n), Ce1(n), Ce2(n), Cres(n), theta(n), thetadot(n), thetadotdot(n), J1(n), J2(n), J1s(n), J2s(n), Hres(n),Chigh1(n),Chigh2(n),states.buoy_up,states.restricted_flow,Cf1arr,Cf2arr] ...
    = butter_churner(n, wh(n), ... known variables
    y(n-1), ydotdot(n-1), yp2(n-1), Fb(n-1), F1(n-1), F2(n-1), Pp1(n-1), Pm1(n-1), Pp2(n-1), Pm2(n-1), ...
    Cf1(n-1), Cf2(n-1), Ce1(n-1), Ce2(n-1), Cres(n-1), theta(n-1), J1(n-1), J2(n-1), J1s(n-1), J2s(n-1), Hres(n-1), Chigh1(n-1), Chigh2(n-1), ... Variables at previous time step
    Cf1arr,Cf2arr,...
    y(n-2), y(n-3), yp2(n-2), yp2(n-3), theta(n-2), theta(n-3), ... Extra old info past original time step
    c, states, module1, module2, bounds1, bounds2, prop, ...
    Pep1(n), Pep2(n));
n;
states;
feed_water_drawn = feed_water_drawn + max(0,yp2dot(n)*c.Ap2/c.delt);
fresh_water_produced = fresh_water_produced + states.mem2_permeate_step * J2(n)*c.Am2/c.delt;
brine_produced = brine_produced + (1-states.mem2_permeate_step) * max(0,(1-c.Yrec)*-ydot(n)*c.Ap1/c.delt) + (1-states.mem1_permeate_step) * max(0,-yp2dot(n)*c.Ap2/c.delt);   
flushcount1 = flushcount1 + (1-states.mem1_permeate_step) * max(0,-yp2dot(n)*c.Ap2/c.delt);
flushcount2 = flushcount2 + (1-states.mem2_permeate_step) * max(0,(1-c.Yrec)*-ydot(n)*c.Ap1/c.delt);

for n = 5:numsteps
    [y(n), ydot(n), ydotdot(n), yp2(n), yp2dot(n), yp2dotdot(n), Fb(n), F1(n), F2(n), Pp1(n), Pm1(n), Pp2(n), Pm2(n), ...
     Cf1(n), Cf2(n), Ce1(n), Ce2(n), Cres(n), theta(n), thetadot(n), thetadotdot(n), J1(n), J2(n), J1s(n), J2s(n), Hres(n),Chigh1(n),Chigh2(n),states.buoy_up,states.restricted_flow,Cf1arr,Cf2arr] ...
    = butter_churner(n, wh(n), ... known variables
    y(n-1), ydotdot(n-1), yp2(n-1), Fb(n-1), F1(n-1), F2(n-1), Pp1(n-1), Pm1(n-1), Pp2(n-1), Pm2(n-1), ...
    Cf1(n-1), Cf2(n-1), Ce1(n-1), Ce2(n-1), Cres(n-1), theta(n-1), J1(n-1), J2(n-1), J1s(n-1), J2s(n-1), Hres(n-1), Chigh1(n-1), Chigh2(n-1), ... Variables at previous time step
    Cf1arr,Cf2arr,...
    y(n-2), y(n-3), yp2(n-2), yp2(n-3), theta(n-2), theta(n-3), ... Extra old info past original time step
    c, states, module1, module2, bounds1, bounds2, prop, ...
    Pep1(n), Pep2(n));
    n
    states;
    feed_water_drawn = feed_water_drawn + max(0,yp2dot(n)*c.Ap2/c.delt);
    fresh_water_produced = fresh_water_produced + states.mem2_permeate_step * J2(n)*c.Am2/c.delt;
    brine_produced = brine_produced + (1-states.mem2_permeate_step) * max(0,(1-c.Yrec)*-ydot(n)*c.Ap1/c.delt) + (1-states.mem1_permeate_step) * max(0,-yp2dot(n)*c.Ap2/c.delt);
    flushcount1 = flushcount1 + (1-states.mem1_permeate_step) * max(0,-yp2dot(n)*c.Ap2/c.delt)
    flushcount2 = flushcount2 + (1-states.mem2_permeate_step) * max(0,(1-c.Yrec)*-ydot(n)*c.Ap1/c.delt);

    %% state logic
    if t(n) > 48
        states.mem1_permeate_step = false;
        states.mem2_permeate_step = false;
    end
    if switching && states.mem1_permeate_step && Cf1(n) > 70
        states.mem1_permeate_step = false;
        flushcount1 = 0;
    end
    if switching && states.mem1_permeate_step == false && flushcount1 > mem1_volume && ydot(n)<0
    %if switching &&  states.mem1_permeate_step == false && Cf1(n) < 40
    
        states.mem1_permeate_step = true;
    end
    if switching && states.mem2_permeate_step && Cf2(n) > 5
        states.mem2_permeate_step = false;
        flushcount2 = 0;
    end
    if switching && states.mem2_permeate_step == false && flushcount2 > mem2_volume && ydot(n)>0
    %if states.mem2_permeate_step == false && Cf2(n) < 0.003
        states.mem2_permeate_step = true;
    end
end
recovery_ratio = fresh_water_produced / feed_water_drawn
brine_produced
fresh_water_produced
feed_water_drawn

%% CALL ALGORITHM 1



%% Figures & debugging
figure (1)
plot(t,wh,'r')
hold on
plot(t,y,'b')
plot(t,yp2,'g')
title("Buoy and Piston Position Plot")
legend("wave height","buoy/piston 1","piston 2")
xlabel("time [s]")
ylabel("height [m]")

figure (2)
plot(t,ydot,'r')
hold on
plot(t,yp2dot,'b')
title("Piston Velocity Plot")
legend("buoy/piston 1","piston 2")
xlabel("time [s]")
ylabel("velocity [m/s]")

figure (3)
plot(t,Fb,'r')
hold on
title("Internal Forces")
plot(t,F1,'b')
plot(t,F2,'g')
legend("Buoy - Piston 1","Piston 1 - Lever","Piston 2 - Lever")
xlabel("time [s]")
ylabel("force [N]")

figure(4)
plot(t,J1,'r')
hold on
title("Membrane Fluxes")
plot(t,J2,'g')
plot(t,J1s,'k')
plot(t,J2s,'b')
legend("membrane 1 water flux","membrane 2 water flux","membrane 1 salt flux","membrane 2 salt flux")
xlabel("time [s]")
ylabel("water flux [m3/m2s], salt flux [kg/m2s]")

figure (5)
plot(t,Cf1,'r')
hold on
plot(t,Ce1,'g')
title("Concentration in Membrane 1")
xlabel("time [s]")
legend("Feed", "Permeate")
ylabel("concentration [g/L]")
ylim([0,80])



figure (6)
plot(t,Cf2,'b')
hold on
plot(t,Ce2,'k')
title("Concentration in Membrane 2")
legend("Feed","Permeate")
ylim([0,0.0025])
xlabel("time [s]")
ylabel("concentration [g/L]")


figure(7)
plot(t,Hres,'r')
title("Reservoir Height")
xlabel("time [s]")
ylabel("height [m]")


figure(8)
plot(t,Cres,'r')
title("Reservoir Concentration")
xlabel("time [s]")
ylabel("concentration [g/L]")

figure(9)
plot(t,Pm1,'r')
hold on
plot(t,Pm2,'b')
plot(t,Pp1,'k')
plot(t,Pp2,'g')
title("Pressure at Various Locations")
legend("membrane 1 inlet","membrane 2 inlet","piston 1 face","piston 2 face")
xlabel("time [s]")
ylabel("pressure [bar]")

