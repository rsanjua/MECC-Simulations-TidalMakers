function [S_carr_next,P_warr,Q_warr,Js_carr,Jw_carr,Sm_carr,Sp_carr,S_warr,Delta_carr]=PFRO_WEC2(S_in,Q_in,S_carr,dt,OperationMode,N,module,bounds,prop,method)
if Q_in~=0
    switch lower(OperationMode)
        case "permeateproduction"
            Q_out = 0;% blocking the brine
            [P_warr,Q_warr,Js_carr,Jw_carr,Sm_carr,Sp_carr,Delta_carr]=DynamicsPermeateProduction(S_carr,Q_in,Q_out,N,module,bounds,prop);
        case "flushing"
            P_out = 0;% during flushing, the RO outlet is at atmospheric pressure
            [P_warr,Q_warr,Js_carr,Jw_carr,Sm_carr,Sp_carr,Delta_carr]=DynamicsFlushing(S_carr,Q_in,P_out,N,module,bounds,prop);
    end
    dx = module.memLength/N;
    Rho_carr = polyval(prop.polyRho_S,S_carr);% [kg/m^3] density
    rho_0 = polyval(prop.polyRho_S,0);% [kg/m^3] pure water density    
    S_warr = zeros(N+1,1);    
    Q_carr = (Q_warr(1:end-1)+Q_warr(2:end))/2;
    switch lower(method)
        case "upwind"
% using first order upwind
            S_warr(1) = S_in;% inlete condition
            S_warr(2) = S_carr(1);
            S_warr(3:end) = S_carr(2:N);
        case "beamwarming"
% using beam warming method
            S_warr(1) = S_in;% inlete condition
            S_warr(2) = 2*S_carr(1)-S_in;
            S_warr(3:end) = (3*S_carr(2:N)-S_carr(1:N-1))/2;
        case "dispersion"
            S_warr(1) = S_in;% inlete condition
            DispCoef= DispersionCoeff(Q_carr,S_carr,module,prop);
            for jj=2:N+1
                S_warr(jj)=S_warr(jj-1)

            end
            S_warr(2:end) = (3*S_carr(2:N)-S_carr(1:N-1))/2;
    end
    dSdx = (S_warr(2:end)-S_warr(1:end-1))/dx;
    S_carr_next=S_carr-dt*((module.memLength/module.Vol_mem).*Q_carr.*dSdx+(module.memArea/module.Vol_mem).*(Js_carr-S_carr.*(Js_carr/1000+rho_0*Jw_carr))./Rho_carr);
else % if we dont have any inlet flow rate:
    S_carr_next = S_carr;
    P_warr = zeros(N+1,1);
    Q_warr = zeros(N+1,1);
    Js_carr = zeros(N,1);
    Jw_carr = zeros(N,1);
    Sm_carr = zeros(N,1);
    Sp_carr = zeros(N,1);
    S_warr = zeros(N+1,1);
    switch lower(method)
        case "upwind"
% using first order upwind
            S_warr(1) = S_in;% inlete condition
            S_warr(2) = S_carr(1);
            S_warr(3:end) = S_carr(2:N);
        case "beamwarming"
% using beam warming method
            S_warr(1) = S_in;% inlete condition
            S_warr(2) = 2*S_carr(1)-S_in;
            S_warr(3:end) = (3*S_carr(2:N)-S_carr(1:N-1))/2;
    end
end
end

function [P_warr,Q_warr,Js_carr,Jw_carr,Sm_carr,Sp_carr,Delta_carr]=DynamicsFlushing(S_carr,Q_in,P_out,N,module,bounds,prop)
    funchandle=@(X) SystemOfEquationsFlushing(X,S_carr,Q_in,P_out,N,module,bounds,prop);
    X0=zeros(N,1);
    X0(1:N) = Unbounding(1,bounds.P_min,bounds.P_max); 
 
    options=optimoptions('fsolve','FunctionTolerance',1e-9,'OptimalityTolerance',1e-9);
    [IndepVar,fval,exitflag,~] = fsolve(funchandle,X0,options);
    P_warr = [Bounding(IndepVar(1:N),bounds.P_min,bounds.P_max) ; P_out];                   %[bar g] pressure (wall based array)
    Q_warr = Q_in*ones(N+1,1); %[m^3/s] flow rate at the walls
    Jw_carr = zeros(N,1);          %[m/s] water flux at the cells 
    Js_carr = zeros(N,1);          %[g/m^2/s] salt flux at the cells 
    Sm_carr = S_carr;
    Sp_carr = S_carr;
    Delta_carr = zeros(N,1);

end
function [P_warr,Q_warr,Js_carr,Jw_carr,Sm_carr,Sp_carr,Delta_carr]=DynamicsPermeateProduction(S_carr,Q_in,Q_out,N,module,bounds,prop)
    funchandle=@(X) SystemOfEquationsPermeateProduction(X,S_carr,Q_in,Q_out,N,module,bounds,prop);
    X0=zeros(4*N,1);
    X0(1:N+1) = Unbounding(polyval(prop.polyPi_S,S_carr(end)),bounds.P_min,bounds.P_max); 
    Q_linear=linspace(Q_in,Q_out,N+1);
    X0(N+2:2*N) = Unbounding(Q_linear(2:end-1),bounds.Q_min,bounds.Q_max);
    X0(2*N+1:3*N) = [Unbounding(Q_in/module.memArea,bounds.Jw_min,bounds.Jw_max)];          %[m/s] water flux at the cells 
    X0(3*N+1:4*N) = [Unbounding(module.B_s*(S_carr.*polyval(prop.polyRho_S,S_carr)),bounds.Js_min,bounds.Js_max)];          %[g/m^2/s] salt flux at the cells 
    options=optimoptions('fsolve','FunctionTolerance',1e-9,'OptimalityTolerance',1e-9);
    [IndepVar,fval,exitflag,~] = fsolve(funchandle,X0,options);
    P_warr = Bounding(IndepVar(1:N+1),bounds.P_min,bounds.P_max);                   %[bar g] pressure (wall based array)
    Q_warr = [Q_in ; Bounding(IndepVar(N+2:2*N),bounds.Q_min,bounds.Q_max) ;Q_out]; %[m^3/s] flow rate at the walls
    Jw_carr = [Bounding(IndepVar(2*N+1:3*N),bounds.Jw_min,bounds.Jw_max)];          %[m/s] water flux at the cells 
    Js_carr = [Bounding(IndepVar(3*N+1:4*N),bounds.Js_min,bounds.Js_max)];          %[g/m^2/s] salt flux at the cells 
    rho_0 = polyval(prop.polyRho_S,0); % [kg/m^3] pure water density
    dx = module.memLength/N;% [m] length of each cell
    x_carr = ((1:N)'-0.5)*dx;
    Q_carr = (Q_warr(2:end)+Q_warr(1:end-1))/2;% [m^3/s] cell value for flow rate
    Sp_carr = Js_carr./(Js_carr/1000+Jw_carr*rho_0); %[g/kg] permeate salinity
    k_carr = MassTransferCoeff(Q_carr,S_carr,x_carr,module,prop);% [m/s] mass transfer coefficient

    D = polyval(prop.polyD_S,S_carr);
    Delta_carr = D./k_carr;% thickness of mass transfer boundary layer. Equation from https://doi.org/10.1021/acs.est.8b02771


    Sm_carr = (S_carr-Sp_carr).*exp(Jw_carr./k_carr)+Sp_carr;% [g/kg] salinity near active layer
end


function EQS=SystemOfEquationsFlushing(IndepVar,S_carr,Q_in,P_out,N,module,bounds,prop)
% this function contains the governing equations for the flushing
% % Input arguments:
% IndepVar: Independent variables for the system of equations
% S_carr : [g/kg] salinity (cell based array)
% N : number of cells
    Q_carr = Q_in*ones(N,1);
    P_warr = [Bounding(IndepVar(1:N),bounds.P_min,bounds.P_max) ; P_out];                   %[bar g] pressure (wall based array)
    dx = module.memLength/N;% [m] length of each cell
    PressureDropEq = (P_warr(2:end)-P_warr(1:end-1))/dx+PressureDropRate(Q_carr,S_carr,module,prop);% N equations
% scaling the equations so that they are in the same order of magnitude
    PressureDropEq = PressureDropEq./PressureDropRate(bounds.Q_max,bounds.S_max,module,prop);
    % EQS = [WaterConservEq(:) ; PressureDropEq(:) ; WaterFluxEq(:) ; SaltFluxEq(:)];
    EQS = PressureDropEq(:);
end
function EQS=SystemOfEquationsPermeateProduction(IndepVar,S_carr,Q_in,Q_out,N,module,bounds,prop)
% this function contains the governing equations for the permeate production
% % Input arguments:
% IndepVar: Independent variables for the system of equations
% S_carr : [g/kg] salinity (cell based array)
% N : number of cells
    P_warr = Bounding(IndepVar(1:N+1),bounds.P_min,bounds.P_max);                   %[bar g] pressure (wall based array)
    Q_warr = [Q_in ; Bounding(IndepVar(N+2:2*N),bounds.Q_min,bounds.Q_max) ;Q_out]; %[m^3/s] flow rate at the walls
    Jw_carr = [Bounding(IndepVar(2*N+1:3*N),bounds.Jw_min,bounds.Jw_max)];          %[m/s] water flux at the cells 
    Js_carr = [Bounding(IndepVar(3*N+1:4*N),bounds.Js_min,bounds.Js_max)];          %[g/m^2/s] salt flux at the cells 
    Rho_carr = polyval(prop.polyRho_S,S_carr);% [kg/m^3] density
    DerRho_carr = polyval(prop.polyDerRho_S,S_carr);% [kg^2/g/m^3]
    rho_0 = polyval(prop.polyRho_S,0); % [kg/m^3] pure water density
    dx = module.memLength/N;% [m] length of each cell
    x_carr = ((1:N)'-0.5)*dx;
    P_carr = (P_warr(2:end)+P_warr(1:end-1))/2;% [bar g] cell value for pressure  
    Q_carr = (Q_warr(2:end)+Q_warr(1:end-1))/2;% [m^3/s] cell value for flow rate
    Sp_carr = Js_carr./(Js_carr/1000+Jw_carr*rho_0); %[g/kg] permeate salinity
    Rho_p_carr = polyval(prop.polyRho_S,Sp_carr);% [kg/m^3] permeate density
    k_carr = MassTransferCoeff(Q_carr,S_carr,x_carr,module,prop);% [m/s] mass transfer coefficient
    D = polyval(prop.polyD_S,S_carr);
    Delta_carr = D./k_carr;% thickness of mass transfer boundary layer. Equation from https://doi.org/10.1021/acs.est.8b02771
    %Sm_carr = (S_carr-Sp_carr).*exp(Jw_carr./k_carr)+Sp_carr;% [g/kg] salinity near active layer
    H_ch = module.H_ch;
    Sm_carr=(exp((Jw_carr.*Delta_carr.*rho_0)./(D.*Rho_carr)).*((H_ch.*S_carr)./2 - (Js_carr.*(3.*Jw_carr.^3.*Delta_carr.^2.*rho_0.^3 + 6.*D.^2.*Jw_carr.*Rho_carr.^2.*rho_0 - 6.*D.^2.*Jw_carr.*Rho_carr.^2.*rho_0.*exp((Jw_carr.*Delta_carr.*rho_0)./(D.*Rho_carr)) + 6.*D.*Jw_carr.^2.*Delta_carr.*Rho_carr.*rho_0.^2))./(H_ch.*Jw_carr.^4.*Rho_carr.*rho_0.^3) + (Js_carr.*(12.*D.^3.*Rho_carr.^3 - 12.*D.^3.*Rho_carr.^3.*exp((Jw_carr.*Delta_carr.*rho_0)./(D.*Rho_carr)) + 2.*Jw_carr.^3.*Delta_carr.^3.*rho_0.^3 + 6.*D.*Jw_carr.^2.*Delta_carr.^2.*Rho_carr.*rho_0.^2 + 12.*D.^2.*Jw_carr.*Delta_carr.*Rho_carr.^2.*rho_0))./(H_ch.^2.*Jw_carr.^4.*Rho_carr.*rho_0.^3)))./(H_ch./2 - (- 2.*Delta_carr.^3 + 3.*H_ch.*Delta_carr.^2)./H_ch.^2 + (6.*D.*Rho_carr.*(2.*D.^2.*Rho_carr.^2 - 2.*D.^2.*Rho_carr.^2.*exp((Jw_carr.*Delta_carr.*rho_0)./(D.*Rho_carr)) + Jw_carr.^2.*Delta_carr.^2.*rho_0.^2 + 2.*D.*Jw_carr.*Delta_carr.*Rho_carr.*rho_0))./(H_ch.^2.*Jw_carr.^3.*rho_0.^3) - (6.*D.*Rho_carr.*(Jw_carr.^2.*Delta_carr.*rho_0.^2 + D.*Jw_carr.*Rho_carr.*rho_0 - D.*Jw_carr.*Rho_carr.*rho_0.*exp((Jw_carr.*Delta_carr.*rho_0)./(D.*Rho_carr))))./(H_ch.*Jw_carr.^3.*rho_0.^3)) - (Js_carr.*(exp((Jw_carr.*Delta_carr.*rho_0)./(D.*Rho_carr)) - 1))./(Jw_carr.*Rho_carr);
    Rho_m_carr = polyval(prop.polyRho_S,Sm_carr);% density near membrane active layer
    Pi_m_carr = polyval(prop.polyPi_S,Sm_carr);% [bar] osmotic pressure near the active layer
    Pi_p_carr = polyval(prop.polyPi_S,Sp_carr); % [bar] osmotic pressure in the permeate
    %WaterConservEq = (Q_warr(2:end)-Q_warr(1:end-1))/dx+Jw_carr*module.memArea/module.memLength;% N equations
    MassConservEq = Rho_carr.*(Q_warr(2:end)-Q_warr(1:end-1))/dx+(module.memArea/module.memLength)*((Js_carr/1000+rho_0*Jw_carr).*(1+S_carr.*DerRho_carr./Rho_carr)-Js_carr.*DerRho_carr./Rho_carr);
    PressureDropEq = (P_warr(2:end)-P_warr(1:end-1))/dx+PressureDropRate(Q_carr,S_carr,module,prop);% N equations
    WaterFluxEq = Jw_carr-module.L_p.*(P_carr - (Pi_m_carr-Pi_p_carr));% N equations
    MassConservEq = MassConservEq./(polyval(prop.polyRho_S,bounds.S_max)*bounds.Q_max./module.memLength);
    SaltFluxEq = Js_carr-module.B_s.*(Rho_m_carr.*Sm_carr-Rho_p_carr.*Sp_carr);% N equations
% scaling the equations so that they are in the same order of magnitude
    %WaterConservEq = WaterConservEq./(bounds.Jw_max.*module.memArea./module.memLength);
    PressureDropEq = PressureDropEq./PressureDropRate(bounds.Q_max,bounds.S_max,module,prop);
    WaterFluxEq = WaterFluxEq./bounds.Jw_max;
    SaltFluxEq = SaltFluxEq./bounds.Js_max;
    % EQS = [WaterConservEq(:) ; PressureDropEq(:) ; WaterFluxEq(:) ; SaltFluxEq(:)];
    EQS = [MassConservEq(:) ; PressureDropEq(:) ; WaterFluxEq(:) ; SaltFluxEq(:)];
    disp(max(WaterFluxEq(:)))
end
function DispCoef =DispersionCoeff(Q,S,module,prop)
Gamma=1;
D = polyval(prop.polyD_S,S);
L_h=0.1;
DispCoef = Gamma.*(2./105).*((module.dHydraulic).^2./(D.*(module.crossArea)^2)).*(Q.^2)+(L_h^2/module.Vol_mem).*Q;
end


function k = MassTransferCoeff(Q,S,x,module,prop)
% calculates the mass transfer coefficient
% correlation from (https://doi.org/10.1016/j.seppur.2022.122121) valid for Reynolds number between 20 and 150. Our Re
% might be lower than 20, but since it is still laminar, we assume the
% correlation holds for Re number less than 20 as well
% input arguments
% S: [g/kg] salinity
% Q: [m^3/s] volumetric flow rate.
% output arguments
% x: [m] position from the module entrance
% k : [m/s] mass transfer coefficient
Re = Reynolds(Q,S,module,prop);
Sc = Schmidt(S,prop);
Sh = 2.401.*((Re.*Sc).^0.297).*(x./module.dHydraulic).^(-0.279);  
Sh(Sh<2*module.dHydraulic/module.H_ch)=2*module.dHydraulic/module.H_ch;% for making sure, the boundary layer thickness is not more than half of the channel height

D = polyval(prop.polyD_S,S);
k = Sh.*D./module.dHydraulic;

end
function dPdx = PressureDropRate(Q,S,module,prop)
% calculates the pressure drop rate
% input arguments
% S: [g/kg] salinity
% Q: [m^3/s] volumetric flow rate.
% output argument:
% dPdx : [bar/m] pressure drop rate
Re=Reynolds(Q,S,module,prop);
f=FrictionFactor(Re);
rho=polyval(prop.polyRho_S,S);
dPdx=f.*rho.*(Q.^2)./(2.*module.dHydraulic.*module.crossArea.^2)./1e5;
end
function f = FrictionFactor(Re)
% calculates the friction factor
% the correlation is from https://doi.org/10.1016/j.jwpe.2019.100820 which
% is for Re<250 and for 28 mil spacer with 45 degree angel between the
% spacer filaments and the from direction
% % Input arguments:
% Re: Reynolds number
f=63.*Re.^(-0.64);% this considers laminar flow
end
function Pe = Peclet(S,Q,prop,module)
Sc = Schmidt(S,prop);
Re = Reynolds(Q,S,module,prop);
Pe = Re.*Sc;
end
function Sc = Schmidt(S,prop)
% calculates the Schmidt number
% % Input arguments:
% S: [g/kg] salinity
Nu=polyval(prop.polyNu_S,S);% [m^2/s] kinematic viscosity
D=polyval(prop.polyD_S,S);% [m^2/s] diffusion coefficient
Sc=Nu./D;
end
function Re = Reynolds(Q,S,module,prop)
% calculates the Reynolds number
% % Input arguments:
% Q: [m^3/s] volumetric flow rate.
% S: [g/kg] salinity
Nu=polyval(prop.polyNu_S,S);% [m^2/s] kinematic viscosity
u=Q./module.crossArea;% [m/s] cross flow velocity
Re=u.*module.dHydraulic./Nu;
end
function y = Bounding(x,lb,ub)
% converts the parameter x without any bounds to the parameter y with lb
% and ub as lower and upper bounds
% % input arguments:
% x: Unbounded variable
% lb: lower bound
% ub: upper bound
    y=lb+(ub-lb).*(tanh(x)+1)./2;
end
function y = Unbounding(x,lb,ub)
% converts parameter x with lb and ub as lower and upper bounds to the parameter y without any bounds
% lb: lower bound
% ub: upper bound
% x: Bounded variable
    y=atanh(2.*(x-lb)./(ub-lb)-1);
end