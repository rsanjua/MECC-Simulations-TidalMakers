function [y, ydot, ydotdot, yp2, yp2dot, yp2dotdot, Fb, F1, F2, Pp1, Pm1, Pp2, Pm2, ...
     Cf1, Cf2, Ce1, Ce2, Cres, theta, thetadot, thetadotdot, J1, J2, J1s, J2s, Hres, Chigh1, Chigh2, buoy_up, flow_restricted,Cf1arr,Cf2arr] ...
    = butter_churner(t, wh, ... known variables
    yn1, ydotdotn1, yp2n1, Fbn1, F1n1, F2n1, Pp1n1, Pm1n1, Pp2n1, Pm2n1, ...
    Cf1n1, Cf2n1, Ce1n1, Ce2n1, Cresn1, thetan1, J1n1, J2n1, J1sn1, J2sn1, Hresn1, Chigh1n1, Chigh2n1, ... Variables at previous time step
    Cf1arrn1,Cf2arrn1,...
    yn2, yn3, yp2n2, yp2n3, thetan2, thetan3, ... Extra old info past original time step
    c, states, module1, module2, bounds1, bounds2, prop,...
    Pep1, Pep2)
% is ryan cool?
ryan_is_cool = true;

% old values for residual and SOR
yold = yn1;
ydotdotold = ydotdotn1;
yp2old = yp2n1;
thetaold = thetan1;
phi2old = asind((c.L2-c.w2*cosd(thetaold))/c.H2);
phi1old = asind((c.L1-c.w1*cosd(thetaold))/c.H1);
Fbold = Fbn1;
F1old = F1n1;
F1xold = F1old*sind(phi1old);
F1yold = F1old*cosd(phi1old);
F2old = F2n1;
F2xold = F2old*sind(phi2old);
F2yold = F2old*sind(phi2old);

Pp1old = Pp1n1;
Pm1old = Pm1n1;
Pp2old = Pp2n1;
Pm2old = Pm2n1;
Cf1old = Cf1n1;
Cf2old = Cf2n1;
Ce1old = Ce1n1;
Ce2old = Ce2n1;
Cresold = Cresn1;
J1old = J1n1;
J1sold = J1sn1;
J2old = J2n1;
J2sold = J2sn1;
Hresold = Hresn1;


residual = 1000;
iterationcount = 0;
residual_cutoff = 0.0000005;
iteration_limit = 500;

buoy_up = states.buoy_up;
flow_restricted = states.restricted_flow;

%% LETS CHURN SOME BUTTER %%
while residual > residual_cutoff && iterationcount < iteration_limit
    residual = 0;
    % if flow_restricted == false
    % 
    %     % 1. Calculate y from buoy equation
    %     y = buoyy(wh,ydotdotold,Fbold,c);
    %     y = SOR(y,yold,c.omega);
    %     yold = y;
    % 
    %     % 2. Calculate all the geometries
    %     [~, ydot, ydotdot, theta, thetadot, thetadotdot, yp2, yp2dot, yp2dotdot,phi1,phi2] = geometry(y,0,yn1,yn2,yn3,0,0,yp2n1,yp2n2,yp2n3,thetan1,thetan2,thetan3,c,"y");
    %     if ydot > 0
    %         buoy_up = true;
    %     else
    %         buoy_up = false;
    %     end
    % else
    %     if buoy_up == true
    % 
    %         % 1. Calculate max velocity from constraint
    %         if states.mem1_permeate_step
    %             yp2dot = -c.vmaxp2;
    %         else
    %             yp2dot = -c.vmaxf2;
    %         end
    % 
    %         % 2. Calculate all the geometries
    %         [y, ydot, ydotdot, theta, thetadot, thetadotdot, yp2, ~, yp2dotdot,phi1,phi2] = geometry(0,0,yn1,yn2,yn3,0,yp2dot,yp2n1,yp2n2,yp2n3,thetan1,thetan2,thetan3,c,"yp2dot");
    %     else
    %         % 1. Calculate max velocity from constraint
    %         if states.mem2_permeate_step
    %             ydot = -c.vmaxp1;
    %         else
    %             ydot = -c.vmaxf1;
    %         end
    % 
    %         % 2. Calculate all the geometries
    %         [y, ~, ydotdot, theta, thetadot, thetadotdot, yp2, yp2dot, yp2dotdot,phi1,phi2] = geometry(0,ydot,yn1,yn2,yn3,0,0,yp2n1,yp2n2,yp2n3,thetan1,thetan2,thetan3,c,"ydot");
    % 
    %     end
    % end
    % 1.  calculate y and y dot
    y = buoyy(wh,ydotdotold,Fbold,c);
    %y = SOR(y,yold,c.omega);
    ydot = (3*y - 4*yn1 + yn2)/(2*c.delt);
    
    % 2. bound y dot and recalculate all geometries
    if states.mem1_permeate_step
        ub = c.vmaxnegp1;
    else
        ub = c.vmaxnegf1;
    end
    if states.mem2_permeate_step
        lb = -c.vmaxp1;
    else
        lb = -c.vmaxf1;
    end
    ydot = Bounding(ydot,lb,ub);
    [y, ~, ydotdot, theta, thetadot, thetadotdot, yp2, yp2dot, yp2dotdot,phi1,phi2] = geometry(0,ydot,yn1,yn2,yn3,0,0,yp2n1,yp2n2,yp2n3,thetan1,thetan2,thetan3,c,"ydot");
    
    if ydot>0
        buoy_up = true;
    else
        buoy_up = false;
    end
    
    %residual = residual + abs((y-yold)/y)
    yold = y;

    % 3. Calculate forces 1 and 2
    Pamb2 = c.rho*c.g*(c.Wdepth-c.heqp2-yp2-c.tp2)*10^-5;
    F2y = c.mp2*yp2dotdot + c.mp2*c.g + Pamb2*c.A2R*10^5 - Pp2old*c.Ap2*10^5;
    F2 = F2y/cosd(phi2);
    F2x = F2*sind(phi2);
    F2x = SOR(F2x,F2xold,c.omega);
    F2y = SOR(F2y,F2yold,c.omega);
    F2 = SOR(F2,F2old,c.omega);
    %residual = residual + abs((F2-F2old)/F2)
    F2old = F2;
    F2xold = F2x;
    F2yold = F2y;

    F1 = ((c.I*thetadotdot) + (c.w2*F2y*cosd(theta)) + (c.w2*F2x*sind(theta)) + (((c.w1+c.w2)/2-c.w1)*c.mlev*c.g*cosd(theta))) / (c.w1*(cosd(theta)*cosd(phi1)-sind(theta)*sind(phi1)));
    F1x = F1*sind(phi1);
    F1y = F1*cosd(phi1);
    F1x = SOR(F1x,F1xold,c.omega);
    F1y = SOR(F1y,F1yold,c.omega);
    F1 = SOR(F1,F1old,c.omega);
    %residual = residual + abs((F1-F1old)/F1)
    F1old = F1;
    F1xold = F1x;
    F1yold = F1y;

    % % 3.1 calculate buoy force if equation hasn't been used (I think i might
    % do this at the end to check if we need a switch back)
    % if flow_restricted
    %     Fb = -c.mb*c.g-c.mb*ydotdot+c.rho*c.g*pi/3*((wh-y+c.yeq+c.r/tand(c.phib))^3 * (tand(c.phib)*tand(c.phib)) - (r^3/tand(c.phib)));
    %     Fb = SOR(Fb,Fbold,c.omega);
    %     Fbold = Fb;
    % end

    % 4. Calculate flow rate

    if (ryan_is_cool)
        % 4. Calculate fluxes
        if states.mem1_permeate_step
            J1 = max(-yp2dot*c.Ap2/c.Am1,0);
            J1 = SOR(J1,J1old,c.omega);
            if buoy_up
                J1s = c.Ps1*(Cf1old-Ce1old);
                J1s = SOR(J1s,J1sold,c.omega);
            else
                J1s = 0;
            end
            %residual = residual + abs((J1-J1old)/J1)
            %residual = residual + abs((J1s-J1sold)/J1s) % theres gonna be some funky zero case here
        else
            J1 = 0;
            J1s = 0;
        end
        if states.mem2_permeate_step
            J2 = max(-ydot*c.Ap1/c.Am2,0);
            J2 = SOR(J2,J2old,c.omega);
            if buoy_up == false
                J2s = c.Ps2*(Cf2old-Ce2old);
                J2s = SOR(J2s,J2sold,c.omega);
            else
                J2s = 0;
            end
            %residual = residual + abs((J2-J2old)/J2)
            %residual = residual + abs((J2s-J2sold)/J2s)
        else
            J2 = 0;
            J2s = 0;
        end
    
    
        J1old = J1;
        J1sold = J1s;
        J2old = J2;
        J2sold = J2s;
    
        % 5. Calculate Pressures (all pressures in bar)
        if buoy_up
            Vm1 = -c.Ap2*yp2dot/c.Am1in;
            if states.mem1_permeate_step
                Pm1 = J1/c.Pw1+Pep1*10^5+ c.iv*(Cf1old-Ce1old)*c.R*c.T;
                Pm1 = Pm1*10^-5;
                Pm1 = SOR(Pm1,Pm1old,c.omega);
            else
                Pm1 = c.ff1*c.rho*Vm1*Vm1*c.Lm1/c.Dm1 + c.Pb1*10^5;
                Pm1 = Pm1*10^-5;
                Pm1 = SOR(Pm1,Pm1old,c.omega);
            end
            if states.mem2_permeate_step
                Pm2 = Pep2*10^5+ c.iv*(Cf2old-Ce2old)*c.R*c.T;
                Pm2 = Pm2*10^-5;
                Pm2 = SOR(Pm2,Pm2old,c.omega);
            else
                Pm2 = c.Pb2*10^5;
                Pm2 = Pm2*10^-5;
                Pm2 = SOR(Pm2,Pm2old,c.omega);
            end
            %residual = residual + abs((Pm1-Pm1old)/Pm1);
            %residual = residual + abs((Pm2-Pm2old)/Pm2);
    
            Pm1old = Pm1;
            Pm2old = Pm2;
    
    
            Vin1 = max(ydot*c.Ap1/c.Apipe, 0);   %%%%%%%%%%% pick up heere
            Rein1 = c.rho*Vin1*c.Dpipe/c.mu + 0.0001; % adding .0001 because of the stability case
            frick = 64/Rein1; % assume laminar flow for now, don't know how to type a moody diagram
            H_Lin1 = (5 + frick*c.Ls3/c.Dpipe)*(Vin1*Vin1/2/c.g);
            Pp1 = c.Pres + c.rho*c.g*(Hresold-H_Lin1-c.alpha3*Vin1*Vin1/2/c.g-c.z1mem1-c.heqp1+c.z1-y)*10^-5;
            Pp1 = SOR(Pp1,Pp1old,c.omega);
    
            Vout2 = -yp2dot*c.Ap2/c.Apipe;
            Reout2 = c.rho*Vout2*c.Dpipe/c.mu;
            frick = 64/Reout2;
            H_Lout2 = (5 + frick*c.Ls2/c.Dpipe) * (Vout2*Vout2/2/c.g);
            Pp2 = Pm1old - c.rho*c.g*(c.alpha2*Vout2*Vout2/2/c.g + c.z2mem1 + c.heqp2 - c.z2 + yp2 - H_Lout2 - c.alpha2*Vm1*Vm1/2/c.g)*10^-5;
            Pp2 = SOR(Pp2,Pp2old,c.omega);
        else
            Vm2 = -c.Ap1*ydot/c.Am2in;
            if states.mem1_permeate_step
                Pm1 = Pep1*10^5+ c.iv*(Cf1old-Ce1old)*c.R*c.T;
                Pm1 = Pm1*10^-5;
                Pm1 = SOR(Pm1,Pm1old,c.omega);
            else
                Pm1 = c.Pb1*10^5;
                Pm1 = Pm1*10^-5;
                Pm1 = SOR(Pm1,Pm1old,c.omega);
            end
            if states.mem2_permeate_step
                Pm2 = J2old/c.Pw2+Pep2*10^5+ c.iv*(Cf2old-Ce2old)*c.R*c.T;
                Pm2 = Pm2*10^-5;
                Pm2 = SOR(Pm2,Pm2old,c.omega);
            else
                Pm2 = c.ff1*c.rho*Vm2*Vm2*c.Lm2/c.Dm2 + c.Pb2*10^5;
                Pm2 = Pm2*10^-5;
                Pm2 = SOR(Pm2,Pm2old,c.omega);
            end
            %residual = residual + abs((Pm1-Pm1old)/Pm1);
            %residual = residual + abs((Pm2-Pm2old)/Pm2);
            Pm1old = Pm1;
            Pm2old = Pm2;
    
            Vout1 = -ydot*c.Ap1/c.Apipe;
            Reout1 = c.rho*Vout1*c.Dpipe/c.mu;
            frick = 64/Reout1;
            H_Lout1 = (5+frick*c.Ls4/c.Dpipe) * (Vout1*Vout1/2/c.g);
            Pp1 = Pm2old - c.rho*c.g*(c.alpha4*Vout1*Vout1/2/c.g + c.z1mem2 + c.heqp1 - c.z1 + y - H_Lout1-c.alpha4*Vm2*Vm2/2/c.g)*10^-5;
            Pp1 = SOR(Pp1,Pp1old,c.omega);
    
            Vin2 = yp2dot*c.Ap2/c.Apipe;
            Rein2 = c.rho*Vin2*c.Dpipe/c.mu;
            frick = 64/Rein2;
            H_Lin2 = (5.3 + frick*c.Ls1/c.Dpipe) * (Vin2*Vin2/2/c.g);
            Pp2 = c.Pinlet + c.rho*c.g*(-H_Lin2-c.alpha1*Vin2*Vin2/2/c.g - c.heqp2 - yp2)*10^-5;
            Pp2 = SOR(Pp2,Pp2old,c.omega);
        end
        %residual = residual + abs((Pp1-Pp1old)/Pp1);
        %residual = residual + abs((Pp2-Pp2old)/Pp2);
    
        Pp1old = Pp1;
        Pp2old = Pp2;
    
        % 6. Calculate all concentrations
        if buoy_up
            if states.mem1_permeate_step
                Ce1 = J1s/J1;
                Ce1 = SOR(Ce1,Ce1old,c.omega);
                Cf1 = -c.Ap2*yp2dot*c.delt*(c.Cseawater-Ce1)/c.Vfeed1 + Cf1n1;
                Cf1 = SOR(Cf1,Cf1old,c.omega);
                Chigh1 = Cf1;
    
                Ce2 = 0;
                Cf2 = Cf2n1;
                Chigh2 = Cf2;
            else
                Ce1 = 0;
                % Ce1 = SOR(Ce1,Ce1old,c.omega);
                Chigh1 = Chigh1n1;
                Cf1 = -c.Ap2*yp2dot*c.delt*(c.Cseawater-Chigh1)/c.Vfeed1 + Cf1n1;
                Cf1 = SOR(Cf1,Cf1old,c.omega);
    
                Ce2 = 0;
                Cf2 = Cf2n1;
                Chigh2 = Cf2;
            end
        else
            if states.mem2_permeate_step
                Ce1 = 0;
                Cf1 = Cf1n1;
                Chigh1 = Cf1;
    
                Ce2 = J2s/J2;
                Ce2 = SOR(Ce2,Ce2n1,c.omega);
                Cf2 = -c.Ap1*ydot*c.delt*(Cresold-Ce2)/c.Vfeed2 + Cf2n1;
                Cf2 = SOR(Cf2,Cf2old,c.omega);
                Chigh2 = Cf2;
            else
                Ce1 = 0;
                Cf1 = Cf1n1;
                Chigh1 = Cf1;
    
                Ce2 = 0;
                Chigh2 = Chigh2n1;
                Cf2 = -c.Ap1*ydot*c.delt*(Cresold-Chigh2)/c.Vfeed2 + Cf2n1;
                Cf2 = SOR(Cf2,Cf2old,c.omega);
    
            end
        end
        %residual = residual + abs((Cf1-Cf1old)/Cf1);
        %residual = residual + abs((Cf2-Cf2old)/Cf2);
    
        Cf1old = Cf1;
        Cf2old = Cf2;
        Cf1arr = ones(1,c.N)*Cf1;
        Cf2arr = ones(1,c.N)*Cf2;

    else
        % 4. Run Ali's code PFRO_WEC2
        if states.mem1_permeate_step
            OperationMode = "permeateproduction";
        else
            OperationMode = "flushing";
        end
        Q_in1 = max(-yp2dot*c.Ap2, 0);
        [Cf1arr,Pm1_warr,~,J1s_carr,J1_carr,~,Ce1_carr,~,~]=PFRO_WEC2(c.Cseawater,Q_in1,Cf1arrn1,c.delt,OperationMode,c.N,module1,bounds1,prop,"upwind");
        if states.mem2_permeate_step
            OperationMode = "permeateproduction";
        else
            OperationMode = "flushing";
        end
        Q_in2 = max(-ydot*c.Ap1,0);
        [Cf2arr,Pm2_warr,~,J2s_carr,J2_carr,~,Ce2_carr,~,~]=PFRO_WEC2(Cresold,Q_in2,Cf2arrn1,c.delt,OperationMode,c.N,module2,bounds2,prop,"upwind");
        
        % 5. Transform arrays into quantities
        Cf1 = sum(Cf1arr,"all")*c.delx1/c.Lm1;
        Cf2 = sum(Cf2arr,"all")*c.delx2/c.Lm2;
        Ce1 = sum(Ce1_carr,"all")*c.delx1/c.Lm1;
        Ce2 = sum(Ce2_carr,"all")*c.delx2/c.Lm2;        
        J1 = sum(J1_carr,"all")*c.delx1/c.Lm1;
        J2 = sum(J2_carr,"all")*c.delx2/c.Lm2;
        J1s = sum(J1s_carr,"all")*c.delx1/c.Lm1;
        J2s = sum(J2s_carr,"all")*c.delx2/c.Lm2;
        Pm1 = sum(Pm1_warr,"all")*c.delx1/c.Lm1;
        Pm2 = sum(Pm2_warr,"all")*c.delx2/c.Lm2;


        % 6. Calculate pressure at the piston faces
        Vm1 = max(-c.Ap2*yp2dot/c.Am1in,0);
        Vm2 = max(-c.Ap1*ydot/c.Am2in,0);
        if buoy_up    

            Vin1 = max(ydot*c.Ap1/c.Apipe, 0);   %%%%%%%%%%% pick up heere
            Rein1 = c.rho*Vin1*c.Dpipe/c.mu + 0.0001; % adding .0001 because of the stability case
            frick = 64/Rein1; % assume laminar flow for now, don't know how to type a moody diagram
            H_Lin1 = (5 + frick*c.Ls3/c.Dpipe)*(Vin1*Vin1/2/c.g);
            Pp1 = c.Pres + c.rho*c.g*(Hresold-H_Lin1-c.alpha3*Vin1*Vin1/2/c.g-c.z1mem1-c.heqp1+c.z1-y)*10^-5;
            Pp1 = SOR(Pp1,Pp1old,c.omega);
    
            Vout2 = max(-yp2dot*c.Ap2/c.Apipe,0);
            Reout2 = c.rho*Vout2*c.Dpipe/c.mu + 0.0001;
            frick = 64/Reout2;
            H_Lout2 = (5 + frick*c.Ls2/c.Dpipe) * (Vout2*Vout2/2/c.g);
            Pp2 = Pm1old - c.rho*c.g*(c.alpha2*Vout2*Vout2/2/c.g + c.z2mem1 + c.heqp2 - c.z2 + yp2 - H_Lout2 - c.alpha2*Vm1*Vm1/2/c.g)*10^-5;
            Pp2 = SOR(Pp2,Pp2old,c.omega);
        else
    
            Vout1 = max(-ydot*c.Ap1/c.Apipe,0);
            Reout1 = c.rho*Vout1*c.Dpipe/c.mu + 0.0001;
            frick = 64/Reout1;
            H_Lout1 = (5+frick*c.Ls4/c.Dpipe) * (Vout1*Vout1/2/c.g);
            Pp1 = Pm2old - c.rho*c.g*(c.alpha4*Vout1*Vout1/2/c.g + c.z1mem2 + c.heqp1 - c.z1 + y - H_Lout1-c.alpha4*Vm2*Vm2/2/c.g)*10^-5;
            Pp1 = SOR(Pp1,Pp1old,c.omega);
    
            Vin2 = max(yp2dot*c.Ap2/c.Apipe,0);
            Rein2 = c.rho*Vin2*c.Dpipe/c.mu+0.0001;
            frick = 64/Rein2;
            H_Lin2 = (5.3 + frick*c.Ls1/c.Dpipe) * (Vin2*Vin2/2/c.g);
            Pp2 = c.Pinlet + c.rho*c.g*(-H_Lin2-c.alpha1*Vin2*Vin2/2/c.g - c.heqp2 - yp2)*10^-5;
            Pp2 = SOR(Pp2,Pp2old,c.omega);
        end
        %residual = residual + abs((Pp1-Pp1old)/Pp1);
        %residual = residual + abs((Pp2-Pp2old)/Pp2);
        Pp1old = Pp1;
        Pp2old = Pp2;
        Chigh1 = Chigh1n1;
        Chigh2 = Chigh2n1;
    end

    % 7. Calculate Reservoir things
    if buoy_up
        if states.mem1_permeate_step
            Hres = c.delt*(J1*c.Am1-ydot*c.Ap1)/c.Ares + Hresn1;
            Hres = SOR(Hres,Hresold,c.omega);
            Cres = c.delt*(J1*c.Am1*Ce1 - ydot*c.Ap1*Cresn1)/c.Ares/Hres + Cresn1;
            Cres = SOR(Cres,Cresold,c.omega);

        else
            Hres = -c.delt*ydot*c.Ap1/c.Ares + Hresn1;
            Hres = SOR(Hres,Hresold,c.omega);
            Cres = -c.delt*ydot*c.Ap1*Cresn1/c.Ares/Hres + Cresn1;
            Cres = SOR(Cres,Cresold,c.omega);

        end
    else
        if states.mem2_permeate_step
            Hres = Hresn1;
            Cres = Cresn1;

        else
            Hres = -c.Yrec*ydot*c.Ap1*c.delt/c.Ares + Hresn1;
            Hres = SOR(Hres,Hresold,c.omega);
            Cres = -c.Yrec*ydot*c.Ap1*Chigh2*c.delt/c.Ares/Hres + Cresn1;
            Cres = SOR(Cres,Cresold,c.omega);

        end
    end
    %residual = residual + abs((Hres-Hresold)/Hres);
    %residual = residual + abs((Cres-Cresold)/Cres);
    Hresold = Hres;
    Cresold = Cres;

    % 8. Calculate buoy force
    Pamb1 = c.rho*c.g*(c.Wdepth-c.heqp1-y-c.tp1)*10^-5;
    Fb = c.mp1*ydotdot-Pp1*c.Ap1*10^5 - F1y + c.mp1*c.g + Pamb1*c.A1R*10^5;
    Fb = SOR(Fb,Fbold,c.omega);
    residual = residual + abs(Fb-Fbold);
    Fbold = Fb;

    % 9. Check the logic
    % if flow_restricted == false
    %     if states.mem1_permeate_step && -yp2dot > c.vmaxp2
    %         flow_restricted = true;
    %     elseif states.mem1_permeate_step == false && -yp2dot > c.vmaxf2
    %         flow_restricted = true;
    %     elseif states.mem2_permeate_step && -ydot > c.vmaxp1
    %         flow_restricted = true;
    %     elseif states.mem2_permeate_step == false && -ydot > c.vmaxf1
    %         flow_restricted = true;
    %     end
    % 
    % else
    % 
    %     ytest = buoyy(wh,ydotdot,Fb,c);
    %     [~, ydottest, ~, ~, ~, ~, ~, yp2dottest, ~,~,~] = geometry(ytest,0,yn1,yn2,yn3,0,0,yp2n1,yp2n2,yp2n3,thetan1,thetan2,thetan3,c,"y");
    % 
    %     if ydottest<0
    %         buoy_up = false;
    %     else
    %         buoy_up = true;
    %     end
    %     if states.mem1_permeate_step && buoy_up && abs(yp2dottest)<c.vmaxp2
    %         flow_restricted = false;
    %     elseif states.mem1_permeate_step == false && buoy_up && abs(yp2dottest)<c.vmaxf2
    %         flow_restricted = false;
    %     elseif states.mem2_permeate_step && buoy_up == false && abs(ydottest) < c.vmaxp1
    %         flow_restricted = false;
    %     elseif states.mem2_permeate_step == false && buoy_up == false && abs(ydottest) < c.vmaxf1
    %         flow_restricted = false;
    %     end
    % 
    % end
iterationcount = iterationcount + 1;
residual;

end % loop end
% ytest
% ydottest
% yp2dottest
end

function [xnew] = SOR(xbar, xold, omega)
% This function implements successive over/underrelaxation to help with
% iterative convergence
    xnew = xbar*omega + (1-omega)*xold;
end

function [y] = buoyy(wh,ydotdot,Fb,c)
(c.mb*ydotdot + c.mb*c.g + Fb);
y = wh + c.yeq + c.r/tand(c.phib) - nthroot(((((c.mb*ydotdot + c.mb*c.g + Fb) / (c.rho*c.g)) + ((pi*c.r^3) / (3*tand(c.phib))))+(3/(pi*tand(c.phib)*tand(c.phib)))),(3));
    
end

function [y, ydot, ydotdot, theta, thetadot, thetadotdot, yp2, yp2dot, yp2dotdot,phi1,phi2] = geometry(y,ydot,yn1,yn2,yn3,yp2,yp2dot,yp2n1,yp2n2,yp2n3,thetan1,thetan2,thetan3,c,known)
% this function calculates the positional variables and their derivatives
% of the entire system based on previous information and one piece of
% current information noted by the variable known
    switch known
        case "y"
            % y = y;
            ydot = (3*y - 4*yn1 + yn2)/(2*c.delt); % ydot
            theta = thetacalcy(y,c.w1, c.Lc1, c.L1); % 17a
            yp2 = c.w2*sind(theta) + c.Lc2*(1-sqrt(1-((c.L2-c.w2*cosd(theta))/c.H2)^2)); % 17b
            yp2dot = (3*yp2 - 4*yp2n1 + yp2n2)/(2*c.delt); % yp2dot
    
        case "ydot"
            y = (2*c.delt*ydot + 4*yn1 - yn2)/3; %ydot
            % ydot = ydot;
            theta = thetacalcy(y,c.w1, c.Lc1, c.L1); % 17a
            yp2 = c.w2*sind(theta) + c.Lc2*(1-sqrt(1-((c.L2-c.w2*cosd(theta))/c.H2)^2)); % 17b
            yp2dot = (3*yp2 - 4*yp2n1 + yp2n2)/(2*c.delt); % yp2dot
        
        case "yp2dot"
            
            yp2 = (2*c.delt*yp2dot + 4*yp2n1 - yp2n2)/3;
            theta = thetacalcy(-yp2,c.w2, c.Lc2, c.L2); % 17a
            y = -1*c.w1*sind(theta) + c.Lc1*(1-sqrt(1-((c.Lc1-c.w1*cosd(theta))/c.H1)^2));
            % yp2dot = (3*yp2 - 4*yp2n1 + yp2n2)/(2*c.delt) % yp2dot
            ydot = (3*y - 4*yn1 + yn2)/(2*c.delt);
            
    end
   
    ydotdot = (2*y - 5*yn1 + 4*yn2 - yn3) / (c.delt^2); % ydotdot
    thetadot = (3*theta - 4*thetan1 + thetan2)/(2*c.delt); % thetadot
    yp2dotdot = (2*yp2 - 5*yp2n1 + 4*yp2n2 - yp2n3) / (c.delt^2); % yp2dotdot
    thetadotdot = (2*theta - 5*thetan1 + 4*thetan2 - thetan3) / (c.delt^2); % thetadotdot
    phi2 = asind((c.L2-c.w2*cosd(theta))/c.H2);
    phi1 = asind((c.L1-c.w1*cosd(theta))/c.H1);

end

function [theta] = thetacalcy(y,w1, Lc1, L1)
%thetacalcy Calculate theta when piston 1 position is known
%   The function uses the bracketing root finding method to calculate theta
%   from y. can also calculate from yp2, but sign of w1 should be negative.
% Calculate theta using root finding method (bracketing)
% thetamax = 45;
% thetamin = -45;
% thetaguess = 0;
% rootval = 100;
% cutoff = 0.000005;
% y;
% while abs(rootval) > cutoff
%     rootval = -w1*sind(thetaguess) + Lc1*(1-sqrt(1-((L1-w1*cosd(thetaguess))/H1)^2)) - y;
%     if rootval < -cutoff
%         thetamax = thetaguess;
%         thetaguess = (thetamax + thetamin) /2;
%     elseif rootval > cutoff
%         thetamin = thetaguess;
%         thetaguess = (thetamax + thetamin) /2;
%     end
% end
% 
% 
% 
% theta = thetaguess;
% end

% Calculate theta using Newton-Raphson method
cutoff = 0.0000005;
iterationcap = 2000;
theta = 0;
oldtheta = theta;
iteration = 0;
residual = cutoff+1;
while residual > cutoff && iteration < iterationcap
    
    ftheta = -w1*sind(theta) + Lc1*(1-sqrt(1-((L1-abs(w1)*cosd(theta))/Lc1)^2))-y;
    fprimetheta = (abs(w1)*sind(theta)*(L1-abs(w1)*cosd(theta)))/(Lc1*sqrt(1-((L1-abs(w1)*cosd(theta))/Lc1)^2))-w1*cosd(theta);
    oldtheta = theta;
    theta = theta - ftheta/fprimetheta;
    residual = abs((theta-oldtheta)/oldtheta);
    iteration = iteration + 1;
end
end


%% CODE BELOW WAS WRITTEN BY ALI NADERI BENI
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
                S_warr(jj)=S_warr(jj-1);

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
    Delta_carr = zeros(N,1);
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
 
    options=optimoptions('fsolve','FunctionTolerance',1e-9,'OptimalityTolerance',1e-9,'MaxIterations',10);
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
    options=optimoptions('fsolve','FunctionTolerance',1e-9,'OptimalityTolerance',1e-9,'MaxIterations',10);
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