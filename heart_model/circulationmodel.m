%
% This function implements the heart model described in 
% Kim, et al, Ann Biomed Eng 2009;37(11):2153-69
% coupled to an RCR circuit (as opposed to a 3D FEM model in Kim's paper).
%

function circulationmodel

global p_la L_av L_va v_lv0 v_d0 R_c R_d C R_AV R_VA

% Heart Parameters (from Kim, et al, Ann Biomed Eng 2009;37(11):2153-69.)
% All parameters in cgs units
R_AV = 0;
R_VA = 0;
L_av = 0.67; % atrioventricular inductance (mitral valve open)
L_va = 0.69; % ventricular-aortic inductance (aortic valve open)
v_lv0 = -33; % ventricular volume at "zero pressure"
% note, above volume not physical, it's the (extraplated) y-intercept of
% the Elastance lines (i.e., linear elastance only hold locally)
p_la = 12*1333.2; % left atrial pressure (1mmHg = 1333.2 dyne/cm^2)
v_lv_i = 130; % assumed initial amount of blood in the LV

% Arterial Parameters (Westerhof, et al, J Appl Physiol 1971;31:776-81)
R_c = 90; % proximal (characteristic) resistance
R_d = 1200; % distal resistance
C = 8e-4; % capacitance
v_d_i = 1000; % assumed initial amount of blood in the arteries
v_d0 = 0.93*v_d_i; % unpressurized arterial blood volume

% Number of cycles for integration
num_cycles = 10;

% Initialize structures to store the solution history and final cycle
soln = struct('t', 0, 'Q_m', 0, 'v_lv', v_lv_i, 'Q_a', 0, 'v_d', v_d_i);
final = struct('t', [], 'Q_m', [], 'v_lv', [], 'Q_a', [], 'v_d', []);

for cycle=1:num_cycles
    
    % Integrate ODE's for isovolumic contraction phase
    options = odeset('Events', @aortic_open, 'Maxstep', 0.01);
    [t_soln, y_soln, t_event, y_event] = ode45( ...
        @(t,y)dydt_isovolumic(t, y), ...
        [soln.t(end) num_cycles], ...
        [0 soln.v_lv(end) 0 soln.v_d(end)], ...
        options);
    soln = UpdateSoln(soln, t_soln, t_event, y_soln, y_event);
    if cycle==num_cycles
        final = FinalSoln(final, t_soln, t_event, y_soln, y_event);
    end
    
    % Integrate ODE's for ejection phase
    options = odeset('Events', @aortic_close, 'Maxstep', 0.01);
    [t_soln, y_soln, t_event, y_event] = ode45( ...
        @(t,y)dydt_ejection(t, y), ...
        [soln.t(end) num_cycles], ...
        [0 soln.v_lv(end) 0 soln.v_d(end)], ...
        options);
    soln = UpdateSoln(soln, t_soln, t_event, y_soln, y_event);
    if cycle==num_cycles
        final = FinalSoln(final, t_soln, t_event, y_soln, y_event);
    end
    
    elas = elastance(t_event);
    fprintf(1,'Aortic valve closing with pLV - pA = %f\n', ...
        elas(1)*(y_event(2) - v_lv0) - 1 / C * (y_event(4) - v_d0) - y_event(3) * R_c);
    
    % Integrate ODE's for isovolumic relaxation phase
    options = odeset('Events', @mitral_open, 'Maxstep', 0.01);
    [t_soln, y_soln, t_event, y_event] = ode45( ...
        @(t,y)dydt_isovolumic(t, y), ...
        [soln.t(end) num_cycles], ...
        [0 soln.v_lv(end) 0 soln.v_d(end)], ...
        options);
    soln = UpdateSoln(soln, t_soln, t_event, y_soln, y_event);
    if cycle==num_cycles
        final = FinalSoln(final, t_soln, t_event, y_soln, y_event);
    end
    
    % Integrate ODE's for filling phase
    options = odeset('Events', @mitral_close, 'Maxstep', 0.01);
    [t_soln, y_soln, t_event, y_event] = ode45( ...
        @(t,y)dydt_filling(t, y), ...
        [soln.t(end) num_cycles], ...
        [0 soln.v_lv(end) 0 soln.v_d(end)], ...
        options);
    soln = UpdateSoln(soln, t_soln, t_event, y_soln, y_event);
    if cycle==num_cycles
        final = FinalSoln(final, t_soln, t_event, y_soln, y_event);
    end
    
    elas = elastance(t_event);
    fprintf(1, 'Mitral valve closing with pLA - pLV = %f\n', ...
        p_la - elas(1)*(y_event(2) - v_lv0));
    
end

p_lv_final = elastance(final.t) .* (final.v_lv - v_lv0) / 1333.2;
p_d_final = 1 / C * (final.v_d - v_d0) / 1333.2;

work = polyarea(final.v_lv, 1333.2*p_lv_final) / 1e7;
display(work);

figure(1);
plot(final.v_lv, p_lv_final); hold on;
title('Pressure vs. volume');
xlabel('volume (ml)'); ylabel('pressure (mmHg)'); hold off;
figure(2);
plot(final.t - final.t(1), p_lv_final); hold on;
plot(final.t - final.t(1), p_d_final, ':');
title('Pressure vs. time (last cycle)');
legend('p_{lv}', 'p_d');
xlabel('time (sec)'); ylabel('pressure (mmHg)'); hold off;
figure(3);
plot(final.t - final.t(1), final.Q_m); hold on;
plot(final.t - final.t(1), final.Q_a, ':');
title('Flow vs. time  (last cycle)');
legend('Q_m', 'Q_a');
xlabel('time (sec)'); ylabel('flow (ml/s)'); hold off;
figure(4);
plot(final.t - final.t(1), elastance(final.t))
title('Elastance');
xlabel('time (sec)'); ylabel('Elastance'); hold off;

function [soln] = UpdateSoln(soln, t_soln, t_event, y_soln, y_event)
% Appends solution from most recent ODE integration to complete history
soln.t    = [soln.t;    t_soln(2:end);   t_event];
soln.Q_m  = [soln.Q_m;  y_soln(2:end,1); y_event(1)];
soln.v_lv = [soln.v_lv; y_soln(2:end,2); y_event(2)];
soln.Q_a  = [soln.Q_a;  y_soln(2:end,3); y_event(3)];
soln.v_d  = [soln.v_d;  y_soln(2:end,4); y_event(4)];

function [final] = FinalSoln(final, t_soln, t_event, y_soln, y_event)
% Saves off solution from ODE integrations of last cycle
final.t    = [final.t;    t_soln; t_event];
final.Q_m  = [final.Q_m;  y_soln(:,1); y_event(1)];
final.v_lv = [final.v_lv; y_soln(:,2); y_event(2)];
final.Q_a  = [final.Q_a;  y_soln(:,3); y_event(3)];
final.v_d  = [final.v_d;  y_soln(:,4); y_event(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following functions test for valve opening and closing conditions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, isterminal, direction] = aortic_open(t, y)
global v_lv0 v_d0 C
p_aortic = 1/C * (y(4)' - v_d0); % (note Q_a = 0)
p_ventricle = elastance(t) .* (y(2)' - v_lv0);
value = p_aortic - p_ventricle;
isterminal = 1;   % stop the integration
direction = -1;   % consider only zero where the event function decreases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, isterminal, direction] = aortic_close(~, y)
% check when flow goes to zero
value = y(3);     % detect height = 0
isterminal = 1;   % stop the integration
direction = -1;   % % consider only zero where the event function decreases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, isterminal, direction] = mitral_open(t, y)
global p_la v_lv0 
% Check for zero crossing of pressure difference between LV and LA
p_lv = elastance(t) .* (y(2) - v_lv0);
value = p_lv - p_la; % detect height = 0
isterminal = 1;   % stop the integration
direction = -1;   % % consider only zero where the event function decreases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, isterminal, direction] = mitral_close(~, y)
% check when flow goes to zero
value = y(1);     % detect height = 0
isterminal = 1;   % stop the integration
direction = -1;   % consider only zero where the event function decreases


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following functions evaluate the derivatives of the state variables.

% dQ_m/dt=1/L_av*(p_la-E_lv*(v_lv-v_lv0))
% dv_lv/dt=Q_m-Q_a
% dQ_a/dt=1/L_va*(E_lv*(v_lv-v_lv0)-1/C*(v_d-v_d0)-Q_a*R_c)
% dv_d/dt=Q_a-(v_d-v_d0)/(R_d*C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = dydt_isovolumic(~,y)
global v_d0 R_d C
dydt=zeros(4,1);
dydt(1)=0;
dydt(2)=0;
dydt(3)=0;
dydt(4)=-(y(4)-v_d0)/(R_d*C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = dydt_ejection(t,y)
global v_lv0 L_va v_d0 R_c R_d C R_VA
dydt=zeros(4,1);
dydt(1)=0;
dydt(2)=-y(3);
dydt(3)=1/L_va*(elastance(t)*(y(2)-v_lv0)-R_VA*y(3)-1/C*(y(4)-v_d0)-y(3)*R_c);
dydt(4)=y(3)-(y(4)-v_d0)/(R_d*C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = dydt_filling(t,y)
global p_la L_av v_lv0 v_d0 R_d C R_AV
dydt=zeros(4,1);
dydt(1)=1/L_av*(p_la-y(1)*R_AV-elastance(t)*(y(2)-v_lv0));
dydt(2)=y(1);
dydt(3)=0;
dydt(4)=-(y(4)-v_d0)/(R_d*C);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [E] = elastance(t)
% Return the elastance value from the normal data in the following paper:
% Single-Beat Estimation of End-Systolic Pressure-Volume Relation in
% Humans, A New Method With the Potential for Noninvasive Application
% Authors: Hideaki Senzaki, MD; Chen-Huan Chen, MD; David A. Kass, MD
% Circulation 1996; 94:2497-2506

% pi/2 was added to this table to cause sin()=1 for the first term.
Senzaki_table_4=[28.38975 pi/2;
    37.58583 .08367674;
    21.02345 -1.486758;
    7.665592 2.865675;
    4.809436 .1677238;
    4.181973 4.630239;
    1.940692 3.088379;
    .5870049 -.3053668;
    1.181256 4.410703;
    .84039 3.181538;
    .02259011 1.242886;
    .3071458 4.156753;
    .3226207 2.946186];
hrt_rate=95; % Table 1 states that normal subjects had heart rates of 95 beats per minute.
E_es=2.0*1333; % Table 1 states that normal subjects had end-systolic elastances of 2.0 mmHg/ml.

period=60/hrt_rate;
t_sin = t/period*2*pi;

E = E_es * sum(10^-2 * ones(size(t)) * Senzaki_table_4(:,1)' .* sin(t_sin*(0:12) + ones(size(t)) * Senzaki_table_4(:,2)'),2);
