clear all; clc; 

nx = 2; ny = 1; 
R1 = 0.02; X1 = 0.0243; Xtr = 0.00557;
KQi = 0.1; KVi = 40; R = R1; X = X1+Xtr; E = 1.0164;
Param = [KQi KVi R X E];

%% EKF
zo = 1*ones(3,1);
ns = 1000;
N = 2*(nx/ny) - 1;
options_steady_states=optimset('Display','iter','maxfunevals', 100000000,'maxiter', 1000); % for steady state
g= @(x0) steady_states(x0,Param); 
[x0_true,~] = fsolve(g,zo,options_steady_states); % finding equilibrium
x0_true = [0.5,0.75,1.02142203324943];
dt_meas = 1/ns;   % Measurement sampling interval (s)
t_meas = 0:dt_meas:1; % Measurement time points
N_meas = length(t_meas); % Number of measurements
dt_ode = dt_meas; 


% Filter's initial state estimate
x0_est =  [x0_true(1);x0_true(2);x0_true(3)];

% Filter's initial covariance
len_filter = length(x0_est);
P0 = 4*eye(len_filter);
Q_scale = 0.25;
Q = Q_scale*eye(3);
Q_filter = Q_scale*eye(len_filter);
R = 0.01;
% Sv = L_Y_theta;



x_true = zeros(3, N_meas);
x_true(:, 1) = x0_true;

% Intermediate time points for ODE solver within each measurement interval
t_interp_steps = (t_meas(2) - t_meas(1)) / dt_ode; % Number of dt_ode steps per dt_meas
if abs(t_interp_steps - round(t_interp_steps)) > 1e-9
    error('dt_meas must be a multiple of dt_ode for consistent simulation steps.');
end
t_interp_steps = round(t_interp_steps); % Ensure integer

current_true_state = x0_true; % State of the true system that evolves continuously

M_DAE2=eye(3);
M_DAE2(3,3)=0;
% M_DAE2(10:12,10:12)=0;
optODE=odeset('Mass',M_DAE2,'RelTol',1e-12);
M_DAE3=eye(12);
M_DAE3(3,3)=0;
M_DAE3(10:12,10:12)=0;
optODE3=odeset('Mass',M_DAE3,'RelTol',1e-12);
rng(10,"threefry"); 
random_process = 5*randn(1)*[3;2;1];
rng(10,"threefry"); 
random_measurement = randn(1, N_meas);
for i = 1:(N_meas - 1)
    t_start_interval = t_meas(i);
    t_end_interval = t_meas(i+1);
    
    % Simulate true system over each small dt_ode step within the interval
    for j = 1:t_interp_steps %This nested loop is to have higher accuracy for the added process noise
        t_current = t_start_interval + (j-1)*dt_ode;
        t_next = t_start_interval + j*dt_ode;
        
        M = [[1 0 0]',eye(3)];
        % Integrate current_true_state over dt_ode
        [~,x_dae_sol] = ode15s(@(t,x_dae_sol) true_dynamics_no_sens(t,x_dae_sol,Param,M),[t_current, t_next],current_true_state,optODE);
        current_true_state = x_dae_sol(end, :)';
        
        % Add continuous process noise at each small dt_ode step
        current_true_state(1:3) = current_true_state(1:3) + sqrt(Q * dt_ode) * random_process;
    end
    
    x_true(:, i+1) = current_true_state; % Store the true state at the measurement point
end

% Generate noisy measurements only at the measurement time points
y_meas = x_true(2,:).*x_true(3,:) + sqrt(R) * random_measurement;

x_est = zeros(len_filter, N_meas);
P_est = zeros(len_filter, len_filter, N_meas);

x_est(:, 1) = x0_est;
P_est(:, :, 1) = P0;
L_k(:, 1) = zeros(len_filter,1);
y_innov_all = zeros(N_meas,1);

% Store predicted state and covariance for plotting
x_pred_plot = zeros(3*4, N_meas);
x_pred_plot(1:3, 1) = x0_est(1:3); 
x_pred_plot(4:end,1) = reshape(eye(3),9,1);
y_pred_plot = zeros(1, N_meas);
y_pred_plot2 = zeros(1, N_meas);
y_pred_plot(1) = x0_true(2);
y_pred_plot2(1) = x0_true(2);
P_pred_plot = zeros(len_filter, len_filter, N_meas);
P_pred_plot(:,:,1) = P0;
C_k_all = zeros(N_meas,3);
Sx2_k_all = zeros(N_meas,3);
Sx2_kplus1_all = zeros(N_meas,3);
Sxy_k_all = zeros(N_meas,3);
Sxy_kplus1_all = zeros(N_meas,3);

x2_est_2 = zeros(1, N_meas); 
x3_est_2 = zeros(1, N_meas);

for k = 1:(N_meas - 1)
    % --- Prediction Step (Continuous-time integration) ---
    t_start = t_meas(k);
    t_end = t_meas(k+1);
    x_k_k = x_pred_plot(:, k);
    x_k_k(1:3) = x_est(1:3, k);
    P_k_k = P_est(:, :, k);
    x_kplus1_k = x_est(:, k);

    M_DAE4=eye(21);
    M_DAE4(3,3)=0; M_DAE4(10:12,10:12)=0;
    optODE4=odeset('Mass',M_DAE4,'RelTol',1e-12);
    zo = x_k_k(3);
    g = @(w0) algebraic_constraint(w0,[Param,x_k_k(2)]); 
    options_algebraic_constraint_solve = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-16, 'TolX', 1e-16);
    [V_consistent,~] = fsolve(g,zo,options_algebraic_constraint_solve);
    x_k_k(3) = V_consistent;
    A = observer_A_function(x_kplus1_k,Param);
    Sx3_k = [-1/A(3,3) * [A(3,1) A(3,2)],0];
    x_k_k(4:end) = reshape(eye(3),9,1);
    x_k_k(10:end) = Sx3_k';
    z0_p0 = [x_k_k;P_k_k(:)];
    [~,Z_P_dae_sol] = ode15s(@(t,z0_p0) observer_dynamics_w_covariances(t,z0_p0,Param,Q_filter,len_filter),[t_start, t_end],z0_p0,optODE4);
    xf_dae_est = Z_P_dae_sol(end,1:4*len_filter);
    P_ode_sol = Z_P_dae_sol(end,4*len_filter+1:end);
    P_kplus1_k = reshape(P_ode_sol(end,:),len_filter,len_filter);  
    
    

    x0_dae_pred = x_pred_plot(:, k);
    x0_dae_pred(4:end) = reshape(eye(3),9,1);
    [~,x_dae_sol_pred] = ode15s(@(t,x_dae_sol) observer_dynamics(t,x_dae_sol,Param,M),[t_start, t_end],x0_dae_pred,optODE3);
    xf_dae_pred = x_dae_sol_pred(end,1:3)';
    
    Sx1_k = x_dae_sol_pred(1,4:6); Sx2_k = x_dae_sol_pred(1,7:9); Sx3_k = x_dae_sol_pred(1,10:12);
    Sx1_kplus1 = x_dae_sol_pred(end,4:6); Sx2_kplus1 = x_dae_sol_pred(end,7:9); Sx3_kplus1 = x_dae_sol_pred(end,10:12);
    
    x_kplus1_k(1:3) = xf_dae_est(end,1:3)';

    x_pred_plot(:, k+1) = x_dae_sol_pred(end,:)';
    P_pred_plot(:, :, k+1) = P_kplus1_k;
    
    

    % End of Prediction step

    % Begin Measurements update

    
    % Get the current measurement
    y_kplus1 = y_meas(k+1);
    % Calculate Innovation 
    y_innov = y_kplus1 - x_kplus1_k(2).*x_kplus1_k(3);
    % Calculate Jacobian-Like of the Measurement function
    Sy1 = Sx2_k*x_kplus1_k(3) + Sx3_k*x_kplus1_k(2);
    second_term = [x_dae_sol_pred(end,4:6)-x_dae_sol_pred(1,4:6);x_dae_sol_pred(end,7:9)-x_dae_sol_pred(1,7:9);x_dae_sol_pred(end,10:12)-x_dae_sol_pred(1,10:12)];
    Sy = Sy1 + Sy1 * second_term;
    Sxy_k_all(k,:) = Sy; 
    C_kplus1 = Sy;

    S_kplus1 = C_kplus1 * P_kplus1_k * C_kplus1' + R ;
    
    % Calculate the observer gain
    L_kplus1 = P_kplus1_k * C_kplus1' / S_kplus1;
    % Using the LSERC test, zero out the appropriate L entries
    [Diff_States_lno,Diff_States_lo,Alg_States_lno,Alg_States_lo] = LSERC(x_kplus1_k(1:3),t_start,t_end,N);
    for i = 1:length(Diff_States_lno)
        L_kplus1(i) = 0;
    end
    for j = 1:length(Alg_States_lno)
        L_kplus1(j) = 0;
    end
    % Update State Estimate    
    x_kplus1_k(1:3) = x_kplus1_k(1:3) + L_kplus1 * y_innov;
    zo = x_kplus1_k(3);
    g = @(w0) algebraic_constraint(w0,[Param,x_kplus1_k(2)]); 
    options_algebraic_constraint_solve = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-16, 'TolX', 1e-16);
    [V_consistent,~] = fsolve(g,zo,options_algebraic_constraint_solve);
    x_kplus1_k(3) = V_consistent;
   
    x_kplus1_kplus1 = x_kplus1_k;
    
    % Update Covariance Estimate
    P_kplus1_kplus1 = (eye(len_filter) - L_kplus1 * C_kplus1) * P_kplus1_k...
         * (eye(len_filter) - L_kplus1 * C_kplus1)' + L_kplus1 * R * L_kplus1';
          
    % Store updated values
    x_est(:, k+1) = x_kplus1_kplus1;
    P_est(:, :, k+1) = P_kplus1_kplus1;
    L_k(:, k+1) = L_kplus1;
    y_innov_all(k+1) = y_innov;
end
%% Plot Results
close all; clc
font_size = 50;
set(0,'DefaultFigureWindowStyle','docked');
fig = figure('DefaultAxesFontSize',font_size,'defaultLineLineWidth',7,'DefaultAxesTitleFontWeight','bold');
% Plot State 1
plot(t_meas, x_true(1,:),  '*', 'color', '#80B3FF',  'LineWidth', 3, 'MarkerSize', 10); hold on;
plot(t_meas, x_est(1,:), 'k-', 'LineWidth', 3); hold on;
plot(t_meas, x_pred_plot(1,:), '-','color','#77AC30', 'LineWidth', 3); 

legend('True', 'Estimated','Predicted', 'Location', 'best','Fontsize',30);
ylabel('$V_{ref}$','interpreter','latex');
grid on;
% Plot State 2
font_size = 50;
fig = figure('DefaultAxesFontSize',font_size,'defaultLineLineWidth',7,'DefaultAxesTitleFontWeight','bold');
plot(t_meas, x_true(2,:),  '*', 'color', '#80B3FF',  'LineWidth', 3, 'MarkerSize', 10); hold on;
plot(t_meas, x_est(2,:), 'k-', 'LineWidth', 3);
plot(t_meas, x_pred_plot(2,:), '-','color','#77AC30', 'LineWidth', 3); 

legend('True', 'Estimated','Predicted', 'Location', 'best','Fontsize',30);
ylabel('$E''''_{q}$','interpreter','latex');
grid on;
% Plot State 3
font_size = 50;
fig = figure('DefaultAxesFontSize',font_size,'defaultLineLineWidth',7,'DefaultAxesTitleFontWeight','bold');
plot(t_meas, x_true(3,:), '*', 'color', '#80B3FF',  'LineWidth', 3, 'MarkerSize', 10); hold on;
plot(t_meas, x_est(3,:), 'k-', 'LineWidth', 3); hold on;
plot(t_meas, x_pred_plot(3,:), '-','color','#77AC30', 'LineWidth', 3); 

legend('True', 'Estimated','Predicted', 'Location', 'best','Fontsize',30);
ylabel('$V$','interpreter','latex');
grid on;

% Plot measurement
font_size = 50;
fig = figure('DefaultAxesFontSize',font_size,'defaultLineLineWidth',7,'DefaultAxesTitleFontWeight','bold');
plot(t_meas, x_true(2,:).*x_true(3,:),  'd', 'color', '#80B3FF',  'LineWidth', 3, 'MarkerSize', 10); hold on;
plot(t_meas, y_meas, 'r*', 'MarkerSize', 2); hold on;
plot(t_meas, x_est(2,:).*x_est(3,:), 'k-', 'LineWidth', 3); hold on;
plot(t_meas, x_pred_plot(2,:).*x_pred_plot(3,:), '-','color','#77AC30', 'LineWidth', 3); 

legend('True', 'Measurements', 'Estimated','Predicted', 'Location', 'best','Fontsize',30);
xlabel('Time (s)');
ylabel('$y$','interpreter','latex');
grid on;

function [Diff_States_lno,Diff_States_lo,Alg_States_lno,Alg_States_lo] = LSERC(z0,t0,tf,N)
    nx = 2; k_directions = 4;
    R1 = 0.02; X1 = 0.0243; Xtr = 0.00557;
    KQi = 0.1; KVi = 40; R = R1; X = X1+Xtr; E = 1.0164;
    Param = [KQi KVi R X E];
    M_DAE=eye(3+3*k_directions);
    M_DAE(3,3)=0;
    for i =(2*k_directions+3)+1: 3+3*k_directions
        M_DAE(i,i)=0;
    end
    optODE=odeset('Mass',M_DAE,'RelTol',1e-12);
    T_f = [];
    T_P_Matrix_f = [];
    Diff_States = [1 2];
    Diff_States_lno = [];
    Diff_States_lo = [];

    Alg_States = [3];
    Alg_States_lno = [];
    Alg_States_lo = [];

    M = [[1 0 0]',eye(3)];
    % Initial conditions
    % z0 = [0.5,0.75,1.02142203324943];
    Vref0 = z0(1); Eq0 = z0(2); V0 = z0(3);
    SVref0 = [1 0 0]*M;
    SEq0 = [0 1 0]*M;
    A = observer_A_function([Vref0,Eq0,V0],Param);
    SV0 = [-1/A(3,3) * [A(3,1) A(3,2)],0]*M;
    Z0 =[Vref0;Eq0;V0;SVref0(:);SEq0(:);SV0(:)]; %vector of states & state sensitivies to be used in ODE solver
    tspan = linspace(t0,tf,N);
    [Time,Z] = ode15s(@(t,Z) true_dynamics(t,Z,Param,M),tspan,Z0(:),optODE);
    StepCount = length(Time);
    % Calculate the the SVD for the LSERC matrix to check identifiability
    h = Z(:,2).*Z(:,3);
    dhdx2 = Z(:,3);
    dhdw = Z(:,2);
    LD_Y_theta = zeros(StepCount,nx+1);
    step = 1;
    for ti=1:StepCount
        Sx1 = [Z(ti,4) Z(ti,5) Z(ti,6)];
        Sx2 = [Z(ti,8) Z(ti,9) Z(ti,10)];
        Sw = [Z(ti,12) Z(ti,13) Z(ti,14)];
        LD_Y_theta(step,:) = dhdx2(ti).*Sw + dhdw(ti).*Sx2;
        step = step+1;
    end
    M = [[1 0]',eye(2)];
    L_Y_theta = LD_Y_theta/M;
    [U_L,S_L,V_L] = svd(L_Y_theta);
    diag_SL = diag(S_L);
    P_sing_vec = [];
    for i=1:length(diag_SL)
        if diag_SL(i)<1e-4
            P_sing_vec = [P_sing_vec V_L(:,i)];
        end
    end
    [R,p] = rref(P_sing_vec',1e-10);
    C = num2cell( p , 1 );

    if find(Diff_States_lo==1)~= 1
        Diff_States_lno = [Diff_States_lno C];
    end
    Diff_States_lo = setdiff(Diff_States, Diff_States_lno);
    % If all differential states are local observable, then all algebraic states
    % are locally observable
    for s = 1:length(Diff_States_lno) 
        j = Diff_States_lno(s);
        Sw_j = Z(:,(2*k_directions+3)+j);
        Psw_w = Sw_j;
        if rank(Psw_w) ~= 0
            Alg_States_lno = [3];
        end
    end
    Alg_States_lo = setdiff(Alg_States, Alg_States_lno);
end

function dZdt = true_dynamics(t,Z,Param,M)
k_directions = size(M,2);
Vref = Z(1);
Eq = Z(2);
V = Z(3);

KQi = Param(1);
KVi = Param(2);
R = Param(3);
X = Param(4);
E = Param(5);
Xeq = 0.8;

Qcmd = 0.195197282638509; %constant
Q = V*Eq/Xeq - V^2/Xeq;
I = 0.780082515717402;
P = V * I;

dVrefdt = KQi*Qcmd - KQi*Q;
dEqdt = KVi*Vref - KVi*V;
dVdt = V^4 - 2*P*R*V^2 - 2*Q*X*V^2 - E^2*V^2 + R^2*P^2 + X^2*P^2 + R^2*Q^2 + X^2*Q^2;

SVref = Z(3+1:3+k_directions);
SEq = Z(3+k_directions+1:3+2*k_directions);
SV = Z(3+2*k_directions+1:3+3*k_directions);

dQdz = [SV(1)*Eq/Xeq+SEq(1)*V/Xeq-2*V*SV(1)/Xeq; SV(2)*Eq/Xeq+SEq(2)*V/Xeq-2*V*SV(2)/Xeq; SV(3)*Eq/Xeq+SEq(3)*V/Xeq-2*V*SV(3)/Xeq; SV(4)*Eq/Xeq+SEq(4)*V/Xeq-2*V*SV(4)/Xeq];

dSVrefdt = [-KQi*dQdz(1); -KQi*dQdz(2); -KQi*dQdz(3); -KQi*dQdz(4)];
dSEqdt = [KVi*SVref(1)-KVi*SV(1); KVi*SVref(2)-KVi*SV(2); KVi*SVref(3)-KVi*SV(3); KVi*SVref(4)-KVi*SV(4)];
dSVdt = [4*V^3*SV(1)-2*P*R*2*V*SV(1)-4*2*Q*X*V*SV(1)-2*dQdz(1)*X*V^2-2*E^2*V*SV(1)+0+0+2*R^2*Q*dQdz(1)+2*X^2*Q*dQdz(1);...
         4*V^3*SV(2)-2*P*R*2*V*SV(2)-4*2*Q*X*V*SV(2)-2*dQdz(2)*X*V^2-2*E^2*V*SV(2)+0+0+2*R^2*Q*dQdz(2)+2*X^2*Q*dQdz(2);...
         4*V^3*SV(3)-2*P*R*2*V*SV(3)-4*2*Q*X*V*SV(3)-2*dQdz(3)*X*V^2-2*E^2*V*SV(3)+0+0+2*R^2*Q*dQdz(3)+2*X^2*Q*dQdz(3); ...
         4*V^3*SV(4)-2*P*R*2*V*SV(4)-4*2*Q*X*V*SV(4)-2*dQdz(4)*X*V^2-2*E^2*V*SV(4)+0+0+2*R^2*Q*dQdz(4)+2*X^2*Q*dQdz(4)];
     
dZdt = [dVrefdt;dEqdt;dVdt;dSVrefdt;dSEqdt;dSVdt];

end

function dZdt = true_dynamics_no_sens(t,Z,Param,M)

k_directions = size(M,2);
Vref = Z(1);
Eq = Z(2);
V = Z(3);

KQi = Param(1);
KVi = Param(2);
R = Param(3);
X = Param(4);
E = Param(5);
Xeq = 0.8;

Qcmd = 0.195197282638509; %constant
Q = V*Eq/Xeq - V^2/Xeq;
I = 0.780082515717402;
P = V * I;

dVrefdt = KQi*Qcmd - KQi*Q;
dEqdt = KVi*Vref - KVi*V;
dVdt = V^4 - 2*P*R*V^2 - 2*Q*X*V^2 - E^2*V^2 + R^2*P^2 + X^2*P^2 + R^2*Q^2 + X^2*Q^2;

dZdt = [dVrefdt;dEqdt;dVdt];

end


function dZdt = observer_dynamics(t,Z,Param,M)
k_directions = 3;
Vref = Z(1);
Eq = Z(2);
V = Z(3);

KQi = Param(1);
KVi = Param(2);
R = Param(3);
X = Param(4);
E = Param(5);
Xeq = 0.8;

Qcmd = 0.195197282638509; %constant
Q = V*Eq/Xeq - V^2/Xeq;
I = 0.780082515717402;
P = V * I;

dVrefdt = KQi*Qcmd - KQi*Q;
dEqdt = KVi*Vref - KVi*V;
dVdt = V^4 - 2*P*R*V^2 - 2*Q*X*V^2 - E^2*V^2 + R^2*P^2 + X^2*P^2 + R^2*Q^2 + X^2*Q^2;

SVref = Z(3+1:3+k_directions);
SEq = Z(3+k_directions+1:3+2*k_directions);
SV = Z(3+2*k_directions+1:3+3*k_directions);

dQdz = [SV(1)*Eq/Xeq+SEq(1)*V/Xeq-2*V*SV(1)/Xeq; SV(2)*Eq/Xeq+SEq(2)*V/Xeq-2*V*SV(2)/Xeq; SV(3)*Eq/Xeq+SEq(3)*V/Xeq-2*V*SV(3)/Xeq];

dSVrefdt = [-KQi*dQdz(1); -KQi*dQdz(2); -KQi*dQdz(3)];
dSEqdt = [KVi*SVref(1)-KVi*SV(1); KVi*SVref(2)-KVi*SV(2); KVi*SVref(3)-KVi*SV(3)];
dSVdt = [4*V^3*SV(1)-2*P*R*2*V*SV(1)-4*2*Q*X*V*SV(1)-2*dQdz(1)*X*V^2-2*E^2*V*SV(1)+0+0+2*R^2*Q*dQdz(1)+2*X^2*Q*dQdz(1);...
         4*V^3*SV(2)-2*P*R*2*V*SV(2)-4*2*Q*X*V*SV(2)-2*dQdz(2)*X*V^2-2*E^2*V*SV(2)+0+0+2*R^2*Q*dQdz(2)+2*X^2*Q*dQdz(2); ...
         4*V^3*SV(3)-2*P*R*2*V*SV(3)-4*2*Q*X*V*SV(3)-2*dQdz(3)*X*V^2-2*E^2*V*SV(3)+0+0+2*R^2*Q*dQdz(3)+2*X^2*Q*dQdz(3)];

dZdt = [dVrefdt;dEqdt;dVdt;dSVrefdt;dSEqdt;dSVdt];
end


function dZdt = observer_dynamics_w_covariances(t,Z,Param,Q_filter,len_filter)
Vref = Z(1);
Eq = Z(2);
V = Z(3);

KQi = Param(1);
KVi = Param(2);
R = Param(3);
X = Param(4);
E = Param(5);
Xeq = 0.8;

Qcmd = 0.195197282638509; %constant
Q = V*Eq/Xeq - V^2/Xeq;
I = 0.780082515717402;
P = V * I;

dVrefdt = KQi*Qcmd - KQi*Q;
dEqdt = KVi*Vref - KVi*V;
dVdt = V^4 - 2*P*R*V^2 - 2*Q*X*V^2 - E^2*V^2 + R^2*P^2 + X^2*P^2 + R^2*Q^2 + X^2*Q^2;

SVref = Z(4:6);
SEq = Z(7:9);
SV = Z(10:12);

dQdz = [SV(1)*Eq/Xeq+SEq(1)*V/Xeq-2*V*SV(1)/Xeq; SV(2)*Eq/Xeq+SEq(2)*V/Xeq-2*V*SV(2)/Xeq; SV(3)*Eq/Xeq+SEq(3)*V/Xeq-2*V*SV(3)/Xeq];

dSVrefdt = [-KQi*dQdz(1); -KQi*dQdz(2); -KQi*dQdz(3)];
dSEqdt = [KVi*SVref(1)-KVi*SV(1); KVi*SVref(2)-KVi*SV(2); KVi*SVref(3)-KVi*SV(3)];
dSVdt = [4*V^3*SV(1)-2*P*R*2*V*SV(1)-4*2*Q*X*V*SV(1)-2*dQdz(1)*X*V^2-2*E^2*V*SV(1)+0+0+2*R^2*Q*dQdz(1)+2*X^2*Q*dQdz(1);...
         4*V^3*SV(2)-2*P*R*2*V*SV(2)-4*2*Q*X*V*SV(2)-2*dQdz(2)*X*V^2-2*E^2*V*SV(2)+0+0+2*R^2*Q*dQdz(2)+2*X^2*Q*dQdz(2); ...
         4*V^3*SV(3)-2*P*R*2*V*SV(3)-4*2*Q*X*V*SV(3)-2*dQdz(3)*X*V^2-2*E^2*V*SV(3)+0+0+2*R^2*Q*dQdz(3)+2*X^2*Q*dQdz(3)];

P_k_k = Z(13:21);
P_k_k_matrix = reshape(P_k_k,len_filter,len_filter);

A = [dSVrefdt';
     dSEqdt';
     dSVdt'];
dPdt_matrix =  A * P_k_k_matrix + P_k_k_matrix * A' + Q_filter;
dPdt = dPdt_matrix(:);   

dZdt = [dVrefdt;dEqdt;dVdt;dSVrefdt;dSEqdt;dSVdt;dPdt];
end

function f = steady_states(x_guess,Param)

Vref = x_guess(1);
Eq = x_guess(2);
V = x_guess(3);

KQi = Param(1);
KVi = Param(2);
R = Param(3);
X = Param(4);
E = Param(5);
Xeq = 0.8;

Qcmd = 0.195197282638509; %constant
Q = V*Eq/Xeq - V^2/Xeq;
I = 0.780082515717402;
P = V * I;
PFE = 10;
Qcmd = tan(PFE)*P;
       
f=[ ...
    KQi*Qcmd - KQi*Q;  
    KVi*Vref - KVi*V;
    V^4 - 2*P*R*V^2 - 2*Q*X*V^2 - E^2*V^2 + R^2*P^2 + X^2*P^2 + R^2*Q^2 + X^2*Q^2];

end

function val = algebraic_constraint(V, Param)
    KQi = Param(1);
    KVi = Param(2);
    R = Param(3);
    X = Param(4);
    E = Param(5);
    Eq = Param(6);
    
    Xeq = 0.8;
    Qcmd = 0.195197282638509; %constant
    Q = V*Eq/Xeq - V^2/Xeq;
    I = 0.780082515717402;
    P = V * I;
    PFE = 10;
    Qcmd = tan(PFE)*P;

    val = V.^4 - (2*(P*R + Q*X) + E.^2).*V.^2 + (R.^2 + X.^2).*(P.^2 + Q.^2);
end

function A = observer_A_function(Z,Param)
Vref = Z(1);
Eq = Z(2);
V = Z(3);

KQi = Param(1);
KVi = Param(2);
R = Param(3);
X = Param(4);
E = Param(5);
Xeq = 0.8;

Qcmd = 0.195197282638509; %constant
Q = V*Eq/Xeq - V^2/Xeq;
I = 0.780082515717402;
P = V * I;

dQdz = [0, V/Xeq, Eq/Xeq-2*V/Xeq];

dSVrefdt = [-KQi*dQdz(1), -KQi*dQdz(2), -KQi*dQdz(3)];
dSEqdt = [KVi, 0, -KVi];
dSVdt = [-2*dQdz(1)*X*V^2+2*R^2*Q*dQdz(1)+2*X^2*Q*dQdz(1),...
         -2*dQdz(2)*X*V^2+2*R^2*Q*dQdz(2)+2*X^2*Q*dQdz(2), ...
         4*V^3-2*P*R*2*V-2*Q*X*2*V-2*E^2*V+R^2*Q*dQdz(3)+2*X^2*Q*dQdz(3)];
     
A = [dSVrefdt;dSEqdt;dSVdt];
end