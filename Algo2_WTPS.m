clear all;clc;close all;
np = 5; k_directions = 6; threshold_percent_max = 0.05;
%Use LHS to construct the reference parameters set
R1 = 0.02; X1 = 0.0243; Xtr = 0.00557;
KQi = 0.1; KVi = 40; R = R1; X = X1+Xtr; E = 1.0164;
P0_true = [KQi KVi R X E];

P_set = {P0_true,(1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true,(1+0.1*rand)*P0_true,...
         P0_true,(1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true,(1+0.1*rand)*P0_true,...
         P0_true,(1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true,(1+0.1*rand)*P0_true,...
         P0_true,(1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true, (1+0.1*rand)*P0_true,(1+0.1*rand)*P0_true,(1+0.1*rand)*P0_true};
nonidentifiable_param_index = [];

M_DAE=eye(3+3*k_directions);
M_DAE(3,3)=0;
for i =(2*k_directions+3)+1: 3+3*k_directions
    M_DAE(i,i)=0;
end
optODE=odeset('Mass',M_DAE,'RelTol',1e-8);
T_f=[];
T_P_Matrix_f=[];



for pp = 1:length(P_set)
    P0 = cell2mat(P_set(pp));
    Param = P0;
    zo = 1*ones(3,1);
    options_steady_states=optimset('Display','iter','maxfunevals', 100000000,'maxiter', 1000); % for steady state
    g= @(x0) steady_states(x0,P0); 
    [z0,~] = fsolve(g,zo,options_steady_states); % finding equilibrium

    % Primary probing directions within the directions matrix P
    D_set = {[[1 0 0 0 0]',eye(5)],[[0 1 0 0 0]',eye(5)],[[0 0 1 0 0]',eye(5)],[[0 0 0 1 0]',eye(5)],[[0 0 0 0 1]',eye(5)]...
             [[-1 0 0 0 0]',eye(5)],[[0 -1 0 0 0]',eye(5)],[[0 0 -1 0 0]',eye(5)],[[0 0 0 -1 0]',eye(5)],[[0 0 0 0 -1]',eye(5)]};
    % Value of epsilon for twin probing (epsilon_twin)
    epsilon_twin_probing = 0.01;
    % Directions of perturbations of epsilon_twin
    M_twin_probing_set = {[[0 1 0 0 0]',zeros(5)],[[1 0 0 0 0]',zeros(5)],[[0 0 0 1 0]',zeros(5)],[[0 0 1 0 0]',zeros(5)]}; 

    % Value of epsilon for singular probing(epsilon_sing)
    epsilon_sing_probing = 0.01;
    D_sing=[];
    V_sing=[];
    Theta = [1 2 3 4 5];
    Theta_lni = [];
    Theta_li = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Primary & twin Probing stages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:length(D_set)
       for k_twin_prob = 1:length(M_twin_probing_set)
        M = cell2mat(D_set(k)); %Get current directions matrix
        M_twin_probing = cell2mat(M_twin_probing_set(k_twin_prob)); %Get current epsilon_twin
        if abs(M(:,1)) == M_twin_probing(:,1)
            continue
        end
        M_twin_probing(:,1) = M_twin_probing(:,1).*epsilon_twin_probing;
        if sum(sign(M(:,1))) == 1 %check the sign of the sum of elements of the first column of M
           M_twin_probing = M_twin_probing+M; 
        elseif sum(sign(M(:,1))) == -1
           M_twin_probing = -M_twin_probing+M;
        end    
        %start solving the forward sensitivity system for every direction in
        %primary and twin probing directions matrices
        % Initial conditions
        Vref0 = z0(1); Eq0 = z0(2); V0 = z0(3);
        SVref0 = [0 0 0 0 0]*M;
        SEq0 = [0 0 0 0 0]*M;
        SV0 = [0 0 0 0 0]*M;
        SVref0_twin_probing = [0 0 0 0 0]*M_twin_probing;
        SEq0_twin_probing = [0 0 0 0 0]*M_twin_probing;
        SV0_twin_probing = [0 0 0 0 0]*M_twin_probing;
        Z0 =[Vref0;Eq0;V0;SVref0(:);SEq0(:);SV0(:)]; %vector of states & state sensitivies to be used in ODE solver
        Z0_twin_probing =[Vref0;Eq0;V0;SVref0_twin_probing(:);SEq0_twin_probing(:);SV0_twin_probing(:)];
        % Time interval/steps
        t0 = 0; tf = 1;
        N = 2*np; 
        tspan = t0:1/(2*N):tf;
        [Time,Z] = ode15s(@(t,Z) dynamics(t,Z,P0,M),tspan,Z0(:),optODE);
        [T_twin_probing,Z_twin_probing] = ode15s(@(t_twin_probing,Z_twin_probing) dynamics(t_twin_probing,Z_twin_probing,P0,M_twin_probing),tspan,Z0_twin_probing(:),optODE);
        StepCount = length(Time);
        % Calculate the the SVD for the LSERC matrix to check identifiability
        h = Z(:,3); %the output is the voltage
        [m,dummy1]=size(h);
        dhdV = 1;
        Yp = zeros(StepCount,np);
        LD_Y_theta = zeros(StepCount,np+1);
        LD_Y_theta_twin_probing = zeros(StepCount,np+1);
        t=tspan(StepCount);
        step = 1;
        p1 = P0(1);
        p2 = P0(2);
        for ti=1:StepCount
            %We need to calculate relative sensitivities  
            xp = [Z(ti,16) Z(ti,17) Z(ti,18) Z(ti,19) Z(ti,20) Z(ti,21)];
            xp_twin_probing = [Z_twin_probing(ti,16) Z_twin_probing(ti,17) Z_twin_probing(ti,18) Z_twin_probing(ti,19) Z_twin_probing(ti,20) Z_twin_probing(ti,21)];
            dhdp = xp;
            dhdp_twin_probing = xp_twin_probing;
            LD_Y_theta(step,:) = dhdV * xp ;
            LD_Y_theta_twin_probing(step,:) = dhdV * xp_twin_probing ;
            step = step+1;
        end
        states_solution = [Z(:,1);Z(:,2)];
        [U_LD,S_LD,V_LD] = svd(LD_Y_theta);
        [U_LD_twin_probing,S_LD_twin_probing,V_LD_twin_probing] = svd(LD_Y_theta_twin_probing);

        L_Y_theta = LD_Y_theta/M;
        L_Y_theta_twin_probing = LD_Y_theta_twin_probing/M_twin_probing;
        [U_L,S_L,V_L] = svd(L_Y_theta);
        diag_SL = diag(S_L);
        P_sing_vec = [];
        for i=1:length(diag_SL)
            if diag_SL(i)<1e-14
                P_sing_vec = [P_sing_vec V_L(:,i)];
            end
        end
        [R,p] = rref(P_sing_vec',1e-10);
        C = num2cell( p , 1 );
        if k==1&&k_twin_prob==1
            Theta_lni = p;
        else
        for i=1:length(C)
            if find(Theta_li==1)~= 1
                Theta_lni = [Theta_lni C];
            end
        end
        end
        [U_L_twin_probing,S_L_twin_probing,V_L_twin_probing] = svd(L_Y_theta_twin_probing);
        %Collect results in tables for easier access
        T.Param = array2table(P0);
        T.LD_Y_theta = array2table(LD_Y_theta);
        T_P_Matrix = array2table(M);
        T.LD_Y_theta_twin_probing = LD_Y_theta_twin_probing;
        T.L_Y_theta = L_Y_theta;
        T.L_Y_theta_twin_probing = L_Y_theta_twin_probing;
        T.S_LD = S_LD;
        T.S_LD_twin_probing = S_LD_twin_probing;
        T.S_L = S_L;
        T.R = R;
        T.p = p;
        T.normalized_sigma = diag(S_L);
        T.V_L = V_L;
        T.S_L_twin_probing = S_L_twin_probing;
        T.time = Time;
        T.states_solution = states_solution;
        T.algebraic = Z(:,3);
        T_f=[T_f;T];
        T_P_Matrix_f = [T_P_Matrix_f;T_P_Matrix];
       end
    end
end
%% 
   clear all;
   np = 5; k_directions = 6;
   R1 = 0.02; X1 = 0.0243; Xtr = 0.00557; %X1+Xtr=0.029870000000000
   KQi = 0.1; KVi = 40; R = R1; X = X1+Xtr; E = 1.0164;
   P0 = [KQi KVi R X E];
   Param_reduced = [P0(3)+0.01*P0(3) P0(4)+0.01*P0(4)];
%    Param_reduced = [P0(3) P0(4)];
   M = [[1 0 0 0 0]',eye(5)];
   M_DAE=eye(3+3*k_directions);
   M_DAE(3,3)=0;
   for i =(2*k_directions+3)+1: 3+3*k_directions
       M_DAE(i,i)=0;
   end
   optODE=odeset('Mass',M_DAE,'RelTol',1e-8);
   zo = 1*ones(3,1);
   options_steady_states=optimset('Display','iter','maxfunevals', 100000000,'maxiter', 1000); % for steady state
   g= @(x0) steady_states(x0,P0); 
   [z0,~] = fsolve(g,zo,options_steady_states); % finding equilibrium
   Vref0 = z0(1); Eq0 = z0(2); V0 = z0(3);
   SVref0 = [0 0 0 0 0]*M;
   SEq0 = [0 0 0 0 0]*M;
   SV0 = [0 0 0 0 0]*M;
   Z0 =[Vref0;Eq0;V0;SVref0(:);SEq0(:);SV0(:)];
   
   t0 = 0; tf = 1;
   N = 1*np; 
   tspan = t0:1/N:tf;
   [Time,Z_true] = ode15s(@(t,Z) dynamics(t,Z,P0,M),tspan,Z0(:),optODE);
   x_true = Z_true(:,1:3);
   
   output_samples = x_true;
   
   p_est = param_est_oc(Param_reduced,P0,tspan,output_samples);
      
   k_directions = 2;
   M_DAE=eye(3+3*k_directions);
   M_DAE(3,3)=0;
   for i =(3+2*k_directions)+1: (3+3*k_directions)
        M_DAE(i,i)=0;
   end
   optODE=odeset('Mass',M_DAE,'RelTol',1e-8);
   zo = 1*ones(3,1);
   options_steady_states=optimset('Display','iter','maxfunevals', 100000000,'maxiter', 1000); % for steady state
   P_steady_states = [P0(1) P0(2) p_est(1) p_est(2) P0(5)];
   g= @(x0) steady_states(x0,P_steady_states); 
   [z0,~] = fsolve(g,zo,options_steady_states); % finding equilibrium
   M = eye(2);
   Vref0 = z0(1); Eq0 = z0(2); V0 = z0(3);
   SVref0 = [0 0]*M;
   SEq0 = [0 0]*M;
   SV0 = [0 0]*M;
   Z0 =[Vref0;Eq0;V0;SVref0(:);SEq0(:);SV0(:)];
   t0 = 0; tf = 1;
   N = 2*np; 
   tspan = t0:1/(2*N):tf;
   [Time,Z_est] = ode15s(@(t,Z) dynamics_reduced(t,Z,p_est,P0),tspan,Z0(:),optODE);
   
   x_est = Z_est(:,1:3);
   
   %% 
   close all;
   set(0,'DefaultFigureWindowStyle','docked');font_size = 40;
   fig = figure('DefaultAxesFontSize',font_size,'defaultLineLineWidth',5,'DefaultAxesTitleFontWeight','bold');
   plot(tspan,x_true,'r-');
   hold on;
   plot(tspan,x_est,'k--');

   hold on;
   plot(tspan,x_p0,'b');
   hold on;
   plot(tspan,x_avg,'Color',[0.4660 0.6740 0.1880]);
   
   xlabel('$t$','Interpreter','latex')
   ylabel('$y(t; \theta_r)$','Interpreter','latex')
   legend('$y(t;\theta_{true})$','$y(t;\theta_{est})$','$y(t;\theta_{initial})$','$y(t;\theta_{avg})$','Interpreter','latex','Location','northwest')
   %% 
function dZdt = dynamics(t,Z,Param,M)
k_directions = size(M,2);
Vref = Z(1);
Eq = Z(2);
V = Z(3);

KQi = Param(1);
KVi = Param(2);
R = Param(3);
X = Param(4);
E = Param(5);
% Xeq = Param(6);
Xeq = 0.8;

% Qcmd = 0.195197282638509; %constant
Q = V*Eq/Xeq - V^2/Xeq;
I = 0.948572674561663;
P = V * I;
PFE = 10;
Qcmd = tan(PFE)*P;

dVrefdt = KQi*Qcmd - KQi*Q;
dEqdt = KVi*Vref - KVi*V;
dVdt = V^4 - 2*P*R*V^2 - 2*Q*X*V^2 - E^2*V^2 + R^2*P^2 + X^2*P^2 + R^2*Q^2 + X^2*Q^2;

SVref = Z(3+1:3+k_directions);
SEq = Z(3+k_directions+1:3+2*k_directions);
SV = Z(3+2*k_directions+1:3+3*k_directions);

dQdParam = [SV(1)*Eq/Xeq+V*SEq(1)/Xeq-2*SV(1)*V/Xeq;SV(2)*Eq/Xeq+V*SEq(2)/Xeq-2*SV(2)*V/Xeq;SV(3)*Eq/Xeq+V*SEq(3)-2*SV(3)*V/Xeq/Xeq;SV(4)*Eq/Xeq+V*SEq(4)/Xeq-2*SV(4)*V/Xeq;SV(5)*Eq/Xeq+V*SEq(5)/Xeq-2*SV(5)*V/Xeq;SV(6)*Eq/Xeq+V*SEq(6)/Xeq-2*SV(6)*V/Xeq];

dSVrefdt = [Qcmd-Q-KQi*dQdParam(1);Qcmd-Q-KQi*dQdParam(2);-KQi*dQdParam(3);-KQi*dQdParam(4);-KQi*dQdParam(5);-KQi*dQdParam(6)];
dSEqdt = [KVi*SVref(1)-KVi*SV(1);KVi*SVref(2)-KVi*SV(2);Vref+KVi*SVref(3)-V-KVi*SV(3);KVi*SVref(4)-KVi*SV(4);KVi*SVref(5)-KVi*SV(5);KVi*SVref(6)-KVi*SV(6)];
dSVdt = [4*V^3*SV(1)-2*P*R*2*V*SV(1)-4*2*Q*X*V*SV(1)-2*dQdParam(1)*X*V^2-2*E^2*V*SV(1)+0+0+2*R^2*Q*dQdParam(1)+2*X^2*Q*dQdParam(1);...
         4*V^3*SV(2)-2*P*R*2*V*SV(2)-4*2*Q*X*V*SV(2)-2*dQdParam(2)*X*V^2-2*E^2*V*SV(2)+0+0+2*R^2*Q*dQdParam(2)+2*X^2*Q*dQdParam(2); ...
         4*V^3*SV(3)-2*P*R*2*V*SV(3)-4*2*Q*X*V*SV(3)-2*dQdParam(3)*X*V^2-2*E^2*V*SV(3)+0+0+2*R^2*Q*dQdParam(3)+2*X^2*Q*dQdParam(3);...
         4*V^3*SV(4)-2*P*V^2-2*P*R*2*V*SV(4)-4*2*Q*X*V*SV(4)-2*dQdParam(4)*X*V^2-2*E^2*V*SV(4)+2*R*P^2+0+2*R*Q^2+2*R^2*Q*dQdParam(4)+2*X^2*Q*dQdParam(4);...
         4*V^3*SV(5)-2*P*R*2*V*SV(5)-2*Q*V^2-4*2*Q*X*V*SV(5)-2*dQdParam(5)*X*V^2-2*E^2*V*SV(5)+0+2*X*P^2+2*R^2*Q*dQdParam(5)+2*X*Q^2+2*X^2*Q*dQdParam(5);...
         4*V^3*SV(6)-2*P*R*2*V*SV(6)-4*2*Q*X*V*SV(6)-2*dQdParam(6)*X*V^2-2*E*V^2-2*E^2*V*SV(6)+0+0+2*R^2*Q*dQdParam(6)+2*X^2*Q*dQdParam(6)];
if t>0.8
    break_here = 'X';
end
dZdt = [dVrefdt;dEqdt;dVdt;dSVrefdt;dSEqdt;dSVdt];
end

function p_est = param_est_oc(Param_reduced,P0,tspan,output_samples)
   Algorithm = 'sqp';
   optNLP = optimset( 'Algorithm',Algorithm,'Hessian','bfgs','LargeScale', 'off', 'GradObj', 'on', 'GradConstr', 'off',...
        'DerivativeCheck', 'off', 'Display', 'iter-detailed', 'TolX', 1e-10,...
        'TolFun', 1e-10, 'TolCon', 1e-10, 'MaxFunEval', 30000, 'Maxiter', 1e+03 );
    zo = 1*ones(3,1);
    options_steady_states=optimset('Display','iter','maxfunevals', 100000000,'maxiter', 1000); % for steady state
    P_steady_states = [P0(1) P0(2) Param_reduced(1) Param_reduced(2) P0(5)];
    g= @(x0) steady_states(x0,P_steady_states); 
    [z0,~] = fsolve(g,zo,options_steady_states); % finding equilibrium
    M = eye(2);
    Vref0 = z0(1); Eq0 = z0(2); V0 = z0(3);
    SVref0 = [0 0]*M;
    SEq0 = [0 0]*M;
    SV0 = [0 0]*M;
    Z0 =[Vref0;Eq0;V0;SVref0(:);SEq0(:);SV0(:)];

% Sequential Approach of Dynamic Optimization
    [ p_opt ] = fmincon( @(Param_reduced)obj(tspan,Z0,Param_reduced,P0,output_samples), Param_reduced, [], [], [], [],...
        [], [], @(Param_reduced)ctr(tspan,Z0,Param_reduced,P0,output_samples), optNLP);
    p_est = p_opt;
end

function [ J, dJ ] = obj(tspan,Z0,Param_reduced,P0,output_samples)
       [f,df] = fun(tspan,Z0,Param_reduced,P0,output_samples);
%        f = fun(tspan,Z0,Param_reduced,P0,output_samples);
        J = f;
        dJ = df;
 end


function [ c, ceq, dc, dceq ] = ctr(tspan,Z0,Param_reduced,P0,output_samples) %no constraints
    if nargout == 2
        f = fun( tspan,Z0,Param_reduced,P0,output_samples);
        ceq = [];
        c = [];
    else
        [f,df] = fun( tspan,Z0,Param_reduced,P0,output_samples);
        ceq = [];
        dceq = [];
        c = [];
        dc = [];
    end
end
function [ f, df ] = fun( tspan,z0,Param_reduced,P0,output_samples)
        k_directions = 2;
        M_DAE=eye(3+3*k_directions);
        M_DAE(3,3)=0;
        for i =(3+2*k_directions)+1: (3+3*k_directions)
            M_DAE(i,i)=0;
        end
        optODE=odeset('Mass',M_DAE,'RelTol',1e-8);
        f = (z0(3)-output_samples(1,3))^2;
        df = (2*(z0(3)-output_samples(1,3))*z0(8:9)')';
        for i=2:size(output_samples,1)
            time = [tspan(i-1) tspan(i)];
            output_sample_state = output_samples(i,3);
            [time_span,Z] = ode15s(@(t,Z) dynamics_reduced(t,Z,Param_reduced,P0),time,z0(:),optODE);
            f = f + (Z(end,3)-output_sample_state)^2 ;
            df = df + 2*(Z(end,3)-output_sample_state)*Z(end,8:9)';
            z0 = Z(end,:);
        end


end

function dZdt = dynamics_reduced(t,Z,Param_reduced,P0)
k_directions = 2;
Vref = Z(1);
Eq = Z(2);
V = Z(3);

KQi = P0(1);
KVi = P0(2);
R = Param_reduced(1);
X = Param_reduced(2);
E = P0(5);
% Xeq = Param(6);
Xeq = 0.8;

Qcmd = 0.195197282638509; %constant
Q = V*Eq/Xeq - V^2/Xeq;
I = 0.948572674561663;
P = V * I;
PFE = 10;
Qcmd = tan(PFE)*P;

dVrefdt = KQi*Qcmd - KQi*Q;
dEqdt = KVi*Vref - KVi*V;
dVdt = V^4 - 2*P*R*V^2 - 2*Q*X*V^2 - E^2*V^2 + R^2*P^2 + X^2*P^2 + R^2*Q^2 + X^2*Q^2;

SVref = Z(3+1:3+k_directions);
SEq = Z(3+k_directions+1:3+2*k_directions);
SV = Z(3+2*k_directions+1:3+3*k_directions);

dQdParam = [SV(1)*Eq/Xeq+V*SEq(1)-2*SV(1)*V/Xeq/Xeq;SV(2)*Eq/Xeq+V*SEq(2)/Xeq]-2*SV(2)*V/Xeq;

dSVrefdt = [-KQi*dQdParam(1);-KQi*dQdParam(2)];
dSEqdt = [KVi*SVref(1)-KVi*SV(1);KVi*SVref(2)-KVi*SV(2)];
dSVdt = [4*V^3*SV(1)-2*P*V^2-2*P*R*2*V*SV(1)-4*2*Q*X*V*SV(1)-2*dQdParam(1)*X*V^2-2*E^2*V*SV(1)+2*R*P^2+0+2*R*Q^2+2*R^2*Q*dQdParam(1)+2*X^2*Q*dQdParam(1);...
         4*V^3*SV(2)-2*P*R*2*V*SV(2)-2*Q*V^2-4*2*Q*X*V*SV(2)-2*dQdParam(2)*X*V^2-2*E^2*V*SV(2)+0+2*X*P^2+2*R^2*Q*dQdParam(2)+2*X*Q^2+2*X^2*Q*dQdParam(2)];

dZdt = [dVrefdt;dEqdt;dVdt;dSVrefdt;dSEqdt;dSVdt];
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
% Xeq = Param(6);
Xeq = 0.8;

Qcmd = 0.195197282638509; %constant
Q = V*Eq/Xeq - V^2/Xeq;
I = 0.948572674561663;
P = V * I;
PFE = 10;
Qcmd = tan(PFE)*P;
       
f=[ ...
    KQi*Qcmd - KQi*Q;  
    KVi*Vref - KVi*V;
    V^4 - 2*P*R*V^2 - 2*Q*X*V^2 - E^2*V^2 + R^2*P^2 + X^2*P^2 + R^2*Q^2 + X^2*Q^2];

end