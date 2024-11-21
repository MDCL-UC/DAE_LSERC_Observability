clear all;clc;close all;
np = 13; k_directions = np+1; threshold_percent_max = 0.05;
%Use LHS to construct the reference parameters set
p1 = 291.2; p2 = 0.0317; p3 = 0.0123; p4 = 4.92e-6; p5 = 20; p6 = 60; p7 = 0.2659; p8 = 357.8; p9 = 60; p10 = 7; p11 = 0.5; p12 = 0.05; gamma = 5.9;
P0_true = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 gamma];
pL = 10^(-1); pU = 10;
P_set = {P0_true,pL*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true, ...
        (pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true, ...
        (pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true, ...
        (pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true, ...
        (pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,(pL+pU*rand)*P0_true,pU*P0_true};

nonidentifiable_param_index = [];
T_f=[];
T_P_Matrix_f=[];

for pp = 1:length(P_set)
    P0 = cell2mat(P_set(pp));
    Param = P0;
    % Primary probing directions within the directions matrix P
    A = eye(np);
    D_set = {[A(:,1),A],[-A(:,1),A],[A(:,2),A],[-A(:,2),A]};
    % Value of epsilon for twin probing (epsilon_twin)
    epsilon_twin_probing = 0.01;
    % Directions of perturbations of epsilon_twin
    M_twin_probing_set = {[A(:,2),A],[A(:,1),A]};
    % Value of epsilon for singular probing(epsilon_sing)
    epsilon_sing_probing = 0.01;
    D_sing=[];
    V_sing=[];
    Theta = [1 2 3 4 5 6 7 8 9 10 11 12 13];
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
        G0 = P0(1); X_state0 = 0; I0 = P0(8)+P0(10);
        SG0 = M(1,:);
        SX0 = zeros(1,14);
        SI0 = M(8,:)+M(10,:);
        Z0 =[G0;X_state0;I0;SG0(:);SX0(:);SI0(:)]; %vector of states & state sensitivies to be used in ODE solver
        Z0_twin_probing =[G0;X_state0;I0;SG0(:);SX0(:);SI0(:)];
        % Time interval/steps
        t0 = 0; tf = 180;
        N = np; 
        tspan = t0:1/(2*N):tf;
        [Time,Z] = ode45(@(t,Z) dynamics(t,Z,P0,M),tspan,Z0(:));
        [T_twin_probing,Z_twin_probing] = ode45(@(t_twin_probing,Z_twin_probing) dynamics(t_twin_probing,Z_twin_probing,P0,M_twin_probing),tspan,Z0_twin_probing(:));
        
        StepCount = length(Time);
        % Calculate the the SVD for the LSERC matrix to check identifiability
        h = Z(:,1); %the output is the voltage
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
            xp = Z(ti,3+k_directions+1:3+2*k_directions);
            xp_twin_probing = Z_twin_probing(ti,3+k_directions+1:3+2*k_directions);
            dhdp = xp;
            dhdp_twin_probing = xp_twin_probing;
            LD_Y_theta(step,:) = dhdV * xp ;
            LD_Y_theta_twin_probing(step,:) = dhdV * xp_twin_probing ;
            step = step+1;
        end
        states_solution = [Z(:,1) Z(:,2) Z(:,3)];
        [U_LD,S_LD,V_LD] = svd(LD_Y_theta);
        [U_LD_twin_probing,S_LD_twin_probing,V_LD_twin_probing] = svd(LD_Y_theta_twin_probing);

        L_Y_theta = LD_Y_theta/M;
        L_Y_theta_twin_probing = LD_Y_theta_twin_probing/M_twin_probing;
        [U_L,S_L,V_L] = svd(L_Y_theta);
        diag_SL = diag(S_L);
        P_sing_vec = [];
        for i=1:length(diag_SL)
            if diag_SL(i)<1e-10
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
np = 13; k_directions = np+1; threshold_percent_max = 0.05;
%Use LHS to construct the reference parameters set
p1 = 291.2; p2 = 0.0317; p3 = 0.0123; p4 = 4.92e-6; p5 = 20; p6 = 60; p7 = 0.2659; p8 = 357.8; p9 = 60; p10 = 7; p11 = 0.5; p12 = 0.05; gamma = 5.9;
P0 = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 gamma];
Param_reduced = [gamma+0.1*gamma];
gamma_avg = 28.5528;
Param_reduced = gamma_avg;
M = eye(np);
M = [M(:,1),M];

G0 = P0(1); X_state0 = 0; I0 = P0(8)+P0(10);
SG0 = M(1,:);
SX0 = zeros(1,14);
SI0 = M(8,:)+M(10,:);
Z0 =[G0;X_state0;I0;SG0(:);SX0(:);SI0(:)];
   
% Time interval/steps
t0 = 0; tf = 180;
N = np; 
tspan = t0:1/(N):tf;
[Time,Z_true] = ode45(@(t,Z) dynamics(t,Z,P0,M),tspan,Z0(:));

output_samples = Z_true(:,3);
   
p_est = param_est_oc(Param_reduced,P0,tspan,output_samples);
         
%% 
function dZdt = dynamics(t,Z,P0,M)
k_directions = size(M,2);

G = Z(1);
X = Z(2);
I = Z(3);

p1 = P0(1);
p2 = P0(2);
p3 = P0(3);
p4 = P0(4);
p5 = P0(5);
p6 = P0(6);
p7 = P0(7);
p8 = P0(8);
p9 = P0(9);
p10 = P0(10);
p11 = P0(11);
p12 = P0(12);
gamma = P0(13);

u = median([0,gamma*(t-p6)/(p5-p6),gamma]);

dGdt = -p2*(G-p9)-G*X + p11*exp(-p12*t);
dXdt = -p3*X + p4*(I-p10);
dIdt = u - p7*(I-p10);

SG = Z(3+1:3+k_directions);
SX = Z(3+k_directions+1:3+2*k_directions);
SI = Z(3+2*k_directions+1:3+3*k_directions);


u_prime = SLmid(zeros(15,1),[gamma*(t-p6)/(p5-p6);(t-p6)/(p5-p6)*M(13,:)' - gamma*(t-p6)/(p5-p6)^2*M(5,:)' + gamma*(t-p5)/(p5-p6)^2*M(6,:)'],...
                [gamma;M(13,:)']);

dSGdt = -(G-p9)*M(2,:)' + p2*M(9,:)' -(p2+X)*SG-G*SX + exp(-p12*t)*(M(11,:)'-p11*t*M(12,:)');
dSXdt = -X*p3 + (I-p10)*M(4,:)' - p4*M(10,:)' - p3*SX + p4*SI;
dSIdt = u_prime - (I-p10)*M(7,:)' + p7*M(10,:)' - p7*SI;

dZdt = [dGdt;dXdt;dIdt;dSGdt;dSXdt;dSIdt];
end

function p_est = param_est_oc(Param_reduced,P0,tspan,output_samples)
   Algorithm = 'sqp';
   optNLP = optimset( 'Algorithm',Algorithm,'Hessian','bfgs','LargeScale', 'off', 'GradObj', 'on', 'GradConstr', 'off',...
        'DerivativeCheck', 'off', 'Display', 'iter-detailed', 'TolX', 1e-14,...
        'TolFun', 1e-14, 'TolCon', 1e-14, 'MaxFunEval', 30000, 'Maxiter', 1e+03 );
    
    M = 1;
    G0 = P0(1); X_state0 = 0; I0 = P0(8)+P0(10);
    SG0 = 0;
    SX0 = 0;
    SI0 = 0;
    Z0 =[G0;X_state0;I0;SG0(:);SX0(:);SI0(:)];

% Sequential Approach of Dynamic Optimization
    [ p_opt ] = fmincon( @(Param_reduced)obj(tspan,Z0,Param_reduced,P0,output_samples), Param_reduced, [], [], [], [],...
        [], [], @(Param_reduced)ctr(tspan,Z0,Param_reduced,P0,output_samples), optNLP);
    p_est = p_opt;
end

function [ J, dJ ] = obj(tspan,Z0,Param_reduced,P0,output_samples)
       [f,df] = fun(tspan,Z0,Param_reduced,P0,output_samples);
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
        k_directions = 1;

        f = (z0(3)-output_samples(1))^2;
        df = (2*(z0(3)-output_samples(1))*z0(6)')';
        for i=2:size(output_samples,1)
            time = [tspan(i-1) tspan(i)];
            output_sample_state = output_samples(i);
            [time_span,Z] = ode45(@(t,Z) dynamics_reduced(t,Z,Param_reduced,P0),time,z0(:));
            f = f + (Z(end,3)-output_sample_state)^2 ;
            df = df + 2*(Z(end,3)-output_sample_state)*Z(end,6)';
            z0 = Z(end,:);
        end
    

end

function dZdt = dynamics_reduced(t,Z,Param_reduced,P0)
k_directions = 1;
G = Z(1);
X = Z(2);
I = Z(3);

p1 = P0(1);
p2 = P0(2);
p3 = P0(3);
p4 = P0(4);
p5 = P0(5);
p6 = P0(6);
p7 = P0(7);
p8 = P0(8);
p9 = P0(9);
p10 = P0(10);
p11 = P0(11);
p12 = P0(12);
gamma = Param_reduced(1);

if t>60
    break_here = 'X';
end
u = median([0,gamma*(t-p6)/(p5-p6),gamma]);

dGdt = -p2*(G-p9)-G*X + p11*exp(-p12*t);
dXdt = -p3*X + p4*(I-p10);
dIdt = u - p7*(I-p10);

SG = Z(3+1:3+k_directions);
SX = Z(3+k_directions+1:3+2*k_directions);
SI = Z(3+2*k_directions+1:3+3*k_directions);


u_prime = SLmid([0;0],[gamma*(t-p6)/(p5-p6);(t-p6)/(p5-p6)*1],...
                [gamma;1]);

dSGdt = 0 + 0 -(p2+X)*SG-G*SX + exp((-p12*t)*0-p11*t*0);
dSXdt = -X*p3 + (I-p10)*0 - p4*0 - p3*SX + p4*SI;
dSIdt = u_prime - (I-p10)*0 + p7*0 - p7*SI;

dZdt = [dGdt;dXdt;dIdt;dSGdt;dSXdt;dSIdt];
end
function SLmid_output = SLmid(xM,yM,zM)
% Shifted Lmid 
% Returns the lexicographic mid of 3 vectors arguments, left-shifted by
% one element
% All Vectors xM,yM and zM must have the same size
 SLmid_output = yM(2:end);
 for k=1:size(yM,2)
    if xM(k)>yM(k)
       SLmid_output = xM(2:end);
       if zM(k)<xM(k)
           SLmid_output = zM(2:end);
       end
       break
    elseif xM(k)<yM(k)
       SLmid_output = yM(2:end);
       if zM(k)<yM(k)
           SLmid_output = zM(2:end);
       end
       break
    end
 end
end