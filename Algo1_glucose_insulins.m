clear all;clc;close all;
np = 13; k_directions = np+1; threshold_percent_max = 0.05;
%Use LHS to construct the reference parameters set


T_f=[];
T_P_Matrix_f=[];

A = eye(np);
% Primary probing directions within the directions matrix P
D_set = {[A(:,1),A],[A(:,2),A],[A(:,3),A],[A(:,4),A],[A(:,5),A],[A(:,6),A],[A(:,7),A],[A(:,8),A],[A(:,9),A],[A(:,10),A],[A(:,11),A],[A(:,12),A],[A(:,13),A],...
         [-A(:,1),A],[-A(:,2),A],[-A(:,3),A],[-A(:,4),A],[-A(:,5),A],[-A(:,6),A],[-A(:,7),A],[-A(:,8),A],[-A(:,9),A],[-A(:,10),A],[-A(:,11),A],[-A(:,12),A],[-A(:,13),A]};
%  D_set = {[A(:,1),A]};
% Value of epsilon for twin probing (epsilon_twin)
epsilon_twin_probing = 0.01;
% Directions of perturbations of epsilon_twin
M_twin_probing_set = {[A(:,2),A],[A(:,1),A],[A(:,4),A],[A(:,3),A],[A(:,6),A],[A(:,5),A],[A(:,8),A],[A(:,7),A],[A(:,10),A],[A(:,9),A],[A(:,11),A],[A(:,12),A],[A(:,13),A],[A(:,12),A]}; 
% M_twin_probing_set = {[A(:,2),A]};
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
        p1 = 291.2; p2 = 0.0317; p3 = 0.0123; p4 = 4.92e-6; p5 = 0.0039; p6 = 79.0353; p7 = 0.2659; p8 = 357.8; p9 = 60; p10 = 7; p11 = 0.5; p12 = 0.05;
        p5 = 10; p6 = 0.0423;
        p5 = 20; p6 = 60; gamma = 5.9;
        P0 = [p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 gamma];
        % Initial conditions
        G0 = p1; X_state0 = 0; I0 = p8+p10;
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
        T_f=[T_f;T];
        T_P_Matrix_f = [T_P_Matrix_f;T_P_Matrix];
       end
    end
%% 
close all;
set(0,'DefaultFigureWindowStyle','docked');font_size = 30;
fig = figure('DefaultAxesFontSize',font_size,'defaultLineLineWidth',5,'DefaultAxesTitleFontWeight','bold');
plot(Time,Z(:,1),'b-');
fig = figure('DefaultAxesFontSize',font_size,'defaultLineLineWidth',5,'DefaultAxesTitleFontWeight','bold');
plot(Time,Z(:,3),'r-');    
fig = figure('DefaultAxesFontSize',font_size,'defaultLineLineWidth',5,'DefaultAxesTitleFontWeight','bold');
plot(Time,(p1)*Z(:,33)./Z(:,3),'b-'); 
hold on
plot(Time,(p2)*Z(:,34)./Z(:,3),'r-');
hold on
plot(Time,(p3)*Z(:,35)./Z(:,3),'g-');
hold on
plot(Time,(p4)*Z(:,36)./Z(:,3),'m-');
hold on
plot(Time,(p5)*Z(:,37)./Z(:,3),'b--');
hold on
plot(Time,(p6)*Z(:,38)./Z(:,3),'r--');
hold on
plot(Time,(p7)*Z(:,39)./Z(:,3),'k--');
hold on
plot(Time,(p8)*Z(:,40)./Z(:,3),'m--');
hold on
plot(Time,(p9)*Z(:,41)./Z(:,3),'b-.');
hold on
plot(Time,(p10)*Z(:,42)./Z(:,3),'r-.');
hold on
plot(Time,(p11)*Z(:,43)./Z(:,3),'g-.');
hold on
plot(Time,(p12)*Z(:,44)./Z(:,3),'m-.');
hold on
plot(Time,gamma*Z(:,45)./Z(:,3),'k-.');
hold on
xline(20,'k-','LineWidth',3)
hold on
xline(60,'k-','LineWidth',3)
ylim([-2 4])
xlim([0 180])
legend('$p_1$','$p_2$','$p_3$','$p_4$','$p_5$','$p_6$','$p_7$','$p_8$','$p_9$','$p_{10}$','$p_{11}$','$p_{12}$','$\gamma$','Interpreter','Latex')
disp('\Theta_{lni}}')
Theta_lni
disp('\Theta_{li}}')
Theta_li = setdiff(Theta,Theta_lni)

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
% u = p5*t*max(0,G-p6);
% u = p5*exp(-p6*t);
u = median([0,gamma*(t-p6)/(p5-p6),gamma]);

dGdt = -p2*(G-p9)-G*X + p11*exp(-p12*t);
dXdt = -p3*X + p4*(I-p10);
dIdt = u - p7*(I-p10);

SG = Z(3+1:3+k_directions);
SX = Z(3+k_directions+1:3+2*k_directions);
SI = Z(3+2*k_directions+1:3+3*k_directions);

% SLmax_output = SLmax(zeros(14,1),[G-p6;SG-M(6,:)']);
% u_prime = t*max(0,G-p6)* M(5,:)'+p5*t*SLmax_output;

% u_prime =  exp(-p6*t)*M(5,:)' - p5*t*exp(-p6*t)*M(6,:)';

u_prime = SLmid(zeros(15,1),[gamma*(t-p6)/(p5-p6);(t-p6)/(p5-p6)*M(13,:)' - gamma*(t-p6)/(p5-p6)^2*M(5,:)' + gamma*(t-p5)/(p5-p6)^2*M(6,:)'],...
                [gamma;M(13,:)']);

dSGdt = -(G-p9)*M(2,:)' + p2*M(9,:)' -(p2+X)*SG-G*SX + exp(-p12*t)*(M(11,:)'-p11*t*M(12,:)');
dSXdt = -X*p3 + (I-p10)*M(4,:)' - p4*M(10,:)' - p3*SX + p4*SI;
dSIdt = u_prime - (I-p10)*M(7,:)' + p7*M(10,:)' - p7*SI;

dZdt = [dGdt;dXdt;dIdt;dSGdt;dSXdt;dSIdt];
end
function SLmax_output = SLmax(xM1,yM2)
% Shifted Lmax 
% Returns the lexicographic maximum of 2 vectors xM1&yM2, left-shifted by
% one element
% Vectors xM1&yM2 must have the same size
 SLmax_output = xM1(2:end);
 for k=1:size(xM1,2)
    if xM1(k)>yM2(k)
       SLmax_output = xM1(2:end);
       break
    elseif xM1(k)<yM2(k)
       SLmax_output = yM2(2:end);
       break
    end
 end
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