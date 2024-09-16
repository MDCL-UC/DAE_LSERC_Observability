clear all;clc;close all;
np = 4; k_directions = 5; threshold_percent_max = 0.05;
%Use LHS to construct the reference parameters set
P0 = [0 30 1 -2];

M_DAE=eye(2+2*k_directions);
M_DAE(2,2)=0;
for i =(k_directions+2)+1: 2+2*k_directions
    M_DAE(i,i)=0;
end
optODE=odeset('Mass',M_DAE,'RelTol',1e-8);
T_f=[];
T_P_Matrix_f=[];

% Primary probing directions within the directions matrix P
D_set = {[[1 0 0 0]',eye(4)],[[0 1 0 0]',eye(4)],[[0 0 1 0]',eye(4)],[[0 0 0 1]',eye(4)]...
         [[-1 0 0 0]',eye(4)],[[0 -1 0 0]',eye(4)],[[0 0 -1 0]',eye(4)],[[0 0 0 -1]',eye(4)]};
% Value of epsilon for twin probing (epsilon_twin)
epsilon_twin_probing = 0.01;
% Directions of perturbations of epsilon_twin
M_twin_probing_set = {[[0 1 0 0]',zeros(4)],[[1 0 0 0]',zeros(4)],[[0 0 0 1]',zeros(4)],[[0 0 1 0]',zeros(4)]}; 
% Number of iterations for the singularity probing stage
q = 1;
% Value of epsilon for singular probing(epsilon_sing)
epsilon_sing_probing = 0.01;
D_sing=[];
V_sing=[];
Theta = [1 2 3 4];
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
        x0 = P0(3)*P0(4); 
        w0 = P0(2)-abs(x0);
        X0 = [0 0 P0(4) P0(3)]*M;
        W0 = [0 1 0 0]*M - fsign(x0,X0)*X0;
        X0_twin_probing = [0 0 P0(4) P0(3)]*M_twin_probing;
        W0_twin_probing = [0 1 0 0]*M_twin_probing - fsign(x0,X0_twin_probing)*X0_twin_probing;
        Z0 =[x0;w0;X0(:);W0(:)]; %vector of states & state sensitivies to be used in ODE solver
        Z0_twin_probing =[x0;w0;X0_twin_probing(:);W0_twin_probing(:)]; %vector of states & state sensitivies to be used in ODE solver
        % Time interval/steps
        t0 = 0; tf = 1;
        N = 2*np; 
        t_star = 1/P0(1)*log(P0(2)/(P0(2)+P0(3)*P0(4)));
        tspan = t0:1/(2*N):tf;
        [Time,Z] = ode15s(@(t,Z) dynamics(t,Z,P0,M),tspan,Z0(:),optODE);
        [T_twin_probing,Z_twin_probing] = ode15s(@(t_twin_probing,Z_twin_probing) dynamics(t_twin_probing,Z_twin_probing,P0,M_twin_probing),tspan,Z0_twin_probing(:),optODE);
        StepCount = length(Time);
        % Calculate the the SVD for the LSERC matrix to check identifiability
        h = Z(:,1);
        [m,dummy1]=size(h);
        dhdx = 1;
        Yp = zeros(StepCount,np);
        LD_Y_theta = zeros(StepCount,np+1);
        LD_Y_theta_twin_probing = zeros(StepCount,np+1);
        t=tspan(StepCount);
        step = 1;
        p1 = P0(1);
        p2 = P0(2);
        for ti=1:StepCount
            %We need to calculate relative sensitivities  
            xp = [Z(ti,3) Z(ti,4) Z(ti,5) Z(ti,6) Z(ti,7)];
            xp_twin_probing = [Z_twin_probing(ti,3) Z_twin_probing(ti,4) Z_twin_probing(ti,5) Z_twin_probing(ti,6) Z_twin_probing(ti,7)];
            dhdp = xp;
            dhdp_twin_probing = xp_twin_probing;
            LD_Y_theta(step,:) = dhdx * xp ;
            LD_Y_theta_twin_probing(step,:) = dhdx * xp_twin_probing ;
            step = step+1;
        end
        solution = Z(:,1);
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
        if q>0
            %Save a record of the directions that resulted in a rank deficient LSERC
            %matrix
            if rank(S_L)<np
                zero_vector_index = find(all(S_L==0));
                d_sing = M(:,1);
                v_sing = V_L(:,zero_vector_index);
                D_sing = [D_sing d_sing];
                V_sing = [V_sing v_sing];
            end
        end
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
        T.solution = solution;
        T.algebraic = Z(:,2);
        T_f=[T_f;T];
        T_P_Matrix_f = [T_P_Matrix_f;T_P_Matrix];
       end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Singularity probing stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while q>0
        check_empty=isempty(V_sing);
        if ne(check_empty,1)
            epsilon_sing_probing = 0.01;
            Param_sing = [P0+epsilon_sing_probing*V_sing(:,1)];
            D_sing = unique(D_sing.','rows').';
            Matrix_sing_probing_set = {[D_sing(:,1),eye(4)]}; %construct the singularity probing directions matrix from directions saved in the previous stages
            T_f_sing=[];
            for k = 1:length(Matrix_sing_probing_set)
               for i=1:length(Param_sing)
                M_Matrix_sing = cell2mat(Matrix_sing_probing_set(k));
                x0 = P0(3)*P0(4); 
                w0 = P0(2)-abs(x0);
                X0 = [0 0 P0(4) P0(3)]*M_Matrix_sing;
                W0 = [0 1 0 0]*M_Matrix_sing - fsign(x0,X0)*X0;
                Z0 =[x0;w0;X0(:);W0(:)]; %vector of states & state sensitivies to be used in ODE solver
                Z0_twin_probing =[x0;w0;X0_twin_probing(:);W0_twin_probing(:)]; %vector of states & state sensitivies to be used in ODE solver
                % Time interval/steps
                t0 = 0; tf = 1;
                N = 2*np; 
                t_star = 1/P0(1)*log(P0(2)/(P0(2)+P0(3)*P0(4)));
                tspan = t0:1/(2*N):tf;
                [Time,Z] = ode15s(@(t,Z) dynamics(t,Z,P0,M_Matrix_sing),tspan,Z0(:),optODE);
                StepCount = length(Time);
                % Calculate the the SVD for the LSERC matrix to check identifiability
                h = Z(:,1);
                [m,dummy1]=size(h);
                dhdx = 1;
                Yp = zeros(StepCount,np);
                LD_Y_theta = zeros(StepCount,np+1);
                t=tspan(StepCount);
                step = 1;
                p1 = P0(1);
                p2 = P0(2);
                for ti=1:StepCount
                    %We need to calculate relative sensitivities  
                    xp = [Z(ti,3) Z(ti,4) Z(ti,5) Z(ti,6) Z(ti,7)];
                    dhdp = xp;
                    LD_Y_theta(step,:) = dhdx * xp ;
                    step = step+1;
                end
                solution = Z(:,1);
                [U_LD,S_LD,V_LD] = svd(LD_Y_theta);

                L_Y_theta = LD_Y_theta/M_Matrix_sing;
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
                for i=1:length(C)
                    if find(Theta_li==1)~= 1
                        Theta_lni = [Theta_lni C];
                    end
                end
               end
            end
                end
               q = q-1;
            end

    
disp('\Theta_{lni}}')
Theta_lni
disp('\Theta_{li}}')
Theta_li = setdiff(Theta,Theta_lni)

function dZdt = dynamics(t,Z,Param,M)
t_star = 1/Param(1)*log(Param(2)/(Param(2)+Param(3)*Param(4)));
k_directions = size(M,2);
x = Z(1);
w = Z(2);

p1 = Param(1);
p2 = Param(2);
p3 = Param(3);
p4 = Param(4);

dxdt = max(p1,0)*w;
dwdt = -x+w-p2;
if t>t_star
    dwdt = -(-x+w-p2);
end
X = Z(2+1:2+k_directions);
W = Z(2+k_directions+1:2+2*k_directions);
xM1 = [0, 0, 0, 0, 0, 0];
yM2 = [p1*w, M(1,:)*w + p1*W(:)'];
dXdt = SLmax(xM1,yM2)';
% dWdt = fsign(x,X(:))*X(:) + W(:) - M(2,:)';
% dXdt = M(1,:)'*w + p1*W(:);
dWdt = -X(:) + W(:) - M(2,:)';
if t>t_star
    dWdt = -(-X(:) + W(:) - M(2,:)');
end
dZdt = [dxdt;dwdt;dXdt;dWdt];
end
