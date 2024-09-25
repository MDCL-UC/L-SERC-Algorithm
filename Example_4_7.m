clear all
% clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
np = 2; %number of parameters theta
Param = zeros(np,1);
Param(1) = 1; %p1
Param(2) = 1; %p2
% States
n = 1; %number of states
x0 = Param(2); %initial conditions
D_sing=[];
V_sing=[];

% Setup primary probing directions within the directions matrix P
P_Matrix_set = {[[1 0]',eye(2)],[[-1 0 ]',eye(2)]...
               ,[[0 1]',eye(2)],[[0 -1]',eye(2)]};
% Setup epsilon for twin probing, epsilon_twin
epsilon_twin_probing = 0.01; %epsilon_twin value
Matrix_twin_probing_set = {[[1 0]',zeros(2)],[[0 1]',zeros(2)]}; %directions of perturbations of epsilon_twin
T_f=[];
T_P_Matrix_f=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Primary & twin Probing stages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(P_Matrix_set)
   for k_twin_prob = 1:length(Matrix_twin_probing_set)
    P_Matrix = cell2mat(P_Matrix_set(k)); %Get current directions matrix
    P_Matrix_twin_probing = cell2mat(Matrix_twin_probing_set(k_twin_prob)); %Get current epsilon_twin
    if abs(P_Matrix(:,1)) == P_Matrix_twin_probing(:,1)
        continue
    end
    P_Matrix_twin_probing(:,1) = P_Matrix_twin_probing(:,1).*epsilon_twin_probing;
    if sum(sign(P_Matrix(:,1))) == 1 %check the sign of the sum of elements of the first column of P
       P_Matrix_twin_probing = P_Matrix_twin_probing+P_Matrix; 
    elseif sum(sign(P_Matrix(:,1))) == -1
       P_Matrix_twin_probing = -P_Matrix_twin_probing+P_Matrix;
    end    
    %start solving the forward sensitivity system for every direction in
    %primary and twin probing directions matrices
    xp0 = [P_Matrix(2,1) P_Matrix(2,2) P_Matrix(2,3)]; %dx/dp initial conditions
    xp0_twin_probing = [P_Matrix_twin_probing(2,1) P_Matrix_twin_probing(2,2) P_Matrix_twin_probing(2,3)]; %dx/dp initial conditions
    X0 =[x0;xp0(:)]; %vector of states & state sensitivies to be used in ODE solver
    X0_twin_probing =[x0;xp0_twin_probing(:)]; %vector of states & state sensitivies to be used in ODE solver
    % Time interval/steps
    N = 2*np; 
    tspan = 0:1/(N):1;
    [Time,X] = ode45(@(t,X) Example_4_7_ODE_solve(t,X,Param,P_Matrix),tspan,X0(:));
    [T_twin_probing,X_twin_probing] = ode45(@(t_twin_probing,X_twin_probing) Example_4_7_ODE_solve(t_twin_probing,X_twin_probing,Param,P_Matrix_twin_probing),tspan,X0_twin_probing(:));
    StepCount = length(Time);
    
    % Calculate the the SVD for the LSERC matrix to check identifiability
    h = X(:,1);
    [m,dummy1]=size(h);
    dhdx = 1;
    Yp = zeros(StepCount,np);
    LD_Y_theta = zeros(StepCount,np+1);
    LD_Y_theta_twin_probing = zeros(StepCount,np+1);
    t=tspan(StepCount);
    step = 1;
    p1 = Param(1);
    p2 = Param(2);
    for ti=1:StepCount
        xp = [X(ti,2) X(ti,3) X(ti,4)];
        xp_twin_probing = [X_twin_probing(ti,2) X_twin_probing(ti,3) X_twin_probing(ti,4)];
        dhdp = xp;
        dhdp_twin_probing = xp_twin_probing;
        LD_Y_theta(step,:) = dhdx * xp ;
        LD_Y_theta_twin_probing(step,:) = dhdx * xp_twin_probing ;
        step = step+1;
    end
    solution = X(:,1);
    [U_LD,S_LD,V_LD] = svd(LD_Y_theta);
    [U_LD_twin_probing,S_LD_twin_probing,V_LD_twin_probing] = svd(LD_Y_theta_twin_probing);

    L_Y_theta = LD_Y_theta/P_Matrix;
    L_Y_theta_twin_probing = LD_Y_theta_twin_probing/P_Matrix_twin_probing;
    [U_L,S_L,V_L] = svd(L_Y_theta);
    [U_L_twin_probing,S_L_twin_probing,V_L_twin_probing] = svd(L_Y_theta_twin_probing);
    %Save a record of the directions that resulted in a rank deficient LSERC
    %matrix
    if rank(S_L)<np
        zero_vector_index = find(all(S_L==0));
        d_sing = P_Matrix(:,1);
        v_sing = V_L(:,zero_vector_index);
        D_sing = [D_sing d_sing];
        V_sing = [V_sing v_sing];
    end
    %Collect results in tables for easier access
    T = array2table(LD_Y_theta);
    T_P_Matrix = array2table(P_Matrix);
    T.LD_Y_theta_twin_probing = LD_Y_theta_twin_probing;
    T.L_Y_theta = L_Y_theta;
    T.L_Y_theta_twin_probing = L_Y_theta_twin_probing;
    T.S_LD = S_LD;
    T.S_LD_twin_probing = S_LD_twin_probing;
    T.S_L = S_L;
    T.S_L_twin_probing = S_L_twin_probing;
    T.solution = solution;
    T_f=[T_f;T];
    T_P_Matrix_f = [T_P_Matrix_f;T_P_Matrix];
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Singularity probing stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check_empty=isempty(V_sing);
if ne(check_empty,1)
    epsilon_sing_probing = 0.01;
    Param_sing = [Param+epsilon_sing_probing*V_sing(:,1)];
    D_sing = unique(D_sing.','rows').';
    Matrix_sing_probing_set = {[D_sing(:,1),eye(2)]}; %construct the singularity probing directions matrix from directions saved in the previous stages
    if size(D_sing,2)>1
    for i=2:size(D_sing,2)
        Matrix_sing_probing_set{end+1} = [D_sing(:,i),eye(2)];
        Param_sing = [Param_sing Param+epsilon_sing_probing*V_sing(:,i) Param-epsilon_sing_probing*V_sing(:,i)];
    end
    end
    Param_sing = unique(Param_sing.','rows').';

    T_f_sing=[];
    for k = 1:length(Matrix_sing_probing_set)
       for i=1:size(Param_sing,2)
        P_Matrix_sing = cell2mat(Matrix_sing_probing_set(k));
        %start solving the forward sensitivity system for every direction in
        %the singularity probing directions matrix
        Param = Param_sing(:,i);
        xp0_sing = [P_Matrix_sing(2,1) P_Matrix_sing(2,2) P_Matrix_sing(2,3)]; %dx/dp initial conditions
        X0_sing =[x0;xp0_sing(:)]; %vector of states & state sensitivies to be used in ODE solver
        % Time interval/steps
        N = 2*np; 
        tspan = 0:1/(N):1;
        [T_sing,X_sing] = ode45(@(t_sing,X_sing) Example_4_7_ODE_solve(t_sing,X_sing,Param,P_Matrix_sing),tspan,X0_sing(:));
        StepCount = length(T_sing);
        % Calculate the the SVD for the LSERC matrix to check identifiability
        h = X(:,1);
        [m,dummy1]=size(h);
        dhdx = 1;
        Yp = zeros(StepCount,np);
        LD_Y_theta_sing = zeros(StepCount,np+1);
        t=tspan(StepCount);
        step = 1;
        p1 = Param(1);
        p2 = Param(2);
        for ti=1:StepCount
            xp_sing = [X_sing(ti,2) X_sing(ti,3) X_sing(ti,4)];
            dhdp = xp_sing;
            LD_Y_theta_sing(step,:) = dhdx * xp_sing ;
            step = step+1;
        end
        solution_sing = X_sing(:,1);
        [U_LD_sing,S_LD_sing,V_LD_sing] = svd(LD_Y_theta_sing);
        L_Y_theta_sing = LD_Y_theta_sing/P_Matrix_sing;
        [U_L_sing,S_L_sing,V_L_sing] = svd(L_Y_theta_sing);
        %Collect results in tables for easier access
        T = array2table(LD_Y_theta_sing);
        T.L_Y_theta_sing = L_Y_theta_sing;
        T.S_LD_sing = S_LD_sing;
        T.S_L_sing = S_L_sing;
        T.solution_sing = solution_sing;
        T_f_sing=[T_f_sing;T];
       end
    end
end
if (rank(S_L)==size(S_L,2)) || (rank(S_L_sing)==size(S_L_sing,2))
    disp('L-SERC identifiable');
else
    disp('not L-SERC identifiable');
end
function dXdt = Example_4_7_ODE_solve(t,X,Param,P_Matrix)
% x value
x = X(1);
% Parameters values
p1 = Param(1);
p2 = Param(2);

f = 0;
g = 1-exp(-1*(x-p1));

Jf_x = 0;
Jf_p = [0 , 0];
Jg_x = exp(p1 - x);
Jg_p = [ -exp(p1 - x), 0 ];


dxdt = max(f,g);

X_1 = X(2);
X_2 = X(3);
X_3 = X(4);
N = [X_1 X_2 X_3];

Slmax_input = [f Jf_p*P_Matrix+Jf_x*N;g Jg_p*P_Matrix+Jg_x*N];

Slmax_output = SLmax(Slmax_input(1,:),Slmax_input(2,:));

% ODE system2
dxpdt1 = Slmax_output(1);
dxpdt2 = Slmax_output(2);
dxpdt3 = Slmax_output(3);

dxpdt = [dxpdt1;dxpdt2;dxpdt3];
dXdt = [dxdt;dxpdt];
end