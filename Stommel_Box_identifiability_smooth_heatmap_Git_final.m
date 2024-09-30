clear
clc
close all
format long
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Initialization section

% Choose high_precision_computation == 'Y' to use HPF toolbox
% (much slower code but allows for smaller ||alpha|| values)
% or choose high_precision_computation == 'N' for standard Matlab precision
% (minimum value allowed for accurate results is ||alpha|| = 1e-6
high_precision_computation = 'N';
% Number of parameters (theta) of the system
ThetaCount = 3;
% Values of reference parameters
Param = zeros(ThetaCount,1);
Param(1) = 3; %theta1
Param(2) = 1.1; %theta2
Param(3) = 0.3; %theta3
% Number of states of the system
n = 2; 
% Initial conditions
T_0=1; 
V_0=2; 
x0 = [T_0;V_0];
% Reference input values
B = 2;
B_hat = 1;
A = B - B_hat;
omega = A/0.05;
u1_ref = B*sin(omega*0);
u2_ref = B_hat*sin(omega*0);
%Primary directions matrix
M_Matrix = [[1 0 0]',eye(3)]; %M=[e_1 I_3]
% Value of epsilon for twin probing (epsilon_twin), turned off in this case
epsilon_twin_probing = 0;
% Number of iterations for the singularity probing stage
q = 0;
% Value of epsilon for singular probing(epsilon_sing)
epsilon_sing_probing = 0;
X0 = zeros(2,4);

%vector of states & state sensitivies to be used in ODE solver
X0_All =[x0;X0(:)];

% Time interval/steps
N = 2*ThetaCount;
multiplier = 1000;
N = N*multiplier;
tspan = 0:1/(N):1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Nonsmooth section
%Solving LD-sensistivity matrix
[Time,X] = ode45(@(t,X) Nonsmooth_Stommel_Box_model_ODE(t,X,Param,M_Matrix),tspan,X0_All(:));

% Ytheta identifiability check matrix
%h=max
m=1;
StepCount = length(Time);
number_of_samples = 8;
time_sample_increment = (length(Time)-1)/number_of_samples;
t_samples = [1];
for s=1:number_of_samples
    t_samples = [t_samples 1+s*time_sample_increment];
end
I = find(Time == 0.75); %nonsmoothness point
t_samples = [t_samples I];
t_samples=sort(t_samples);
t_samples = ceil(t_samples);

Y_theta_case_id = zeros(StepCount*m,ThetaCount+1);

Y_theta_samples_case_id = [];

t=tspan(StepCount);
step = 1;
step_sample = 1;
Time_subset_plot=[];
alpha_vec_norm = [];
%This loop to solve for y in each time step
for ti=1:StepCount
    T_min_case_id = 1.5;
    h_case_id = max(X(ti,1),T_min_case_id);
    y_output_case_id(step:step+m-1) = h_case_id;
    dh_case_id = SLmax([X(ti,1) X(ti,3) X(ti,5) X(ti,7) X(ti,9)],[T_min_case_id 0 0 0 0]);
    Y_theta_case_id(step:step+m-1,1:end) = dh_case_id;
    if ismember(ti, t_samples)
        Y_theta_samples_case_id = [Y_theta_samples_case_id;dh_case_id];
        Time_subset_plot = [Time_subset_plot Time(ti)];
        step_sample = step_sample+m;
    end
    step = step+m;
end
%Get SY and SVD identifiability
[U_LD_case_id,S_LD_case_id,V_LD_case_id] = svd(Y_theta_samples_case_id);
Sy_case_id = Y_theta_case_id*[0 0 0;eye(3)];
Sy_samples_case_id = Y_theta_samples_case_id*[0 0 0;eye(3)];
[U_L_case_id,S_L_case_id,V_L_case_id] = svd(Sy_samples_case_id);
%State variables
T_Sol = X(:,1);
V_Sol = X(:,2);
Diff_Sol = abs(T_Sol - V_Sol);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Smoothing results section
total_samples_smooth = 20;
alpha_vals = [1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,0];
choosen_alpha = 1e-1;
% alpha_vals = [1e-9,1e-10];
Rank_table = zeros(length(alpha_vals),total_samples_smooth);
alpha_vec = zeros(2,length(alpha_vals));
for ii=1:length(alpha_vals)
% Smoothing parameters
alpha1 = alpha_vals(ii)/sqrt(2);
alpha2 = alpha_vals(ii)/sqrt(2);
alpha_vec(:,ii) = [alpha1;alpha2];
alpha_vec_norm(ii) = vecnorm(alpha_vec(:,ii));
%Solving LD-sensistivity matrix
[Time_smooth,X_smooth] = ode45(@(t,X) Stommel_Box_model_ODE_smoothing(t,X,Param,M_Matrix,alpha1),tspan,X0_All(:));

% Ytheta identifiability check matrix
%h=max
m=1;
StepCount_smooth = length(Time_smooth);
for jj = 1:total_samples_smooth
number_of_samples_smooth = jj;
time_sample_increment_smooth = (length(Time_smooth)-1)/total_samples_smooth;
t_samples_smooth = [];
for s=1:number_of_samples_smooth
    t_samples_smooth = [t_samples_smooth 1+s*time_sample_increment_smooth];
end

t_samples_smooth=sort(t_samples_smooth);
t_samples_smooth = ceil(t_samples_smooth);
Y_theta_case_smooth = zeros(StepCount_smooth*m,ThetaCount+1);
Y_theta_samples_case_smooth = [];
t_smooth=tspan(StepCount_smooth);
step_smooth = 1;
step_sample_smooth = 1;
Time_subset_plot_smooth=[];

%This loop to solve for y in each time step
 for ti=1:StepCount_smooth
    T_min = 1.5;
    if alpha2 == 0
        h_case_smooth = max(X_smooth(ti,1),T_min);
        dh_case_smooth = SLmax([X_smooth(ti,1) X_smooth(ti,3) X_smooth(ti,5) X_smooth(ti,7) X_smooth(ti,9)],[T_min 0 0 0 0]);
    elseif  (alpha_vec_norm(ii) == 1e-9 || alpha_vec_norm(ii) == 1e-10) && (high_precision_computation == 'Y')
        h_case_smooth = hpf(1/2*(hpf(X_smooth(ti,1)+T_min+sqrt(hpf((X_smooth(ti,1)-T_min)^2)+hpf(alpha2^2)))));
        dh_case_smooth = [hpf(1/2*X_smooth(ti,3)*(1/(2*(X_smooth(ti,1)-T_min)))+0+1/2*1/sqrt(hpf((X_smooth(ti,1)-T_min)^2)+hpf(alpha2^2)))*(2*(X_smooth(ti,1)-T_min))...
                   1/2*X_smooth(ti,5)*(1/(2*(X_smooth(ti,1)-T_min))+0+1/2*1/sqrt(hpf((X_smooth(ti,1)-T_min)^2)+hpf(alpha2^2)))*(2*(X_smooth(ti,1)-T_min))...
                   1/2*X_smooth(ti,7)*(1/(2*(X_smooth(ti,1)-T_min))+0+1/2*1/sqrt(hpf((X_smooth(ti,1)-T_min)^2)+hpf(alpha2^2)))*(2*(X_smooth(ti,1)-T_min))...
                   1/2*X_smooth(ti,9)*(1/(2*(X_smooth(ti,1)-T_min))+0+1/2*1/sqrt(hpf((X_smooth(ti,1)-T_min)^2)+hpf(alpha2^2)))*(2*(X_smooth(ti,1)-T_min))];
    else
        h_case_smooth = 1/2*(X_smooth(ti,1)+T_min+sqrt((X_smooth(ti,1)-T_min)^2+alpha2^2));
        dh_case_smooth = [1/2*X_smooth(ti,3)*(1/(2*(X_smooth(ti,1)-T_min))+0+1/2*1/sqrt((X_smooth(ti,1)-T_min)^2+alpha2^2))*(2*(X_smooth(ti,1)-T_min))...
                   1/2*X_smooth(ti,5)*(1/(2*(X_smooth(ti,1)-T_min))+0+1/2*1/sqrt((X_smooth(ti,1)-T_min)^2+alpha2^2))*(2*(X_smooth(ti,1)-T_min))...
                   1/2*X_smooth(ti,7)*(1/(2*(X_smooth(ti,1)-T_min))+0+1/2*1/sqrt((X_smooth(ti,1)-T_min)^2+alpha2^2))*(2*(X_smooth(ti,1)-T_min))...
                   1/2*X_smooth(ti,9)*(1/(2*(X_smooth(ti,1)-T_min))+0+1/2*1/sqrt((X_smooth(ti,1)-T_min)^2+alpha2^2))*(2*(X_smooth(ti,1)-T_min))];

    end   
    y_output_case2(step_smooth:step_smooth+m-1) = h_case_smooth;
    Y_theta_case_smooth(step_smooth:step_smooth+m-1,1:end) = dh_case_smooth;

     if ismember(ti, t_samples_smooth)
        Y_theta_samples_case_smooth = [Y_theta_samples_case_smooth;dh_case_smooth];
        Time_subset_plot_smooth = [Time_subset_plot_smooth Time_smooth(ti)];
        step_sample_smooth = step_sample_smooth+m;
     end
    step_smooth = step_smooth+m;
end  


%Get SY and SVD case2
[U_LD_smooth,S_LD_smooth,V_LD_smooth] = svd(double(Y_theta_samples_case_smooth));
Sy_smooth = Y_theta_case_smooth*[0 0 0;eye(3)];
Sy1_smooth_all(:,ii) =  Sy_smooth(:,1);
Sy_samples_smooth = double(Y_theta_samples_case_smooth)*[0 0 0;eye(3)];
[U_L_smooth,S_L_smooth,V_L_smooth] = svd(Sy_samples_smooth);
Rank_table(ii,jj)=rank(S_L_smooth,1e-20000);

if alpha_vec_norm(ii)== choosen_alpha
    Sy_smooth_sens_plot = Sy_smooth;
    alpha_vec_norm_plot = alpha_vec_norm(ii);
end
end
end

Rank_table_logical=Rank_table>=1;
Rank_table_logical=double(Rank_table_logical);
Rank_table_logical(Rank_table_logical==1)=NaN;
[rows,cols]=find(Rank_table_logical==0);
Rank_table_logical(1:rows-1,cols)=1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% 
% Plotting section
close all
set(0,'DefaultFigureWindowStyle','docked')
% %1- Plot solution vs output
% %a)solution
figure_combined = figure('Name','Solution vs Output','DefaultAxesFontSize',24,'defaultLineLineWidth',4);
sol_vs_output_subplot=subplot(2,2,[1 2],'Parent',figure_combined);
plot(Time,T_Sol,'color',[0.9290 0.6940 0.1250],'Parent',sol_vs_output_subplot);
hold on
plot(Time,V_Sol,'b','Parent',sol_vs_output_subplot);
%b)output
hold on
plot(Time,y_output_case_id(1:m:end),'k--','LineWidth',5,'Parent',sol_vs_output_subplot);
xlabel('Time')
lgnd_sol_vs_output=legend('$T(t)$ ','$V(t)$ ','$y(t)$','Location','best','Interpreter','latex','FontSize',24);
ylim([0.9 2.2]);
xline(0.75,'k-.','LineWidth',2,'HandleVisibility','off');
xline(0.4,'k-.','LineWidth',2,'HandleVisibility','off');


%2- Plot S^L_y
sensitivity_subplot=subplot(2,2,3,'Parent',figure_combined);
subplot_title = append('$||\alpha|| = $ ',num2str(choosen_alpha));
plot2=plot(Time,Sy_case_id(1:m:end,1),'b',Time,Sy_case_id(1:m:end,2),'k',Time,Sy_case_id(1:m:end,3),'r' ...
    ,Time_smooth,Sy_smooth_sens_plot(:,1),'b:',Time_smooth,Sy_smooth_sens_plot(:,2),'k:',Time_smooth,Sy_smooth_sens_plot(:,3),'r:'...
    ,Time_subset_plot,Sy_samples_case_id(1:m:end,1),'b*',Time_subset_plot,Sy_samples_case_id(1:m:end,2),'k*',Time_subset_plot,Sy_samples_case_id(1:m:end,3),'r*'...
    ,'Parent',sensitivity_subplot,'MarkerSize',20);
xlabel('Time');
formatSpec = '%s';
% txt1 = '${\mathrm{col}}_1(S_y)$ ' + num2str(alpha_vec_norm_plot,formatSpec);
lgnd_S_Y_case_id=legend('${\mathrm{col}}_1(S^L_y)$','${\mathrm{col}}_2(S^L_y)$ ','${\mathrm{col}}_3(S^L_y)$'...
                       ,'${\mathrm{col}}_1(S_y)$','${\mathrm{col}}_2(S_y)$ ','${\mathrm{col}}_3(S_y)$'...
                       ,'Location','northwest','Interpreter','latex','FontSize',20);
ylim([-0.5 1]);
xline(0.75,'k-.','LineWidth',2,'HandleVisibility','off');
xline(0.4,'k-.','LineWidth',2,'HandleVisibility','off');
xticks(0:0.2:1)
title(subplot_title,'Interpreter','latex') 
% 3- Heatmap
heatmap_subplot=subplot(2,2,4,'Parent',figure_combined);
h = heatmap(Rank_table_logical);
map = [1 0 0;
       0.9290 0.6940 0.1250];
h.Colormap = map;
h.MissingDataColor = [0.4660 0.6740 0.1880];
h.CellLabelColor = 'none';
h.XLabel = 'Time';
h.YLabel = '||\alpha||';
h.FontSize = 24;
h.ColorbarVisible = 'off'; h.GridVisible='off';
alpha_vec_norm = vecnorm(alpha_vec);
alpha_vec_norm_cell = num2cell(alpha_vec_norm');
alpha_vec_norm_cell(2) = {'1e-01'};alpha_vec_norm_cell(3) = {'1e-02'};alpha_vec_norm_cell(4) = {'1e-03'};alpha_vec_norm_cell(5) = {'1e-04'};
h.XDisplayLabels(total_samples_smooth/5:total_samples_smooth/5:end) = num2cell(str2double(h.XDisplayLabels(total_samples_smooth/5:total_samples_smooth/5:end))./total_samples_smooth);
for ii = 1:length(h.XDisplayLabels)
    if str2double(h.XDisplayLabels(ii))>1
        h.XDisplayLabels(ii) = {''};
    end
end
h.XDisplayLabels(1) = {'0'};
h.YDisplayLabels = alpha_vec_norm_cell;
h.YDisplayLabels(end)={'L-SERC'};
s = struct(h);
s.XAxis.TickLabelRotation = 0;

function dXdt = Nonsmooth_Stommel_Box_model_ODE(t,X,Param,M_Matrix)
%states
T = X(1);
V = X(2);
% V = abs(T-V);
% Parameters values
theta1 = Param(1);
theta2 = Param(2);
theta3 = Param(3);
%Control inputs
%A/omega = 0.05 figure 1 dynamic tipping, values of B and B_hat, figure 6 
B = 2;
B_hat = 1;
A = B - B_hat;
omega = A/0.05;
u1 = B*sin(omega*t);
u2 = B_hat*sin(omega*t);
% ODE system1 (f-function)
dTdt = theta1 + u1 - T - T*abs(T-V); %f1
dVdt = theta2 + u2 - V*theta3 - V*abs(T-V); %f2

%S_x equations
X11 = X(3);
X21 = X(4);
X12 = X(5);
X22 = X(6);
X13 = X(7);
X23 = X(8);
X14 = X(9);
X24 = X(10);

%Directional derivative with respect to thetas
df1 = [1 0 0]*M_Matrix;
df2 = [0 1 0]*M_Matrix-[0 0 V]*M_Matrix;

%absolute function derivative
u = T-V;
U = [X11-X21 X12-X22 X13-X23 X14-X24];
fsign_u = fsign(u,U);

% ODE system2
X11_dt = df1(1) - X11 - abs(T-V)*X11 - T*fsign_u*(X11-X21);
X12_dt = df1(2) - X12 - abs(T-V)*X12 - T*fsign_u*(X12-X22);
X13_dt = df1(3) - X13 - abs(T-V)*X13 - T*fsign_u*(X13-X23);
X14_dt = df1(4) - X14 - abs(T-V)*X14 - T*fsign_u*(X14-X24);
X21_dt = df2(1) - theta3*X21 - abs(T-V)*X21 - V*fsign_u*(X11-X21);
X22_dt = df2(2) - theta3*X22 - abs(T-V)*X22 - V*fsign_u*(X12-X22);
X23_dt = df2(3) - theta3*X23 - abs(T-V)*X23 - V*fsign_u*(X13-X23);
X24_dt = df2(4) - theta3*X24 - abs(T-V)*X24 - V*fsign_u*(X14-X24);

dXdt = [dTdt;dVdt;X11_dt;X21_dt;X12_dt;X22_dt;X13_dt;X23_dt;X14_dt;X24_dt];
end
function dXdt = Stommel_Box_model_ODE_smoothing(t,X,Param,M_Matrix,alpha)
%states
T = X(1);
V = X(2);
% V = abs(T-V);
% Parameters values
theta1 = Param(1);
theta2 = Param(2);
theta3 = Param(3);
%Control inputs
%A/omega = 0.05 figure 1 dynamic tipping, values of B and B_hat, figure 6 
B = 2;
B_hat = 1;
A = B - B_hat;
omega = A/0.05;
u1 = B*sin(omega*t);
u2 = B_hat*sin(omega*t);

% Smoothing abs
smooth_abs = sqrt((T-V)^2+alpha^2);

% ODE system1 (f-function)
dTdt = theta1 + u1 - T - T*smooth_abs; %f1
dVdt = theta2 + u2 - V*theta3 - V*smooth_abs; %f2

%S_x equations
X11 = X(3);
X21 = X(4);
X12 = X(5);
X22 = X(6);
X13 = X(7);
X23 = X(8);
X14 = X(9);
X24 = X(10);

%Directional derivative with respect to thetas
df1 = [1 0 0]*M_Matrix;
df2 = [0 1 0]*M_Matrix-[0 0 V]*M_Matrix;

%smooth abs function derivative
grad_smooth_abs = (T-V)*((T-V)^2+alpha^2)^(-1/2);

% ODE system2
X11_dt = df1(1) - X11 - smooth_abs*X11 - T*grad_smooth_abs*(X11-X21);
X12_dt = df1(2) - X12 - smooth_abs*X12 - T*grad_smooth_abs*(X12-X22);
X13_dt = df1(3) - X13 - smooth_abs*X13 - T*grad_smooth_abs*(X13-X23);
X14_dt = df1(4) - X14 - smooth_abs*X14 - T*grad_smooth_abs*(X14-X24);
X21_dt = df2(1) - theta3*X21 - smooth_abs*X21 - V*grad_smooth_abs*(X11-X21);
X22_dt = df2(2) - theta3*X22 - smooth_abs*X22 - V*grad_smooth_abs*(X12-X22);
X23_dt = df2(3) - theta3*X23 - smooth_abs*X23 - V*grad_smooth_abs*(X13-X23);
X24_dt = df2(4) - theta3*X24 - smooth_abs*X24 - V*grad_smooth_abs*(X14-X24);

dXdt = [dTdt;dVdt;X11_dt;X21_dt;X12_dt;X22_dt;X13_dt;X23_dt;X14_dt;X24_dt];
end
function fsign_x = fsign(x,X)
% % M should be row vector of size 1xk
 column_count = size(X,2);
 first_non_zero = 0; 
 if sum(x) ~= 0
  first_non_zero = sign(x); 
 else    
  for k=1:column_count
       S = sum(X(:,k));
       if S~=0
          first_non_zero =  sign(X(:,k));
          break
       end
  end
 end
 fsign_x = first_non_zero';
end