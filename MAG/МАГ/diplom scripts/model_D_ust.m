close all
clear all
clc

N = 5000;

sigma_ksi_V = 0.3;
sigma_ksi_alpha = 0.02;

sigma_n_x = 0.2;
sigma_n_y = 0.2;
sigma_n_omega = 0.01;

% sigma_ksi_V = 0.5;
% sigma_ksi_alpha = 0.1;
% 
% sigma_n_x = 0.1;
% sigma_n_y = 0.1;
% sigma_n_omega = 0.1;
Vx_true_0 = (18.33 - 4.2)/(24.02-12.02);
Vy_true_0 = (9.04 - 1.91)/(24.02-12.02);
V_true_0 = sqrt(Vx_true_0^2 + Vy_true_0^2);

T = 0.1;
% mu = 0.01/T;
% nu = 0.02/T;
t = 5:T:30;
% N = (24.02-12.02)/T;
% N1 = (12.02-5)/T;
% N2  = (30-24.02)/T;
% N_all = (30-5)/T;

N1 = fix(12.02/T);
N2 = fix(24.02/T);
N_all = fix(30/T);

% for k = 1:N + 100
for k = 1:N_all    
    if k >= 1 && k <= N1
        x_true(k) = 4.2;
        y_true(k) = 1.91;
        V_true(k) = 0;
        Vx_true(k) = 0;
        Vy_true(k) = 0;
        alpha_true(k) = atan2((9.04 - 1.91),(18.33 - 4.2));
        omega_true(k) = 0;
        if k == 1
            t(k) = 0;
        else
            t(k) = t(k-1) + T;
        end
        if k == N1
           Vx_true(k) = Vx_true_0;
           Vy_true(k) = Vy_true_0;
           V_true(k) = V_true_0;
        end
        
    elseif k >= N2
        x_true(k) = x_true(k - 1);
        y_true(k) = y_true(k - 1);
        %V_true(k) = 0;
        Vx_true(k) = 0;
        Vy_true(k) = 0;
        V_true(k) = 0;
        alpha_true(k) = alpha_true(k - 1);
        omega_true(k) = 0;
        t(k) = t(k-1) + T;
    else
        %V_true(k) = V_true(k - 1) + T*sigma_ksi_V*randn;
        Vx_true(k) = Vx_true(k - 1);
        Vy_true(k) = Vy_true(k - 1);
         V_true(k) =  V_true(k-1);
%         x_true(k) = x_true(k - 1) + cos(alpha_true(k-1))*V_true(k-1)*T;
%         y_true(k) = y_true(k - 1) + sin(alpha_true(k-1))*V_true(k-1)*T;
        x_true(k) = x_true(k - 1) + Vx_true(k-1)*T;
        y_true(k) = y_true(k - 1) + Vy_true(k-1)*T;
                %alpha_true(k) = alpha_true(k-1)*(1-mu*T) + mu*T*sigma_ksi_alpha*randn;
        %alpha_true(k) = alpha_true(k-1) + T*sigma_ksi_alpha*randn;
        %omega_true(k) = (alpha_true(k) - alpha_true(k - 1))/T;
        alpha_true(k) = alpha_true(k-1);
        omega_true(k) = 0;
        t(k) = t(k-1) + T;
        end
    %V_true(k) = sqrt(Vx_true(k)^2 + Vy_true(k)^2);
    y_meas(:,k) = [x_true(k) + sigma_n_x*randn; y_true(k) + sigma_n_y*randn; omega_true(k) + sigma_n_omega*randn];
end


% figure(1)
% plot(V_true)
% title('Измерения скорости','FontSize',18)
% grid on

% figure(2)
% title('Х изм, У изм, Х оценка, У оценка','FontSize',18)
% plot(x_true)
% hold on
% plot(y_true)
% plot(y_meas(1,:))
% plot(y_meas(2,:))
% grid on
% 
% figure(3)
% title('Курс измер, скорость курса измер, Оценка курса','FontSize',18)
% plot(alpha_true)
% hold on
% plot(omega_true)
% plot(y_meas(3,:))
% grid on

%x_est = [3.5; 1; 0.1; 0.5; 0.1];
x_est = [4.2; 1.91; 0; atan2((9.04 - 1.91),(18.33 - 4.2)); 0];
%atan2((9.04 - 1.91),(18.33 - 4.2))
Dx_est = eye(5);
Dx_est = [2.90561546557470e-05,1.41089044394337e-05,7.13209782458592e-05,-4.16965147678531e-08,-4.51819501562224e-20;1.41089044394337e-05,8.22202230580691e-06,3.59947696376184e-05,8.26098113612398e-08,8.95268561245708e-20;7.13209782458592e-05,3.59947696376184e-05,0.000452838260593956,-6.15569761402546e-12,-1.75041578815690e-24;-4.16965147678531e-08,8.26098113612398e-08,-6.15569761411971e-12,2.22054326705367e-08,3.99680319548197e-13;-4.51819501562147e-20,8.95268561245526e-20,-1.75041578115088e-24,3.99680319548269e-13,9.99600319680544e-09];
%sigma_ksi_V_f = 1e0*sqrt(sigma_ksi_Vx^2 + sigma_ksi_Vy^2);
%sigma_ksi_alpha_f = 1e0*sqrt(sigma_ksi_Vx^2 + sigma_ksi_Vy^2);

D_ksi = [sigma_ksi_V^2         0 ;
                0      (sigma_ksi_alpha)^2];

Dn = [sigma_n_x^2 0 0; 0 sigma_n_y^2 0; 0 0 sigma_n_omega^2];            
%Dn = [sigma_n_x^2 0 0; 0 sigma_n_y^2 0; 0 0 sig];     

%mu = 0.01/T;
            
G = [0 0; 0 0; T 0; 0 0; 0 1];            
C = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 0 0 1];
% 
% omega_meas = y_meas(3,:);
% y_meas(3,:) = [];
 
for k = 2:N_all%N + 100
    x_ext(1,k) = x_est(1,k-1) + x_est(3,k-1)*cos(x_est(4,k-1))*T;
    x_ext(2,k) = x_est(2,k-1) + x_est(3,k-1)*sin(x_est(4,k-1))*T;
    x_ext(3,k) = x_est(3,k-1);
    x_ext(4,k) = x_est(4,k-1) + x_est(5, k-1)*T;
    x_ext(5,k) = x_est(5,k-1);
    
    dFdx = [1 0 cos(x_ext(4,k))*T  -x_ext(3,k)*sin(x_ext(4,k))*T 0;
            0 1 sin(x_ext(4,k))*T   x_ext(3,k)*cos(x_ext(4,k))*T 0;
            0 0       1                       0           0;
            0 0       0                       1           T;
            0 0       0                       0           1];
    
    Dx_ext = dFdx*Dx_est*dFdx' + G*D_ksi*G';
%     Dx_est = (Dx_ext^-1 + C'*Dn*C)^-1;
    K = Dx_ext*C'/(C*Dx_ext*C'+Dn);
    x_est(:,k) = x_ext(:,k) + K*(y_meas(:,k) - C*x_ext(:,k));
    Dx_est = (eye(5)-K*C)*Dx_ext;
    Mx(1,k) = sqrt(Dx_est(1,1));
    My(1,k) = sqrt(Dx_est(2,2));
    MV(1,k) = sqrt(Dx_est(3,3));
    Malp(1,k) = sqrt(Dx_est(4,4));
    Momeg(1,k) = sqrt(Dx_est(5,5));
end

for q = 1:length(omega_true)
   if q == 1
      omega_tr_int(q) = omega_true(q)*T + alpha_true(1);
      omega_ms_int(q) = y_meas(3,q)*T + alpha_true(1);
   else
      omega_tr_int(q) = omega_tr_int(q - 1) + omega_true(q)*T;
      omega_ms_int(q) = omega_ms_int(q - 1) + y_meas(3,q)*T;
   end
end    

%% Static before moving
x_error_stat1 = x_est(1,1:N1) - x_true(1:N1);
y_error_stat1 = x_est(2,1:N1) - y_true(1:N1);
alpha_error_stat1 = x_est(4,1:N1) - alpha_true(1:N1);

mean_x_stat1 = mean(x_error_stat1);
mean_y_stat1 = mean(y_error_stat1);
mean_alpha_stat1 = mean(alpha_error_stat1);

std_x_stat1 = std(x_error_stat1);
std_y_stat1 = std(y_error_stat1);
std_alpha_stat1 = std(alpha_error_stat1);

drmse_stat1 = 2*sqrt(std_x_stat1^2 + std_y_stat1^2);
%% Moving
x_error_mov = x_est(1,N1:N2) - x_true(N1:N2);
y_error_mov = x_est(2,N1:N2) - y_true(N1:N2);
alpha_error_mov = x_est(4,N1:N2) - alpha_true(N1:N2);

mean_x_mov = mean(x_error_mov);
mean_y_mov = mean(y_error_mov);
mean_alpha_mov = mean(alpha_error_mov);

std_x_mov = std(x_error_mov);
std_y_mov = std(y_error_mov);
std_alpha_mov = std(alpha_error_mov);


mean_X_meas_mov = mean(x_true(N1:N2));
std_X_meas_mov = std(x_true(N1:N2));

mean_Y_meas_mov = mean(y_true(N1:N2));
std_Y_meas_mov = std(y_true(N1:N2));

mean_alpha_meas_mov = mean(alpha_true(N1:N2));
std_alpha_mov = std(alpha_true(N1:N2));

drmse_mov = 2*sqrt(std_X_meas_mov^2 + std_Y_meas_mov^2);
%% Static after moving
x_error_stat2 = x_est(1,N2:N_all) - x_true(N2:N_all);
y_error_stat2 = x_est(2,N2:N_all) - y_true(N2:N_all);
alpha_error_stat2 = x_est(4,N2:N_all) - alpha_true(N2:N_all);

mean_x_stat2 = mean(x_error_stat2);
mean_y_stat2 = mean(y_error_stat2);
mean_alpha_stat2 = mean(alpha_error_stat2);

std_x_stat2 = std(x_error_stat2);
std_y_stat2 = std(y_error_stat2);
std_alpha_stat2 = std(alpha_error_stat2);


mean_X_meas_stat_2 = mean(x_true(N2:N_all));
std_X_meas_stat_2 = std(x_true(N2:N_all));

mean_Y_meas_stat_2 = mean(y_true(N2:N_all));
std_Y_meas_stat_2 = std(y_true(N2:N_all));

mean_alpha_meas_stat_2 = mean(alpha_true(N2:N_all));
std_alpha_stat_2 = std(alpha_true(N2:N_all));

drmse_stat2 = 2*sqrt(std_X_meas_stat_2^2 + std_Y_meas_stat_2^2);
%%
velocity_error = x_est(3,:) - sqrt(Vx_true.^2 + Vy_true.^2);

%%
%график траектории
figure(3)
plot(x_true, y_true,'b-', 'linewidth',2)
hold on
%plot(y_meas(1,:), y_meas(2,:),'k-', 'linewidth',1)
grid on
xlabel('$$X, m$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$Y, m$$','Interpreter','LaTeX', 'FontSize',24)
legend({"Истинная траектория","Смоделированная траектория"},'FontSize',10)
set(gca, 'FontSize',17)

%Угол курса
figure(4)
plot(t, alpha_true(1,:),'b-', 'linewidth',2)
hold on
plot(t, x_est(4,:), 'r-','linewidth',2)
grid on
ylabel('$$ \alpha, rad $$', 'Interpreter','LaTeX','FontSize',24)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
%title('Истинный угол курса и его оценка','FontSize',18)
legend({"Истинный угол курса","Оценка фильтра СШП+ИНС"},'FontSize',10)
set(gca, 'FontSize',17)

%производная угла курса
figure(7)
plot(t, omega_true(1,:),'b-', 'linewidth',2)
hold on
plot(t, y_meas(3,:),'k-', 'linewidth',0.5)
plot(t, x_est(5,:), 'r-','linewidth',2)
grid on
ylabel('$$ \dot{\alpha}, rad/s $$', 'Interpreter','LaTeX','FontSize',24)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
%title('Истинная скорость угола курса, смоделированное измерение и оценка','FontSize',18)
legend({"Истинная скорость угла курса", "Cмоделированное измерение","Оценка фильтра СШП+ИНС"},'FontSize',10)
set(gca, 'FontSize',17)

%Координата Х
figure(5)
% title('','FontSize',18)
plot(t, x_true(1,:),'b-', 'linewidth',2)
hold on
grid on
plot(t, y_meas(1,:),'k-', 'linewidth',0.5)
plot(t, x_est(1,:),'r-', 'linewidth',2)
ylabel('$$X, m$$','Interpreter','LaTeX', 'FontSize',24)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
%title('Зависимость истинной координаты Х, смоделированного измерения и оценки от времени','FontSize',18)
legend({"Истинная координата Х", "Cмоделированное измерение", "Оценка фильтра СШП+ИНС"},'FontSize',10)
set(gca, 'FontSize',17)

%Координата У
figure(6)
plot(t, y_true(1,:),'b-', 'linewidth',2)
hold on
grid on
plot(t, y_meas(2,:),'k-', 'linewidth',0.5)
plot(t, x_est(2,:),'r-', 'linewidth',2)
ylabel('$$Y, m$$','Interpreter','LaTeX', 'FontSize',24)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
%title('Зависимость истинной координаты Y, смоделированного измерения и оценки от времени','FontSize',18)
legend({"Истинная координата Y","Cмоделированное измерение", "Оценка фильтра СШП+ИНС"},'FontSize',10)
set(gca, 'FontSize',17)

%Оценка траектории
figure(9)
plot(x_true, y_true,'b-', 'linewidth',2)
hold on
plot(y_meas(1,:), y_meas(2,:),'k-', 'linewidth',0.5)
plot(x_est(1,:),x_est(2,:),'r-', 'linewidth',2)
grid on
%title('Истинная траектория, ее модель и еe оценка','FontSize',18)
xlabel('$$X, m$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$Y, m$$','Interpreter','LaTeX', 'FontSize',24)
legend({"Истинная траектория","Модель радиоизмерений", "Оценка фильтра СШП+ИНС"},'FontSize',10)
set(gca, 'FontSize',17)


%Ошибка оценки угла курса
figure(14)
plot(t(1:N1),alpha_error_stat1,'b','linewidth',2)
hold on
plot(t,3*Malp,'r--','linewidth',2)
plot(t,3*0.096,'m.')
plot(t(N1:N2),alpha_error_mov,'b','linewidth',2)
plot(t(N2:end),alpha_error_stat2,'b','linewidth',2)
plot(t,-3*Malp,'r--','linewidth',2)
plot(t,-3*0.096,'m.','linewidth',2)
grid on
title('Ошибка оценки угла курса','FontSize',18)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \widehat{\alpha} - \alpha, rad $$', 'Interpreter','LaTeX','FontSize',24)
set(gca, 'FontSize',17)
%legend("Ошибка оценки угла курса","Предельная граница ошибки фильтрации, рассчитанная по оценке матрицы дисперсии","Предельная граница ошибки фильтрации, рассчитанная аналитически",'FontSize',10)

%ошибка оценкии координаты Х
figure(15)
plot(t(1:N1),x_error_stat1,'b-','linewidth',2)
hold on
plot(t(N1:N2),x_error_mov,'b-','linewidth',2)
plot(t(N2:end),x_error_stat2,'b-','linewidth',2)
plot(t,3*Mx,'r--','linewidth',2)
plot(t,-3*Mx,'r--','linewidth',2)
plot(t,3*0.0939,'m.','linewidth',2)
plot(t,-3*0.0939,'m.','linewidth',2)
grid on
title('Ошибка оценкии координаты Х','FontSize',18)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \widehat{X} - X, m $$','Interpreter','LaTeX', 'FontSize',24)
set(gca, 'FontSize',17)
%legend("Ошибка оценки угла курса","Предельная граница ошибки фильтрации, рассчитанная по оценке матрицы дисперсии","Предельная граница ошибки фильтрации, рассчитанная аналитически",'FontSize',10)

%ошибка оценкии координаты Y
figure(16)
plot(t(1:N1),y_error_stat1,'b-','linewidth',2)
hold on
plot(t(N1:N2),y_error_mov,'b-','linewidth',2)
plot(t(N2:end),y_error_stat2,'b-','linewidth',2)
plot(t,3*My,'r--','linewidth',2)
plot(t,-3*My,'r--','linewidth',2)
plot(t,3*0.0939,'m.','linewidth',2)
plot(t,-3*0.0939,'m.','linewidth',2)
grid on
title('Ошибка оценкии координаты Y','FontSize',18)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \widehat{Y} - Y, m $$','Interpreter','LaTeX', 'FontSize',24)
set(gca, 'FontSize',17)
% legend("Смоделированное измерение","Оценка фильтра СШП+ИНС")

%Ошибка оценкии скорости
figure(17)
plot(t, velocity_error,'b-','linewidth',2)
hold on
plot(t,3*MV,'r--','linewidth',2)
plot(t,-3*MV,'r--','linewidth',2)
plot(t,3*0.1798,'m.','linewidth',2)
plot(t,-3*0.1798,'m.','linewidth',2)
grid on
title('Ошибка оценкии скорости','FontSize',18)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \widehat{V} - V, m/s $$','Interpreter','LaTeX','FontSize',24)
set(gca, 'FontSize',17)
% legend("Смоделированное измерение","Оценка фильтра СШП+ИНС")

% Оценка скорости
figure(18)
plot(t, V_true(1,:),'b-', 'linewidth',2)
hold on
plot(t, x_est(3,:), 'r-','linewidth',2)
grid on
ylabel('$$ V, m/s $$', 'Interpreter','LaTeX','FontSize',24)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
%title('Истинный угол курса и его оценка','FontSize',18)
legend({"Истинная скорость потребителя","Оценка фильтра СШП+ИНС"},'FontSize',10)
set(gca, 'FontSize',17)

