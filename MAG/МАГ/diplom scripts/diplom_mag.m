clear all
close all
clc

N = 5000;

sigma_ksi_V = 0.3;
sigma_ksi_alpha = 0.1;

sigma_n_x = 0.2;
sigma_n_y = 0.2;
sigma_n_omega = 0.1;

% sigma_ksi_V = 0.5;
% sigma_ksi_alpha = 0.1;
% 
% sigma_n_x = 0.1;
% sigma_n_y = 0.1;
% sigma_n_omega = 0.1;
Vx_true_0 = (18.33 - 4.2)/(24.02-12.02);
Vy_true_0 = (9.04 - 1.91)/(24.02-12.02);

T = 0.1;
% mu = 0.01/T;
% nu = 0.02/T;
N = (24.02-12.02)/T;
for k = 1:N
    if k == 1
        x_true(k) = 4.2;
        y_true(k) = 1.91;
        %V_true(k) = 0;
        Vx_true(k) = Vx_true_0;
        Vy_true(k) = Vy_true_0;
        alpha_true(k) = atan2((9.04 - 1.91),(18.33 - 4.2));
        omega_true(k) = 0;
        t(k) = 0;
    else
        %V_true(k) = V_true(k - 1) + T*sigma_ksi_V*randn;
        Vx_true(k) = Vx_true(k - 1);
        Vy_true(k) = Vy_true(k - 1);
%         x_true(k) = x_true(k - 1) + cos(alpha_true(k-1))*V_true(k-1)*T;
%         y_true(k) = y_true(k - 1) + sin(alpha_true(k-1))*V_true(k-1)*T;
        x_true(k) = x_true(k - 1) + Vx_true(k-1)*T;
        y_true(k) = y_true(k - 1) + Vy_true(k-1)*T;
        %alpha_true(k) = alpha_true(k-1)*(1-mu*T) + mu*T*sigma_ksi_alpha*randn;
        alpha_true(k) = alpha_true(k-1) + T*sigma_ksi_alpha*randn;
        omega_true(k) = (alpha_true(k) - alpha_true(k - 1))/T;
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

x_est = [4.2; 1.91; 0; atan2((9.04 - 1.91),(18.33 - 4.2)); 0];
Dx_est = eye(5);

%sigma_ksi_V_f = 1e0*sqrt(sigma_ksi_Vx^2 + sigma_ksi_Vy^2);
%sigma_ksi_alpha_f = 1e0*sqrt(sigma_ksi_Vx^2 + sigma_ksi_Vy^2);

D_ksi = [sigma_ksi_V^2         0 ;
                0      (sigma_ksi_alpha)^2];

Dn = [sigma_n_x^2 0 0; 0 sigma_n_y^2 0; 0 0 sigma_n_omega^2];            
%Dn = [sigma_n_x^2 0 0; 0 sigma_n_y^2 0; 0 0 sig];     

%mu = 0.01/T;
            
G = [0 0; 0 0; T 0; 0 0; 0 T];            
C = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 0 0 1];
% 
% omega_meas = y_meas(3,:);
% y_meas(3,:) = [];
 
for k = 2:N
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
    M(1,k) = sqrt(Dx_est(2,2));
end

figure(4)

plot(t, x_est(4,:), 'linewidth',2)
hold on
plot(t, alpha_true(1,:), 'linewidth',2)
grid on
ylabel('$$ \alpha, rad/s $$', 'Interpreter','LaTeX','FontSize',24)
xlabel("t, s", 'FontSize',24)
title('Истинный угол курса и его оценка','FontSize',18)
legend("Measurements","UWB+INS")

figure(5)
% title('','FontSize',18)
plot(t, x_est(1,:), 'linewidth',2)
hold on
plot(t, x_true(1,:), 'linewidth',2)
grid on
ylabel("Х, m", 'FontSize',24)
xlabel("t, s", 'FontSize',24)
title('Зависимость истинной координаты Х и ее оценки от времени','FontSize',18)

figure(6)
plot(t, x_est(2,:), 'linewidth',2)
hold on
plot(t, y_true(1,:), 'linewidth',2)
grid on
ylabel("Y, m", 'FontSize',24)
xlabel("t, s", 'FontSize',24)
title('Зависимость истинной координаты Y и ее оценки от времени','FontSize',18)

% figure(7)
% plot((y_true(1,:)- x_est(2,:)))
% hold on
% plot(3*M(1,:))
% grid on

% figure(8)
% plot(t, x_true)
% hold on
% plot(t, y_true)
% grid on

figure(9)
plot(x_true, y_true, 'linewidth',2)
hold on
plot(x_est(1,:),x_est(2,:), 'linewidth',2)
grid on
title('Истинная траектория и еe оценка','FontSize',18)
xlabel("Х, m", 'FontSize',24)
ylabel("Y, m", 'FontSize',24)
legend("Measurements","UWB+INS")
