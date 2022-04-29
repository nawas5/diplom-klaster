close all
clear all
clc
config = Config();
config.sigma_ksi = 10;

[uwb, gyro, acc, mag, euler] = log_reader_1_0();

%only for vertical
gyro_new = interp1(gyro(1,:),gyro(4,:),uwb(1,:));
config.sigma_n = 0.1;
config.sigma_phi = 0.02;
config.sigma_ksi_v = 0.4;
config.sigma_ksi_phi = 0.1;
config.T = 0.1;
y = [uwb(2,:);uwb(3,:);gyro_new];
%y = y(:,80:180);
X_true = [y(1,1);y(2,1);1.5;atan2((9.04 - 1.91),(18.33 - 4.2));0];

X1(:,1) = X_true(:,1);
Fil1 = EKF_2(X1(:,1), y(1:2,1), config);
Fil1.Dx = [0.0194226457970021,0.000626771108315620,0.00185336046539954,0.0196281539455231,0.0125972524041403;0.000626771108315621,0.0150882599905291,-0.0135512284324082,0.00289161121067866,0.00173121499410728;0.00185336046539954,-0.0135512284324082,0.0274808198289386,-0.000502840435743853,-0.000630486841949936;0.0196281539455231,0.00289161121067866,-0.000502840435743853,0.0338676136820128,0.0297902233229618;0.0125972524041403,0.00173121499410728,-0.000630486841949936,0.0297902233229618,0.0398697835143881];
X2(:,1) = X_true(:,1);
Fil2 = EKF_3(X2(:,1), y(:,1), config);
Fil2.Dx = [0.0151800844461052,4.93275007504195e-06,0.00124987254400736,0.0111338301507617,0.00497369040407111;4.93275007504188e-06,0.0150183918521298,-0.0136244025988374,0.000987562519275644,0.000392877668340553;0.00124987254400736,-0.0136244025988374,0.0274400544325498,-7.91283192913706e-06,-8.88421339111869e-05;0.0111338301507617,0.000987562519275644,-7.91283192913627e-06,0.0160616350133531,0.0127434095644531;0.00497369040407111,0.000392877668340553,-8.88421339111869e-05,0.0127434095644531,0.0216707105804582];
D_x(:,1) = Fil2.Dx(1,1);
D_y(:,1) = Fil2.Dx(2,2);
D_V(:,1) = Fil2.Dx(3,3);
D_alp(:,1) = Fil2.Dx(4,4);


for i = 1:length(y)
    if i > 1
       
        
        Fil1 = Fil1.Update(y(1:2,i), config.T, config);
        X1(:,i) = Fil1.X;
        
        Fil2 = Fil2.Update(y(:,i), config.T, config);
        X2(:,i) = Fil2.X;
        M_x(:,i) = sqrt(Fil2.Dx(1,1));
        M_y(:,i) = sqrt(Fil2.Dx(2,2));
        M_V(:,i) = sqrt(Fil2.Dx(3,3));
        M_alp(:,i) = sqrt(Fil2.Dx(4,4));
        M_omeg(:,i) = sqrt(Fil2.Dx(4,4));
    end
end

%% Reference
N1 = 63;
N2 = 175;
N_all = length(uwb(1,:));

Vx_true_0 = (18.33 - 4.2)/(24.02-12.02);
Vy_true_0 = (9.04 - 1.91)/(24.02-12.02);

% for k = 1:N + 100
for k = 1:N_all    
    if k >= 1 && k <= N1
        x_true(k) = 3.0802;
        y_true(k) = 1.3319;
        %V_true(k) = 0;
        Vx_true(k) = 0;
        Vy_true(k) = 0;
        alpha_true(k) = atan2((9.04 - 1.91),(18.33 - 4.2));
        omega_true(k) = 0;
        if k == 1
            T = 0.1;
        else
            T = uwb(1,k) - uwb(1,k-1);
        end
        if k == N1
           Vx_true(k) = Vx_true_0;
           Vy_true(k) = Vy_true_0;            
        end
        
    elseif k >= N2
        x_true(k) = x_true(k - 1);
        y_true(k) = y_true(k - 1);
        %V_true(k) = 0;
        Vx_true(k) = 0;
        Vy_true(k) = 0;
        alpha_true(k) = alpha_true(k - 1);
        omega_true(k) = 0;
      
    else
        T = uwb(1,k) - uwb(1,k-1);
        %V_true(k) = V_true(k - 1) + T*sigma_ksi_V*randn;
        Vx_true(k) = Vx_true(k - 1);
        Vy_true(k) = Vy_true(k - 1);
%         x_true(k) = x_true(k - 1) + cos(alpha_true(k-1))*V_true(k-1)*T;
%         y_true(k) = y_true(k - 1) + sin(alpha_true(k-1))*V_true(k-1)*T;
        x_true(k) = x_true(k - 1) + Vx_true(k-1)*T;
        y_true(k) = y_true(k - 1) + Vy_true(k-1)*T;
        %alpha_true(k) = alpha_true(k-1)*(1-mu*T) + mu*T*sigma_ksi_alpha*randn;
        %alpha_true(k) = alpha_true(k-1) + T*sigma_ksi_alpha*randn;
        %omega_true(k) = (alpha_true(k) - alpha_true(k - 1))/T;
        alpha_true(k) = alpha_true(k-1);
        omega_true(k) = 0;
        
        end
 end


%% Static before moving
x_error_stat1 = X2(1,1:N1) - x_true(1:N1);
y_error_stat1 = X2(2,1:N1) - y_true(1:N1);
alpha_error_stat1 = X2(4,1:N1) - alpha_true(1:N1);

mean_x_stat1 = mean(x_error_stat1);
mean_y_stat1 = mean(y_error_stat1);
mean_alpha_stat1 = mean(alpha_error_stat1);

std_x_stat1 = std(x_error_stat1);
std_y_stat1 = std(y_error_stat1);
std_alpha_stat1 = std(alpha_error_stat1);

drmse_stat1 = 2*sqrt(std_x_stat1^2 + std_y_stat1^2);

%% Moving
x_error_mov = X2(1,N1:N2) - x_true(N1:N2);
y_error_mov = X2(2,N1:N2) - y_true(N1:N2);
alpha_error_mov = X2(4,N1:N2) - alpha_true(N1:N2);

mean_x_mov = mean(x_error_mov);
mean_y_mov = mean(y_error_mov);
mean_alpha_mov = mean(alpha_error_mov);

std_x_mov = std(x_error_mov);
std_y_mov = std(y_error_mov);
std_alpha_mov = std(alpha_error_mov);

drmse_mov = 2*sqrt(std_x_mov^2 + std_y_mov^2);
%% Static after moving
x_error_stat2 = X2(1,N2:200) - x_true(N2:200);
y_error_stat2 = X2(2,N2:200) - y_true(N2:200);
alpha_error_stat2 = X2(4,N2:200) - alpha_true(N2:200);

mean_x_stat2 = mean(x_error_stat2);
mean_y_stat2 = mean(y_error_stat2);
mean_alpha_stat2 = mean(alpha_error_stat2);

std_x_stat2 = std(x_error_stat2);
std_y_stat2 = std(y_error_stat2);
std_alpha_stat2 = std(alpha_error_stat2);


velocity_error = X2(3,:) - sqrt(Vx_true.^2 + Vy_true.^2);

drmse_st2 = 2*sqrt(std_x_stat2^2 + std_y_stat2^2);

D_coord = 0.0939;
%%
figure(2)
% title("x")
plot(uwb(1,:),x_true(1,:),'k--','linewidth',2)
hold on
plot(uwb(1,:),y(1,:),'b','linewidth',2)
grid on
% plot(X1(1,:),'b','linewidth',2)
plot(uwb(1,:),X2(1,:),'r','linewidth',2)
ylabel('$$X, m$$','Interpreter','LaTeX', 'FontSize',24)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
%title('Зависимость истинной координаты Х и ее оценки от времени','FontSize',18)
legend({"Референс", "Измерение", "Оценка фильтра СШП+ИНС"},'FontSize',10)
set(gca, 'FontSize',17)

figure(3)
plot(uwb(1,:),y_true(1,:),'k--','linewidth',2)
hold on
plot(uwb(1,:), y(2,:),'b','linewidth',2)
% plot(X1(2,:),'b','linewidth',2)
plot(uwb(1,:),X2(2,:),'r','linewidth',2)
grid on
ylabel('$$Y, m$$','Interpreter','LaTeX', 'FontSize',24)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
legend({"Референс", "Измерение", "Оценка фильтра СШП+ИНС"},'FontSize',10)
set(gca, 'FontSize',17)
%title('Зависимость истинной координаты Y и ее оценки от времени','FontSize',18)

% figure(3)
% title("V")
% grid on
% hold on
% plot(X1(3,:),'b','linewidth',2)
% plot(X2(3,:),'r','linewidth',2)
% legend("-","+")

% figure(4)
% title("phi")
% 
% grid on
% hold on
% plot(X1(4,:),'b','linewidth',2)
% plot(X2(4,:),'r','linewidth',2)
% legend("-","+")

figure(4)
plot(uwb(1,:),omega_true(1,:),'k--','linewidth',2)
hold on
plot(uwb(1,:),y(3,:),'b','linewidth',2)
grid on
% plot(X1(5,:),'b','linewidth',2)
plot(uwb(1,:),X2(5,:),'r','linewidth',2)
xlim([uwb(1,1) uwb(1,end)])
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \dot{\alpha}, rad/s $$', 'Interpreter','LaTeX', 'FontSize',24)
legend({"Референс", "Измерение", "Оценка фильтра СШП+ИНС"},'FontSize',10)
%title('оценка производной курса от времени','FontSize',18)
set(gca, 'FontSize',17)

figure(5)
plot(uwb(1,:),alpha_true(1,:),'k--','linewidth',2)
grid on
hold on
% plot(X1(5,:),'b','linewidth',2)
plot(uwb(1,:),X2(4,:),'r','linewidth',2)
xlim([uwb(1,1) uwb(1,end)])
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \alpha, rad $$', 'Interpreter','LaTeX', 'FontSize',24)
legend({"Референс", "Оценка фильтра СШП+ИНС"},'FontSize',10)
%title('оценка курса от времени','FontSize',18)
set(gca, 'FontSize',17)

figure(6)
plot(x_true(1,:),y_true(1,:),'k--','linewidth',2)
hold on
plot(y(1,:),y(2,:),'b','linewidth',2)
grid on 
% plot(X1(1,:),X1(2,:),'b','linewidth',1)
plot(X2(1,:),X2(2,:),'r','linewidth',2)
%title('Истинная траектория и еe оценка','FontSize',18)
xlabel('$$X, m$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$Y, m$$','Interpreter','LaTeX', 'FontSize',24)
legend({"Референс", "Измерение", "Оценка фильтра СШП+ИНС"},'FontSize',10)
%title('оценка производной курса от времени','FontSize',18)
set(gca, 'FontSize',17)

%%
figure(7)
%plot(y(1,:) - uwb(2,:),'k')
%plot(X1(1,:) - uwb(2,:),'b','linewidth',2)
plot(uwb(1,:),X2(1,:) - uwb(2,:),'r','linewidth',2)
hold on
plot(uwb(1,:),3*D_coord,'m.','linewidth',8)
plot(uwb(1,:),-3*D_coord,'m.','linewidth',8)

plot(uwb(1,:),4*M_x,'b--','linewidth',2)
plot(uwb(1,:),-4*M_x,'b--','linewidth',2)

grid on
xlim([uwb(1,1) uwb(1,end)])
title('Ошибка оценки координаты Х','FontSize',18)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \widehat{X} - X, m $$','Interpreter','LaTeX', 'FontSize',24)
set(gca, 'FontSize',17)


figure(8)
title("y")
% plot(y(2,:) - uwb(3,:),'k')
%plot(X1(2,:) - uwb(3,:),'b','linewidth',2)
plot(uwb(1,:),X2(2,:) - uwb(3,:),'r','linewidth',2)
grid on
hold on
plot(uwb(1,:),4*M_y,'b--','linewidth',2)
plot(uwb(1,:),-4*M_y,'b--','linewidth',2)

plot(uwb(1,:),3*0.0939,'m.','linewidth',4)
plot(uwb(1,:),-3*0.0939,'m.','linewidth',4)

xlim([uwb(1,1) uwb(1,end)])

title('Ошибка оценки координаты Y','FontSize',18)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \widehat{Y} - Y, m $$','Interpreter','LaTeX', 'FontSize',24)
set(gca, 'FontSize',17)

figure(9)
title("V")
%plot(X1(3,:) - X_true(3,:),'b','linewidth',2)
%plot(uwb(1,:),X2(3,:) - X_true(3,:),'r','linewidth',2)
plot(uwb(1,:),5*M_V,'b--','linewidth',2)
grid on
hold on
plot(uwb(1,:),-5*M_V,'b--','linewidth',2)

plot(uwb(1,:),velocity_error,'r-','linewidth',2)

plot(uwb(1,:),3*0.1798,'m.','linewidth',4)
plot(uwb(1,:),-3*0.1798,'m.','linewidth',4)

xlim([uwb(1,1) uwb(1,end)])
title('Ошибка оценкии скорости','FontSize',18)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \widehat{V} - V, m/s $$','Interpreter','LaTeX','FontSize',24)
set(gca, 'FontSize',17)

figure(10)
title("phi")

%plot(X1(4,:) - X_true(4,:),'b','linewidth',2)
plot(uwb(1,:),X2(4,:) - X_true(4,:),'r','linewidth',2)
grid on
hold on
plot(uwb(1,:),3*M_alp,'b--','linewidth',2)
plot(uwb(1,:),-3*M_alp,'b--','linewidth',2)

plot(uwb(1,:),3*0.096,'m.','linewidth',4)
plot(uwb(1,:),-3*0.096,'m.','linewidth',4)

xlim([uwb(1,1) uwb(1,end)])
title('Ошибка оценки угла курса','FontSize',18)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \widehat{\alpha} - \alpha, rad $$', 'Interpreter','LaTeX','FontSize',24)
set(gca, 'FontSize',17)

figure(11)
% title("dphi")
%plot(y(3,:) - X_true(5,:),'k*')
plot(uwb(1,:),3*M_omeg,'b--','linewidth',2)
grid on
hold on
plot(uwb(1,:),-3*M_omeg,'b--','linewidth',2)
%plot(X1(5,:) - X_true(5,:),'b','linewidth',2)
plot(uwb(1,:),X2(5,:) - X_true(5,:),'r','linewidth',2)
xlim([uwb(1,1) uwb(1,end)])
title('Ошибка оценки скорости угла курса','FontSize',18)
xlabel('$$t, s$$','Interpreter','LaTeX', 'FontSize',24)
ylabel('$$ \widehat{\dot{\alpha}} - \dot{\alpha}, rad/s $$', 'Interpreter','LaTeX','FontSize',24)
set(gca, 'FontSize',17)

