clear all
close all
clc

[uwb, gyro, acc, mag, euler] = log_reader_1_0('121243.txt');

% A = [1 3 4 7; 5 2 3 1; 9 8 3 1];
% B = [2 6 8 9; 9 8 3 1; 5 2 3 1];
% 
% C1 = [A(1,:) B(1,:)];
% [D,I] = sort(C1);
% 
% C = [A(2:3,:) B(2:3,:)];
% C = C(:,I);
% C = [D; C];

%C = unique(A,B);
% last_ind_A = 1;
% last_ind_B = 1;
% maxT = 9;
% C = [];
% k = 1;
% while 1
%     if isempty(C)
%         if A(1,last_ind_A) < B(1,last_ind_B)
%            C = A(:,last_ind_A);
%            if last_ind_A ~= size(A,2)
%              last_ind_A = last_ind_A + 1;
%            end
%         else
%            C = B(:,last_ind_B);
%            if last_ind_B ~= size(B,2)
%               last_ind_B = last_ind_B + 1;
%            end
%         end
%     else
%         if A(1,last_ind_A) < B(1,last_ind_B)
%            C = [C A(:,last_ind_A)];
%            if last_ind_A ~= size(A,2)
%               last_ind_A = last_ind_A + 1;
%            end
%         else
%            C = [C B(:,last_ind_B)];
%            if last_ind_B ~= size(B,2)
%               last_ind_B = last_ind_B + 1;
%            end
%         end
%     end
%     if C(1,end) == maxT
%        break
%     end
%     k = k + 1;
% end

C1 = [uwb(1,:) gyro(1,:)];
[D,I] = sort(C1);

DataMeas = [uwb(2:4,:) gyro(2:4,:)];
DataMeas = DataMeas(:,I);
DataMeas = [D; DataMeas];

sigma_ksi_V = 0.5;
sigma_ksi_alpha = 0.1;

sigma_n_xy = sqrt(std(uwb(2,9:55))^2 + std(uwb(3,9:55))^2);

sigma_n_x = sigma_n_xy;
sigma_n_y = sigma_n_xy;

sigma_n_wxyz = sqrt(std(gyro(2,10:430))^2 + std(gyro(3,10:430))^2 + std(gyro(4,10:430))^2);
bias_omega = mean(gyro(4,50:300));
sigma_n_omega = sigma_n_wxyz;
%sigma_n_omega = 10;

T = 1;
% x_est = [5.05; 2.86; 0; 0];
% x_est_stor = [];
% Dx_est = eye(4);
% D_ksi = [sigma_ksi_V^2         0;
%                 0      (sigma_ksi_alpha)^2];          
% Dn = [sigma_n_x^2 0; 0 sigma_n_y^2];        
% G = [0 0; 0 0; T 0; 0 T];            
% C = [1 0 0 0;
%      0 1 0 0];
sigma_ksi_Valpha = 0.1;
x_est = [5.22; 2.71; 0; 0.8543; 0];
x_est_stor = [];
Dx_est = eye(5);
D_ksi = [sigma_ksi_V^2         0  ;
                0      (sigma_ksi_Valpha)^2];          
        
G = [0 0; 0 0; T 0; 0 0; 0 T];  

% omega_meas = y_meas(3,:);
% y_meas(3,:) = [];
 
% for k = 2:length(DataMeas)    
%     if DataMeas(4,k) ~= 0 
%         T = DataMeas(1,k) - DataMeas(1,k-1);
%         x_ext(1,1) = x_est(1) + x_est(3)*cos(x_est(4))*T;
%     	x_ext(2,1) = x_est(2) + x_est(3)*sin(x_est(4))*T;
%         x_ext(3,1) = x_est(3);
%         x_ext(4,1) = x_est(4) + (DataMeas(3,k) - bias_omega)*T;
%     
%         dFdx = [1 0 cos(x_ext(4))*T  -x_ext(3)*sin(x_ext(4))*T;
%                 0 1 sin(x_ext(4))*T   x_ext(3)*cos(x_ext(4))*T;
%                 0 0       1                       0           ;
%                 0 0       0                       1];
%     
%         Dx_ext = dFdx*Dx_est*dFdx' + G*D_ksi*G';
%         
%         x_est = x_ext;
%         Dx_est = Dx_ext;
%         x_est_stor = [x_est_stor x_est];
%     else
%         y_meas = DataMeas(2:3,k);
%         x_ext = x_est;
%         Dx_ext = Dx_est;
% %       Dx_est = (Dx_ext^-1 + C'*Dn*C)^-1;
%         K = Dx_ext*C'/(C*Dx_ext*C'+Dn);
%         x_est = x_ext + K*(y_meas - C*x_ext);
%         Dx_est = (eye(4)-K*C)*Dx_ext;
%         x_est_stor = [x_est_stor x_est];
%         %M(1,k) = sqrt(Dx_est(2,2));
%     end
% end

for k = 2:length(DataMeas)      
        T = DataMeas(1,k) - DataMeas(1,k-1);
        x_ext(1,1) = x_est(1) + x_est(3)*cos(x_est(4))*T;
    	x_ext(2,1) = x_est(2) + x_est(3)*sin(x_est(4))*T;
        x_ext(3,1) = x_est(3);
        x_ext(4,1) = x_est(4) + x_est(5)*T;
        x_ext(5,1) = x_est(5);
    
        dFdx = [1 0 cos(x_ext(4))*T  -x_ext(3)*sin(x_ext(4))*T 0;
                0 1 sin(x_ext(4))*T   x_ext(3)*cos(x_ext(4))*T 0;
                0 0       1                       0            0;
                0 0       0                       1            T;
                0 0       0                       0            1];
    
        Dx_ext = dFdx*Dx_est*dFdx' + G*D_ksi*G';
    if DataMeas(4,k) == 0   
       y_meas = DataMeas(2:3,k);
       C = [1 0 0 0 0;
            0 1 0 0 0];
       Dn = [sigma_n_x^2 0;
              0 sigma_n_y^2];
    else
       y_meas = DataMeas(4,k);
       C = [0 0 0 0 1]; 
       Dn = sigma_n_omega^2;
    end
       K = Dx_ext*C'/(C*Dx_ext*C'+Dn);
       x_est = x_ext + K*(y_meas - C*x_ext);
       Dx_est = (eye(5)-K*C)*Dx_ext;
       x_est_stor = [x_est_stor x_est];
       %M(1,k) = sqrt(Dx_est(2,2));
  
end

alpha_meas = 0.8543;
for k = 2:length(gyro)
    Ts = gyro(1,k) -  gyro(1,k-1);
    alpha_meas(k) = alpha_meas(k - 1) + (gyro(4,k)-bias_omega)*Ts;
end

figure
plot(gyro(1,:),alpha_meas)
grid on

figure
plot(gyro(1,:),gyro(4,:)-bias_omega)
hold on
plot(gyro(1,:),gyro(2,:))
plot(gyro(1,:),gyro(3,:))
grid on

figure
plot(x_est_stor(1,:))
hold on
plot(x_est_stor(2,:))
grid on

figure
plot(x_est_stor(1,:),x_est_stor(2,:))
hold on
plot(uwb(2,:),uwb(3,:))
grid on

figure
plot(DataMeas(1,2:end),x_est_stor(4,:))
hold on
plot(gyro(1,:),alpha_meas)
grid on
% figure(4)
% 
% plot(x_est(4,:))
% hold on
% plot(alpha_true(1,:))
% grid on
% title('Измер скорости курса, оценка скорости курса','FontSize',18)
% 
% figure(5)
% % title('','FontSize',18)
% plot(x_est(1,:))
% hold on
% plot(x_true(1,:))
% grid on
% title('Измер Х, оценка Х','FontSize',18)
% 
% figure(6)
% plot(x_est(2,:))
% hold on
% plot(y_true(1,:))
% grid on
% title('Измер У, оценка У','FontSize',18)
% 
% figure(7)
% plot((y_true(1,:)- x_est(2,:)))
% hold on
% plot(3*M(1,:))
