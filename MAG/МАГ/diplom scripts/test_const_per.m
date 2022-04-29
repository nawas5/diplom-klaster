close all 
clear all
clc

config = Config();
config.sigma_ksi = 10;


trace = ExtendedTrace(config, clock, 3);
trace.V = 1;
trace.N_periods = 0;
trace = trace.Initialization();
trace.X(1,1) = 0;
trace.X(2,1) = 1;
trace.X(4,1) = 0;
trace.X(5,1) = 0;
trace = trace.Go;

StateVector = trace.StateVector;
for i = 1:length(StateVector)
    X_true(1,i) = StateVector(1,i);
    X_true(2,i) = StateVector(4,i);
    X_true(3,i) = norm(StateVector([2 5],i));
    X_true(4,i) = atan2(StateVector(5,i),StateVector(2,i));
    if i > 1
        X_true(5,i-1) = (X_true(4,i) - X_true(4,i-1))/config.T;
    end
end

y(1,:) = X_true(1,:) + normrnd(0, config.sigma_n, [1 length(X_true)]);
y(2,:) = X_true(2,:) + normrnd(0, config.sigma_n, [1 length(X_true)]);
y(3,:) = X_true(5,:) + normrnd(0, config.sigma_phi, [1 length(X_true)]);

X1(:,1) = X_true(:,1);
Fil1 = EKF_2(X1(:,1), y(1:2,1), config);

X2(:,1) = X_true(:,1);
Fil2 = EKF_3(X2(:,1), y(:,1), config);
%Fil2.Dx = [0.00706180925246408,-4.20359534236868e-06,-0.000225402413394873,0.000180157820434087,5.62005599106602e-11;-4.20359534236902e-06,0.00706363047701007,-0.000179983314592340,-0.000225614976795541,-6.98159815859314e-11;-0.000225402413394873,-0.000179983314592340,2.45518573342991e-05,-1.99513246705830e-10,1.25368553905550e-15;0.000180157820434088,-0.000225614976795541,-1.99513246707884e-10,2.45343768610111e-05,9.80454848922937e-09;5.62005599106594e-11,-6.98159815859314e-11,1.25368553906216e-15,9.80454848922916e-09,9.90195135927764e-07];




for i = 1:length(y)
    if i > 1
       
        
        Fil1 = Fil1.Update(y(1:2,i), config.T, config);
        X1(:,i) = Fil1.X;
        
        Fil2 = Fil2.Update(y(:,i), config.T, config);
        X2(:,i) = Fil2.X;
        
       
        
    end
end

STD_uwb_coord_x = sqrt((1/(length(X1)-1))*sum((X1(1,:)-X_true(1,:)).^2));
STD_uwb_coord_y = sqrt((1/(length(X1)-1))*sum((X1(2,:)-X_true(2,:)).^2));
STD_uwb_sum = sqrt(STD_uwb_coord_x^2 + STD_uwb_coord_y^2);
M_uwb_X = mean(X1(1,:));
M_uwb_Y = mean(X1(2,:));

STD_uwb_ins_coord_x = sqrt((1/(length(X2)-1))*sum((X2(1,:)-X_true(1,:)).^2));
STD_uwb_ins_coord_y = sqrt((1/(length(X2)-1))*sum((X2(2,:)-X_true(2,:)).^2));
STD_uwb_ins_sum = sqrt(STD_uwb_ins_coord_x^2 + STD_uwb_ins_coord_y^2);
M_uwb_ins_X = mean(X2(1,:));
M_uwb_ins_Y = mean(X2(2,:));

Gain = sqrt((M_uwb_X^2+M_uwb_Y^2+STD_uwb_coord_x^2+STD_uwb_coord_y^2)/(M_uwb_ins_X^2+M_uwb_ins_Y^2+STD_uwb_ins_coord_x^2+STD_uwb_ins_coord_y^2));


figure(1)
title("Координата Х")
plot(y(1,:) - X_true(1,:),'k')
grid on
hold on
plot(X1(1,:) - X_true(1,:),'b','linewidth',2)
plot(X2(1,:) - X_true(1,:),'r','linewidth',2)
xlabel("t, s")
ylabel("X, m")
legend("Measurements","UWB","UWB+INS")

figure(2)
title("y")
plot(y(2,:) - X_true(2,:),'k')
grid on
hold on
plot(X1(2,:) - X_true(2,:),'b','linewidth',2)
plot(X2(2,:) - X_true(2,:),'r','linewidth',2)
xlabel("t, s")
ylabel("Y, m")
legend("Measurements","UWB","UWB+INS")

figure(3)
title("V")
grid on
hold on
plot(X1(3,:) - X_true(3,:),'b','linewidth',2)
plot(X2(3,:) - X_true(3,:),'r','linewidth',2)
xlabel("t, s")
ylabel("V, m/s")
legend("UWB","UWB+INS")

figure(4)
title("phi")
grid on
hold on
plot(X1(4,:) - X_true(4,:),'b','linewidth',2)
plot(X2(4,:) - X_true(4,:),'r','linewidth',2)
xlabel("t, s")
ylabel('$$ \alpha, rad/s $$', 'Interpreter','LaTeX')
legend("UWB","UWB+INS")

figure(5)
title("dphi")
plot(y(3,:) - X_true(5,:),'k*')
grid on
hold on
plot(X1(5,:) - X_true(5,:),'b','linewidth',2)
plot(X2(5,:) - X_true(5,:),'r','linewidth',2)
xlabel("t, s")
ylabel('$$ \dot{\alpha}, rad/s $$', 'Interpreter','LaTeX')
legend("Measurements","UWB","UWB+INS")

figure(6)
plot(X_true(1,:),X_true(2,:),'k','linewidth',2)
grid on 
hold on
plot(X1(1,:),X1(2,:),'b','linewidth',1)
plot(X2(1,:),X2(2,:),'r','linewidth',1)
xlabel("Х, m")
ylabel("Y, m")
legend("Measurements","UWB","UWB+INS")

figure(7)
plot(X_true(1,:),X_true(2,:),'k','linewidth',2)
xlabel("Х, m")
ylabel("Y, m")
grid on 


