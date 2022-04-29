classdef EKF_2
    
    properties (Access = public) % атрибуты класса
        X; % вектор состояния
        Dx; % матрица дисперсий ошибок фильтрации
        sigma_n; % СКО шума наблюдения ToA отметок
        sigma_phi;
        sigma_ksi_v; % СКО формирующего шума по скорости
        sigma_ksi_phi;
        c; % скорость света
        D_ksi; % матрица формирующих шумов
    end
    
    methods (Access = public)
        function obj = EKF_2(X0, y, config) % конструктор
           obj.X = X0;
           obj.Dx = eye(length(X0));
           obj.sigma_n = config.sigma_n; % м
           obj.sigma_phi = config.sigma_phi; % м
           obj.sigma_ksi_v = config.sigma_ksi_v; % м/с
           obj.sigma_ksi_phi = config.sigma_ksi_phi; % м/с
           obj.D_ksi = eye(2);
           obj.D_ksi(1,1) = obj.sigma_ksi_v^2;
           obj.D_ksi(2,2) = obj.sigma_ksi_phi^2;
        end
 
        
        function obj = Update(obj, y, dt, config) % загрузка новых измерений = обновление фильтра
            
            X_prev = obj.X;
            Dx = obj.Dx;
             
            G = [0 0; 0 0; dt 0; 0 0; 0 dt];
            
            D_n = config.sigma_n^2 * eye(size(y,1));
            
            F = obj.make_F(X_prev, config, dt);
            dF = obj.make_dF(X_prev, config, dt);
            
            X_ext = F;

            D_x_ext = dF * Dx * dF' + G * obj.D_ksi * G';
            dS = obj.make_dS(X_ext, config);
            S = dS*X_ext;
            K = D_x_ext * dS' * inv(dS*D_x_ext*dS' + D_n);
            Dx = D_x_ext - K * dS * D_x_ext;
            X_prev = X_ext + K*(y - S);
            
            obj.X = X_prev;
            obj.Dx = Dx;
            
        end
        
        function [ S ] = make_S( obj, X, config )
            for i = 1:size(config.posts,2)
                S(i,1) = sqrt((X(1,1) - config.posts(1,i))^2 + (X(3,1) - config.posts(2,i))^2 + (config.hei - config.posts(3,i))^2) + X(5,1);
            end
        end
        
        function [ F ] = make_F( obj, X, config, dt )
            F(1,1) = X(1,1) + X(3,1) * cos(X(4,1))*dt;
            F(2,1) = X(2,1) + X(3,1) * sin(X(4,1))*dt;
            F(3,1) = X(3,1);
            F(4,1) = X(4,1) + X(5,1)*dt;
            F(5,1) = X(5,1);
        end
        
        function [ dF ] = make_dF( obj, X, config, dt )
            dF(1,:) = [1 0 cos(X(4,1))*dt -X(3,1) * sin(X(4,1))*dt 0];
            dF(2,:) = [0 1 sin(X(4,1))*dt X(3,1) * cos(X(4,1))*dt 0];
            dF(3,:) = [0 0 1 0 0];
            dF(4,:) = [0 0 0 1 dt];
            dF(5,:) = [0 0 0 0 1];
                
        end


        function [ dS ] = make_dS( obj, X, config )
            dS = [1 0 0 0 0;
                  0 1 0 0 0;];
        end

        function [ delta ] = make_delta( obj, y, X, config )
            X = obj.X;
            for i = 1:length(y)
                x = config.posts(:,i);
                delta(i,1) = y(i)/(1 + (X(2)*(X(1) - x(1)) + X(4)*(X(3) - x(2)))/(config.c*sqrt((X(1) - x(1))^2 + (X(3) - x(2))^2 + (config.hei - x(3))^2)));
            end
        end
        
    end
    
    
    
end



