classdef ExtendedTrace
    
%     StateVector as follows X = [x Vx ax y Vy ay h d T]'
    
    properties (Access = public)
        posts;                  % координаты постов
        posts_number;           % количество постов
        sigma_n;                % шум наблюдений
        X;                      % вектор состояния
        ToA;                    % измерения, сек
        ranges;                 % дальности
        mode;                   % режим полета 1-V-const, 2-V-const/acc-const/V-const, 3-V-const/V_ug-const/V-const
        start_time_str;         % начальное время
        start_time_sec;         % начальное время, сек
        T_period = 1;        % период излучения (минимальный), сек
        pause = 0;              % возможный интервал паузы начала излучения, сек
        N_periods = 0;        % число периодов, которые могут быть пропущены
        sigma_T = 0;           % ско девиации периодов, сек
        max_coord;              % пределы по рабочей зоне
        max_V;                  % пределы по скорости
        max_acc;                % пределы по ускорениям
        V = 1;                % модуль скорости, м/с
        acc=1;                  % ускорение
        V_ug = 0.01;           % угловой скорости при повороте
        angle;
        N_measurements = 200;   % число измерений
        c;                      % скорость света, м/c
        hei;                    % высота ИРИ, м
        Measurements;           % массив измерений на выдачу
        StateVector;            % вектор состояния на выдачу
        fligth_time = 1000;       % время полета, сек
        sigma_h;            % формирующий шум по высоте
        R_mode4;            % радиус полета по кругу
        w_mode4;        % угловая скорость полета по кругу
        angle_mode4;    % угол при полете по кругу
    end
    
    methods (Access = public)
        function obj = ExtendedTrace(config, start_time,mode) % конструктор
           obj.mode=mode;
           obj.start_time_str = datestr(datetime(start_time),'dd.MM.yyyy HH:mm:ss.FFF');
           obj.start_time_sec = start_time(4)*3600 + start_time(5)*60 + start_time(6);
           obj.posts = config.posts;
           obj.posts_number = config.posts_number;
           obj.hei = config.hei;
           obj.c = config.c;
           obj.max_coord = config.max_coord;
           obj.max_V = config.max_V;
           obj.max_acc = config.max_acc;
           obj.sigma_n = config.sigma_t;
           obj.sigma_h = config.sigma_h;
           if obj.mode == 4
               obj.R_mode4 = 50e3;
               
           end
%            obj.Measurements = zeros(obj.posts_number,obj.N_measurements);
           obj = Initialization(obj);
%            obj.StateVector = zeros(length(obj.X),obj.N_measurements);
           
        end
        
        function obj = Go(obj)
            stop_time = obj.start_time_sec + obj.fligth_time;
            i = 1;
            time=obj.X(9);
            switch obj.mode
                case 1
                    
                while obj.X(9) < stop_time
                obj = obj.UpdateTrace_mode1;
                obj.StateVector(:,i) = obj.X;
                obj.Measurements(:,i) = obj.ToA;
                if abs(obj.X(1)) > (obj.max_coord) || abs(obj.X(4)) > (obj.max_coord)
                    break
                end
                i = i + 1;
                end
                
                case 2
                    
                    while obj.X(9)< stop_time
                if obj.X(9) < time + round(2*(stop_time-time)/3) && obj.X(9) > time + round((stop_time-time)/3)
%                 if obj.X(6) < time + round(1*(stop_time-time)/3) 
                    obj = obj.UpdateTrace_mode2;
                else
                    obj = obj.UpdateTrace_mode1;
                 end
                obj.StateVector(:,i) = obj.X;
                obj.Measurements(:,i) = obj.ToA;
                if abs(obj.X(1)) > (obj.max_coord) || abs(obj.X(4)) > (obj.max_coord)
                    break
                end
                i = i + 1;
                    end
                    
                case 3
                while obj.X(9)< stop_time
                if obj.X(9) < time+round(2*(stop_time-time)/3) && obj.X(9) > time + round((stop_time-time)/3)
                    obj = obj.UpdateTrace_mode3;
                else
                    obj = obj.UpdateTrace_mode1;
                 end
                obj.StateVector(:,i) = obj.X;
                obj.Measurements(:,i) = obj.ToA;
                if abs(obj.X(1)) > (obj.max_coord) || abs(obj.X(4)) > (obj.max_coord)
                    break
                end
                i = i + 1;
                end
                
                case 4
                
                while obj.X(9) < stop_time
                obj = obj.UpdateTrace_mode4;
                obj.StateVector(:,i) = obj.X;
                obj.Measurements(:,i) = obj.ToA;
                if abs(obj.X(1)) > (obj.max_coord) || abs(obj.X(4)) > (obj.max_coord)
                    break
                end
                i = i + 1;
                end
                    
                
                
                
            end
            
            
        end
        
        function show(obj)
            plot(obj.posts(1,:),obj.posts(2,:),'^')
            hold on
            plot(obj.StateVector(1,:),obj.StateVector(4,:),'.-')
            grid on
            axis([-obj.max_coord obj.max_coord -obj.max_coord obj.max_coord])
            daspect([1 1 1])
        end
        
    end
    
    methods (Access = public)
        
        function obj = Initialization(obj)
            obj.X(1,1) = randi([-obj.max_coord obj.max_coord]);
            obj.X(4,1) = randi([-obj.max_coord obj.max_coord]);
            obj.angle = randi([0 2*314])/100;
            obj.X(2,1) = obj.V * cos(obj.angle);
            obj.X(5,1) = obj.V * sin(obj.angle);
            obj.X(3,1) = 0;
            obj.X(6,1) = 0;
            obj.X(7,1) = obj.hei;
            obj.X(8,1) = obj.T_period;
            obj.X(9,1) = obj.start_time_sec + randi([0 obj.pause]);
            
            if obj.mode == 4
                obj.w_mode4 = obj.V/obj.R_mode4;
                obj.angle_mode4 = randi([0 2*314])/100;
                obj.X(1,1) = obj.R_mode4 * cos(obj.angle_mode4);
                obj.X(4,1) = obj.R_mode4 * sin(obj.angle_mode4);
                obj.X(2,1) = 0;
                obj.X(5,1) = 0;
                obj.X(3,1) = 0;
                obj.X(6,1) = 0;
                obj.X(7,1) = obj.hei;
                obj.X(8,1) = obj.T_period;
                obj.X(9,1) = obj.start_time_sec + randi([0 obj.pause]);
            end
        end
        
        function obj = UpdateTrace_mode1(obj)
            T = obj.X(8,1);
            F = [1 T T^2/2 0 0 0 0 0 0;
                 0 1 T 0 0 0 0 0 0;
                 0 0 1 0 0 0 0 0 0;
                 0 0 0 1 T T^2/2 0 0 0;
                 0 0 0 0 1 T 0 0 0;
                 0 0 0 0 0 1 0 0 0;
                 0 0 0 0 0 0 1 0 0;
                 0 0 0 0 0 0 0 0 0;
                 0 0 0 0 0 0 0 1 1];
             
            obj.X = F*obj.X;
            obj.X(3) = 0;
            obj.X(6) = 0;
            obj.X(8) = obj.T_period * (1 + randi([0 obj.N_periods])) + normrnd(0,obj.sigma_T);
            obj.X(7) = obj.X(7) + T*normrnd(0,obj.sigma_h);
            
            for i = 1:obj.posts_number
                obj.ranges(i,1) = norm(obj.posts(:,i) - [obj.X(1); obj.X(4); obj.X(7)]);
                obj.ToA(i,1) = obj.X(9) + obj.ranges(i,1)/obj.c + normrnd(0,obj.sigma_n);
            end
        end
        
        function obj = UpdateTrace_mode2(obj)
            T = obj.X(8,1);
            F = [1 T T^2/2 0 0 0 0 0 0;
                 0 1 T 0 0 0 0 0 0;
                 0 0 1 0 0 0 0 0 0;
                 0 0 0 1 T T^2/2 0 0 0;
                 0 0 0 0 1 T 0 0 0;
                 0 0 0 0 0 1 0 0 0;
                 0 0 0 0 0 0 1 0 0;
                 0 0 0 0 0 0 0 0 0;
                 0 0 0 0 0 0 0 1 1];
            obj.X = F*obj.X;
            obj.X(8) = obj.T_period * (1 + randi([0 obj.N_periods])) + normrnd(0,obj.sigma_T);
            obj.X(7) = obj.X(7) + T*normrnd(0,obj.sigma_h);
            
            %adding acceleration
            obj.X(3) = obj.acc*cos(obj.angle);
            obj.X(6) = obj.acc*sin(obj.angle);
            
            
            for i = 1:obj.posts_number
                obj.ranges(i,1) = norm(obj.posts(:,i) - [obj.X(1); obj.X(4); obj.X(7)]);
                obj.ToA(i,1) = obj.X(9) + obj.ranges(i,1)/obj.c + normrnd(0,obj.sigma_n);
            end
        end
        
        function obj = UpdateTrace_mode3(obj)
            T = obj.X(8,1);
            %increasing angle
            obj.angle=obj.angle+obj.V_ug;

             obj.X(1,1) = obj.X(1,1)+obj.X(2,1)*T;
             obj.X(2,1) = obj.V * cos(obj.angle);
             obj.X(3,1) = 0;
             obj.X(4,1) = obj.X(4,1)+obj.X(5,1)*T;
             obj.X(5,1) = obj.V * sin(obj.angle);
             obj.X(6,1) = 0;
             obj.X(7,1) = obj.X(7) + T*normrnd(0,obj.sigma_h);
             obj.X(9,1) = obj.X(8,1)+ obj.X(9,1);
             obj.X(8,1) = obj.T_period * (1 + randi([0 obj.N_periods])) + normrnd(0,obj.sigma_T);
             
            for i = 1:obj.posts_number
                obj.ranges(i,1) = norm(obj.posts(:,i) - [obj.X(1); obj.X(4); obj.X(7)]);
                obj.ToA(i,1) = obj.X(9) + obj.ranges(i,1)/obj.c + normrnd(0,obj.sigma_n);
            end
        end
        
        
        function obj = UpdateTrace_mode4(obj)
            T = obj.X(8,1);
            %increasing angle
            obj.angle_mode4 = obj.angle_mode4 + obj.w_mode4 * T;

             obj.X(1,1) = obj.R_mode4 * cos(obj.angle_mode4);
             obj.X(4,1) = obj.R_mode4 * sin(obj.angle_mode4);
             obj.X(3,1) = 0;
             obj.X(2,1) = 0;
             obj.X(5,1) = 0;
             obj.X(6,1) = 0;
             obj.X(7,1) = obj.X(7) + T*normrnd(0,obj.sigma_h);
             obj.X(9,1) = obj.X(8,1)+ obj.X(9,1);
             obj.X(8,1) = obj.T_period * (1 + randi([0 obj.N_periods])) + normrnd(0,obj.sigma_T);
             
            for i = 1:obj.posts_number
                obj.ranges(i,1) = norm(obj.posts(:,i) - [obj.X(1); obj.X(4); obj.X(7)]);
                obj.ToA(i,1) = obj.X(9) + obj.ranges(i,1)/obj.c + normrnd(0,obj.sigma_n);
            end
        end
        
        
    end
    
end

