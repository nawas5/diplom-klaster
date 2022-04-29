function [uwb, gyro, acc, mag, euler] = log_reader_1_0(filename)

if nargin == 0
    [file, path] = uigetfile('*.*');
    filename = fullfile(path,file);  
end

    f = fopen(filename);

    uwb = [];
    gyro = [];
    acc = [];
    mag = [];
    euler = [];
    
    
    while feof(f)==0 
        s=fgetl(f);
        
        if contains(s,"GYR")
            k = size(gyro,2) + 1;
            S = split(s);
            gyro(:,k) = [str2num(S{2,1})/1e9; str2num(S{3,1}); str2num(S{4,1}); str2num(S{5,1})];
        end
        
        if contains(s,"EYR")
            k = size(euler,2) + 1;
            S = split(s);
            euler(:,k) = [str2num(S{2,1})/1e9; str2num(S{3,1}); str2num(S{4,1}); str2num(S{5,1})];
        end
        
        if contains(s,"ACC")
            k = size(acc,2) + 1;
            S = split(s);
            acc(:,k) = [str2num(S{2,1})/1e9; str2num(S{3,1}); str2num(S{4,1}); str2num(S{5,1})];
        end
        
        if contains(s,"MAG")
            k = size(mag,2) + 1;
            S = split(s);
            mag(:,k) = [str2num(S{2,1})/1e9; str2num(S{3,1}); str2num(S{4,1}); str2num(S{5,1})];
        end
        
        if contains(s,"server")
            k = size(uwb,2) + 1;
            S = split(s);
            uwb(:,k) = [str2num(S{2,1})/1e9; str2num(S{4,1}); str2num(S{5,1}); str2num(S{6,1})];
        end
        
    end
    
    if size(uwb,2) > 0
        T_min = min([acc(1,1) gyro(1,1) mag(1,1) uwb(1,1) euler(1,1)]);
    else
        T_min = min([acc(1,1) gyro(1,1) mag(1,1) euler(1,1)]);
    end
    
    acc(1,:) = acc(1,:) - T_min;
    gyro(1,:) = gyro(1,:) - T_min;
    mag(1,:) = mag(1,:) - T_min;
    euler(1,:) = euler(1,:) - T_min;
    
    if size(uwb,2) > 0
        uwb(1,:) = uwb(1,:) - T_min; 
        subplot(411)
        title("UWB")
        plot(uwb(1,:),uwb(2:3,:),'linewidth',2)
        grid on
        xlabel('t, sec')
        ylabel('meter')
    end
    
    subplot(412)
    title("Euler")
    plot(euler(1,:),euler(2:4,:),'linewidth',2)
    grid on
    xlabel('t, sec')
    ylabel('rad')
    
    subplot(413)
    title("acc")
    plot(acc(1,:),acc(2:4,:),'linewidth',2)
    grid on
    xlabel('t, sec')
    ylabel('хз')
    
    subplot(414)
    title("gyro")
    plot(gyro(1,:),gyro(2:4,:),'linewidth',2)
    grid on
    xlabel('t, sec')
    ylabel('хз')
    
    fclose(f);
    
end

