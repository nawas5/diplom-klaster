clear all
clc
close all
% Начало замера тестов
tic
% -------------------------------------------------------------------------
% TODO Инициализация данных.
rng(10)
rr = 100; % Число реализаций
stat_matrix = nan(rr,3); % Матрица статистики всех реализаций
time = datestr(now,'ddmmyyhhMMss');
fid = fopen([time '_logs.txt'], 'wt'); % Открыли файл для записи
freq = [1.09 1.5 5.48 9.8 16]*1e9; % Несущая частота
dur = [50 100 500 20000 65000]*1e-9; % Длительность импульса
T_min = 2e-6;
for ll=1:rr % Прогоняем много реализаций
N_signal1 = 3; % Число импульсов в паттерне 1.
N_signal2 = 7; % Число импульсов в паттерне 2.
N_signal3 = 10; % Число импульсов в паттерне 3.
N_signal4 = 15;
N_signal5 = 20;
% Формирование 1-го паттерна.
prev_T = 0;
pattern_freq = freq(randi(length(freq)));
pattern_dur = dur(randi(length(dur)));
pattern_signal1 = nan(4,N_signal1); %4 параметра импульса.
for i = 1:N_signal1
 pattern_signal1(:,i) = make_impulse( pattern_freq, pattern_dur, T_min, prev_T );
 prev_T = pattern_signal1(1,i);
end
% Формирование 2-го паттерна.
prev_T = 0;
pattern_freq = freq(randi(length(freq)));
pattern_dur = dur(randi(length(dur)));
pattern_signal2 = nan(4,N_signal2); % 4 параметра импульса.
for i = 1:N_signal2
 pattern_signal2(:,i) = make_impulse( pattern_freq, pattern_dur, T_min, prev_T );
 prev_T = pattern_signal2(1,i);
end
%Формирование 3-го паттерна.
prev_T = 0;
pattern_freq = freq(randi(length(freq)));
pattern_dur = dur(randi(length(dur)));
pattern_signal3 = nan(4,N_signal3); % 4 параметра импульса.
for i = 1:N_signal3
 pattern_signal3(:,i) = make_impulse( pattern_freq, pattern_dur, T_min, prev_T );
 prev_T = pattern_signal3(1,i);
end
%Формирование 4-го паттерна.
prev_T = 0;
pattern_freq = freq(randi(length(freq)));
pattern_dur = dur(randi(length(dur)));
pattern_signal4 = nan(4,N_signal4); % 4 параметра импульса.
for i = 1:N_signal4
 pattern_signal4(:,i) = make_impulse( pattern_freq, pattern_dur, T_min, prev_T );
 prev_T = pattern_signal4(1,i);
end
%Формирование 5-го паттерна.
prev_T = 0;
pattern_freq = freq(randi(length(freq)));
pattern_dur = dur(randi(length(dur)));
pattern_signal5 = nan(4,N_signal5); % 4 параметра импульса.
for i = 1:N_signal5
 pattern_signal5(:,i) = make_impulse( pattern_freq, pattern_dur, T_min, prev_T );
 prev_T = pattern_signal5(1,i);
end
% prev_T = clock;
% prev_T = prev_T(6);
prev_T = 19.3680000000000;
N = 10000;
% Закидывание паттерна во всю выборку.
k1=0;
k2=0;
k3=0;
k4=0;
k5=0;
imp = nan(4,N);
positions1=0;
positions2=0;
positions3=0;
positions4=0;
positions5=0;
i = 1;
% TODO Добавление паттеринов в матрицу импульсов.
while i < N+1
 imp(:,i) = make_impulse( freq, dur, T_min, prev_T );
 if (randi(1000) == 300) && (i<(N-N_signal1))
 k1=k1+1;
 % i
 positions1(k1)= i;
 for j = 1:N_signal1
 imp(:,i) = pattern_signal1(:,j);
 imp(1,i) = imp(1,i) + prev_T;
 i = i + 1;
 end
 prev_T = imp(1,i-1); 
 
 elseif (randi(1000) == 100) && (i<(N-N_signal2))
 k2=k2+1;
 % i
 positions2(k2)= i;
 for j = 1:N_signal2
 imp(:,i) = pattern_signal2(:,j);
 imp(1,i) = imp(1,i) + prev_T;
 i = i + 1;
 end
 prev_T = imp(1,i-1);
 elseif (randi(1000) == 10) && (i<(N-N_signal3))
 k3=k3+1;
 % i
 positions3(k3)= i;
 for j = 1:N_signal3
 imp(:,i) = pattern_signal3(:,j);
 imp(1,i) = imp(1,i) + prev_T;
 i = i + 1;
 end
 prev_T = imp(1,i-1);
 elseif (randi(1000) == 10) && (i<(N-N_signal4))
 k4=k4+1;
 % i
 positions4(k4)= i;
 for j = 1:N_signal4
 imp(:,i) = pattern_signal4(:,j);
 imp(1,i) = imp(1,i) + prev_T;
 i = i + 1;
 end
 prev_T = imp(1,i-1);
 elseif (randi(1000) == 10) && (i<(N-N_signal5))
 k5=k5+1;
 % i
 positions5(k5)= i;
 for j = 1:N_signal5
 imp(:,i) = pattern_signal5(:,j);
 imp(1,i) = imp(1,i) + prev_T;
 i = i + 1;
 end
 prev_T = imp(1,i-1);
 else
 prev_T = imp(1,i);
 i = i + 1;
 end
 
end
% Истинные позиции паттернов.
true_positions.position1 = positions1';
true_positions.position2 = positions2';
true_positions.position3 = positions3';
true_positions.position4 = positions4';
true_positions.position5 = positions5';
vector_k = [k1 k2 k3 k4 k5];
% vector_k = [k1 k2 k3];
% Добавление шума.
imp(1,:) = imp(1,:) + normrnd(0, 50e-9,[1, N]);
imp(2,:) = imp(2,:) + normrnd(0, 1e3,[1, N]);
imp(3,:) = imp(3,:) + normrnd(0, 10e-9,[1, N]);
imp(4,1) = 0;
for i=2:N
 imp(4,i) = imp(1,i) - imp(1,i-1);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%TODO графики пространства признаков.
%пространство признаков до нормировки

% figure(1)
% scatter(imp(2,:)*1e-9,imp(3,:)*1e6,50,'filled')
% legend('Принятые импульсы');
% grid on
% title('Импульсы в пространстве признаков (f,\deltat)')
% xlabel('f,ГГц','FontSize',20,...
% 'FontWeight','bold')
% ylabel('\deltat,мкc','FontSize',20,...
% 'FontWeight','bold')
% set(gca, 'FontSize', 20)
% 
% 
% figure(2)
% scatter3(imp(1,:),imp(2,:)*1e-9,imp(3,:)*1e6,50,'filled')
% grid on
% legend('Принятые импульсы');
% title('Импульсы в пространстве признаков (ToA,f,\deltat)')
% ylabel('f,ГГц','FontSize',20,...
% 'FontWeight','bold')
% zlabel('\deltat,мкc','FontSize',20,...
% 'FontWeight','bold')
% xlabel('t,c','FontSize',20,...
% 'FontWeight','bold')
% set(gca, 'FontSize', 20)
% 
% figure(22)
% scatter3(imp(4,:)*1e6,imp(2,:)*1e-9,imp(3,:)*1e6,50,'filled')
% grid on
% legend('Принятые импульсы');
% title('Импульсы в пространстве признаков (ToA,f,\deltat)')
% ylabel('f,ГГц','FontSize',20,...
% 'FontWeight','bold')
% zlabel('\deltat,мкc','FontSize',20,...
% 'FontWeight','bold')
% xlabel('T,мкc','FontSize',20,...
% 'FontWeight','bold')
% set(gca, 'FontSize', 20)
% 
% 
% %построим распределения параметров импульсов
% figure(23)
% subplot(3,1,1)
% histogram(imp(2,:)*1e-9,'Normalization','probability')
% xlabel('f,ГГц','FontSize',20,...
% 'FontWeight','bold')
% subplot(3,1,2)
% histogram(imp(3,:)*1e6,'Normalization','probability')
% xlabel('\deltat,мкc','FontSize',20,...
% 'FontWeight','bold')
% subplot(3,1,3)
% histogram(imp(4,:)*1e6,'Normalization','probability')
% xlabel('T,мкc','FontSize',20,...
% 'FontWeight','bold')
% grid on

% Другой способ нормировки (z-стандартизация).
mean_z = mean(imp,2);
std_z = [std(imp(1,:)) std(imp(2,:)) std(imp(3,:)) std(imp(4,:)) ]';
imp_norm = (imp - mean_z)./std_z;
imp_norm = imp_norm';
%imp_norm = imp_norm(:,2:3);
imp_norm = imp_norm(:,2:4); %с учетом периода.
% -------------------------------------------------------------------------
% Проход окном.
answer_DBSCAN = struct();
% Окно в два паттерна.
%vector_N = [N_signal1 N_signal2];
% Окно в один паттерн.
vector_N = [N_signal1];
for p=1:length(vector_N)
 N_signal = vector_N(p);
 IMP = nan(1,N_signal*length(imp_norm(1,:)));
 for i = 1:(N-N_signal+1)
 iii = [];
 for j = 1:N_signal
 iii = [iii imp_norm((i-1) + j,:)];
 end
 IMP(i,:) = iii;
 end
 % -------------------------------------------------------------------------
 % -------------------------------------------------------------------------
 % Предобработка перед кластеризацией (удаление начального значения периода в
 % строках
 
 IMP(:,3) = [];
 
 %$crutch$
 minpts = min(vector_k); %число соседей (1. берет равный или больше числа feature+1) (2. число паттернов)
 kD = pdist2(IMP,IMP,'euc','Smallest',minpts); % The minpts smallest pairwise 

 
 kD(minpts+1,:) = sort(kD(end,:));
 kD(minpts+2,1) = 0;
 for z=2:length(kD)
 kD(minpts+2,z) = kD(minpts+1,z)-kD(minpts+1,z-1);
 end
 % -------------------------------------------------------------------------
 % -------------------------------------------------------------------------
 % Графики
 % figure(1)
 % subplot(2,1,1)
 % plot(kD(end-2,:))
 % ylabel('Euc','FontSize',20,...
 % 'FontWeight','bold')
 % xlabel('номер паттерна,k','FontSize',20,...
 % 'FontWeight','bold') 
 % grid on
 % subplot(2,1,2)
 % plot(kD(end-2,:),'.','MarkerSize',10)
 % ylabel('Euc','FontSize',20,...
 % 'FontWeight','bold')
 % xlabel('номер паттерна,k','FontSize',20,...
 % 'FontWeight','bold')
 % title('Взаимные расстояния до минимального числа соседей в кластере')
 % set(gca, 'FontSize', 20)
 % grid on
 
 % Графики
% figure(2)
% subplot(2,1,1)
% plot(kD(end-1,:))
% ylabel('Euc','FontSize',20,...
% 'FontWeight','bold')
% xlabel('номер паттерна,k','FontSize',20,...
% 'FontWeight','bold')
% title('Отсортированные взаимные расстояния до минимального числа соседей в кластере') 
% grid on
% set(gca, 'FontSize', 20)
% subplot(2,1,2)
% plot(kD(end,:))
% ylabel("Euc'",'FontSize',20,...
% 'FontWeight','bold')
% xlabel('номер паттерна,k','FontSize',20,...
% 'FontWeight','bold')
% title('Производная верхнего графика') 
% grid on
% set(gca, 'FontSize', 20)
 [val,idx] = max(kD(end,1:round(end/2)));
 % epsilon = 0.01; % Равный колену графика (чем меньше espsilon, тем больше вероятность появляения нескольких кластеров или выбросов)
 % epsilon = max(kD(end,1:(end-1)))/2; % end-1, так как end точно не подходит
 % TODO добавил / 5
 epsilon = kD(end-1,idx)/5; %$crutch$
 
 if ll == 17
 aaaaaaa = 5;
 ffff = aaaaaaa;
 end
 labels = dbscan(IMP,epsilon,minpts); %класстеризуем через DBSCAN
 
 % disp(['Число кластеров - ' num2str(max(labels))])
 curr_str = ['=== Номер реализации -- ', num2str(ll), ' ====']; 
 fprintf(fid,'%s\n',curr_str);
 curr_str = 'Кластеризация DBSCAN:'; 
 fprintf(fid,'%s\n',curr_str);
 for v=1:max(labels)
 % disp(['Позиции ' num2str(v) '-ого кластера - ' num2str(find(labels==v)')])
 evalc(['answer_DBSCAN.claster', num2str(p), num2str(v), ' = find(labels==v)']);
 curr_str = ['Позиции ' num2str(v) '-ого кластера - ' num2str(find(labels==v)')]; 
 fprintf(fid,'%s\n',curr_str);
 end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Сопоставление паттернов со сдвигами в отдельные массивы. 
fieldnames_ans = fieldnames(answer_DBSCAN);
data_set = [];
for v=1:length(fieldnames(answer_DBSCAN))
 work_field_name = fieldnames_ans(v);
 work_field=eval(strcat('answer_DBSCAN.',work_field_name{1,1}));
 data_set = [data_set work_field'];
end
data_set = sort(data_set);
post_proc1 = struct();
temp_v = 0;
index=1;
index_internal = 2;
temp_v(1) = data_set(1);
for i=2:length(data_set)
 left = data_set(i-1);
 right = data_set(i);
 if left==right-1
 temp_v(index_internal) = right;
 index_internal = index_internal + 1;
 if data_set(i)==data_set(end)
 evalc(['post_proc1.index', num2str(index), ' = temp_v']); 
 end
 else
 evalc(['post_proc1.index', num2str(index), ' = temp_v']);
 index = index + 1;
 index_internal = 2;
 temp_v = 0;
 temp_v(1) = data_set(i);
 if data_set(i)==data_set(end)
 evalc(['post_proc1.index', num2str(index), ' = temp_v']); 
 end
 end
end
% -------------------------------------------------------------------------
if ll==17
 aaa=5;
end
% -------------------------------------------------------------------------
% Новая пост-обработка.
new_post_proc = struct();
threshold_freq = 2 * 20*1e9 * 0.001 * 6; %$crutch$
threshold_dur = 2 * 0.1*1e-6; %$crutch$
fieldnames_ans = fieldnames(post_proc1);
count = 1;
for i=1:length(fieldnames(post_proc1))
 temp_post_proc = struct();
 work_field_name = fieldnames_ans(i);
 work_field = eval(strcat('post_proc1.',work_field_name{1,1}));
 
 if size(work_field, 2) ~= 1
 for ii = 1:size(work_field, 2)
 flag = true;
 fieldnames_ans_in_temp_post_proc = fieldnames(temp_post_proc);
 for iii = 1:numel(fieldnames_ans_in_temp_post_proc)
 work_field_name_new = fieldnames_ans_in_temp_post_proc(iii);
 left_ind = eval(strcat('temp_post_proc.',work_field_name_new{1,1}));
 right_ind = work_field(ii);
 
 left_freq = mean(imp(2, left_ind));
 right_freq = imp(2, right_ind);
 left_dur = mean(imp(3, left_ind));
 right_dur = imp(3, right_ind);
 if i == 4
 ffffff = 5;
 end
 if abs(left_freq - right_freq) < threshold_freq && abs(left_dur -right_dur) < threshold_dur
 str = strcat('temp_post_proc.', work_field_name_new{1,1});
temp = [left_ind, right_ind];
evalc([str ' = temp']);
flag = false;
break;
 end
 
 end
 if flag
 temp_v = work_field(ii);
 evalc(['temp_post_proc.index', num2str(count), ' = temp_v']); 
 count = count + 1;
 end
 
 end

 else
 temp_v = work_field(1);
 evalc(['temp_post_proc.index', num2str(count), ' = temp_v']); 
 count = count + 1;
 end
 fieldnames_ans_in_new_post_proc = fieldnames(temp_post_proc);
 for ii=1:length(fieldnames_ans_in_new_post_proc)
 work_field_name_temp_post_proc = fieldnames_ans_in_new_post_proc(ii);
 work_field_temp_post_proc = eval(strcat('temp_post_proc.',work_field_name_temp_post_proc{1,1}));
 evalc(['new_post_proc.', work_field_name_temp_post_proc{1,1}, ' = work_field_temp_post_proc']);
 end
end
post_proc1 = new_post_proc;
% -------------------------------------------------------------------------
%$crutch$ сортировка структуры по длине каждого поля
% -------------------------------------------------------------------------
post_proc2 = struct();
massiv_lengths = nan(length(fieldnames(post_proc1)),3);
fieldnames_ans = fieldnames(post_proc1);
for i=1:length(fieldnames(post_proc1))
 work_field_name = fieldnames_ans(i);
 work_field = eval(strcat('post_proc1.',work_field_name{1,1}));
 temp_var = length(work_field);
 massiv_lengths(i,1) = temp_var;
end
[sort_seq,ss_idx] = sort(massiv_lengths(:,1));
massiv_lengths(:,2) = sort_seq;
massiv_lengths(:,3) = ss_idx;
j=1;
% for i=drange(massiv_lengths(:,3))
for i=1:length(fieldnames(post_proc1))
 num = massiv_lengths(i,3);
 work_field_name = fieldnames_ans(num);
 work_field = eval(strcat('post_proc1.',work_field_name{1,1}));
 evalc(['post_proc2.index', num2str(j), ' = work_field']);
 j=j+1;
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Избавление от сдвигов.
post_proc3 = struct();
if ll==17
 aaa=5;
end
fieldnames_ans = fieldnames(post_proc2);
for v=1:length(fieldnames(post_proc2))
 work_field_name = fieldnames_ans(v);
 work_field = eval(strcat('post_proc2.',work_field_name{1,1}));
 len = work_field(length(work_field))-work_field(1);
 if isfield(post_proc3, ['index' num2str(len-1)]) %$crutch$
 temp=eval(strcat('post_proc3.',['index' num2str(len-1)]));
 temp = [temp work_field(1)];
 evalc(['post_proc3.index', num2str(len-1), ' = temp']);
 elseif isfield(post_proc3, ['index' num2str(len+1)]) %$crutch$
 temp=eval(strcat('post_proc3.',['index' num2str(len+1)]));
 temp = [temp work_field(1)];
 evalc(['post_proc3.index', num2str(len+1), ' = temp']);
 elseif ~isfield(post_proc3, ['index' num2str(len)])
 evalc(['post_proc3.index', num2str(len), ' = work_field(1)']);
 else 
 temp = eval(strcat('post_proc3.',['index' num2str(len)]));
 temp = [temp work_field(1)];
 evalc(['post_proc3.index', num2str(len), ' = temp']);
 end 
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Защита от паттернов одинаковой длинны.
% Откоментить, если не нужен нижележащий блок.
% answer_postproc = post_proc3;
% -------------------------------------------------------------------------
% Закоментить этот блок, если он не нужен.
count = 1;
answer_postproc = struct();
% TODO выбираем threshold.
threshold_freq = 2 * 20*1e9 * 0.001 * 6; %$crutch$
threshold_dur = 2 * 0.1*1e-6;
fieldnames_ans = fieldnames(post_proc3);
for v=1:length(fieldnames_ans)
 time_answer = struct();
 work_field_name = fieldnames_ans(v);
 
 work_field = eval(strcat('post_proc3.',work_field_name{1,1}));
 c = 1;
 for k=1:size(work_field,2)
 new_i = work_field(k);
% for new_i = drange(work_field)
 flag = true;
 if ~isempty(time_answer)
 fieldnames_ans_1 = fieldnames(time_answer);
 for i_in_time_answer = 1:length(fieldnames_ans_1)
 work_field_in_time_answer = eval(strcat('time_answer.', fieldnames_ans_1{i_in_time_answer}));
 left_freq = imp(2, new_i);
 left_dur = imp(3, new_i);
 index = strsplit(work_field_name{1,1}, 'x');
 index = str2double(index(2));
 pattern_freq_temp = imp(2, work_field_in_time_answer(1):(work_field_in_time_answer(1) + vector_N(1) + index-1));
 pattern_dur_temp = imp(3, work_field_in_time_answer(1):(work_field_in_time_answer(1) + vector_N(1) + index-1));
 
 right_freq = mean(pattern_freq_temp); 
 right_dur = mean(pattern_dur_temp);
 
 if abs(right_freq - left_freq) <= threshold_freq && abs(right_dur - left_dur) <= threshold_dur
 temp = [work_field_in_time_answer new_i];
evalc([strcat('time_answer.', fieldnames_ans_1{i_in_time_answer}) ' = temp']);
 flag = false;
break;
 end
 end
 end
 if flag 
 str = strcat('time_answer.', work_field_name{1,1}, '_', num2str(count));
 temp = new_i;
 evalc([str ' = temp']);
 count = count + 1;
 end 
 c = c + 1;
 end
 fieldnames_ans_1 = fieldnames(time_answer);
 count = 1;
 for i_in_time_answer = 1:length(fieldnames_ans_1)
 work_field_name = fieldnames_ans_1(i_in_time_answer); 
 work_field = eval(strcat('time_answer.',work_field_name{1,1}));
 str = strcat('answer_postproc.', work_field_name{1,1});
 evalc([str ' = work_field']);
 end
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
 % Обратная нормировка.
 IMP_out(:,1) = IMP(:,1).* std_z(2)' + mean_z(2)';
 IMP_out(:,2) = (IMP(:,2).* std_z(3)' + mean_z(3)').*1e6 ;
 for j=1:(length(IMP(1,:))-2)
 if mod(j,3)==1 %частота
 IMP_out(:,2+j) = IMP(:,2+j).* std_z(2)' + mean_z(2)'; 
 end
 if mod(j,3)==2 %длительность импульса
 IMP_out(:,2+j) = ( IMP(:,2+j).* std_z(3)' + mean_z(3)' ).*1e6; %переводим в мкс
 end
 if mod(j,3)==0 %период
 IMP_out(:,2+j) = ( IMP(:,2+j).* std_z(4)' + mean_z(4)' ).*1e6; %переводим в мкс
 end
 end
 
figure(1)
scatter(IMP_out(:,2)*1e-9,IMP_out(:,3)*1e6,50,'filled')
legend('Принятые импульсы');
grid on
title('Импульсы в пространстве признаков (f,\deltat)')
xlabel('f,ГГц','FontSize',20,...
'FontWeight','bold')
ylabel('\deltat,мкc','FontSize',20,...
'FontWeight','bold')
set(gca, 'FontSize', 20)


% figure(2)
% scatter3(imp(1,:),imp(2,:)*1e-9,imp(3,:)*1e6,50,'filled')
% grid on
% legend('Принятые импульсы');
% title('Импульсы в пространстве признаков (ToA,f,\deltat)')
% ylabel('f,ГГц','FontSize',20,...
% 'FontWeight','bold')
% zlabel('\deltat,мкc','FontSize',20,...
% 'FontWeight','bold')
% xlabel('t,c','FontSize',20,...
% 'FontWeight','bold')
% set(gca, 'FontSize', 20)
% 
% figure(22)
% scatter3(imp(4,:)*1e6,imp(2,:)*1e-9,imp(3,:)*1e6,50,'filled')
% grid on
% legend('Принятые импульсы');
% title('Импульсы в пространстве признаков (ToA,f,\deltat)')
% ylabel('f,ГГц','FontSize',20,...
% 'FontWeight','bold')
% zlabel('\deltat,мкc','FontSize',20,...
% 'FontWeight','bold')
% xlabel('T,мкc','FontSize',20,...
% 'FontWeight','bold')
% set(gca, 'FontSize', 20)
% 
% 
% %построим распределения параметров импульсов
% figure(23)
% subplot(3,1,1)
% histogram(imp(2,:)*1e-9,'Normalization','probability')
% xlabel('f,ГГц','FontSize',20,...
% 'FontWeight','bold')
% subplot(3,1,2)
% histogram(imp(3,:)*1e6,'Normalization','probability')
% xlabel('\deltat,мкc','FontSize',20,...'FontWeight','bold')
% subplot(3,1,3)
% histogram(imp(4,:)*1e6,'Normalization','probability')
% xlabel('T,мкc','FontSize',20,...
% 'FontWeight','bold')
% grid on
 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
 % Вывод результатов в виде наглядных паттернов в кластерах
 total_proc = struct();
 imp=imp';
 for z =1:numel(fieldnames(answer_postproc))
 fieldnames_ans = fieldnames(answer_postproc);
 work_field_name = fieldnames_ans(z);
 index = strsplit(work_field_name{1,1}, 'x');
 index = strsplit(index{1,2}, '_');
 index = str2double(index(1));
 if index == 0 %паттерн из 3 импульсов (8 параметров)
 work_field=eval(strcat('answer_postproc.',work_field_name{1,1}));
 j=1;
 input_data_new=nan(length(work_field),length(IMP_out(1,:))+2);
 for i=1:length(work_field)
 input_data_new(j,1) = work_field(i);
 input_data_new(j,2) = imp(work_field(i),1);
 input_data_new(j,3:end) = IMP_out(work_field(i),:);
 j=j+1;
 end
 evalc(['total_proc.claster', num2str(z) ' = input_data_new']);
 else %больше 3 импульсов
 work_field=eval(strcat('answer_postproc.',work_field_name{1,1}));
 j=1;
 input_data_new=nan(length(work_field) ,length(IMP_out(1,:)) + index*length(imp_norm(1,:))+2);
 for i=1:length(work_field)
 input_data_new(j,1) = work_field(i);
 input_data_new(j,2) = imp(work_field(i),1);
 input_data_new(j,3:length(IMP_out(1,:))+2) = IMP_out(work_field(i),:);
 result = IMP_out((work_field(i)+1:work_field(i)+index),(length(IMP_out(1,:))-length(imp_norm(1,:))+1:length(IMP_out(1,:))));
 result = reshape(result.',1,[]);
 input_data_new(j,length(IMP_out(1,:))+1+1+1:end) = result;
 j=j+1;
 end
 evalc(['total_proc.claster', num2str(z) ' = input_data_new']);
 end
 end
 imp=imp';
%--------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Формирование ответа
total = struct();
answer_postproc = struct();
fieldnames_ans = fieldnames(total_proc);
count = 1;
for v=1:length(fieldnames_ans)
 work_field_name = fieldnames_ans(v);
 work_field = eval(strcat('total_proc.',work_field_name{1,1}));
 work_field = sortrows(work_field, 1);
 
 
 %$crutch$
 if size(work_field, 1) >= 4
 work_field_answer_postproc = work_field(:, 1);
 work_field_total = work_field(:, 2:end);
 evalc(['answer_postproc.claster', num2str(count) ' = work_field_answer_postproc']);
 evalc(['total.claster', num2str(count) ' = work_field_total']);
 count = count + 1;
 end
end
%--------------------------------------------------------------------------
fieldnames_ans = fieldnames(answer_postproc);
count = 1;
for v=1:length(fieldnames_ans)
 work_field_name = fieldnames_ans(v);
 work_field = eval(strcat('answer_postproc.',work_field_name{1,1}));
 work_field_answer_postproc = work_field(:, 1);
 work_field_total = work_field(:, 2:end);
 evalc(['answer_postproc.claster', num2str(count) ' = work_field_answer_postproc']);
 evalc(['total.claster', num2str(count) ' = work_field_total']);
 count = count + 1;
end
%--------------------------------------------------------------------------
% Вывод позиций обработанных кластеров
fprintf(fid,'%s\n','');
% Показываем ответ алгоритма.
total_answer_mas = 0;
fieldnames_ans = fieldnames(answer_postproc);
disp(['Номер реализации - ' num2str(ll)])
disp('Результат с выхода алгоритма(DBSCAN+постобработка): ')
curr_str1 = 'Результат с выхода алгоритма(DBSCAN+постобработка): ';
disp(['Число кластеров - ' num2str(length(fieldnames_ans))])
curr_str2 = ['Число кластеров - ' num2str(length(fieldnames_ans))];
fprintf(fid,'%s\n',curr_str1);
fprintf(fid,'%s\n',curr_str2); 
 for v=1:length(fieldnames(answer_postproc))
 work_field_name = fieldnames_ans(v);
 work_field=eval(strcat('answer_postproc.',work_field_name{1,1}))';
 total_answer_mas(v,1:length(work_field)) = work_field;
 disp(['Позиции ' num2str(v) '-ого кластера - ' num2str(work_field)])
 curr_str = ['Позиции ' num2str(v) '-ого кластера - ' num2str(work_field)];
 fprintf(fid,'%s\n',curr_str);
 end
% Показывает истинный ответ.
true_positions_mas = 0;
fieldnames_ans = fieldnames(true_positions);
disp('Истиное расположение паттернов в кластерах: ')
curr_str1 = 'Истиное расположение паттернов в кластерах: ';
disp(['Число кластеров - ' num2str(length(fieldnames_ans))])
curr_str2 = ['Число кластеров - ' num2str(length(fieldnames_ans))];
fprintf(fid,'%s\n',curr_str1);
fprintf(fid,'%s\n',curr_str2);
 for v=1:length(fieldnames(true_positions))
 work_field_name = fieldnames_ans(v);
 work_field=eval(strcat('true_positions.',work_field_name{1,1}));
 true_positions_mas(v,1:length(work_field)) = work_field;
 disp(['Позиции ' num2str(v) '-ого кластера - ' num2str(work_field') ])
 curr_str = ['Позиции ' num2str(v) '-ого кластера - ' num2str(work_field') ];
 fprintf(fid,'%s\n',curr_str);
 end
fprintf(fid,'%s\n','');
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Блок оценки точности результатов. 
 fprintf("\n");
 fprintf(fid,'%s\n',''); 
 disp(['Реальных кластеров - ' num2str(length(true_positions_mas(:,1))) ' Алгоритм выявил - ' num2str(length(total_answer_mas(:,1)))])
 curr_str = ['Реальных кластеров - ' num2str(length(true_positions_mas(:,1))) ' Алгоритм выявил - ' num2str(length(total_answer_mas(:,1)))];
 fprintf(fid,'%s\n',curr_str);
 stat_matrix(ll,3) = length(total_answer_mas(:,1));
 
 total_answer_mas = sort(reshape(total_answer_mas,1,[]));
 true_positions_mas = sort(reshape(true_positions_mas,1,[]));
 j=1;
 % Блок обработки (исключение нулей).
 for i=1:length(total_answer_mas)
 if total_answer_mas(j) == 0
 total_answer_mas(j)=[];
 if j>length(total_answer_mas)
 break
 end
 else
 j=j+1;
 if j>length(total_answer_mas)
 break
 end
 end
 end
 j=1;
 
 % Блок обработки (исключение нулей).
 for i=1:length(true_positions_mas)
 if true_positions_mas(j) == 0
 true_positions_mas(j)=[];
 if j>length(true_positions_mas)
 break
 end
 else
 j=j+1;
 if j>length(true_positions_mas)
 break
 end
 end
 end
 
 C1 = intersect(total_answer_mas,true_positions_mas);
 if length(total_answer_mas)>length(true_positions_mas)
 C2 = setdiff(total_answer_mas,true_positions_mas);
 else
 C2 = setdiff(true_positions_mas,total_answer_mas);
 end 
 
 FAl = 100*length(C2)/(N-max(vector_N)+1-sum(vector_k)); %False Alarm (не совсем правильно)
 accuracy = length(C1)/length(true_positions_mas)*100;
 disp(['Точность распознавания паттернов алгоритмом на данной реализации составила -' num2str(accuracy) '%'])
 disp(['Процент определенных ложных паттернов - ' num2str(FAl) '%'])
 curr_str1 = ['Точность распознавания паттернов алгоритмом на данной реализации составила - ' num2str(accuracy) '%'];
 curr_str2 = ['Процент определенных ложных паттернов - ' num2str(FAl) '%'];
 fprintf(fid,'%s\n',curr_str1);
 fprintf(fid,'%s\n',curr_str2);
 
 stat_matrix(ll,1) = accuracy;
 stat_matrix(ll,2) = FAl;
 fprintf("\n");
 fprintf(fid,'%s\n','');
 fprintf(fid,'%s\n','');
 fprintf(fid,'%s\n','');
 
end 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Общая статистика по всем реализациям.
Mean_accuracy = sum(stat_matrix(:,1))/rr;
Mean_FAl = sum(stat_matrix(:,2))/rr;
Mean_clusters = sum(stat_matrix(:,3))/rr;
fprintf("\n");
fprintf(fid,'%s\n','');
disp(['Результат кластеризации на ' num2str(rr) ' случаях, средняя точность определения паттернов составила - ' num2str(Mean_accuracy) '%']) 
disp(['Средний процент определенных ложных паттернов составил - ' num2str(Mean_FAl) 
'%']) 
disp(['Среднее число выявляемых кластеров составило - ' num2str(Mean_clusters)])
curr_str1 = ['Результат кластеризации на ' num2str(rr) ' случаях, средняя точность определения паттернов составила - ' num2str(Mean_accuracy) '%'];
curr_str2 = ['Средний процент определенных ложных паттернов составил - ' 
num2str(Mean_FAl) '%'];
curr_str3 = ['Среднее число выявляемых кластеров составило - ' num2str(Mean_clusters)];
fprintf(fid,'%s\n',curr_str1);
fprintf(fid,'%s\n',curr_str2);
fprintf(fid,'%s\n',curr_str3);
%--------------------------------------------------------------------------
% Закрытие файлов.
fclose(fid); 
fclose all;
% Конец замеров тестов
toc
%--------------------------------------------------------------------------
% Функция формирования имплуьсов 
function [ imp ] = make_impulse( freq, dur, T_min, prev_T )
 
 delta = T_min * randi(10);
 ToA = prev_T + delta; % Формирование времени прихода импульса.
 
 freq1 = freq(randi(length(freq))); % Формирование несущей импульса.
 f = freq1 + randi([0 1])*(-1)^randi([0 1]) * 0.001 * randi(6) * freq1;
 
 dur1 = dur(randi(length(dur))); % Формирование длительности импульса.
 d = dur1 + 0*randi([0 1])*(-1)^randi([0 1]) * 0.1 * randi(6) * dur1;
 imp = [ToA; f; d; delta];
end