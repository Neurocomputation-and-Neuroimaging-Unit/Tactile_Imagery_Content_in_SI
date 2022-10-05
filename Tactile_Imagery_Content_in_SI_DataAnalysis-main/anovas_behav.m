%% ANOVAs
%1st: eine Variable (Press Flutt Vibro der Imag), 3 Stufen -> Performance?
    % load data (see "plot_responses")
    % take imag-part -> anova???
%2nd: eine Variable, 3 Stufen, Messwiederholung (6 Runs) -> Performance?

%% load

clear all

%data source directory
src_dir      = 'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Data';

%subject identifiers
cd(src_dir)
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

sbjs = [1 2 4 5 6 7 8 10 12 13 14 16 17 18 19 20 21 22 23 24 26 27 28 29 30 31 32]; 

% directory with original logfiles
src_dir='C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Data';
all_responses = [];
all_data = [];
all_data3 = [];

for s = sbjs
    
    sbj_data2 = []; % 3 conds, 3 values
    sbj_data3 = [];
        
    log_dir = strcat(src_dir, filesep, SJs{s}, filesep, "logs");
    cd(log_dir);
    
    for r = [1 2 3 4 5 6]   
        
        name_run_file = strcat('responses_', num2str(SJs{s}),"_run0", num2str(r), '.txt');
        
        if exist(name_run_file) == 2
            
%             display(name_run_file);
            
            responses = load(strcat('responses_', num2str(SJs{s}),"_run0", num2str(r), '.txt'));
            
            % slice idx [4 5 6]
            imag_resp = responses(:, [4 5 6]); % press flutt vibro
            count_resp = [];
            edges = [-0.5 0.5 1.5 2.5];
            
            for i = [1 2 3]
                count_resp = [count_resp; histcounts(imag_resp(:,i), edges)];
            end
            
        else
            
            display('###################################################################')
            display(['################## ' SJs{s} ' run: ' num2str(r) ', ERROR ###################'])
            display('###################################################################')
    
        end
        
        sbj_data2 = [sbj_data2; count_resp(:,3)'];
        sbj_data3 = [sbj_data3 count_resp(:,3)];
        
    end
    
    all_data = [all_data; sbj_data2];
    all_data3 = [all_data3; sbj_data3];

end

group = {'Press' 'Flutt' 'Vibro'};
[anv other] = anova1(all_data, group)

freq = {'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'; 'Press'; 'Flutt'; 'Vibro'};
t = table(freq, all_data3(:, 1), ...
                all_data3(:, 2), ...
                all_data3(:, 3), ...
                all_data3(:, 4), ...
                all_data3(:, 5), ...
                all_data3(:, 6), ...
    'VariableNames', {'freq', 'run1', 'run2', 'run3', 'run4', 'run5', 'run6'});
time = table([1 2 3 4 5 6]', 'VariableNames', {'runs'});
rm = fitrm(t, 'run1-run6~freq', 'WithinDesign', time);
ranova(rm)
close all

