
function [average] = averageResults(src_dir, sj, folder_names, outputDir)

cd(src_dir)
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

sjDir = [src_dir '\' SJs{sj}];

relevant_folders = folder_names;

input = zeros(8,1);

for f = 1:size(relevant_folders, 2)
    
    this_input = load([sjDir cell2mat(relevant_folders(:,f)) '\res_accuracy_minus_chance.mat']);
    input = [input + this_input(1,1).results.accuracy_minus_chance.output];
    
end

average = input/size(relevant_folders,2);

writematrix(average, [cell2mat(outputDir) '.txt']);