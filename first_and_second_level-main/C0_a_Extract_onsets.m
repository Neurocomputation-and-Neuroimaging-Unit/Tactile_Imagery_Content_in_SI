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
log_dir='C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\OriginalLogs';

condnames = {'StimPress','StimFlutt','StimVibro','ImagPress','ImagFlutt','ImagVibro','Null_1','Null_2','Att_1','Att_2','Motion'};  % condition names

for s = sbjs %1:length(SJs)
    cd(log_dir)
    
    outputdir=[src_dir filesep SJs{s} filesep 'logs'];
    if ~exist(outputdir, 'dir')
        mkdir(outputdir)
    end
    
    runs=dir(['FMRI_log_' num2str(s) '_*.mat']);
    
    onsets={};
    if ~isempty(runs)
        for r=1:length(runs)
            load(runs(r).name)
            temprun = ['run0' num2str(r)];
            
             %% AUSANHMSWEISE (für worked/failed)
%             % cond 1: stimulation pressure (2 Hz)
%             onsets{r,1} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==1) & (log_ExPra19.Design(4,:)==2))/1000;
%             % cond 2: stimulation flutter (30 Hz)
%             onsets{r,2} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==1) & (log_ExPra19.Design(4,:)==3))/1000;
%             % cond 3: stimulation vibration (150 Hz)
%             onsets{r,3} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==1) & (log_ExPra19.Design(4,:)==1))/1000;
%             
%             % cond 4: imagination press successfull (answer = 2)
%             onsets{r,4} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==2) & (log_ExPra19.Design(4,:)==2) & (log_ExPra19.responses(1,:)==2))/1000;
%             % cond 5: imagination flutter successfull (answer = 2)
%             onsets{r,5} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==2) & (log_ExPra19.Design(4,:)==3) & (log_ExPra19.responses(1,:)==2))/1000;
%             % cond 6: imagination vibro successfull (answer = 2)
%             onsets{r,6} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==2) & (log_ExPra19.Design(4,:)==1) & (log_ExPra19.responses(1,:)==2))/1000;
%             
%             % cond 7: null
%             onsets{r,7} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==0) & (log_ExPra19.Design(4,:)==0))/1000;
%             
%             % cond 8: attention
%             onsets{r,8} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==3) & (log_ExPra19.Design(4,:)>0))/1000;
%             
%             % cond 9: motion regressor
%             onsets{r,9} = (log_ExPra19.Design(1,(log_ExPra19.Design(3,:)>0) & (log_ExPra19.responses(1,:)>0))/1000)+4;
%             %onset der Reaktion per Knopfdruck: beginn trial + prep-time(1s) +
%             %stim+time(3s) + answer-time(2s) -> zu ENDE der answer-time???
%             
%             tempOnsData = { condnames{1}, num2str(onsets{r,1});
%                             condnames{2}, num2str(onsets{r,2});
%                             condnames{3}, num2str(onsets{r,3});
%                             condnames{4}, num2str(onsets{r,4});
%                             condnames{5}, num2str(onsets{r,5});
%                             condnames{6}, num2str(onsets{r,6});
%                             condnames{7}, num2str(onsets{r,7});
%                             condnames{8}, num2str(onsets{r,8});
%                             condnames{9}, num2str(onsets{r,9})};

            %% Gründe: Null trials random in 2 Gruppen, als bessere Baseline für Conjunction (??)     
            
            % cond 1: stimulation pressure (2 Hz)
            onsets{r,1} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==1) & (log_ExPra19.Design(4,:)==2))/1000;
            % cond 2: stimulation flutter (30 Hz)
            onsets{r,2} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==1) & (log_ExPra19.Design(4,:)==3))/1000;
            % cond 3: stimulation vibration (150 Hz)
            onsets{r,3} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==1) & (log_ExPra19.Design(4,:)==1))/1000;
            
            % cond 4: imagination pressure (2 Hz)
            onsets{r,4} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==2) & (log_ExPra19.Design(4,:)==2))/1000;
            % cond 5: imagination flutter (30 Hz)
            onsets{r,5} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==2) & (log_ExPra19.Design(4,:)==3))/1000;
            % cond 6: imagination vibration (150 Hz)
            onsets{r,6} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==2) & (log_ExPra19.Design(4,:)==1))/1000;
            
            % cond 7: null
            rando = randperm(6);
            nulls = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==0) & (log_ExPra19.Design(4,:)==0))/1000;
            onsets{r,7} = sort(nulls(rando([1 2 3])));
            onsets{r,8} = sort(nulls(rando([4 5 6])));
            
            % cond 8: attention
%             onsets{r,9} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==3) & (log_ExPra19.Design(4,:)>0))/1000;
            rando2 = randperm(6);
            atts = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==3) & (log_ExPra19.Design(4,:)>0))/1000;
            onsets{r,9} = sort(atts(rando2([1 2 3])));
            onsets{r,10} = sort(atts(rando2([4 5 6])));
            
            % cond 9: motion regressor
            onsets{r,11} = (log_ExPra19.Design(1,(log_ExPra19.Design(3,:)>0) & (log_ExPra19.responses(1,:)>0))/1000)+4;
            %onset der Reaktion per Knopfdruck: beginn trial + prep-time(1s) +
            %stim+time(3s) + answer-time(2s) -> zu ENDE der answer-time???
            
            tempOnsData = { condnames{1}, num2str(onsets{r,1});
                            condnames{2}, num2str(onsets{r,2});
                            condnames{3}, num2str(onsets{r,3});
                            condnames{4}, num2str(onsets{r,4});
                            condnames{5}, num2str(onsets{r,5});
                            condnames{6}, num2str(onsets{r,6});
                            condnames{7}, num2str(onsets{r,7});
                            condnames{8}, num2str(onsets{r,8});
                            condnames{9}, num2str(onsets{r,9});
                            condnames{10}, num2str(onsets{r,10});
                            condnames{11}, num2str(onsets{r,11})};
                        
             %% NORMALERWEISE      

%             % cond 1: stimulation pressure (2 Hz)
%             onsets{r,1} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==1) & (log_ExPra19.Design(4,:)==2))/1000;
%             % cond 2: stimulation flutter (30 Hz)
%             onsets{r,2} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==1) & (log_ExPra19.Design(4,:)==3))/1000;
%             % cond 3: stimulation vibration (150 Hz)
%             onsets{r,3} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==1) & (log_ExPra19.Design(4,:)==1))/1000;
%             
%             % cond 4: imagination pressure (2 Hz)
%             onsets{r,4} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==2) & (log_ExPra19.Design(4,:)==2))/1000;
%             % cond 5: imagination flutter (30 Hz)
%             onsets{r,5} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==2) & (log_ExPra19.Design(4,:)==3))/1000;
%             % cond 6: imagination vibration (150 Hz)
%             onsets{r,6} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==2) & (log_ExPra19.Design(4,:)==1))/1000;
%             
%             % cond 7: null
%             onsets{r,7} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==0) & (log_ExPra19.Design(4,:)==0))/1000;
%             
%             % cond 8: attention
%             onsets{r,8} = log_ExPra19.Design(1,(log_ExPra19.Design(3,:)==3) & (log_ExPra19.Design(4,:)>0))/1000;
%             
%             % cond 9: motion regressor
%             onsets{r,9} = (log_ExPra19.Design(1,(log_ExPra19.Design(3,:)>0) & (log_ExPra19.responses(1,:)>0))/1000)+4;
%             %onset der Reaktion per Knopfdruck: beginn trial + prep-time(1s) +
%             %stim+time(3s) + answer-time(2s) -> zu ENDE der answer-time???
%             
%             tempOnsData = { condnames{1}, num2str(onsets{r,1});
%                             condnames{2}, num2str(onsets{r,2});
%                             condnames{3}, num2str(onsets{r,3});
%                             condnames{4}, num2str(onsets{r,4});
%                             condnames{5}, num2str(onsets{r,5});
%                             condnames{6}, num2str(onsets{r,6});
%                             condnames{7}, num2str(onsets{r,7});
%                             condnames{8}, num2str(onsets{r,8});
%                             condnames{9}, num2str(onsets{r,9})};
            
            tempxlsname = [outputdir filesep 'Onsets4_' SJs{s} '_' temprun '.xls'];
%             xlswrite(tempxlsname, tempOnsData);
% 
%             delete(tempxlsname);
            writecell(tempOnsData, tempxlsname);
            
            clear temp*
            
        end 
        cd(outputdir)
        eval(['save all_onsets4_' SJs{s} '.mat onsets condnames'])
        clear onsets runs outputdir
    end %end if
end
