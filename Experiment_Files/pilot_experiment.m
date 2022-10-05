%% ========================================================================

%  ACUPUNCTURE PILOT EXPERIMENT

% =========================================================================
% Prepare experiment
% =========================================================================
% Basics

% cd /Users/Ringelblume/Desktop/Dropbox/IMACU/Pilotexperiment;
sca;
clc;
close all;
clear all;
clearvars;

% Default settings for setting up Psychtoolbox
% PsychDefaultSetup(2);

% ========================================================================
% Adjustables

fmri = 1;
stimul = 1;
par = input('Participant number: ')
startrun = input('Start with run: ')

balancecue = mod(par, 6); % for balancing cue symbol
balancecolour = mod(par, 2); % for balancing cue colour
numruns = 7; % number of blocks (including practice block)
pract = 48; % number of practice trials for trial indices
stimdur = 3;
predur = 1;
ansdur = 2;
anskeys = {'1!', '2@'};
ansbalance = {'ja               nein', 'nein               ja'};
rl = [0 1];
flutter_f = 30;
vibro_f = 150;
base_f = 2;
isi = [2 3 4 5 6]; % inter-stimulus interval in seconds
winheight = 600; % window height
winwidth = 800; % window width

prj_dir = pwd;

% =========================================================================
% Configure quaerosys

if stimul
    amp = 4095;
    fade_out = 400;
    res = 1; %0.5;  % resoultion of signal in ms (QueroSys Stimulator has a maximal resolution of 0.5ms)
    t = [res:res:stimdur*1000]/1000;

%     ParAddr=(hex2dec('378')); %Parallel Port Adress
%     lptwrite(ParAddr,0);  %what'that?

%     loadlibrary('C:\Jan\TacMoDo\Quaerobox\QuaeroSys\stimlib0.dll', 'C:\Jan\TacMoDo\Quaerobox\QuaeroSys\stimlibrel.h')
    loadlibrary('stimlib0.dll', 'stimlibrel.h')
    calllib('stimlib0', 'setProperty', 'local_buffer_size', 1000000)
    calllib('stimlib0', 'initStimulator', '98J9007703J8RJHA9J27J0B053KBKR'); % Licence oktober 16
    
    % Creating basic pulse and frequency modulations (at half frequency
    % to allow for abs() later on)
    base_pulse = sin(pi*base_f*t);
    vibro_modul = (amp/2).*sin(pi*vibro_f*t);
    flutter_modul = (amp/2).*sin(pi*flutter_f*t);
    base_modul = (amp/2);

    % Create a fade out of the signal
%     damping_win = ones(1,length(t));
%     damping_win(length(damping_win)-fade_out:length(damping_win)) = 1:-1/fade_out:0;
%     damping_win(1:fade_out+1) = 0:1/fade_out:1;

    vibro_stim = abs(vibro_modul.*base_pulse).*2;
    
    flutter_stim = abs(flutter_modul.*base_pulse).*2;
    
    base_stim = abs(base_modul.*base_pulse).*2;

%     base_stim = [0:amp/500:amp repmat(amp, 1, 500) ...
%         amp:-1*amp/500:0 0:amp/500:amp repmat(amp, 1, 500) amp:-1*amp/500:0];
%     base_stim_att = base_stim/2 + amp_correction/2;
end

if balancecue == 1
    vibro3cue = 'star.png';
    vibro1cue = 'circle.png';
    vibro2cue = 'pentagon.png';
elseif balancecue == 2
    vibro3cue = 'circle.png';
    vibro1cue = 'star.png';
    vibro2cue = 'pentagon.png';
elseif balancecue == 3
    vibro3cue = 'star.png';
    vibro1cue = 'pentagon.png';
    vibro2cue = 'circle.png';
elseif balancecue == 4
    vibro3cue = 'pentagon.png';4
    vibro1cue = 'star.png';
    vibro2cue = 'circle.png';
elseif balancecue == 5
    vibro3cue = 'pentagon.png';
    vibro1cue = 'circle.png';
    vibro2cue = 'star.png';
elseif balancecue == 0
    vibro3cue = 'circle.png';
    vibro1cue = 'pentagon.png';
    vibro2cue = 'star.png';
end

if balancecolour == 1
    stimcolour = 'blue';
    imagcolour = 'yellow';
elseif balancecolour == 0
    stimcolour = 'yellow';
    imagcolour = 'blue';
end

% =========================================================================

% Read in stimulus from csv file

trialpool = readtable('trialpool.csv', 'Delimiter', ';', 'ReadVariableNames', true);
blocklen = height(trialpool);


% Set text strings for greeting, instructions etc.
greet = ['Herzlich Willkommen zu unserem Experiment! \n\n' ...
    'Um fortzufahren, drücken Sie bitte \n\n' ...
    'eine beliebige Taste.'];
instr1 = ['In den folgenden ' num2str(numruns) ' Blöcken \n\n' ...
    'wird Ihr Zeigefinger abwechselnd auf verschiedene \n\n' ...
    'Arten stimuliert, beziehungsweise Sie werden gebeten, \n\n' ...
    'sich die Stimulation intensiv vorzustellen. \n\n' ...
    'Mehr dazu auf der nächsten Seite. \n\n' ...
    'Weiter mit beliebiger Taste.'];
instr2 = ['Vor jedem Durchgang sehen Sie ein kleines Symbol, \n\n' ...
    'das Ihnen anzeigt, welche Stimulation als nächstes folgt \n\n' ...
    'bzw. welche Stimulation Sie sich vorstellen sollen. \n\n' ...
    'Sie haben ' num2str(predur) ' Sekunde um sich vorzubereiten, \n\n' ...
    'dann wird das Symbol größer und \n\n' ... 
    'die Stimulation/Imagination beginnt. \n\n' ...
    'Mehr zu den Symbolen auf der nächsten Seite. \n\n' ...
    'Weiter mit beliebiger Taste.'];
instr3 = ['Vor einem Durchgang mit STIMULATION \n\n' ...
    'erscheint eins der folgenden Symbole: \n\n' ...
    'Flutterpuls:                        \n\n' ...
    'Vibrationspuls:                  \n\n' ...
    'Druckpuls:                         \n\n' ...
    'Weiter mit beliebiger Taste.'];
instr4 = ['Vor einem Durchgang mit IMAGINATION \n\n' ...
    'erscheint eins der folgenden Symbole: \n\n' ...
    'Flutterpuls:                        \n\n' ...
    'Vibrationspuls:                  \n\n' ...
    'Druckpuls:                         \n\n' ...
    'Weiter mit beliebiger Taste.'];
instr5 = ['In der Stimulationsbedingung erfolgt nicht bei jedem \n\n' ...
    'Durchgang eine Stimulation. Nach jedem Durchgang \n\n' ...
    'werden Sie gebeten, mittels der Knöpfe \n\n' ...
    'zu entscheiden, ob Sie etwas gespürt haben oder nicht. \n\n' ...
    'Die Antwortoptionen erscheinen auf dem Bildschirm, \n\n' ...
    'bitte drücken Sie entsprechend rechts oder links. \n\n' ...
    'Starte mit beliebiger Taste!'];
stim_quest = ['Gab es eine Stimulation?'];
imag_quest = ['Ist Ihnen die Imagination gelungen?'];
paus = ['Pause\n\n'...
    'Um fortzufahren, bitte eine beliebige Taste drücken.'];
bye = ['Vielen Dank für Ihre Teilnahme!'];

% Import cue images and make textures with transparent backgrounds
[cross, ~, alpha] = imread('Cues/cross.png');
cross(:, :, 4) = alpha;
[imagvibro3symb, ~, alpha] = imread(['Cues/' imagcolour '_' vibro3cue]);
imagvibro3symb(:, :, 4) = alpha;
[imagvibro1symb, ~, alpha] = imread(['Cues/' imagcolour '_' vibro1cue]);
imagvibro1symb(:, :, 4) = alpha;
[imagvibro2symb, ~, alpha] = imread(['Cues/' imagcolour '_' vibro2cue]);
imagvibro2symb(:, :, 4) = alpha;
[stimvibro3symb, ~, alpha] = imread(['Cues/' stimcolour '_' vibro3cue]);
stimvibro3symb(:, :, 4) = alpha;
[stimvibro1symb, ~, alpha] = imread(['Cues/' stimcolour '_' vibro1cue]);
stimvibro1symb(:, :, 4) = alpha;
[stimvibro2symb, ~, alpha] = imread(['Cues/' stimcolour '_' vibro2cue]);
stimvibro2symb(:, :, 4) = alpha;


% Configure parallel port for scanner trigger

if fmri
    disp('Running experiment with scanner triggers')
    daqregister('C:\Programme\MATLAB\R2009a\toolbox\daq\daq\private\mwparallel.dll')
    reg_timing = [];
    po = daqhwinfo; 
    po.InstalledAdaptors
    dio = digitalio('parallel', 'LPT1');
    in_line = addline(dio, 11, 'In');
    dat = 0;
else disp('Running experiment without scanner')
end

%% START COGENT (just for creating logs because they do it this way in Berlin)

% START COGENT
% start_cogent;
addpath(strcat(prj_dir, filesep, 'logs'));

%% ========================================================================
% Start experiment
% =========================================================================

% Start screen

[win, rect] = Screen('OpenWindow', 2, [150 150 150], [-1*winwidth 0 0 winheight]);
%[win, rect] = Screen('OpenWindow', 0, [150 150 150], [0 0 winwidth winheight]);
precuesize = [0 0 winwidth/28 winwidth/28];
[preCue, xOffsetsigS, yOffsetsigS] = CenterRect(precuesize, rect);
maincuesize = [0 0 winwidth/16 winwidth/16];
[mainCue, xOffsetsigS, yOffsetsigS] = CenterRect(maincuesize, rect);
instcuepos1 = [winwidth*15/20 winheight*9/20 winwidth*15/20+winwidth/28 winheight*9/20+winwidth/28];
instcuepos2 = [winwidth*15/20 winheight*10.5/20 winwidth*15/20+winwidth/28 winheight*10.5/20+winwidth/28];
instcuepos3 = [winwidth*15/20 winheight*12/20 winwidth*15/20+winwidth/28 winheight*12/20+winwidth/28];

% Select specific text font, style and size and prepare cue textures  

Screen('TextFont', win, 'Arial');
Screen('TextSize', win, 20);
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

imagvibro3tex = Screen('MakeTexture', win, imagvibro3symb);
imagvibro1tex = Screen('MakeTexture', win, imagvibro1symb);
imagvibro2tex = Screen('MakeTexture', win, imagvibro2symb);
stimvibro3tex = Screen('MakeTexture', win, stimvibro3symb);
stimvibro1tex = Screen('MakeTexture', win, stimvibro1symb);
stimvibro2tex = Screen('MakeTexture', win, stimvibro2symb);
crosstex      = Screen('MakeTexture', win, cross);

HideCursor([0]);
KbQueueCreate(-1);
% KbQueueFlush;
% KbQueueStart;

% Show greeting and instruction text

DrawFormattedText(win, greet, 'center', 'center', [0 0 0]);
Screen(win,'flip');
KbQueueWait(-1);

DrawFormattedText(win, '', 'center', 'center', [0 0 0]);
Screen(win,'flip');
WaitSecs(0.5);

DrawFormattedText(win, instr1, 'center', 'center', [0 0 0]);
Screen(win,'flip');
KbQueueWait(-1);

DrawFormattedText(win, instr2, 'center', 'center', [0 0 0]);
Screen(win,'flip');
KbQueueWait(-1);

Screen('DrawTexture', win, stimvibro3tex, [], instcuepos1);
Screen('DrawTexture', win, stimvibro1tex, [], instcuepos2);
Screen('DrawTexture', win, stimvibro2tex, [], instcuepos3);
DrawFormattedText(win, instr3, 'center', 'center', [0 0 0]);
Screen(win,'flip');
KbQueueWait(-1);

Screen('DrawTexture', win, imagvibro3tex, [], instcuepos1);
Screen('DrawTexture', win, imagvibro1tex, [], instcuepos2);
Screen('DrawTexture', win, imagvibro2tex, [], instcuepos3);
DrawFormattedText(win, instr4, 'center', 'center', [0 0 0]);
Screen(win,'flip');
KbQueueWait(-1);

DrawFormattedText(win, instr5, 'center', 'center', [0 0 0]);
Screen(win,'flip');
KbQueueWait(-1);
DrawFormattedText(win, '', 'center', 'center', [0 0 0]);
Screen(win,'flip');

% =========================================================================
% Loop over blocks plus practice block

KbQueueCreate(-1);

for m = startrun : numruns
    
    %============================ Initiate log files ======================
    % To avoid overwriting of Logfiles the start time is included in the Log-file-Names
    uhr = clock;    % clock returns a six element date vector containing the current time
                    % and date in decimal form: [year month day hour
                    % minute
                    % seconds]
    zeit  = ['_t_' num2str(uhr(4)) '_' num2str(uhr(5))];
    if fmri == 1 && m >= 2
        log_path =       [prj_dir filesep 'logs' filesep 'FMRI_log_'  num2str(par)  '_'   num2str(m)  zeit  '.mat'];
        config_log(strcat(prj_dir, filesep, 'logs', filesep, 'FMRI_log_', num2str(par), '_',  num2str(m), zeit, '.txt'));
    else
        log_path =       [prj_dir filesep 'logs' filesep 'TRAIN_log_'  num2str(par)  '_'   num2str(m)  zeit  '.mat'];
        config_log(strcat(prj_dir, filesep, 'logs', filesep, 'TRAIN_log_', num2str(par), '_',  num2str(m), zeit, '.txt'));
    end
    
    start_cogent;
    
    % COGENT-LOGFILE
    logstring(strcat('Subject: ', num2str(par) , '  Run Number: ', num2str(m)));
    logstring(strcat('Start-time:  ', num2str(uhr(4)), ':', num2str(uhr(5)), 'Uhr  -' , ' Date: ', date ));
    
    if m == 1
        sequence = trialpool;
        sequence.Properties.VariableNames = trialpool.Properties.VariableNames;
    else
        inds = randperm(length([1:height(trialpool)]));
        sequence = trialpool(inds,:);
        sequence.Properties.VariableNames = trialpool.Properties.VariableNames;
    end
    
    rl2 = Shuffle(repmat(rl, 1, ceil(height(sequence)/length(rl))));
    rl2 = rl2(1:blocklen);
    
    log_ExPra19.subject_name    = {''};
    log_ExPra19.subject_no      = num2str(par);
    log_ExPra19.symbols.vibro       = vibro1cue;
    log_ExPra19.symbols.pressure    = vibro2cue;
    log_ExPra19.symbols.flutter     = vibro3cue;
    log_ExPra19.colours.stim        = stimcolour;
    log_ExPra19.colours.imag        = imagcolour;
    if fmri == 1 && m >= 2
        log_ExPra19.fmri            = fmri;
    else
        log_ExPra19.fmri            = 0;
    end
    log_ExPra19.run             = m;
    log_ExPra19.date            = date;
    log_ExPra19.time            = clock;
    log_ExPra19.stimdur         = stimdur;
    log_ExPra19.Design          = zeros(5, blocklen);
    isitemp                     = Shuffle(repmat(isi, 1, ceil(blocklen/length(isi))));
    isitemp                     = isitemp(1:blocklen);  
    log_ExPra19.Design(5,:)     = isitemp;
    isi2                        = isitemp + (predur + stimdur + ansdur);
    log_ExPra19.Design(1,:)     = (cumsum(isi2)-(predur + stimdur + ansdur))*1000;
    log_ExPra19.Design(2,:)     = log_ExPra19.Design(1,:)/2000;

    log_ExPra19.responses       = zeros(3, blocklen);
    log_ExPra19.responses(3,:)  = rl2;
    save(log_path, 'log_ExPra19');
    
    % EXPLANATION OF THE DESIGN MATRIX
    % 1: Trial onset (beginning of proper trial) according to plan
    % 2: Volume at pre-trial onset
    % 3: Stimulus condition (1="stim", 2="imag", 3="att")
    % 4: Stimulus modality (1="vibration", 2="pressure", 3="flutter")
    % 5: ISI before this trial
    % 6: Trial onset (beginning of pre-trial) according to plan
    
    % responses
    % 1: answer given (2 = participant answered yes, 1 = participant answered no,
    %       0 = participant did not press either button)
    % 2: correctness (for stimulation & detection conditions, 2 = correct, 
    %  1 = incorrect, 0 = no answer or imagination trial)
    % 3: right-left answer coding per trial (0 = yes on left button, 1 = yes on right button)
    
    if fmri && m >= 2 % Wait for scanner trigger
        disp('Waiting for Trigger...')
        while 1
            dat = getvalue(dio);
            if dat == 1;
                logstring('FIRST TRIGGER from the MRI scanner');
                disp('FIRST TRIGGER from the MRI scanner');
                t_first = GetSecs; % time counter
            break;
            end
        end
    else
        t_first = GetSecs; % time counter
        logstring('START TIME (In fMRI: FIRST TRIGGER from the MRI scanner)');
        disp('START TIME (In fMRI: FIRST TRIGGER from the MRI scanner)');
    end
    

    % =========================================================================
    % Loop over stimuli within current block

    for i = 1 : blocklen
        % Show fixation cross
        Screen('DrawTexture', win, crosstex, [], mainCue);
        Screen('Flip', win);
        logstring(['Beginning of ISI']);
        disp(['Beginning of ISI']);

        actiontime = log_ExPra19.Design(1,i)/1000;
        logstring(['fMRI-ActionTime' num2str(i) '  ' num2str(actiontime)]);
        logstring(['Onset in volume: ' num2str(log_ExPra19.Design(2,i))]);
        
        KbQueueFlush;
        KbQueueStart;

        % Show stimulus
        if strcmp(sequence.Condition(i), 'stim')
            log_ExPra19.Design(3,i) = 1;
            if strcmp(sequence.Stimulation(i), 'vibro1')
                log_ExPra19.Design(4,i) = 1;
                if stimul    
                    calllib('stimlib0', 'setPinBlock', 0, 5, 1, 1, 1, 1, 1, 1, 1, 1);
                    calllib('stimlib0', 'setPinBlock', 1, 5, 1, 1, 1, 1, 1, 1, 1, 1);
    %                 calllib('stimlib0', 'waitForTrigger', 16, 0);
                    for s=1:length(vibro_stim)
                        calllib('stimlib0', 'setDAC', 1, round(vibro_stim(s)));
                        calllib('stimlib0', 'wait', 5, res*2);
                    end
                end
                symbol = [stimcolour '_' vibro1cue];
                logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
                while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
                end
                disp(['Pre-cue onset']);
                Screen('DrawTexture', win, stimvibro1tex, [], preCue);
                Screen('Flip', win);
                logstring(['Pre-cue: ' symbol]);
                tic
                if stimul
                    WaitSecs(predur-0.51); % the "startStimulation" command takes 510ms to execute in this condition
                    calllib('stimlib0', 'startStimulation');
                    timing_control(i)=toc;
                    calllib('stimlib0', 'stopStimulation');
                    Screen('DrawTexture', win, stimvibro1tex, [], mainCue);
                    Screen('Flip', win);
                    logstring(['Cue: ' symbol]);
                    WaitSecs(stimdur);
                    calllib('stimlib0', 'setDAC', 1, 0);
                    calllib('stimlib0', 'wait', 5, 1);
                    calllib('stimlib0', 'startStimulation');
                    calllib('stimlib0', 'stopStimulation');
                else
                    WaitSecs(predur);
                    Screen('DrawTexture', win, stimvibro1tex, [], mainCue);
                    Screen('Flip', win);
                    logstring(['Cue: ' symbol]);
                    WaitSecs(stimdur);
                end
            elseif strcmp(sequence.Stimulation(i), 'vibro2')
                log_ExPra19.Design(4,i) = 2;
                if stimul    
                    calllib('stimlib0', 'setPinBlock', 0, 5, 1, 1, 1, 1, 1, 1, 1, 1);
                    calllib('stimlib0', 'setPinBlock', 1, 5, 1, 1, 1, 1, 1, 1, 1, 1);
    %                 calllib('stimlib0', 'waitForTrigger', 16, 0);
                    for s=1:length(base_stim)
                        calllib('stimlib0', 'setDAC', 1, round(base_stim(s)));
                        calllib('stimlib0', 'wait', 5, res*2);
                    end
                end
                symbol = [stimcolour '_' vibro2cue];
                logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
                while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
                end
                disp(['Pre-cue onset']);
                Screen('DrawTexture', win, stimvibro2tex, [], preCue);
                Screen('Flip', win);
                logstring(['Pre-cue: ' symbol]);
                tic
                if stimul
                    WaitSecs(predur-0.51); % the "startStimulation" command takes 420ms to execute in this condition
                    calllib('stimlib0', 'startStimulation');
                    timing_control(i)=toc;
                    calllib('stimlib0', 'stopStimulation');
                    Screen('DrawTexture', win, stimvibro2tex, [], mainCue);
                    Screen('Flip', win);
                    logstring(['Cue: ' symbol]);
                    WaitSecs(stimdur);
                    calllib('stimlib0', 'setDAC', 1, 0);
                    calllib('stimlib0', 'wait', 5, 1);
                    calllib('stimlib0', 'startStimulation');
                    calllib('stimlib0', 'stopStimulation');
                else
                    WaitSecs(predur);
                    Screen('DrawTexture', win, stimvibro2tex, [], mainCue);
                    Screen('Flip', win);
                    logstring(['Cue: ' symbol]);
                    WaitSecs(stimdur);
                end
            elseif strcmp(sequence.Stimulation(i), 'vibro3')
                log_ExPra19.Design(4,i) = 3;
                if stimul    
                    calllib('stimlib0', 'setPinBlock', 0, 5, 1, 1, 1, 1, 1, 1, 1, 1);
                    calllib('stimlib0', 'setPinBlock', 1, 5, 1, 1, 1, 1, 1, 1, 1, 1);
    %                 calllib('stimlib0', 'waitForTrigger', 16, 0);
                    for s=1:length(flutter_stim)
                        calllib('stimlib0', 'setDAC', 1, round(flutter_stim(s)));
                        calllib('stimlib0', 'wait', 5, res*2);
                    end
                end
                symbol = [stimcolour '_' vibro3cue];
                logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
                while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
                end
                disp(['Pre-cue onset']);
                Screen('DrawTexture', win, stimvibro3tex, [], preCue);
                Screen('Flip', win);
                logstring(['Pre-cue: ' symbol]);
                tic
                if stimul
                    WaitSecs(predur-0.51); % the "startStimulation" command takes 510ms to execute in this condition
                    calllib('stimlib0', 'startStimulation');
                    timing_control(i)=toc;
                    calllib('stimlib0', 'stopStimulation');
                    Screen('DrawTexture', win, stimvibro3tex, [], mainCue);
                    Screen('Flip', win);
                    logstring(['Cue: ' symbol]);
                    WaitSecs(stimdur);
                    calllib('stimlib0', 'setDAC', 1, 0);
                    calllib('stimlib0', 'wait', 5, 1);
                    calllib('stimlib0', 'startStimulation');
                    calllib('stimlib0', 'stopStimulation');
                else
                    WaitSecs(predur);
                    Screen('DrawTexture', win, stimvibro3tex, [], mainCue);
                    Screen('Flip', win);
                    logstring(['Cue: ' symbol]);
                    WaitSecs(stimdur);
                end    
            end
            quest = char(strcat(stim_quest, '\n\n\n', ansbalance(rl2(i)+1)));
            DrawFormattedText(win, quest, 'center', 'center', [0 0 0]);
            Screen(win,'flip');
            logstring('Waiting for response');
            WaitSecs(ansdur);
        elseif strcmp(sequence.Condition(i), 'imag')
            log_ExPra19.Design(3,i) = 2;
            if strcmp(sequence.Stimulation(i), 'vibro1')
                log_ExPra19.Design(4,i) = 1;
                symbol = [imagcolour '_' vibro1cue];
                logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
                while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
                end
                disp(['Pre-cue onset']);
                Screen('DrawTexture', win, imagvibro1tex, [], preCue);
                Screen('Flip', win);
                logstring(['Pre-cue: ' symbol]);
                tic
                WaitSecs(predur);
                Screen('DrawTexture', win, imagvibro1tex, [], mainCue);
                Screen('Flip', win);
                logstring(['Cue: ' symbol]);
                WaitSecs(stimdur);
            elseif strcmp(sequence.Stimulation(i), 'vibro2')
                log_ExPra19.Design(4,i) = 2;
                symbol = [imagcolour '_' vibro2cue];
                logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
                while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
                end
                disp(['Pre-cue onset']);
                Screen('DrawTexture', win, imagvibro2tex, [], preCue);
                Screen('Flip', win);
                logstring(['Pre-cue: ' symbol]);
                tic
                WaitSecs(predur);
                Screen('DrawTexture', win, imagvibro2tex, [], mainCue);
                Screen('Flip', win);
                logstring(['Cue: ' symbol]);
                WaitSecs(stimdur);
            elseif strcmp(sequence.Stimulation(i), 'vibro3')
                log_ExPra19.Design(4,i) = 3;
                symbol = [imagcolour '_' vibro3cue];
                logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
                while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
                end
                disp(['Pre-cue onset']);
                Screen('DrawTexture', win, imagvibro3tex, [], preCue);
                Screen('Flip', win);
                logstring(['Pre-cue: ' symbol]);
                tic
                WaitSecs(predur);
                Screen('DrawTexture', win, imagvibro3tex, [], mainCue);
                Screen('Flip', win);
                logstring(['Cue: ' symbol]);
                WaitSecs(stimdur);
            end
            quest = char(strcat(imag_quest, '\n\n\n', ansbalance(rl2(i)+1)));
            DrawFormattedText(win, quest, 'center', 'center', [0 0 0]);
            Screen(win,'flip');
            logstring('Waiting for response');
            WaitSecs(ansdur);
        elseif strcmp(sequence.Condition(i), 'att')
            log_ExPra19.Design(3,i) = 3;
            if strcmp(sequence.Stimulation(i), 'vibro1')
                log_ExPra19.Design(4,i) = 1;
                symbol = [stimcolour '_' vibro1cue];
                logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
                while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
                end
                disp(['Pre-cue onset']);
                Screen('DrawTexture', win, stimvibro1tex, [], preCue);
                Screen('Flip', win);
                logstring(['Pre-cue: ' symbol]);
                tic
                WaitSecs(predur);
                Screen('DrawTexture', win, stimvibro1tex, [], mainCue);
                Screen('Flip', win);
                logstring(['Cue: ' symbol]);
                WaitSecs(stimdur);
            elseif strcmp(sequence.Stimulation(i), 'vibro2')
                log_ExPra19.Design(4,i) = 2;
                symbol = [stimcolour '_' vibro2cue];
                logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
                while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
                end
                disp(['Pre-cue onset']);
                Screen('DrawTexture', win, stimvibro2tex, [], preCue);
                Screen('Flip', win);
                logstring(['Pre-cue: ' symbol]);
                tic
                WaitSecs(predur);
                Screen('DrawTexture', win, stimvibro2tex, [], mainCue);
                Screen('Flip', win);
                logstring(['Cue: ' symbol]);
                WaitSecs(stimdur);
            elseif strcmp(sequence.Stimulation(i), 'vibro3')
                log_ExPra19.Design(4,i) = 3;
                symbol = [stimcolour '_' vibro3cue];
                logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
                while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
                end
                disp(['Pre-cue onset']);
                Screen('DrawTexture', win, stimvibro3tex, [], preCue);
                Screen('Flip', win);
                logstring(['Pre-cue: ' symbol]);
                tic
                WaitSecs(predur);
                Screen('DrawTexture', win, stimvibro3tex, [], mainCue);
                Screen('Flip', win);
                logstring(['Cue: ' symbol]);
                WaitSecs(stimdur);
            end
            quest = char(strcat(stim_quest, '\n\n\n', ansbalance(rl2(i)+1)));
            DrawFormattedText(win, quest, 'center', 'center', [0 0 0]);
            Screen(win,'flip');
            logstring('Waiting for response');
            WaitSecs(ansdur);
        elseif strcmp(sequence.Condition(i), 'null')
            logstring(['Wait for ActionTime  =' num2str(log_ExPra19.Design(1,i))]);
            while GetSecs - t_first < actiontime
                    GetSecs;
                    [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
                    thisResp = KbName(lastPress);
                    if iscell(thisResp)
                        thisResp = thisResp{1};
                    end
                    if strcmp(thisResp, 'esc')
                        sca;
                        break;
                    elseif strcmp(thisResp, 'q') && m == 1
                        break;
                    end
            end
            disp(['Pre-cue onset']);
            symbol = ['cross.png'];
            Screen('DrawTexture', win, crosstex, [], mainCue);
            Screen('Flip', win);
            logstring(['Cue: ' symbol]);
            WaitSecs(predur);
            WaitSecs(stimdur+ansdur);
        end

        % Write out rating
        [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
        thisResp = KbName(lastPress);
        if iscell(thisResp)
            thisResp = thisResp{1};
        end
        if strcmp(thisResp, 'esc')
            sca;
            break;
        elseif strcmp(thisResp, 'q') && m == 1
            break;
        end

        if strcmp(thisResp, char(anskeys(rl2(i)+1)))
            log_ExPra19.responses(1,i) = 2;
        elseif strcmp(thisResp, char(anskeys(2-rl2(i))))
            log_ExPra19.responses(1,i) = 1;
        end
        if log_ExPra19.responses(1,i) == 2 && strcmp(sequence.Condition(i), 'stim')
            log_ExPra19.responses(2,i) = 2;
        elseif log_ExPra19.responses(1,i) == 1 && strcmp(sequence.Condition(i), 'att')
            log_ExPra19.responses(2,i) = 2;
        elseif log_ExPra19.responses(1,i) == 2 && strcmp(sequence.Condition(i), 'att')
            log_ExPra19.responses(2,i) = 1;
        elseif log_ExPra19.responses(1,i) == 1 && strcmp(sequence.Condition(i), 'stim')
            log_ExPra19.responses(2,i) = 1;
        end

        logstring(['Pressed: ', thisResp, ' - Meaning: ', num2str(log_ExPra19.responses(1,i)), ' - Correct: ', num2str(log_ExPra19.responses(2,i))]);
        save(log_path, 'log_ExPra19');
        
        if stimul
            % QuaeroSys Check
            if checkQuaeroRemote
                logstring('QuaeroReset'); 
            end
        end
    end


% =========================================================================
% Save log file with run stats, show pause screen
        
    rundur = GetSecs - t_first;
    log_ExPra19.QuickStats.Missed_Trials = sum(log_ExPra19.responses(1,:)==0);
    log_ExPra19.QuickStats.rundur = rundur;
    save(log_path, 'log_ExPra19');
    stop_cogent;
    disp(['Missed trials: ' num2str(sum(log_ExPra19.responses(1,:)==0))]);
    disp(['Correct perc trials: ' num2str(sum(log_ExPra19.responses(2,:)==2))]);

    if m < numruns
        Screen('DrawTexture', win, crosstex, [], mainCue);
        Screen('Flip', win);
        WaitSecs(4);
        DrawFormattedText(win, paus, 'center', 'center', [0 0 0]);
        Screen(win,'flip');
        logstring('PAUSE: waiting for key press');
        KbQueueWait(-1);
        [pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(-1);
        thisResp = KbName(lastPress);
        if iscell(thisResp)
            thisResp = thisResp{1};
        end
        if strcmp(thisResp, 'esc')
            sca;
            break;
        end
        DrawFormattedText(win, '', 'center', 'center', [0 0 0]);
        Screen(win,'flip');
        WaitSecs(0.5);
    end
end

% =========================================================================
% Show goodbye screen

DrawFormattedText(win, bye, 'center', 'center', [0 0 0]);
Screen(win,'flip');
WaitSecs(5);
pause(0.4); 
sca
KbQueueStop;

fclose('all');

%CLOSE COGENT
% cgshut;


