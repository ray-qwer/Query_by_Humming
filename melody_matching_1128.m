% load dataset from dataset.mat
dataset = load('dataset.mat');
dataset = dataset.title_melody_pair;
titles = [dataset{:}];
titles = [titles.title];

% load humming record sound
song = input("which song","s");
if song == "n"
    return
end
song = string(song);
title = strsplit(song,"_");
title = title(2);
song = "D:\NTUEE\4_grade\first_semester\topic\djj\TestSongs\"+song +".wav";

% get midi of the humming record
[humming_midi,beats] = getMidiNBeatByHumming(song);

% compare with each song
distOfEachSong = zeros(1,length(dataset));
for i = 1:length(dataset)
    distOfEachSong(i) = getSongMatching(dataset{i}.melody,humming_midi,dataset{i}.beats,beats);
end
[score,index] = min(distOfEachSong);
ori_index = find(titles == title);
ori_score = distOfEachSong(ori_index);
disp(titles(index));
disp(score);
disp(titles(ori_index));
disp(ori_score);
least5 = sort(distOfEachSong);
least5 = least5(1:5)

function score = getSongMatching(t_melody,q_melody,t_beats,q_beats)
    target_diff = getDiffofMIDI(t_melody);
    query_diff = getDiffofMIDI(q_melody);
    target_rel_beats = getRelativeBeats(t_beats);
    query_rel_beats = getRelativeBeats(q_beats);
    lambda = 0.5;
    penalty = 5;
    DP = zeros(length(q_melody),length(t_melody));
    DP(:,1) = [0:length(query_diff)]*penalty;
    for i = 2:length(q_melody)
        for j = 2:length(t_melody)
            diff = abs(target_diff(j-1)-query_diff(i-1));
            diff = getDist(diff);
            beat_diff = getBeatDiff(target_rel_beats(j-1),query_rel_beats(i-1));
            diff = diff + beat_diff * lambda;
            if i ~= length(q_melody)
                DP(i,j) = min([DP(i-1,j)+penalty,DP(i,j-1)+penalty,DP(i-1,j-1)+diff]);
            else
                DP(i,j) = min([DP(i-1,j)+penalty,DP(i,j-1),DP(i-1,j-1)+diff]);
            end
            
        end
    end
    score = DP(end,end);
end
function beat_diff = getBeatDiff(t,q)
    beat_diff = abs(log2(t) - log2(q));
end
function relBeats = getRelativeBeats(beats)
    relBeats = zeros(1,length(beats)-1);
    for i = 1:length(beats)-1
        relBeats(i) = beats(i+1)/beats(i);
    end
end
function dist = getDist(diff)
    b = mod(diff+6, 12)-6;
    a = floor((diff+6)/ 12);
    dist = abs(a) + abs(b);
end
function diff = getDiffofMIDI(target)
    diff = conv(target,[1,-1]);
    diff = diff(2:end-1);
end

function [midi,beats] = getMidiNBeatByHumming(song)
    [x,fs] = audioread(song);
    % global variables: hyperparameters
    slot = round(fs/50); % 20 ms
    shift = round(fs/100); % 10ms
    tau = 0.04; % denoise
    amp_max = 0.8;
    onset_gap = 2; 
    quiet_threshold = 0.4;
    onset_threshold = 4;
    %
    env_A = getSignalEnvelope(x,amp_max,tau,shift,slot);
    onset_set = getOnset(env_A,quiet_threshold,onset_gap,onset_threshold);
    % fft testing one first
    % parameters
    % task 1&2: choose segment
    midi = getMidi(x,fs,onset_set,50,shift);
    beats = getBeat(onset_set);
end
function beats = getBeat(onset_set)
    inteval = zeros(1,length(onset_set)-1);
    for i = 1:length(onset_set)-1
        inteval(i) = onset_set(i+1)-onset_set(i);
    end
    std_beat = median(inteval);
    beats = 2.^(round(log2(inteval/std_beat)));
end
function midi = getMidi(x,fs,onset_set,smoother_width,shift)
    smoother = genSmoother(smoother_width);
    low_freq_threshold = 80;
    base_freq = zeros(1,length(onset_set)-1);
    for i = 1:length(onset_set)-1
        seg = x(onset_set(i)*shift:onset_set(i+1)*shift-1);
        fft_x = fft(seg);
        seg_len = length(seg);
        P2 = abs(fft_x/seg_len);
        P1 = P2(1:floor(seg_len/2)+1);
        P1(2:end-1) = P1(2:end-1)*2;
        low_threshold_of_segment = floor(low_freq_threshold*seg_len/fs);
        P1 = conv(P1,smoother);
        P1 = P1(ceil(length(smoother)/2):end-floor(length(smoother)/2));
        max_freq_magnitude = max(P1(low_threshold_of_segment:end)); % should i start from low_freq_threshold?
        flag = 0;
        for j = low_threshold_of_segment:length(P1)
            if flag == 1
                if P1(j) < P1(j-1) && P1(j-1) > max_freq_magnitude*0.2
                    base_freq(i) = (j-1)*fs/seg_len;
                    break;
                end
            end
            if P1(j) > P1(j-1)
                flag = 1;
            else 
                flag = 0;
            end
        end
    end
    midi = round(48 + 12*log2(base_freq/261.63));
end
function smoother = genSmoother(maxValue)
    if maxValue <=0
        smoother = [0];
        return
    end
    a = [1:maxValue];
    b = [maxValue-1:-1:1];
    smoother = [a,b];
    smoother = smoother/sum(smoother);
end
function gini_diff = giniMethod(segment)
    gini_threshold = 0.1;
    sort_arr = sort(segment,'descend');
    summation = sum(sort_arr);
    Xs(length(sort_arr)) = 0;
    accum_num = 0;
    interval = 1/length(sort_arr);
    for i = 1:length(sort_arr)
        accum_num = accum_num + sort_arr(i);
        Xs(i) = accum_num;
    end
    Xs = Xs / summation;
    gini_t_index = round(gini_threshold*length(sort_arr));
    reimann = sum(Xs(1:gini_t_index));
    reimann_sum = reimann * interval;
    gini_diff = reimann_sum - gini_threshold^2/2;
end
function plot_wave(x,y,order)
    subplot(2,1,order);
    plot(x,y);
    xlim([x(1),x(end)]);
end
function envA = getSignalEnvelope(x,amp_max,tau,shift,slot)
    max_in_signal = max(abs(x));
    % envelope_signal = zeros(1,floor(size(x,1)/shift));
    slot_num = 0;
    envelope = zeros(1,ceil(size(x,1)/shift));
    for index_start = 1:shift:size(x,1)
        slot_num = slot_num + 1;
        index_end = index_start + slot -1;
        if index_end >= size(x,1)
            index_end = size(x,1);
        end
        interval = abs(x(index_start:index_end));
        Ain = max(max(interval)-tau,0);
        Ain = Ain*amp_max/max_in_signal;
        if Ain-tau<0
            Ain = 0;
        end
        envelope(slot_num) = Ain;
    end
    %slot_count = [0:length(envelope)-1];
    % normalize
    mean_of_Ak = mean(envelope);
    envA = envelope/(0.2+0.1*mean_of_Ak);
end
function onset_set = getOnset(env_A,quiet_threshold,onset_gap,onset_threshold) 
    % unchanged variables
    fractional_power = 0.7;
    match_filter = [3,3,4,4,-1,-1,-2,-2,-2,-2,-2,-2];
    
    % take fractional power
    env_A = env_A.^fractional_power;
    % perform convolution with match_filter
    onset_output = conv(env_A, match_filter);
    onset_output = onset_output(1:(end-11));

    flag = 1;
    onset_candidate = 0;
    onset_prev = 0;
    over_threshold_count = 0;
    onset_set = [];
    gap_too_long = false;
    while 1
        for i = 1:length(onset_output) 
            if flag == 2
                if env_A(i) > quiet_threshold
                    over_threshold_count = over_threshold_count + 1;
                    if over_threshold_count > 0.1/0.01
                        if onset_candidate - onset_prev > onset_gap/0.01 && onset_threshold >3
                            gap_too_long = true;
                            break
                        end
                        onset_set(end+1) = onset_candidate;
                        onset_prev = onset_candidate;
                        flag = 0;
                    end
                else
                    flag = 0;
                end
            elseif flag == 0
                flag = 1;
                over_threshold_count = 0;
            end
            if onset_output(i) > onset_threshold
                if flag == 1 && i - onset_prev > 0.2/0.01
                    onset_candidate = i;
                    flag = 2;
                end
            end
        end
        if  gap_too_long
            if onset_threshold <= 3 % onset_threshold is less than denoise
                break
            end
            gap_too_long = false;
            flag = 0;
            over_threshold_count = 0;
            onset_set = [];
            onset_candidate = 0;
            onset_prev = 0;
            onset_threshold = onset_threshold - 0.5;
        else
            break
        end
    end
    onset_set = onset_set - 3;
    onset_set(end+1) = length(onset_output);
end