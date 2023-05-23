clear all

addpath(fullfile("functions"))

if(~exist("results","dir"))
    mkdir results
end

%% specifiy your variables here

% specifiy path to dataset
path = "C:\Users\cschr\Desktop\radar_dataset_small_0311";

% choose all IDs
%ID_list = get_ID_list();
% or choose an ID from ID_list
ID_list = "28";

W_MOE = 600; % number of elements (max 1096) inside of the window for Model Order Estimation
W_MUSIC = 700; % number of elements (max 1096) inside of the window for MUSIC Algorithm 
st_windowlength = 200; % slow time window length in samples
N_cov = 10; % number of calculated covariance matrices per window
alpha_RD = 3; % tuning parameter for RD
use_ref_model_order = false; % use reference model order / number of persons in the scene
MODEL_ORDER_MUSIC = 15; % model order used for MUSIC algorithm

% turn plots on
plot_range_profile = false;
plot_on = true;
plot_MOE = false;
plot_RD = false;
plot_spatial = false;

%% radar constants
c = physconst('LightSpeed');
delta = 0.02; % antenna spacing
Bw = 1.7e9; % bandwidth
f0 = 6.3e9; % start frequency

%% main
num_runs = length(ID_list);

for run = 1:num_runs

    %% load radar data
    data1 = load(fullfile(path, "m" + ID_list(run) + "_radar.mat"));
    data_radar = data1.data_radar;
    meta_data_radar = data1.meta_data_radar;
    clear data1
    
    % get data which belongs to rx antennas 2, 6, 10 and 14
    % and concatenate data for tx antenna 1 and 17 to create virtual antenna array
    S = permute(cat(3,data_radar.tx_1(:,:,[1,3,5,7]), data_radar.tx_17(:,:,[1,3,5,7])), [1,3,2]);
    clear data_radar
    
    time_axis_st = meta_data_radar.time_axis_st;
    
    %fs_st = meta_data_radar.fs_st_median; % slow time sampling frequency
    [K, M, L] = size(S); % K: fast time, M: antennas, L: slow time
    
    k = 0:K-1; % number of requency steps
    deltaf = Bw/K; % frequency step
    fc = (K-1)/2*deltaf+f0; % center frequency
    
    %% load reference data
    
    if(exist(fullfile(path, "m" + ID_list(run) + "_reference.mat"), "file"))
        data2 = load(fullfile(path, "m" + ID_list(run) + "_reference.mat"));
        data_ref = data2.data_ref;
        meta_data_ref = data2.meta_data_ref;
        clear data2
        
        ref_timeaxis_st = meta_data_ref.time_axis;
        fs_ref = meta_data_ref.fs_ref;
    
        ref_pulse = data_ref.ecg;
        ref_respiration = data_ref.respiration - movmean(data_ref.respiration,5*fs_ref,1);
    
        ref_range = meta_data_ref.range;
        ref_doa = meta_data_ref.doa;
        person_ID = meta_data_ref.person_ID;
    
        ref_P = length(person_ID);

        obstacle = meta_data_ref.obstacle;
    
        clear data_ref    
    else
        ref_P = 0;
        ref_pulse = [];
        ref_respiration = [];
    
        ref_range = [];
        ref_doa = [];
        person_ID = [];

        ref_param = get_reference_parameter(ID_list(run));
        obstacle = ref_param.obstacle;

        clear ref_param
    end
    
    %% calculate 2D MUSIC spectrum
   
    % grid size of MUSIC spectrum
    delta_d = 0.05;
    delta_theta = 0.05;
    max_d = 4.5;
    max_theta = 0.45*pi;
    d = 0:delta_d:max_d;
    theta = -max_theta:delta_theta:max_theta;
    
    st_windowlength = 200; % slow time window length in samples

    if(ID_list(run) == "0" || ID_list(run) == "0-1" || ID_list(run) == "0-2") % these measurements only have a length of 200
        st_windowlength = 100;
    end
    
    num_windows = floor(L/st_windowlength);
    
    
    P_MUSIC = zeros(length(d), length(theta), num_windows);
    P_MUSIC_sum = zeros(length(d), length(theta), num_windows);
    
    est_P_music = nan(num_windows,1);
    est_P = nan(num_windows,1);
    est_P_med = nan(num_windows,1);
    fs_st = nan(num_windows,1);
    
    S_filt = zeros(K,M,st_windowlength,num_windows);
    
    for win = 1:num_windows
    
        wind = (win-1)*st_windowlength+1:((win-1)*st_windowlength+st_windowlength);
        S_sec = S(:,:,wind);
        time_axis_st_sec = meta_data_radar.time_axis_st(wind);
    
        fs_st(win) = 1/median(diff(time_axis_st_sec)); %slow time sampling frequency for every window
    
        %% clutter removal
        Prefilter_length = 5; % seconds
        S_filt(:,:,:,win) = S_sec - movmean(S_sec,Prefilter_length*fs_st(win),3);

        % alternative SVD clutter removal
        %S_filt(:,:,:,win) = clutter_SVD(S_sec, 1, 3);

       
        %% plot range profile over slow time
        if plot_range_profile
            ant_range = 1; % exemplary range profile for antenna 1
            % the range profile is obtained by taking the ifft along fast time of S

            [range_profile, range_vec] = get_Range_Profile(squeeze(S_filt(:,ant_range,:,win)),deltaf,K);
    
            figure
            imagesc(time_axis_st_sec, range_vec, range_profile)
            set(gca,'YDir','normal')
            for i = 1:ref_P
                yline(ref_range(i),'r-',"P" + num2str(person_ID(i)))
            end
            xlabel('Slow Time (s)')
            ylabel('Distance (m)')
            ylim([0 4])
            title("Range Profile")
        
            % a range profile with finer granularity can be achieved by zero padding N > K
            [range_profile, range_vec] = get_Range_Profile(squeeze(S_filt(:,ant_range,:,win)),deltaf,8182);

            figure
            imagesc(time_axis_st_sec, range_vec, range_profile)
            set(gca,'YDir','normal')
            for i = 1:ref_P
                yline(ref_range(i),'r-',"P" + num2str(person_ID(i)))
            end
            xlabel('Slow Time (s)')
            ylabel('Distance (m)')
            ylim([0 4])
            title("HRRP")
        end
        
        %% 2D MUSIC with constant model order
    
        %S_mean = mean(S_filt,3);
        %S_mean = S_filt(:,:,1,win);
    
        for i = 1:N_cov
            [R_mat(:,:,i), W_K, W_M, W_KM] = SpacialSmoothing(S_filt(:,:,1+(i-1)*floor(st_windowlength/N_cov),win), W_MUSIC);
        end
        R_tilde = mean(R_mat,3);
        clear R_mat
    
        % forward backward averaging
        J = fliplr(eye(W_K*W_M));
        R = (R_tilde + J * R_tilde' * J) ./ 2;
        
        [V, EV] = svd(R);
        
        if use_ref_model_order
            est_P_music(win) = ref_P;
        else
            est_P_music(win) = MODEL_ORDER_MUSIC;
        end
    
        % steering vector
        a_steer = @(d,theta) exp(-1i*2*pi*(f0+k(1:W_K)*deltaf).'*1/c*(2.*d+sin(theta)*delta*(0:W_M-1)));
                    
        P_MUSIC(:,:,win) = Music2DSpectrum(V, est_P_music(win), [delta_d, delta_theta, max_d, max_theta], a_steer);
    
        if(win == 1)
            P_MUSIC_sum(:,:,win) = P_MUSIC(:,:,win);
        else
            P_MUSIC_sum(:,:,win) = P_MUSIC_sum(:,:,win-1) + P_MUSIC(:,:,win);
        end
      
    
        %% estimate number of persons
    
        if W_MOE ~= W_MUSIC
            for i = 1:N_cov
                [R_mat(:,:,i), W_K, W_M, W_KM] = SpacialSmoothing(S_filt(:,:,1+(i-1)*floor(st_windowlength/N_cov),win), W_MOE);
            end
            R_tilde = mean(R_mat,3);
            clear R_mat
        
            % forward backward averaging
            J = fliplr(eye(W_K*W_M));
            R = (R_tilde + J * R_tilde' * J) ./ 2;
    
            [V, EV] = svd(R);
        end
    
        lambda = diag(EV); % eigenvalues of R
        
        D = min(2*W_KM, W_M*W_K); % position of drop
    
        est_P(win) = moe_RD(lambda(1:D), alpha_RD, plot_RD);
    
        est_P_med(win) = floor(median(est_P, "omitnan"));
    
        if(isnan(est_P_med(win)))
            est_P_med(win) = 0;
        end
    
        if plot_MOE
            figure
            plot(log10(lambda(1:30)))
            xline(ref_P,'r-','Ref Model Order')
            xline(est_P(win),'r--','Est Model Order', "LabelVerticalAlignment","bottom")
            xlabel('Sorted Eigenvalues')
            ylabel('Manitude (dB)')
            grid on
        end
    end
    
    % figure
    % imagesc(rad2deg(theta),d,abs(P_MUSIC_sum(:,:,win)))
    % set(gca,'YDir','normal')
    % hold on
    % colorbar
    % title("Sum " + num2str(win))
    
        
    %% Location Estimates
    
    target_group_distance = 0.3; % Maximum distance in m beteen peaks to belong to the same person
    
    d_est = cell(1, num_windows);
    theta_est = cell(1, num_windows);
    
    d_est_group = cell(1, num_windows);
    theta_est_group = cell(1, num_windows);
    
    d_est_group_only = cell(1, num_windows);
    theta_est_group_only = cell(1, num_windows);
    
    for win = 1:num_windows
        % finding peaks on 2D Music Spectrum
        Pmusic = P_MUSIC_sum(:,:,win);
        Pmusic_T = Pmusic.';
        [~, Loc] = findpeaks(Pmusic(:));
        [~, Loc_T] = findpeaks(Pmusic_T(:));
        [Col, Row] = ind2sub(size(Pmusic),Loc);
        [Ind] = sub2ind(size(Pmusic_T),Row,Col);
        Ind_intersect = intersect(Ind,Loc_T);
        [theta_index, d_index] = ind2sub(size(Pmusic_T), Ind_intersect);
        
        % sorting peaks with respect to height
        Pmusic_value = zeros(1,length(theta_index));
        for i = 1:length(theta_index)
            Pmusic_value(i) = Pmusic(d_index(i),theta_index(i));
        end
        [~,Index_sorted] = sort(Pmusic_value,'descend');
        theta_p{win} = theta(theta_index(Index_sorted.')).';
        d_p{win} = d(d_index(Index_sorted.')).';

%         figure
%         imagesc(rad2deg(theta),d,abs(P_MUSIC_sum(:,:,win)))
%         set(gca,'YDir','normal')
%         hold on
%         scatter(rad2deg(theta_p{win}), d_p{win},"w")
    
        % less peaks found as persons estimated, set estimated # of persons equal found peaks
        if(est_P_med(win) > length(theta_p{:,win}))
            est_P_med(win) = length(theta_p{:,win});
        end

        % Grouping a minimal number of higest peaks to form as many groups as the estimated model order.
        for P_test = est_P_med(win):length(theta_p{:,win})
            % identifying peaks that are located close to each other
            theta_test = theta_p{win}(1:P_test);
            d_test = d_p{win}(1:P_test);
            [X,Y] = pol2cart(theta_test,d_test);
            closeLoc = [];
            for nn = 1:P_test
                for jj = 1:P_test
                    if sqrt((X(nn)-X(jj))^2 + (Y(nn)-Y(jj))^2) < target_group_distance
    
                        if(nn == 1 && jj == 1)
                            closeLoc = [closeLoc; nn, jj];
                        end
    
                        % dismiss peaks that are already close to a higher peak
                        if(isempty(intersect(closeLoc(:,2), jj)))
                            closeLoc = [closeLoc; nn, jj];
                        end
                    end
                end
            end
            % sorting peaks into the same groups if they are close to each other
            PersonGroup = {};
            for i = 1:size(closeLoc,1)
                found = false;
                for j = 1:length(PersonGroup)
                    % adding peak to an existing group
                    if intersect(PersonGroup{j},closeLoc(i,:))
                        PersonGroup{j} = [PersonGroup{j} closeLoc(i,:)];
                        found = true;
                        break;
                    end
                end
                % adding a new group for a peak that is not near another group
                if ~found
                    PersonGroup{end+1} = closeLoc(i,:);
                end
            end
            % The loop is repeated with an aditional peak if there are less groups than the estimated model order.
            if length(PersonGroup)>=est_P_med(win) || length(PersonGroup)>=5
                break;
            end
        end
    
        d_est_group_only{win} = cell(length(PersonGroup),1);
        theta_est_group_only{win} = cell(length(PersonGroup),1);    
        for j = 1:length(PersonGroup)
            % The indeces in each group contain duplicates that need to be filtered out. 
            PersonGroup{j} = unique(PersonGroup{j});
            [X,Y] = pol2cart(theta_p{win}(PersonGroup{j}),d_p{win}(PersonGroup{j}));
            X = mean(X);
            Y = mean(Y);
            [theta_est_temp, d_est_temp] = cart2pol(X,Y);
            d_est{win}(j,:) = d_est_temp;
            theta_est{win}(j,:) = theta_est_temp;
    
            d_est_group{win}{j,:} = d_p{win}(PersonGroup{j}).';
            theta_est_group{win}{j,:} = theta_p{win}(PersonGroup{j}).';
    
            if(length(d_est_group{win}{j,:}) > 1)
                d_est_group_only{win}{j,:} = d_p{win}(PersonGroup{j}).';
                theta_est_group_only{win}{j,:} = theta_p{win}(PersonGroup{j}).';
            end
        end
    
        est_P_med_group(win) = length(d_est{win});
    end
        
    %% Spacial Filtering
    K_max = 90; % floor(c/(2*delta*deltaf)-f0/deltaf)

    k_a = 0:K_max-1;
    a = @(d,theta) exp(-1i*2*pi*(f0+k_a*deltaf).'*1/c*(2.*d+sin(theta)*delta*([0:M-1])));
    
    for win = 1:num_windows
        w{win} = get_spatial_filter(a,d_est{win},theta_est{win});
        weight_mat = blackman(K_max)*blackman(M).';
        w{win} = weight_mat(:).*w{win};
    end
    
    for win = num_windows
        filter_response = [];
        for j = 1:size(w{win},2)
            for i = 1:length(theta)
                for ii = 1:length(d)
                    a_test = a(d(ii),theta(i));
                    a_test = a_test(1:K_max,:);
                    filter_response(ii,i,j) = w{win}(:,j)'*a_test(:);
                end
            end
        end
        if plot_spatial
            figure
            hold on
            subplot(2,1,1)
            plot(squeeze(sum(abs(filter_response).^2,1)))
            subplot(2,1,2)
            plot(squeeze(sum(abs(filter_response).^2,2)))
        end
    end
    
    %% Tracking
    %track estimated locations
    
    tracking_distance = 0.2;
    
    clear track
    for win = 1:num_windows
        track{win} = [];
        if(win == 1)
            track{win} = 1:est_P_med_group(win);
        else
            for i = 1:est_P_med_group(win)
                [X_this,Y_this] = pol2cart(theta_est{win}(i),d_est{win}(i));
    
                for win_prev = 1:win-1
    
                    for j = 1:est_P_med_group(win_prev)
                        [X_prev,Y_prev] = pol2cart(theta_est{win_prev}(j), d_est{win_prev}(j));
                
                        err = sqrt((X_prev - X_this)^2+(Y_prev - Y_this)^2);
        
                        if err < tracking_distance
                            track{win}(i) = track{win_prev}(j);
                        end
                    end
                end
    
                if(length(track{win}) < i)
                    %catch if track is empty
                    if(isempty(track{win}))
                         track_temp = track;
                         track_temp{win} = nan;
                         temp = max(cellfun(@(x) max(x), track_temp)) + 1;
                    else
                        temp = max(cellfun(@(x) max(x), track)) + 1;
                    end
                    track{win} = [track{win} temp];
                    clear track_temp temp
                end
            end
    
            if(length(track{win}) < est_P_med_group(win))
                track{win} = [track{win} (max(cellfun(@(x) max(x), track)) + [1:(est_P_med_group(win) - length(track{win}))])];
            end
        end
        if(isempty(track{win}))
            track{win} = 0;
        end
    end
    
    if plot_on 
        numSub = numSubplots(num_windows);
        figure
        for win = 1:num_windows
            % data for circles around reference positions
            [th1, r1] = get_Ref_Circles(ref_doa, ref_range, target_group_distance, ref_P);
            
            subplot(numSub(1), numSub(2),win)
    
            imagesc(rad2deg(theta),d,abs(P_MUSIC_sum(:,:,win)))
            set(gca,'YDir','normal')
            hold on
    
            for i = 1:ref_P
                text(ref_doa(i),ref_range(i), "P" + num2str(person_ID(i)),"Color",'w', 'HorizontalAlignment','center')       
            end
            plot(th1*180/pi, r1,'white','LineWidth',1);
        
            cc = ["o", "+", "*", "x", "^", "v", ">", "<", ".", "square", "diamond", "pentagram", "hexagram","o", "+", "*", "x", "^", "v", ">", "<", "."];
            for i = 1:est_P_med_group(win)
                plot(rad2deg(theta_est_group_only{win}{i}),d_est_group_only{win}{i},cc(i),'LineWidth',1,'Color','g')
                text(rad2deg(theta_est{win}(i)), d_est{win}(i), num2str(track{win}(i)),"Color",'r', 'HorizontalAlignment','center')  
            end
            xlabel('DoA (deg)')
            ylabel('Range (m)')
            title("M" + ID_list(run) + ", W" + num2str(win))
        end


        %% save data for polar plot 
        ii = 1;
        for i = 1:length(d)
            for j = 1:length(theta)
                plot_data(ii,:) = [d(i) 180-(rad2deg(theta(j))+90) P_MUSIC_sum(i,j,num_windows)];
                ii = ii + 1;
            end
        end

        table_titles = ["radius", "angle", "value"];

        T = array2table(plot_data);
        T.Properties.VariableNames = table_titles;
        writetable(T, fullfile("results", "musicSpectrum-ID" + num2str(ID_list(run))+ ".csv"))

    end
    
    
    
    %% Vital Sign Estimation
    
    for win = 1:num_windows
        if(est_P_med_group(win) > 0)
            %y = zeros(slowtime_windowlength,P_targets);
            vecS = zeros(M*K_max,st_windowlength);
            
            for i = 1:st_windowlength
                temp = squeeze(S_filt(1:K_max,:,i,win));
                vecS = temp(:);
                for p = 1:est_P_med_group(win)
                    y{win}(i,p) = w{win}(:,p)'*vecS;
                end
            end


            % phase unwrapping based on DACM algorithm from https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6583987
            % seems to work, but not a large difference
%             eta_dacm = zeros(st_windowlength-1, est_P_med_group(win));
%             for i = 2:st_windowlength
%                 eta_dacm_temp = zeros(1, est_P_med_group(win));
%                 for ii = 2:i
%                     eta_dacm_temp = eta_dacm_temp + (real(y{win}(ii,:)).*(imag(y{win}(ii,:)) - imag(y{win}(ii-1,:))) - ... 
%                                     (real(y{win}(ii,:)) - real(y{win}(ii-1,:))).*imag(y{win}(ii,:))) ...
%                                     ./ (real(y{win}(ii,:)).^2 + imag(y{win}(ii,:)).^2);
%                 end
%                 eta_dacm(i,:) = -c/4/pi/f0*eta_dacm_temp;
%             end
%             eta{win} = eta_dacm - movmean(eta_dacm,50);

            % dumb approach
            eta_temp = -c/4/pi/f0*unwrap(angle(y{win}));
            eta{win} = eta_temp - movmean(eta_temp,50);
        else
            eta{win} = zeros(st_windowlength,1);
        end
    end
    
    %% plot VS for all estimated persons
    
    P_vs = max(cellfun(@(x) max(x), track));
    
    eta_persons = zeros(st_windowlength,num_windows, P_vs);
    
    for j = 1:P_vs
        j_found = cellfun(@(x) x == j, track, 'UniformOutput', false);
        for win = 1:num_windows
            if(sum(j_found{win}) > 0)
                eta_persons(:,win,j) = mean(eta{win}(:,j_found{win}), 2);
            end
        end
    end
    
    vs_t = nan(num_windows, st_windowlength);
    for win = 1:num_windows
        wind = (win-1)*st_windowlength+1:((win-1)*st_windowlength+st_windowlength);
        vs_t(win,:) = time_axis_st(wind);
    end
    
    if plot_on
        numSub = numSubplots(P_vs);
        figure
         for j = 1:P_vs
            subplot(numSub(1), numSub(2), j)
        
            plot(vs_t.', eta_persons(:,:,j))
            title("est-P" + num2str(j))
            xlim([0 max(vs_t,[],"all")])
            grid on
            xlabel("Time in seconds")
         end
         sgtitle("VS of estimated persons")
    end
    
    
    %% frequency estimation for all estimated persons
    
    N_fft_est = 2048;
    wind_fft = blackman(st_windowlength);
    A = sum(abs(wind_fft./st_windowlength).^2);
    
    clear freq_vector
    % sampling frequency per window can fluctuate, hence freq vector for every window is different
    for win = 1:num_windows
        freq_vector(:,win) = (-N_fft_est/2:N_fft_est/2-1)*fs_st(win)/N_fft_est;
    end
    
    I = abs(fftshift(fft(wind_fft.*eta_persons,N_fft_est,1),1)).^2 ./ N_fft_est ./ A;
    
    I = I(N_fft_est/2+1:end,:,:);
    freq_vector = freq_vector(N_fft_est/2+1:end,:);
    
    % interpolate all periodograms to be on same grid
    freq_all = mean(freq_vector(:,max(freq_vector) == min(max(freq_vector))), 2); % mean catch same minmax
    I_interp = nan(N_fft_est/2, num_windows, P_vs);
    for win = 1:num_windows
        for j = 1:P_vs
            I_interp(:,win,j) = interp1(freq_vector(:,win), I(:,win,j), freq_all);
        end
    end
    
    I_mean = squeeze(mean(I_interp,2));
    
    % frequency peak based on individual spectra
    [~, freq_locs] = max(I_interp, [], 1);
    
    est_respiration_freq = squeeze(freq_all(freq_locs));
    est_respiration_freq(est_respiration_freq == 0) = nan;
    
    est_respiration_freq_persons = median(est_respiration_freq, "omitnan");
    
    if plot_on
        numSub = numSubplots(P_vs);
        figure
        for j = 1:P_vs
            subplot(numSub(1), numSub(2), j)
            imagesc(1:num_windows,freq_all,I_interp(:,:,j))
            set(gca,'YDir','normal')
            hold on
            plot(est_respiration_freq(:,j), "xr")
            plot(est_respiration_freq_persons(j)*ones(num_windows,1), "g")
            %line([est_respiration_freq_persons(j) est_respiration_freq_persons(j)], [0 max(I(:,:,j),[],"all")], "Linewidth", 1.5)
            ylim([0 2])
            title("est-P" + num2str(j))
            grid on
            xlabel("windows")
            ylabel("Frequency in Hz")
        end
        sgtitle("Spectrum of VS of estimated persons")
        legend("Track max of Estimated", "Median of Estimated")
        
        
        numSub = numSubplots(P_vs);
        figure
        for j = 1:P_vs
            subplot(numSub(1), numSub(2), j)
            semilogy(freq_all*60, I_mean(:,j))
            hold on
            title("est-P" + num2str(j))
            xlim([0 100])
            xticks([0 20 40 60 80 100])
            grid on
            xlabel("frequency (acts/min)")
        end
        sgtitle("Mean-Spectrum over all windows of VS of estimated persons")
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Statring from here reference parameters are used 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Matching
    % Match estimated locations with reference locations
    
    loc_error = NaN(ref_P, num_windows);
    
    match = cell(1, num_windows);
    match_stats = cell(1, num_windows);

    for win = 1:num_windows
    
        match{win} = - ones(length(d_est{win}), 1);
        match_stats{win} = - ones(length(d_est{win}), 1);
        
        for i = 1:length(d_est{win})
        
            [X_est,Y_est] = pol2cart(theta_est{win}(i),d_est{win}(i));
        
            for j = 1:ref_P       
                [X_ref,Y_ref] = pol2cart(deg2rad(ref_doa(j)),ref_range(j));
        
                err = sqrt((X_ref-X_est)^2+(Y_ref-Y_est)^2);
    
                if (err < target_group_distance && sum(match{win} == j) == 0) % avoid matching two positions to the same person
                    match{win}(i) = j;
                    loc_error(j,win) = err; % size ref_P x win, -1 if no matching location found
                end
                % for stats
                if (err < target_group_distance) % avoid matching two positions to the same person
                    match_stats{win}(i) = j;
                end
            end
        end
    end
    
    missed_detections = ref_P - cellfun(@(x) sum(unique(x) ~= -1), match_stats).';
    false_detections = cellfun(@(x) sum(x == -1), match_stats).';
    loc_error_mean = mean(loc_error, 1, "omitnan").';
    TPP = (cellfun(@(x) sum(unique(x) ~= -1), match_stats)/ref_P).'; % true positive proportion = # detections / ref. # of persons
    FDP = false_detections./est_P_med; % false discovery proportion = # false detections / est. # of persons
    
%     disp("Statistics of last window, NOT averaged over windows!")
%     disp("mean loc error: " + num2str(loc_error_mean(num_windows)))
%     disp("missed detections: " + num2str(missed_detections(num_windows)))
%     disp("false detections: " + num2str(false_detections(num_windows)))
%     disp("TPP: " + num2str(TPP(num_windows)))
%     disp("FDP: " + num2str(FDP(num_windows)))

    %% plot errors
    if plot_on
        if(ref_P  > 0)
            for j = 1:ref_P
                leg(j) = "P" + j;
            end
            
            figure
            plot(loc_error.')
            legend(leg)
            xlabel("windows")
            ylabel("loc error in m")
            ylim([0 target_group_distance])
        end
    end

    %% plot vitals of reference Persons
    vss = zeros(num_windows, st_windowlength, ref_P);
    
    for j = 1:ref_P
        for win = 1:num_windows
            if(sum(match{win} == j) ~= 0)
                vss(win,:,j) = eta{win}(:,match{win} == j);
            end
        end
       
        % correct sign changes between windows
        vs(1,:,j) = vss(1,:,j);
        for win = 2:num_windows
            if(sign(vs(win-1,end,j)) ~= sign(vss(win,1,j)))
                vs(win,:,j) = -vss(win,:,j);
            else
                vs(win,:,j) = vss(win,:,j);
            end
        end
    end
    
    % for win = 1:num_windows
    %     wind = (win-1)*st_windowlength+1:((win-1)*st_windowlength+st_windowlength);
    %     vs_t(win,:) = time_axis_st(wind);
    % end
    
    if plot_on
        numSub = numSubplots(ref_P);
        figure
        for j = 1:ref_P
            subplot(numSub(1), numSub(2), j)
            plot(ref_timeaxis_st, (ref_respiration(:,j))./max(ref_respiration(:,j)), "b")
            hold on
            plot(vs_t.', vs(:,:,j).'./max(vs(:,:,j).'))
            title("ref-P" + num2str(j))
            xlim([0 max(vs_t,[],"all")])
            grid on
        end
        sgtitle("VS of reference persons")
    end

    %% plot spectrum of reference persons
    N_fft_ref = 4096;

    %respiration
    I_ref = nan(N_fft_ref,num_windows,ref_P);
    I_ref_resp_half = nan(N_fft_ref/2,num_windows,ref_P);

    if(ref_P > 0) 
        % downsample by factor 10
        down_factor = 10;
        ref_respiration_down = ref_respiration(1:down_factor:end,:);
        N_ref = length(ref_respiration_down);
        
        ref_timeaxis_st_down = ref_timeaxis_st(1:down_factor:end,:);
       
        freq_vector_resp_ref = (-N_fft_ref/2:N_fft_ref/2-1)*fs_ref/down_factor/N_fft_ref;
        
        for win = 1:num_windows
        
            ref_respiration_down_win{win} = ref_respiration_down(ref_timeaxis_st_down <= vs_t(win,end) & ref_timeaxis_st_down >= vs_t(win,1),:);
            N_ref = length(ref_respiration_down_win{win});
            wind_fft = blackman(N_ref);
            A = sum(abs(wind_fft./N_ref).^2);
        
            I_ref(:,win,:) = abs(fftshift(fft(wind_fft.*ref_respiration_down_win{win},N_fft_ref,1),1)).^2 ./ N_fft_ref ./ A;
        
        end
        
        I_ref_resp_half = I_ref(N_fft_ref/2+1:end,:,:);
        freq_vector_ref_resp_half = freq_vector_resp_ref(N_fft_ref/2+1:end).';
    end

    % pulse
    I_ref_pulse = nan(N_fft_ref,num_windows,ref_P);
    I_ref_pulse_half = nan(N_fft_ref/2,num_windows,ref_P);

    if(ref_P > 0) 
        % downsample by factor 10
        down_factor = 1;
        ref_pulse_down = ref_pulse(1:down_factor:end,:);
        N_ref = length(ref_pulse_down);
        
        ref_timeaxis_st_down = ref_timeaxis_st(1:down_factor:end,:);
       
        freq_vector_pulse_ref = (-N_fft_ref/2:N_fft_ref/2-1)*fs_ref/down_factor/N_fft_ref;
        
        for win = 1:num_windows
        
            ref_pulse_down_win{win} = ref_pulse_down(ref_timeaxis_st_down <= vs_t(win,end) & ref_timeaxis_st_down >= vs_t(win,1),:);
            N_ref = length(ref_pulse_down_win{win});
            wind_fft = blackman(N_ref);
            A = sum(abs(wind_fft./N_ref).^2);
        
            I_ref_pulse(:,win,:) = abs(fftshift(fft(wind_fft.*ref_pulse_down_win{win},N_fft_ref,1),1)).^2 ./ N_fft_ref ./ A;
        
        end
        
        I_ref_pulse_half = I_ref_pulse(N_fft_ref/2+1:end,:,:);
        freq_vector_ref_pulse_half = freq_vector_pulse_ref(N_fft_ref/2+1:end).';
    end



    % match estimated and reference signal
    est_respiration_freq_P_ref = nan(num_windows, ref_P);
    I_interp_P_ref = nan(N_fft_est/2, num_windows, ref_P);
    for j = 1:ref_P
        for win = 1:num_windows
            if(sum((match{win} == j).' .* track{win}) ~= 0)
                est_respiration_freq_P_ref(:,j) = est_respiration_freq(:, sum((match{win} == j).' .* track{win}));
                I_interp_P_ref(:,:,j) = I_interp(:,:,sum((match{win} == j).' .* track{win}));
            end
        end
    end

    if plot_on
        numSub = numSubplots(ref_P);
        figure
        for j = 1:ref_P
            subplot(numSub(1), numSub(2), j)
            imagesc(1:num_windows,freq_vector_ref_resp_half,I_ref_resp_half(:,:,j))
            set(gca,'YDir','normal')
            hold on
            plot(est_respiration_freq_P_ref(:,j), "xr")
            ylim([0 2])
            title("ref-P" + num2str(j))
            grid on
            xlabel("windows")
            ylabel("Frequency in Hz")
        end
        sgtitle("Spectrum of VS of reference persons")
        legend("Max of estimated persons")


        I_ref_resp_half_mean = squeeze(mean(I_ref_resp_half, 2));
        I_ref_pulse_half_mean = squeeze(mean(I_ref_pulse_half, 2));
        I_interp_P_ref_mean = squeeze(mean(I_interp_P_ref, 2));

        numSub = numSubplots(ref_P);
        figure
        for j = 1:ref_P
            subplot(numSub(1), numSub(2), j)
            semilogy(freq_all*60, I_interp_P_ref_mean(:,j)./max(I_interp_P_ref_mean(:,j)))
            hold on
            semilogy(freq_vector_ref_resp_half*60, I_ref_resp_half_mean(:,j)./max(I_ref_resp_half_mean(:,j)))
            semilogy(freq_vector_ref_pulse_half*60, I_ref_pulse_half_mean(:,j)./max(I_ref_pulse_half_mean(:,j)))
            title("P" + num2str(j))
            xlim([0 120])
            xticks([0 20 40 60 80 100 120])
            grid on
            xlabel("frequency (acts/min)")
        end
        legend("est", "ref-breath", "ref-pulse", "Location","southwest")
        sgtitle("Mean Periodogram over all windows")
    end
    
    %% calculate error in breathing frequency

    if(~isempty(I_ref_resp_half))
        [~, freq_locs] = max(I_ref_resp_half, [], 1);
        
        ref_respiration_freq = squeeze(freq_vector_ref_resp_half(freq_locs));
        ref_respiration_freq(ref_respiration_freq == 0) = nan;
    
        resp_error_person_abs = sqrt(median((est_respiration_freq_P_ref - ref_respiration_freq).^2,1, "omitnan"));
        resp_error_abs = mean(resp_error_person_abs, "omitnan");
    
        resp_error_person_rel = sqrt(median(((est_respiration_freq_P_ref - ref_respiration_freq)./ref_respiration_freq).^2,1, "omitnan"));
        resp_error_rel = mean(resp_error_person_rel, "omitnan");

        resp_error_person_abs_abs = median(est_respiration_freq_P_ref - ref_respiration_freq, "omitnan");
        resp_error_person_rel_abs = median((est_respiration_freq_P_ref - ref_respiration_freq)./ref_respiration_freq, "omitnan");

    else
        resp_error_person_abs_abs = nan;
        resp_error_person_rel_abs = nan;
        resp_error_person_abs = nan;
        resp_error_person_rel = nan;
        resp_error_abs = nan;
        resp_error_rel = nan;
    end


    if plot_on
        numSub = numSubplots(ref_P);
        figure
        for j = 1:ref_P
            subplot(numSub(1), numSub(2), j)
            plot(est_respiration_freq_P_ref(:,j))
            hold on
            plot(ref_respiration_freq(:,j))
            ylim([0 1])
            title("ref-P" + num2str(j))
            grid on
            xlabel("windows")
            ylabel("Frequency in Hz")
        end
        sgtitle("Spectrum of VS of reference persons")
        legend("est", "ref")

        %save data for plot
        table_titles = ["x", append(repmat("refP",1,ref_P), string(person_ID)), append(repmat("estP",1,ref_P), string(person_ID))];

        T = array2table([(1:num_windows).', ref_respiration_freq, est_respiration_freq_P_ref]);
        T.Properties.VariableNames = table_titles;
        writetable(T, fullfile("results", "estVSrefBreath-ID" + num2str(ID_list(run))+ ".csv"))
    end
    
    %% create table with results
    
    eval{run, 1} = str2double(replace(ID_list(run), "-", ""));
    eval{run, 2} = mean(loc_error, "all", "omitnan"); % mean location error over all locations and windows
    eval{run, 3} = mean(loc_error(:,end), "omitnan"); % mean location error over all locations and last window
    eval{run, 4} = mean(loc_error, 2, "omitnan");% mean location error per person
    eval{run, 5} = mean(missed_detections); % mean missed detections over all windows
    eval{run, 6} = missed_detections(end); % missed detections of last windows
    eval{run, 7} = mean(false_detections); % mean false detections over all windows
    eval{run, 8} = false_detections(end); % false detections of last windows
    eval{run, 9} = mean(TPP); % mean TPP over all windows
    eval{run, 10} = TPP(end); % TPP of last window
    eval{run, 11} = mean(FDP); % mean FDP over all window
    eval{run, 12} = FDP(end); % FDP of last window
    eval{run, 13} = obstacle; %free space (0), door (1), wall (2)
    eval{run, 14} = length(person_ID); % number of persons
    eval{run, 15} = person_ID; % person IDs
    eval{run, 16} = resp_error_person_abs; % absolute error respiration per person
    eval{run, 17} = resp_error_person_rel; % relative error respiration per person
    eval{run, 18} = resp_error_person_abs_abs; % absolute error respiration per person
    eval{run, 19} = resp_error_person_rel_abs; % relative error respiration per person

    save(fullfile("results","VS_M"+ID_list(run)+".mat"), ...
        "ref_timeaxis_st","ref_respiration","vs_t","vs", ...
        "I_mean","freq_all","I_ref_resp_half","freq_vector_ref_resp_half");

    disp("ID: " + ID_list(run) + ", W_MOE: " + num2str(W_MOE) + ", W_MUSIC: " + num2str(W_MUSIC) + ", alpha_RD: " + num2str(alpha_RD))
end

loc_table_titles = ["ID", "MeanLocError", "LastLocError", "PersonLocError", "MeanMissedDetections", "LastMissedDetections","MeanFalseDetections",...
                "LastFalseDetections", "MeanTPP", "LastTPP", "MeanFDP", "LastFDP", "obstacle", "numPersons", "Persons", ...
                "MedianSquaredErrorRespAbs", "MedianSquaredErrorRespRel", "ErrorRespAbsLast", "ErrorRespRelLast"];

T = array2table(eval);
T.Properties.VariableNames = loc_table_titles;
save(fullfile("results", "eval_numWin-" + num2str(num_windows) + "_Wmoe-" + num2str(W_MOE) + "_Wmusic-"+ num2str(W_MUSIC)...
               + "_Ncov-"+ num2str(N_cov) + "_alphaRD-" + num2str(alpha_RD) + ".mat"), "T")


