function Pacticle_Filter_2d
format compact
%% PREDEFINES
% Data
load ('data20171107.mat');
% Control signal
ut = [5,0];
% DEFINE
dT = 0.025;
wB = 4;
% Max moves
Nstep = 625;
% Noise
% Standard deviation of measurement range, bearing
SIGmeas_range = 0.2;                        
SIGmeas_bearing = 2/180*pi;
% Standard deviation of process v,tht
SIGproc_v = 0.5;                   
SIGproc_tht = 3/180*pi;           
SIGproc = [SIGproc_v,SIGproc_tht]; 
% Initial deviation
SIGstate = 10;
% Number of particles (the more, the better the particle
%% MAIN
%% INITIALIZATION
NPARTICLES = 300;
CHI = [SIGstate*rand(NPARTICLES,3),ones(NPARTICLES,1)/NPARTICLES];
% Clear all figure
clf
% lw,sw,aw are used to plot
% and cum_rms is used to calculate rms 
lw = [];%array to hold largest weight particle
cum_rms_lw = 0; %lw cumulative error
sw = [];%array to hold smallest weight particle
cum_rms_sw = 0; %sw cumulative error
aw = [];%array to hold average weight particle
cum_rms_aw = 0; %aw cumulative error
cum_rms_nf = 0; %no filter cumulative error
cum_rms_xodo = 0;%cumulative error between XODO and XTRUE
%% Main loop
for mv = 1:Nstep
%% laser data available at step == 8,16,.. 
    % Count numbers of detected landmark
    NLANDMARKS = countNotNAN(Z(:,:,mv));
    if(NLANDMARKS > 0)
        % Choose the first found landmark
        idx = Z(3,1,mv);
        % Choose the best landmark
        % idx = choose_landmark(CHI,mv,Z,NLANDMARKS,lm);
        
        N = 0; %normalization factor for cumsum function
        CHIbar_t = zeros(NPARTICLES,4);
        CHI_t = zeros(NPARTICLES,4);
        % FIRST SAMPLE and WEIGHT CALCULATE
        for i=1:NPARTICLES
            xtm1 = [CHI(i,1),CHI(i,2),CHI(i,3)];
            xt = sample_motion_model_2d(VG(:,mv),xtm1,SIGproc,dT,wB);
            [dist,b_ang] = xyphi2zPredict(xt,lm(:,idx));
            % Compute weight
            % w1: weight of dist, w2: weight of bearing angle
            w1 = likelihood_calculate(Z(1,1,mv),dist,SIGmeas_range);
            w2 = 1/abs(Z(2,1,mv)-b_ang);
            wt = w1*w2;
            CHIbar_t(i,:) = [xt,wt];
            N = N+wt;
        end
        % Find largest and smallest weight particle
        % and calculate rms error
        [next_lw,next_sw,next_aw] = lsa(CHIbar_t);
        lw = [lw next_lw];
        cum_rms_lw = cal_cumrms(cum_rms_lw,next_lw,XTRUE(:,mv));
        sw = [sw next_sw];
        cum_rms_sw = cal_cumrms(cum_rms_sw,next_sw,XTRUE(:,mv));
        aw = [aw next_aw];
        cum_rms_aw = cal_cumrms(cum_rms_aw,next_aw,XTRUE(:,mv));

        % RESAMPLE 
        cumwt = cumsum(CHIbar_t(:,4))/N;
        for i=1:NPARTICLES
            index = [];
            while(isempty(index))
                swt = rand;
                index = find( cumwt>= swt,1,'first');
            end
            xt = CHIbar_t(index,1:3);
            CHI_t(i,:) = [xt,1/NPARTICLES];
        end
        CHI = CHI_t;
    else
    %% NO DATA 
        for i=1:NPARTICLES
            CHI_t = zeros(NPARTICLES,4);
            xtm1 = [CHI(i,1),CHI(i,2),CHI(i,3)];
            xt = sample_motion_model_2d(VG(:,mv),xtm1,SIGproc,dT,wB);
            CHI_t(i,:) = [xt,1/NPARTICLES];
        end
        CHI = CHI_t;
    end
    % Start point
    if(mv == 8)
        plot(CHI(1,1),CHI(1,2),'hr','LineWidth',3);
        hold on;
    % End point
    elseif(mv == Nstep-1)
        plot(CHI(1,1),CHI(1,2),'hy','LineWidth',3);
    end
    cum_rms_xodo = cal_cumrms(cum_rms_xodo,XODO(:,mv),XTRUE(:,mv));
end
%% Calculate rms error of paths
lw_rms = cum_rms_lw/length(lw)
sw_rms = cum_rms_sw/length(sw)
aw_rms = cum_rms_aw/length(aw)
xodo_rms = cum_rms_xodo/Nstep
%% Plot
plot(lm(1,:),lm(2,:),'^','LineWidth',3);
plot(XTRUE(1,:),XTRUE(2,:),'-r','LineWidth',3);
plot(XODO(1,:),XODO(2,:),'-k','LineWidth',3);
plot(lw(1,:),lw(2,:),'-b','LineWidth',2);
plot(sw(1,:),sw(2,:),'-g','LineWidth',2);
plot(aw(1,:),aw(2,:),'-m','LineWidth',2);
legend('Start',...
       'End',...
       'Landmark',...
       'Real Pose',...
       'Predict Pose',...
       'Largest Weight',...
       'Smallest Weight',...
       'Average Weight')
   for i = 1:10
        text(lm(1,i)+1,lm(2,i),sprintf('lm %d',i));
   end
end

%% UTIL FUNCTIONS
%% countNotNAN
function cnt = countNotNAN(Z)
%OBJECT: count numbers of landmark detected
Z = isnan(Z);
cnt = length(find(Z ~=1))/3;
end

%% initialise_particles
function p = initialise_particles(np)
for i=1:np
    p(i).w = 1/np;
    p(i).x = [0;0;0];
    p(i).z = [0;0];
end
end
%% likelihood_calculate
function wt = likelihood_calculate(zt, xt, SIGmeas)
% evaluate likelihood of measurement zt given prior with mean xt and std SIGmeas
wt = 1/(2*pi*SIGmeas^2)^(1/2)*exp(-1/2*(xt-zt)^2/SIGmeas^2);
end
%% lsa
function [lw,sw,aw] = lsa(CHIbar_t)
% @OBJECT: find largest,smallest,average value in CHIbar_t
weights = CHIbar_t(:,4);
x_bar = CHIbar_t(:,1);
y_bar = CHIbar_t(:,2);

[l_val,l_idx] = max(weights);
lw = [x_bar(l_idx);y_bar(l_idx)];

[s_val,s_idx] = min(weights);
sw = [x_bar(s_idx);y_bar(s_idx)];

a_val = median(weights);
a_idx = find(abs(weights-a_val)<0.005,1,'first');
aw = [x_bar(a_idx);y_bar(a_idx)];
end
%% choose_landmark
function idx = choose_landmark(CHI,mv,Z,NLANDMARKS,lm)
SIGmeas_range = 0.2;
SIGmeas_bearing = 2/180*pi;
likelihood = 0;
idx = Z(3,1,mv); % landmark index in landmark table
for lm_idx = 1:NLANDMARKS
    tmp_x = CHI(1,1);
    tmp_y = CHI(1,2);
    tmp_phi = CHI(1,3);
    dist = sqrt((tmp_x-lm(1,lm_idx)^2+(tmp_y-lm(2,lm_idx))^2)); % distance to lm_idx landmark
    w_dist = likelihood_calculate(Z(1,lm_idx,mv),dist,SIGmeas_range);
    tmp_lm_x = lm(1,Z(3,lm_idx,mv));
    tmp_lm_y = lm(2,Z(3,lm_idx,mv));
    tmp_dx = tmp_lm_x - tmp_x;
    tmp_dy = tmp_lm_y - tmp_y;
    angle = atan2(tmp_dy,tmp_dx) - tmp_phi;
    w_angle = likelihood_calculate(Z(2,lm_idx,mv),angle,SIGmeas_bearing);
    wt = w_dist+w_angle;
    if(wt>likelihood)
        likelihood = wt;
        idx = Z(3,lm_idx,mv);
    end
end
end
