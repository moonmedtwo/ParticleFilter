function Pacticle_Filter_2d
format compact
%% PREDEFINES
% Data
load ('data.mat');
% Control signal
ut = [5,0];
% DEFINE
dT = 0.025;
wB = 4;
% Max moves
Nstep = 625;
% Noise
% Standard deviation of measurement range, bearing
SIGmeas_range = 1;                        
SIGmeas_bearing = 2/180*pi;
% Standard deviation of process v,tht
SIGproc_v = 0.5;                   
SIGproc_tht = 3/180*pi;           
SIGproc = [SIGproc_v,SIGproc_tht]; 
% Number of particles (the more, the better the particle
NPARTICLES = 100;
%% INITIALIZATION
 CHI = [zeros(NPARTICLES,3),ones(NPARTICLES,1)/NPARTICLES];
%% MAIN
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
% Main loop
for mv = 1:Nstep
    % Compute true data: output XTRUE DONE
    % laser data available at step == 8,16,.. 
    if(mod(mv,8)==0)
        % Count numbers of detected landmark
        NLANDMARKS = countNotNAN(Z(:,:,mv));
        if(NLANDMARKS > 0)
            % Choose the first found landmark
            lm_idx = Z(3,1,mv);
            N = 0; %normalization factor for cumsum function
            CHIbar_t = zeros(NPARTICLES,4);
            CHI_t = zeros(NPARTICLES,4);
            % FIRST SAMPLE and WEIGHT CALCULATE
            for i=1:NPARTICLES
                % No given ut or waypoint, cannot sample
                % Instead using XODO 
%                 xtm1 = [XODO(1,mv),XODO(2,mv),XODO(3,mv)];
%                 xt = sample_motion_model_2d(ut,xtm1,SIGproc,dT,wB);
                xt = [XODO(1,mv)+rand(1),XODO(2,mv)+rand(1),XODO(3,mv)+3/pi*180*rand(1)];
                cum_rms_nf = cum_rms_nf +...
                    sqrt((xt(1)-XTRUE(1,mv))^2+(xt(2)-XTRUE(2,mv))^2);
                
                [dist,b_ang] = xyphi2zPredict(xt,lm(:,lm_idx));
                % Compute weight
                % w1: weight of dist, w2: weigth of bearing angle
                w1 = measurement_model(Z(1,1,mv),dist,SIGmeas_range);
                w2 = measurement_model(Z(2,1,mv),b_ang,SIGmeas_bearing);
                wt = 0.5*w1+0.5*w2;
                CHIbar_t(i,:) = [xt,wt];
                N = N+wt;
            end
            % RESAMPLE 
            cumwt = cumsum(CHIbar_t(:,4))/N;
            for i=1:NPARTICLES
                swt = rand;
                index = find( cumwt>= swt,1,'first');
                xt = CHIbar_t(index,1:3);
                CHI_t(i,:) = [xt,1/NPARTICLES];
            end
            CHI = CHI_t;
            
            % Find largest and smallest weight particle
            weights = CHIbar_t(:,4);
            x_bar = CHIbar_t(:,1);
            y_bar = CHIbar_t(:,2);
            [l_val,l_idx] = max(weights);
            lw = [lw [x_bar(l_idx);y_bar(l_idx)]];
            cum_rms_lw = cum_rms_lw + ...
            sqrt((x_bar(l_idx)-XTRUE(1,mv))^2+(y_bar(l_idx)-XTRUE(2,mv))^2);
            
            [s_val,s_idx] = min(weights);
            sw = [sw [x_bar(s_idx);y_bar(s_idx)]];
            cum_rms_sw = cum_rms_sw + ...
                sqrt((x_bar(s_idx)-XTRUE(1,mv))^2+(y_bar(s_idx)-XTRUE(2,mv))^2);
        
            a_val = median(weights);
            a_idx = find(abs(weights-a_val)<0.005,1,'first');
            aw = [aw [x_bar(a_idx);y_bar(a_idx)]];
            cum_rms_aw = cum_rms_aw + ...
                sqrt((x_bar(a_idx)-XTRUE(1,mv))^2+(y_bar(a_idx)-XTRUE(2,mv))^2);
            
        end
    end
end
%% Calculate rms error of paths
lw_rms = cum_rms_lw/length(lw)
sw_rms = cum_rms_sw/length(sw)
aw_rms = cum_rms_aw/length(aw)
nf_rms = cum_rms_nf/length(aw)
%% Plot
plot(lm(1,:),lm(2,:),'^','LineWidth',3);
hold on;
plot(XTRUE(1,:),XTRUE(2,:),'-r','LineWidth',3);
plot(XODO(1,:),XODO(2,:),'-k','LineWidth',3);
plot(lw(1,:),lw(2,:),'-b','LineWidth',2);
plot(sw(1,:),sw(2,:),'-g','LineWidth',2);
plot(aw(1,:),aw(2,:),'-m','LineWidth',2);
legend('Landmark',...
       'Real Pose',...
       'Predict Pose',...
       'Largest Weight',...
       'Smallest Weight',...
       'Average Weight')
end

%% UTIL FUNCTIONS
function cnt = countNotNAN(Z)
%OBJECT: count numbers of landmark detected
Z = isnan(Z);
cnt = length(find(Z ~=1))/3;
end

function p = initialise_particles(np)
for i=1:np
    p(i).w = 1/np;
    p(i).x = [0;0;0];
    p(i).z = [0;0];
end
end

function wt = measurement_model(zt, xt, SIGmeas)
% evaluate likelihood of measurement zt given prior with mean xt and std SIGmeas
wt = 1/(2*pi*SIGmeas^2)^(1/2)*exp(-1/2*(xt-zt)^2/SIGmeas^2);
end
