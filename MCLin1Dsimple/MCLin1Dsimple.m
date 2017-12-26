% Particle filter example in 1D. "MCLin1Dsimple" implements Monte-Carlo
% Localization (MCL) in one Dimension.
%
% Simulates particle filtering in 1D using MCL algorithm (table 8.2) from
% "Probabilistic Robotics",  (Intelligent Robotics and Autonomous Agents
% series) Hardcover ? August 19, 2005 by Sebastian Thrun, Wolfram Burgard,
% Dieter Fox.
%
% This simple implementation is designed for teaching.  It has pauses, and
% you must press a key to step through the simulation.  Displays the
% particles before moving in blue (with a histogram), the actual position
% after moving with a red line, the measurement with a magenta line, the
% particles propagated forward by the motion model and weighted according
% to the measurement in black, and the resampled particles in red (with a
% histogram).
%
% Aaron T. Becker, University of Houston
% April 20, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MCLin1Dsimple
format compact
%%% PLOT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;
set(f1,'name', 'Monte Carlo Localization')
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SIGmeas = 1;  %standard deviation of measurement noise
SIGproc = 0.5; %standard deviation of process noise
SIGinit = 4; %standard deviation of initial position

M = 100;  %number of particles (the more, the better the particle
%Probability Mass Function (PMF) matches the true Probability Distribution
% Function (PDF).

CHI = [SIGinit*randn(M,1), ones(M,1)/M];  %array of particles and associated weights
xACTt = SIGinit*randn(1); %true position of the robot (unknown to particles)
ut = 5; %control input (constant motion in x-direction)
max_moves = 10;
for mv = 1:max_moves
    % move ACTUAL robot
    xACTt = sample_motion_model(ut, xACTt, SIGproc);
    %take measurement
    zt = take_measurement(xACTt,SIGmeas);
    % apply particle filter
    CHI_prev = CHI; %saving the old partcles (optional)
    [CHI,CHIbar] = AlgorithmMCL( CHI, ut, zt, SIGproc, SIGmeas);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Plot the results (you can skip these) %%%%%%%%%%%%%%%%%%%%%%%%%
    clf
    plot( CHI_prev(:,1),0.5*ones(size(CHI(:,2))),'b.');
    hold on
    plot( [xACTt,xACTt], [0,1],'-r' );
    plot( [zt,zt], [0,1],'-m' );
    plot( CHIbar(:,1),CHIbar(:,2),'k.');
    plot( CHI(:,1),0.55*ones(size(CHI(:,2))),'r.');
    legend('before move','actual position','measurement','after move,weighted','resampled','location','best')
    
    [N,X] = hist(CHI(:,1));
    hCHIhist = bar(X,N/sum(N),'r');
    uistack(hCHIhist,'bottom')
    
    [N,X] = hist(CHI_prev(:,1));
    hCHIm1hist = bar(X,N/sum(N),'b');
    uistack(hCHIm1hist,'bottom')
    %xlabel('position')
    ylabel('cumulative mass')
    axis tight
    a = axis;
    axis([a(1),a(2),0,1]);
    a = axis;
    
    xlabel('position')
    axis(a);
    title('<Paused>  Press a key to continue')
    pause
    %%%%%%%%%%%%%%%%%%%%%%%%% END Plot the results %%%%%%%%%%%%%%%%%%%%%%%%%
end
end

function [CHIt,CHIbart] =  AlgorithmMCL( CHItm1, ut,zt,SIGproc,SIGmeas)
    %Implements Monte Carlo Localization (MCL) algorithm, Table 8.2, page 252,
    %"Probabilistic Robotics"
    M = size(CHItm1,1);
    CHIbart = zeros(size(CHItm1));
    CHIt = zeros(size(CHItm1));
    N = 0; %normalization factor
    for m = 1:M
        xtm1 = CHItm1(m,1);
        xt = sample_motion_model(ut, xtm1,SIGproc);
        wt = measurement_model(zt, xt, SIGmeas);
        CHIbart(m,:) = [xt,wt];
        N = N+wt;
    end

    % normalize
    cumwt = cumsum(CHIbart(:,2))/N;
    for m = 1:M %Resampling step
        %draw m-th sample with probability proportional to wt with two methods:
        swt = rand; %1.) randomly sample
        %swt = cumwt(end)/M*(m-1/2); %2.) swt steps through CMF, and selects
        %the associated particle.  This is a low-variance sampling method.
        index = find( cumwt>= swt,1,'first');
        xt = CHIbart(index,1);
        CHIt(m,:) = [xt,1/M]; %add particle to CHI
    end
end

function xt = sample_motion_model(ut, xtm1, SIGproc)
%sample with 1D noise
xt = xtm1 + ut+ SIGproc*randn(1);
end

function wt = measurement_model(zt, xt, SIGmeas)
% evaluate likelihood of measurement zt given prior with mean xt and std SIGmeas
wt = 1/(2*pi*SIGmeas^2)^(1/2)*exp(-1/2*(xt-zt)^2/SIGmeas^2);
end

function zt = take_measurement(xACTt,SIGmeas)
% actual state pertubed by Gaussian noise
zt = xACTt+randn(1)*SIGmeas;
end
