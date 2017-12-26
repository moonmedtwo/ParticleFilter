function xtrue = sample_motion_model_2d(ut,xtm1,SIGproc,dT,wB)
% @OBJECT calculate robot pose each sampling interval
% @input
%   ut control signal
%   xtm1 last robot pose
%   SIGproc process noise
%   dT sampling interval
%   wB wheelBase
v = ut(1) + SIGproc(1)*randn(1);
tht = ut(2) + SIGproc(2)*randn(1);

x1 = xtm1(1);
y1 = xtm1(2);
phi1 = xtm1(3);

dX = [v*dT*cos(tht +phi1),v*dT*sin(tht+phi1),v*dT*sin(tht)/wB];

xtrue = xtm1+dX;
end