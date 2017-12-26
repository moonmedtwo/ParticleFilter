function [d,b] = xyphi2zPredict(curXODO,lm)
% OBJECT: Convert x,y,phi to distance and bearing angle seen from vehicle
% OBJECT: to a specific land mark
% d: distance predicted
% b: bearing angle predicted
% input
% curXODO: current predicted x,y,phi
% lm: landmark specified

lmx = lm(1);
lmy = lm(2);

curX = curXODO(1);
curY = curXODO(2);
curPhi = curXODO(3);

dx = lmx - curX;
dy = lmy - curY;

d = sqrt((curX-lmx)^2+(curY-lmy)^2);
b = atan2(dy,dx) - curPhi;
