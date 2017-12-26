%% cum_rms
function cum_rms = cal_cumrms(cum_rms,val,trueval)
% @OBJECT calculate cumulative root mean square error
cum_rms = cum_rms + ...
    sqrt((val(1)-trueval(1))^2 + (val(2)-trueval(2))^2);
end

