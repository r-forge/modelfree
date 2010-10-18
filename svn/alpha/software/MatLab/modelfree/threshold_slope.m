function  [x_th,slope] = threshold_slope(pfit,xfit,thresh)
%
% Finds the approximate value of x (=x_th) for which the value of the
% psychometric function is equal thresh and the approximate value of slope
% in x_th  
% 
% INPUT
% 
% pfit - estimated values of the psychometric function
% xfitr - stimulus levels in which the function was estimated
%
% OPTIONAL INPUT
% 
% thresh - value for which to estimate threshold; default is 0.5
%
% OUTPUT
% 
% x_th - estimated threshold
% slope - estimated values of slope, i.e. discrete derivative of pfit at x_th

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PROGRAM

%%%% CHECK INPUT PARAMETERS
% First 2 arguments are mandatory
if (nargin<2)
    error('Check input. First 2 arguments are mandatory');
end

%%%%
%%%% DEFAULT VALUES
if (nargin<3), thresh = .5; end

%%%%
%%%% INITIAL VALUES
% data in column vectors
tmp1(:,1) = pfit;
pfit = tmp1;

tmp1(:,1) = xfit;
xfit = tmp1;

% threshold
x_th = xfit(abs(pfit-thresh)==min(abs(pfit-thresh)));

if (length(x_th)>1), 
    % if there are many point for the same threshold value, then function
    % is flat in this point and slope=0
    slope = 0;
    x_th = mean(x_th);
else
    % slope
    ind = find(xfit == x_th);
    slope = (pfit(min(ind+1,end))-pfit(max(ind-1,1)))/(xfit(min(ind+1,end))-xfit(max(ind-1,1)));
end