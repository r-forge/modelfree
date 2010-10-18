function [sd,th0] = bootstrap_sd_th(TH,r,m,x,N,h0,X,link,guessing,lapsing,...
                        K,p,ker,maxiter,tol)
%
% Finds a bootstrap estimate of the standard deviation of the estimated
% threshold for the local polynomial estimate of the psychometric function
% with guessing and lapsing rates specified in lims.
%
% INPUT
%
% TH - required threshold level
% r - number of successes in points x
% m - number of trials in points x 
% x - stimulus levels 
% N - number of bootstrap replications; N should be at least 200 for 
%   reliable results  
% h0 - bandwidth
%
% OPTIONAL INPUT
%
% X - set of value for which to calculate the estimates of the psychometric
% function for the threshold estimation; if not given 1000 equally spaced
% points from mininmum to maximum of x are used
% link - name of the link function to be used; default is 'logit'
% guessing - guessing rate; default is 0
% lapsing - lapsing rate; default is 0
% K - power parameter for Weibull and reverse Weibull link; default is
% 2
% p - degree of the polynomial; default p = 1
% ker - kernel function for weights; default 'normpdf'
% maxiter - maximum number of iterations in Fisher scoring; default is 50
% tol - tolerance level at which to stop Fisher scoring; default is 1e-6
%
% OUTPUT
% 
% sd - bootstrap estimate of the standard deviation of the threshold
% estimator
% th0 - threshold estimate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PROGRAM

%%%% CHECK INPUT PARAMETERS
% First 6 arguments are mandatory
if (nargin<6)
    error('Check input. First 6 arguments are mandatory');
end

%%%%
%%%% DEFAULTS
if (nargin<7)
    X = linspace(min(x),max(x),1000)';
end

if (nargin<8)
    link = 'logit';
    disp('default link function is ''logit''');
end

if (nargin<9)
    guessing = 0;
    disp('default guessing rate is zero');
end

if (nargin<10)
    lapsing = 0;
    disp('default lapsing rate is zero');
end

if (nargin<11)
    K = 2;
    if strcmp(link, 'weibull')
        disp('default exponent for Weibull link function is 2');
    elseif strcmp(link, 'revweibull')
        disp('default exponent for reverse Weibull link function is 2');
    end
end

if (nargin<12)
    p = 1;
    disp('degree of the polynomial to be fitted on the linear scale is 1');
end

if (nargin<13)
    ker = 'normpdf';
    disp('default kernel is ''normpdf''');
end

if (nargin<14)
    maxiter=200;
    disp('default maximum number of iterations is 200');
end

if (nargin<15)
    tol = 1e-6;
    disp('default tolerance is 1e-6');
end

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
if ~isscalar( TH )
    error( 'threshold level should be a scalar' );
elseif TH <= 0 || TH >= 1
    error( 'threshold level should be between 0 and 1' );
end

clear data;
data(1).content = x;
data(2).content = r;
data(3).content = m;
checkinput( 'psychometricdata2', data );
clear data;
if ~isscalar( N )
    error( 'number of bootstrap samples should be a scalar' );
elseif N <= 0 || round( N ) ~= N
    error( 'number of bootstrap samples have to be a positive integer' );
elseif N < 200,
    warning('smoothpsych:smallN', 'number of bootstrap samples should be greater than 200 \n otherwise the results might be unreliable' )
end
[rX, cX] = size( X );
if cX > 1
    error('X (values where to estimate the PF) has to be column vector');
end
if rX < 2
    error('At least 2 values needed for vector X');
end
checkinput( 'bandwidth', h0 );
checkinput( 'linkfunction', link );
if length( guessing ) > 1
    error( 'guessing rate must be scalar' );
end
if length( lapsing ) > 1
    error( 'lapsing rate must be scalar' );
end
checkinput( 'guessingandlapsing', [ guessing 1 - lapsing ] );
checkinput( 'exponentk', K );
checkinput( 'degreepolynomial', p );
checkinput( 'kernel', ker );
checkinput( 'maxiter', maxiter );
checkinput( 'tolerance', tol );

%%%%% INITIAL VALUES
% data in column vectors
tmp1(:,1) = r;
r = tmp1;

tmp1(:,1) = m;
m = tmp1;

tmp1(:,1) = x;
x = tmp1;

n = length(x);

%%%% INITIAL ESTIMATE
% initial estimates with bandiwdth h0
f = locglmfit(x,r,m,x,h0,link,guessing,lapsing,K,p,ker,...
    maxiter,tol);

% dense version for estimation of the threshold
F = locglmfit(X,r,m,x,h0,link,guessing,lapsing,K,p,ker,maxiter,tol);

%%%%%%%%%
%%%% THRESHOLD ESTIMATE
th0 = threshold_slope(F,X,TH);

%%%%%%%%%
%%%% BOOTSTRAP SAMPLING

% re-sampling
M = repmat(m,1,N);
samp = binornd(M,repmat(f,1,N));

% exclude "degenerate samples" if min(M)>1
if (min(M)>1),
    for i = 1:N,
        ok(i) = (length(unique(samp(:,i)))>3);
    end

    while (min(ok)==0),
        Lok = sum(ok==0);
        samp(:,ok==0) = binornd(repmat(m,1,Lok),repmat(f,1,Lok));
        findok = find(ok==0);
        for i = findok,
            ok(i) = (length(unique(samp(:,i)))>3);
        end
    end
end

%%% INITIATE VARABLE IN WHICH DATA ARE STORED
th_boot = zeros(1,N);

%%%%%%%%%%%%%%%%%%%%%
%%%%%%% BOOTSTRAP ESTIMATES OF THE THRESHOLD
for i = 1:N,

    ftmp = locglmfit(X,samp(:,i),m,x,h0,link,guessing,lapsing,...
                K,p,ker,maxiter,tol);
    th_boot(i) = threshold_slope(ftmp,X,TH);

end

sd = sqrt(var(th_boot));