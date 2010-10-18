function h = bandwidth_bootstrap(r,m,x,H,N,h0,link,guessing,lapsing,...
                        K,p,ker,maxiter,tol,method)
%
% Finds bootstrap estimate of the optimal bandwidth h for a local
% polynomial estimate of the psychometric function with specified guessing
% and lapsing rates.
%
% INPUT
%
% r - number of successes in points x
% m - number of trials in points x 
% x - stimulus levels 
% H - minimum and maximum values of bandwidth to be considered
% N - number of bootstrap replications
%
% OPTIONAL INPUT
%
% h0 - pilot bandwidth; if not specified, then the scaled plug-in bandwidth
% is used
% link - name of the link function to be used; default is 'logit'
% guessing - guessing rate; default is 0
% lapsing - lapsing rate; default is 0
% K - power parameter for Weibull and reverse Weibull link; default is
% 2
% p - degree of the polynomial; default p=1
% ker - kernel function for weights; default 'normpdf'
% maxiter - maximum number of iterations in Fisher scoring; default is 50
% tol - tolerance level at which to stop Fisher scoring; default is 1e-6
% method - loss function to be used in bootstrap: choose from:
% 'ISEeta', 'ISE', 'deviance'; by default all possible values are
% calculated
%
% OUTPUT
% 
% h - bootstrap bandwidth for the chosen "method"; if no method was
% specified, then it is three row vector with entries corresponding to the
% estimated bandwidths on p-scale, on eta-scale and for mean deviance 

% First 5 arguments are mandatory
if (nargin<5)
    error('Check input. First 5 arguments are mandatory');
end

%%%%
%%%% DEFAULTS
if (nargin<6)
    h0 = [];
end

if (nargin<7)
    link = 'logit';
    disp('default link function is ''logit''');
end

if (nargin<8)
    guessing = 0;
    disp('default guessing rate is zero');
end

if (nargin<9)
    lapsing = 0;
    disp('default lapsing rate is zero');
end

if (nargin<10)
    K = 2;
    if strcmp(link, 'weibull')
        disp('default exponent for Weibull link function is 2');
    elseif strcmp(link, 'revweibull')
        disp('default exponent for reverse Weibull link function is 2');
    end
end

if (nargin<11)
    p = 1;
    disp('degree of the polynomial to be fitted on the linear scale is 1');
end

if (nargin<12)
    ker = 'normpdf';
    disp('default kernel is ''normpdf''');
end

if (nargin<13)
    maxiter=200;
    disp('default maximum number of iterations is 200');
end

if (nargin<14)
    tol = 1e-6;
    disp('default tolerance is 1e-6');
end

if (nargin<15)
    disp('bootstrap bandwidth calculated for the 3 methods');
    method = 'all';
end

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
clear data;
data(1).content = x;
data(2).content = r;
data(3).content = m;
checkinput( 'psychometricdata2', data );
clear data
checkinput( 'minmaxbandwidth', H );
checkinput( 'bootstrapreplications', N );
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
checkinput( 'method', method );


%%%%%
%%%%% INITIAL VALUES

% data in column vectors
tmp1(:,1) = r;
r = tmp1;

tmp1(:,1) = m;
m = tmp1;

tmp1(:,1) = x;
x = tmp1;

n = length(x);

%%%% OBTAIN INITIAL BANDWIDTH, IF NOT GIVEN
if (isempty(h0)), h0 = (1.5 * n^.1) * bandwidth_plugin(r,m,x,p,str2func(ker),link); end

%%%% INITIAL ESTIMATE
% pilot over-smoothed estimate with bandiwdth h0
[f,eta] = locglmfit(x,r,m,x,h0,link,guessing,lapsing,K,p,ker,...
    maxiter,tol);

%%%% SAMPLING

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
lH =length(H);

f1 = repmat(f,1,N);
eta1 = repmat(eta,1,N);


%%%% INITIAL VALUES FOR OPTIMISATION
options = optimset('MaxIter', 50); 

%%%%%%%%%%%% BANDWIDTH
% MINIMA FOR MEAN LOSS FUNCTIONS

if strcmp(method,'ISE'),
    % MISE
    h = fminbnd(@get_mise_p,H(1),H(end),options);

elseif strcmp(method,'ISEeta'),
    % MISE ETA
    h = fminbnd(@get_mise,H(1),H(end),options);

elseif strcmp(method,'deviance'),
    % DEVINACE
    h = fminbnd(@get_dev,H(1),H(end),options);

else
    % MISE
    h(1,:) = fminbnd(@get_mise_p,H(1),H(end),options);
    
    % MISE ETA
    h(2,:) = fminbnd(@get_mise,H(1),H(end),options);
    
    % DEVINACE
    h(3,:) = fminbnd(@get_dev,H(1),H(end),options);
end


% -------------------------------------------------------------------------
	function mise = get_mise_p(h)
		% get mise for this value of h
		% this is a nested function, so it shares variables!
		for i = 1:N,
			ftmp = locglmfit(x,samp(:,i),m,x,h,link,...
                            guessing,lapsing,K,p,ker,maxiter,tol);
            fest(:,i) = ftmp;
        end
		mise = mean(ISE(f1, fest)); % return MISE for this h
	end
% -------------------------------------------------------------------------
	function mise = get_mise(h)
		% get mise on eta scale on for this value of h
		% this is a nested function, so it shares variables!
		for i = 1:N,
            [ftmp,etatmp] = locglmfit(x,samp(:,i),m,x,h,link,...
                            guessing,lapsing,K,p,ker,maxiter,tol);
            fest(:,i) = etatmp;
        end
		mise = mean(ISE(eta1, fest)); % return MISE for this h
    end
% -------------------------------------------------------------------------
	function mD = get_dev(h)
		% get mean devinace for this value of h
		% this is a nested function, so it shares variables!
		for i = 1:N,
			ftmp = locglmfit(x,samp(:,i),m,x,h,link,...
                            guessing,lapsing,K,p,ker,maxiter,tol);
            D(i) = deviance(r, m, ftmp);
        end
		mD = mean(D); % return MISE for this h
	end
% -------------------------------------------------------------------------


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% INTERNAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% LOSS FUNCTION

function y = ISE(f1,f2)

y = sum((f1 - f2).^2);
end


function D = deviance(Y,f)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PROGRAM

% retrieve data
k = Y(:,1);
m = Y(:,2);

% adjustment to aviod degenerate values
k(k>=m) = k(k>=m)-.001;
k(k<=0) = .001;

f(f>=1) = 1-.001;
f(f<=0) = .001;   


% deviance
D = 2*sum((k .* log(k./(m.*f)) + (m-k) .* log((m-k)./(m-m.*f))));

end