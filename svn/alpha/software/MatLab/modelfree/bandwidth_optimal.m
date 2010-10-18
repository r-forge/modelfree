function h = bandwidth_optimal(ptrue,r,m,x,H,link,guessing,lapsing,...
                            K,p,ker,maxiter,tol,method)
%
% Finds the optimal bandwidth for a local polynomial estimate
% of the psychometric function with specified guessing and lapsing rates. The
% difference between this function and bandwidth_cross_validation is that
% here the true psychometric function is known.
%
% INPUT
% 
% ptrue - the true function. Vector with the value of ptrue at each
% stimulus level x
% r - number of successes in points x
% m - number of trials in points x 
% x - stimulus levels 
% H - minimum and maximum values of bandwidth to be considered
%
% OPTIONAL INPUT
%
% link - name of the link function to be used; default is 'logit'
% guessing - guessing rate; default is 0
% lapsing - lapsing rate; default is 0
% K - power parameter for Weibull and reverse Weibull link; default is
% 2
% p - degree of the polynomial; default p=1
% ker - kernel function for weights; default 'normpdf'
% maxiter - maximum number of iterations in Fisher scoring; default is 50
% tol - tolerance level at which to stop Fisher scoring; default is 1e-6
% method - loss function to be used in cross-validation: choose from:
% 'ISEeta', 'ISE', 'deviance'; by default all possible values are
% calculated
%
% OUTPUT
% 
% h - optimal bandwidth for the chosen "method"; if no method was
% specified, then it is three row vector with entries corresponding to the
% estimated bandwidths on p-scale, on eta-scale and for mean deviance 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PROGRAM
% First 5 arguments are mandatory
if (nargin<5)
    error('Check input. First 5 arguments are mandatory');
end

%%%%
%%%% DEFAULTS
if (nargin<6)
    link = 'logit';
    disp('default link function is ''logit''');
end

if (nargin<7)
    guessing = 0;
    disp('default guessing rate is zero');
end

if (nargin<8)
    lapsing = 0;
    disp('default lapsing rate is zero');
end

if (nargin<9)
    K = 2;
    if strcmp(link, 'weibull')
        disp('default exponent for Weibull link function is 2');
    elseif strcmp(link, 'revweibull')
        disp('default exponent for reverse Weibull link function is 2');
    end
end

if (nargin<10)
    p = 1;
    disp('degree of the polynomial to be fitted on the linear scale is 1');
end

if (nargin<11)
    ker = 'normpdf';
    disp('default kernel is ''normpdf''');
end

if (nargin<12)
    maxiter=200;
    disp('default maximum number of iterations is 200');
end

if (nargin<13)
    tol = 1e-6;
    disp('default tolerance is 1e-6');
end

if (nargin<14)
    disp('optimal bandwidth calculated for the 3 methods');
    method = 'all';
end

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS

clear data;
data(1).content = x;
data(2).content = r;
data(3).content = m;
checkinput( 'psychometricdata2', data );
clear data
% The vector ptrue with values of the true psychometric functions should be a
% column vector with the same length as x. Furthermore, the elements of
% ptrue should be in [0,1]
[ rptrue cptrue ] = size( ptrue );
if cptrue ~= 1 || rptrue ~= length( x )
    error( 'ptrue should be a column vector of the same length as x' );
end
if any( ptrue < 0  | ptrue > 1 )
    error( 'elements of ptrue are probabilities and therefore between 0 and 1' );
end
checkinput( 'minmaxbandwidth', H );
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

% data in column vectors
tmp1(:,1) = r;
r = tmp1;

tmp1(:,1) = m;
m = tmp1;

tmp1(:,1) = x;
x = tmp1;


% initial values for the loops
LH = length(H);
Lx = length(x);

options = optimset('MaxIter', 50); 

%%%%% CROSS-VALIDATION
%%%%% BANDIWDTH

if strcmp(method,'ISE'),
    %%% p-scale
    h = fminbnd(@get_ise_p,H(1),H(end),options);
    
elseif strcmp(method,'ISEeta'),
    %%%%% eta scale
    h = fminbnd(@get_ise,H(1),H(end),options);
    
elseif strcmp(method,'deviance'),
    %%%%% DEVINACE
    h = fminbnd(@get_dev,H(1),H(end),options);
    
else
    %%% p-scale
    h(1,:) = fminbnd(@get_ise_p,H(1),H(end),options);
    %%%%% eta scale
    h(2,:) = fminbnd(@get_ise,H(1),H(end),options);
    %%%%% DEVINACE
    h(3,:) = fminbnd(@get_dev,H(1),H(end),options);    
end


% -------------------------------------------------------------------------
	function ise = get_ise_p(h)
		% get ise for this value of h
		% this is a nested function, so it shares variables!
		fest = locglmfit(x,r,m,x,h,link,guessing,lapsing,...
            K,p,ker,maxiter,tol);
            
		ise = ISE(ptrue, fest); % return MISE for this h
 
	end
% -------------------------------------------------------------------------
	function ise = get_ise(h)
		% get ise on eta scale on for this value of h
		% this is a nested function, so it shares variables!
        [ftmp,etaest] = locglmfit(x,r,m,x,h,link,guessing,...
            lapsing,K,p,ker,maxiter,tol);

        fit_eta = log(ptrue./(1-ptrue));
      
        ise = ISE(fit_eta, etaest); % return ISE for this h
    end
% -------------------------------------------------------------------------
	function D = get_dev(h)
		% get devinace for this value of h
		% this is a nested function, so it shares variables!
    	ftmp = locglmfit(x,r,m,x,h,link,guessing,lapsing,...
            K,p,ker,maxiter,tol);
    
        D = deviance(ptrue.*m, m,ftmp); % return MISE for this h
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
