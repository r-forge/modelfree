function [pfit,etafit,H] = locglmfit(xfit,r,m,x,h,link,guessing,lapsing,...
    K,p,ker,maxiter,tol)
%
% Local polynomial estimator for the psychometric function (PF) and eta
% function (PF transformed by link) for binomial data; also returns the hat
% matrix.  
%
% Actual calculations are done in LOCGLMFIT_PRIVATE or
% LOCGLMFIT_SPARSE_PRIVATE depending on the size of the data set. Here the
% data are split into several parts to speed up the calculations.  
%
%
%INPUT
%
% xfit - points in which to calculate the estimate
% r - number of successes in points x
% m - number of trials in points x 
% x - stimulus levels 
% h - bandwidth(s)
%
% OPTIONAL INPUT
%
% link - name of the link function to be used; default is 'logit'
% guessing - guessing rate; default is 0
% lapsing - lapsing rate; default is 0
% K - power parameter for Weibull and reverse Weibull link; default is
% 2
% p - degree of the polynomial; default p = 1
% ker - kernel function for weights; default 'normpdf'
% maxiter - maximum number of iterations in Fisher scoring; default is 200
% tol - tolerance level at which to stop Fisher scoring; default is 1e-6
%
% OUTPUT
%
% pfit - value of the local polynomial estimate in points xfit
% etafit - estimate of eta (link of pfit)
% H - hat matrix (OPTIONAL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PROGRAM

%%%% CHECK INPUT PARAMETERS
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

%%%% CHECK ROBUSTNESS OF INPUT PARAMETERS
[ rxfit, cxfit ] = size( xfit );
if cxfit ~= 1
    error('vector of x points should be column vector');
end
clear data;
data(1).content = x;
data(2).content = r;
data(3).content = m;
checkinput( 'psychometricdata2', data );
clear data
[ rh, ch ] = size( h );
if ch ~= 1
    error('vector of bandwidths should be column vector');
end
if ~( length( h ) == 1 || rh == rxfit )
    error( 'bandwidth h must be either a scalar or a vector with the same number of elements as x' );
end
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

%%%%%
%%%%% INITIAL VALUES
split = 20;

Lxfit = length(xfit);
Lx = length(x);

if (Lx > 15),
    % big data
    fun_estim = 'locglmfit_sparse_private';
else
    % small data
    fun_estim = 'locglmfit_private';
end

%%%% SPLIT AND EVALUATION
if length( h ) == 1
    % with Hat matrix
    if nargout == 3,

        if (Lxfit<=split),

            % small x
            [pfit,etafit,H] = feval(fun_estim,xfit,r,m,h,x,link,guessing,lapsing,...
                K,p,ker,maxiter,tol);
        else
            % large x
            % number of parts into which the fitting is divided
            fLx = floor(Lxfit/split);

            % initialise output
            pfit = [];
            etafit = [];
            H = [];   

  
                for i = 0:(fLx-1),
                    % part of the fit
                    [pfit1,etafit1,H1] = feval(fun_estim,xfit(i*split + (1:split)),r,m,...
                        h,x,link,guessing,lapsing,K,p,ker,maxiter,tol);
                    % put the fits together
                    pfit = [pfit;pfit1];
                    etafit = [etafit;etafit1];
                    H = [H;H1];   
                end
                %final part of the fit
                if ((split*fLx)<Lxfit),
                    [pfit1,etafit1,H1] = feval(fun_estim,xfit((1+split*fLx):Lxfit),r,m,...
                        h,x,link,guessing,lapsing,K,p,ker,maxiter,tol);

                    % put the fits together
                    pfit = [pfit;pfit1];
                    etafit = [etafit;etafit1];
                    H = [H;H1];   
                end

 
        end

    else % no Hat matrix

        if (Lxfit<=split),

            % small x
            [pfit,etafit] = feval(fun_estim,xfit,r,m,h,x,link,guessing,lapsing,...
                K,p,ker,maxiter,tol);
        else
            % large x
            % number of parts into which the fitting is divided
            fLx = floor(Lxfit/split);

            % initialise output
            pfit = [];
            etafit = [];

            for i = 0:(fLx-1),
                % part of the fit
                [pfit1,etafit1] = feval(fun_estim,xfit(i*split + (1:split)),r,m,h,x,...
                    link,guessing,lapsing,K,p,ker,maxiter,tol);
                % put the fits together
                pfit = [pfit;pfit1];
                etafit = [etafit;etafit1];
            end
            %final part of the fit
            if ((split*fLx)<Lxfit),
                [pfit1,etafit1] = feval(fun_estim,xfit((1+split*fLx):Lxfit),r,m,h,x,...
                    link,guessing,lapsing,K,p,ker,maxiter,tol);

                % put the fits together
                pfit = [pfit;pfit1];
                etafit = [etafit;etafit1];
            end

        end
    
    end % if nargout == 3
else % if length( h ) == 1
    % with Hat matrix
    if nargout == 3,

        if (Lxfit<=split),

            % small x
            [pfit,etafit,H] = feval(fun_estim,xfit,r,m,h,x,link,guessing,lapsing,...
                K,p,ker,maxiter,tol);
        else
            % large x
            % number of parts into which the fitting is divided
            fLx = floor(Lxfit/split);

            % initialise output
            pfit = [];
            etafit = [];
            H = [];   

                for i = 0:(fLx-1),
                 % part of the fit
                    [pfit1,etafit1,H1] = feval(fun_estim,xfit(i*split + (1:split)),r,m,...
                         h(i*split + (1:split)),x,link,guessing,lapsing,K,p,ker,maxiter,tol);
                    % put the fits together
                    pfit = [pfit;pfit1];
                    etafit = [etafit;etafit1];
                    H = [H;H1];   
                end
                %final part of the fit
                if ((split*fLx)<Lxfit),
                    [pfit1,etafit1,H1] = feval(fun_estim,xfit((1+split*fLx):Lxfit),r,m,...
                        h((1+split*fLx):Lxfit),x,link,guessing,lapsing,K,p,ker,maxiter,tol);
                    % put the fits together
                    pfit = [pfit;pfit1];
                    etafit = [etafit;etafit1];
                    H = [H;H1];   
                end
        end
    else % no Hat matrix
        if (Lxfit<=split),

            % small x
            [pfit,etafit] = feval(fun_estim,xfit,r,m,h,x,link,guessing,lapsing,...
                K,p,ker,maxiter,tol);
        else
            % large x

            % number of parts into which the fitting is divided
            fLx = floor(Lxfit/split);

            % initialise output
            pfit = [];
            etafit = [];

            for i = 0:(fLx-1),
                % part of the fit
                [pfit1,etafit1] = feval(fun_estim,xfit(i*split + (1:split)),r,m,h(i*split + (1:split)),x,...
                    link,guessing,lapsing,K,p,ker,maxiter,tol);
                % put the fits together
                pfit = [pfit;pfit1];
                etafit = [etafit;etafit1];
            end
            %final part of the fit
            if ((split*fLx)<Lxfit),
                [pfit1,etafit1] = feval(fun_estim,xfit((1+split*fLx):Lxfit),r,m,h((1+split*fLx):Lxfit),x,...
                    link,guessing,lapsing,K,p,ker,maxiter,tol);

                % put the fits together
                pfit = [pfit;pfit1];
                etafit = [etafit;etafit1];
            end

        end
    
    end % if nargout == 3
end % if length( h ) == 1