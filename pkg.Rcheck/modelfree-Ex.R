pkgname <- "modelfree"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('modelfree')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("bandwidth_bootstrap")
### * bandwidth_bootstrap

flush(stderr()); flush(stdout())

### Name: bandwidth_bootstrap
### Title: Bootstrap estimate of bandwidth
### Aliases: bandwidth_bootstrap
### Keywords: nonparametric models regression nonlinear

### ** Examples

data("01_Miranda")
h<- bandwidth_bootstrap( example01$r, example01$m, example01$x, c( 0.1, 10 ), 10 )



cleanEx()
nameEx("bandwidth_cross_validation")
### * bandwidth_cross_validation

flush(stderr()); flush(stdout())

### Name: bandwidth_cross_validation
### Title: Cross-validation estimate of bandwidth
### Aliases: bandwidth_cross_validation
### Keywords: nonparametric models regression nonlinear

### ** Examples

data("01_Miranda")
h<- bandwidth_cross_validation( example01$r, example01$m, example01$x, c( 0.1, 10 ) )



cleanEx()
nameEx("bandwidth_plugin")
### * bandwidth_plugin

flush(stderr()); flush(stdout())

### Name: bandwidth_plugin
### Title: Plug in estimation of Bandwidth
### Aliases: bandwidth_plugin
### Keywords: nonparametric models regression nonlinear

### ** Examples

data("01_Miranda")
h<-bandwidth_plugin( example01$r, example01$m, example01$x )



cleanEx()
nameEx("binom_g")
### * binom_g

flush(stderr()); flush(stdout())

### Name: binom_g
### Title: Psychometric function with guessing rate
### Aliases: binom_g
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
value <-binom_g( example01$r, example01$m, example01$x, "logit", 1, 2, 0.01 )



cleanEx()
nameEx("binom_gl")
### * binom_gl

flush(stderr()); flush(stdout())

### Name: binom_gl
### Title: Psychometric function with guessing and lapsing rates
### Aliases: binom_gl
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" );
value <-binom_gl( example01$r, example01$m, example01$x, "logit", 1, 2, c( 0.01, 0.01 ) );



cleanEx()
nameEx("binom_l")
### * binom_l

flush(stderr()); flush(stdout())

### Name: binom_l
### Title: Psychometric function with lapsing rate
### Aliases: binom_l
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
value <-binom_l( example01$r, example01$m, example01$x, "logit", 1, 2, 0.01 )



cleanEx()
nameEx("binom_lims")
### * binom_lims

flush(stderr()); flush(stdout())

### Name: binom_lims
### Title: Psychometric function with guessing and lapsing rates
### Aliases: binom_lims
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
value <-binom_lims( example01$r, example01$m, example01$x )



cleanEx()
nameEx("binom_revweib")
### * binom_revweib

flush(stderr()); flush(stdout())

### Name: binom_revweib
### Title: Psychometric function fitting for reverse Weibull link function
### Aliases: binom_revweib
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
value <- binom_revweib( example01$r, example01$m, example01$x )



cleanEx()
nameEx("binom_weib")
### * binom_weib

flush(stderr()); flush(stdout())

### Name: binom_weib
### Title: Psychometric function fitting for Weibull link function
### Aliases: binom_weib
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
value <- binom_revweib( example01$r, example01$m, example01$x )



cleanEx()
nameEx("binomfit_lims")
### * binomfit_lims

flush(stderr()); flush(stdout())

### Name: binomfit_lims
### Title: Generalized linear model fit with guessing and lapsing rates
### Aliases: binomfit_lims
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
value <- binomfit_lims( example01$r, example01$m, example01$x )



cleanEx()
nameEx("bootstrap_ci_sl")
### * bootstrap_ci_sl

flush(stderr()); flush(stdout())

### Name: bootstrap_ci_sl
### Title: Bootstrap estimate of confidence interval for slope estimation
### Aliases: bootstrap_ci_sl
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
bwd <- 0.2959
value <- bootstrap_ci_sl( 0.5, example01$r, example01$m, example01$x, 10, bwd )



cleanEx()
nameEx("bootstrap_ci_th")
### * bootstrap_ci_th

flush(stderr()); flush(stdout())

### Name: bootstrap_ci_th
### Title: Bootstrap estimate of confidence interval for threshold
###   estimation
### Aliases: bootstrap_ci_th
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
bwd <- 0.2959;
value <- bootstrap_ci_th( 0.5, example01$r, example01$m, example01$x, 10, bwd );



cleanEx()
nameEx("bootstrap_sd_sl")
### * bootstrap_sd_sl

flush(stderr()); flush(stdout())

### Name: bootstrap_sd_sl
### Title: Bootstrap estimate the standard deviation of slope estimation
### Aliases: bootstrap_sd_sl
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
bwd <- 0.2959
value <- bootstrap_sd_sl( 0.5, example01$r, example01$m, example01$x, 10, bwd )



cleanEx()
nameEx("bootstrap_sd_th")
### * bootstrap_sd_th

flush(stderr()); flush(stdout())

### Name: bootstrap_sd_th
### Title: Bootstrap estimate the standard deviation of threshold
###   estimation
### Aliases: bootstrap_sd_th
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
bwd <- 0.2959
value <- bootstrap_sd_th( 0.5, example01$r, example01$m, example01$x, 10, bwd )



cleanEx()
nameEx("comploglog_link")
### * comploglog_link

flush(stderr()); flush(stdout())

### Name: comploglog_link
### Title: Complementary log-log link function with guessing and lapsing
###   rates
### Aliases: comploglog_link
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-comploglog_link( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("comploglog_link_private")
### * comploglog_link_private

flush(stderr()); flush(stdout())

### Name: comploglog_link_private
### Title: Complementary log-log link function with guessing and lapsing
###   rates
### Aliases: comploglog_link_private
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-comploglog_link_private( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("deviance2")
### * deviance2

flush(stderr()); flush(stdout())

### Name: deviance2
### Title: Deviance between data and fitted function
### Aliases: deviance2
### Keywords: nonparametric

### ** Examples

data( "01_Miranda" )
h = 0.2959
fit <- locglmfit( example01$x, example01$r, example01$m, example01$x, h )
Dev <- deviance2( example01$r, example01$m, fit$fitval )



cleanEx()
nameEx("locglmfit")
### * locglmfit

flush(stderr()); flush(stdout())

### Name: locglmfit
### Title: Local generalized linear fitting
### Aliases: locglmfit
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
xnew = 1.2 * (0:99)/99+0.1
h <- 0.2959
fit <- locglmfit( xnew, example01$r, example01$m, example01$x, h )



cleanEx()
nameEx("locglmfit_private")
### * locglmfit_private

flush(stderr()); flush(stdout())

### Name: locglmfit_private
### Title: Local generalized linear fitting with usual (non-sparse)
###   matrices
### Aliases: locglmfit_private
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
xnew = 1.2 * (0:99)/99+0.1
h <- 0.2959
fit <- locglmfit_private( xnew,  example01$r,  example01$m, example01$x, h, FALSE, "logit_link", 0, 0, 2, 1, "dnorm", 50, 1e-6)



cleanEx()
nameEx("locglmfit_sparse_private")
### * locglmfit_sparse_private

flush(stderr()); flush(stdout())

### Name: locglmfit_sparse_private
### Title: Local generalized linear fitting with sparse matrices
### Aliases: locglmfit_sparse_private
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
xnew = 1.2 * (0:99)/99+0.1
h <- 0.2959
fit <- locglmfit_sparse_private( xnew,  example01$r,  example01$m, example01$x, h, FALSE, "logit_link", 0, 0, 2, 1, "dnorm", 50, 1e-6)



cleanEx()
nameEx("logit_link")
### * logit_link

flush(stderr()); flush(stdout())

### Name: logit_link
### Title: Logit link function with guessing and lapsing rates
### Aliases: logit_link
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-logit_link( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("logit_link_private")
### * logit_link_private

flush(stderr()); flush(stdout())

### Name: logit_link_private
### Title: Logit link function with guessing and lapsing rates
### Aliases: logit_link_private
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-logit_link_private( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("loglog_link")
### * loglog_link

flush(stderr()); flush(stdout())

### Name: loglog_link
### Title: Log-log link function with guessing and lapsing rates
### Aliases: loglog_link
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-loglog_link( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("loglog_link_private")
### * loglog_link_private

flush(stderr()); flush(stdout())

### Name: loglog_link_private
### Title: Log-log link function with guessing and lapsing rates
### Aliases: loglog_link_private
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-loglog_link_private( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("probit_link")
### * probit_link

flush(stderr()); flush(stdout())

### Name: probit_link
### Title: Probit link function with guessing and lapsing rates
### Aliases: probit_link
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-probit_link( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("probit_link_private")
### * probit_link_private

flush(stderr()); flush(stdout())

### Name: probit_link_private
### Title: Probit link function with guessing and lapsing rates
### Aliases: probit_link_private
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-probit_link_private( 0.1, 0.1 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("revweibull_link")
### * revweibull_link

flush(stderr()); flush(stdout())

### Name: revweibull_link
### Title: Reverse Weibull link function with guessing and lapsing rates
### Aliases: revweibull_link
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-revweibull_link( 20 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("revweibull_link_private")
### * revweibull_link_private

flush(stderr()); flush(stdout())

### Name: revweibull_link_private
### Title: Reverse Weibull link function with guessing and lapsing rates
### Aliases: revweibull_link_private
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-revweibull_link_private( 20, 0, 0 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("threshold_slope")
### * threshold_slope

flush(stderr()); flush(stdout())

### Name: threshold_slope
### Title: Threshold and slope of estimated psychometric function
### Aliases: threshold_slope
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
xnew = 1.2 * (0:999)/999+0.1
h = 0.2959
fit <- locglmfit( xnew, example01$r, example01$m, example01$x, h )
value <- threshold_slope( fit$pfit , xnew )



cleanEx()
nameEx("weibull_link")
### * weibull_link

flush(stderr()); flush(stdout())

### Name: weibull_link
### Title: Weibull link function with guessing and lapsing rates
### Aliases: weibull_link
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-weibull_link( 20 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



cleanEx()
nameEx("weibull_link_private")
### * weibull_link_private

flush(stderr()); flush(stdout())

### Name: weibull_link_private
### Title: Weibull link function with guessing and lapsing rates
### Aliases: weibull_link_private
### Keywords: nonparametric models regression nonlinear

### ** Examples

data( "01_Miranda" )
x <- example01$x
r <- example01$r
m <- example01$m
glmdata <- data.frame( cbind( r/m ,m , x ) )
names( glmdata ) <- c( "resp", "m", "x" )
glmformula <- c( "resp ~ x" )
userlink<-weibull_link_private( 20, 0, 0 )
fit <- glm( glmformula, data = glmdata, weights = m, family = binomial( userlink ) )



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
