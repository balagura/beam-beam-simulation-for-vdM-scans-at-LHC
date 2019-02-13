## Comparison of the new and old beam-beam corrections for round beams
##
## Author V. Balagura, balagura@cern.ch (Jan 2019)

## minimal part of old beam-beam correction code from CERN gitlab, 
## ------------------------------------------------------------
##
## Author: Rosen Matev (Rosen.Matev@cern.ch)
##
library(plyr)
library(stringr)
library(reshape2)

splitIdColumn <- function(df, pattern, groups, variable.name='variable', keep.var=F) {
  vars <- str_match(df[, variable.name], pattern)[, groups, drop=F]
  ids <- data.frame(alply(vars, 2, type.convert, as.is=T), stringsAsFactors=F)
  names(ids) <- names(groups)
  newdf <- cbind(df, ids)
  if (!keep.var) newdf[variable.name] <- NULL
  newdf
}

C_BBP <- 1.44e-12/(2*pi)  # alpha*hbar*c/(2*pi) in GeV.um
M_P <- 0.938  # proton mass in GeV

dynbetaInfo <- read.table('old_dynbeta.txt', header=T, quote='', nrows=1)
.tmp <- read.table('old_dynbeta.txt', header=T, quote='', skip=2)
dynbetaTable <- dcast(splitIdColumn(melt(.tmp, id.vars=1), '(.+)\\.([XY])', groups=c(var=2, dir=3)), ...~var)
dynbetaRefTable <- splitIdColumn(melt(.tmp, id.vars=1, value.name='reldb'), 'reldb([xy])_ref\\.([XY])', groups=c(plane=2, dir=3))
rm(.tmp)

quadroot <- function(c, b, a) {
  D <- b^2 - 4*a*c + 0i
  x1 <- (-b-sqrt(D))/(2*a)
  x2 <- (-b+sqrt(D))/(2*a)
  rbind(x1, x2)
}

xiFromRelBetaChange <- function(reldb, tune, method=c('exact', 'approx')) {
  method <- match.arg(method)
  if (method == 'exact') {
    #roots <- polyroot(c((1+df$reldb)^(-2)-1, -4*pi/tan(2*pi*df$Q_ref), 4*pi^2))
    roots <- quadroot((1+reldb)^(-2)-1, -4*pi/tan(2*pi*tune), 4*pi^2)
    if (any(Im(roots) != 0)) stop('Complex roots!')
    roots <- Re(roots)
    if (any(roots[1,] > roots[2,])) stop('Unexpected roots!')
    roots[2,]
  } else if (method == 'approx') {
    -tan(2*pi*tune)/(2*pi)*reldb
  }
}

dynbetaRefTable <- merge(dynbetaRefTable, data.frame(plane=c('x', 'y'), Q_ref=with(dynbetaInfo, c(Qx, Qy))))
dynbetaRefTable <- mutate(dynbetaRefTable,
                          #scale = -sqrt(2)*tan(2*pi*Q_ref)*with(dynbetaInfo, mpc2*eps_ref/n_ref), # [scale] = GeV.um
                          #reldb_inv = scale*reldb, # [reldb_inv] = GeV.um
                          xi = xiFromRelBetaChange(reldb, Q_ref),
                          xi_inv = xi * with(dynbetaInfo, 2*eps_ref/n_ref/(C_BBP/M_P)) # [xi_inv] = 1
)


relBetaChangeFromXi <- function(xi, tune) {
  1/sqrt(1 + 4*pi*xi/tan(2*pi*tune) - 4*pi^2*xi^2) - 1
}

dynbetaXiInvSlice <- function(nsep, plane, dir) {
  d <- dir
  p <- plane
  data <- subset(dynbetaRefTable, plane==p & dir==d)
  approx(data$beamSigmaSep, data$xi_inv, abs(nsep), method='linear', rule=2)$y
}

dynbetaXiInv <- function(nsepx, nsepy, plane) {
  # Return xi^(inv)_{`plane`}(`nsepx`, `nsepy`) in GeV.um
  g_00 <- c(dynbetaXiInvSlice(0, plane, 'X'), dynbetaXiInvSlice(0, plane, 'Y'))
  stopifnot(g_00[1] == g_00[2])
  g_00 <- mean(g_00)
  g_x0 <- dynbetaXiInvSlice(nsepx, plane, 'X')
  g_0y <- dynbetaXiInvSlice(nsepy, plane, 'Y')
  g_x0*g_0y/g_00
}

dynbetaXi <- function(nsepx, nsepy, plane, energy, beta, n_op, sigx_op, sigy_op) {
  # Return the beam-beam parameters at normalized separation (nsepx, nsepy) for a single plane and single beam
  #   nsepx = sepx/sigx; nsepy = ...
  #   energy = beam energy in GeV
  #   beta = beam nominal beta in plane `plane` in m
  #   n_op = opposite beam charge
  #   sigx(y)_op = opposite beam sigma_x(y) in um
  #sig_op <- ifelse(plane=='x', sigx_op, ifelse(plane=='y', sigy_op, NA))
  stopifnot(length(plane)==1)
  sig_op <- if (plane=='x') sigx_op else if (plane=='y') sigy_op else stop('Plane must be either x or y')
  scale <- C_BBP*n_op*(beta*1e6)/(energy*sig_op*(sigx_op+sigy_op))
  xi_inv <- dynbetaXiInv(nsepx, nsepy, plane)
  xi_inv*scale
}

dynbetaRelBetaChange <- function(nsepx, nsepy, plane, energy, beta, tune, n_op, sigx_op, sigy_op, ref_nsep=c(10,0)) {
  # Return \Delta\Beta/\Beta at normalized separation (nsepx, nsepy) for a single plane and single beam
  #   nsepx = sepx/sigx; nsepy = ...
  #   energy = beam energy in GeV
  #   beta = beam nominal beta in plane `plane` in m
  #   tune = beam nominal tune in plane `plane`
  #   n_op = opposite beam charge
  #   sigx(y)_op = opposite beam sigma_x(y) in um
  #   ref_nsep is the separation of the reference beta (in units of (sigx, sigy))

  xi_xy <- dynbetaXi(nsepx, nsepy, plane, energy, beta, n_op, sigx_op, sigy_op)
  rdb_xy <- relBetaChangeFromXi(xi_xy, tune)
  if (is.null(ref_nsep)) {
    xi_ref <- rdb_ref <- NA
    reldb <- rdb_xy
  } else {
    xi_ref <- dynbetaXi(ref_nsep[1], ref_nsep[2], plane, energy, beta, n_op, sigx_op, sigy_op)
    rdb_ref <- relBetaChangeFromXi(xi_ref, tune)
    reldb <- (rdb_xy+1)/(rdb_ref+1)-1  # or approx rdb_xy-rdb_ref
  }
  data.frame(xi=xi_xy, xi_ref, rdb=rdb_xy, rdb_ref, reldb)
}

dynbetaRelLumiChange <- function(sepx, sepy, sigx1, sigx2, sigy1, sigy2, n1, n2, b1params, b2params, ref_nsep=c(10,0), method='diff', output='short') {
  # Return \Delta L/L at separation (sepx, sepy)
  #   sepx, sepy - separation in um
  #   sigx1, sigx2, sigy1, sigy2 - single beam widths in um
  #   n1, n2 - charge in bunches in units of e
  #   b1params, b2params - beam machine parameters: energy (GeV), betax(y) (m), tunex(y)
  #   ref_nsep is the separation of the reference beta (in units of (sigx, sigy))

  px1 <- dynbetaRelBetaChange(sepx/sigx1, sepy/sigy1, 'x', b1params$energy, b1params$betax, b1params$tunex, n2, sigx2, sigy2, ref_nsep) # X plane, beam 1
  px2 <- dynbetaRelBetaChange(sepx/sigx2, sepy/sigy2, 'x', b2params$energy, b2params$betax, b2params$tunex, n1, sigx1, sigy1, ref_nsep) # X plane, beam 2
  py1 <- dynbetaRelBetaChange(sepx/sigx1, sepy/sigy1, 'y', b1params$energy, b1params$betay, b1params$tuney, n2, sigx2, sigy2, ref_nsep) # Y plane, beam 1
  py2 <- dynbetaRelBetaChange(sepx/sigx2, sepy/sigy2, 'y', b2params$energy, b2params$betay, b2params$tuney, n1, sigx1, sigy1, ref_nsep) # Y plane, beam 2

  csigx <- sqrt(sigx1^2+sigx2^2)
  csigy <- sqrt(sigy1^2+sigy2^2)
  #d_csigx <- (px1$reldb*sigx1^2 + px2$reldb*sigx2^2)/(2*csigx)
  #d_csigy <- (py1$reldb*sigy1^2 + py2$reldb*sigy2^2)/(2*csigy)
  #csigx_mod <- csigx + d_csigx
  #csigy_mod <- csigy + d_csigy
  csigx_mod <- sqrt((1+px1$reldb)*sigx1^2 + (1+px2$reldb)*sigx2^2)
  csigy_mod <- sqrt((1+py1$reldb)*sigy1^2 + (1+py2$reldb)*sigy2^2)
  d_csigx <- csigx_mod - csigx
  d_csigy <- csigy_mod - csigy

  rdLdcsig <- function(sep, csig) ((sep/csig)^2 - 1)
  if (method == 'deriv') {
    # per plane effect
    dLL <- rdLdcsig(sepx, csigx)/csigx*d_csigx + rdLdcsig(sepy, csigy)/csigy*d_csigy
  } else if (method == 'diff') {
    L <- function(csigx, csigy) exp(-0.5*(sepx/csigx)^2 - 0.5*(sepy/csigy)^2)/(csigx*csigy)
    dLL <- L(csigx_mod, csigy_mod)/L(csigx, csigy) - 1
  } else if (method == 'wrong') {
    L <- function(csigx, csigy) 1/(csigx*csigy)
    dLL <- L(csigx_mod, csigy_mod)/L(csigx, csigy) - 1
  }

  if (output == 'full') {
    data.frame(px1=px1, px2=px2, py1=py1, py2=py2,
               csigx, csigy, d_csigx, d_csigy,
               rdLdcsigx=rdLdcsig(sepx, csigx), rdLdcsigy=rdLdcsig(sepy, csigy),
               dLL)
  } else {
    dLL
  }
}
## -------------------- End of minimal old beam-beam correction code from CERN gitlab --------------------

library(data.table, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)

dir <- commandArgs(TRUE)[1]
print(dir)
config <- read.table(file.path(dir, 'config.txt'), sep='=', col.names=c('var','value'),
                     stringsAsFactors=FALSE)
config <- as.data.table(config)[,var := trimws(var)]
setkey(config, var)
config.scalar <- function(name) as.numeric(config[J(name)]$value)
beta <- config.scalar('beta') # in m
p <- config.scalar('p') # in GeV
sig1x <- config.scalar('sig1.x') #  in um
sig1y <- config.scalar('sig1.y')
sig2x <- config.scalar('sig2.x')
sig2y <- config.scalar('sig2.y')
N1 <- config.scalar('N1')
N2 <- config.scalar('N2')
Qx <- config.scalar('Qx')
Qy <- config.scalar('Qy')
##
stopifnot(sig1x == sig1y & sig1x == sig2x & sig1x == sig2y)
sig <- sig1x

stopifnot(N1 == N2)
N <- N1

vdm.sig <- sqrt(2) * sig

orbitRelLumiChange <- function(sepx, sepy, sig, n1, n2, beta, Qx, Qy) {
    ## beta in m
    ## sepx,y, sig in um
    alpha <- 1/137.035
    hbar <- 0.197327e-15 # in Gev * m
    beta0 <- 1
    Z1 <- Z2 <- 1
    k <- 2 * Z1 * Z2 * (n1 + n2) * alpha * hbar / beta0 / p * beta * 1e12 # in um^2
    z <- sepx + 1i * sepy
    r2 <- Re(z * Conj(z))
    kick = k / 2 * z / r2 * (1 - exp(-r2 / 2 / vdm.sig^2))
    kick[z == 0] <- 0
    dx <- Re(kick) / tan(pi * Qx)
    dy <- Im(kick) / tan(pi * Qy)
    dLLx <- exp((-(sepx + dx)^2 + sepx^2) / 2 / vdm.sig^2)
    dLLy <- exp((-(sepy + dy)^2 + sepy^2) / 2 / vdm.sig^2)
    dLLx * dLLy
}

old.lumi.w.over.wo.fun <- function(sepx.range = seq(0, 5, by=0.05) * sig) {
    bparams <- list(energy=p,betax=beta,betay=beta,tunex=Qx,tuney=Qy)
    ## Take eg. the x-scan, the y-scan is the same as the beams are round
    beta.cor <- dynbetaRelLumiChange(sepx=sepx.range, sepy=0,
                                     sigx1=sig, sigx2=sig, sigy1=sig, sigy2=sig,
                                     n1=N,
                                     n2=N,
                                     b1params = bparams,
                                     b2params = bparams,
                                     ref_nsep=c(10,0))
    orbit.cor <- orbitRelLumiChange(sepx.range,0,sig,N,N,beta,Qx,Qy)
    old <- data.table(separation = sepx.range,
                      beta = 1 + beta.cor,
                      orbit = orbit.cor)
    old[, beta.and.orbit := beta * orbit]
    old[, vdm.profile := {
        . <- 2/sqrt(2*pi)/vdm.sig * exp(-separation^2 / 2 / vdm.sig^2)
    }]
}
old.lumi.w.over.wo <- old.lumi.w.over.wo.fun()

new.lumi.w.over.wo <- {
    . <- fread(file.path(dir, 'summary.txt'),
               col.names=c('step','x2','y2','int.to.int0.correction','int0.analytic',
                           'int0','int0.to.analytic', 'int0.rel.err',
                           'x','y',
                           'kick.x.analytic','kick.y.analytic',
                           'x.analytic','y.analytic'),
               colClasses=c('integer', rep('double',13)))
    .[, vdm.profile := {
        . <- 2/sqrt(2*pi)/vdm.sig * exp(-.$x2^2 / 2 / vdm.sig^2)
    }]
    . <- .[,list(separation = x2,
                 cor = 1 + 2 * (int.to.int0.correction - 1),
                 ## correction is doubled, as both bunches are affected by
                 ## int.to.int0.correction
                 vdm.profile)]
}

xsec <- {
    integ.to.ref <- function(x, y, cor) {
        ## trapezoidal integration is better than the sum which is biased, as
        ## it is equivalent to an integration by rectangles, and the left half
        ## of the first rectangle with the center at zero is in the negative
        ## region which in fact should be excluded.        
        trapezoidal.integ <- pracma::trapz(x, y * cor)
        reference <- pracma::trapz(x, y)
        trapezoidal.integ / reference
    }
    ## formula for x-section: integ^2 / mu0:
    xsec.cor <- function(x, y, cor) {
        int.w.divided.by.wo <- integ.to.ref(x, y, cor)
        int.w.divided.by.wo^2 / cor[1]
    }
    old.xsec <- with(old.lumi.w.over.wo, xsec.cor(separation, vdm.profile, beta.and.orbit))
    new.xsec <- with(new.lumi.w.over.wo, xsec.cor(separation, vdm.profile, cor))
    list(old = old.xsec, new = new.xsec, new.minus.old = new.xsec - old.xsec)
}

gaussian.fit.estimator <- function(h) {
    nev <- sum(h$y)
    mean <- sum(h$y * h$x) / nev
    sd <- sqrt(sum(h$y * (h$x-mean)^2) / nev)
    list(nev=nev, mean=mean, sd=sd)
}
## Fit to single Gaussian.
## If not given, initial "start" values are calculated from the number of events, mean() and sd() 
hfit.g <- function(h, start = gaussian.fit.estimator(h), ...) {
    bin.width <- h$breaks[2] - h$breaks[1]
    error.squared <- pmax(1, h$y) # assign error=1 to bins with content<1 (eg. zero)
    x <- h$x
    y <- h$y
    suppressWarnings(
        nls(y ~ bin.width * nev * dnorm(x, mean, sd), weights=1/error.squared,
            start=start, ...))
}

old.reco.from.new <- {
    bparams <- list(energy=p,betax=beta,betay=beta,tunex=Qx,tuney=Qy)
    n <- copy(new.lumi.w.over.wo)
    nbins <- nrow(n)
    h <- with(n, {
        bin <- diff(separation[1:2])
        y <- cor * vdm.profile * 1e10
        y <- c(y[nbins:1],y[2:nbins])
        list(breaks=c(-separation[nbins:1]-bin/2, separation+bin/2),
             x=c(-separation[nbins:1],separation[2:nbins]),
             y=rpois(2*nbins-1, y))
    })
    f <- hfit.g(h, control = list(warnOnly = TRUE))
    ## Warnings possible but convergence is still good enough
    ##    coef(summary(f))[3,1:2]
    ## without beam-beam:
    h0 <- with(n, {
        bin <- diff(separation[1:2])
        y <- vdm.profile * 1e10
        y <- c(y[nbins:1],y[2:nbins])
        list(breaks=c(-separation[nbins:1]-bin/2, separation+bin/2),
             x=c(-separation[nbins:1],separation[2:nbins]),
             y=rpois(2*nbins-1, y))
    })
    f0 <- hfit.g(h0, control = list(warnOnly = TRUE))
    ##    coef(summary(f0))[3,1:2]
    sig.new <- as.numeric(sig * coef(f)[3] / coef(f0)[3])
    ## Take eg. the x-scan, the y-scan is the same as the beams are round
    beta.cor <- dynbetaRelLumiChange(sepx=n$separation, sepy=0,
                                     sigx1=sig.new, sigx2=sig.new, sigy1=sig.new, sigy2=sig.new,
                                     n1=N,
                                     n2=N,
                                     b1params = bparams,
                                     b2params = bparams,
                                     ref_nsep=c(10,0))
    orbit.cor <- orbitRelLumiChange(n$separation,0,sig.new,N,N,beta,Qx,Qy)
    old <- n[, `:=`(beta = 1 + beta.cor,
                    orbit = orbit.cor)]
    old[, beta.and.orbit := beta * orbit]
    old[, reco := cor / beta.and.orbit]
    list(old,
         sig.reco.to.sig = sig.new / sig)
}

xsec.old.from.new <- {
    integ.to.ref <- function(x, y, cor) {
        ## trapezoidal integration is better than the sum which is biased, as
        ## it is equivalent to an integration by rectangles, and the left half
        ## of the first rectangle with the center at zero is in the negative
        ## region which in fact should be excluded.        
        trapezoidal.integ <- pracma::trapz(x, y * cor)
        reference <- pracma::trapz(x, y)
        trapezoidal.integ / reference
    }
    ## formula for x-section: integ^2 / mu0:
    xsec.cor <- function(x, y, cor) {
        int.w.divided.by.wo <- integ.to.ref(x, y, cor)
        int.w.divided.by.wo^2 / cor[1]
    }
    with(old.reco.from.new[][[1]], xsec.cor(separation, vdm.profile, reco))
}

integ <- fread(file.path(dir, 'integrals.txt.gz'),
          col.names=c('step','x2','y2','turn','int'))
integ0 <- integ[turn < 1000, list(int0 = mean(int)), by = step]

summary <- fread(file.path(dir, 'summary.txt'),
                 col.names=c('step','x2','y2','int.to.int0.correction','int0.analytic',
                             'int0','int0.to.analytic', 'int0.rel.err',
                             'x','y',
                             'kick.x.analytic','kick.y.analytic',
                             'x.analytic','y.analytic'),
                 colClasses=c('integer', rep('double',13)))
summary[, correction := 200 * (int.to.int0.correction - 1)]

center <- fread(file.path(dir, 'centers.txt.gz'),
                col.names=c('step','x2','y2','turn','x','y'))

rxy.weight <- fread(file.path(dir, 'rx_ry_weights.txt.gz'),
                     col.names=c('step','point','rx','ry','w'))

## -------------------- Plots --------------------
{
    txt.16 <- element_text(size = 16)
    bw1 <- theme_bw() + theme(title = txt.16, axis.title = txt.16,
                              legend.title = txt.16, legend.text = txt.16)
    ## use black-white + font size 16 theme by default
    theme_set(theme_get() + bw1 + theme(axis.text = txt.16))
}
gg.cor <- function() {
    old <- melt(old.lumi.w.over.wo[,list(separation, beta, orbit, beta.and.orbit)], id.vars = 'separation')
    combined <- rbind(old,
                      new.lumi.w.over.wo[, list(separation,
                                                variable = 'new',
                                                value = cor)])
     qplot(data = combined, separation, value, geom=c('point','line'), color=variable) +
        geom_hline(aes(yintercept = 1), linetype = 'dashed') +
        labs(x = 'Beams separation, um',
             y = 'L(with bb) / L(wo bb)')
}

gg.vdm.profiles <- function() {
    . <- old.lumi.w.over.wo[, {
        list(separation = separation,
             old = vdm.profile * beta.and.orbit,
             wo.beam.beam = vdm.profile)
    }]
    combined <- rbind(melt(., id.vars = 'separation'),
                      new.lumi.w.over.wo[, list(separation,
                                                value = vdm.profile * cor,
                                                variable = 'new')])
     qplot(data = combined, separation, value, geom=c('point','line'), color=variable) +
        labs(x = 'Beams separation, um',
             y = 'vdM profiles, normalized without beam-beam')
}

gg.vdm.profile.change <- function() {
    combined <- rbind(old.lumi.w.over.wo[, list(separation = separation,
                                                  y = vdm.profile * (beta.and.orbit - 1),
                                                  variable = 'old')],
                      new.lumi.w.over.wo[, list(separation,
                                                y = vdm.profile * (cor - 1),
                                                variable = 'new')])
     qplot(data = combined, separation, y, geom=c('point','line'), color=variable) +
        geom_hline(aes(yintercept = 0), linetype = 'dashed') +
        labs(x = 'Beams separation, um',
             y = 'vdM profiles (with - normalized wo) beam-beam')
}

gg.old.from.new <- function() {
    . <- melt(old.reco.from.new[[1]], id.vars=c('separation','vdm.profile'))
    pl1 <- qplot(data=., separation, value, color=variable,geom=c('line','point')) +
        geom_line(data = .[variable == 'reco'],
                  aes(separation, value, color=variable), size=2) +
        labs(x = 'Bem separation, um', y = 'L(with bb) / L(wo bb)')
    pl2 <- qplot(data=., separation, vdm.profile*(value - 1), color=variable,geom=c('line','point')) +
        geom_line(data = .[variable == 'reco'],
                  aes(separation, vdm.profile*(value - 1), color=variable), size=2) +
        labs(x = 'Bem separation, um', y = 'vdM profiles (with - normalized wo) beam-beam')
    list(pl1, pl2)
}

gg.integ.per.turn <- function() {
    . <- merge(integ, integ0, by='step')[, correction := int / int0]
    ## integrals per turn
    qplot(data=., turn, correction, size = I(0.01)) +
        facet_wrap(~ reorder(paste0(x2,',',y2),
                             sort(as.numeric(paste0(x2,y2))))
                   ) +
        geom_hline(aes(yintercept = 1), linetype = 'dashed')
}

gg.integ.per.100.turns <- function() {
    . <- copy(integ)[, turn := (turn %/% 100) * 100]
    . <- .[, list(int = mean(int)), by = .(step, x2, y2, turn)]
    . <- merge(., integ0, by='step')
    .[, correction := int / int0]
    ## integrals averaged per 100 turns
    qplot(data=., turn, correction) +
        facet_wrap(~ reorder(paste0(x2,',',y2),
                             sort(as.numeric(paste0(x2,y2))))
                   ) +
        geom_hline(aes(yintercept = 1))
}

gg.correction <- function() {
    qplot(data=summary, reorder(paste0(x2,',',y2),
                                  sort(as.numeric(paste0(x2,y2)))),
          correction, geom = c('point','line'), group=1) +
        geom_hline(aes(yintercept = 0), linetype = 'dashed') +
        labs(x = 'Beam separation, um', y = '(L(w bb) / L(wo bb) - 1) * 100%') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

gg.center <- function() {
    . <- melt(summary, id.vars=c('x2','y2'),
              measure.vars=c('x','x.analytic','y','y.analytic'))
    .[, axis := ifelse(grepl('^x',variable), 'x', 'y')]
    .[, method := ifelse(grepl('analytic',variable), 'analytic', 'numeric')]
    ggplot() +
        geom_point(data=.[method=='numeric'],
                   aes(x=reorder(paste0(x2,',',y2),
                                 sort(as.numeric(paste0(x2,y2)))),
                       y=value)) +
        geom_line(data=.,
                  aes(x=reorder(paste0(x2,',',y2),
                                sort(as.numeric(paste0(x2,y2)))),
                      y=value,
                      group=interaction(axis,method), color=axis,
                      linetype = method)) +
        labs(x = 'Beam separation, um', y = '<X>, um') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot.dir <- file.path(dir, 'pdf')
if (! file.exists( plot.dir )) {
    if (! dir.create(plot.dir, recursive = TRUE) ) {
        stop('Can not create subdirectory ', plot.dir)
    }
}

cat('------------------------------------------------------------\n')
cat('PDF files will be stored in', plot.dir, '\n')
cat('------------------------------------------------------------\n')
i.plot <- 0
save.plot <- function(plot.function.name) {
    pdf(file.path(plot.dir, paste0(plot.function.name, '.pdf')), width=12, height=8)
    . <- get(plot.function.name)()
    if ( is.list(.) & all(
                          sapply(., function(x) 'ggplot' %in% class(x))
                      )) {
        lapply(., print)
    } else
        print(.)
    i.plot <<- i.plot + 1
    cat('PDF file', i.plot, '\n')
    invisible(dev.off())
}
save.plot('gg.cor')
save.plot('gg.vdm.profiles')
save.plot('gg.vdm.profile.change')
save.plot('gg.old.from.new')
save.plot('gg.integ.per.turn')
save.plot('gg.integ.per.100.turns')
save.plot('gg.correction')
save.plot('gg.center')

data.dir <- file.path(dir, 'Rdata')
if (! file.exists( data.dir )) {
    if (! dir.create(data.dir, recursive = TRUE) ) {
        stop('Can not create subdirectory ', data.dir)
    }
}
save(old.lumi.w.over.wo,
     new.lumi.w.over.wo,
     xsec,
     old.reco.from.new,
     xsec.old.from.new,
     center,
     summary,
     integ,
     file = file.path(data.dir, 'results.rds'))

print('Cross-section corrections: Xsec(beam-beam)/Xsec(NO beam-beam)')
print(xsec)
print('Cross-section when applying old correction on top of new data')
print(xsec.old.from.new)
