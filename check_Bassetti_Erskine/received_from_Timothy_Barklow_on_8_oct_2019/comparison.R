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
library(RcppFaddeeva)

splitIdColumn <- function(df, pattern, groups, variable.name='variable', keep.var=F) {
  vars <- str_match(df[, variable.name], pattern)[, groups, drop=F]
  ids <- data.frame(alply(vars, 2, type.convert, as.is=T), stringsAsFactors=F)
  names(ids) <- names(groups)
  newdf <- cbind(df, ids)
  if (!keep.var) newdf[variable.name] <- NULL
  newdf
}

C_BBP <- 1.44e-12/(2*pi)  # alpha*hbar*c/(2*pi) in GeV.um
M_P <- 0.938272  # proton mass in GeV

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

args <- commandArgs(TRUE)

useAveragePoints <- length(args) > 1
ylimEnabled <- length(args) > 3

cat(' useAveragePoints= ', useAveragePoints, '\n')

useBaseErskForOrbit <- TRUE

cat(' useBaseErskForOrbit= ', useBaseErskForOrbit, '\n')


#dir <- commandArgs(TRUE)[1]
dir <- args[1]
print(dir)
dirPoints <- ifelse(useAveragePoints,args[2],"")
ylimMin <- ifelse(ylimEnabled,as.double(args[3]),0.97)
ylimMax <- ifelse(ylimEnabled,as.double(args[4]),1.03)
print(dirPoints)
config <- read.table(file.path(dir, 'config.txt'), sep='=', col.names=c('var','value'),
                     stringsAsFactors=FALSE)
config <- as.data.table(config)[,var := trimws(var)]
setkey(config, var)
config.scalar <- function(name) as.numeric(config[J(name)]$value)
config.vector <- function(name) as.numeric(unlist(strsplit(trimws(config[J(name)]$value,which="left"), split=" ")))
beta <- config.scalar('beta') # in m
kvbbd <- beta*1e6
p <- config.scalar('p') # in GeV
gamma <- p/M_P
rProton <- 1.5346955e-18
rpOverGamma <- rProton/gamma
sig1x <- config.scalar('sig1.x') #  in um
sig1y <- config.scalar('sig1.y')
sig2x <- config.scalar('sig2.x')
sig2y <- config.scalar('sig2.y')
umRoundBeamOffset <- 1e-8
capsig1x <- sqrt(2)*sig1x+umRoundBeamOffset
capsig1y <- sqrt(2)*sig1y
capsig2x <- sqrt(2)*sig2x+umRoundBeamOffset
capsig2y <- sqrt(2)*sig2y
Z1 <- config.scalar('Z1')
Z2 <- config.scalar('Z2')
N1 <- config.scalar('N1')
N2 <- config.scalar('N2')
Qx <- config.scalar('Qx')
Qy <- config.scalar('Qy')
config.scalar('sig1.x') #  in um
config.scalar('sig1.y')
config.scalar('sig2.x')
config.scalar('sig2.y')
bsigx2 <- config.vector('bsigx2')
bsigy2 <- config.vector('bsigy2')
x2config <- config.vector('x2')
y2config <- config.vector('y2')
if(is.na(x2config[1])) {
   x2config <- bsigx2*sig2x
   y2config <- bsigy2*sig2y
} else {
   bsigx2 <- x2config/sig2x
   bsigy2 <- y2config/sig2y
}
cat(' initialization bsigx2= ',bsigx2, '\n')
cat(' initialization x2config= ',x2config, '\n')
cat(' initialization bsigy2= ',bsigy2, '\n')
cat(' initialization y2config= ',y2config, '\n')
#cat(' is.na(bsigx2[1])= ', is.na(bsigx2[1]), '\n')
#cat(' is.na(x2config[1])= ', is.na(x2config[1]), '\n')
#cat(' is.na(bsigy2[1])= ', is.na(bsigy2[1]), '\n')
#cat(' is.na(y2config[1])= ', is.na(y2config[1]), '\n')
bsigx2 > bsigy2
isSepx <- sum(bsigx2 > bsigy2) > 0
cat(' initialization isSepx= ',isSepx, '\n')
#useDefaultSepxRange <- TRUE
useDefaultSepxRange <- FALSE
cat(' useDefaultSepxRange= ',useDefaultSepxRange, '\n')

##
alpha <- 1/137.035
hbar <- 0.197327e-15 # hbarc in Gev * m
beta0 <- 1
##
#stopifnot(sig1x == sig1y & sig1x == sig2x & sig1x == sig2y)
stopifnot(sig1x == sig2x & sig1y == sig2y)
sig <- 0.5*(sig2x+sig2y)

stopifnot(N1 == N2)
N <- N1


BaseErsk <- function(sepx,sepy) {
# convert all dimensions from um to meters and take abs value of deflections
  sepxm <- abs(sepx*1e-6)
  sepym <- abs(sepy*1e-6)
  sigmaxm <- capsig2x*1e-6
  sigmaym <- capsig2y*1e-6
  eps0 <- 1
  sep1 <- sepxm
  sep2 <- sepym
  sigma1 <- sigmaxm
  sigma2 <- sigmaym
  if(sigmaxm < sigmaym) {
    sep1 <- sepym
    sep2 <- sepxm
    sigma1 <- sigmaym
    sigma2 <- sigmaxm
  }
  SS <- sqrt(2*(sigma1*sigma1-sigma2*sigma2))
  factBE <- 2*sqrt(pi)/(2*eps0*SS)
  etaBE <- complex(real=sep1*sigma2/sigma1, imaginary=sep2*sigma1/sigma2)
  zetaBE <- complex(real=sep1, imaginary=sep2)
  w1 <- Faddeeva_w(zetaBE/SS)
  w2 <- Faddeeva_w(etaBE/SS)
  valBE <- factBE*(w1-exp(-sep1*sep1/(2.*sigma1*sigma1)-sep2*sep2/(2.*sigma2*sigma2))*w2)
  Ex <- 0
  Ey <- 0
  if(sigmaxm > sigmaym) {
    Ex <- abs(Im(valBE))*sign(sepx)
    Ey <- abs(Re(valBE))*sign(sepy)
  } else {
    Ey <- abs(Im(valBE))*sign(sepy)
    Ex <- abs(Re(valBE))*sign(sepx)
  }
  complex(real=Ex, imaginary=Ey)
}

EllipDeflectAngles <- function(sepx,sepy) {
  KK <- 2*rpOverGamma
  Efield <- BaseErsk(sepx,sepy)
  Ex <- Re(Efield)
  Ey <- Im(Efield)
  dfangx <- KK*N2*Ex
  dfangy <- KK*N2*Ey
  complex(real=dfangx, imaginary=dfangy)
}

bbpoints <- function(dirFunc,sepx,sepy) {
  . <- fread(file.path(dirFunc, 'summary.txt'),
             col.names=c('step','x2','y2','int.to.int0.correction','int0.analytic',
                         'int0','int0.to.analytic', 'int0.rel.err',
                         'x','y',
                         'kick.x.analytic','kick.y.analytic',
                         'x.analytic','y.analytic'),
             colClasses=c('integer', rep('double',13)))
  if(isSepx) {
    x2vals <- .[,x2]
    int.to.int0vals <- .[,int.to.int0.correction]
#    cat(' x2vals= ',x2vals,'\n')
#    cat(' int.to.int0vals= ',int.to.int0vals,'\n')
#    cat(' cor= ',1 + 2 * (int.to.int0vals - 1),'\n')
                    ## correction is doubled, as both bunches are affected by
                    ## int.to.int0.correction
    orbFun <- approxfun(x2vals,1 + 2 * (int.to.int0vals - 1))
    orbFun(sepx)
  } else {
    y2vals <- .[,y2]
    int.to.int0vals <- .[,int.to.int0.correction]
#    cat(' y2vals= ',y2vals,'\n')
#    cat(' int.to.int0vals= ',int.to.int0vals,'\n')
#    cat(' cor= ',1 + 2 * (int.to.int0vals - 1),'\n')
                    ## correction is doubled, as both bunches are affected by
                    ## int.to.int0.correction
    orbFun <- approxfun(y2vals,1 + 2 * (int.to.int0vals - 1))
    orbFun(sepy)
  }
}

opticalDistortions <- function(sepx,sepy) {
  bbpoints(dir,sepx,sepy)-bbpoints(dirPoints,sepx,sepy)+1
}



if(useAveragePoints) {
  capSigma <- {
    . <- fread(file.path(dir, 'summary.txt'),
               col.names=c('step','x2','y2','int.to.int0.correction','int0.analytic',
                           'int0','int0.to.analytic', 'int0.rel.err',
                           'x','y',
                           'kick.x.analytic','kick.y.analytic',
                           'x.analytic','y.analytic'),
               colClasses=c('integer', rep('double',13)))
    if(isSepx) {
      x2vals0 <- .[,x2]
      x2vals <- seq(x2vals0[1],tail(x2vals0,n=1),length.out=100000)
      L0vals <- exp(-0.25*(x2vals/sig2x)^2)
      capSig0 <- sum(L0vals)
      cat(' capSig0= ',capSig0,'\n')
      corvalsBB <- bbpoints(dir,x2vals,0)
      capSigBB <- sum(corvalsBB*L0vals)
      capSigRatioBB <- capSigBB/capSig0/corvalsBB[1]
      cat(' capSigBB= ',capSigBB,'\n')
      cat(' capSigRatioBB = ',capSigRatioBB,'\n')
      cat(' capSigRatioBB percent change= ',(capSigRatioBB-1)*100,'\n')
      muPeakRatioBB=corvalsBB[1]
      cat(' muPeakRatioBB = ',muPeakRatioBB,'\n')
      cat(' muPeakRatioBB percent change= ',(muPeakRatioBB-1)*100,'\n')
      corvalsOrbit <- bbpoints(dirPoints,x2vals,0)
      capSigOrbit <- sum(corvalsOrbit*L0vals)
      capSigRatioOrbit <- capSigOrbit/capSig0/corvalsOrbit[1]
      cat(' capSigOrbit= ',capSigOrbit,'\n')
      cat(' capSigRatioOrbit = ',capSigRatioOrbit,'\n')
      cat(' capSigRatioOrbit percent change= ',(capSigRatioOrbit-1)*100,'\n')
      muPeakRatioOrbit=corvalsOrbit[1]
      cat(' muPeakRatioOrbit = ',muPeakRatioOrbit,'\n')
      cat(' muPeakRatioOrbit percent change= ',(muPeakRatioOrbit-1)*100,'\n')
      corvalsOptic <- opticalDistortions(x2vals,0)
      capSigOptic <- sum(corvalsOptic*L0vals)
      capSigRatioOptic <- capSigOptic/capSig0/corvalsOptic[1]
      cat(' capSigOptic= ',capSigOptic,'\n')
      cat(' capSigRatioOptic = ',capSigRatioOptic,'\n')
      cat(' capSigRatioOptic percent change= ',(capSigRatioOptic-1)*100,'\n')
      muPeakRatioOptic=corvalsOptic[1]
      cat(' muPeakRatioOptic = ',muPeakRatioOptic,'\n')
      cat(' muPeakRatioOptic percent change= ',(muPeakRatioOptic-1)*100,'\n')
    } else {
      y2vals0 <- .[,y2]
      y2vals <- seq(y2vals0[1],tail(y2vals0,n=1),length.out=100000)
      L0vals <- exp(-0.25*(y2vals/sig2y)^2)
      capSig0 <- sum(L0vals)
      cat(' capSig0= ',capSig0,'\n')
      corvalsBB <- bbpoints(dir,0,y2vals)
      capSigBB <- sum(corvalsBB*L0vals)
      capSigRatioBB <- capSigBB/capSig0/corvalsBB[1]
      cat(' capSigBB= ',capSigBB,'\n')
      cat(' capSigRatioBB = ',capSigRatioBB,'\n')
      cat(' capSigRatioBB percent change= ',(capSigRatioBB-1)*100,'\n')
      muPeakRatioBB=corvalsBB[1]
      cat(' muPeakRatioBB = ',muPeakRatioBB,'\n')
      cat(' muPeakRatioBB percent change= ',(muPeakRatioBB-1)*100,'\n')
      corvalsOrbit <- bbpoints(dirPoints,0,y2vals)
      capSigOrbit <- sum(corvalsOrbit*L0vals)
      capSigRatioOrbit <- capSigOrbit/capSig0/corvalsOrbit[1]
      cat(' capSigOrbit= ',capSigOrbit,'\n')
      cat(' capSigRatioOrbit = ',capSigRatioOrbit,'\n')
      cat(' capSigRatioOrbit percent change= ',(capSigRatioOrbit-1)*100,'\n')
      muPeakRatioOrbit=corvalsOrbit[1]
      cat(' muPeakRatioOrbit = ',muPeakRatioOrbit,'\n')
      cat(' muPeakRatioOrbit percent change= ',(muPeakRatioOrbit-1)*100,'\n')
      corvalsOptic <- opticalDistortions(0,y2vals)
      capSigOptic <- sum(corvalsOptic*L0vals)
      capSigRatioOptic <- capSigOptic/capSig0/corvalsOptic[1]
      cat(' capSigOptic= ',capSigOptic,'\n')
      cat(' capSigRatioOptic = ',capSigRatioOptic,'\n')
      cat(' capSigRatioOptic percent change= ',(capSigRatioOptic-1)*100,'\n')
      muPeakRatioOptic=corvalsOptic[1]
      cat(' muPeakRatioOptic = ',muPeakRatioOptic,'\n')
      cat(' muPeakRatioOptic percent change= ',(muPeakRatioOptic-1)*100,'\n')
    }
  }
}

orbitRelLumiChange <- function(sepx, sepy, sig, n1, n2, beta, Qx, Qy) {
  if(useAveragePoints) {
    bbpoints(dirPoints,sepx,sepy)
  } else {
    ## beta in m
    ## sepx,y, sig in um
    vdm.s <- sqrt(2) * sig
    vdm.sx <- sqrt(2) * sig2x
    vdm.sy <- sqrt(2) * sig2y
    k <- 2 * Z1 * Z2 * (n1 + n2) * alpha * hbar / beta0 / p * beta * 1e12 # in um^2
    z <- sepx + 1i * sepy
    r2 <- Re(z * Conj(z))
    kick = k / 2 * z / r2 * (1 - exp(-r2 / 2 / vdm.s^2))
    kick[z == 0] <- 0
    kickBE = kvbbd*EllipDeflectAngles(Re(z),Im(z))
    if(useBaseErskForOrbit) {
      kick <- kickBE
    }
    cat(' Re(kick)= ', Re(kick), '\n')
    cat(' Re(kickBE)= ', Re(kickBE), '\n')
    cat(' Im(kick)= ', Im(kick), '\n')
    cat(' Im(kickBE)= ', Im(kickBE), '\n')
    dx <- Re(kick) / tan(pi * Qx)
    dy <- Im(kick) / tan(pi * Qy)
    dLLx <- exp((-(sepx + dx)^2 + sepx^2) / 2 / vdm.sx^2)
    dLLy <- exp((-(sepy + dy)^2 + sepy^2) / 2 / vdm.sy^2)
    dLLx * dLLy
  }
}

vdm.sig <- sqrt(2) * sig
## quadrupole kick approximation:
k.angle <- 2 * Z1 * Z2 * N * alpha * hbar / beta0 / p * 1e12

beta.linearized.rel.lumi.change <- function(x) {
    kick.x.deriv<- function(x) {
        ifelse(x==0,
               k.angle/2/sig^2,
               k.angle * (-1 / x^2 + (1 / sig^2 + 1 / x^2) *exp(-x^2 / 2/sig^2)))
    }
    kick.y.deriv <- function(x) {
        ifelse(x==0,
               k.angle/2/sig^2,
               k.angle / x^2 *(1 - exp(-x^2 / 2/sig^2)))
    }
    rel.beta.x <- function(x) {
        beta * kick.x.deriv(x) / 2 / tan(2*pi*Qx)
    }
    rel.beta.y <- function(x) {
        beta * kick.y.deriv(x) / 2 / tan(2*pi*Qy)
    }
    vdm.sigx_mod <- (1 + rel.beta.x(x)/2) * vdm.sig
    vdm.sigy_mod <- (1 + rel.beta.y(x)/2) * vdm.sig
    L <- function(vdm.sigx, vdm.sigy) exp(-0.5 * (x/vdm.sigx)^2) / vdm.sigx / vdm.sigy
    L(vdm.sigx_mod, vdm.sigy_mod)/L(vdm.sig, vdm.sig) - 1
}
## ----------------------------------
old.lumi.w.over.wo.fun <- function(sepx.range = seq(0, 5, by=0.05) * sig) {
    sepxRange <- `if`(useDefaultSepxRange,`if`(isSepx,sepx.range,0),x2config)
    sepyRange <- `if`(useDefaultSepxRange,`if`(isSepx,0,sepx.range),y2config)
    bparams <- list(energy=p,betax=beta,betay=beta,tunex=Qx,tuney=Qy)
    ## Take eg. the x-scan, the y-scan is the same as the beams are round
#    beta.cor <- dynbetaRelLumiChange(sepx=sepx.range, sepy=0,
    beta.cor <- dynbetaRelLumiChange(sepx=sepxRange, sepy=sepyRange,
                                     sigx1=sig, sigx2=sig, sigy1=sig, sigy2=sig,
                                     n1=N,
                                     n2=N,
                                     b1params = bparams,
                                     b2params = bparams,
                                     ref_nsep=c(10,0))
#    orbit.cor <- orbitRelLumiChange(sepx.range,0,sig,N,N,beta,Qx,Qy)
    orbit.cor <- orbitRelLumiChange(sepxRange,sepyRange,sig,N,N,beta,Qx,Qy)
    if(useAveragePoints) {
      optic.distort <- opticalDistortions(sepxRange,sepyRange)
      old <- data.table(separation = `if`(isSepx,sepxRange,sepyRange),
                        beta.MADX = 1 + beta.cor,
                        orbit = orbit.cor,
  		      optic = optic.distort)
  
      old[, beta.linear := 1 + beta.linearized.rel.lumi.change(separation)]
      old[, beta.MADX.orbit := beta.MADX * orbit]
      old[, vdm.profile := {
          . <- 2/sqrt(2*pi)/vdm.sig * exp(-separation^2 / 2 / vdm.sig^2)
      }]
    } else {
      old <- data.table(separation = `if`(isSepx,sepxRange,sepyRange),
                        beta.MADX = 1 + beta.cor,
                        orbit = orbit.cor)
  
      old[, beta.linear := 1 + beta.linearized.rel.lumi.change(separation)]
      old[, beta.MADX.orbit := beta.MADX * orbit]
      old[, vdm.profile := {
          . <- 2/sqrt(2*pi)/vdm.sig * exp(-separation^2 / 2 / vdm.sig^2)
      }]
    }
}
old.lumi.w.over.wo <- old.lumi.w.over.wo.fun()

#old.lumi.w.over.wo.alt.range.fun <- function(sepx.range = seq(0, 432, by=21.6) * 1) {
old.lumi.w.over.wo.alt.range.fun <- function(sepx.range=c(0, 22, 43, 65, 86, 108,  432)) {
    sepxRangeAlt <- `if`(isSepx,sepx.range,0)
    sepyRangeAlt <- `if`(isSepx,0,sepx.range)
    bparams <- list(energy=p,betax=beta,betay=beta,tunex=Qx,tuney=Qy)
    ## Take eg. the x-scan, the y-scan is the same as the beams are round
    beta.cor <- dynbetaRelLumiChange(sepx=sepxRangeAlt, sepy=sepyRangeAlt,
                                     sigx1=sig, sigx2=sig, sigy1=sig, sigy2=sig,
                                     n1=N,
                                     n2=N,
                                     b1params = bparams,
                                     b2params = bparams,
                                     ref_nsep=c(10,0))
    orbit.cor <- orbitRelLumiChange(sepxRangeAlt,sepyRangeAlt,sig,N,N,beta,Qx,Qy)
    old <- data.table(separation = `if`(isSepx,sepxRangeAlt,sepyRangeAlt),
                      beta.MADX = 1 + beta.cor,
                      orbit = orbit.cor)

    old[, beta.linear := 1 + beta.linearized.rel.lumi.change(separation)]
    old[, beta.MADX.orbit := beta.MADX * orbit]
    old[, vdm.profile := {
        . <- 2/sqrt(2*pi)/vdm.sig * exp(-separation^2 / 2 / vdm.sig^2)
    }]
}
old.lumi.w.over.wo.alt.range <- old.lumi.w.over.wo.alt.range.fun()

new.lumi.w.over.wo <- {
    . <- fread(file.path(dir, 'summary.txt'),
               col.names=c('step','x2','y2','int.to.int0.correction','int0.analytic',
                           'int0','int0.to.analytic', 'int0.rel.err',
                           'x','y',
                           'kick.x.analytic','kick.y.analytic',
                           'x.analytic','y.analytic'),
               colClasses=c('integer', rep('double',13)))
    if(isSepx) {
       .[, vdm.profile := {
           . <- 2/sqrt(2*pi)/vdm.sig * exp(-.$x2^2 / 2 / vdm.sig^2)
       }]
       . <- .[,list(separation = x2,
                    cor = 1 + 2 * (int.to.int0.correction - 1),
                    ## correction is doubled, as both bunches are affected by
                    ## int.to.int0.correction
                    vdm.profile)]
    } else {
       .[, vdm.profile := {
           . <- 2/sqrt(2*pi)/vdm.sig * exp(-.$y2^2 / 2 / vdm.sig^2)
       }]
       . <- .[,list(separation = y2,
                    cor = 1 + 2 * (int.to.int0.correction - 1),
                    ## correction is doubled, as both bunches are affected by
                    ## int.to.int0.correction
                    vdm.profile)]
    }
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
    old.xsec <- with(old.lumi.w.over.wo, xsec.cor(separation, vdm.profile, beta.MADX.orbit))
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
#    cat(' old.reco.from.new  separation= ', separation, '\n')
    cat(' old.reco.from.new  nbins= ', nbins, '\n')
    h <- with(n, {
        bin <- diff(separation[1:2])
        y <- cor * vdm.profile * 1e10
        y <- c(y[nbins:1],y[2:nbins])
        list(breaks=c(-separation[nbins:1]-bin/2, separation+bin/2),
             x=c(-separation[nbins:1],separation[2:nbins]),
             y=rpois(2*nbins-1, y))
    })
    cat(' old.reco.from.new point 001 ','\n')
    f <- hfit.g(h, control = list(warnOnly = TRUE))
    cat(' old.reco.from.new point 002 ','\n')
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
    cat(' old.reco.from.new point 003 ','\n')
    f0 <- hfit.g(h0, control = list(warnOnly = TRUE))
    cat(' old.reco.from.new point 004 ','\n')
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
    old <- n[, `:=`(beta.MADX = 1 + beta.cor,
                    orbit = orbit.cor)]
    old[, beta.MADX.orbit := beta.MADX * orbit]
    old[, reco := cor / beta.MADX.orbit]
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
  if(useAveragePoints) {
    old <- melt(old.lumi.w.over.wo[,list(separation, orbit, optic)],
                id.vars = 'separation')
    combined <- rbind(old,
                  new.lumi.w.over.wo[, list(separation,
                                            variable = 'new',
                                            value = cor)])
    if(ylimEnabled) {
      ggp <- ggplot() +
#        ylim(0.97,1.03) +
        ylim(ylimMin,ylimMax) +
        geom_line (data = combined[variable!='new'], aes(separation, value, color=variable)) +
        geom_point(data = combined[variable=='new'], aes(separation, value, color=variable), size=2) +
        geom_line(data = combined[variable=='new'], aes(separation, value, color=variable), linetype = 'dotted') +
        geom_hline(aes(yintercept = 1), linetype = 'dashed') +
        scale_colour_manual(values = c(new='black', orbit='green', optic='purple')) +
        labs(x = 'Beams separation, um',
             y = 'L(with bb) / L(wo bb)')
      pggp <- ggplot_build(ggp)
    } else {
      ggp <- ggplot() +
        geom_line (data = combined[variable!='new'], aes(separation, value, color=variable)) +
        geom_point(data = combined[variable=='new'], aes(separation, value, color=variable), size=2) +
        geom_line(data = combined[variable=='new'], aes(separation, value, color=variable), linetype = 'dotted') +
        geom_hline(aes(yintercept = 1), linetype = 'dashed') +
        scale_colour_manual(values = c(new='black', orbit='green', optic='purple')) +
        labs(x = 'Beams separation, um',
             y = 'L(with bb) / L(wo bb)')
      pggp <- ggplot_build(ggp)
    }
    
  } else {
    old <- melt(old.lumi.w.over.wo[,list(separation, beta.MADX, orbit, beta.MADX.orbit, beta.linear)],
                id.vars = 'separation')
    combined <- rbind(old,
                  new.lumi.w.over.wo[, list(separation,
                                            variable = 'new',
                                            value = cor)])
#    print(getOption("max.print"))					    
#    print('from inside gg.cor print combined')
#    print(combined,max=999)
    ggp <- ggplot() +
        geom_line (data = combined[variable!='new'], aes(separation, value, color=variable)) +
        geom_point(data = combined[variable=='new'], aes(separation, value, color=variable), size=2) +
        geom_line(data = combined[variable=='new'], aes(separation, value, color=variable), linetype = 'dotted') +
        geom_hline(aes(yintercept = 1), linetype = 'dashed') +
        scale_colour_manual(values = c(new='black', beta.MADX.orbit='red', beta.MADX='blue',
                                       orbit='green', beta.linear='purple')) +
        labs(x = 'Beams separation, um',
             y = 'L(with bb) / L(wo bb)')
    pggp <- ggplot_build(ggp)
  }
}

gg.vdm.profiles <- function() {
    . <- old.lumi.w.over.wo[, {
        list(separation = separation,
             old = vdm.profile * beta.MADX.orbit,
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
gg.vdm.profiles <- function() {
    . <- old.lumi.w.over.wo[, {
        list(separation = separation,
             old = vdm.profile * beta.MADX.orbit,
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
                                                  y = vdm.profile * (beta.MADX.orbit - 1),
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
dirName <- gsub('^.*/','',dir)
cat(' dirName= ',dirName, '\n')
i.plot <- 0
save.plot <- function(plot.function.name) {
#    pdf(file.path(plot.dir, paste0(plot.function.name, '.pdf')), width=12, height=8)
    pdf(file.path(plot.dir, paste0(plot.function.name, '_', dirName, '.pdf')), width=12, height=8)
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



#save.plot(paste('gg.cor_',dirName,sep=''))
#save.plot(paste('gg.vdm.profiles_',dirName,sep=''))
#save.plot(paste('gg.vdm.profile.change_',dirName,sep=''))
#save.plot(paste('gg.old.from.new_',dirName,sep=''))
#save.plot(paste('gg.integ.per.turn_',dirName,sep=''))
#save.plot(paste('gg.integ.per.100.turns_',dirName,sep=''))
#save.plot(paste('gg.correction_',dirName,sep=''))
#save.plot(paste('gg.center_',dirName,sep=''))

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
print('print old')
print(old)
print('new.lumi.w.over.wo')
print(new.lumi.w.over.wo)
print('old.lumi.w.over.wo')
print(old.lumi.w.over.wo)
print('gg.cor')
print(gg.cor)
print('old.lumi.w.over.wo.alt.range')
print(old.lumi.w.over.wo.alt.range)
