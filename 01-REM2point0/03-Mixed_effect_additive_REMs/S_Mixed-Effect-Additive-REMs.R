## ----setup, include=FALSE-----------------------------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
if (!require("mgcv", quietly = TRUE)) {
  # If not installed, install it
  install.packages("mgcv")
  # Load the package
  library("mgcv")
} else {
  if (!require("splines", quietly = TRUE)) {
    install.packages("splines")
    library("splines")
  } else {
    if (!require("ggplot2", quietly = TRUE)) {
      install.packages("ggplot2")
      library("ggplot2")
    } else {
      if (!require("tidyverse", quietly = TRUE)) {
        install.packages("tidyverse")
        library("tidyverse")
      } else {
        if (!require("RColorBrewer", quietly = TRUE)){
          install.packages("RColorBrewer")
        } else {
          if (!require("mgcViz", quietly = TRUE)){
            install.packages("mgcViz")
          } else {
            library("mgcv")
            library("splines")
            library("ggplot2")
            library("tidyverse")
            library("RColorBrewer")
            library("mgcViz")
          }
        }
      }
    }
  }
}

## -----------------------------------------------------------------------------
pal.blue <- brewer.pal(9, "Blues")
pal.rose <- brewer.pal(9, "RdPu")
colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#999999", "#66C2A5", "#FC8D62")

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("02-Images/insects.pdf")

## -----------------------------------------------------------------------------
load(file="01-Data/01-Inputs/input00.RData")

## -----------------------------------------------------------------------------
head(FR[,c("LifeForm", "Taxon",
           "Region", "FirstRecord",
           "Source")])

## -----------------------------------------------------------------------------
native[sample(1:nrow(native), 10),c("species", "region")]

## -----------------------------------------------------------------------------
head(first_records[,c("year","lf",
           "species", 
           "region")])

## -----------------------------------------------------------------------------
invaded.regions <- function(sp.n, r.n, y, native, first_records){
  
  # Convert input arguments to numeric type if not already
  sp.n <- as.numeric(sp.n)
  r.n <- as.numeric(r.n)
  y <- as.numeric(y)
  
  # Get unique combinations of species number and region number from native data
  t <- unique(as.vector(subset(native, sp.num == sp.n, r.num)))
  
  # Get region numbers from first records data where species number matches
  pr <- as.vector(subset(first_records, sp.num == sp.n, r.num))
  
  # If the invasion is both present in first records data and in native range
  # consider the former as actual piece of information
  t <- setdiff(t, pr)
  
  # Find indices of first records occurring before end date for the species
  set.sp <- which(first_records$sp.num == sp.n & first_records$year < y)
  
  # Combine regions in native range with regions in first records before date
  t <- na.omit(c(t, first_records$r.num[set.sp]))
  
  # Do not consider the involved region
  inv <- unlist(setdiff(t, r.n))
  
  # Return invaded regions
  return(inv)
}

## -----------------------------------------------------------------------------
# Define sender set - species
spec <- unique(first_records$species)
(s <- length(spec))
# Define receiver set - regions
reg.lf <- unique(first_records$region)
(r <- length(reg.lf))

## -----------------------------------------------------------------------------
reg <- colnames(data_distance)

## -----------------------------------------------------------------------------
creating_case_control_dataset <- function(first_records,
                                 spec, reg.lf, reg,
                                 seed=1234
                                 ){
  
  reg.lf.num <- match(reg.lf, reg)
  s <- length(spec)
  r <- length(reg.lf)
  
  ## POSSIBLE (S,R) INTERACTIONS ####
  # Initialize a vector to keep track of the risk set size over time
  at.risk <- NULL
  # Create a matrix to indicate possible (species, region) interactions
  alien.occ <- matrix(0, nrow = s, ncol = r) 
  rownames(alien.occ) <- spec
  colnames(alien.occ) <- reg.lf
  for (n.sp in 1:s){
    # Identify native regions for each species
    nat.id <- unique(native$r.num[native$sp.num == n.sp]) 
    # Identify regions where the species is not native
    possible.to <- setdiff(reg.lf.num,nat.id)
    # Mark these regions as possible invasion sites for the species
    alien.occ[n.sp,reg[possible.to]] <- 1
  }
  
  ## COLLECTING INFORMATION ####
  dat.gam <- data.frame(matrix(NA, nrow=nrow(first_records),ncol=6))
  colnames(dat.gam) <- c("y", "year", 
                         "sp1", "r1", 
                         "sp2", "r2")
  # The response is fixed and equal to 1
  dat.gam[,1] <- rep(1, nrow(first_records))
  
  set.seed(seed)
  # For each FR:
  for (i in 1:nrow(first_records)){
    
    ### INFORMATION CONCERNING THE EVENT ####
    # year of the invasion event
    dat.gam[i,2] <- year <- first_records[i,"year"]
    # invading species
    dat.gam[i,3] <- s.ev <- first_records[i,"species"]
    # invaded country
    dat.gam[i,4] <- r.ev <- first_records[i,"region"]
    
    ### POSSIBLE EVENTS ####
    # Events occurred at the same time of the considered event
    # are removed from the risk set
    sub_stp <- first_records[first_records$year==year,
                             c("species", "region")]
    ni <- nrow(sub_stp)
    for (j in 1:ni){
      # Mark these (species, region) pairs as not at risk
      alien.occ[sub_stp[j,1],sub_stp[j,2]] <- 0
    }
    at.risk <- c(at.risk, sum(alien.occ==1))
    
    ### SAMPLING THE NON-EVENT ####
    sr.nv<-sample(which(alien.occ!=0),1)
    # species non-event
    dat.gam[i,5] <- s.nv <- spec[(sr.nv-1)%%s+1]
    # region non-event
    dat.gam[i,6] <- r.nv <- reg.lf[(sr.nv-1)%/%s+1]
  }
  
  return(dat.gam)
}

## -----------------------------------------------------------------------------
dat.gam <- creating_case_control_dataset(first_records, 
                                         spec, reg.lf, reg)

## -----------------------------------------------------------------------------
dat.gam$sp1.num <- match(dat.gam$sp1, spec)
dat.gam$r1.num <- match(dat.gam$r1, reg)
dat.gam$sp2.num <- match(dat.gam$sp2, spec)
dat.gam$r2.num <- match(dat.gam$r2, reg)

## -----------------------------------------------------------------------------
head(dat.gam[c("year","sp1", "r1","sp2", "r2")])
head(dat.gam[c("year","sp1.num", "r1.num","sp2.num", "r2.num")])

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("02-Images/trade-idea.pdf")

## -----------------------------------------------------------------------------
load("01-Data/01-Inputs/input01.RData")

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("02-Images/trade-covariate-computation.pdf")

## -----------------------------------------------------------------------------
t <- which(data_trade$transfer < 0)
data_trade$transfer[t] <- 0

## -----------------------------------------------------------------------------
trade.funct <- function(inv, r.n, y, reg, data_trade){
  
  # Convert input arguments to numeric type if not already
  r.n <- as.numeric(r.n)
  y <- as.numeric(y)
  
  # Check if there are already invaded countries
  if(length(inv)!=0){
    
    # Find rows that involve invaded countries as sending trade
    u <- which(data_trade$FromRegion %in% reg[inv])
    # Find rows that involve region of interest as receiving trade
    v <- which(data_trade$ToRegion == reg[r.n])
    # Find the intersection of the two sets
    w <- intersect(u,v)
    x <- data_trade[w,]
    # Consider the trade instances occurred before or at the time of interest
    x <- x[x$year<=y,]
    trade_value <- NULL
    # If there are rows in the filtered dataset
    if(nrow(x)>0){
      # For each invaded country, the maximum year is recorded
      o <- aggregate(x$year, list(x$FromRegion), FUN=max)
      # For each of them, the corresponding transfer is stored
      for (o.i in 1:nrow(o)){
        trade_value <- c(trade_value, 
                         x$transfer[x$FromRegion==o[o.i,1] & 
                                      x$year==o[o.i,2]])}
    }
  } else { 
  # If there are not already invaded countries, trade is set equal to 0
    trade_value <- 0
  }
  # Compute the log-transformed sum of trade values (with an added constant 1)
  log_trade.value <- ifelse(length(trade_value)>0, 
                            log(sum(trade_value, na.rm =T)+1),0)
  
  # Return the computed log-transformed trade value
  return(log_trade.value)
}

## -----------------------------------------------------------------------------
log_trade <- function(sp.n, r.n, y, native, first_records, reg, data_trade){
  
  inv <- invaded.regions(sp.n = sp.n, 
                         r.n = r.n, 
                         y = y, 
                         native = native, 
                         first_records = first_records)
  
  log_trade.value <- ifelse(r.n==match("USACanada", reg), 
                            mean(trade.funct(inv = inv,
                                             r.n = match("United States",reg),
                                             y = y, 
                                             reg = reg,
                                             data_trade = data_trade),
                                 trade.funct(inv = inv,
                                             r.n = match("Canada",reg),
                                             y = y, 
                                             reg = reg,
                                             data_trade = data_trade)),
                            trade.funct(inv = inv,
                                        r.n = r.n,
                                        y = y,
                                        reg = reg, 
                                        data_trade = data_trade))
  
  return(log_trade.value)
}

## -----------------------------------------------------------------------------
dat.gam$tr1 <- apply(dat.gam[,c("sp1.num", "r1.num", "year")], 1, 
                            function(x) log_trade(sp.n = x[1], 
                                                  r.n = x[2], 
                                                  y = x[3], 
                                                  native = native, 
                                                  first_records = 
                                                    first_records,
                                                  reg = reg,
                                                  data_trade = data_trade))
dat.gam$tr2 <- apply(dat.gam[,c("sp2.num", "r2.num", "year")], 1, 
                            function(x) log_trade(x[1], x[2], x[3], 
                                                  native = native, 
                                                     first_records = 
                                                       first_records,
                                                  reg = reg,
                                                  data_trade = data_trade))
dat.gam$tr = dat.gam$tr1 - dat.gam$tr2

## -----------------------------------------------------------------------------
q = 10

## -----------------------------------------------------------------------------
bspline <- function(x, k, i, m = 2) {
  # ith B-spline basis function of order m at the values in x
  # given knot locations in k
  
  if (m == -1) {
    # Base case of the recursion: 
    # when m is -1, we are at the lowest order basis function
    res <- as.numeric(x < k[i + 1] & x >= k[i])  
    # Returns 1 if x is within the interval [k[i], k[i+1])
    
  } else {
    
    # Recursive case: 
    # B-spline basis function from lower order basis functions
    # Calculate the first term's coefficient
    z0 <- (x - k[i]) / (k[i + m + 1] - k[i])
    # Calculate the second term's coefficient
    z1 <- (k[i + m + 2] - x) / (k[i + m + 2] - k[i + 1])
    
    # Recursive calls to the lower order basis functions
    res <- z0 * bspline(x, k, i, m - 1) + z1 * bspline(x, k, i + 1, m - 1)
  }
  return(res)  # Return the evaluated B-spline basis function
}

## -----------------------------------------------------------------------------
# Equally spaced knots from 1880 to 2005
n_knots = 14.
knots = seq(from = min(dat.gam$year),
            to = max(dat.gam$year),
            length.out = n_knots)

## -----------------------------------------------------------------------------
m = 2
basis = matrix(0, nrow = length(dat.gam$year), ncol = q)
for (j in 1:q) {
  basis[, j] = bspline(dat.gam$year, k = knots, i = j)
}

plot(y = basis[, 1],
     x = dat.gam$year,
     ylab = 'b(t) - basis function of time', xlab = 't - time',
     col = 0,
     cex = 0.2)
for (j in 1:10) {
  points(y = basis[, j],
         x = dat.gam$year,
         cex = 0.2,
         col = colors[j])
}
for (k in knots) {
  abline(v=k, lty=2, lwd=0.5)
}

## -----------------------------------------------------------------------------
x.ev <- dat.gam$tr1
x.nv <- dat.gam$tr2
x <- x.ev - x.nv

## -----------------------------------------------------------------------------
stp <- dat.gam$year

## -----------------------------------------------------------------------------
gam_tr.only <- gam(y ~ s(stp, by=x) - 1,
    family="binomial"(link = 'logit'),
    method="REML", data=dat.gam)

## -----------------------------------------------------------------------------
plot(gam_tr.only)

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("02-Images/distance-idea.pdf")

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("02-Images/distance-covariate-computation.pdf")

## -----------------------------------------------------------------------------
load(file="01-Data/02-Gam-Fits/gam_d.only.RData")

## -----------------------------------------------------------------------------
plot(gam_d.only)

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("02-Images/insect-idea.pdf")

## -----------------------------------------------------------------------------
load("01-Data/01-Inputs/input03.RData")

## -----------------------------------------------------------------------------
sp1 <- dat.gam$sp1
sp2 <- dat.gam$sp2
sp <- factor(c(sp1,sp2))
dim(sp) <- c(length(sp1),2)

## -----------------------------------------------------------------------------
load(file="01-Data/02-Gam-Fits/gam_sp.only.RData")

## -----------------------------------------------------------------------------
re.species <- coefficients(gam_sp.only)
names(re.species) <- levels(sp)

## -----------------------------------------------------------------------------
sort(re.species, decreasing = TRUE)[1:5]

## -----------------------------------------------------------------------------
sort(re.species)[1:5]

## -----------------------------------------------------------------------------
load("01-Data/01-Inputs/input04.RData")

## -----------------------------------------------------------------------------
unit <- rep(1, nrow(dat.gam))
stp = dat.gam$year
X = cbind(dat.gam$d1,dat.gam$d2)
I = cbind(unit,-unit)		

## -----------------------------------------------------------------------------
gam_complete <- gam(y ~ dt + 
                      s(stp, by=tr) +
                      s(X, by=I) +
                      s(sp, by=I, bs="re") - 1,
    family="binomial"(link = 'logit'),
    method="REML", data=dat.gam)

## -----------------------------------------------------------------------------
load(file="01-Data/02-Gam-Fits/gam_dt.only.RData")
load(file="01-Data/02-Gam-Fits/gam_tr.only.RData")
load(file="01-Data/02-Gam-Fits/gam_d.only.RData")
load(file="01-Data/02-Gam-Fits/gam_sp.only.RData")

## -----------------------------------------------------------------------------
AIC(gam_dt.only)

## -----------------------------------------------------------------------------
AIC(gam_tr.only)

## -----------------------------------------------------------------------------
AIC(gam_d.only)

## -----------------------------------------------------------------------------
AIC(gam_sp.only)

## -----------------------------------------------------------------------------
AIC(gam_complete)

## -----------------------------------------------------------------------------
sort(re.species, decreasing = TRUE)[1:5]
sort(re.species)[1:5]

## -----------------------------------------------------------------------------
re.species_complete <- coefficients(gam_complete)[21:length(coefficients(gam_complete))]
names(re.species_complete) <- levels(sp)
sort(re.species_complete, decreasing = TRUE)[1:5]
sort(re.species_complete)[1:5]

## -----------------------------------------------------------------------------
sort(re.species)[1:3]
sort(re.species_complete)[1:3]

