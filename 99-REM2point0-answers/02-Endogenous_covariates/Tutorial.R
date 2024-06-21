## ----setup, include=FALSE---------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE-------------------------
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

## ---------------------------------------------------------
pal.blue <- brewer.pal(9, "Blues")
pal.rose <- brewer.pal(9, "RdPu")
colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7", "#999999", "#66C2A5", "#FC8D62")

## ----echo=FALSE-------------------------------------------
knitr::include_graphics("02-Images/class-network.pdf")

## ---------------------------------------------------------
load("01-Data/01-Inputs/input01.RData")

## ---------------------------------------------------------
head(dat.gam)

## ---------------------------------------------------------
receiver.teacher <- class$info$teacher[dat.gam$receiver] 
non.receiver.teacher <- class$info$teacher[dat.gam$non.receiver]
receiver.teacher.diff <- receiver.teacher - non.receiver.teacher
dat.gam <- cbind(dat.gam,
                 receiver.teacher,non.receiver.teacher,
                 receiver.teacher.diff, y=1)

## ---------------------------------------------------------
summary(gam(y ~ receiver.teacher.diff - 1, 
            data=dat.gam, 
            family="binomial"))

## ---------------------------------------------------------
head(dat.gam)

## ---------------------------------------------------------
receiver_popularity <- function(dat.gam, k, ev=TRUE) {
  
  # Check if k is within bounds
  if (k < 1 || k > nrow(dat.gam)) {
    stop("k is out of bounds.")
  }
  
  # If k is 1, there's no history to compare, return 0
  if (k == 1) {
    return(0)
  } else {
    
    # Calculate the history of receivers up to the (k-1)th element
    history_receiver <- dat.gam$receiver[1:(k-1)]
    
    # Ensure that dat.gam has the required columns
    if (!"receiver" %in% names(dat.gam)) {
      stop("dat.gam does not contain the column 'receiver'.")
    }
    if (!"non.receiver" %in% names(dat.gam) && !ev) {
      stop("dat.gam does not contain the column 'non.receiver'.")
    }
    
    # Calculate popularity based on the value of ev
    if (ev) {
      return(sum(history_receiver == dat.gam$receiver[k]))
    } else {
      return(sum(history_receiver == dat.gam$non.receiver[k]))
    }
  }
}


## ---------------------------------------------------------
dat.gam$endPop1 <- sapply(1:nrow(dat.gam), 
                          function(x) receiver_popularity(dat.gam, x, ev=TRUE))
dat.gam$endPop2 <- sapply(1:nrow(dat.gam), 
                          function(x) receiver_popularity(dat.gam, x, ev=FALSE))
endPop_diff = dat.gam$endPop1 - dat.gam$endPop2

## ---------------------------------------------------------
summary(gam(y ~ endPop_diff - 1, 
            data=dat.gam, 
            family="binomial"))

## ---------------------------------------------------------
rec_cyc <- function(dat.gam, k, ev=TRUE) {
  
  # Determine which column indices to use based on the 'ev' flag
  if (ev) {
    s <- dat.gam[k, 2]   # 's' is the second column (event sender)
    r <- dat.gam[k, 3]   # 'r' is the third column (event receiver)
  } else {
    s <- dat.gam[k, 4]   # 's' is the fourth column (event sender)
    r <- dat.gam[k, 5]   # 'r' is the fifth column (event receiver)
  }
  
  # Check if k is within bounds of dat.gam
  if (k < 1 || k > nrow(dat.gam)) {
    stop("k is out of bounds.")
  }
  
  # If k is 1, there's no history to compare, return (0, 0)
  if (k == 1) {
    return(c(0, 0))
  } else {
    # Initialize variables to store reciprocity and cyclic closure indicators
    rec_id <- 0
    cyc_id <- 0
    
    # Determine the starting point in the past to begin examining
    prev <- which(dat.gam[1:(k-1), 2] == s & dat.gam[1:(k-1), 3] == r)
    if (length(prev) > 0) {
      prev <- prev[length(prev)] + 1
    } else {
      prev <- 1
    }
    
    # Check if there are third parties interacting towards 's'
    t <- which(dat.gam[prev:k, 3] == s)
    t <- prev + t - 1
    l <- dat.gam[t, 2]
    common.nodes <- unique(l)
    
    if (length(common.nodes) != 0) {
      # Iterate over each common node 'a'
      for (a in common.nodes) {
        
        # Check if 'a' is the receiver 'r' (reciprocity check)
        if (a == r) {
          rec_id <- max(rec_id, 1)
        } else {
          # Find the latest time 'a' interacted towards 's'
          stop_t <- dat.gam$time[t[which(dat.gam[t, 2] == a & 
                                         dat.gam[t, 3] == s)]]
          stop_t <- max(stop_t)
          
          # Check if 'r' interacted towards 'a' before 'stop_t'
          t2 <- which(dat.gam[prev:max(which(dat.gam$time == stop_t)), 2] == r &
                      dat.gam[prev:max(which(dat.gam$time == stop_t)), 3] == a)
          
          # If at least one such event is found, set cyclic closure indicator
          if (length(t2) > 0) {
            cyc_id <- max(cyc_id, 1)
          } else {
            # Check if 'r' interacted towards 'a' before 's' interacted towards 'r'
            t3 <- which(dat.gam[1:prev, 2] == r & dat.gam[1:prev, 3] == a)
            if (length(t3) > 0) {
              cyc_id <- max(cyc_id, 1)
            }
          }
        }
      }
    }
  }
  
  # Return the reciprocity and cyclic closure indicators
  return(c(rec_id, cyc_id))
}

## ---------------------------------------------------------
events <- t(sapply(1:nrow(dat.gam), 
                 function(x) rec_cyc(dat.gam, x, ev=TRUE)))
non_events <- t(sapply(1:nrow(dat.gam), 
                function(x) rec_cyc(dat.gam, x, ev=FALSE)))

dat.gam <- cbind(dat.gam, 
                 rec1 = events[,1],
                 cyc1 = events[,2],
                 rec2 = non_events[,1],
                 cyc2 = non_events[,2])
dat.gam$rec_diff = dat.gam$rec1 - dat.gam$rec2
dat.gam$cyc_diff = dat.gam$cyc1 - dat.gam$cyc2

## ---------------------------------------------------------
events_alt <- n_effects(dat.gam[,c(2,3,1)], dat.gam[,c(4,5,1)])[[1]]
non_events_alt <- n_effects(dat.gam[,c(2,3,1)], dat.gam[,c(4,5,1)])[[2]]

## ---------------------------------------------------------
# events_alt$r1==dat.gam$rec1

## ---------------------------------------------------------
# non_events_alt$r1==dat.gam$rec2

## ---------------------------------------------------------
summary(gam(y ~ endPop_diff + rec_diff - 1, data=dat.gam, 
            family="binomial"))

## ---------------------------------------------------------
summary(gam(y ~ endPop_diff + rec_diff + cyc_diff - 1, data=dat.gam, 
            family="binomial"))

## ---------------------------------------------------------
transitive_closure_slt <- function(dat.gam, k, alpha, ev=TRUE) {
  
  current <- as.numeric(dat.gam[k, 1])
  
  # Determine which column indices to use based on the 'ev' flag
  if (ev) {
    s <- dat.gam[k, 2]   # 's' is the second column (event sender)
    r <- dat.gam[k, 3]   # 'r' is the third column (event receiver)
  } else {
    s <- dat.gam[k, 4]   # 's' is the fourth column (event sender)
    r <- dat.gam[k, 5]   # 'r' is the fifth column (event receiver)
  }
  
  # Check if k is within bounds of dat.gam
  if (k < 1 || k > nrow(dat.gam)) {
    stop("k is out of bounds.")
  }
  
  # If k is 1, there's no history to compare, return (0, 0)
  if (k == 1) {
    return(c(0, 0, 0))
  } else {
    
    # Determine the starting point in the past to begin examining
    prev <- which(dat.gam[1:(k-1), 2] == s & dat.gam[1:(k-1), 3] == r)
    if (length(prev) > 0) {
      prev <- prev[length(prev)] + 1
    } else {
      prev <- 1
    }
    
    trs_rec <- trs_long <- trs_id <- 0
    
    # Check if there are third parties interacting towards 'r'
    t <- which(dat.gam[prev:k,3]==r)
    t <- prev + t - 1
    l <- dat.gam[t,2]
    common.nodes <- unique(l)
    
    if(length(common.nodes)!=0){
      
      # Iterate over each common node 'a'
      for(a in common.nodes){
      
      # Find the latest time 'a' interacted towards 'r'
      stop_t <- dat.gam$time[t[which(dat.gam[t,2]==a & dat.gam[t,3] == r)]]
      stop_t <- max(stop_t)
      
      # Check if 's' interacted towards 'a' before 'stop_t'
      t2 <- which(dat.gam[prev:max(which(dat.gam$time==stop_t)),2]==s &
                  dat.gam[prev:max(which(dat.gam$time==stop_t)),3]==a)
      
      # If at least one such event is found, 
      # check if latest time 'a' interacted towards 'r' 
      # is larger or smaller than (current-alpha)
      if (length(t2)>0){
        if (stop_t>=(current-alpha)){
          trs_rec <- max(trs_rec,1)
        } else {
          trs_long <- max(trs_long, 1)
        }
        trs_id <- max(trs_id, 1)
      } else {
        
        # Check if 's' interacted towards 'a' before 's' interacted towards 'r'
        t3 <- which(dat.gam[1:prev,2]==s & dat.gam[1:prev,3]==a)
        if (length(t3) > 0) {
          if (stop_t>=(current-alpha)){
            trs_rec <- max(trs_rec,1)
            } else {
              trs_long <- max(trs_long, 1)
              }
          trs_id <- max(trs_id, 1)
        }
      }
      }
    }
  }
  # Return the reciprocity and cyclic closure indicators
  return(c(trs_rec, trs_long, trs_id))
}

## ----echo=FALSE-------------------------------------------
knitr::include_graphics("02-Images/timeline.pdf")

## ---------------------------------------------------------
# rem_dat <- matrix(c(1, "A", "B", 2, "A", "B",
#                     3, "A", "B", 4, "B", "A",
#                     5, "A", "B", 6, "B", "C",
#                     7, "C", "B", 8, "B", "A",
#                     9, "A", "C", 10, "C", "A",
#                     11, "C", "B", 12, "A", "B",
#                     13, "A", "B", 14, "A", "B",
#                     15, "B", "C", 16, "A", "B",
#                     17, "A", "C", 18, "C", "B",
#                     19, "B", "C", 20, "A", "B"),
#                   ncol=3, byrow=TRUE)
# colnames(rem_dat) <- c("time", "sender", "receiver")
# rem_dat <- data.frame(rem_dat)
# rem_dat$time <- as.numeric(rem_dat$time)

## ---------------------------------------------------------
# save.image(file="input02.RData")

## ---------------------------------------------------------
load(file="01-Data/01-Inputs/input02.RData")

## ---------------------------------------------------------
rec_cyc_slt <- function(dat.gam, k, alpha, ev=TRUE) {
  
  current <- as.numeric(dat.gam[k, 1])
  
  # Determine which column indices to use based on the 'ev' flag
  if (ev) {
    s <- dat.gam[k, 2]   # 's' is the second column (event sender)
    r <- dat.gam[k, 3]   # 'r' is the third column (event receiver)
  } else {
    s <- dat.gam[k, 4]   # 's' is the fourth column (event sender)
    r <- dat.gam[k, 5]   # 'r' is the fifth column (event receiver)
  }
  
  # Check if k is within bounds of dat.gam
  if (k < 1 || k > nrow(dat.gam)) {
    stop("k is out of bounds.")
  }
  
  # If k is 1, there's no history to compare, return (0, 0)
  if (k == 1) {
    return(c(0, 0, 0, 0, 0, 0))
  } else {
    # Initialize variables to store reciprocity and cyclic closure indicators
    rec_rec <- rec_long <- rec_id <- 0
    cyc_rec <- cyc_long <- cyc_id <- 0
    
    # Determine the starting point in the past to begin examining
    prev <- which(dat.gam[1:(k-1), 2] == s & dat.gam[1:(k-1), 3] == r)
    if (length(prev) > 0) {
      prev <- prev[length(prev)] + 1
    } else {
      prev <- 1
    }
    
    # Check if there are third parties interacting towards 's'
    t <- which(dat.gam[prev:k, 3] == s)
    t <- prev + t - 1
    l <- dat.gam[t, 2]
    common.nodes <- unique(l)
    
    if (length(common.nodes) != 0) {
      # Iterate over each common node 'a'
      for (a in common.nodes) {
        
        # Check if 'a' is the receiver 'r' (reciprocity check)
        if (a == r) {
          stop_t <- dat.gam$time[t[which(dat.gam[t, 2] == r & 
                                         dat.gam[t, 3] == s)]]
          stop_t <- max(stop_t)
          if (stop_t>=(current-alpha)){
            rec_rec <- max(rec_rec,1)
            } else {
            rec_long <- max(rec_long, 1)
          }
          rec_id <- max(rec_id, 1)
        } else {
          # Find the latest time 'a' interacted towards 's'
          stop_t <- dat.gam$time[t[which(dat.gam[t, 2] == a & 
                                         dat.gam[t, 3] == s)]]
          stop_t <- max(stop_t)
          
          # Check if 'r' interacted towards 'a' before 'stop_t'
          t2 <- which(dat.gam[prev:max(which(dat.gam$time == stop_t)), 2] == r &
                      dat.gam[prev:max(which(dat.gam$time == stop_t)), 3] == a)
          
          # If at least one such event is found, set cyclic closure indicator
          if (length(t2) > 0) {
            if (stop_t>=(current-alpha)){
              cyc_rec <- max(cyc_rec,1)
            } else {
              cyc_long <- max(cyc_long, 1)
            }
            cyc_id <- max(cyc_id, 1)
          } else {
            # Check if 'r' interacted towards 'a' before 's' interacted towards 'r'
            t3 <- which(dat.gam[1:prev, 2] == r & dat.gam[1:prev, 3] == a)
            if (length(t3) > 0) {
              if (stop_t>=(current-alpha)){
                cyc_rec <- max(cyc_rec,1)
              } else {
                cyc_long <- max(cyc_long, 1)
              }
              cyc_id <- max(cyc_id, 1)
            }
          }
        }
      }
    }
  }
  
  # Return the reciprocity and cyclic closure indicators
  return_object <- c(rec_rec, rec_long, rec_id, cyc_rec, cyc_long, cyc_id)
  names(return_object) <- c("rec_short", 
                            "rec_long",
                            "rec_id",
                            "cyc_short",
                            "cyc_long",
                            "cyc_id")
  return(return_object)
}

## ---------------------------------------------------------
knitr::include_graphics("02-Images/dynamics-simulated.pdf")

## ---------------------------------------------------------
frame_1 <- t(sapply(1:20, function(x) rec_cyc_slt(dat.gam = rem_dat, k=x,
                                      alpha = 3)))
colnames(frame_1) <- c("rec_short", 
                       "rec_long",
                       "rec_id",
                       "cyc_short",
                       "cyc_long",
                       "cyc_id")
frame_1

## ---------------------------------------------------------
frame_2 <- t(sapply(1:20, function(x) transitive_closure_slt(dat.gam = rem_dat, k=x,
                                      alpha = 3)))
colnames(frame_2) <- c("trs_short", "trs_long", "trs_id")
frame_2

## ---------------------------------------------------------
repetition_slt <- function(dat.gam, k, alpha, ev=TRUE, current=NULL) {
  
  if(is.null(current)){
    current <- as.numeric(dat.gam[k, 1])
  }
  
  # Determine which column indices to use based on the 'ev' flag
  if (ev) {
    s <- dat.gam[k, 2]   # 's' is the second column (event sender)
    r <- dat.gam[k, 3]   # 'r' is the third column (event receiver)
  } else {
    s <- dat.gam[k, 4]   # 's' is the fourth column (event sender)
    r <- dat.gam[k, 5]   # 'r' is the fifth column (event receiver)
  }
  
  # Check if k is within bounds of dat.gam
  if (k < 1 || k > nrow(dat.gam)) {
    stop("k is out of bounds.")
  }
  
  # If k is 1, there's no history to compare, return 0
  if (k == 1) {
    return(c(0, 0, 0))
  } else {
    # Initialize variables to store reciprocity and cyclic closure indicators
    rep_rec <- rep_long <- rep_id <- 0
    
    # Determine the starting point in the past to begin examining
    prev <- which(dat.gam[1:(k-1), 2] == s & dat.gam[1:(k-1), 3] == r)
    if (length(prev) > 0) {
      stop_t <- max(dat.gam[prev,"time"])
      if (stop_t>=(current-alpha)){
        rep_rec <- max(rep_rec,1)
      } else {
        rep_long <- max(rep_long, 1)
        }
      rep_id <- 1
    }
  }
  return(c(rep_rec, rep_long, rep_id))
}

## ---------------------------------------------------------
frame_3 <- t(sapply(1:20, function(x) repetition_slt(dat.gam = rem_dat, k=x,
                                      alpha = 3)))
colnames(frame_3) <- c("rep_short", "rep_long", "rep_id")
frame_3

## ----echo=FALSE-------------------------------------------
knitr::include_graphics("02-Images/exo-endo.pdf")

## ---------------------------------------------------------
load("01-Data/01-Inputs/input03.RData")

## ---------------------------------------------------------
M = gam(one ~ ## OUT-DEGREE + OUT-ACTIVITY + IN-DEGREE + IN-ACTIVITY
          RecOutDeg + RecNROutAct + RecInDeg + RecNRInAct +
          ## NODE BALANCE
          RecNetBal +
          ## TIME FROM LAST EVENT
          RecLastTime +
          ## ASSORTATIVITY BY DEGREE AND BY INTENSITY
          SdrOutDeg_RecInDeg + SdrNROutAct_RecNRInAct +
          ## REPETITION AND RECIPROCITY
          EdTrt + EdRec +
          ## DYADIC BALANCE
          EdBalDif + 
          ## TRANSITIVE CLOSURE + CYCLIC CLOSURE + SHARED LENDERS + SHARED BORROWERS
          EdTran + EdClos + EdShrBors + EdShrLens - 1,

        data=dat.gam, 
        family="binomial")

## ---------------------------------------------------------
# save(M, file="01-Data/02-Gam-Fits/baseline-model.RData")
summary(M)

## ---------------------------------------------------------
M_SLT = gam(one ~  
              ## OUT-DEGREE + OUT-ACTIVITY
              RecOutDeg + ShortRecOutDeg + 
              RecNROutAct+ ShortRecNROutAct + 
              ## IN-DEGREE + IN-ACTIVITY
              RecInDeg + ShortRecInDeg + 
              RecNRInAct+ ShortRecNRInAct + 
              ## NODE BALANCE
              RecNetBal +
              ## NODE TIME FROM LAST EVENT
              RecLastTime +
              ## ASSORTATIVITY BY DEGREE AND BY INTENSITY
              SdrOutDeg_RecInDeg + SdrShortOutDeg_ShortRecInDeg + 
              SdrNROutAct_RecNRInAct + SdrShortNROutAct_ShortRecNRInAct +
              ## REPETITION
              EdTrt + ShortEdTrt +
              ## RECIPROCITY
              EdRec + ShortEdRec + 
              ## DYADIC BALANCE
              EdBalDif + 
              ## TRANSITIVE CLOSURE
              EdTran + ShortEdTran + 
              ## CYCLIC CLOSURE
              EdClos + ShortEdClos +
              # ## SHARED BORROWERS
              EdShrBors + ShortEdShrBors + 
              ## SHARED LENDERS
              EdShrLens + ShortEdShrLens -1,

            family = "binomial", 
            data = dat.gam)

## ---------------------------------------------------------
# save(M_SLT, file="01-Data/02-Gam-Fits/short-long-term-model.RData")
summary(M_SLT)

## ---------------------------------------------------------
M_TW = gam(one ~  
             ## OUT-DEGREE + OUT-ACTIVITY + IN-DEGREE + IN-ACTIVITY
             RecOutDeg + RecNROutAct + RecInDeg + RecNRInAct +
             ## NODE BALANCE
             RecNetBal +
             ## TIME FROM LAST EVENT
             RecLastTime +
             ## ASSORTATIVITY BY DEGREE AND BY INTENSITY
             SdrOutDeg_RecInDeg + SdrNROutAct_RecNRInAct +
             ## REPETITION AND RECIPROCITY
             EdTWTrt + EdTWRec +
             ## DYADIC BALANCE
             EdBalDif + 
             ## TRANSITIVE CLOSURE + CYCLIC CLOSURE + SHARED BORROWERS + SHARED LENDERS
             EdTWTran + EdTWClos + EdTWShrBors + EdShrLens -1,

           family = "binomial", 
           data = dat.gam)

## ---------------------------------------------------------
# save(M_TW, file="01-Data/02-Gam-Fits/time-weighted-model.RData")
summary(M_TW)

## ---------------------------------------------------------
M_EW1 = gam(one ~  
              ## OUT-DEGREE + OUT-ACTIVITY + IN-DEGREE + IN-ACTIVITY
              RecOutDeg + RecNROutAct + RecInDeg + RecNRInAct +
              ## NODE BALANCE
              RecNetBal +
              ## TIME FROM LAST EVENT
              RecLastTime +
              ## ASSORTATIVITY BY DEGREE AND BY INTENSITY
              SdrOutDeg_RecInDeg + SdrNROutAct_RecNRInAct +
              ## REPETITION AND RECIPROCITY
              # EdTrt + EdTrtEvtWei + EdRec + EdRecEvtWei +
              EdTrt + EdWeiTrt + EdRec + EdWeiRec +
              ## DYADIC BALANCE
              EdBalDif + 
              ## TRANSITIVE CLOSURE + CYCLIC CLOSURE + SHARED BORROWERS + SHARED LENDERS
              EdTran + EdClos + EdShrBors + EdShrLens -1,

            family = "binomial", 
            data = dat.gam)

summary(M_EW1)

## ---------------------------------------------------------
M_EW2 = gam(one ~
              ## OUT-DEGREE + OUT-ACTIVITY + IN-DEGREE + IN-ACTIVITY
              RecOutDeg + RecNROutAct + RecInDeg + RecNRInAct +
              ## NODE BALANCE
              RecNetBal +
              ## TIME FROM LAST EVENT
              RecLastTime +
              ## ASSORTATIVITY BY DEGREE AND BY INTENSITY
              SdrOutDeg_RecInDeg + SdrNROutAct_RecNRInAct +
              ## REPETITION AND RECIPROCITY
              EdTrt + EdNORWeiTrt + EdRec + EdNORWeiRec + 
              ## DYADIC BALANCE
              EdBalDif + 
              ## TRANSITIVE CLOSURE + CYCLIC CLOSURE + SHARED BORROWERS + SHARED LENDERS
              EdTran + EdClos + EdShrBors + EdShrLens -1,              
            family = "binomial", 
            data = dat.gam)

summary(M_EW2)

## ---------------------------------------------------------
M_EW3 = gam(one ~
              ## OUT-DEGREE + OUT-ACTIVITY + IN-DEGREE + IN-ACTIVITY
              RecOutDeg + RecNROutAct + RecInDeg + RecNRInAct +
              ## NODE BALANCE
              RecNetBal +
              ## TIME FROM LAST EVENT
              RecLastTime +
              ## ASSORTATIVITY BY DEGREE AND BY INTENSITY
              SdrOutDeg_RecInDeg + SdrNROutAct_RecNRInAct +
              ## REPETITION AND RECIPROCITY
              EdTrt + EdTrtEvtWei + EdRec + EdRecEvtWei + 
              ## DYADIC BALANCE
              EdBalDif + 
              ## TRANSITIVE CLOSURE + CYCLIC CLOSURE + SHARED BORROWERS + SHARED LENDERS
              EdTran + EdClos + EdShrBors + EdShrLens -1,              

            family = "binomial", 
            data = dat.gam)

summary(M_EW3)

## ---------------------------------------------------------
# save(M_EW1, M_EW2, M_EW3, file="01-Data/02-Gam-Fits/event-weighted-model.RData")

## ---------------------------------------------------------
save.image("01-Data/03-Output/output.RData")

