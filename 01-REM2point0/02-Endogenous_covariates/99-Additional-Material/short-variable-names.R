short_term_covariates <- setdiff(setdiff(colnames(M_SLT$model), "one"),
                                 setdiff(colnames(M$model), "one"))
short_term_covariates

colnames(dat.gam)[match(short_term_covariates, colnames(dat.gam))]<- c("ShortRecOutDeg", 
                             "ShortRecNROutAct",
                             "ShortRecInDeg",
                             "ShortRecNRInAct",
                             "SdrShortOutDeg_ShortRecInDeg",
                             "SdrShortNROutAct_ShortRecNRInAct",
                             "ShortEdTrt",
                             "ShortEdRec", 
                             "ShortEdTran",
                             "ShortEdClos",
                             "ShortEdShrBors",
                             "ShortEdShrLens")

save(dat.gam, file="input03.RData")
