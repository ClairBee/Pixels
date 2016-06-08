
library(msm)

####################################################################################################

# 2 - FITTING MULTI-STATE MODELS                                                                ####

# by default, the data are assumed to represent snapshots of the process at arbitrary times.

# pairwise transition matrix (state table)
statetable.msm(state, PTNUM, data = cav)

# define allowed state transitions in Q matrix format (diagonals are trivial)
Q <- rbind (c(0, 0.25, 0, 0.25),
            c(0.166, 0, 0.166, 0.166),
            c(0, 0.25, 0, 0.25),
            c(0, 0, 0, 0) )

# model 1: simple bidirectional model
{
    # estimate values of initial Q matrix before attempting maximum likelihood estimation
    Q.crude <- crudeinits.msm(state ~ years, PTNUM, data = cav, qmatrix = Q)
    
    # fit the model (using initial user-estimated Q matrix)
    cav.msm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = Q, deathexact = 4)
    
    # alternatively, get crude estimates of initial values as part of fitting procedure
    cav.msm.crude <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = Q, deathexact = 4, gen.inits = T)
    
    cav.msm; cav.msm.crude
    # actually very similar.
    # From the estimated intensities, we see patients are three times as likely to develop symptoms 
    # than die without symptoms (transitions from state 1). After disease onset (state 2), progression 
    # to severe symptoms (state 3) is 50% more likely than recovery. Once in the severe state, death is
    # more likely than recovery, and a mean of 1 / -0.44 = 2.3 years is spent in state 3 before death or recovery.
}


# model 2: sex as a covariate
{
    cavsex.msm <- msm( state ~ years, subject = PTNUM, data = cav, qmatrix = Q, deathexact = 4, 
                       covariates = ~ sex)
    
    cavsex.msm
    
    # The sizes of the confidence intervals for some of the hazard ratios suggests there is no information
    # in the data about the corresponding covariate effects, leading to a likelihood that is a flat function
    # of these parameters, and this model should be simplified. The first column shown in the output is
    # the estimated transition intensity matrix, with the covariate z set to its mean
    # value in the data. This represents an average intensity matrix for the population of 535 male and
    # 87 female patients. To extract separate intensity matrices for male and female patients (z = 0 and 1
    # respectively), use the function qmatrix.msm, as shown below. 
    
    qmatrix.msm(cavsex.msm, covariates = list(sex = 0)) # male
    qmatrix.msm(cavsex.msm, covariates = list(sex = 1)) # female
}

# model 2a: transition-specific covariates
{
    # different transition rates may be easily modelled on different covariates
    # by specifying a named list of formulae as the covariates argument.
    
    # In the model below, the transition rate from state 1 to state 2 and
    # the rate from state 1 to state 4 are each modelled on sex as a covariate, but no other intensities have
    # covariates on them.
    
    cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                       qmatrix = Q, deathexact = 4,
                       covariates = list("1-2" = ~ sex, "1-4" = ~sex) )
    
}


# model 3: constrained covariate effects
{
    # We may also want to constrain the effect of a covariate to be equal for certain transition rates, to
    # reduce the number of parameters in the model, or to investigate hypotheses on the covariate effects. A
    # constraint argument can be used to indicate which of the transition rates have common covariate
    # effects.
    
    cav3.msm <- msm( state ~ years, subject = PTNUM, data = cav,
                     qmatrix = Q, deathexact = 4,
                     covariates = ~ sex,
                     constraint = list(sex = c(1,2,3,1,2,3,2)) )
    
    # This constrains the effect of sex to be equal for the progression rates q12 , q23, equal for the
    # death rates q14 , q24 , q34, and equal for the recovery rates q21 , q32 . The intensity parameters are
    # assumed to be ordered by reading across the rows of the transition matrix, starting at the first row:
    # (q12 , q14 , q21, q23 , q24 , q32, q34 ), giving constraint indicators (1,2,3,1,2,3,2). Any vector of
    # increasing numbers can be used for the indicators. Negative entries can be used to indicate that some
    # effects are equal to minus others: (1,2,3,-1,2,3,2) sets the fourth effect to be minus the first.
}


# model 4: fixed parameters
{
    # For exploratory purposes we may want to fit a model assuming that some parameters are fixed, and
    # estimate the remaining parameters. 
        
    cav4.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     qmatrix = Q, deathexact = 4,
                     control = list(trace=2, REPORT=1),
                     fixedpars = c(6, 7) )
    
}

qmatrix.msm(cav.msm)                # intensity matrix
pmatrix.msm(cav.msm)                # transition probability matrix
pmatrix.msm(cav.msm, t = 10)        # transition probability matrix after time t

# confidence intervals are computationally intensive but can be obtained using
pmatrix.msm(cav.msm, t = 10, ci = "norm")       # faster (~1sec)
pmatrix.msm(cav.msm, t = 10, ci = "boot")       # takes FOR EVER. Stopped at 9 minutes.

    # Most of msm’s output functions present confidence intervals based on asymptotic standard errors 
    # calculated from the Hessian, or transformations of these using the delta method. The asymptotic
    # standard errors are expected to be underestimates of the true standard errors (Cramer-Rao lower bound).


pnext.msm(cav.msm)                  # probability that each state is next (more intuitive parameterisation?)
sojourn.msm(cav.msm)                # mean sojourn times in each state
totlos.msm(cav.msm)                 # total length of stay in each state

efpt.msm(cav.msm, tostate = 4)      # expected first passage time (hitting time)
envisits.msm(cav.msm)               # expected number of visits

qratio.msm(cav.msm, ind1=c(2,1), ind2=c(1,2))       # ratio of transition intensities
                                                    # recovery is 1.8 times more likely than progression

hazard.msm(cavsex.msm)

# 2 - MODEL ASSESSMENT                                                                          ####    

# survival plot
plot(cav.msm, legend.pos = c(8, 1))

    # This shows that the 10-year survival probability with severe CAV is approximately 0.1, as opposed
    # to 0.3 with mild CAV and 0.5 without CAV. With severe CAV the survival probability diminishes very
    # quickly to around 0.3 in the first five years after transplant.


# bootstrapping - takes ages (stopped after ~8m)
system.time(q.list <- boot.msm(cav.msm, stat = function(x){qmatrix.msm(x)$estimates}))


# contribution of each subject to log-likelihood
logLik.msm(cav.msm, by.subject=TRUE)

# If only a few subjects give an infinite log-likelihood, then you can check whether their state 
# histories are particularly unusual and conflict with the model. For example, they might appear to
# make unusually large jumps between states in short periods of time. For models with misclassification,
# note that the default true initial state distribution initprobs puts all individuals in
# true state 1 at their first observation. If someone starts in a much higher state, this may result
# in an infinite log-likelihood, and changing initprobs would be sensible.

# To compare the relative fit of two nested models, it is easy to compare their likelihoods. However 
# it is not always easy to determine how well a fitted multi-state model describes an irregularly-observed
# process. Ideally we would like to compare observed data with fitted or expected data under the model. 
# If there were times at which all individuals were observed then the fit of the expected numbers 
# in each state (or prevalences) can be assessed directly at those times.
prevalence.msm(cav.msm, times=seq(0,20,2))
plot.prevalence.msm(cav.msm, mintime=0, maxtime=20)

# A set of expected counts can be produced if the process begins at a common time for all individuals.
    
pearson.msm(cav.msm, timegroups = 2, transitions = c(1,2,3,4,5,6,7,8,9,9,9,10))

# model 5: multi-state model with misclassification
{
    Qm <- rbind(c(0, 0.148, 0, 0.0171),
                c(0, 0, 0.202, 0.081),
                c(0, 0, 0, 0.126),
                c(0, 0, 0, 0))
    ematrix <- rbind(c(0, 0.1, 0, 0),
                     c(0.1, 0, 0.1, 0),
                     c(0, 0.1, 0, 0),
                     c(0, 0, 0, 0))
    cavmisc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                       qmatrix = Qm, ematrix = ematrix, deathexact = 4,
                       obstrue = firstobs)
    cavmisc.msm
    plot.prevalence.msm(cavmisc.msm, mintime=0, maxtime=20)
    
}

# model 6: misclassification model with misclassifications modelled on sex
{
    cavmiscsex.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                          qmatrix = Qm, ematrix = ematrix,
                          deathexact = 4, misccovariates = ~sex,
                          obstrue=firstobs)
    cavmiscsex.msm
    
    # The large confidence interval for the odds ratio for 1/2 misclassification suggests there is no informa-
    # tion in the data about the difference between genders in the false positive rates for angiography. On
    # the other hand, women have slightly more false negatives.
    
    ematrix.msm(cavmiscsex.msm, covariates = list(sex = 0))     # misclassifications for male
    ematrix.msm(cavmiscsex.msm, covariates = list(sex = 1))     # misclassifications for female

    # odds ratio for misclassification    
    odds.msm(cavmiscsex.msm)
    
}

vit <- viterbi.msm(cavmisc.msm)     # reconstruct most likely state path
vit[vit$subject==100103,]           # state path (observed & most likely) for single patient

# The results for patient 100103 are shown, who appeared to ‘recover’ to a less severe state of disease while
# in state 3. We assume this is not biologically possible for the true states, so we expect that either the
# observation of state 3 at time 4.98 was an erroneous observation of state 2, or their apparent state
# 2 at time 5.94 was actually state 3. According to the expected path constructed using the Viterbi
# algorithm, it is the observation at time 5.94 which is most probably misclassified.

