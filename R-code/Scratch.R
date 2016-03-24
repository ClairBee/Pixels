

find.freq <- function(x)
{
    n <- length(x)
    spec <- spec.ar(c(x),plot=FALSE)
    if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
    {
        period <- round(1/spec$freq[which.max(spec$spec)])
        if(period==Inf) # Find next local maximum
        {
            j <- which(diff(spec$spec)>0)
            if(length(j)>0)
            {
                nextmax <- j[1] + which.max(spec$spec[j[1]:500])
                period <- round(1/spec$freq[nextmax])
            }
            else
                period <- 1
        }
    }
    else
        period <- 1
    return(period)
}

find.freq(sm[1200, 993:1950])
find.freq(sm[1800, 993:1950])
find.freq(sm[1500, 993:1950])
zz <- apply(sm[,993:1950], 1, find.freq)

o.plot(zz)
table(zz)
hist(zz, breaks = "fd", xlim = c(0,50))
