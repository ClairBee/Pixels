
# PARAMETRIC MODEL & SD USED TO IDENTIFY BAD PIXELS WITH JOHNSON QUANTILES

library("IO.Pixels")
load.pixel.means(); load.pixel.sds()

####################################################################################################
# # change ordering on 'dead' & 'hot': dim/bright takes priority (more stringent condition)

# Also need to find way to identify lines of bright pixels, which may not be picked up here
####################################################################################################

# BAD PIXELS - 14.10.09                                                                         ####

# fit parametric model & use Johnson quantiles to get bad pixel map (dead/dim/normal/bright/hot)

# use Johnson quantiles of SD to pick up further bad pixels (static/normal/erratic)

# compare to official map

####################################################################################################



# BAD PIXELS - 14.11.18                                                                         ####

# fit parametric model & use Johnson quantiles to get bad pixel map (dead/dim/normal/bright/hot)

# use Johnson quantiles of SD to pick up further bad pixels (static/normal/erratic)

# compare to official map

# compare to previous map

# examine new bad pixels in detail

####################################################################################################


