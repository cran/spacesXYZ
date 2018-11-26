

#   temperature a vector of temperatures in the interval [4000,25000] K
#
#   returns a matrix with 2 columns:  x, y
#   if any temperature is outside the valid range x,y are set to NA
#

daylightLocus <-  function( temperature, space=1931 )
    {
    ok  = is.numeric(temperature) && 0<length(temperature)
    if( ! ok )
        {
        log.string( ERROR, "temperature='%s' is invalid.  It must be a numeric vector with positive length.", as.character(temperature) )
        return(NULL)
        }
    
    if( ! match(space,c(1960,1976,1931),nomatch=FALSE) )
        {
        log.string( ERROR, "space='%s' is invalid.",  as.character(space[1]) )
        return(NULL)
        }
        
    x = rep( NA_real_, length(temperature) )
    
    t_inv = 1.e3 / temperature    
    
    idx = ( 4000 <= temperature  &  temperature < 7000 )

    x[ idx ] =  0.244063  +  0.09911 * t_inv[idx]  +  2.9678 * t_inv[idx]^2  -  4.6070 * t_inv[idx]^3

    
    idx = ( 7000 <= temperature  &  temperature <= 25000 )
    
    x[ idx ] = 0.237040  +  0.24748 * t_inv[idx]  +  1.9018 * t_inv[idx]^2  -  2.0064 * t_inv[idx]^3

    y   = -3 * x^2  + 2.870 * x  - 0.275
        
    out = cbind( x=x, y=y )

    rnames  = names(temperature)
    if( is.null(rnames) )   rnames = sprintf( "%dK", round(temperature) )
    
    rownames(out)   = rnames
    
    if( space != 1931 )
        out = uvfromxy( out, space=space )
    
    return( out )
    }
    