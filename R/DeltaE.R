

#   Lab1        Nx3 matrix,
#   Lab2        Nx3 matrix,
#   metric      which metric to use
#
#   returns     numeric N-vector of pairwise differences

DeltaE  <-  function( Lab1, Lab2, metric=1976 )
    {
    Lab1 = prepareNxM( Lab1 )
    if( is.null(Lab1) )  return(NULL)

    Lab2 = prepareNxM( Lab2 )
    if( is.null(Lab2) )  return(NULL)

    if( nrow(Lab1)==1  &&  1<nrow(Lab2) )
        #   replicate Lab1
        Lab1    = matrix(Lab1,nrow(Lab2),3,byrow=TRUE)

    if( 1<nrow(Lab1)  &&  nrow(Lab2)==1 )
        #   replicate Lab2
        Lab2    = matrix(Lab2,nrow(Lab1),3,byrow=TRUE)

    if( nrow(Lab1) != nrow(Lab2) )
        {
        log_level( ERROR, "nrow(Lab1) = %d != %d = nrow(Lab2).",
                            nrow(Lab1), nrow(Lab2) )
        return(NULL)
        }

    #   Lab1 and Lab2 passed all checks

    #   find names if possible
    rnames  = NULL
    if( ! is.null(rownames(Lab1) ) )
        rnames  = rownames(Lab1)
    else if( ! is.null(rownames(Lab2) ) )
        rnames  = rownames(Lab2)

    metricname  = c('1976','1994','2000')
    metric  = as.character(metric)
    idx.metric  = match( metric, metricname, nomatch=0 )
    if( length(idx.metric)==0  ||  any(idx.metric==0) )
        {
        metric  = metric[ idx.metric==0 ]
        log_level( ERROR, "metric='%s' is invalid.", metric )
        return( NULL )
        }

    n   =   nrow(Lab1)  # same as Lab2
    out = matrix( NA_real_, n, length(idx.metric) )
    rownames(out)   = rnames

    if( 1 < length(idx.metric) )
        colnames(out)   = sprintf( "DeltaE.%s", metricname[idx.metric] )

    for( j in 1:length(idx.metric) )
        {
        idx = idx.metric[j]

        if( idx == 1 )
            #   1976.  do the whole column in one shot - vectorized
            out[ ,j] = DeltaE.1976( Lab1, Lab2 )
        else if( idx == 2 )
            {
            #   1994.  one pair at a time
            for( i in 1:n )
                out[i,j] = DeltaE.1994( Lab1[i,], Lab2[i,] )
            }
        else if( idx == 3 )
            {
            #   2000.  one pair at a time
            for( i in 1:n )
                out[i,j] = DeltaE.2000( Lab1[i,], Lab2[i,] )
            }
        }

    if( length(idx.metric) == 1 )
        {
        #   change nx1 matrix to just a plain vector, and copy rownames to names
        dim(out)    = NULL
        names(out)  = rnames
        }

    return(out)
    }


#   Lab1    Nx3 matrix
#   Lab2    Nx3 matrix
DeltaE.1976  <-  function( Lab1, Lab2 )
    {
    delta   = Lab1 - Lab2

    sqrt( rowSums(delta*delta) )
    }


#   Lab1           a 3-vector
#   Lab2           a 3-vector
#   .KL .KC .KH    weighting factors
#
#   returns        color difference by CIED2000  - insanely complicated !

DeltaE.2000  <-  function( Lab1, Lab2,  .KL=1, .KC=1, .KH=1 )
    {
    Lbar = (Lab1[1] + Lab2[1]) / 2

    C1 = sqrt( Lab1[2]^2  + Lab1[3]^2 )
    C2 = sqrt( Lab2[2]^2  + Lab2[3]^2 )
    Cbar = (C1 + C2) / 2

    g   = sqrt( Cbar^7 / (Cbar^7 + 25^7) )
    G   = (1 - g) / 2

    a1 = Lab1[2] * (1 + G)
    a2 = Lab2[2] * (1 + G)


    #   recompute chromas, drop the primes
    C1 = sqrt( a1^2  + Lab1[3]^2 )
    C2 = sqrt( a2^2  + Lab2[3]^2 )
    Cbar = (C1 + C2) / 2

    #mess = sprintf( "C1=%g   C2=%g\n", C1,  C2 )
    #cat( mess )

    h1 = atan2( Lab1[3], a1 ) %% (2*pi)
    h2 = atan2( Lab2[3], a2 ) %% (2*pi)

    #mess = sprintf( "h1=%g   h2=%g\n", h1,  h2 )
    #cat( mess )

    Hbar    = (h1 + h2)/2
    hdelta  = h2 - h1
    if( pi < abs(h1 - h2) )
        {
        Hbar    = Hbar + pi
        hdelta  = hdelta  - sign(hdelta)*2*pi
        }

    deg_to_rad  = pi/180

    T = 1 - 0.17*cos(Hbar - 30*deg_to_rad) + 0.24*cos(2*Hbar) + 0.32*cos(3*Hbar + 6*deg_to_rad) - 0.20*cos(4*Hbar - 63*deg_to_rad)

    Ldelta  = Lab2[1] - Lab1[1]
    Cdelta  = C2 - C1
    Hdelta  = 2 * sqrt(C1*C2) * sin( hdelta / 2 )

    #mess = sprintf( "hdelta=%g   Hdelta=%g\n", hdelta,  Hdelta )
    #cat( mess )

    SL  = (Lbar - 50)^2
    SL  = 1 + 0.015 * SL / sqrt(20 + SL)
    SC  = 1 + 0.045 * Cbar
    SH  = 1 + 0.015 * Cbar * T

    theta_delta = 30*deg_to_rad * exp( -((Hbar*180/pi - 275)/25)^2 )    # in radians

    g   = sqrt( Cbar^7 / (Cbar^7 + 25^7) )  #   recall Cbar is really Cbar'
    RC  = 2*g
    RT  = -RC * sin( 2*theta_delta )

    Lterm   = Ldelta / (.KL * SL)
    Cterm   = Cdelta / (.KC * SC)
    Hterm   = Hdelta / (.KH * SH)

    out = sqrt( Lterm^2 + Cterm^2 + Hterm^2  +  RT*Cterm*Hterm )

    return( out )
    }


DeltaE.1994  <-  function( Lab1, Lab2,  .KL=1, .K1=0.045, .K2=0.015 )
    {
    Delta   = Lab1  - Lab2
    C1      = sqrt( sum(Lab1[2:3]^2) )
    C2      = sqrt( sum(Lab2[2:3]^2) )
    Delta.C = C1 - C2
    
    #   Dr. Ben Jann found that the next line was missing the 'sum', and so return value very wrong.
    #   Bug was fixed in v 1.2-x
    Delta.H = sqrt( max( sum(Delta[2:3]^2) - Delta.C^2, 0) )  # the max() is necessary because of non-zero numerical precision

    SL  = 1
    #   SC  = 1 + .K1*C1    asymmetric variant
    #   SH  = 1 + .K2*C1    asymmetric variant

    C12 = sqrt(C1*C2)   # symmetric variant, from Hunt - The Reproduction of Colour - 6th edition
    SC  = 1 + .K1*C12   # symmetric variant
    SH  = 1 + .K2*C12   # symmetric variant

    tmp = c( Delta[1], Delta.C, Delta.H ) / c( .KL*SL, SC, SH )

    return( sqrt( sum(tmp^2) ) )
    }


