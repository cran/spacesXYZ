
    
############        CCTfrom***()     ###################    
    
#   XYZ     an Nx3 matrix, or a vector that can be converted to such a matrix
#   returns CCT, or NA if outside valid range

CCTfromXYZ    <- function( XYZ, isotherms='robertson', locus='robertson', strict=FALSE  )
    {
    uv  = uvfromXYZ( XYZ, space=1960 )        
    
    if( is.null(uv) )   return(NULL)
    
    out = CCTfromuv( uv, isotherms=isotherms, locus=locus, strict=strict  ) 
        
    return( out )
    }


#   xy  an Nx2 matrix, or a vector that can be converted to such a matrix
    
CCTfromxy  <- function( xy, isotherms='robertson', locus='robertson', strict=FALSE )
    {
    uv  = uvfromxy( xy, space=1960 )        
    
    if( is.null(uv) )   return(NULL)
    
    out = CCTfromuv( uv, isotherms=isotherms, locus=locus, strict=strict  ) 
        
    return(out)
    }


#   uv  an Nx2 matrix, or a vector that can be converted to such a matrix
    
CCTfromuv  <- function( uv, isotherms='robertson', locus='robertson', strict=FALSE )
    {
    uv  = prepareNxM( uv, M=2 )
    if( is.null(uv) )   return(NULL)

    isofull = c( 'native.Rob', 'Robertson', 'McCamy' )
        
    idx.isotherms   = pmatch( tolower(isotherms), tolower(isofull), nomatch=0, duplicates.ok=TRUE )
    if( length(idx.isotherms)==0  ||  any( idx.isotherms==0 ) )
        {
        log.string( ERROR, "isotherms='%s' is invalid.", as.character(isotherms) )
        return( NULL )
        }            
        
    idx.locus   = pmatch( tolower(locus), c('robertson'), nomatch=0 )
    if( idx.locus == 0 )
        {
        log.string( ERROR, "locus='%s' is invalid.", as.character(locus) )
        return( NULL )
        }               

    n   = nrow(uv)
    out = matrix( NA_real_, n, length(idx.isotherms) )
    rownames(out)   = rownames(uv)
    colnames(out)   = isofull[ idx.isotherms ]

    #   assign matrix one element at a time
    for( j in 1:length(idx.isotherms) )
        {
        #   now on column j        
        idx = idx.isotherms[j]
        
        if( idx == 1 )
            {
            #   native Robertson spline
            for( i in 1:n )
                out[i,j]  = CCTfromuv_native( uv[i, ], locus=locus, strict=strict  )
            }
        else if( idx == 2 )
            {
            #   Robertson isotherms
            for( i in 1:n )
                out[i,j]  = CCTfromuv_Robertson( uv[i, ], locus=locus, strict=strict  )
            }        
        else if( idx == 3 )
            {
            #   McCamy
            for( i in 1:n )
                {
                #   McCamy, so convert from uv to xy.
                denom   = uv[i,1] - 4*uv[i,2] + 2
                
                if( is.na(denom) ||  denom < 1.e-16 )   next

                xy      = c( 1.5*uv[i,1], uv[i,2] ) / denom            
                out[i,j]  = CCTfromxy_McCamy( xy, locus=locus, strict=strict )
                }
            }
        }
        
    if( length(idx.isotherms) == 1 )
        {
        #   change nx1 matrix to just a plain vector, and copy rownames to names
        rnames      = rownames(out)
        dim(out)    = NULL
        names(out)  = rnames
        }

    return(out)
    }

    
#   uv      point, usually off locus
#   locus   ignored
#   strict  check distance to locus
#
#   computes where the native isotherm through uv intersects the locus,
#   and returns the native temperature there.
#   It works by searching along the locus, while 'sweeping around' the normal line
#   to find the normal line passing through uv.

CCTfromuv_native  <- function( uv, locus, strict )
    {
    if( any( is.na(uv) ) )  return(NA_real_)
    
    if( TRUE )
        {
        #   not sure whether this special case is worth the time
        #   see whether uv is actually in the LUT !   When it happens, the result is impressive and skips all below.
        idx = which( (uv[1] == p.dataCCT$u) & (uv[2] == p.dataCCT$v) )
        if( length(idx) == 1 )
            #   we have a hit
            return( 1.e6 / p.dataCCT$mired[idx] )
        }
    
    
    
    myfun <- function( mir, uv )
        {
        #   mir         = 1.e6 / temp
        
        uv.locus    = c( p.uvfromMired[[1]]( mir ), p.uvfromMired[[2]]( mir ) )
        
        tangent     = c( p.uvfromMired[[1]]( mir, deriv=1 ), p.uvfromMired[[2]]( mir, deriv=1 ) )
        
        return( sum( tangent * (uv.locus - uv) ) )  # relevant dot product
        }
    
    intervalMired    = range( p.dataCCT$mired )
    
    #   check endpoint values
    f1  = myfun( intervalMired[1], uv )
    f2  = myfun( intervalMired[2], uv )
    
    if( 0 < f1 * f2 )    
        {
        #   same sign, uv is invalid
        log.string( WARN, "uv=%g,%g is in invalid region.  CCT cannot be calculated.", uv[1], uv[2] )
        return( NA_real_ )
        }
        
    if( f1 == 0 )
        mired.end = intervalMired[1]
    else if( f2 == 0 )
        mired.end = intervalMired[2]
    else 
        {
        #log.string( DEBUG, "stats::uniroot() on interval [%g,%g],    endpoint values %g and %g.", 
        #                        intervalMired[1], intervalMired[2],  f1, f2 )
        
        #   reduced the tolerance here, but still only takes about 8 iterations
        res = try( stats::uniroot( myfun, interval=intervalMired, uv=uv, tol=.Machine$double.eps^0.5 ),  silent=FALSE )
        
        if( class(res) == "try-error" )    
            {
            cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
            return( NA_real_ )
            }
            
        mired.end   = res$root   
        
        #log.string( DEBUG, "stats::uniroot() successful after %d iterations.  mired.end=%g", res$iter, mired.end )
        }

    #   temp.end   = res$root

    if( strict )
        {
        #   mired.end   = 1.e6 / temp.end
        uv.locus    = c( p.uvfromMired[[1]]( mired.end ), p.uvfromMired[[2]]( mired.end ) )        
        resid       = uv - uv.locus
        dist        = sqrt( sum(resid^2) )
        
        if( 0.05 < dist )
            {
            log.string( WARN, "uv=%g,%g is invalid, because its distance to the Planckian locus = %.7f > 0.05.  (mired=%g)",
                                uv[1], uv[2], dist, mired.end )
            return( NA_real_ )
            }        
        }

    return( 1.e6 / mired.end )
    }
    
    
#   uv          a 2-vector with uv 1960 
#   locus       only 'robertson' is available, only used if strict=TRUE
#   strict      also compute distance to locus, and rejects if too far
#
#   returns CCT, or NA if strict is TRUE and too far from locus
#
#   requires private data frame p.dataCCT,  which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()
    
CCTfromuv_Robertson  <- function( uv, locus, strict )
    {
    if( any( is.na(uv) ) )  return(NA_real_)
    
    idx.locus  = pmatch( tolower(locus), c('robertson'), nomatch=0 )
    if( idx.locus == 0 )
        {
        log.string( ERROR, "locus='%s' is invalid.", as.character(locus) )
        return(  NA_real_ )
        }

    di = (uv[2] - p.dataCCT$v) - p.dataCCT$t * (uv[1] - p.dataCCT$u)
    #   print( di )
    
    n = length(di)
    
    #   di should be decreasing, and with exactly 1 zero crossing
    i   = 0
    
    #   check for an exact 0
    idx = which( di == 0 )
    
    if( 0 < length(idx) )
        {
        #   exactly 1 is OK
        if( length(idx) == 1 )  i = max( idx[1], 2 )
        }
    else
        {
        #   check for a true crossing
        idx = which( di[1:(n-1)] * di[2:n] < 0 ) + 1
        
        #   exactly 1 is OK
        if( length(idx) == 1 )  i = idx[1]      #     ||  i==1 ) ?
        }
        
    if( i == 0 )
        {
        log.string( WARN, "uv=%g,%g is in invalid region. Cannot find unique zero crossing. CCT cannot be calculated.", uv[1], uv[2] )
        return( NA_real_ )
        }
        
    d0  = di[i]   / sqrt( 1 + p.dataCCT$t[i]^2 )
    dm  = di[i-1] / sqrt( 1 + p.dataCCT$t[i-1]^2 )
    p   = dm / (dm - d0)
    
    mired   = (1-p)*p.dataCCT$mired[i-1]  +  p*p.dataCCT$mired[i]
    
    #   offset   = (1-p) * c(p.dataCCT$u[i-1],p.dataCCT$v[i-1])  +  p * c(p.dataCCT$u[i],p.dataCCT$v[i])  - uv
    
    CCT = 1.e6 / mired                
    
    if( strict )
        {
        res = nativeFromRobertson( CCT, locus )
        
        if( is.null(res) )  return( NA_real_ )
        
        offset  = uv - res$uv   # res$uv is on the locus
        
        test    = sqrt( sum(offset*offset) )
        if( 0.05 < test )
            {
            log.string( WARN, "uv=%g,%g is invalid, because its distance to the Planckian locus = %g > 0.05.  (mired=%g)",
                                uv[1], uv[2], test, mired )
            return( NA_real_ )
            }
        }        

    return( CCT )
    }
    
    
#   CCT     temperature defining a Robertson isotherm, in K.   CCT=Inf is OK, mired is then 0.
#   locus   currently ignored. Using the Robertson spline
#
#   find where the isotherm intersects the locus, using uniroot()
#   in general, we are looking at the spline locus between isotherm i and isotherm i+1; a sector.
#
#   return  a list with:
#           uv      the point of intersection on the locus
#           CCT     the native parameterization CCT in K
#           normal  unit vector along the isotherm, and approximate normal to locus
#   Or NULL in case of error.

nativeFromRobertson <- function( CCT, locus )
    {    
    mired   = 1.e6 / CCT
    
    i   = findInterval( mired, p.dataCCT$mired )
    
    if( i == 0 )
        {
        log.string( WARN, "CCT=%g is out of the Robertson LUT range. mired < %g", CCT, p.dataCCT$mired[1] )
        return(NULL)
        }
            
    n   = length( p.dataCCT$mired )
    if( i == n )
        {
        #   an extra check. 
        #   although I have found that 1.e6 / (1.e6 / 600) == 600 exactly on my PC, this might not be TRUE on other hardware,
        #   so use an epsilon just in case.
        epsilon = 1.e-12
        if( p.dataCCT$mired[n] + epsilon  <  mired )
            {
            log.string( WARN, "CCT=%g is outside the Robertson LUT range. %g < mired", CCT, p.dataCCT$mired[n] )
            return(NULL)
            }     

        #   otherwise set mired to what it probably should be (600),
        #   and keeping going with i == n, and the next special if block will take care of it !
        mired = p.dataCCT$mired[n] 
        }

    out = list()
    
    if( p.dataCCT$mired[i] == mired )    
        {
        #   special case, CCT is at one of the defining isotherms.  root at endpoint so avoid using uniroot()
        #log.string( DEBUG, "mired=%g exactly. isotherm %d.  uniroot unnecessary.", mired, i )
        
        out$uv  = c( p.dataCCT$u[i], p.dataCCT$v[i] )
        out$CCT = CCT      
        
        tani        = c( 1, p.dataCCT$t[i] )    # tangent vector to isotherm
        len         = sqrt( sum(tani^2) )
        out$normal  = -tani / len               # normal to locus (approx)
        
        return(out)
        }
        
    #   we are looking at the spline locus between isotherm i and isotherm i+1. a sector.
    #   See Fig. 2(3.11) in W&S. 
        
    #   extract the 2 points
    uvi     = c( p.dataCCT$u[i], p.dataCCT$v[i] )
    uvi1    = c( p.dataCCT$u[i+1], p.dataCCT$v[i+1] )
      
    #   precompute unit normals to the isotherms at both ends.  Maybe these should be moved into the LUT someday.
    #   these are very close to the unit tangents to the locus.
    tani    = c( 1, p.dataCCT$t[i] )            # tangent vector to isotherm
    len     = sqrt( sum(tani^2) )
    normi   = c( -tani[2], tani[1] ) / len      # unit normal to isotherm at knot i
    
    tani1   = c( 1, p.dataCCT$t[i+1] )          # tangent vector to isotherm
    len     = sqrt( sum(tani1^2) )
    normi1  = c( -tani1[2], tani1[1] ) / len    # unit normal to isotherm at knot i+1
        
    myfun   <- function( mired, CCT )
        {
        uv  = c( p.uvfromMired[[1]]( mired ), p.uvfromMired[[2]]( mired ) )
        
        deltai  = uv - uvi
        di      = sum( deltai*normi )
        deltai1 = uv - uvi1
        di1     = sum( deltai1*normi1 )
        
        denom   = di - di1
        
        ok  = (0 <= di)  &&  (di1 <= 0) && (0 < denom)
        if( ! ok )
            log.string( ERROR, "Bad signs. di=%g and di1=%g. mired=%g. i=%d", di, di1, mired, i )
    
        alpha   = di / denom   # this must be in closed interval [0,1]
    
        mir = (1-alpha)*p.dataCCT$mired[i]  +  alpha*p.dataCCT$mired[i+1]
            
        return( 1.e6/mir - CCT )
        }
    
    rangeMired    = c( p.dataCCT$mired[i], p.dataCCT$mired[i+1] )

    #log.string( DEBUG, "stats::uniroot() on interval [%g,%g], isotherms=%d and %d.  endpoint values %g and %g.", 
    #                        rangeMired[1], rangeMired[2], i, i+1, myfun(rangeMired[1],CCT), myfun(rangeMired[2],CCT)  )

    res = try( stats::uniroot( myfun, interval=rangeMired, CCT ),  silent=FALSE )
    
    if( class(res) == "try-error" )    
        {
        cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
        return( NULL )
        }        
        
    #log.string( DEBUG, "stats::uniroot() successful after %d iterations.  isotherm %d", res$iter, i )

    out$uv  = c( p.uvfromMired[[1]]( res$root ), p.uvfromMired[[2]]( res$root ) )
    out$CCT = 1.e6 / res$root
    
    #   for the unit normal go through distance calculation and interpolation again !
    deltai  = out$uv - uvi
    di      = sum( deltai*normi )
    deltai1 = out$uv - uvi1
    di1     = sum( deltai1*normi1 )
    denom   = di - di1    
    alpha   = di / denom   # this must be in closed interval [0,1]
    tangent = (1-alpha)*normi  +  alpha*normi1
    normal  = c( -tangent[2], tangent[1] )      # tangent to locus
    len     = sqrt( sum(normal^2) )
    
    out$normal  = normal / len      #; print( out$normal )
            
    return( out )
    }
        
    
#   C. S. McCamy
#   Correlated color temperature as an explicit function of chromaticity coordinates.
#   Color Research & Application
#   Volume 17, Issue 2, pages 142-144, 
#   April 1992
#
#   'nocheck' means do not check that xy is finite, and do not check distance from locus

CCTfromxy_McCamy_nocheck  <- function( xy )
    {
    topbot  = xy - c(0.3320,0.1858)     # the 2nd term is the point where all the isotherms intersect
    
    if( topbot[2] <= 0 )  return( NA_real_ )
    
    n   = topbot[1]/topbot[2]
    
    out = -449*n^3 + 3525*n^2 - 6823.3*n + 5520.33
    
    if( out <= 0 )  return( NA_real_ )  # this can happen
    
    return( out )
    }
    
CCTfromxy_McCamy  <- function( xy, locus, strict )
    {
    if( any( is.na(xy) ) )  return( NA_real_ )
    
    CCT = CCTfromxy_McCamy_nocheck( xy )
    
    if( is.finite(CCT)  &&  strict )
        {
        #   find intersection of McCamy isotherm and the locus
        res = nativeFromMcCamy( CCT, locus )
        
        if( is.null(res) )  return( NA_real_ )
        
        uv  = c( 4*xy[1] , 6*xy[2] ) / ( -2*xy[1]  + 12*xy[2]  + 3 )
        
        offset  = uv - res$uv
        dist    = sqrt( sum(offset^2) )
        if( 0.05 < dist )
            {
            log.string( WARN, "uv=%g,%g is invalid, because its distance to the Planckian locus = %g > 0.05.  CCT=%g",
                                uv[1], uv[2], dist, CCT  )
            return( NA_real_ )
            }
        }
    
    return(CCT)
    }
    
    
############        planckLocus()     ###################    

#   temperature a N-vector of temperatures, in K.  Inf is OK
#   locus       must match 'robertson'
#   param       'native, 'robertson', or 'mccamy'
#   delta       offset from the locus, in top direction.  Always along the corresponding isotherm. Distance in uv-1960 plane.
#   space       year of chromaticity space
#
#   returns     an Nx2 matrix of uv's, or NULL in case of argument error
#
planckLocus  <- function( temperature, locus='robertson', param='robertson', delta=0, space=1960 )
    {
    ok  = is.numeric(temperature)  &&  0<length(temperature)
    if( ! ok )
        {
        log.string( ERROR, "Argument temperature is invalid.  It must be a numeric vector with positive length."  )
        return(NULL)
        }

    idx.locus   = pmatch( tolower(locus), c('robertson'), nomatch=0 )
    if( idx.locus == 0 )
        {
        log.string( ERROR, "locus='%s' is invalid.", as.character(locus) )
        return( NULL )
        }        
        
    idx.param  = pmatch( tolower(param), c('native','robertson','mccamy'), nomatch=0 )
    if( idx.param == 0 )
        {
        log.string( ERROR, "param='%s' is invalid.", as.character(param) )
        return( NULL )
        }
        
    n   = length(temperature)
    
    ok  = is.numeric(delta)  &&  length(delta) %in% c(1,n)
    if( ! ok )
        {
        log.string( ERROR, "Argument delta is invalid.  It must be a numeric vector with length 1 or %d.", n  )
        return(NULL)
        }
    if( length(delta) == 1 )    delta = rep( delta[1], n )

    
    if( ! match(space,c(1960,1976,1931),nomatch=FALSE) )
        {
        log.string( ERROR, "space='%s' is invalid.",  as.character(space[1]) )
        return(NULL)
        }    
    
    uv  = matrix( NA_real_, n, 2 )
    
    rnames  = names(temperature)
    if( is.null(rnames) )   rnames = sprintf( "%gK", round(temperature) )
        
    rownames(uv)    = rnames
    colnames(uv)    = c('u','v')        
        
    # amazingly, 600 seems to involute perfectly with 1.e6, though not with 1        
    temperaturerange    = c( 1.e6 / p.dataCCT$mired[ length(p.dataCCT$mired) ], Inf )        # 33334  100000 )
    
    ok  =   temperaturerange[1]<=temperature  &  temperature<=temperaturerange[2]
    temperature[ ! ok ] = NA_real_

    
    if( idx.param == 1 )
        {
        #   the easy one - plain 'native' spline
        mired   = 1.e6 / temperature            # temperature==Inf is OK here, mired is then 0
        uv[ ,1]   = p.uvfromMired[[1]]( mired )
        uv[ ,2]   = p.uvfromMired[[2]]( mired )
        
        if( any( delta!=0 ) )
            {
            #   handle the delta offset
            for( i in 1:n )
                {
                if( delta[i] == 0  ||  is.na(mired[i]) ) next
                
                #   compute unit normal to the locus at point i
                normal  = c( -p.uvfromMired[[2]]( mired[i], deriv=1 ), p.uvfromMired[[1]]( mired[i], deriv=1 ) )
                len     = sqrt( sum(normal^2) )
                if( len == 0 )  next    # should never happen
                normal  = normal / len
                
                uv[i, ] = uv[i, ] + delta[i]*normal
                }
            }
        }    
    else
        {
        for( i in 1:n )
            {
            if( is.na(temperature[i]) ) next
            
            #if( is.infinite(temperature[i]) )
            #    {
            #    #   special case
            #    uv[ ,1]   = p.uvfromMired[[1]]( 0 )
            #    uv[ ,2]   = p.uvfromMired[[2]]( 0 )
            #    next
            #    }
            
            if( idx.param == 3 )
                {
                res = nativeFromMcCamy( temperature[i], locus )
                }
            else
                {
                #   has to be Robertson
                res = nativeFromRobertson( temperature[i], locus )     # also works if temperature[i] is Inf
                }
                
            if( is.null(res) ) next
            
            uv[i, ] = res$uv
            
            #   handle the delta offset
            if( is.finite(delta[i])  &&  delta[i]!=0 )
                {
                uv[i, ] = uv[i, ] + delta[i] * res$normal
                }
            }
        }
        
    #   now handle the space
    if( space == 1931 )   
        {
        xy  = uv    # get the size and rownames right
        colnames(xy) = c('x','y')
        
        denom   = uv[ ,1] - 4*uv[ ,2] + 2
        
        xy[ ,1] = 1.5*uv[ ,1] / denom
        xy[ ,2] = uv[ ,2] / denom        
        
        return( xy )
        }
    else if( space == 1976 )
        {
        colnames(uv) = c("u'","v'")
        uv[ ,2] = 1.5 * uv[ ,2]
        }

    return(uv)
    }
    
    
#   CCT     temperature defining a McCamy isotherm, in K
#   locus   currently ignored
#
#   find where the isotherm intersects the (Robertson) locus, using uniroot()
#
#   return  a list with:
#           uv      the point of intersection on the locus
#           CCT     the native parameterization CCT in K
#           normal  to the locus, always along the McCamy isotherm for CCT
#       or NULL in case of error.
nativeFromMcCamy <- function( CCT, locus )
    {
    #   the following limits were computed from the endpoints of the Robertson lookup table
    #   and rounded to the appropriate integer
    ok  = 1702 <= CCT  &&  CCT <= 34539
    
    if( ! ok )  return(NULL)    # too small or too big for McCamy           33333
       
    myfun <- function( mired, CCT )
        {
        uv  = c( p.uvfromMired[[1]](mired), p.uvfromMired[[2]](mired) )
        
        xy  = c( 1.5*uv[1], uv[2] ) / ( uv[1] - 4*uv[2] + 2 )    
        
        return( CCTfromxy_McCamy_nocheck(xy) - CCT )
        }
        
    rangeMired    = range( p.dataCCT$mired )
        
    #   log.string( DEBUG, "myfun() at endpoints: %g and %g.", myfun(rangeMired[1],CCT), myfun(rangeMired[2],CCT) )
    
    res = try( stats::uniroot( myfun, interval=rangeMired, CCT ),  silent=FALSE )
    
    if( class(res) == "try-error" )    
        {
        cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
        return( NULL )
        }        
    # log.string( DEBUG, "uniroot() successful after %d iterations, for McCamy.", res$iter )
    
    uv  = c( p.uvfromMired[[1]]( res$root ), p.uvfromMired[[2]]( res$root ) )
    
    out = list()
    out$uv  = uv
    out$CCT = 1.e6 /  res$root
    
    # normal of locus, along the isotherm
    xy  = c( 1.5*uv[1], uv[2] ) / ( uv[1] - 4*uv[2] + 2 )    
        
    meet    = c(0.3320,0.1858) #  the point where all the isotherms meet
    topbot  = xy - meet
    
    alpha       = topbot[1] / topbot[2]
    C           = sum( c(1,-alpha) * meet )
    normal      = c( alpha - 4*C, 1.5 - C )
    len         = sqrt( sum(normal^2) )
    out$normal  = normal / len          #; print( out$normal )
    
    return(out)
    }
    
    
