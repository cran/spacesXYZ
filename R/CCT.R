
    
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

    #   process isotherms    
    if( length(isotherms) == 0 )
        {
        log.string( ERROR, "isotherms is invalid, because length(isotherms)=0."  )
        return( NULL )
        }      
        
    isofull = c( 'native', 'Robertson', 'McCamy' )
        
    isotherms   =  as.character(isotherms)
    
    isotherms[ is.na(isotherms) ]   = 'native'
        
    idx.isotherms   = pmatch( tolower(isotherms), tolower(isofull), nomatch=0, duplicates.ok=TRUE )
    if( any( idx.isotherms==0 ) )
        {
        log.string( ERROR, "isotherms='%s' is invalid.", isotherms[idx.isotherms==0] )
        return( NULL )
        }            
        
    #   process locus
    ok  = is.character(locus)  &&  length(locus)==1
    if( ! ok )
        {
        log.string( ERROR, "Argument locus is invalid.  It must be a character vector with length 1." )
        return(NULL)
        }

    locusfull = c( 'Robertson', 'precision' )        
                     
    idx.locus   = pmatch( tolower(locus), tolower(locusfull), nomatch=0 )
    if( idx.locus == 0 )
        {
        log.string( ERROR, "locus='%s' is invalid.", as.character(locus) )
        return( NULL )
        }            

    if( idx.locus == 1 )
        locus.list  = p.uvCubicsfromMired
    else
        locus.list  = p.uvQuinticsfromMired
        

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
                out[i,j]  = CCTfromuv_native( uv[i, ], locus.list, strict  )
            }
        else if( idx == 2 )
            {
            #   Robertson isotherms
            for( i in 1:n )
                out[i,j]  = CCTfromuv_Robertson( uv[i, ], locus.list, strict  )
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
                out[i,j]  = CCTfromxy_McCamy( xy, locus.list, strict )
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

    
#   uv          point, usually off locus
#   locus.list  ufun, vfun, and miredInterval
#   strict      check distance to locus
#
#   computes where the native isotherm through uv intersects the locus,
#   and returns the native temperature there.
#   It works by searching along the locus, while 'sweeping around' the normal line
#   to find the normal line passing through uv.

CCTfromuv_native  <- function( uv, locus.list, strict )
    {
    if( any( is.na(uv) ) )  return(NA_real_)
    
    #   log.string( TRACE, "uv=%g,%g", uv[1], uv[2] )
    #   print( str(locus.list) )
    
    
    if( FALSE )
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
        
        uv.locus    = c( locus.list$ufun( mir ), locus.list$vfun( mir ) )
        
        tangent     = c( locus.list$ufun( mir, deriv=1 ), locus.list$vfun( mir, deriv=1 ) )
        
        return( sum( tangent * (uv.locus - uv) ) )  # relevant dot product
        }
    
    miredInterval    = locus.list$miredInterval           # range( p.dataCCT$mired )
    
    #   check endpoint values
    f1  = myfun( miredInterval[1], uv )
    f2  = myfun( miredInterval[2], uv )
    
    if( 0 < f1 * f2 )    
        {
        #   same sign, uv is invalid
        log.string( WARN, "uv=%g,%g is in invalid region.  CCT cannot be calculated.", uv[1], uv[2] )
        return( NA_real_ )
        }
        
    if( f1 == 0 )
        mired.end = miredInterval[1]
    else if( f2 == 0 )
        mired.end = miredInterval[2]
    else 
        {
        #log.string( DEBUG, "stats::uniroot() on interval [%g,%g],    endpoint values %g and %g.", 
        #                        miredInterval[1], miredInterval[2],  f1, f2 )
        
        #   reduced the tolerance here, but still only takes about 8 iterations
        res = try( stats::uniroot( myfun, interval=miredInterval, uv=uv, tol=.Machine$double.eps^0.5 ),  silent=FALSE )
        
        if( inherits(res,"try-error" ) )    # class(res) == "try-error" )    
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
        uv.locus    = c( locus.list$ufun( mired.end ),locus.list$vfun( mired.end ) )        
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
    
CCTfromuv_Robertson  <- function( uv, locus.list, strict )
    {
    if( any( is.na(uv) ) )  return(NA_real_)
    
    mired   = miredfromuv_Robertson_nocheck( uv, extrap=FALSE )
    
    if( is.na(mired) )  return(NA_real_)
    
    CCT = 1.e6 / mired                
    
    if( strict )
        {
        res = nativeFromRobertson( CCT, locus.list )
        
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
    
    
  
#   uv          a 2-vector with uv 1960 
#   extrap      TRUE means extrapolate past mired=0 and 600.  FALSE means return NA_real
#   requires private data frame p.dataCCT,  which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()
#   'nocheck' means do not that input is NA and do not check locus

miredfromuv_Robertson_nocheck  <- function( uv, extrap=FALSE )
    {
    di = (uv[2] - p.dataCCT$v) - p.dataCCT$t * (uv[1] - p.dataCCT$u)
    #   print( di )
    
    n = length(di)
    
    #   di should be decreasing, and with exactly 1 zero crossing
    #   i   = 0
    
    #   check for an exact 0
    idx = which( di == 0 )
    
    if( 0 < length(idx) )
        {
        if( length(idx) == 1 )  
            #   exactly 1 is OK        
            return( p.dataCCT$mired[idx] )
        
        #   more than 1 should not happen
        log.string( WARN, "uv=%g,%g lies on more than 1 isotherm. Mired cannot be calculated.", uv[1], uv[2] )
        return( NA_real_ )
        }

    #   check for a true crossing
    idx = which( di[1:(n-1)] * di[2:n] < 0 )
    
    if( length(idx) == 0 )    
        {
        #   all di must be pos or neg
        if( ! extrap )
            {
            log.string( WARN, "uv=%g,%g is in invalid region. Found no zero-crossings. Mired cannot be calculated.", uv[1], uv[2] )
            return( NA_real_ )
            }        
            
        #   do a very crude extrapolation, only useful for root-finding
        if( di[1] < 0 )
            return( p.dataCCT$mired[1] + 100 * di[1] )
        else
            return( p.dataCCT$mired[n] + 30 * di[n] )
        }
    else if( 2 <= length(idx) )        #     ||  i==1 ) ?
        {
        log.string( WARN, "uv=%g,%g is in invalid region. Found multiple zero-crossings. Mired cannot be calculated.", uv[1], uv[2] )
        return( NA_real_ )
        }        
        
    #   exactly 1 zero-crossing is OK.  This is the normal situation        
    i = idx[1] + 1
    
    d0  = di[i]   / sqrt( 1 + p.dataCCT$t[i]^2 )
    dm  = di[i-1] / sqrt( 1 + p.dataCCT$t[i-1]^2 )
    p   = dm / (dm - d0)
    
    mired   = (1-p)*p.dataCCT$mired[i-1]  +  p*p.dataCCT$mired[i]
        
    return( mired )
    }
    
    
#   CCT         temperature defining a Robertson isotherm, in K.   CCT=Inf is OK, mired is then 0.
#   locus.list  ufun, vfun, and miredInterval
#
#   find where the isotherm intersects the locus, using uniroot()
#   in general, we are looking at the spline locus between isotherm i and isotherm i+1; a sector.
#
#   return  a list with:
#           uv      the point of intersection on the locus
#           CCT     the native parameterization CCT in K
#           normal  unit vector along the isotherm, and approximate normal to locus
#                   this is used to translate the point uv a distance delta away from the locus
#   Or NULL in case of error.

nativeFromRobertson <- function( CCT, locus.list )
    {    
    mired   = 1.e6 / CCT
    
    j   = findInterval( mired, p.dataCCT$mired )
    
    if( j == 0 )
        {
        log.string( WARN, "CCT=%g is out of the Robertson LUT range. mired < %g", CCT, p.dataCCT$mired[1] )
        return(NULL)
        }
            
    n   = length( p.dataCCT$mired )
    if( j == n )
        {
        #   an extra check. 
        #   although I have found that 1.e6 / (1.e6 / 600) == 600 exactly on my PC, this might not be TRUE on other hardware,
        #   so use an epsilon just in case.
        epsilon = 1.e-12
        if( p.dataCCT$mired[n] + epsilon  <  mired )
            {
            log.string( WARN, "CCT=%g is outside the Robertson LUT range. %g < %g mired", 
                            CCT, p.dataCCT$mired[n], mired )
            return(NULL)
            }     

        #   otherwise set mired to what it probably should be (600),
        #   and keeping going with j == n, and the next special if block will take care of it !
        mired = p.dataCCT$mired[n] 
        }

    out = list()
    
    if( p.dataCCT$mired[j] == mired )    
        {
        #   special case, CCT is at one of the defining isotherms.
        #   The root is at the endpoint, so avoid using uniroot().
        #log.string( DEBUG, "mired=%g exactly. isotherm %d.  uniroot unnecessary.", mired, j )
        
        out$uv  = c( p.dataCCT$u[j], p.dataCCT$v[j] )
        out$CCT = CCT      
        
        tanj        = c( 1, p.dataCCT$t[j] )    # tangent vector to isotherm
        len         = sqrt( sum(tanj^2) )
        out$normal  = -tanj / len               # normal to locus (approx)
        
        return(out)
        }
        
    #   we are looking at the spline locus between isotherm j and isotherm j+1. a sector.
    #   See Fig. 2(3.11) in W&S. 
    
    #   compute fraction of the way between j and j+1
    alpha   =   (mired - p.dataCCT$mired[j]) / (p.dataCCT$mired[j+1] - p.dataCCT$mired[j])
        
    #   extract the 2 points on isotherms j and j+1
    pj      = c( p.dataCCT$u[j], p.dataCCT$v[j] )
    pj1     = c( p.dataCCT$u[j+1], p.dataCCT$v[j+1] )
      
    #   compute unit normals to the isotherms at both ends.  Maybe these should be moved into the LUT someday.
    #   these are very close to the unit tangents to the locus.
    tanj    = c( 1, p.dataCCT$t[j] )            # tangent vector to isotherm
    len     = sqrt( sum(tanj^2) )
    normj   = c( -tanj[2], tanj[1] ) / len      # unit normal to isotherm j
    
    tanj1   = c( 1, p.dataCCT$t[j+1] )          # tangent vector to isotherm
    len     = sqrt( sum(tanj1^2) )
    normj1  = c( -tanj1[2], tanj1[1] ) / len    # unit normal to isotherm j+1

    #   compute the normal for the intermediate isotherm
    norm    = (1-alpha)*normj  +  alpha*normj1
    #   norm    = norm / sqrt( sum(norm^2) )        do not make it a unit normal
    
    #   compute point p on the intermediate isotherm, on the segment [pj,pj1]
    delta   = pj1 - pj
    s       = alpha * sum( delta * normj1 ) / sum( delta * norm )
    p       = pj  +  s*delta
    
    #   from now on, mired refers to the parameter on locus.list, and not the Robertson mired isotherm.
    
    myfun   <- function( mired )
        {
        uv  = c( locus.list$ufun( mired ), locus.list$vfun( mired ) )
        
        return( sum( (uv - p)*norm ) )
        }
    
    miredInterval   = locus.list$miredInterval         #c( p.dataCCT$mired[i], p.dataCCT$mired[i+1] )

    #   check endpoint values
    f1  = myfun( miredInterval[1] )
    f2  = myfun( miredInterval[2] )
    
    if( 0 < f1 * f2 )    
        {
        #   same sign, should not happen
        log.string( WARN, "CCT=%g  mired=%g  p=%g,%g.  norm=%g,%g. Test function has the same sign at endpoints [%g,%g]. %g and %g CCT cannot be calculated.", 
                                CCT, mired, p[1], p[2], norm[1], norm[2], miredInterval[1], miredInterval[2], f1, f2  )
        return( NULL )
        }
        
    if( f1 == 0 )
        mired.end = miredInterval[1]
    else if( f2 == 0 )
        mired.end = miredInterval[2]    
    else
        {
        #log.string( DEBUG, "stats::uniroot() on interval [%g,%g], isotherms=%d and %d.  endpoint values %g and %g.", 
        #                        rangeMired[1], rangeMired[2], i, i+1, myfun(rangeMired[1],CCT), myfun(rangeMired[2],CCT)  )

        res = try( stats::uniroot( myfun, interval=miredInterval ),  silent=FALSE )
        
        if( inherits(res,"try-error" ) )    # class(res) == "try-error" )    
            {
            cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
            return( NULL )
            }     

        mired.end = res$root
        }
        
    #log.string( DEBUG, "stats::uniroot() successful after %d iterations.  isotherm %d", res$iter, i )

    out$uv      = c( locus.list$ufun( mired.end ), locus.list$vfun( mired.end ) )
    out$CCT     = 1.e6 / mired.end
    norm        = norm / sqrt( sum(norm^2) )    # make intermediate norm a unit vector     
    out$normal  = c( -norm[2], norm[1] )        # normal / len      #; print( out$normal )
            
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
    
    w   = topbot[1]/topbot[2]
    
    out = ((-449*w + 3525)*w - 6823.3)*w  +  5520.33
    
    if( out <= 0 )  return( NA_real_ )  # this *can* happen if w is too big
    
    return( out )
    }
    
CCTfromxy_McCamy  <- function( xy, locus.list, strict )
    {
    if( any( is.na(xy) ) )  return( NA_real_ )
    
    CCT = CCTfromxy_McCamy_nocheck( xy )
    
    if( is.finite(CCT)  &&  strict )
        {
        #   find intersection of McCamy isotherm and the locus
        res = nativeFromMcCamy( CCT, locus.list )
        
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
    
    
#   CCT         temperature defining a McCamy isotherm, in K
#   locus.list  ufun, vfun, and miredInterval
#
#   find where the isotherm intersects the locus, using uniroot()
#
#   return  a list with:
#           uv      the point of intersection on the locus
#           CCT     the native parameterization CCT in K
#           normal  to the locus, always along the McCamy isotherm for CCT
#       or NULL in case of error.
nativeFromMcCamy <- function( CCT, locus.list )
    {
    if( 34530 < CCT )   # the McCamy isotherm just below CCT=Inf.  34539
        return(NULL)
        
    if( CCT < 1621 )    # the cubic polynomial 'turns around' at about 1620.245
        return(NULL)
        
    ifun    <- function( w )    { ((-449*w + 3525)*w - 6823.3)*w + 5520.33 - CCT }
        
    #   the following interval comes from playing around with the cubic
    res = try( stats::uniroot( ifun, interval=c(-1.91,1.28), tol=.Machine$double.eps^0.33 ), silent=FALSE )
    
    if( inherits(res,"try-error" ) )    # class(res) == "try-error" )    
        {
        cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
        return( NULL )
        }        
    # log.string( DEBUG, "stats::uniroot() inverted McCamy cubic successful after %d iterations.", res$iter )
    
    #   find equation of the McCamy isotherm for this CCT
    alpha   = res$root        
    meet    = c(0.3320,0.1858)          #  the (x,y) 1931 where all the isotherms meet    
    C       = sum( c(1,-alpha) * meet )
    normal  = c( 1.5 - C, 4*C - alpha ) # normal to the isotherm
    uv0     = uvfromxy( meet, 1960 )    # all isotherms pass through this point.  uv0 = (0.2908709, 0.2441738)    
        
    myfun <- function( mired )
        {
        uv  = c( locus.list$ufun(mired), locus.list$vfun(mired) )
        
        return( sum( (uv-uv0)*normal ) )              # or ( sum(uv*normal) - 2*C )
        
        #xy  = c( 1.5*uv[1], uv[2] ) / ( uv[1] - 4*uv[2] + 2 )    
        #return( CCTfromxy_McCamy_nocheck(xy) - CCT )
        }
        
    miredInterval    = locus.list$miredInterval    # range( p.dataCCT$mired )
        
    #   log.string( DEBUG, "myfun() at endpoints: %g and %g.", myfun(miredInterval[1],CCT), myfun(miredInterval[2],CCT) )
    
    #   check endpoint values
    f1  = myfun( miredInterval[1] )
    f2  = myfun( miredInterval[2] )
    
    if( 0 < f1 * f2 )    
        {
        #   same sign, should not happen
        log.string( WARN, "CCT=%g.  Test function has the same sign at endpoints [%g,%g]. %g and %g. Intersection of isotherm and locus cannot be calculated.", 
                                CCT,  miredInterval[1], miredInterval[2], f1, f2  )
        return( NULL )
        }
        
    if( f1 == 0 )
        mired.end = miredInterval[1]
    else if( f2 == 0 )
        mired.end = miredInterval[2]    
    else
        {    
        res = try( stats::uniroot( myfun, interval=miredInterval, tol=.Machine$double.eps^0.33 ),  silent=FALSE )
        
        if( inherits(res,"try-error" ) )    # class(res) == "try-error" )    
            {
            cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
            return( NULL )
            }        
        # log.string( DEBUG, "stats::uniroot() found mired.end after %d iterations, for McCamy.", res$iter )
        
        mired.end   = res$root
        }
        
    uv  = c( locus.list$ufun( mired.end ), locus.list$vfun( mired.end ) )
    
    out = list()
    out$uv  = uv
    out$CCT = 1.e6 / mired.end
    out$normal  = c(-normal[2],normal[1]) / sqrt( sum(normal^2) )          #  approximate normal to the locus
    
    return(out)
    }
    
    
    
    
    
############        planckLocus()     ###############################################################

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
    #   process temperature
    ok  = is.numeric(temperature)  &&  0<length(temperature)
    if( ! ok )
        {
        log.string( ERROR, "Argument temperature is invalid.  It must be a numeric vector with positive length."  )
        return(NULL)
        }

    #   process locus
    ok  = is.character(locus)  &&  length(locus)==1
    if( ! ok )
        {
        log.string( ERROR, "Argument locus is invalid.  It must be a character vector with length 1." )
        return(NULL)
        }

    locusfull = c( 'Robertson', 'precision' )        
               
    idx.locus   = pmatch( tolower(locus), tolower(locusfull), nomatch=0 )
    if( idx.locus == 0 )
        {
        log.string( ERROR, "locus='%s' is invalid.", as.character(locus) )
        return( NULL )
        }    

    if( idx.locus == 1 )
        locus.list  = p.uvCubicsfromMired
    else
        locus.list  = p.uvQuinticsfromMired

        
    #   process param
    if( length(param) != 1 )
        {
        log.string( ERROR, "param is invalid, because it has length=%d != 1.", length(param)  )
        return( NULL )
        }        
        
    param   = as.character(param)
    param[ is.na(param) ] = 'native'
        
    paramfull = c( 'native', 'Robertson', 'McCamy' )        
        
    idx.param  = pmatch( tolower(param), tolower(paramfull), nomatch=0 )
    if( idx.param == 0 )
        {
        log.string( ERROR, "param='%s' is invalid.", param  )
        return( NULL )
        }
    
    #   process delta
    n   = length(temperature)
    
    ok  = is.numeric(delta)  &&  length(delta) %in% c(1,n)
    if( ! ok )
        {
        log.string( ERROR, "Argument delta is invalid.  It must be a numeric vector with length 1 or %d.", n  )
        return(NULL)
        }
    if( length(delta) == 1 )    delta = rep( delta[1], n )

    #   process space    
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
        
    if( FALSE )     # now let temperatures out of range be taken care of by lower functions
        {
        # amazingly, 600 seems to involute perfectly with 1.e6
        temperaturerange    = range( 1.e6 / locus.list$miredInterval )        # 33334  100000 )
        
        ok  =   temperaturerange[1]<=temperature  &  temperature<=temperaturerange[2]
        temperature[ ! ok ] = NA_real_
        }
    
    if( idx.param == 1 )
        {
        #   the easy one - plain 'native' 
        mired   = 1.e6 / temperature            # temperature==Inf is OK here, mired is then 0
        uv[ ,1]   = locus.list$ufun( mired )
        uv[ ,2]   = locus.list$vfun( mired )
        
        if( any( delta!=0 ) )
            {
            #   handle the delta offset
            for( i in 1:n )
                {
                if( delta[i] == 0  ||  is.na(mired[i]) ) next
                
                #   compute unit normal to the locus at point i
                normal  = c( -locus.list$vfun( mired[i], deriv=1 ), locus.list$ufun( mired[i], deriv=1 ) )
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
            #    uv[ ,1]   = p.uvCubicsfromMired[[1]]( 0 )
            #    uv[ ,2]   = p.uvCubicsfromMired[[2]]( 0 )
            #    next
            #    }
            
            if( idx.param == 3 )
                {
                res = nativeFromMcCamy( temperature[i], locus.list )
                }
            else
                {
                #   has to be Robertson
                res = nativeFromRobertson( temperature[i], locus.list )     # also works if temperature[i] is Inf
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
        
    #   now handle the output space
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
    
