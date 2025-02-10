

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


#   uv          an Nx2 matrix, or a vector that can be converted to such a matrix
#   isotherms   a character vector of isotherms, there can be multiple isotherms
#   locus       a single character string

CCTfromuv  <- function( uv, isotherms='robertson', locus='robertson', strict=FALSE )
    {
    uv  = prepareNxM( uv, M=2 )
    if( is.null(uv) )   return(NULL)

    #   process isotherms
    if( length(isotherms) == 0 )
        {
        log_level( ERROR, "isotherms is invalid, because length(isotherms)=0."  )
        return( NULL )
        }

    isofull = c( 'native', 'Robertson', 'McCamy' )

    isotherms   =  as.character(isotherms)

    #   isotherms[ is.na(isotherms) ]   = 'native'  kind of stupid option

    idx.isotherms   = pmatch( tolower(isotherms), tolower(isofull), nomatch=0, duplicates.ok=TRUE )
    if( any( idx.isotherms==0 ) )
        {
        log_level( ERROR, "isotherms='%s' is invalid.", isotherms[idx.isotherms==0] )
        return( NULL )
        }

    #   process locus
    ok  = is.character(locus)  &&  length(locus)==1
    if( ! ok )
        {
        log_level( ERROR, "Argument locus is invalid.  It must be a character vector with length 1." )
        return(NULL)
        }

    locusfull = c( 'Robertson', 'precision' )

    idx.locus   = pmatch( tolower(locus), tolower(locusfull), nomatch=0 )
    if( idx.locus == 0 )
        {
        log_level( ERROR, "locus='%s' is invalid.", as.character(locus) )
        return( NULL )
        }

    if( idx.locus == 1 )
        #   the Robertson locus (cubic spline)
        locus.list  = p.uvCubicsfromMired
    else
        #   the precision locus (quintic spline)
        locus.list  = p.uvQuinticsfromMired


    n   = nrow(uv)

    out = matrix( NA_real_, n, length(idx.isotherms) )
    rownames(out)   = rownames(uv)
    colnames(out)   = isofull[ idx.isotherms ]

    Duv = matrix( NA_real_, n, length(idx.isotherms) )
    rownames(Duv)   = rownames(uv)
    colnames(Duv)   = isofull[ idx.isotherms ]

    #   assign matrix one element at a time
    
    for( i in 1:n )
        {
        for( j in 1:length(idx.isotherms) )
            {
            #   now on column j
            idx = idx.isotherms[j]

            if( idx == 1 )
                {
                #   native isotherms, from either Robertson locus or the precision locus
                #for( i in 1:n )
                #    out[i,j]  = CCTfromuv_native( uv[i, ], locus.list, strict  )
                res = fullCCTfromuv_native( uv[i, ], locus.list )
                }
            else if( idx == 2 )
                {
                #   Robertson isotherms
                res = fullCCTfromuv_Robertson( uv[i, ], locus.list )                
                }
            else if( idx == 3 )
                {
                #   McCamy isotherms
                res = fullCCTfromuv_McCamy( uv[i, ], locus.list )                
                }
                
            if( is.null(res) )  next

            out[i,j]    = res$CCT
            Duv[i,j]    = res$Duv
            }
        }

    if( strict )
        {
        #   if abs(Duv) is too large, or if Duv is NA,  change CCT to NA
        out[ 0.05<abs(Duv)  |  is.na(Duv) ]  = NA_real_
        }

    if( length(idx.isotherms) == 1 )
        {
        #   change nx1 matrix to just a plain vector, and copy rownames to names
        rnames      = rownames(out)
        dim(out)    = NULL
        names(out)  = rnames
        
        dim(Duv)    = NULL
        names(Duv)  = rnames
        }

    attr( out, "Duv" )  = Duv

    return( out )
    }





#   uv          point, usually off locus
#   locus.list  ufun, vfun, and miredInterval

#   computes the native isotherm of the given locus that passes through uv
#   It works by searching along the locus, while 'sweeping around' the normal line
#   to find the normal line passing through uv.  It uses stats::uniroot()

#   returns a list with these items:
#       uv0     the point on the locus
#       CCT     the locus parameter at uv0
#       normal  unit normal to the locus at uv0, upward pointing
#       Duv     signed distance from uv to uv0; it's along the locus normal


fullCCTfromuv_native  <- function( uv, locus.list )
    {
    if( any( is.na(uv) ) )  return(NULL)

    myfun <- function( mir )
        {
        #   mir         = 1.e6 / temp

        uv.locus    = c( locus.list$ufun( mir ), locus.list$vfun( mir ) )

        tangent     = c( locus.list$ufun( mir, deriv=1 ), locus.list$vfun( mir, deriv=1 ) )

        return( sum( tangent * (uv.locus - uv) ) )  # relevant dot product
        }

    miredInterval    = locus.list$miredInterval           # range( p.dataCCT$mired )

    #   check endpoint values
    f1  = myfun( miredInterval[1] )
    f2  = myfun( miredInterval[2] )

    if( 0 < f1 * f2 )
        {
        #   same sign, uv is invalid
        log_level( WARN, "uv=%g,%g is in invalid region.  CCT cannot be calculated.", uv[1], uv[2] )
        return( NULL )
        }

    if( f1 == 0 )
        mired.root  = miredInterval[1]
    else if( f2 == 0 )
        mired.root  = miredInterval[2]
    else
        {
        #   reduced the tolerance here, but still only takes about 8 iterations
        res = try( stats::uniroot( myfun, interval=miredInterval, tol=.Machine$double.eps^0.5 ),  silent=FALSE )

        if( inherits(res,"try-error" ) )    # class(res) == "try-error" )
            {
            cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
            return( NULL )
            }

        mired.root   = res$root

        #log_level( DEBUG, "stats::uniroot() successful after %d iterations.  mired.end=%g", res$iter, mired.end )
        }

    out = list()

    out$uv0 = c( locus.list$ufun( mired.root ),locus.list$vfun( mired.root ) )
    out$CCT = 1.e6 / mired.root

    tangent = c( locus.list$ufun( mired.root, deriv=1 ), locus.list$vfun( mired.root, deriv=1 ) )

    out$normal  = c( -tangent[2], tangent[1] ) / sqrt( sum(tangent^2) )

    offset  = uv - out$uv0

    out$Duv     = sum( offset * out$normal )

    #   to determine the sign, use offset[2] which is v = the vertical direction in uv plane
    #out$Duv = sign(offset[2]) * sqrt( sum(offset^2) )   # signed distance

    return( out )
    }



#   uv          a 2-vector with uv 1960
#   locus.list  both Robertson and precision are OK 
#
#   returns a list with
#       CCT     of the isotherm that passes through uv
#       Duv     signed distance of uv from locus, along the isotherm    #, or NA if strict is TRUE and too far from locus
#
#   requires private data frame p.dataCCT,  which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()

fullCCTfromuv_Robertson  <- function( uv, locus.list )
    {
    if( any( is.na(uv) ) )  return( NULL )

    data    = datafromuv_Robertson_nocheck( uv, extrap=FALSE )
    
    mired   = data$mired

    if( is.na(mired) )  return( NULL )

    CCT = 1.e6 / mired

    res = nativeFromRobertson( CCT, locus.list )

    if( is.null(res) )  return( NULL )

    offset  = uv - res$uv   # res$uv is on the locus

    # test    = sqrt( sum(offset*offset) )
    
    #   offset and data$tangent are parallel, 
    #   so their dot product is signed distance of offset along tangent
    Duv = sum( offset * data$tangent )
    
    if( FALSE )
        {
        #   sanity check, compare data$tangent and res$normal
        #   the difference should be very small
        dif = data$tangent - res$normal
        
        cat( "dif = ", dif, '\n' )
        }

    return( list( CCT=CCT, Duv=Duv ) )
    }



#   uv          a 2-vector with uv 1960
#   extrap      TRUE means extrapolate past mired=0 and 600.  FALSE means return NA_real
#   requires private data frame p.dataCCT,  which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()
#   'nocheck' means do not that input is NA and do not check locus
#
#   returns a list with:
#       mired   of the Robertson isotherm that passes through uv
#       tangent unit tangent vector along the isotherm, pointing UP
#
#   uniroot() is not used, and there is no locus involved

datafromuv_Robertson_nocheck  <- function( uv, extrap=FALSE )
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
        if( 2 <= length(idx) )
            {
            #   more than 1 should not happen
            log_level( WARN, "uv=%g,%g lies on more than 1 isotherm. Mired cannot be calculated.", uv[1], uv[2] )
            return( list( mired=NA_real_, tangent=NA_real_ ) )
            }

        #   exactly 1 is easy
        mired   = p.dataCCT$mired[idx]

        tangent = unitize( -c( 1, p.dataCCT$t[idx] ) )    #   take negative so it points UP

        return( list( mired=mired, tangent=tangent ) )
        }

    #   look for zero crossings
    idx = which( di[1:(n-1)] * di[2:n] < 0 )

    if( length(idx) == 0 )
        {
        #   all di must be pos or neg
        if( ! extrap )
            {
            log_level( WARN, "uv=%g,%g is in invalid region. Found no zero-crossings. Mired cannot be calculated.", uv[1], uv[2] )
            return( list( mired=NA_real_, tangent=NA_real_ ) )
            }

        #   do a very crude extrapolation, only useful for root-finding
        if( di[1] < 0 )
            return( p.dataCCT$mired[1] + 100 * di[1] )
        else
            return( p.dataCCT$mired[n] + 30 * di[n] )
        }
    else if( 2 <= length(idx) )        #     ||  i==1 ) ?
        {
        log_level( WARN, "uv=%g,%g is in invalid region. Found multiple zero-crossings. Mired cannot be calculated.", uv[1], uv[2] )

        return( list( mired=NA_real_, tangent=NA_real_ ) )
        }

    #   exactly 1 zero-crossing is OK.  This is the normal situation
    i = idx[1] + 1

    #   find the unit tangent vectors along both isotherms
    tan0    = -c( 1, p.dataCCT$t[i] )       #   take negative so it points UP
    tan1    = -c( 1, p.dataCCT$t[i-1] )     #   take negative so it points UP

    norm0   = sqrt( sum(tan0^2) )
    norm1   = sqrt( sum(tan1^2) )

    tan0    = tan0 / norm0
    tan1    = tan1 / norm1

    d0  = di[i]   / norm0
    d1  = di[i-1] / norm1
    rho = d1 / (d1 - d0)        #  interpolation parameter

    mired   = (1-rho)*p.dataCCT$mired[i-1]  +  rho*p.dataCCT$mired[i]
    tangent = (1-rho)*tan1  +  rho*tan0

    out = list()
    out$mired   = mired
    out$tangent = unitize( tangent )

    return( out )
    }


#   CCT         temperature defining a Robertson isotherm, in K.   CCT=Inf is OK, mired is then 0.
#   locus.list  ufun, vfun, and miredInterval
#
#   find where the given isotherm intersects the given locus, using uniroot()

#   in general, we are looking at the spline locus between isotherm i and isotherm i+1; a sector.
#
#   return  a list with:
#           uv      the point of intersection on the given locus
#           CCT     the native parameter of uv in K; not necessarily the same as the given CCT
#           normal  unit vector along the Robertson isotherm, pointing up, and approximately normal to locus
#                   this is used to translate the point uv a distance delta away from the locus
#   Or NULL in case of error.

nativeFromRobertson <- function( CCT, locus.list )
    {
    mired   = 1.e6 / CCT

    j   = findInterval( mired, p.dataCCT$mired )

    if( j == 0 )
        {
        log_level( WARN, "CCT=%g is out of the Robertson LUT range. mired < %g", CCT, p.dataCCT$mired[1] )
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
            log_level( WARN, "CCT=%g is outside the Robertson LUT range. %g < %g mired",
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
        #log_level( DEBUG, "mired=%g exactly. isotherm %d.  uniroot unnecessary.", mired, j )

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

    #   compute the normal for the intermediate isotherm at temperature CCT
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
        log_level( WARN, "CCT=%g  mired=%g  p=%g,%g.  norm=%g,%g. Test function has the same sign at endpoints [%g,%g]. %g and %g CCT cannot be calculated.",
                                CCT, mired, p[1], p[2], norm[1], norm[2], miredInterval[1], miredInterval[2], f1, f2  )
        return( NULL )
        }

    if( f1 == 0 )
        mired.root  = miredInterval[1]
    else if( f2 == 0 )
        mired.root  = miredInterval[2]
    else
        {
        #log_level( DEBUG, "stats::uniroot() on interval [%g,%g], isotherms=%d and %d.  endpoint values %g and %g.",
        #                        rangeMired[1], rangeMired[2], i, i+1, myfun(rangeMired[1],CCT), myfun(rangeMired[2],CCT)  )

        res = try( stats::uniroot( myfun, interval=miredInterval ),  silent=FALSE )

        if( inherits(res,"try-error" ) )    # class(res) == "try-error" )
            {
            cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
            return( NULL )
            }

        mired.root = res$root
        }

    #log_level( DEBUG, "stats::uniroot() successful after %d iterations.  isotherm %d", res$iter, i )

    out$uv      = c( locus.list$ufun( mired.root ), locus.list$vfun( mired.root ) )
    out$CCT     = 1.e6 / mired.root
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
#
#   returns CCT of the McCamy isotherm that passes through xy

#   returns a list with:
#       CCT     of the McCamy isotherm that passes through uv
#       tangent unit tangent vector along the isotherm, pointing UP

#   uniroot() is not used, and there is no locus involved

datafromuv_McCamy_nocheck  <- function( uv )
    {
    xy  = xyfromuv( uv )
    
    xy0     = c(0.3320,0.1858)      # the point where all the isotherms intersect, in xy plane
    
    topbot  = xy - xy0    

    if( topbot[2] <= 0 )  return( NULL )

    w   = topbot[1]/topbot[2]

    CCT = ((-449*w + 3525)*w - 6823.3)*w  +  5520.33

    if( CCT <= 0 )  return( NULL )  # this *can* happen if w is too big

    uv0 = uvfromxy( xy0, space=1960 )       # the point where all the isotherms intersect, in uv plane

    tangent = unitize( uv - uv0 )
    
    out = list( CCT=CCT, tangent=tangent )

    return( out )
    }


#   returns a list with
#       CCT     of the isotherm that passes through uv
#       Duv     signed distance of uv from locus, along the isotherm    #, or NA if strict is TRUE and too far from locus

fullCCTfromuv_McCamy  <- function( uv, locus.list )
    {
    if( any( is.na(uv) ) )  return( NULL )
    
    data    = datafromuv_McCamy_nocheck( uv )
    
    if( is.null(data) )  return( NULL )
    
    CCT = data$CCT
    
    #   find intersection of McCamy isotherm and the locus
    res = nativeFromMcCamy( CCT, locus.list )
       
    if( ! is.null(res) )
        {
        Duv = sum( (uv - res$uv) * data$tangent )
        
        if( FALSE )
            {
            #   sanity check, compare data$tangent and res$normal
            #   the difference should be very small
            dif = data$tangent - res$normal
            
            cat( "fullCCTfromuv_McCamy().  dif = ", dif, '\n' )
            }            
        }
    else
        Duv = NA_real_
        
    out = list( CCT=CCT, Duv=Duv )
    
    return( out )
    }


#   CCT         temperature defining a McCamy isotherm, in K
#   locus.list  ufun, vfun, and miredInterval
#
#   find where the given isotherm intersects the given locus, using uniroot()
#
#   return  a list with:
#           uv      the point of intersection of the given isotherm, and the given locus
#           CCT     the native parameter for uv in K, not necessarily the same as the given CCT
#           normal  unit vector along the McCamy isotherm, pointing up, approximately normal to the locus
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
    # log_level( DEBUG, "stats::uniroot() inverted McCamy cubic successful after %d iterations.", res$iter )

    #   find equation of the McCamy isotherm for this CCT
    alpha   = res$root
    xy0     = c(0.3320,0.1858)              #  the (x,y) 1931 where all the isotherms meet
    C       = sum( c(1,-alpha) * xy0 )
    normal  = c( 1.5 - C, 4*C - alpha )     # normal to the isotherm
    uv0     = uvfromxy( xy0, space=1960 )   # all isotherms pass through this point.  uv0 = (0.2908709, 0.2441738)

    myfun <- function( mired )
        {
        uv  = c( locus.list$ufun(mired), locus.list$vfun(mired) )

        return( sum( (uv-uv0)*normal ) )              # or ( sum(uv*normal) - 2*C )

        #xy  = c( 1.5*uv[1], uv[2] ) / ( uv[1] - 4*uv[2] + 2 )
        #return( CCTfromxy_McCamy_nocheck(xy) - CCT )
        }

    miredInterval    = locus.list$miredInterval    # range( p.dataCCT$mired )

    #   log_level( DEBUG, "myfun() at endpoints: %g and %g.", myfun(miredInterval[1],CCT), myfun(miredInterval[2],CCT) )

    #   check endpoint values
    f1  = myfun( miredInterval[1] )
    f2  = myfun( miredInterval[2] )

    if( 0 < f1 * f2 )
        {
        #   same sign, should not happen
        log_level( WARN, "CCT=%g.  Test function has the same sign at endpoints [%g,%g]. %g and %g. Intersection of isotherm and locus cannot be calculated.",
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
        # log_level( DEBUG, "stats::uniroot() found mired.end after %d iterations, for McCamy.", res$iter )

        mired.end   = res$root
        }

    uv  = c( locus.list$ufun( mired.end ), locus.list$vfun( mired.end ) )

    out = list()
    out$uv  = uv
    out$CCT = 1.e6 / mired.end
    out$normal  = c(-normal[2],normal[1]) / sqrt( sum(normal^2) )   #  unit tangent to isotherm, and approximate normal to the locus

    return(out)
    }





############        planckLocus()     ###############################################################

#   This function is exported, for use by the user.

#   temperature a N-vector of temperatures, in K.  Inf is OK
#   locus       'robertson' or 'precision'
#   param       'native, 'robertson', or 'mccamy'.  It determines the isotherm family.
#   Duv         offset from the locus, in top direction.  Always along the corresponding isotherm. Distance in uv-1960 plane.
#   space       year of chromaticity space
#
#   returns     an Nx2 matrix of uv's, or NULL in case of argument error
#
planckLocus  <- function( temperature, locus='robertson', param='robertson', Duv=0, space=1960 )
    {
    #   process temperature
    ok  = is.numeric(temperature)  &&  0<length(temperature)
    if( ! ok )
        {
        log_level( ERROR, "Argument temperature is invalid.  It must be a numeric vector with positive length."  )
        return(NULL)
        }

    #   process locus
    ok  = is.character(locus)  &&  length(locus)==1
    if( ! ok )
        {
        log_level( ERROR, "Argument locus is invalid.  It must be a character vector with length 1." )
        return(NULL)
        }

    locusfull = c( 'Robertson', 'precision' )

    idx.locus   = pmatch( tolower(locus), tolower(locusfull), nomatch=0 )
    if( idx.locus == 0 )
        {
        log_level( ERROR, "locus='%s' is invalid.", as.character(locus) )
        return( NULL )
        }

    if( idx.locus == 1 )
        locus.list  = p.uvCubicsfromMired
    else
        locus.list  = p.uvQuinticsfromMired


    #   process param
    if( length(param) != 1 )
        {
        log_level( ERROR, "param is invalid, because it has length=%d != 1.", length(param)  )
        return( NULL )
        }

    param   = as.character(param)
    
    # param[ is.na(param) ] = 'native'  this replacement is kind of stupid, so removed

    paramfull = c( 'native', 'Robertson', 'McCamy' )

    idx.param  = pmatch( tolower(param), tolower(paramfull), nomatch=0 )
    if( idx.param == 0 )
        {
        log_level( ERROR, "param='%s' is invalid.", param  )
        return( NULL )
        }

    #   process Duv
    n   = length(temperature)

    ok  = is.numeric(Duv)  &&  length(Duv) %in% c(1,n)
    if( ! ok )
        {
        log_level( ERROR, "Argument Duv is invalid.  It must be a numeric vector with length 1 or %d.", n  )
        return(NULL)
        }
    if( length(Duv) == 1 )    Duv = rep( Duv, n )

    #   process space
    if( ! match(space,c(1960,1976,1931),nomatch=FALSE) )
        {
        log_level( ERROR, "space='%s' is invalid.",  as.character(space[1]) )
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

        if( any( Duv!=0 ) )
            {
            #   handle the Duv offset
            for( i in 1:n )
                {
                if( Duv[i] == 0  ||  is.na(mired[i]) ) next

                #   compute unit normal to the locus at point i
                normal  = c( -locus.list$vfun( mired[i], deriv=1 ), locus.list$ufun( mired[i], deriv=1 ) )
                len     = sqrt( sum(normal^2) )
                if( len == 0 )  next    # should never happen
                normal  = normal / len

                uv[i, ] = uv[i, ] + Duv[i]*normal
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
                #   in the next line, we get:
                #       res$uv      where the given temperature[i] isotherm intersects the given locus
                #       res$CCT
                #       res$normal  unit vector along the McCamy isotherm
                res = nativeFromMcCamy( temperature[i], locus.list )
                }
            else
                {
                #   has to be Robertson
                #       res$uv      where the given temperature[i] isotherm intersects the given locus
                #       res$CCT     parameter of res$uv, and not nec. equal to temperature[l]
                #       res$normal  unit vector along the Robertson isotherm
                res = nativeFromRobertson( temperature[i], locus.list )     # also works if temperature[i] is Inf
                }

            if( is.null(res) ) next

            uv[i, ] = res$uv

            #   handle the Duv offset
            if( is.finite(Duv[i])  &&  Duv[i]!=0 )
                {
                uv[i, ] = uv[i, ] + Duv[i] * res$normal
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


########################        deadwood below      #################################


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

    #   log_level( TRACE, "uv=%g,%g", uv[1], uv[2] )
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
        log_level( WARN, "uv=%g,%g is in invalid region.  CCT cannot be calculated.", uv[1], uv[2] )
        return( NA_real_ )
        }

    if( f1 == 0 )
        mired.end = miredInterval[1]
    else if( f2 == 0 )
        mired.end = miredInterval[2]
    else
        {
        #log_level( DEBUG, "stats::uniroot() on interval [%g,%g],    endpoint values %g and %g.",
        #                        miredInterval[1], miredInterval[2],  f1, f2 )

        #   reduced the tolerance here, but still only takes about 8 iterations
        res = try( stats::uniroot( myfun, interval=miredInterval, uv=uv, tol=.Machine$double.eps^0.5 ),  silent=FALSE )

        if( inherits(res,"try-error" ) )    # class(res) == "try-error" )
            {
            cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
            return( NA_real_ )
            }

        mired.end   = res$root

        #log_level( DEBUG, "stats::uniroot() successful after %d iterations.  mired.end=%g", res$iter, mired.end )
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
            log_level( WARN, "uv=%g,%g is invalid, because its distance to the Planckian locus = %.7f > 0.05.  (mired=%g)",
                                uv[1], uv[2], dist, mired.end )
            return( NA_real_ )
            }
        }

    return( 1.e6 / mired.end )
    }



