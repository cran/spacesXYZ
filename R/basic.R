

    
################      XYZ  <-->  xyY    #######################

xyYfromXYZ <- function( XYZ )
    {
    XYZ = prepareNxM(XYZ,3)
    if( is.null(XYZ) )  return(NULL)

    xyY = cbind(NA_real_, NA_real_, XYZ[ ,2])
    rownames(xyY) = rownames(XYZ)
    colnames(xyY) = c('x','y','Y')

    #   .rowSums is faster than rowSums()
    denom   = .rowSums( XYZ, nrow(XYZ), 3 )

    #  denom0 is a logical vector, which may have logical NAs
    denom0  = denom == 0   

    #   denom0 may have logical NAs, but if so the corresponding assignment has no effect at that index
    #   which means denom[] was already NA_real_
    denom[ denom0 ] = NA_real_    
    
    #   NA_real values in denom will propagate into x and y
    #   in this matrix division, we are using the fact that denom is duplicated into the y column
    xyY[ ,1:2]  = XYZ[ ,1:2] / denom 

    warn    = denom0
    if( any(warn,na.rm=TRUE) )
        {
        #   but do not warn if X=Y=Z=0.  This is the pure-black special case.
        warn = warn  &  0<.rowSums( XYZ*XYZ, nrow(XYZ), 3 )
        if( any(warn,na.rm=TRUE) )
            log.string( WARN, "%d of %d XYZ vectors could not be transformed, because X+Y+Z==0.",    
                                    sum(warn,na.rm=TRUE), length(warn) )
        }
        
    return( xyY )
    }

XYZfromxyY <- function( xyY )
    {
    xyY = prepareNxM(xyY)
    if( is.null(xyY) )  return(NULL)

    XYZ <- cbind( NA_real_, xyY[,3], NA_real_)
    rownames(XYZ) = rownames(xyY)
    colnames(XYZ) = c('X','Y','Z')
    
    y   = xyY[ ,2]
    y0  = (y == 0)

    #   y0 is a logical vector that may contain NAs
    #   if y==0 then y0 is TRUE and this NA_real_ will propagate into mult and then into X and Z.  
    #   If y0 is NA, then this has no effect, but y is NA there anyway.
    y[ y0 ] = NA_real_        
    
    mult    = xyY[ ,3] / y
    XYZ[ ,1] = mult * xyY[ ,1]
    XYZ[ ,3] = mult * (1 - xyY[ ,1] - y)
        
    # treat Y=0 as a special case, set all XYZ=0
    #   xy can be anything, even NA, we don't care.
    Y0   = xyY[ ,3] == 0
    XYZ[Y0, ]  = 0    
            
    warn    = is.finite(xyY[ ,1]) &  y0  &  ! Y0
    if( any(warn,na.rm=TRUE) )
        log.string( WARN, "%d of %d xyY vectors could not be transformed because y==0.", 
                    sum(warn,na.rm=TRUE), length(warn) )

    return( XYZ )
    }


################      Lightness  <-->  Y    #######################

Lightness_from_linear  <-  function( Y )
    {
    thresh = (24/116)^3
    
    return( ifelse( Y < thresh, (116/12)^3 *Y, 116*Y^(1/3) - 16 ) )
    }
    
linear_from_Lightness  <-  function( L )
    {
    return( ifelse( L < 8, (12/116)^3 * L, ((L + 16)/116)^3 ) )
    }
    
    
    
    
################      XYZ  <-->  Lab    #######################

    
LabfromXYZ <- function( XYZ, white )
    {
    XYZ = prepareNxM(XYZ)
    if( is.null(XYZ) )  return(NULL)
    
    if( is.character(white) )
        white = standardXYZ( white[1] )

    ok  = is.numeric(white)  &&  length(white)==3  &&  all( 0 < white )
    if( ! ok )
        {
        log.string( ERROR, "white='%s' is invalid.", as.character(white) )
        return(NULL)
        }            

    L   = Lightness_from_linear( XYZ[ ,2]/white[2] )
    a   = (500/116) * (Lightness_from_linear( XYZ[ ,1]/white[1] ) - L)
    b   = (200/116) * (L - Lightness_from_linear( XYZ[ ,3]/white[3] ))

    #   special treatment for near neutrals
    #   if they are this close to gray, then snap to exact gray
    tol = 5.e-12
    neutral = (abs(a) + abs(b) < tol)
    neutral[ is.na(neutral) ]   = FALSE
    if( any(neutral) )
        {
        a[neutral] = 0
        b[neutral] = 0
        }

    Lab = cbind(L, a, b)

    rownames(Lab) = rownames(XYZ)
    colnames(Lab) = c('L','a','b')
    
    Lab
    }
    
    
XYZfromLab    <- function( Lab, white )
    {
    Lab = prepareNxM(Lab)
    if( is.null(Lab) )  return(NULL)
    
    if( is.character(white) )
        white = standardXYZ( white[1] )

    ok  = is.numeric(white)  &&  length(white)==3  &&  all( 0 < white )
    if( ! ok )
        {
        log.string( ERROR, "white='%s' is invalid.", as.character(white) )
        return(NULL)
        }            
    
    X   = white[1] * linear_from_Lightness( (116/500) * Lab[ ,2]  +  Lab[ ,1])
    Y   = white[2] * linear_from_Lightness( Lab[ ,1] ) 
    Z   = white[3] * linear_from_Lightness(-(116/200) * Lab[ ,3]  +  Lab[ ,1])

    XYZ = cbind(X,Y,Z)
    rownames(XYZ) = rownames(Lab)
    colnames(XYZ) = c('X','Y','Z')

    #   treat L==0 as a special case, pure black
    mask    = Lab[ ,1] == 0
    XYZ[mask, ] = 0

    XYZ
    }

    
    
################      XYZ  <-->  Luv    #######################
    
LuvfromXYZ <- function( XYZ, white )
    {
    #   Charles Poynton. Digital Video and HD. p. 281.
    XYZ = prepareNxM(XYZ)
    if( is.null(XYZ) )  return(NULL)
    

    if( is.character(white) )
        white = standardXYZ( white[1] )

    ok  = is.numeric(white)  &&  length(white)==3  &&  all( 0 < white )
    if( ! ok )
        {
        log.string( ERROR, "white='%s' is invalid.", as.character(white) )
        return(NULL)
        }            
        
    #   chromaticity of XYZs  (1976 UCS)
    uv  = uvfromXYZ( XYZ, space=1976 )

    #   chromaticity of white
    uv.ref  = uvfromXYZ( white, space=1976 )

    L = Lightness_from_linear( XYZ[ ,2]/white[2] )
    u = 13 * L * ( uv[ ,1] - uv.ref[1] )
    v = 13 * L * ( uv[ ,2] - uv.ref[2] )
    
    #   special treatment for near neutrals
    #   if they are this close to gray, then snap to exact gray
    tol = 5.e-12
    neutral = (abs(u) + abs(v) < tol)
    neutral[ is.na(neutral) ]   = FALSE
    if( any(neutral) )
        {
        u[neutral] = 0
        v[neutral] = 0
        }    

    Luv = cbind(L,u,v)
    rownames(Luv) = rownames(uv)
    colnames(Luv) = c('L','u','v')

    Luv
    }
    
    
XYZfromLuv <- function( Luv, white )
    {
    #   Charles Poynton. Digital Video and HD. p. 281.
    Luv = prepareNxM(Luv)
    if( is.null(Luv) )  return(NULL)

    #   verify white
    if( is.character(white) )
        white = standardXYZ( white[1] )

    ok  = is.numeric(white)  &&  length(white)==3  &&  all( 0 < white )
    if( ! ok )
        {
        log.string( ERROR, "white='%s' is invalid.", as.character(white) )
        return(NULL)
        }            
        
    #   chromaticity of white
    uv.ref  = uvfromXYZ( white, space=1976 )        
        
    u   = Luv[,2] / ( 13 * Luv[ ,1] )  +  uv.ref[1]
    v   = Luv[,3] / ( 13 * Luv[ ,1] )  +  uv.ref[2]

    Y   = white[2] * linear_from_Lightness( Luv[ ,1] ) 
    X   = Y * (9*u)/(4*v)
    Z   = Y * (12 - 3*u - 20*v)/(4*v)

    XYZ = cbind(X,Y,Z)
    rownames(XYZ) = rownames(Luv)
    colnames(XYZ) = c('X','Y','Z')

    #   treat L==0 as a special case, pure black. 
    #   Without this check, X and Z would be NaN
    mask    = Luv[ ,1] == 0
    XYZ[mask, ] = 0

    XYZ
    }    
        
    
    
################      uv  <---  XYZ    #######################

uvfromXYZ <- function( XYZ, space=1976 )
    {
    XYZ = prepareNxM(XYZ)
    if( is.null(XYZ) )  return(NULL)
    
    if( ! match(space,c(1960,1976),nomatch=FALSE) )
        {
        log.string( ERROR, "space='%s' is invalid.",  as.character(space[1]) )
        return(NULL)
        }
    
        
    denom = XYZ[ ,1]  +  15 * XYZ[ ,2]  +  3 * XYZ[ ,3]
    
    mask    =   denom <= 0
    
    if( any(mask,na.rm=T) )
        {
        log.string( WARN, "%d of %d XYZ vectors could not be transformed, because X + 15Y + 3Z <= 0.",
                                sum(mask,na.rm=T), length(mask) )
        denom[ mask ] = NA_real_                                
        }

    if( space == 1976 )
        {
        uv  = cbind( 4*XYZ[ ,1], 9*XYZ[ ,2] ) / denom
        colnames(uv)   = c("u'","v'")
        }
    else if( space == 1960 )
        {
        uv = cbind( 4*XYZ[ ,1], 6*XYZ[ ,2] ) / denom
        colnames(uv)   = c('u','v')        
        }
        
    rownames(uv)   = rownames(XYZ)
        
    return( uv )
    }    
    
    
uvfromxy <- function( xy, space=1976 )
    {
    xy = prepareNxM( xy, M=2 )
    if( is.null(xy) )  return(NULL)
    
    if( ! match(space,c(1960,1976,1931),nomatch=FALSE) )
        {
        log.string( ERROR, "space='%s' is invalid.",  as.character(space[1]) )
        return(NULL)
        }

    denom = -2 * xy[ ,1]  +  12 * xy[ ,2]  +  3
    
    mask    = denom <= 0
    
    if( any(mask,na.rm=T) )
        {
        log.string( WARN, "%d of %d xy vectors could not be transformed, because -2x + 12y + 3 <= 0.",
                                sum(mask,na.rm=T), length(mask) )
        denom[ mask ] = NA_real_                                
        }

    if( space == 1976 )
        {
        uv  = cbind( 4*xy[ ,1], 9*xy[ ,2] ) / denom
        colnames(uv)   = c("u'","v'")
        }
    else if( space == 1960 )
        {
        uv = cbind( 4*xy[ ,1], 6*xy[ ,2] ) / denom
        colnames(uv)   = c('u','v')        
        }
    else if( space == 1931 )
        {
        #   undocumented 'feature'
        colnames(xy) = c('x','y')
        return(xy)
        }
        
    rownames(uv)   = rownames(xy)
        
    return( uv )
    }    
    

####   requires private data frame p.dataIlluminants, which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()  ##

illumsubset <- function( name )
    {
    if( is.null(name) || (is.character(name) && 0<length(name) && name[1]=='*') )
        # just return everything !
        idx = 1:nrow(p.dataIlluminants)
    else
        #   this idx might have length = 0, but that's OK
        idx = pmatch( toupper(name), rownames(p.dataIlluminants) )
        
    return(idx)
    }
    
#   pmatch() assumes that the rownames of p.dataIlluminants are upper case
    
standardXYZ <- function( name )
    {
    idx = illumsubset( name )

    return( p.dataIlluminants$XYZ[ idx,  ,drop=F ] )
    }
    
standardxy <- function( name )
    {
    idx = illumsubset( name )

    return( p.dataIlluminants$xy[ idx,  ,drop=F ] )
    }
    
    
    