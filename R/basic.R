

    
################      XYZ  <-->  xyY    #######################

xyYfromXYZ <- function( XYZ )
    {
    XYZ = prepareNxM(XYZ)
    if( is.null(XYZ) )  return(NULL)

    xyY = cbind(NA_real_, NA_real_, XYZ[ ,2])
    rownames(xyY) = rownames(XYZ)
    colnames(xyY) = c('x','y','Y')

    denom   = rowSums( XYZ )
    mask    = 0<denom  &  0<=XYZ[ ,2]

    xyY[mask,1]    = XYZ[mask,1] / denom[mask]    # (XYZ[w,1]+XYZ[w,2]+XYZ[w,3])
    xyY[mask,2]    = XYZ[mask,2] / denom[mask]    # (XYZ[w,1]+XYZ[w,2]+XYZ[w,3])
    
    if( sum(mask) < length(mask) )
        log.string( WARN, "%d of %d XYZ vectors could not be transformed.", 
                    length(mask)-sum(mask), length(mask) )

    xyY
    }

XYZfromxyY <- function( xyY )
    {
    xyY = prepareNxM(xyY)
    if( is.null(xyY) )  return(NULL)

    XYZ <- cbind( NA_real_, xyY[,3], NA_real_)
    rownames(XYZ) = rownames(xyY)
    colnames(XYZ) = c('X','Y','Z')
    
    #   treat Y=0 as a special case - pure black
    mask0   = xyY[ ,3] == 0
    XYZ[mask0, ]  = 0    

    mask    = xyY[,2] != 0  &  ! mask0
    mask[ is.na(mask) ] = FALSE
    
    if( any(mask) )
        {
        xyY_sub     = xyY[mask, ,drop=FALSE]
        mult        = xyY_sub[ ,3] / xyY_sub[ ,2]
        XYZ[mask,1] = mult * xyY_sub[ ,1]
        XYZ[mask,3] = mult * (1-xyY_sub[ ,1]-xyY_sub[ ,2])
        }
        
    #   mask[ mask0 ] = TRUE

    if( sum(mask | mask0) < length(mask) )
        log.string( WARN, "%d of %d xyY vectors could not be transformed.", 
                    length(mask)-sum(mask|mask0), length(mask) )

    XYZ
    }


################      Lightness  <-->  Y    #######################

Lightness_from_linear  <-  function( Y )
    {
    thresh = (24/116)^3
    
    out = ifelse( Y < thresh, (116/12)^3 *Y, 116*Y^(1/3) - 16 )
    
    out
    }
    
linear_from_Lightness  <-  function( L )
    {
    out = ifelse( L < 8, (12/116)^3 * L, ((L + 16)/116)^3 )
    
    out
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
    uv  = uvfromXYZ( XYZ, version=1976 )

    #   chromaticity of white
    uv.ref  = uvfromXYZ( white, version=1976 )

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
    uv.ref  = uvfromXYZ( white, version=1976 )        
        
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

uvfromXYZ <- function( XYZ, version=1976 )
    {
    XYZ = prepareNxM(XYZ)
    if( is.null(XYZ) )  return(NULL)
        
    denom = XYZ[ ,1]  +  15 * XYZ[ ,2]  +  3 * XYZ[ ,3]
    
    mask    = denom <= 0
    
    if( any(mask) )
        {
        denom[ mask ] = NA_real_
        log.string( WARN, "%d of %d XYZ vectors could not be transformed, because X + 15Y + 3Z <= 0.",
                                sum(mask), length(mask) )
        }

    if( version == 1976 )
        {
        uv  = cbind( 4*XYZ[ ,1], 9*XYZ[ ,2] ) / denom
        colnames(uv)   = c("u'","v'")
        }
    else if( version == 1960 )
        {
        uv = cbind( 4*XYZ[ ,1], 6*XYZ[ ,2] ) / denom
        colnames(uv)   = c('u','v')        
        }
    else
        {
        log.string( WARN, "version='%s' is invalid.",  as.character(version[1]) )
        return(NULL)
        }
        
    rownames(uv)   = rownames(XYZ)
        
    return( uv )
    }    
    

#   requires private data frame p.dataIlluminants, which is lazy-loaded from sysdata.rda;  see savePrivateDatasets()
standardXYZ <- function( name )
    {
    #   XYZ             = p.dataIlluminants$XYZ 
    #   rownames(XYZ)   = p.dataIlluminants$Name
    
    idx = pmatch( toupper(name), rownames(p.dataIlluminants$XYZ) )

    return( p.dataIlluminants$XYZ[ idx,  ,drop=F ] )
    }
    

    
    