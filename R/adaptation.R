


#   constructor for CAT - Chromatic Adaptation Transform
#
CAT <- function( source.XYZ, target.XYZ, method="Bradford" )
    {
    if( is.character(source.XYZ) )
        source.XYZ = standardXYZ( source.XYZ[1] )
    
    ok  = is.numeric(source.XYZ)  &&  length(source.XYZ)==3  &&  all( is.finite(source.XYZ) )
    if( ! ok )
        {
        log.string( ERROR, "source.XYZ='%s' is invalid.", as.character(source.XYZ[1]) )
        return(NULL)
        }    
        
    if( is.character(target.XYZ) )
        target.XYZ = standardXYZ( target.XYZ[1] )
    
    ok  = is.numeric(target.XYZ)  &&  length(target.XYZ)==3  &&  all( is.finite(target.XYZ) )
    if( ! ok )
        {
        log.string( ERROR, "target.XYZ='%s' is invalid.", as.character(target.XYZ[1]) )
        return(NULL)
        }   
    
    
    full    = names(p.Ma)   # c( "Bradford", "VonKries", "MCAT02", "scaling" )
    
    idx     = pmatch( tolower(method[1]), tolower(full) )
    if( is.na(idx) )
        {
        log.string( ERROR, "method='%s' is invalid.", method )
        return(NULL)
        }    

    Ma = p.Ma[[ idx ]]

    #rownames(Ma)    = c('L','M','S')
    #colnames(Ma)    = c('X','Y','Z')
    
    
    if( identical( as.numeric(source.XYZ), as.numeric(target.XYZ) ) )
        #   trivial case
        M   = diag(3)
    else
        {
        crm_src = coneResponseMatrix( Ma, source.XYZ )
        crm_tgt = coneResponseMatrix( Ma, target.XYZ )
        
        if( is.null(crm_src)  ||  is.null(crm_tgt) )   return(NULL)
        
        M   = solve(crm_tgt) %*% crm_src
        }

    rownames(M) = c('X','Y','Z')
    colnames(M) = c('X','Y','Z')

    
    out = list()
    
    out$method  = full[idx]
    out$Ma      = Ma
    out$source.XYZ  = source.XYZ
    out$source.xyY  = specialxyY( source.XYZ )
    
    out$target.XYZ  = target.XYZ    
    out$target.xyY  = specialxyY( target.XYZ )
    
    out$M       = M

    class(out)  = c( "CAT", class(out) )
    
    return( out )
    }
    

#   x           the CAT
#   XYZ.src     XYZ in the source viewing environment, a 3-vector or an Nx3 matrix
#
#   returns     XYZ adapted to target white
adaptXYZ.CAT  <-  function( x, XYZ.src )
    {
    XYZ.src = prepareNxM( XYZ.src )    
    if( is.null(XYZ.src) )  return(NULL)

    out = tcrossprod( XYZ.src, x$M )
    
    #rownames(out)   = rownames(XYZ.src)
    #colnames(out)   = c('X','Y','Z')
    
    return( out )
    }


#   xyY.src     xyY in the source viewing environment, a 3-vector or an nx3 matrix
#
#   returns  xyY adapted to white.tgt
adaptxyY.CAT  <-  function( x, xyY.src )
    {
    xyY.src = prepareNxM( xyY.src )    
    if( is.null(xyY.src) )  return(NULL)
    
    #   convert to XYZ    
    XYZ.src = XYZfromxyY( xyY.src )
    if( is.null(XYZ.src) ) return(NULL)
    
    #   adapt the XYZ
    XYZ.tgt    = adaptXYZ( x, XYZ.src )
    if( is.null(XYZ.tgt) )  return(NULL)
    
    #   now back to xyY
    xyY.tgt     = xyYfromXYZ( XYZ.tgt )
    
    #   do an additional test for neutrals
    #   to correct possible roundoff error
    neutral = xyY.src[ ,1] == x$source.xyY[1]   &    xyY.src[ ,2] == x$source.xyY[2] 
    
    neutral[ is.na(neutral) ]   = FALSE    
    #   mask    = (xyY.src[ ,1] == white.src[1])  &  (xyY.src[ ,2] == white.src[2])
    
    if( any(neutral) )
        {
        #   for some rows, xyY.src has the exact chromaticity of source.XYZ
        #   therefore, xyY.tgt should have exact chromaticity of target.XYZ
        xyY.tgt[neutral,1] = x$target.xyY[1]
        xyY.tgt[neutral,2] = x$target.xyY[2]
        
        #   log.string( DEBUG, "%d neutrals (of %d) given special treatment.",   sum(neutral), length(neutral) )
        }
    
    return( xyY.tgt )
    }
    
    
    
#   Lab.src     Lab in the source viewing environment, a 3-vector or an Nx3 matrix
#
#   returns  Lab adapted to white.tgt
adaptLab.CAT  <-  function( x, Lab.src )
    {
    Lab.src = prepareNxM( Lab.src )    
    if( is.null(Lab.src) )  return(NULL)
    
    #   convert to XYZ    
    XYZ.src = XYZfromLab( Lab.src, x$source.XYZ )
    if( is.null(XYZ.src) ) return(NULL)
    
    #   adapt the XYZ
    XYZ.tgt    = adaptXYZ( x, XYZ.src )
    if( is.null(XYZ.tgt) )  return(NULL)
    
    #   now back to Lab
    Lab.tgt     = LabfromXYZ( XYZ.tgt, x$target.XYZ )
    
    return( Lab.tgt )
    }
    
    
#   Luv.src     Luv in the source viewing environment, a 3-vector or an Nx3 matrix
#
#   returns  Luv adapted to white.tgt
adaptLuv.CAT  <-  function( x, Luv.src )
    {
    Luv.src = prepareNxM( Luv.src )    
    if( is.null(Luv.src) )  return(NULL)
    
    #   convert to XYZ    
    XYZ.src = XYZfromLuv( Luv.src, x$source.XYZ )
    if( is.null(XYZ.src) ) return(NULL)
    
    #   adapt the XYZ
    XYZ.tgt    = adaptXYZ( x, XYZ.src )
    if( is.null(XYZ.tgt) )  return(NULL)
    
    #   now back to Luv
    Luv.tgt     = LuvfromXYZ( XYZ.tgt, x$target.XYZ )
    
    return( Luv.tgt )
    }
    
    
    
    
    


#   .Ma     the matrix for the method
#   .white  the reference white    
#   returns a matrix that maps .white to (1,1,1)
coneResponseMatrix  <-  function( .Ma, .white )
    {   
    lms = .Ma %*% as.numeric(.white)
    
    if( any( lms<=0 ) )
        {
        log.string( ERROR, "white XYZ=(%g,%g,%g) is too far from neutral.", 
                        .white[1], .white[2], .white[3] )
        return(NULL)
        }
    
    #  return( diag( as.double(1/lms) ) %*% .Ma )
    
    return( as.double(1/lms)  * .Ma )   # a trick which is equivalent, but a little more efficient  
    }
        
        
specialxyY <- function( XYZ, digits=5 )
    {
    xyY  = xyYfromXYZ(XYZ)
    
    xy.rounded  = round(xyY[1:2],digits)
    
    tol = 5.e-12
    
    if( sum( abs( xy.rounded - xyY[1:2] ) ) < tol )
        # XYZ must have been computed from a standard illuminant
        # it cannot be a coincidence
        xyY[1:2] = xy.rounded
        
    xyY        
    }

    
        
#--------       UseMethod() calls           --------------#    
    
adaptXYZ <- function( x, XYZ.src ) 
    {
    UseMethod("adaptXYZ")
    }            
    
adaptxyY <- function( x, xyY.src ) 
    {
    UseMethod("adaptxyY")
    }                
    
adaptLab <- function( x, Lab.src ) 
    {
    UseMethod("adaptLab")
    }                
    
adaptLuv <- function( x, Luv.src ) 
    {
    UseMethod("adaptLuv")
    }                    