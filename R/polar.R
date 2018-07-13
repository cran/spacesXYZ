
################      LCHab  <-->  Lab    #######################

LCHabfromLab    <- function( Lab )
    {
    Lab = prepareNxM(Lab)
    if( is.null(Lab) )  return(NULL)
    
    C   = sqrt( Lab[ ,2]^2 +  Lab[ ,3]^2 )                  #   Chroma
    H   = (atan2( Lab[ ,3], Lab[ ,2] ) * 180/pi) %% 360     #   Hue
    
    LCH = cbind(Lab[ ,1],C,H)
    rownames(LCH) = rownames(Lab)
    colnames(LCH) = c('L','Cab','Hab')
    
    return( LCH )
    }
    
    
LabfromLCHab    <- function( LCHab )
    {
    LCHab = prepareNxM(LCHab)
    if( is.null(LCHab) )  return(NULL)
    
    theta   = LCHab[ ,3] * pi/180
    
    Lab = cbind(LCHab[ ,1], LCHab[ ,2]*cos(theta), LCHab[ ,2]*sin(theta) )
    rownames(Lab) = rownames(LCHab)
    colnames(Lab) = c('L','a','b')
    
    return( Lab )
    }
    
    
################      LCHuv  <-->  Luv    #######################

LCHuvfromLuv    <- function( Luv )
    {
    Luv = prepareNxM(Luv)
    if( is.null(Luv) )  return(NULL)
    
    C   = sqrt( Luv[ ,2]^2 +  Luv[ ,3]^2 )                  #   Chroma
    H   = (atan2( Luv[ ,3], Luv[ ,2] ) * 180/pi) %% 360     #   Hue
    
    LCH = cbind(Luv[ ,1],C,H)
    rownames(LCH) = rownames(Luv)
    colnames(LCH) = c('L','Cuv','Huv')
    
    return( LCH )
    }
    
    
LuvfromLCHuv    <- function( LCHuv )
    {
    LCHuv = prepareNxM(LCHuv)
    if( is.null(LCHuv) )  return(NULL)
    
    theta   = LCHuv[ ,3] * pi/180
    
    Luv = cbind(LCHuv[ ,1], LCHuv[ ,2]*cos(theta), LCHuv[ ,2]*sin(theta) )
    rownames(Luv) = rownames(LCHuv)
    colnames(Luv) = c('L','a','b')
    
    return( Luv )
    }
        