

#   Lab0        Nx3 or 1x3 matrix,  
#   Lab1        Nx3 or 1x3 matrix,  
#
#   returns     numeric N-vector of pairwise differences

DeltaE  <-  function( Lab0, Lab1 )
    {
    Lab0 = prepareNxM( Lab0 )
    if( is.null(Lab0) )  return(NULL)
        
    Lab1 = prepareNxM( Lab1 )
    if( is.null(Lab1) )  return(NULL)
    
    if( nrow(Lab0)==1  &&  1<nrow(Lab1) )
        #   replicate Lab0
        Lab0    = matrix(Lab0,nrow(Lab1),3,byrow=TRUE)
        
    if( 1<nrow(Lab0)  &&  nrow(Lab1)==1 )
        #   replicate Lab1
        Lab1    = matrix(Lab1,nrow(Lab0),3,byrow=TRUE)
        
        
    if( nrow(Lab0) != nrow(Lab1) )
        {
        log.string( ERROR, "nrow(Lab0) = %d != %d = nrow(Lab1).",
                            nrow(Lab0), nrow(Lab1) )
        return(NULL)
        }
        
    #   Lab0 and Lab1 passed all checks
    out = DeltaE.1976( Lab0, Lab1 )
    
    return(out)
    }
        
        
DeltaE.1976  <-  function( Lab0, Lab1 )
    {        
    delta   = Lab0 - Lab1
    
    sqrt( rowSums(delta*delta) )
    }
    
    