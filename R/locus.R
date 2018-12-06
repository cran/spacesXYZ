
#   mired   desired vector of mired (aka mirek) values
#   xyz     path to xyz 1931 responsivities
#   c2      for Planck's formula  nm*K
#
#   returns data.frame with these columns
#       mired   aka beta
#       u
#       v
#       up      1st deriv, wrt beta
#       vp      1st deriv, wrt beta
#       upp     2nd deriv, wrt beta
#       vpp     2nd deriv, wrt beta
#       slope   for comparison with Robertson table

makePrecisionLocus <- function( mired=c( seq(0,390,by=10), seq(400,1000,by=25) ),  xyz="../inst/extdata/ciexyz31_1.csv",  c2=1.4388e7 )
    {
    xyz = read.table( xyz, header=T, sep=',' )
    
    #   split into wavelength and matrix
    lambda  = as.numeric( xyz$Wavelength )  # usually 360 to 830
    xyz     = as.matrix( xyz[ ,2:4] )       #; print( str(xyz) )
    
    out             = matrix( NA_real_, length(mired), 1 + 2*3 )
    colnames(out)   = c('mired','u','v','up','vp','upp','vpp')
    out[ ,1]        = mired
    
    u.num   =   xyz %*% c(4,0,0)
    v.num   =   xyz %*% c(0,6,0)
    denom   =   xyz %*% c(1,15,3)       # since denominator is the same for both u and v
    
    for( i in 1:length(mired) )
        {
        beta    = mired[i]
    
        spec012 = planckFormula( beta, lambda, c2 )     # radiation at every wavelength, plus 2 derivatives wrt beta
        
        # denom012 has 3 numbers, first the value and then 2 derivatives wrt beta
        denom012    = crossprod( spec012, denom )   # same as   t(spec012) %*% denom  
        
        for( cc in c('u','v') )
            {
            if( cc == 'u' )
                numer   = u.num
            else
                numer   = v.num
            
            #   id = either c('u','up','upp')   or   c('v','vp','vpp') 
            id      = paste( cc, c('','p','pp'), sep='' )
            
            numer012        = crossprod( spec012, numer )   # same as t(spec012) %*% numer
            out[i,id[1]]    = numer012[1] / denom012[1]
            out[i,id[2]]    = (numer012[2] - out[i,id[1]]*denom012[2]) / denom012[1]
            out[i,id[3]]    = (numer012[3] - 2*out[i,id[2]]*denom012[2] - out[i,id[1]]*denom012[3]) / denom012[1]
            }
        }
        
    out = as.data.frame(out)
    
    if( FALSE )
        #   append slope, to compare with Robertson table
        out = cbind( out, slope = -out$up/out$vp )
    
    return(out)
    }
    
    
#   beta    mired, 0 is OK
#   lambda  vector of wavelengths, in nm.  Length is 471
#   c2      in Planck's formula
#
#   returns a 3x471 matrix with columns
#       radiation at lambda
#       drad/dbeta      1st deriv, wrt beta
#       d2rad/d2beta    2nd deriv, wrt beta
    
planckFormula <- function( beta, lambda, c2 )
    {
    out = matrix( NA_real_, length(lambda), 3 )
    colnames(out)   = c( 'r', 'rp', 'rpp' )
    
    a       = 1.e-6*c2 / lambda  
    
    if( beta == 0 ) 
        {
        #   special case, the planck function is undefined when beta=0.
        #   But we get a nice analalytic function by using the function beta*planck instead.
        #   Since this amounts to scaling both numerator and denominator by the same amount, it cancels out.
        #   The derivatives at 0 are easily taken from the Taylor series computed in Maple.
        #   We could have used beta*planck for all beta, but the 2 derivatives are more complicated.
        out[ ,1]    = 1/a
        out[ ,2]    = -1/2
        out[ ,3]    = a/6
        }
    else
        {
        e           = exp( a * beta )
        out[ ,1]    =   (e - 1)^(-1)
        out[ ,2]    =   -a * e * (e - 1)^(-2)
        out[ ,3]    =   a^2 * e * (e + 1) * (e - 1)^(-3)
        }
        
        
    out = lambda^(-5) * out
    
    return( out )
    }
