
#   x       in strictly increasing order
#   y       y(x)
#   yp      y'(x)
#   ypp     y''(x)
#   All of these must be the same length

#   returns a function fun( x, deriv=0 ), where deriv can be 0,1,2.
#   The function is class C^2.
#   This function is designed to have the same argument interface as stats::splinefun().
#   For x outside the range, returns NA

quinticfun <- function( x, y, yp, ypp )
    {
    n   = length(x)
    
    ok  = 2<=n  &&  length(y)==n  &&  length(yp)==n  &&  length(ypp)==n
    if( ! ok )
        {
        log.string( ERROR, "Invalid lengths." )
        return(NULL)
        }
        
    ok  = all( 0 < diff(x) )
    if( ! ok )
        {
        log.string( ERROR, "x is not strictly increasing." )
        return(NULL)
        }
        
    #   all quintics have no contant term, only 5 coefficients
    a   = numeric( 5 )
    coeffmat = matrix( NA_real_, n-1, length(a) )

    for( i in 1:(n-1) )
        {
        a[1]    = yp[i] 
        a[2]    = 0.5 * ypp[i]
        
        b       = x[i+1] - x[i]
        bvec    = b ^ (1:5)
        
        mat = matrix( c( bvec[3:5], c(3,4,5)*bvec[2:4], c(6,12,20)*bvec[1:3] ), 3, 3, byrow=TRUE )
        rhs = c( y[i+1]-y[i] - bvec[2]*a[2] - b*a[1], yp[i+1] - b*2*a[2] - a[1], ypp[i+1] - 2*a[2] )
        sol = base::solve( mat, rhs )
        
        if( inherits(sol,"try-error" ) )    #class(sol) == "class-error")
            {
            log.string( ERROR, "Failure of 3x3 solve()." )
            return(NULL)
            }
        
        a[3:5]  = sol
        coeffmat[i, ]   = a
        }
    
    outfun  <- function( xin, deriv=0 )
        {
        ok  = is.numeric(xin)  &&  0<length(xin)
        if( ! ok )
            {
            log.string( ERROR, "xin is invalid." )
            return(NULL)
            }      

        ok  = deriv %in% c(0,1,2)
        if( ! ok )
            {
            log.string( ERROR, "deriv is invalid." )
            return(NULL)
            }      
        
        ivec    = base::findInterval( xin, x )
        
        m   = length(xin)
        out = rep( NA_real_, m )
        
        for( k in 1:m )
            {
            i   = ivec[k]
            
            if( i == 0 )    next    # reject silently, out of range

            if( i == n )
                {
                #   examine more closely
                epsilon = 1.e-12
                if( x[n] + epsilon < xin[k] )   next     # reject silently, out of range

                xin[k] = x[n]   # override value, and next if block will take care of it
                }
                
            if( xin[k] == x[i] )
                {
                #   exactly on a knot !
                out[k] = c( y[i], yp[i], ypp[i] )[ deriv+1 ]                
                next
                }
            
            #   xin[k] is in interval i
            
            u   = xin[k] - x[i] # normalize to interval x[i] to x[i+1]
            a   = coeffmat[i, ]
            
            if( deriv == 0 )
                out[k]    = ( ( ( (a[5]*u + a[4])*u + a[3])*u + a[2])*u + a[1])*u  +  y[i]
            else if( deriv == 1 )
                out[k]    = ( ( (5*a[5]*u + 4*a[4])*u + 3*a[3])*u + 2*a[2])*u + a[1]
            else
                out[k]    = ( (20*a[5]*u + 12*a[4])*u + 6*a[3])*u + 2*a[2]
            }
            
        return( out )
        }
        
    return( outfun )
    }
    