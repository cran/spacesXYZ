

#   returns time in seconds, from an arbitrary origin
gettime <- function()
    {
    if( p.microbenchmark )
        return( microbenchmark::get_nanotime() * 1.e-9 )
    else
        return( as.double( base::Sys.time() ) )
    }
    


    
    
###########     argument processing     ##############
#
#   A   a non-empty numeric NxM matrix, or something that can be converted to be one
#
#   returns such a matrix, or NULL in case of error
#
prepareNxM  <-  function( A, M=3 )
    {    
    ok  = is.numeric(A) &&  0<length(A)  &&  (length(dim(A))<=2)  # &&  (0<M) 
    
    ok  = ok  &&  ifelse( is.matrix(A), ncol(A)==M, ((length(A) %% M)==0)  )
    
    if( ! ok )
        {
        #print( "prepareNx3" )
        #print( sys.frames() )
        mess    = substr( as.character(A)[1], 1, 10 )
        #arglist = list( ERROR, "A must be a non-empty numeric Nx3 matrix (with N>0). A='%s...'", mess )
        #do.call( log_level, arglist, envir=parent.frame(n=3) )
        #myfun   = log_level
        #environment(myfun) = parent.frame(3)
        
        Aname = deparse(substitute(A))        
        
        #   notice hack to make log_level() print name of parent function, and *NOT* prepareNxM()        
        log_level( ERROR, "Argument '%s' must be a non-empty numeric Nx%d matrix (with N>0). %s='%s...'", 
                                    Aname, M, Aname, mess, .topcall=sys.call(-1L) )
        return(NULL)
        }
    
    if( ! is.matrix(A) )
        A = matrix( A, ncol=M, byrow=TRUE )
        
    return( A )
    }
        