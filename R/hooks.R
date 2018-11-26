

#   private global variables.  The initial 'p.' means private

#   This group is read from sysdata.rda
#   p.dataCCT       the Robertson table for CCT.  It does not have to be unlocked

#   This group is assigned during .onAttach()
p.microbenchmark    = FALSE     # logical value, whether the package microbenchmark is loaded.  Once assigned it need never change.
p.uvfromMired       = NULL      # a pair of splines. Once created it need never change.
    
.onLoad <- function( libname, pkgname )
    {   
    #   at this point all globals seem to be unlocked
        
    #packageStartupMessage( ".onLoad() environment ", environmentIsLocked(asNamespace('spacesXYZ')), '\n' )   #  a paradox the environment is not locked !
    #packageStartupMessage( ".onLoad() p.uvfromMired ", bindingIsLocked( "p.uvfromMired", asNamespace('spacesXYZ') ), '\n'  )
    
    p.microbenchmark    <<- base::requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "p.microbenchmark=", p.microbenchmark, '\n' )
           
    p.uvfromMired       <<- list()
    p.uvfromMired[[1]]  <<- stats::splinefun( p.dataCCT$mired, p.dataCCT$u, method='fmm' )     # for u CIE 1960
    p.uvfromMired[[2]]  <<- stats::splinefun( p.dataCCT$mired, p.dataCCT$v, method='fmm' )     # for v CIE 1960
    }
    
    
.onAttach <- function( libname, pkgname )
    {
    #   at this point all globals seem to be locked
    
    #packageStartupMessage( ".onAttach() ", utils::str(p.uvfromMired), '\n' )     shows 2 functions OK
        
    #packageStartupMessage( ".onAttach() ", p.microbenchmark , '\n' )             shows TRUE OK
        
        
    # packageStartupMessage( ".onAttach() ", environmentIsLocked(asNamespace('spacesXYZ')), '\n', collapse=' ' )  a paradox the environment is locked !
    
    #unlockBinding( "p.microbenchmark", asNamespace('spacesXYZ') )     # asNamespace(pkgname) here generates a NOTE ! 
    #p.microbenchmark    <<- requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "p.microbenchmark=", p.microbenchmark, '\n' )
    
    #   make 2 splinefuns here, because it is unsafe to save them in sysdata.rda
    #unlockBinding( "p.uvfromMired", asNamespace('spacesXYZ') )            # asNamespace(pkgname) here generates a NOTE !     
    #p.uvfromMired       <<- list()
    #p.uvfromMired[[1]]  <<- splinefun( p.dataCCT$mired, p.dataCCT$u, method='fmm' )     # for u CIE 1960
    #p.uvfromMired[[2]]  <<- splinefun( p.dataCCT$mired, p.dataCCT$v, method='fmm' )     # for v CIE 1960
    }

    
    