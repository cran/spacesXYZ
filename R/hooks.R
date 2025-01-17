

#   private global variables.  The initial 'p.' means private


#   This group is assigned during .onLoad()
p.microbenchmark        = FALSE     # logical value, whether the package microbenchmark is loaded.  Once assigned it need never change.
p.uvCubicsfromMired     = NULL      # a pair of cubic splines. Once created it need never change.
p.uvQuinticsfromMired   = NULL      # a pair of quintic splines. Once created it need never change.


.onLoad <- function( libname, pkgname )
    {
    #   at this point all globals seem to be unlocked

    #packageStartupMessage( ".onLoad() environment ", environmentIsLocked(asNamespace('spacesXYZ')), '\n' )   #  a paradox the environment is not locked !
    #packageStartupMessage( ".onLoad() p.uvCubicsfromMired ", bindingIsLocked( "p.uvCubicsfromMired", asNamespace('spacesXYZ') ), '\n'  )

    p.microbenchmark    <<- base::requireNamespace( 'microbenchmark', quietly=TRUE )

    if( requireNamespace( 'logger', quietly=FALSE ) )
        {
        #   log_formatter( formatter_mine )
        #   layout_mine and appender_mine are defined in logger.R
        log_formatter( logger::formatter_sprintf, namespace=pkgname )   # force sprintf(), even if glue is installed
        log_layout( layout_mine, namespace=pkgname )                    # put fn() between timestamp and the msg
        log_appender( appender_mine, namespace=pkgname )                # maybe stop on ERROR or FATAL
        log_threshold( WARN, namespace=pkgname )                        # default is INFO
        }



    #   make functions ufun and vfun.  These are all class C^2, but cubics and quintics.
    p.uvCubicsfromMired         <<- list()
    p.uvCubicsfromMired$ufun    <<- stats::splinefun( spacesXYZ::RobertsonLocus$mired, spacesXYZ::RobertsonLocus$u, method='fmm' )     # for u CIE 1960
    p.uvCubicsfromMired$vfun    <<- stats::splinefun( spacesXYZ::RobertsonLocus$mired, spacesXYZ::RobertsonLocus$v, method='fmm' )     # for v CIE 1960
    p.uvCubicsfromMired$miredInterval <<- range( spacesXYZ::RobertsonLocus$mired )

    p.uvQuinticsfromMired       <<- list()
    p.uvQuinticsfromMired$ufun  <<- quinticfun( spacesXYZ::PrecisionLocus$mired, spacesXYZ::PrecisionLocus$u, spacesXYZ::PrecisionLocus$up, spacesXYZ::PrecisionLocus$upp )
    p.uvQuinticsfromMired$vfun  <<- quinticfun( spacesXYZ::PrecisionLocus$mired, spacesXYZ::PrecisionLocus$v, spacesXYZ::PrecisionLocus$vp, spacesXYZ::PrecisionLocus$vpp )
    p.uvQuinticsfromMired$miredInterval <<- range( spacesXYZ::PrecisionLocus$mired )

    #   packageStartupMessage( ".onLoad().  Made quintic pair." )
    }


.onAttach <- function( libname, pkgname )
    {
    #   at this point all globals seem to be locked

    #packageStartupMessage( ".onAttach() ", utils::str(p.uvCubicsfromMired), '\n' )     shows 2 functions OK

    #packageStartupMessage( ".onAttach() ", p.microbenchmark , '\n' )             shows TRUE OK


    # packageStartupMessage( ".onAttach() ", environmentIsLocked(asNamespace('spacesXYZ')), '\n', collapse=' ' )  a paradox the environment is locked !

    #unlockBinding( "p.microbenchmark", asNamespace('spacesXYZ') )     # asNamespace(pkgname) here generates a NOTE !
    #p.microbenchmark    <<- requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "p.microbenchmark=", p.microbenchmark, '\n' )

    #   make 2 splinefuns here, because it is unsafe to save them in sysdata.rda
    #unlockBinding( "p.uvCubicsfromMired", asNamespace('spacesXYZ') )            # asNamespace(pkgname) here generates a NOTE !
    #p.uvCubicsfromMired       <<- list()
    #p.uvCubicsfromMired[[1]]  <<- splinefun( p.dataCCT$mired, p.dataCCT$u, method='fmm' )     # for u CIE 1960
    #p.uvCubicsfromMired[[2]]  <<- splinefun( p.dataCCT$mired, p.dataCCT$v, method='fmm' )     # for v CIE 1960
    }


