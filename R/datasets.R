  

############    private data   ###########
#   p.dataIlluminants       standard illuminants, with 'standardized' XYZ
#   p.Ma                    3x3 adaptation matrices
    
#   an advantage of the private data in "sysdata.rda" is that these
#   do not have to be exposed, and therefore documented
savePrivateDatasets  <- function( .path="sysdata.rda" )
    {
    savevec = character(0)
        

    ##---------------       illuminants table    ------------------------##
    path    = "../inst/extdata/illuminants.txt"
    p.dataIlluminants = read.table( path, sep='\t', header=T, stringsAsFactors=F )
    
    p.dataIlluminants$XYZ = as.matrix( p.dataIlluminants[  ,c('X','Y','Z') ] )
    p.dataIlluminants[ c('X','Y','Z') ] = NULL
    n   = ncol( p.dataIlluminants )
    p.dataIlluminants   = p.dataIlluminants[ ,c(n,1:(n-1)) ]
    attr(p.dataIlluminants,"description") = readComments( path )
    savevec = c( savevec, "p.dataIlluminants" )
    
    #   list of adaptation matrices
    p.Ma    = list()
    p.Ma[[ "Bradford" ]]    = matrix( c(0.8951,0.2664,-0.1614,  -0.7502,1.7135,0.0367,  0.0389,-0.0685,1.0296), 3, 3, byrow=T )
    p.Ma[[ "VonKries" ]]    = matrix( c(0.40024,0.7076,-0.08081,  -0.2263,1.16532,0.0457,  0,0,0.91822), 3, 3, byrow=T )
    p.Ma[[ "MCAT02" ]]      = matrix( c( 0.7328, 0.4296, -0.1624,  -0.7036, 1.6975, 0.0061, 0.0030, 0.0136, 0.9834 ), 3, 3, byrow=T )
    p.Ma[[ "scaling" ]]     = diag(3)
    
    for( k in 1:length(p.Ma) )
        {
        rownames( p.Ma[[k]] )   = c('L','M','S')
        colnames( p.Ma[[k]] )   = c('X','Y','Z')
        }
        
    savevec = c( savevec, "p.Ma" )
    

    ##  finally ready to save it
    save( list=savevec, file=.path, compress='xz' )   #     'xz'  'gzip'  FALSE
    
    return( invisible(TRUE) )
    }    
    
readComments <- function( .path )
    {
    line    = readLines( .path, n=1024 )
    
    out = line[ grepl( "^[ \t]*#", line ) ]
    
    if( length(out) == 0 )  out = NULL
    
    return( out )
    }    