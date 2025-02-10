

p.MCAT02        = matrix( c(0.7328, 0.4296, -0.1624,  -0.7036, 1.6975, 0.0061,  0.0030, 0.0136, 0.9834), 3, 3, byrow=TRUE )
p.MCAT02_inv    = base::solve( p.MCAT02 )

p.MHPE          = matrix( c(0.38971, 0.68898, -0.07868,  -0.22981, 1.18340, 0.04641, 0, 0, 1), 3, 3, byrow=TRUE )


#   surround viewing conditions
p.surroundconditions = rbind(
    'Average'   = c( 0.69, 1, 1 ),
    'Dim'       = c( 0.59, 0.9, 0.9 ),
    'Dark'      = c( 0.525,  0.8, 0.8 )
    )
colnames( p.surroundconditions ) = c( 'c', 'N_c', 'F' )


#   Table 16.2 Data for conversion from hue angle to hue quadrature
p.HUE_DATA    = rbind(
    'h_i'   = c( 20.14, 90.00, 164.25, 237.53, 380.14 ),
    'e_i'   = c( 0.8, 0.7, 1.0, 1.2, 0.8 ),
    'H_i'   = c( 0.0, 100.0, 200.0, 300.0, 400.0 )
    )
colnames( p.HUE_DATA )    = c( 'Red', 'Yellow', 'Green', 'Blue', 'Red' )


CIECAM02fromXYZ   <- function( XYZ, XYZ_w, L_A=100, Y_b=20, surround='Average', discount=TRUE )
    {
    XYZ = prepareNxM( XYZ, M=3 )
    if( is.null(XYZ) )  return(NULL)
    
    if( is.character(XYZ_w) )
        {
        white   = XYZ_w[1]
        
        XYZ_w   = 100 * standardXYZ( white )
        }
    else
        {
        white   = NA_character_
        
        if( length(XYZ_w) == 3 )
            {
            mat = 100 * standardXYZ( '*' )
            
            mat_w   = matrix( XYZ_w, nrow=nrow(mat), ncol=3, byrow=TRUE )
            
            delta   = abs(mat - mat_w)
            
            delta   = rowSums(delta)
            
            idx     = which( delta < 5.e-12 )
            
            if( 0 < length(idx) )
                white = rownames(mat)[ idx[1] ]
            }
        }
        
    ok = is.numeric(XYZ_w)  &&  length(XYZ_w)==3  && all(is.finite(XYZ_w))  &&   all(0 < XYZ_w)        
    if( ! ok )
        {
        log_level( ERROR, "white='%s' is invalid.", as.character(XYZ_w) )
        return(NULL)
        }                

    #   do not want a 1x3 matrix here
    dim(XYZ_w)      = NULL
    names(XYZ_w)    = c('X','Y','Z')
    

    #   get surround conditions parameters
    i   = pmatch( tolower(surround), tolower(rownames(p.surroundconditions)) )
    if( is.na(i) )
        {
        log_level( ERROR, "surround='%s' is invalid.", surround )
        return(NULL)
        }

    #   load surr with surr$c,  surr$N_c,  surr$F
    surr    = as.list( p.surroundconditions[i, ] )  #; print(surr)

    #   put viewing conditions in viewing
    viewing     = viewing_condition_dependent_parameters( XYZ_w[2], L_A, Y_b )  #; print(viewing)

    D   = ifelse( ! discount, surr$F * (1 - (1 / 3.6) * exp((-L_A - 42) / 92)), 1 ) #; cat( "D=", D, '\n' )

    #   compute RGB for white, post_adaptation
    RGB_w   = p.MCAT02  %*%  XYZ_w

    RGB_wc  = full_chromatic_adaptation_forward( RGB_w, RGB_w, XYZ_w[2], D )

    RGB_wp  = HPE_fundamentals( RGB_wc )

    RGB_wa  = post_adaptation_non_linear_response_compression_forward( RGB_wp, viewing$F_L )
    
    #cat( "RGB_w =", RGB_w, '\n' )
    #cat( "RGB_wc =", RGB_wc, '\n' )
    #cat( "RGB_wp =", RGB_wp, '\n' )
    #cat( "RGB_wa =", RGB_wa, '\n' )


    m   = nrow(XYZ)

    h   = rep(NA_real_,m)
    J   = rep(NA_real_,m)
    Q   = rep(NA_real_,m)    
    C   = rep(NA_real_,m)
    M   = rep(NA_real_,m)
    
    for( i in 1:m )
        {
        RGB = p.MCAT02  %*%  XYZ[i, ]

        #   compute RGB for white, post_adaptation

        RGB_c   = full_chromatic_adaptation_forward( RGB, RGB_w, XYZ_w[2], D )

        RGB_p   = HPE_fundamentals( RGB_c )

        RGB_a   = post_adaptation_non_linear_response_compression_forward( RGB_p, viewing$F_L )
        
        #cat( "RGB =", RGB, '\n' )
        #cat( "RGB_c =", RGB_c, '\n' )
        #cat( "RGB_p =", RGB_p, '\n' )
        #cat( "RGB_a =", RGB_a, '\n' )

        res = CIECAM02_from_RGBa( RGB_a, RGB_wa, L_A, Y_b, viewing, surr )

        if( is.null(res) )  next

        h[i]    = res$h
        J[i]    = res$J
        Q[i]    = res$Q                
        C[i]    = res$C
        M[i]    = res$M
        }

    rnames  = row.names(XYZ)
    if( is.null(rnames) )   rnames = 1:m
    
    out = data.frame( row.names=rnames )
    
    colnames(XYZ)   = c('X','Y','Z')
    
    out$XYZ = XYZ
    out$h   = h    
    #   out$ab  = cbind( cos( h * pi/180 ), sin( h * pi/180 ) )
    out$J   = J
    out$Q   = Q    
    out$C   = C
    out$M   = M
    
    out$Jp  = 1700 * J / (1000 + 7*J)
    
    Mp      = log( 1 + 0.0228*M ) / 0.0228      # not log10()
    out$abp = Mp * cbind( cos( h * pi/180 ), sin( h * pi/180 ) )
    colnames(out$abp)   = c('a','b')

    attr( out, 'white' )    = white
    attr( out, 'XYZ_w' )    = XYZ_w
    attr( out, 'L_A' )      = L_A
    attr( out, 'Y_b' )      = Y_b
    attr( out, 'surround' ) = surround
    attr( out, 'discount' ) = discount
    attr( out, 'D' )        = D
    
    return( out )
    }





CIECAM02_from_RGBa    <- function( RGB_a, RGB_wa, L_A, Y_b, viewing, surround )
    {
    ab  = opponent_colour_dimensions_forward( RGB_a )     #; cat( "ab=", ab, '\n' )

    h   = hue_angle( ab )  #; cat( "h=", h, '\n' )

    e_t = eccentricity_factor( h )    #;     cat( "e_t = ", e_t, '\n' )

    A   = achromatic_response_forward( RGB_a, viewing$N_bb )    #;     cat( "A = ", A, '\n' )

    A_w = achromatic_response_forward( RGB_wa, viewing$N_bb )   #;     cat( "A_w = ", A_w, '\n' )
    


    J   = lightness_correlate( A, A_w, surround$c, viewing$z )

    Q   = brightness_correlate( surround$c, J, A_w, viewing$F_L )

    C   = chroma_correlate( J, viewing$n, surround$N_c, viewing$N_cb, e_t, ab[1], ab[2], RGB_a )

    M   = colourfulness_correlate( C, viewing$F_L )

    out = list( h=h, J=J, Q=Q, C=C, M=M )

    return( out )
    }


#   CAM_1, CAM_2    data frames returned from CIECAM02fromXYZ()

DeltaE_CAM02    <- function( CAM_1, CAM_2 )
    {
    n1  = nrow(CAM_1)
    n2  = nrow(CAM_2)
    
    UCS_1   = cbind( CAM_1$Jp, CAM_1$abp )      # UCS_CAM02( CAM_1 )
    UCS_2   = cbind( CAM_2$Jp, CAM_2$abp )      # UCS_CAM02( CAM_2 )
    
    # replicate matrix rows, if necessary
    if( n1 == 1  &&  1 < n2 )
        UCS_1   = matrix( UCS_1, nrow=n2, ncol=3, byrow=TRUE )
    else if( 1 < n1  &&  n2 == 1 )
        UCS_2   = matrix( UCS_2, nrow=n1, ncol=3, byrow=TRUE )
    
    if( nrow(UCS_1) != nrow(UCS_2) )
        {
        log_level( ERROR, "nrow(CAM_1) = %d   and   nrow(CAM_2) = %d, and this is invalid.", n1, n2 )
        return(NULL)
        }
    
    return( sqrt( rowSums( (UCS_1 - UCS_2)^2 ) ) )
    }
    
    



viewing_condition_dependent_parameters <- function( Y_w, L_A, Y_b )
    {
    k   = 1 / (5*L_A + 1)

    F_L = 0.2*k^4*(5*L_A)  +  0.1*(1-k^4)^2*((5*L_A)^(1/3))

    n   = Y_b / Y_w

    N_bb = N_cb = 0.725 * (1 / n) ^ 0.2

    z   = 1.48 + sqrt(n)

    out = list( n=n, F_L=F_L, N_bb=N_bb, N_cb=N_cb, z=z )

    return( out )
    }

full_chromatic_adaptation_forward <- function( RGB, RGB_w, Y_w, D )
    {
    myfun <-  function( x_xw )  {
        x   = x_xw[1]
        xw  = x_xw[2]
        return( ((Y_w * D / xw) + 1 - D) * x )
        }

    apply( cbind(RGB,RGB_w), 1, myfun )
    }

HPE_fundamentals    <- function( RGB_c )
    {
    p.MHPE %*% (p.MCAT02_inv %*% RGB_c)
    }

post_adaptation_non_linear_response_compression_forward <- function( RGB_p, F_L )
    {
    ((((400 * (F_L * RGB_p / 100) ** 0.42) /  (27.13 + (F_L * RGB_p / 100) ** 0.42))) + 0.1)       #  (16.15)
    }

opponent_colour_dimensions_forward  <- function( RGB_a )
    {
    R   = RGB_a[1]
    G   = RGB_a[2]
    B   = RGB_a[3]

    ab = c( R - 12 * G / 11 + B / 11, (R + G - 2 * B) / 9 )

    return( ab )
    }


#   rectangular to degrees
hue_angle   <- function( ab )
    {
    h   = (atan2(ab[2], ab[1]) * 180/pi ) %% 360

    return( h )
    }

eccentricity_factor <- function( h )
    {
    e_t = (cos(2 + h * pi/180) + 3.8) / 4

    return( e_t )
    }



achromatic_response_forward <- function( RGB_a, N_bb )
    {
    R   = RGB_a[1]
    G   = RGB_a[2]
    B   = RGB_a[3]

    A = (2 * R + G + (1 / 20) * B - 0.305) * N_bb   # (16.23)

    return( A )
    }

lightness_correlate <- function( A, A_w, c, z )
    {
    J = 100 * (A / A_w) ^ (c * z)     # (16.24)

    return( J )
    }

brightness_correlate    <- function( c, J, A_w, F_L )
    {
    Q = (4 / c) * sqrt(J / 100) * (A_w + 4) * F_L^0.25      # (16.25)

    return( Q )
    }


chroma_correlate    <- function( J, n, N_c, N_cb, e_t, a, b, RGB_a )
    {
    Ra  = RGB_a[1]
    Ga  = RGB_a[2]
    Ba  = RGB_a[3]
    
    #   temporary_magnitude_quantity_forward(N_c, N_cb, e_t, a, b, RGB_a)
    t   = ( (50000/13) * N_c * N_cb) * (e_t * sqrt(a^2 + b^2)) / (Ra + Ga + (21/20)*Ba)    # (16.26)

    C   = t^0.9 * (J / 100)^0.5 * (1.64 - 0.29^n) ^ 0.73        # (16.27)

    return( C )
    }


colourfulness_correlate <- function( C, F_L )
    {
    M   = C * F_L^0.25

    return( M )
    }



##################      deadwood below  #######################
    
#   CAM     data frame returned from CIECAM02fromXYZ()
#           only 3 columns are used:  J, M, h
#
#   returns matrix with 3 columns:  Jp, ap, bp

UCS_CAM02   <- function( CAM )
    {
    Jp  = 1700 * CAM$J / (1000 + 7*CAM$J)
    
    Mp  = log10( 1 + 0.0228*CAM$M ) / 0.0228
    
    ap  = Mp*cos( CAM$h * pi/180 )
    bp  = Mp*sin( CAM$h * pi/180 )
   
    return( cbind(Jp,ap,bp) )
    }