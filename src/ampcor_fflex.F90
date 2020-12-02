 
SUBROUTINE ampcor_fflex(c_refimg_map,c_srchimg_map, &
        nx_im1, ny_im1, nx_im2, ny_im2, & !!!size original image 
        x0_shift, y0_shift, nx_win, ny_win, nx_search, ny_search, & !!! init shift, subwindow size, search distance
        nx_step, ny_step, x_start, y_start, & !!! steps, position in slc where to start for processing
        jstart_master, jstart_slave, & !!! position where portion of slc was chunked (c_refimg_map,c_srchimg_map)
        searchmap_nx, searchmap_ny, x_st_searchmap, y_st_searchmap, rgsp_searchmap, azsp_searchmap, USE_SEARCHMAP_FLAG, &
        msk, x_st_msk, y_st_msk, rgsp_msk, azsp_msk, USE_MASK_FLAG, &
        shiftm_x, shiftm_y, x_st_shiftm, y_st_shiftm, rgsp_shiftm, azsp_shiftm, USE_SHIFTMAP_FLAG, &
        x0_polyshift, y0_polyshift, USE_POLYSHIFT_FLAG, &
        nx_off, ny_off, &
        xoff_out_map,yoff_out_map,r_snr_map, r_cov_map1, r_cov_map2, r_cov_map3, &
        nx_ref,ny_ref, nx_srch,ny_srch, & ! hidden size of master and slave slc chunks
        nx_ma, ny_ma, nx_se, ny_se, nx_sh, ny_sh) ! hidden

    !f2py complex(4) intent(in) :: c_refimg_map
    !f2py integer(kind=4) intent(hide),depend(c_refimg_map) :: nx_ref=shape(c_refimg_map,0), ny_ref=shape(c_refimg_map,1)

    !f2py complex(4) intent(in) :: c_srchimg_map
    !f2py integer(kind=4) intent(hide),depend(c_srchimg_map) :: nx_srch=shape(c_srchimg_map,0), ny_srch=shape(c_srchimg_map,1)

    !f2py integer(kind=4) intent(in) :: nx_im1, ny_im1, nx_im2, ny_im2
    !f2py integer(kind=4) intent(in) :: x0_shift, y0_shift, nx_win, ny_win, nx_search, ny_search
    !f2py integer(kind=4) intent(in) :: nx_step, ny_step, x_start, y_start
    !f2py integer(kind=4) intent(in) :: jstart_master, jstart_slave

    !f2py intent(in) :: searchmap_nx, searchmap_ny
    !f2py integer(kind=4) intent(hide),depend(searchmap_nx, searchmap_ny) :: nx_se=shape(searchmap_nx,0), ny_se=shape(searchmap_nx,1) 
    !f2py integer intent(in) :: x_st_searchmap, y_st_searchmap, rgsp_searchmap, azsp_searchmap
    !f2py logical intent(in) :: USE_SEARCHMAP_FLAG

    !f2py integer intent(in) :: msk
    !f2py integer(kind=4) intent(hide),depend(msk) :: nx_ma=shape(msk,0), ny_ma=shape(msk,1)
    !f2py integer intent(in) :: x_st_msk, y_st_msk, rgsp_msk, azsp_msk
    !f2py logical intent(in) ::USE_MASK_FLAG

    !f2py real*8 intent(in) :: shiftm_y, shiftm_x
    !f2py integer(kind=4) intent(hide),depend(shiftm_y, shiftm_x) :: nx_sh=shape(shiftm_x,0), ny_sh=shape(shiftm_x,1)
    !f2py integer intent(in) :: x_st_shiftm, y_st_shiftm, rgsp_shiftm, azsp_shiftm
    !f2py logical intent(in) :: USE_SHIFTMAP_FLAG

    !f2py real*8 dimension(3), intent(in) :: x0_polyshift, y0_polyshift
    !f2py logical intent(in) :: USE_POLYSHIFT_FLAG
    !f2py integer(kind=4) intent(in) :: nx_off, ny_off

    !f2py real*8 intent(out) :: xoff_out_map,yoff_out_map,r_snr_map,r_cov_map1,r_cov_map2, r_cov_map3

    implicit none

    !
    integer i_ovs, i_idx, i_idy
    parameter(i_ovs=2)
    parameter(i_idx=512) ! was 2048
    parameter(i_idy=512) ! was 2018

    ! INPUT :
    integer(kind=4) ny_ma,nx_ma,nx_se, ny_se, ny_sh,nx_sh,nx_off,ny_off ! size for mask, searchmap, shiftmap, offsetmap
    integer(kind=4) x0_shift, y0_shift, nx_im1, ny_im1, nx_im2, ny_im2
    integer nx_search, ny_search
    integer nx_win, ny_win
    integer nx_step, ny_step
    integer(kind=4) nx_ref,ny_ref, nx_srch,ny_srch 
    integer(kind=4) x_start, y_start
    integer(kind=4) jstart_master, jstart_slave
    complex(4) c_refimg_map
    dimension c_refimg_map(nx_ref,ny_ref)
    complex(4) c_srchimg_map
    dimension c_srchimg_map(nx_srch,ny_srch)

    integer msk
    dimension msk(nx_ma, ny_ma)
    integer y_st_msk, x_st_msk, azsp_msk, rgsp_msk
    logical USE_MASK_FLAG

    real*8 searchmap_nx, searchmap_ny
    dimension searchmap_nx(nx_se, ny_se)
    dimension searchmap_ny(nx_se, ny_se)
    integer y_st_searchmap, x_st_searchmap, azsp_searchmap, rgsp_searchmap
    logical USE_SHIFTMAP_FLAG

    real*8 shiftm_y, shiftm_x
    dimension shiftm_y(nx_sh, ny_sh)
    dimension shiftm_x(nx_sh, ny_sh)
    integer y_st_shiftm, x_st_shiftm, azsp_shiftm, rgsp_shiftm
    logical USE_SEARCHMAP_FLAG

    real*8 y0_polyshift(3), x0_polyshift(3)
    logical USE_POLYSHIFT_FLAG

    ! OUTPUT:
    real*8 xoff_out_map,yoff_out_map
    dimension xoff_out_map(nx_off,ny_off)
    dimension yoff_out_map(nx_off,ny_off)
    real*8 r_snr_map
    dimension r_snr_map(nx_off,ny_off)
    real*8 r_cov_map1, r_cov_map2, r_cov_map3
    dimension r_cov_map1(nx_off,ny_off)
    dimension r_cov_map2(nx_off,ny_off)
    dimension r_cov_map3(nx_off,ny_off)

    !LOCAL
    real*8 xoff_out, yoff_out, r_snr, r_cov(3)
    complex(4) c_refimg
    dimension c_refimg(nx_win,ny_win)

    complex(4), dimension (:,:), allocatable :: c_srchimg
    
    integer(kind=4) j_off, j_im, i_off, i_im, j_msk, i_msk, i,k
    integer(kind=4) min_y0_shift, max_y0_shift, jmin_im1, jmax_im1, jmin_im2, jmax_im2
    logical process_this_line, process_this_point
    integer(kind=4) i_shift, j_shift, i_search, j_search, max_ny_search
    integer nx_m, ny_m, nx_s, ny_s
    integer(kind=4) i0m, j0m, i0s, j0s
    integer(kind=4) x_shift, y_shift

    complex c_osref(i_ovs*i_idx*i_ovs*i_idy)

    xoff_out_map(:,:)=0.0
    yoff_out_map(:,:)=0.0
    r_snr_map(:,:)=0.0
    r_cov_map1(:,:)=0.0
    r_cov_map2(:,:)=0.0
    r_cov_map3(:,:)=0.0

    !print *,nx_ref,ny_ref, nx_srch,ny_srch,x_start, y_start,jstart_master, jstart_slave

    !print *,'USE_POLYSHIFT_FLAG',USE_POLYSHIFT_FLAG

    !print *,'USE_SEARCHMAP_FLAG',USE_SEARCHMAP_FLAG
    !print *,'USE_MASK_FLAG',USE_MASK_FLAG
    !print *,'USE_SHIFTMAP_FLAG',USE_SHIFTMAP_FLAG

    ! initialize plan for fftw
    do i=3,14
        k=2**i
        call cfft1d_jpl(k,c_osref,0)
    end do

    DO j_off=1,ny_off

    j_im = (j_off-1)*ny_step + y_start
    !print *,j_off,j_im

    IF ( (j_im+y0_shift .GT. 2*ny_search+ny_win/2) .AND. &
        (j_im+y0_shift .LT. ny_im2-ny_win/2-2*ny_search) .AND. &
        (j_im .GT. ny_win/2) .AND. &
        (j_im .LT. ny_im1-ny_win/2) ) THEN

        IF (USE_MASK_FLAG) THEN

            j_msk = INT( (j_im - y_st_msk) / azsp_msk + 1) ! conversion from im y-coord to msk y-coord

            IF (SUM(msk(:,j_msk)) == 0) THEN
                process_this_line = .FALSE.
            ELSE
                process_this_line = .TRUE.
            ENDIF
        ELSE
            process_this_line = .TRUE.
        ENDIF
        ! process this line {{{
        IF (process_this_line) THEN
            jmin_im1 = INT(j_im-ny_win/2)
            jmax_im1 = INT(j_im+ny_win/2)

            min_y0_shift = y0_shift
            max_y0_shift = y0_shift

            IF (USE_POLYSHIFT_FLAG) THEN
                IF (y0_polyshift(1) .GT. 0) THEN
                    min_y0_shift = INT( y0_polyshift(1) + j_im*y0_polyshift(3))
                    max_y0_shift = INT( y0_polyshift(1) + nx_im1*y0_polyshift(2) + j_im*y0_polyshift(3))
                ELSE
                    min_y0_shift = INT( y0_polyshift(1) + nx_im1*y0_polyshift(2) + j_im*y0_polyshift(3) )
                    max_y0_shift = INT( y0_polyshift(1) + j_im*y0_polyshift(3))
                ENDIF
            ENDIF

            IF (USE_SHIFTMAP_FLAG) THEN
                j_shift = INT( (j_im - y_st_shiftm) / azsp_shiftm + 1 ) ! conversion from image y-coord to offsetmap y-corrd
                min_y0_shift = INT( MINVAL(shiftm_y(:,j_shift)) + min_y0_shift)
                max_y0_shift = INT( MAXVAL(shiftm_y(:,j_shift)) + min_y0_shift)
            ENDIF

            IF (USE_SEARCHMAP_FLAG) THEN
                j_search = INT( (j_im - y_st_searchmap) / azsp_searchmap + 1)
                max_ny_search = INT( MAXVAL( searchmap_ny(:,j_search)))
            ELSE
                max_ny_search = ny_search
            ENDIF

            jmin_im2 = INT( j_im + min_y0_shift - ny_win/2 - 2*max_ny_search )
            jmax_im2 = INT( j_im + max_y0_shift + ny_win/2 + 2*max_ny_search )

            IF ((jmin_im2 .LT. 0) .OR. (jmax_im2 .GT. ny_im2)) THEN
                process_this_line = .False.
            ENDIF

            ! loop over x (range) {{{
            DO i_off=1,nx_off ! loop over x-direction (range for SAR)
            i_im = INT( (i_off-1)*nx_step + x_start) ! from output offset x-coord to im x-coord

            IF (USE_MASK_FLAG) THEN
                j_msk = INT( (j_im - y_st_msk) / azsp_msk) ! conversion from im y-coord to msk y-coord
                i_msk = INT( (i_im - x_st_msk) / rgsp_msk) 

                IF (msk(i_msk,j_msk) .EQ. 1) THEN
                    process_this_point = .TRUE.
                ELSE
                    process_this_point = .FALSE.
                ENDIF
            ELSE
                process_this_point = .TRUE.

            ENDIF
            ! process this point {{{
            IF (process_this_point) THEN

                x_shift = x0_shift
                y_shift = y0_shift

                IF (USE_POLYSHIFT_FLAG) THEN
                    x_shift = INT(x0_polyshift(1) + i_im*x0_polyshift(2) + j_im*x0_polyshift(3))
                    y_shift = INT(y0_polyshift(1) + i_im*y0_polyshift(2) + j_im*y0_polyshift(3))
                ENDIF

                IF (USE_SHIFTMAP_FLAG) THEN
                    j_shift = INT( (j_im - y_st_shiftm) / azsp_shiftm + 1)
                    i_shift = INT( (i_im - x_st_shiftm) / rgsp_shiftm + 1)
                    x_shift = INT( shiftm_x(i_shift,j_shift) + x_shift)
                    y_shift = INT( shiftm_y(i_shift,j_shift) + y_shift)
                ENDIF

                IF (USE_SEARCHMAP_FLAG) THEN
                    j_search = INT( (j_im - y_st_searchmap) / azsp_searchmap + 1)
                    i_search = INT( (i_im - x_st_searchmap) / rgsp_searchmap + 1)
                    nx_search = INT( searchmap_nx(i_search,j_search) )
                    ny_search = INT( searchmap_ny(i_search,j_search) )
                ENDIF

                ! call ampcor_flex {{{
                IF ( (i_im+x_shift .GT. 2*nx_search + nx_win/2) .AND. &
                    ( i_im + x_shift .LT. nx_im2 - nx_win/2 - 2*nx_search) .AND. &
                    ( i_im .GT. nx_win/2) .AND. &
                    ( i_im .LT. nx_im1 - nx_win/2) ) THEN

                    ! # extract ref and srch sub-images from lines previously read
                    nx_m = nx_win
                    ny_m = ny_win
                    nx_s = nx_win + 4*nx_search
                    ny_s = ny_win + 4*ny_search

                    i0m = INT(i_im - nx_win/2 + 1)
                    j0m = INT(j_im - ny_win/2 - jstart_master + 1)

                    i0s = INT(i_im + x_shift - nx_win/2 - 2*nx_search +1 )
                    j0s = INT(j_im + y_shift - ny_win/2 - 2*ny_search - jstart_slave + 1)

                    IF ( (i0m .GT. 0) .AND. (j0m .GT. 0) .AND. (i0s .GT. 0) .AND. &
                        (j0s .GT. 0) .AND. ((i0m+nx_m) .LT. nx_ref) .AND. &
                        ((j0m+ny_m) .LT. ny_ref) .AND. ((i0s+nx_s) .LT. nx_srch) .AND. &
                        ((j0s+ny_s) .LT. ny_srch) ) THEN

                        allocate( c_srchimg(nx_s,ny_s) )
                        
                        c_refimg = c_refimg_map(i0m:i0m+nx_m-1, j0m:j0m+ny_m-1 )
                        
                        c_srchimg = c_srchimg_map(i0s:i0s+nx_s-1, j0s:j0s+ny_s-1 )
                        call ampcor_flex(c_refimg,c_srchimg,nx_search,ny_search, &
                            xoff_out,yoff_out,r_snr,r_cov, &
                            nx_m,ny_m,nx_s,ny_s)
                        
                        xoff_out_map(i_off,j_off) = xoff_out+x_shift
                        yoff_out_map(i_off,j_off) = yoff_out+y_shift
                        r_snr_map(i_off,j_off) = r_snr
                        r_cov_map1(i_off,j_off) = r_cov(1)
                        r_cov_map2(i_off,j_off) = r_cov(2)
                        r_cov_map3(i_off,j_off) = r_cov(3) 
                        deallocate(c_srchimg)
                        r_snr=0
                        xoff_out=0
                        yoff_out=0
                        r_cov(1:3)=0
                    ENDIF
                ENDIF
                !}}}
            ENDIF
            ! }}}
            END DO
            ! }}}
        ENDIF
        ! }}}
    ENDIF
    END DO

END

SUBROUTINE ampcor_flex(c_refimg,c_srchimg,nx_search,ny_search, &
        xoff_out,yoff_out,r_snr,r_cov, &
        nx_m,ny_m,nx_s,ny_s) ! hidden

    ! c_refimg = image complex (nx_m,ny_m)
    ! c_srchimg = image complex (nx_m+4*nx_search, ny_m+4*ny_search)
    ! ny_m = size x c_refimg
    ! nx_m = size y c_refimg
    ! ny_s = size x c_srchimg
    ! nx_s = size y c_srchimg 
   
    ! INPUT VARIABLES :
    !------------------

    integer ny_m,nx_m,ny_s,nx_s,nx_search,ny_search
    complex(4) c_refimg,c_srchimg
    dimension c_refimg(nx_m,ny_m)
    dimension c_srchimg(nx_s,ny_s)
    integer i_index,i_indexi

    ! LOCAL VARIABLE :
    !-----------------

    integer*4 i_srchp,nx_sp,ny_sp,i_n2wsxi,i_n2wsyi
    integer*4 i_n2wsxj,i_n2wsyj,i_xx,i_yy
    integer*4 i_covs,i_cw 
    integer*4 i_idx,i_idy 
    real*8 r_covth,r_snrth
    integer*4 i_ovs,i_covsm,i_cwm,i_srchpp
    parameter(i_idx=512) ! was 512
    parameter(i_idy=512) ! was 512
    parameter(i_ovs=2)
    parameter(i_srchpp=4)
    parameter(i_covsm=64)
    parameter(i_cwm=16)

    integer i_sinc_fourier,i_sinc,i_fourier
    parameter(i_sinc=1,i_fourier=2) !i_dump_images=1 means dump debug feature is on
    parameter(i_sinc_fourier=i_sinc)

    integer i_sinc_window
    parameter(i_sinc_window=2)

    integer*4 MAXINTKERLGH      ! maximum interpolation kernel length
    parameter (MAXINTKERLGH=256)

    integer*4 MAXDECFACTOR      ! maximum lags in interpolation kernels
    parameter(MAXDECFACTOR=4096)

    integer*4 MAXINTLGH         ! maximum interpolation kernel array size
    parameter (MAXINTLGH=MAXINTKERLGH*MAXDECFACTOR)
    
    integer*4 i_decfactor       ! Range migration decimation Factor
    integer i_intplength      ! Range migration interpolation kernel length
    real*8  r_fdelay          ! Range migration filter delay
    real*8 r_fintp(0:MAXINTLGH) ! interpolation kernel values
    real*8 r_relfiltlen,r_pedestal,r_beta
    real*4 r_imgi,r_imgj,r_imgc
    real*8 r_snr,r_outside
    dimension r_imgi(i_idx,i_idy)
    dimension r_imgj(i_idx,i_idy)
    dimension r_imgc(i_idx,i_idy)
    complex(4) c_chipref(i_idx*i_idy)
    complex(4) c_chipsch(i_idx*i_idy)
    real*8    r_corr(i_covsm*i_cwm,i_covsm*i_cwm)
    complex(4) c_corr(i_covsm*i_cwm*i_covsm*i_cwm)
    complex(4) c_corrt(i_cwm*i_cwm)
    complex(4) c_ossch(i_ovs*i_idx*i_ovs*i_idy)
    complex(4) c_osref(i_ovs*i_idx*i_ovs*i_idy)

    integer nx_mos,ny_mos,nx_sos,ny_sos,i_wsox,i_wsoy
    integer i_ovss
    integer*4 i_nn(2)
    integer i_dem,i_dir,i_shiftxos,i_shiftyos
    integer i_cpeak(2),i_px,i_py
    real*8 r_max,r_oscoroff(2)

    ! FUNCTION STATEMENT:
    !--------------------

    integer nextpower
    ! OUTPUT VARIABLE:
    !-----------------

    integer i_shiftx,i_shifty
    real*8 r_shfty,r_shftx,r_peak,r_meani,r_meanj
    real*8 r_stdvi,r_stdvj,r_noise,r_cov(3),r_eval1
    real*8 r_eval2,r_evec1(2),r_evec2(2)
    integer i_flag,i_edge(2)
    real*8 r_peakos,r_shftxos,r_shftyos,r_covos(3),r_snros
    real*8 xoff_out,yoff_out,r_mean_cor

    i_covs = 32
    i_cw = 16
    
    ! sinc interploation kernel

    i_decfactor = 4096
    i_weight = 1
    r_pedestal = 0.0
    r_beta = .75
    r_relfiltlen = 6.0

    !call fill_sinc(r_beta,r_relfiltlen,i_decfactor,i_weight,r_pedestal,i_intplength,r_fdelay,r_fintp)

    if (nx_s .NE. nx_m+4*nx_search) THEN
        print*,'nx_s:',nx_s
        print*,'nx_m:',nx_m
        print*,'nx_search:',nx_search
        print*,'nx_m+2*nx_search:',nx_m+4*nx_search
        print*,'ERROR: nx_s =/ nx_m+4*nx_search'
        RETURN
    endif

    if (ny_s .NE. ny_m+4*ny_search) THEN
        print*,'ny_s:',ny_s
        print*,'ny_m:',ny_m
        print*,'ny_search:',ny_search
        print*,'ny_m+2*ny_search:',ny_m+4*ny_search
        print*,'ERROR:ny_s =/ ny_m+4*ny_search'
        RETURN
    endif

    i_srchp = min(ny_search,nx_search,i_srchpp)
    
    nx_sp = nx_m + 2*nx_search
    ny_sp = ny_m + 2*ny_search

    i_n2wsxi = 2**(nextpower(nx_m))
    i_n2wsyi = 2**(nextpower(ny_m))

    i_n2wsxj = 2**(nextpower(nx_sp))
    i_n2wsyj = 2**(nextpower(ny_sp))
    
    r_snrth = 0.0
    r_covth = 1.e10
    r_covth = min(r_covth,999.999998)
    ! get the reference image and search images
    
    r_imgi(1:nx_m,1:ny_m) = cabs(c_refimg(1:nx_m,1:ny_m))
    r_imgj(1:nx_sp,1:ny_sp) = cabs(c_srchimg(1+nx_search:nx_sp+nx_search,1+ny_search:ny_sp+ny_search))

    call correlate(r_imgi,r_imgj,nx_m,ny_m,nx_sp, &
        ny_sp,1,r_meani,r_stdvi,r_meanj, &
        r_stdvj,r_peak,r_noise,r_cov,r_eval1, &
        r_eval2,r_evec1,r_evec2,r_imgc,i_shiftx,i_shifty,i_edge, &
        i_flag)

     r_shftx = float(i_shiftx) - nx_search
     r_shfty = float(i_shifty) - ny_search
    !decide with points are good matches and print out the match values
    IF (i_flag .eq. 0 .and. i_edge(1) .eq. 0 .and. i_edge(2) .eq. 0) THEN  !{{{
        !found a potentially good data point 
        !compute the "snr" {{{
        r_outside = 0.0
        i_cnta = 0
        DO l=max(i_shifty-9,1),min(i_shifty+11,ny_s-ny_m)
            DO k=max(i_shiftx-9,1),min(i_shiftx+11,nx_s-nx_m)
                i_cnta = i_cnta + 1
                r_outside = r_outside + r_imgc(k,l)**2
            ENDDO
        ENDDO
        r_outside = r_outside - r_peak**2
        r_outside = r_outside/(i_cnta-1)

        r_snr = r_peak**2/max(r_outside,1.e-10)
        ! }}}

        IF (r_snr .gt. r_snrth .and. r_cov(1) .lt. r_covth .and. r_cov(2) .lt. r_covth)THEN !{{{
        ! oversample the region around the peak 2 to 1 to estimate the fractional offset {{{

            ! write the reference image and search image around the peak into arrays {{{
            DO i_yy=1,ny_m                        
                DO i_xx=1,nx_m
                    i_index = (i_yy-1)*i_n2wsxi + i_xx
                    IF (i_xx .ge. 1 .and. i_xx .le. nx_m) THEN
                        c_chipref(i_index) = c_refimg(i_xx,i_yy)
                    ELSE
                        c_chipref(i_index) = cmplx(0.0,0.0)
                    ENDIF
                enddo
            enddo
                  
            DO i_yy=1,ny_m
                DO i_xx=nx_m+1,i_n2wsxi
                    i_index = (i_yy-1)*i_n2wsxi + i_xx
                    c_chipref(i_index) = cmplx(0.0,0.0)
                enddo
            enddo
                  
            DO i_yy=ny_m+1,i_n2wsyi
                DO i_xx=1,i_n2wsxi
                    i_index = (i_yy-1)*i_n2wsxi + i_xx
                    c_chipref(i_index) = cmplx(0.0,0.0)
                ENDDO
            ENDDO
                  
           !     now the search image

            DO i_yy=1,ny_sp
                do i_xx=1,nx_sp
                    i_index = (i_yy-1)*i_n2wsxj + i_xx
                    if(i_xx+r_shftx+nx_search .ge. 1 .and. &
                        i_xx+r_shftx+nx_search .le. nx_s .and. &
                        i_yy+r_shfty+ny_search .ge. 1 .and. &
                        i_yy+r_shfty+ny_search .le. ny_s) THEN
                        c_chipsch(i_index) = c_srchimg(NINT(i_xx+r_shftx+nx_search), &
                            NINT(i_yy+r_shfty+ny_search))
                    else
                        c_chipsch(i_index) = cmplx(0.0,0.0)
                    endif
                enddo
            enddo
                  
            do i_yy=1,ny_sp
                do i_xx=nx_sp+1,i_n2wsxj
                    i_index = (i_yy-1)*i_n2wsxj + i_xx
                    c_chipsch(i_index) = cmplx(0.0,0.0)
                enddo
            enddo
                  
            do i_yy=ny_sp+1,i_n2wsyj                        
                do i_xx=1,i_n2wsxj
                    i_index = (i_yy-1)*i_n2wsxj + i_xx
                    c_chipsch(i_index) = cmplx(0.0,0.0)
                enddo
            enddo
            ! }}}

            !Deramp data prior to FFT     
            call derampc(c_chipref,i_n2wsxi,i_n2wsyi)
            call derampc(c_chipsch,i_n2wsxj,i_n2wsyj)
                  
            !forward fft the data
            i_nn(1) = i_n2wsxj
            i_nn(2) = i_n2wsyj

            i_dem = 2
            i_dir = 1
            
            call fourn(c_chipsch,i_nn,i_dir)
            
            i_nn(1) = i_n2wsxi
            i_nn(2) = i_n2wsyi

            call fourn(c_chipref,i_nn,i_dir)

            !spread the spectral data out for inverse transforms
            
            i_nn(1) = i_n2wsxi*i_ovs
            i_nn(2) = i_n2wsyi*i_ovs

            i_dem = 2
            i_dir = -1

            do k=1,i_nn(2)                        
                do l=1,i_nn(1)
                i_index = l + (k-1)*i_nn(1)
                c_osref(i_index) = cmplx(0.0,0.0)
                enddo
            enddo

            do k=1,i_n2wsyi/2                  
                do l=1,i_n2wsxi/2
                    i_index = (k-1)*i_nn(1) + l
                    i_indexi = (k-1)*i_n2wsxi + l
                    c_osref(i_index) = c_chipref(i_indexi)
                    i_index = (i_nn(2) - i_n2wsyi/2 + k - 1)*i_nn(1) + l
                    i_indexi = (k + i_n2wsyi/2 - 1)*i_n2wsxi + l
                    c_osref(i_index) = c_chipref(i_indexi)
                    i_index = (k-1)*i_nn(1) + i_nn(1) - i_n2wsxi/2 + l
                    i_indexi = (k-1)*i_n2wsxi + i_n2wsxi/2 + l
                    c_osref(i_index) = c_chipref(i_indexi)
                    i_index = (i_nn(2) - i_n2wsyi/2 + k - 1)*i_nn(1) + i_nn(1) - i_n2wsxi/2 + l
                    i_indexi = (k + i_n2wsyi/2 - 1)*i_n2wsxi + l + i_n2wsxi/2
                    c_osref(i_index) = c_chipref(i_indexi)
                enddo
            enddo
            
            call fourn(c_osref,i_nn,i_dir)

            i_nn(1) = i_n2wsxj*i_ovs
            i_nn(2) = i_n2wsyj*i_ovs
            i_dem = 2
            i_dir = -1

            do l=1,i_nn(1)
                do k=1,i_nn(2)
                i_index = l + (k-1)*i_nn(1)
                c_ossch(i_index) = cmplx(0.0,0.0)
                enddo
            enddo
                  
            do k=1,i_n2wsyj/2                  
            do l=1,i_n2wsxj/2
            i_index = (k-1)*i_nn(1) + l
            i_indexi = (k-1)*i_n2wsxj + l
            c_ossch(i_index) = c_chipsch(i_indexi)
            i_index = (i_nn(2) - i_n2wsyj/2 + k - 1)*i_nn(1) + l
            i_indexi = (k + i_n2wsyj/2 - 1)*i_n2wsxj + l
            c_ossch(i_index) = c_chipsch(i_indexi)
            i_index = (k-1)*i_nn(1) + i_nn(1) - i_n2wsxj/2 + l
            i_indexi = (k-1)*i_n2wsxj + i_n2wsxj/2 + l
            c_ossch(i_index) = c_chipsch(i_indexi)
            i_index = (i_nn(2) - i_n2wsyj/2 + k - 1)*i_nn(1) + i_nn(1) - i_n2wsxj/2 + l
            i_indexi = (k + i_n2wsyj/2 - 1)*i_n2wsxj + l + i_n2wsxj/2
            c_ossch(i_index) = c_chipsch(i_indexi)
            enddo
            enddo

            !     inverse transform
                  
            call fourn(c_ossch,i_nn,i_dir)
            
                  
            !     detect images and put into correlation arrays
                  
            do i_yy=1,ny_m*i_ovs
                     do i_xx=1,nx_m*i_ovs
                        i_index = i_xx + (i_yy-1)*i_n2wsxi*i_ovs
                        r_imgi(i_xx,i_yy) = cabs(c_osref(i_index)/(i_n2wsxi*i_n2wsyi))
                     enddo
            enddo
                  
            do i_yy=1,ny_sp*i_ovs                        
                     do i_xx=1,nx_sp*i_ovs
                        i_index = i_xx + (i_yy-1)*i_n2wsxj*i_ovs
                        r_imgj(i_xx,i_yy) = cabs(c_ossch(i_index))/(i_n2wsxj*i_n2wsyj)
                enddo
            enddo
            ! }}}
                  
            !     correlate the oversampled chips {{{
            nx_mos = nx_m*i_ovs
            ny_mos = ny_m*i_ovs
            nx_sos = nx_sp*i_ovs 
            ny_sos = ny_sp*i_ovs 
            i_wsox = nx_sos - (nx_mos-1)
            i_wsoy = ny_sos - (ny_mos-1)

            i_ovss = 1

        
            call correlate(r_imgi,r_imgj,nx_mos,ny_mos, &
                nx_sos,ny_sos,i_ovss,r_meani,r_stdvi, &
                r_meanj,r_stdvj,r_peakos, &
                r_noise,r_covos,r_eval1,r_eval2,r_evec1,r_evec2, &
                r_imgc,i_shiftxos,i_shiftyos,i_edge,i_flag)

            r_shftxos = float(i_shiftxos)/i_ovs - float((i_wsox-1)/2)/i_ovs + r_shftx
            r_shftyos = float(i_shiftyos)/i_ovs - float((i_wsoy-1)/2)/i_ovs + r_shfty

            r_outside = 0.0
            i_cnta = 0
            do l=max(i_shiftyos-9,1),min(i_shiftyos+11,i_wsoy)
                do k=max(i_shiftxos-9,1),min(i_shiftxos+11,i_wsox)
                    i_cnta = i_cnta + 1
                    r_outside = r_outside + r_imgc(k,l)**2
                enddo
            enddo

            r_outside = r_outside - r_peakos**2
            r_outside = r_outside/(i_cnta-1)
            r_snros = r_peakos**2/min(r_outside,1.e10)

            r_snros = 10.
            r_covos(1) = 0. 
            r_covos(2) = 0. 
            ! }}}

            if(r_snros .gt. r_snrth .and. r_covos(1) .lt. r_covth .and. r_covos(2) .lt. r_covth) THEN ! {{{
                !print *,'oversample the oversampled correlation surface'
                ! oversample the oversampled correlation surface {{{

                     r_mean_cor = 0.0
                     i_cnta = 0
                     i_px = i_shiftxos+1
                     i_py = i_shiftyos+1
                     r_max = 0.0
                     !print *,'i_px,i_py',i_px,i_py
                     !print *,'maxval ',maxval(r_imgc)
                     !print *,'maxloc ',maxloc(r_imgc)
                     do i_yy=-i_cw/2,i_cw/2-1
                         
                        do i_xx=-i_cw/2,i_cw/2-1
                           
                           i_index = (i_yy+i_cw/2)*i_cw + i_xx + i_cw/2 + 1
                            
                           if (i_xx+i_px .ge. 1 .and. i_xx+i_px .le. nx_mos .and. &
                               i_yy+i_py .ge. 1 .and. i_yy+i_py .le. ny_mos) then

                               c_corrt(i_index) = cmplx(real(abs(r_imgc(i_xx+i_px,i_yy+i_py)/r_peakos)),real(0.))
                               r_mean_cor = r_mean_cor + cabs(c_corrt(i_index))
                               i_cnta = i_cnta + 1
                           else
                               c_corrt(i_index) = cmplx(0.0, 0.0)
                           endif

                        enddo
                        
                     enddo

                !     substract off the mean

                     r_mean_cor = r_mean_cor/max(i_cnta,1)
                     r_mean_cor = 0.0
                     do i_yy=-i_cw/2,i_cw/2-1
                        do i_xx=-i_cw/2,i_cw/2-1
                           i_index = (i_yy+i_cw/2)*i_cw + i_xx + i_cw/2 + 1
                           c_corrt(i_index) = c_corrt(i_index) - cmplx(real(r_mean_cor),0.0)
                        enddo
                     enddo

                !     oversample the correlation surface
                     
                        !     oversample via Fourier transforms {{{
                        !     forward fft the data
                       !print*,'oversample via Fourier transforms' 
                        i_nn(1) = i_cw
                        i_nn(2) = i_cw
                        i_dem = 2
                        i_dir = 1
                        
                        call fourn(c_corrt,i_nn,i_dir)
                        
                        !     spread the spectral data out for inverse transforms
                        
                        i_nn(1) = i_cw*i_covs
                        i_nn(2) = i_cw*i_covs
                        i_dem = 2
                        i_dir = -1
                        
                        do k=1,i_nn(2)                           
                           do l=1,i_nn(1)
                              i_index = (k-1)*i_nn(1) + l
                              c_corr(i_index) = 0.0
                           enddo
                        enddo
                        
                        do l=1,i_cw/2
                           do k=1,i_cw/2
                              i_index = (k-1)*i_nn(1) + l
                              i_indexi = (k-1)*i_cw + l
                              c_corr(i_index) = c_corrt(i_indexi) 
                              i_index = l + (i_nn(2)-i_cw/2+k-1)*i_nn(1)
                              i_indexi = l + (k+i_cw/2-1)*i_cw
                              c_corr(i_index) = c_corrt(i_indexi) 
                              i_index = i_nn(1)-i_cw/2+l + (k-1)*i_nn(2)
                              i_indexi = l+i_cw/2 + (k-1)*i_cw
                              c_corr(i_index) = c_corrt(i_indexi) 
                              i_index = i_nn(1)-i_cw/2+l + (i_nn(2)-i_cw/2+k-1)*i_nn(1)
                              i_indexi = l+i_cw/2 + (k+i_cw/2-1)*i_cw
                              c_corr(i_index) = c_corrt(i_indexi) 
                           enddo
                        enddo
                        
                        !     inverse transform
                        
                        call fourn(c_corr,i_nn,i_dir)
                        
                       ! }}} 
                    ! }}} 
                    !     detect the peak {{{
                     
                     r_max=0.
                     i_cpeak(1)=0
                     i_cpeak(2)=0
                     do i_yy=1,i_cw*i_covs
                        do i_xx=1,i_cw*i_covs
                           i_index = (i_yy-1)*i_cw*i_covs + i_xx
                           r_corr(i_xx,i_yy) = cabs(c_corr(i_index))/((i_cw**2)*(i_cw*i_covs)**2)
                           
                           if (abs(i_xx-i_cw*i_covs/2) .le. i_covs .and. &
                               abs(i_yy-i_cw*i_covs/2) .le. i_covs) then
                               if (r_corr(i_xx,i_yy) .ge. r_max) then
                                   r_max = r_corr(i_xx,i_yy)
                                   i_cpeak(1) = i_xx - i_cw/2*i_covs
                                   i_cpeak(2) = i_yy - i_cw/2*i_covs
                               endif
                           endif
                           enddo
                     enddo
                    
                     r_oscoroff(1) = float(i_cpeak(1)-1)/float(i_covs) 
                     r_oscoroff(2) = float(i_cpeak(2)-1)/float(i_covs)

                     xoff_out = r_oscoroff(1)/i_ovs + r_shftxos
                     yoff_out = r_oscoroff(2)/i_ovs + r_shftyos
                     !print*,'xoff_out,yoff_out:',xoff_out,yoff_out
                    ! }}} 
                     
            else 
                     xoff_out = 0
                     yoff_out = 0

                     ! Bad match at level 2
                     
            endif         !thresholds second pass }}}
                  
        else 
                  
                  ! Bad match at level 1
                  
        endif            !thresholds }}}
               
    endif               !not edge point or no data point }}}

    RETURN    
end

      subroutine correlate(r_imgi,r_imgj,nx_m,ny_m,nx_s, & ! {{{
              ny_s,i_ovs,r_meani,r_stdvi,r_meanj,r_stdvj, &
              r_peak,r_noise,r_cov,r_eval1,r_eval2, &
              r_evec1,r_evec2,r_imgc,i_shftx,i_shfty,i_edge,i_flag)

          !****************************************************************
!**   
!**   FILE NAME: correlate.f
!**   
!**   DATE WRITTEN: /10/10/92
!**   
!**   PROGRAMMER:Scott Hensley / Scott Shaffer
!**   
!**   FUNCTIONAL DESCRIPTION: This routine will do amplitude correlation
!     
!**   on two specified input files.
!**   
!**   ROUTINES CALLED:none
!**   
!**   NOTES: none
!**   
!**   UPDATE LOG:
!**   
!**   Date      Description                              Person
!**   ----      -----------                              ------
!**   /12/12/94   Modified to work with real data.         SH
!**   /02/22/95   Modified to work oversampled data.      SS/SH
!**   
!*****************************************************************
      
      implicit none
      
!     INPUT VARIABLES:
      integer*4 i_idx,i_idy
      parameter(i_idx=512)
      parameter(i_idy=512)

      integer ny_m,nx_m,ny_s,nx_s
      integer i_ovs
      integer i_wsayi,i_wsaxi
      integer i_wsayj,i_wsaxj,i_wsaxyi,i_wsaxyj

      real*4 r_imi,r_imgi,r_imgc
      real*4 r_imj,r_imgj
      dimension r_imi(i_idx,i_idy)
      dimension r_imgc(i_idx,i_idy)
      dimension r_imj(i_idx,i_idy)
      dimension r_imgi(i_idx,i_idy)
      dimension r_imgj(i_idx,i_idy)      

!     OUTPUT VARIABLES:
      real*8 r_shfty,r_shftx,r_peak,r_shrp,r_meani,r_meanj
      real*8 r_stdvi,r_stdvj,r_noise,r_cov(3),r_eval1
      real*8 r_eval2,r_evec1(2),r_evec2(2)
      
!     LOCAL VARIABLES:
      integer i,j,m,n,ix,iy,i_shfty,i_shftx,io
      integer i_cnti,i_cntj,i_cntai,i_cntaj,i_edge(2),i_flag

      real*8 r_sumc,r_sumi,r_smqi
      real*8 r_sumj(0:i_idx,0:i_idy)
      real*8 r_smqj(0:i_idx,0:i_idy)
      real*8 r_crpd(0:i_idx,0:i_idy)
      real*8 r_corr(0:i_idx,0:i_idy)
      real*8 r_corn(0:i_idx,0:i_idy)
      real*8 r_denom

      real*8 r_dxx,r_dyy,r_dxy,r_n2,r_n4,r_u,r_u2

      logical l_init
      
!     DATA STATEMENTS:
      data l_init /.false./
      
!     FUNCTION STATEMENTS:
      
!     PROCESSING STEPS:

      i_edge(1)=0
      i_edge(2)=0
      i_wsayi=ny_m
      i_wsaxi=nx_m
      i_wsayj=ny_s
      i_wsaxj=nx_s
      i_wsaxyi=i_wsayi*i_wsaxi
      i_wsaxyj=i_wsayj*i_wsaxj/i_ovs

      r_cov(1)=0.
      r_cov(2)=0.
      r_cov(3)=0. 

!     compute mean and standard deviations on blocks 

      i_cntai = 0
      i_cntaj = 0
      r_sumi = 0.d0
      r_smqi = 0.d0

      do iy=1,i_wsayj
         do ix=1,i_wsaxj
            r_imgc(ix,iy) = 0.
            r_imi(ix,iy) = 0.
            r_imj(ix,iy) = 0.
            i_cnti=0
            i_cntj=0
            
            r_imj(ix,iy) = r_imgj(ix,iy)
            if(ix .le. nx_m .and. iy .le. ny_m)then
                r_imi(ix,iy) = r_imgi(ix,iy)
                if(r_imi(ix,iy) .ne. 0)then
                    i_cntai = i_cntai+1
                    r_sumi = r_sumi + dble(r_imi(ix,iy))
                    r_smqi = r_smqi + dble(r_imi(ix,iy))**2
                endif
            endif
            if(r_imj(ix,iy) .ne. 0)then
                i_cntaj = i_cntaj+1
            endif   
         enddo
      enddo


      if ( i_cntai .ne. 0 ) then
         r_meani = r_sumi/i_cntai
         r_stdvi = sqrt((r_smqi/i_cntai)-r_meani**2)
      else
         r_meani = 0.
      endif
      if (i_cntai .ge. 0.9*i_wsaxyi .and. &
          i_cntaj .ge. 0.9*i_wsaxyj ) then !have enough real estate

         do iy=0,i_wsayj-1
            r_sumj(0,iy) = 0.
            r_smqj(0,iy) = 0.
            do io = 1,i_ovs
               r_sumj(io,iy) = 0.
               r_smqj(io,iy) = 0.
               do ix=0,(i_wsaxi-1)*i_ovs,i_ovs
                  r_sumj(io,iy) = r_sumj(io,iy) + r_imj(ix+io,iy+1)
                  r_smqj(io,iy) = r_smqj(io,iy) + r_imj(ix+io,iy+1)**2
               enddo
            enddo

            do ix=i_ovs+1,i_wsaxj - (i_wsaxi-1)*i_ovs
               r_sumj(ix,iy) = r_sumj(ix-i_ovs,iy) - r_imj(ix-i_ovs,iy+1 &
                   ) +r_imj(ix+(i_wsaxi-1)*i_ovs,iy+1)
               r_smqj(ix,iy) = r_smqj(ix-i_ovs,iy) - r_imj(ix-i_ovs,iy+1 &
                   )**2 +r_imj(ix+(i_wsaxi-1)*i_ovs,iy+1)**2
            enddo
         enddo

         do ix=0,i_wsaxj - (i_wsaxi-1)*i_ovs-1
            do io=1,i_ovs
               r_sumj(ix,io-1)=0.
               r_smqj(ix,io-1)=0.
               do iy=0,(i_wsayi-1)*i_ovs,i_ovs
                  r_sumj(ix,io-1) = r_sumj(ix,io-1)+r_sumj(ix+1,iy+io-1)
                  r_smqj(ix,io-1) = r_smqj(ix,io-1)+r_smqj(ix+1,iy+io-1)
               enddo
            enddo

            do iy=i_ovs,i_wsayj - (i_wsayi-1)*i_ovs-1
               r_sumj(ix,iy) = r_sumj(ix,iy-i_ovs) - r_sumj(ix+1,iy &
                   -i_ovs)+r_sumj(ix+1,iy+(i_wsayi-1)*i_ovs)
               r_smqj(ix,iy) = r_smqj(ix,iy-i_ovs) - r_smqj(ix+1,iy &
                   -i_ovs)+r_smqj(ix+1,iy+(i_wsayi-1)*i_ovs)
            enddo
         enddo

!         type *,' '
!         do ix=0,i_wsaxj - (i_wsaxi-1)*i_ovs-1
!            do iy=0,i_wsayj - (i_wsayi-1)*i_ovs-1
!               r_sum=0.
!               do ixx=ix+1,ix+i_wsaxi*i_ovs,i_ovs
!                  do iyy=iy+1,iy+i_wsayi*i_ovs,i_ovs
!                     r_sum=r_sum+r_imj(ixx,iyy)
!                  enddo
!               enddo
!               type *,ix,iy,r_sumj(ix,iy),r_sum,r_sumj(ix,iy)-r_sum
!            enddo
!         enddo

         i_shftx = 0
         i_shfty = 0
         r_peak = -9.e27
         do m=0,i_wsaxj - (i_wsaxi-1)*i_ovs-1
            do n=0,i_wsayj - (i_wsayi-1)*i_ovs-1
               r_sumc = 0.
               do j=1,i_wsayi
                  do i=1,i_wsaxi
                     r_sumc = r_sumc + r_imi(i,j)*r_imj((i-1)*i_ovs+m+1,(j-1)*i_ovs+n+1)
                  enddo
               enddo
               r_crpd(m,n) = r_sumc
               r_corr(m,n) = r_sumc - r_meani*r_sumj(m,n)
               r_denom = (r_stdvi*sqrt((r_smqj(m,n)*i_wsaxyi)- &
                   (r_sumj(m,n))**2))
               if ( r_denom .gt. 0. ) then
                  r_corn(m,n) = r_corr(m,n)/r_denom
               else
                  r_corn(m,n) = 0.
               endif
               r_imgc(m+1,n+1) = real(r_corn(m,n))
!               if(nx_m .eq. 112)then
!                  type*, 'r_c = ',m,n,r_corn(m,n),r_crpd(m,n),r_meani*r_sumj(m,n),
!     +                 r_crpd(m,n)-r_meani*r_sumj(m,n),r_sumj(m,n),r_denom
!               endif
               if ( r_peak .lt. r_corn(m,n)) then
                  r_peak = r_corn(m,n)
                  i_shftx = m
                  i_shfty = n
               endif
            enddo
         enddo
!     commpute the curvature of the corrrelation surface to estimate the
!     goodness of the match

         if ( r_peak .gt. 0. ) then

             ix = i_shftx
             iy = i_shfty
             if ( iy .eq. 0 .or. iy .eq. i_wsayj - (i_wsayi-1)*i_ovs-1 ) &
                 i_edge(1)=1
             if ( ix .eq. 0 .or. ix .eq. i_wsaxj - (i_wsaxi-1)*i_ovs-1 ) &
                 i_edge(2)=1
             r_shftx = float(ix)/i_ovs
             r_shfty = float(iy)/i_ovs
             r_meanj = r_sumj(ix,iy)/i_wsaxyi
             r_stdvj = sqrt((r_smqj(ix,iy)/i_wsaxyi)-r_meanj**2)
             r_shrp = (r_peak-(r_corn(max(ix-1,1),iy)+ &
                 r_corn(min(ix+1,i_wsaxj - (i_wsaxi-1)*i_ovs-1),iy))/2.)
             i_flag = 0

             if ( ix .eq. 0 ) then
                 if ( iy .eq. 0 ) then
                     r_dxx = -(r_corn(ix+1,iy)+r_corn(ix+1,iy)- &
                         2*r_corn(ix,iy))
                     r_dyy = -(r_corn(ix,iy+1)+r_corn(ix,iy+1)- &
                         2*r_corn(ix,iy))
                     r_dxy = 0.
                     r_dxx = r_dxx/4 ! added emperically
                     r_dyy = r_dyy/4
                     r_dxy = r_dxy/4
                     r_peak = r_peak/4
                 else if ( iy .eq. i_wsayj - (i_wsayi-1)*i_ovs-1 ) then
                     r_dxx = -(r_corn(ix+1,iy)+r_corn(ix+1,iy)- &
                         2*r_corn(ix,iy))
                     r_dyy = -(r_corn(ix,iy-1)+r_corn(ix,iy-1)- &
                         2*r_corn(ix,iy))
                     r_dxy = 0
                     r_dxx = r_dxx/4 ! added emperically
                     r_dyy = r_dyy/4
                     r_dxy = r_dxy/4
                     r_peak = r_peak/4
                 else
                     r_dxx = -(r_corn(ix+1,iy)+r_corn(ix+1,iy)- &
                         2*r_corn(ix,iy))
                     r_dyy = -(r_corn(ix,iy+1)+r_corn(ix,iy-1)- &
                         2*r_corn(ix,iy))
                     r_dxy = 2*(r_corn(ix+1,iy+1)- &
                         r_corn(ix+1,iy-1))/4
                     r_dxx = r_dxx/2 ! added emperically
                     r_dyy = r_dyy/2
                     r_dxy = r_dxy/2
                     r_peak = r_peak/2
                 endif
             else if ( ix .eq. i_wsaxj - (i_wsaxi-1)*i_ovs-1 ) then
                 if ( iy .eq. 0 ) then
                     r_dxx = -(r_corn(ix-1,iy)+r_corn(ix-1,iy)- &
                         2*r_corn(ix,iy))
                     r_dyy = -(r_corn(ix,iy+1)+r_corn(ix,iy+1)- &
                         2*r_corn(ix,iy))
                     r_dxy = 0
                     r_dxx = r_dxx/4 ! added emperically
                     r_dyy = r_dyy/4
                     r_dxy = r_dxy/4
                     r_peak = r_peak/4
                 else if ( iy .eq. i_wsayj - (i_wsayi-1)*i_ovs-1 ) then
                     r_dxx = -(r_corn(ix-1,iy)+r_corn(ix-1,iy)- &
                         2*r_corn(ix,iy))
                     r_dyy = -(r_corn(ix,iy-1)+r_corn(ix,iy-1)- &
                         2*r_corn(ix,iy))
                     r_dxy = 0
                     r_dxx = r_dxx/4 ! added emperically
                     r_dyy = r_dyy/4
                     r_dxy = r_dxy/4
                     r_peak = r_peak/4
                 else
                     r_dxx = -(r_corn(ix-1,iy)+r_corn(ix-1,iy)- &
                         2*r_corn(ix,iy))
                     r_dyy = -(r_corn(ix,iy+1)+r_corn(ix,iy-1)- &
                         2*r_corn(ix,iy))
                     r_dxy = 2*(r_corn(ix-1,iy-1)- &
                         r_corn(ix-1,iy+1))/4
                     r_dxx = r_dxx/2 ! added emperically
                     r_dyy = r_dyy/2
                     r_dxy = r_dxy/2
                     r_peak = r_peak/2
                 endif
             else if ( iy .eq. 0 ) then
                 r_dxx = -(r_corn(ix+1,iy)+r_corn(ix-1,iy)- &
                     2*r_corn(ix,iy))
                 r_dyy = -(r_corn(ix,iy+1)+r_corn(ix,iy+1)- &
                     2*r_corn(ix,iy))
                 r_dxy = 2*(r_corn(ix+1,iy+1)- &
                     r_corn(ix-1,iy+1))/4
                 r_dxx = r_dxx/2  ! added emperically
                 r_dyy = r_dyy/2
                 r_dxy = r_dxy/2
                 r_peak = r_peak/2
             else if ( iy .eq. i_wsayj - (i_wsayi-1)*i_ovs-1 ) then
                 r_dxx = -(r_corn(ix+1,iy)+r_corn(ix-1,iy)- &
                     2*r_corn(ix,iy))
                 r_dyy = -(r_corn(ix,iy-1)+r_corn(ix,iy-1)- &
                     2*r_corn(ix,iy))
                 r_dxy = 2*(r_corn(ix-1,iy-1)- &
                     r_corn(ix+1,iy-1))/4
                 r_dxx = r_dxx/2  ! added emperically
                 r_dyy = r_dyy/2
                 r_dxy = r_dxy/2
                 r_peak = r_peak/2
             else
                 r_dxx = -(r_corn(ix+1,iy)+r_corn(ix-1,iy)- & 
                     2*r_corn(ix,iy))
                 r_dyy = -(r_corn(ix,iy+1)+r_corn(ix,iy-1)- &
                     2*r_corn(ix,iy))
                 r_dxy = (r_corn(ix+1,iy+1)+ &
                     r_corn(ix-1,iy-1)-r_corn(ix+1,iy-1)- &
                     r_corn(ix-1,iy+1))/4
             endif

            r_n2 = max(1.-r_peak,0.e0)
            r_noise = sqrt(r_n2)
            r_dxx = r_dxx*i_wsaxyi
            r_dyy = r_dyy*i_wsaxyi
            r_dxy = r_dxy*i_wsaxyi

            r_n4 = r_n2**2
            r_n2 = r_n2*2
            r_n4 = r_n4*.5*i_wsaxyi

            r_u = r_dxy**2-r_dxx*r_dyy
            r_u2 = r_u**2       !                    *i_avgx*i_avgy/i_wsaxyi
            if ( r_u .eq. 0 ) then
               r_cov(1)=99.
               r_cov(2)=99.
               r_cov(3)=0.
               i_flag=1
            else
                r_cov(1)=(-r_n2*r_u*r_dyy+r_n4*(r_dyy**2+r_dxy**2)) &
                    /r_u2
                r_cov(2)=(-r_n2*r_u*r_dxx+r_n4*(r_dxx**2+r_dxy**2)) &
                    /r_u2
                r_cov(3)=((r_n2*r_u      -r_n4*(r_dxx+r_dyy))*r_dxy) &
                    /r_u2
            endif
            r_u=sqrt((r_cov(1)+r_cov(2))**2.-4.*(r_cov(1)*r_cov(2)- &
                r_cov(3)**2))
            r_eval1=(r_cov(1)+r_cov(2)+r_u)/2.
            r_eval2=(r_cov(1)+r_cov(2)-r_u)/2.
            if ( r_eval1 .le. 0 .or. r_eval2 .le. 0 ) then
            endif
            
            if ( r_cov(3) .eq. 0 ) then
               if ( r_cov(1) .ge. r_cov(2) ) then
                  r_evec1(1)=1.
                  r_evec1(2)=0.
                  r_evec2(1)=0.
                  r_evec2(2)=1.
               else
                  r_evec1(1)=0.
                  r_evec1(2)=1.
                  r_evec2(1)=1.
                  r_evec2(2)=0.
               endif
            else
               if ( r_cov(1)-r_eval1 .ne. 0. ) then
                  r_evec1(1)=-r_cov(3)/(r_cov(1)-r_eval1)
               else
                  !write(6,*) 'e vector 1 error'
                  r_evec1(1)=999.
               endif
               r_evec1(2)=1.
               r_u=sqrt(r_evec1(1)**2+r_evec1(2)**2)
               r_evec1(1)=r_evec1(1)/r_u
               r_evec1(2)=r_evec1(2)/r_u

               if ( r_cov(1)-r_eval2 .ne. 0. ) then
                  r_evec2(1)=-r_cov(3)/(r_cov(1)-r_eval2)
               else
                  !write(6,*) 'e vector 2 error'
                  r_evec2(1)=999.
               endif
               r_evec2(2)=1.
               r_u=sqrt(r_evec2(1)**2+r_evec2(2)**2)
               r_evec2(1)=r_evec2(1)/r_u
               r_evec2(2)=r_evec2(2)/r_u
            endif

            r_evec1(1)=r_evec1(1)*sqrt(abs(r_eval1)) 
            r_evec1(2)=r_evec1(2)*sqrt(abs(r_eval1)) 
            r_evec2(1)=r_evec2(1)*sqrt(abs(r_eval2)) 
            r_evec2(2)=r_evec2(2)*sqrt(abs(r_eval2)) 

         else

            r_shfty=0
            r_shftx=0
            r_shrp=0.
            i_flag=1
            !write(6,*) 'correlation error'

         endif

      else

         r_shfty=0
         r_shftx=0
         r_shrp=0.
         i_flag=1

      endif


      return

      end ! }}}
      subroutine derampc(c_img,i_dimx,i_dimy)  ! {{{
 
      implicit none
      integer i_dimx,i_dimy,i,j
      complex(4) c_img(i_dimx,i_dimy),c_phdn,c_phac
      real*8 r_phac,r_phdn
 
      c_phdn = cmplx(0.,0.)
      c_phac = cmplx(0.,0.)
 
      do i=1,i_dimx-1
         do j=1,i_dimy
            c_phac = c_phac + c_img(i,j)*conjg(c_img(i+1,j))
         enddo
      enddo
 
      do i=1,i_dimx
         do j=1,i_dimy-1
            c_phdn = c_phdn + c_img(i,j)*conjg(c_img(i,j+1))
         enddo
      enddo
 
      if(cabs(c_phdn) .eq. 0)then
         r_phdn = 0.0
      else
         r_phdn = atan2(aimag(c_phdn),real(c_phdn))
      endif
 
      if(cabs(c_phac) .eq. 0)then
         r_phac = 0.0
      else
         r_phac = atan2(aimag(c_phac),real(c_phac))
      endif
 
!       write(6,*) 'Phase across, down = ',r_phac,r_phdn
      
      do i=1,i_dimx
         do j=1,i_dimy
            c_img(i,j) = c_img(i,j)*cmplx( real(cos(r_phac*i+r_phdn*j)), &
                real(sin(r_phac*i+r_phdn*j)))
         enddo
      enddo
 
      end ! }}}
      subroutine fourn(data,nn,isign) ! {{{

      complex(4) data(*), d(16384)
      integer*4 nn(2)
      integer n
      integer isign,is

      is = -isign
      n = nn(1)
      do i = 1,nn(2)
         call cfft1d_jpl(nn(1),data(1+nn(1)*(i-1)),is)
      end do

      do i = 1,nn(1)

         do j = 1,nn(2)
            d(j) = data(i+nn(1)*(j-1))
         end do
         call cfft1d_jpl(nn(2),d,is)

         do j = 1 , nn(2)
            if(is .eq. 1)then
               d(j) = d(j)*nn(1)*nn(2)
           endif
           data(i+nn(1)*(j-1)) = d(j)
         end do

      end do

      return
      end ! }}}
      integer function nextpower(i_num) ! {{{

!****************************************************************
!**     
!**   FILE NAME: nextpower.f
!**     
!**   DATE WRITTEN: 6/1/97
!**     
!**   PROGRAMMER: Scott Hensley
!**     
!**   FUNCTIONAL DESCRIPTION: Computes the closest number which is a 
!**   power of two and returns the exponent of two for the number that 
!**   is the first power of two exceeding the input number.
!**     
!**   ROUTINES CALLED:
!**     
!**   NOTES: 
!**     
!**   UPDATE LOG:
!**
!**   Date Changed        Reason Changed                  CR # and Version #
!**   ------------       ----------------                 -----------------
!**     
!*****************************************************************

      implicit none

!     INCLUDE FILES:

!     PARAMETER STATEMENTS:

!     INPUT VARIABLES:

      integer i_num
        
!     OUTPUT VARIABLES:

!     LOCAL VARIABLES:

      real*8 r_num,r_log2,r_log2numm1

!     COMMON BLOCKS:

!     EQUIVALENCE STATEMENTS:

!     DATA STATEMENTS:

      data r_log2 /.301029995664d0/

!     FUNCTION STATEMENTS:

!     SAVE STATEMENTS:

      save r_log2

!     PROCESSING STEPS:

      r_num = i_num

      r_log2numm1 = dlog10(r_num - .5d0)/r_log2

      nextpower = int(r_log2numm1)+1
        
      end ! }}} 
      subroutine fill_sinc(r_beta,r_relfiltlen,i_decfactor,i_weight, & ! {{{
              r_pedestal,i_intplength,r_fdelay,r_fintp)

!****************************************************************
!**     
!**   FILE NAME: fill_sinc.f
!**     
!**   DATE WRITTEN: 2/2/98
!**     
!**   PROGRAMMER: Scott Hensley
!**     
!**   FUNCTIONAL DESCRIPTION: This routine computes the sinc interpolation
!**   coefficients needed by the processor for various range and azimuth
!**   interpolations.
!**     
!**   ROUTINES CALLED:
!**     
!**   NOTES: 
!**     
!**   UPDATE LOG:
!**
!**   Date Changed        Reason Changed                  CR # and Version #
!**   ------------       ----------------                 -----------------
!**     
!*****************************************************************

      implicit none

!     INCLUDE FILES:

!     PARAMETER STATEMENTS:
      
      integer*4 MAXDECFACTOR      ! maximum lags in interpolation kernels
      parameter(MAXDECFACTOR=4096)                        
      
      integer*4 MAXINTKERLGH      ! maximum interpolation kernel length
      parameter (MAXINTKERLGH=256)
      
      integer*4 MAXINTLGH         ! maximum interpolation kernel array size
      parameter (MAXINTLGH=MAXINTKERLGH*MAXDECFACTOR)

!     INPUT VARIABLES:

      integer*4 i_decfactor,i_weight
      real*8 r_beta,r_relfiltlen,r_pedestal
        
!     OUTPUT VARIABLES:

      integer i_intplength      ! Range migration interpolation kernel length
      real*8  r_fdelay          ! Range migration filter delay
      real*8 r_fintp(0:MAXINTLGH) ! interpolation kernel values

!     LOCAL VARIABLES:

      real*8 r_filter(0:MAXINTLGH)
      integer*4 i,j,i_filtercoef

!     COMMON BLOCKS:

!     EQUIVALENCE STATEMENTS:

!     DATA STATEMENTS:

!     FUNCTION STATEMENTS:

!     SAVE STATEMENTS:

!     PROCESSING STEPS:

!     get sinc 
      
      call sinc_coef(r_beta,r_relfiltlen,i_decfactor,r_pedestal, &
          i_weight,i_intplength,i_filtercoef,r_filter(0))
      
      r_fdelay = i_intplength/2.0
      
      do i = 0 , i_intplength - 1
         do j = 0 , i_decfactor - 1
            r_fintp(i+j*i_intplength) = real(r_filter(j+i*i_decfactor))
         enddo
      enddo
      
      end ! }}} 
      subroutine sinc_coef(r_beta,r_relfiltlen,i_decfactor,r_pedestal, & ! {{{
          i_weight,i_intplength,i_filtercoef,r_filter)

!****************************************************************
!**     
!**   FILE NAME: sinc_coef.f
!**     
!**   DATE WRITTEN: 10/15/97
!**     
!**   PROGRAMMER: Scott Hensley
!**     
!**   FUNCTIONAL DESCRIPTION: The number of data values in the array 
!**   will always be the interpolation length * the decimation factor, 
!**   so this is not returned separately by the function.
!**     
!**   ROUTINES CALLED:
!**     
!**   NOTES: 
!**     
!**   UPDATE LOG:
!**
!**   Date Changed        Reason Changed                  CR # and Version #
!**   ------------       ----------------                 -----------------
!**     
!*****************************************************************

      implicit none

!     INPUT VARIABLES:

      real*8 r_beta             !the "beta" for the filter
      real*8 r_relfiltlen       !relative filter length
      integer i_decfactor       !the decimation factor
      real*8 r_pedestal         !pedestal height
      integer i_weight          !0 = no weight , 1=weight
        
!     OUTPUT VARIABLES:
      
      integer i_intplength      !the interpolation length
      integer*4 i_filtercoef      !number of coefficients
      real*8 r_filter(*)        !an array of data values 

!     LOCAL VARIABLES:

      real*8 pi,r_wgt,r_s,r_fct,r_wgthgt,r_soff,r_wa
      integer i

!     COMMON BLOCKS:

!     EQUIVALENCE STATEMENTS:

!     DATA STATEMENTS:

!     FUNCTION STATEMENTS:

!     PROCESSING STEPS:

      pi = 4.d0*atan(1.d0)

!     number of coefficients

      i_intplength = nint(r_relfiltlen/r_beta)
      i_filtercoef = i_intplength*i_decfactor
      r_wgthgt = (1.d0 - r_pedestal)/2.d0
      r_soff = (i_filtercoef - 1.d0)/2.d0
      
      do i=0,i_filtercoef-1
         r_wa = i - r_soff
         r_wgt = (1.d0 - r_wgthgt) + r_wgthgt*cos((pi*r_wa)/r_soff)
         r_s = r_wa*r_beta/dble(i_decfactor)
         if(r_s .ne. 0.0)then
            r_fct = sin(pi*r_s)/(pi*r_s)
         else
            r_fct = 1.0
         endif
         if(i_weight .eq. 1)then
            r_filter(i+1) = r_fct*r_wgt
         else
            r_filter(i+1) = r_fct
         endif
      enddo
      
      end ! }}}
      subroutine cfft1d_jpl(n,c,dir) ! {{{

    integer*4  n, dir, ier
    complex*8    c(*)
    integer nmax
    parameter (nmax = 32768)

#include <fftw3.f>

    integer*8 plani(16),planf(16)
    complex*8 in(nmax)
    integer i
    save plani,planf
    
    if(dir.eq.0)then
         do i=3,15
            if(2**i.eq.n)go to 1
         end do
         write(6,*) 'fftw: length unsupported:: ',n
         stop
 1       do j= 1 , n
            in(j) = cmplx(0.,0.)
         end do
         call sfftw_plan_dft_1d(planf(i),n,in,in,FFTW_FORWARD,FFTW_MEASURE)
         call sfftw_plan_dft_1d(plani(i),n,in,in,FFTW_BACKWARD,FFTW_MEASURE)
         return    
    endif

    if(dir.eq.-1)then
        if(n.eq.8)call sfftw_execute_dft(planf(3),c,c)
        if(n.eq.16)call sfftw_execute_dft(planf(4),c,c)
        if(n.eq.32)call sfftw_execute_dft(planf(5),c,c)
        if(n.eq.64)call sfftw_execute_dft(planf(6),c,c)
        if(n.eq.128)call sfftw_execute_dft(planf(7),c,c)
        if(n.eq.256)call sfftw_execute_dft(planf(8),c,c)
        if(n.eq.512)call sfftw_execute_dft(planf(9),c,c)
        if(n.eq.1024)call sfftw_execute_dft(planf(10),c,c)
        if(n.eq.2048)call sfftw_execute_dft(planf(11),c,c)
        if(n.eq.4096)call sfftw_execute_dft(planf(12),c,c)
        if(n.eq.8192)call sfftw_execute_dft(planf(13),c,c)
        if(n.eq.16384)call sfftw_execute_dft(planf(14),c,c)
        if(n.eq.32768)call sfftw_execute_dft(planf(15),c,c)
    end if
    if(dir.eq. 1)then
        if(n.eq.8)call sfftw_execute_dft(plani(3),c,c)
        if(n.eq.16)call sfftw_execute_dft(plani(4),c,c)
        if(n.eq.32)call sfftw_execute_dft(plani(5),c,c)
        if(n.eq.64)call sfftw_execute_dft(plani(6),c,c)
        if(n.eq.128)call sfftw_execute_dft(plani(7),c,c)
        if(n.eq.256)call sfftw_execute_dft(plani(8),c,c)
        if(n.eq.512)call sfftw_execute_dft(plani(9),c,c)
        if(n.eq.1024)call sfftw_execute_dft(plani(10),c,c)
        if(n.eq.2048)call sfftw_execute_dft(plani(11),c,c)
        if(n.eq.4096)call sfftw_execute_dft(plani(12),c,c)
        if(n.eq.8192)call sfftw_execute_dft(plani(13),c,c)
        if(n.eq.16384)call sfftw_execute_dft(plani(14),c,c)
        if(n.eq.32768)call sfftw_execute_dft(plani(15),c,c)
    end if       


    !  if(dir.eq.-1)then
    !      call sfftw_plan_dft_1d(plan,n,c,cout,FFTW_FORWARD,FFTW_ESTIMATE)
    !  end if
    !  if(dir.eq. 1)then
    !      call sfftw_plan_dft_1d(plan,n,c,cout,FFTW_BACKWARD,FFTW_ESTIMATE)
    !  end if
    !  call sfftw_execute_dft(plan,c,cout)
    !  call sfftw_destroy_plan(plan)

     ! c=cout

      return
      end
      ! }}}
