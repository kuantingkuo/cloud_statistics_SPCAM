program cldsize
    use netcdf
    implicit none
!    include "netcdf.inc"
    Integer             , parameter :: NX=64 , NY=1, NZ=28 , NT=720
    real(kind=4), parameter         :: dx=2000. ![m]
    Integer             , parameter :: GRIDSIZE   = NX*NZ
!    Character ( len=* ) , parameter :: FILE_NAME  = "Q_120-20_0001-01.nc"
!    Character ( len=* ) , parameter :: SIZE_FILE  = "Q_120-20_0001-01.txt"
    Logical             , parameter :: WRITE_BACK = .True.
    Integer             , parameter :: CONNECT    = 4 ![4|8] ways connetion
    Integer             , parameter :: C_DIM      = 2
    Character ( len=: ), allocatable :: FILE_NAME, SIZE_FILE
    Integer                         :: Time, x, label, i, j, temp, k, t, num
    Integer                         :: label_data(NX,NZ) , cloudflag(NX,NZ)
    Real(kind=4)                    :: cld_inc_size( GRIDSIZE )
    real(kind=4), dimension(NX,NZ,NT) :: qcalldata, qralldata
    Real(kind=4), dimension(NX,NZ)    :: qc3d_data, qr3d_data, size_data
    real(kind=4), dimension(NZ)       :: height, depth
    integer(kind=1), dimension(NX,NZ) :: cltype
    integer :: ncid, timeid, zid, xid, qid, rid, tsize, labid, sizid, hid, typid
    integer :: latid, lonid, status, length
    real(kind=4) :: lon, lat, hsfc
    character(2), dimension(GRIDSIZE) :: cld_one_type
    logical, parameter :: debug = .false.

    call get_command_argument(1, length=length)
    allocate(character(length) :: FILE_NAME)
    call get_command_argument(1, value=FILE_NAME)
    length = length + 1
    allocate(character(length) :: SIZE_FILE)
    SIZE_FILE = FILE_NAME(1:length-3)//"txt"
    print*, "output file to ", trim(SIZE_FILE)
    open( unit=4567 , file=SIZE_FILE, status="unknown" )

!    open( unit=6789 , file=FILE_NAME , status="old"  , &
!          access="direct" , form='unformatted' , &
!          recl = NX  )

    print*, "read data..."

    if (WRITE_BACK) then
      call check_nf90( nf90_open(file_name, NF90_WRITE, ncid) )
    else
      call check_nf90( nf90_open(file_name, NF90_NOWRITE, ncid) )
    endif
    call check_nf90( nf90_inq_dimid(ncid, "time", timeid) )
    call check_nf90( nf90_inq_dimid(ncid, "crm_nz", zid) )
    call check_nf90( nf90_inq_dimid(ncid, "crm_nx", xid) )
    call check_nf90( nf90_inquire_dimension(ncid, timeid, len=tsize) )
    call check_nf90( nf90_inq_varid(ncid, "qci", qid) )
    call check_nf90( nf90_get_var(ncid, qid, qcalldata(:,:,1:tsize)) )
    call check_nf90( nf90_inq_varid(ncid, "qpr", rid) )
    call check_nf90( nf90_get_var(ncid, rid, qralldata(:,:,1:tsize)) )
    call check_nf90( nf90_inq_varid(ncid, "crm_nz", hid) )
    call check_nf90( nf90_get_var(ncid, hid, height) )
    call check_nf90( nf90_inq_varid(ncid, "lat", latid) )
    call check_nf90( nf90_get_var(ncid, latid, lat) )
    call check_nf90( nf90_inq_varid(ncid, "lon", lonid) )
    call check_nf90( nf90_get_var(ncid, lonid, lon) )

    call surface_height( lon, lat, hsfc )

    do k=2,NZ-1
       depth(k) = (height(k+1) - height(k-1))/2.
    enddo
    depth(1) = height(1)-hsfc + (height(2)-height(1))/2.
    depth(NZ) = height(NZ) - height(NZ-1)

    if (WRITE_BACK) then
      print*, "preparing output variables..."

      status = nf90_def_var(ncid, 'label', NF90_UBYTE, (/xid, zid, timeid/), labid)
      if (status /= NF90_ENAMEINUSE) then
        call check_nf90( nf90_def_var(ncid, 'size', NF90_FLOAT, &
                                      (/xid, zid, timeid/), sizid) )
        call check_nf90( nf90_def_var(ncid, 'type', NF90_UBYTE, &
                                      (/xid, zid, timeid/), typid) )
        call check_nf90( nf90_put_att(ncid, labid, 'long_name', &
                                                   'number of cloud') )
        call check_nf90( nf90_put_att(ncid, labid, 'units', '#') )
        call check_nf90( nf90_put_att(ncid, sizid, 'long_name', 'cloud size') )
        call check_nf90( nf90_put_att(ncid, sizid, 'units', 'km^2') )
        call check_nf90( nf90_put_att(ncid, typid, 'long_name', &
             'cloud types (1:High cloud 2:Altostratus 3:Altocumulus 4:Stratus 5:Stratocumulus 6:Cumulus 7:Nimbostratus 8:Deep convective clouds 10:undefined)') )
        call check_nf90( nf90_enddef(ncid) )
      else
        call check_nf90( nf90_inq_varid(ncid, 'label', labid) )
        call check_nf90( nf90_inq_varid(ncid, 'size', sizid) )
        call check_nf90( nf90_inq_varid(ncid, 'type', typid) )
      endif
    endif

    print*, "start identifying cloud objects..."
    do Time = 1 , tsize
       qc3d_data = qcalldata(:,:,time)
       call find_cloud( qc3d_data, dx, depth, label_data, cld_inc_size, label )
       cld_inc_size = cld_inc_size * 1.e-6 ! m^2 -> km^2
       qr3d_data = qralldata(:,:,time)
       if(debug) print*,'Time:',time
       call cloud_type( label_data, qc3d_data, qr3d_data, height, depth, hsfc, cltype, cld_one_type )

       if ( WRITE_BACK ) then
           call refill_data( label_data , cld_inc_size , size_data)
!            call check ( nf_put_vara_int( ncid(1) , varcldsize , nc_start , nc_count ,label_data ) ) 
           call check_nf90( nf90_put_var(ncid, labid, label_data, &
                                         (/1, 1, time/), (/NX, NZ, 1/)) )
           call check_nf90( nf90_put_var(ncid, sizid, size_data, &
                                         (/1, 1, time/), (/NX, NZ, 1/)) )
           call check_nf90( nf90_put_var(ncid, typid, cltype, &
                                         (/1, 1, time/), (/NX, NZ, 1/)) )

!           do k=1,NZ
!            do j=1,NY  
!             do i=1,NX        
!                sizedata(Time,i,j,k)=label_data(i,j,k)
!                flagdata(Time,i,j,k)=cloudflag(i,j,k)
!             enddo
!            enddo
!           enddo
     
        endif

        temp = 0
        do x  = 1 , label-1
            if ( cld_inc_size(x) .GE. temp ) then
                temp = cld_inc_size(x)
            endif
            write ( 4567,'(I5,I9,F9.3,1X,A2)' ) Time, x, cld_inc_size(x), cld_one_type(x)
        enddo

        !write ( 4567,'(A2,I5,A6,I9)' ) 'T=' , Time , ',size=' , temp 

    end do

    call check_nf90( nf90_close(ncid) )
    close(4567)

contains
    subroutine refill_data( label_data , cld_inc_size , size_data)
        Integer    , intent (in) :: label_data(NX,NZ)
        real(kind=4), intent(in) :: cld_inc_size(GRIDSIZE)
        real(kind=4), dimension(NX,NZ), intent(out) :: size_data
        Integer                     :: x,z
        size_data(:,:) = 0.
        do x=1,NX
           do z=1,NZ
               if ( label_data(x,z) .EQ. 0 ) cycle
               size_data(x,z) = cld_inc_size( label_data(x,z) )
           enddo
        enddo
    end subroutine

    subroutine find_cloud( qc_data, dx, depth, label_data, cld_inc_size, label )
        real       , intent (in)    :: qc_data(NX,NZ) 
        real(kind=4), intent(in)    :: dx, depth(NZ)
        integer    , intent (out)   :: label_data(NX,NZ)
        real(kind=4), intent(out)   :: cld_inc_size(GRIDSIZE)
        integer    , intent (out)   :: label 
        integer                     :: x , z 
        real(kind=4)                :: inc
        integer                     :: search_flag(NX,NZ)
        integer                     :: stackX(GRIDSIZE) , stackZ(GRIDSIZE)
        integer                     :: move_x_3d(6) =  (/ -1,  0, +1,  0,  0,  0 /),&
                                       move_z_3d(6) =  (/  0,  0,  0,  0, +1, -1 /),&
                                       move_x_2d(8)  = (/ -1,  0, +1,  0, -1, -1, +1, +1 /),&
                                       move_z_2d(8)  = (/  0, -1,  0, +1, -1, +1, -1, +1 /)
        integer                     :: cu_move , nex , nez , myx , myy ,myz 
        integer                     :: stack_ptr
        
        label = 1
        label_data(:,:)  = 0
        search_flag(:,:) = 0
        stack_ptr = 0
        stackX(:) = 0
        stackZ(:) = 0

        do x = 1 , NX
           do z =  1 , NZ
               if ( search_flag(x,z) .EQ. 1 ) cycle

               stackX(1) = x
               stackZ(1) = z
               search_flag(x,z) = 1
               stack_ptr = 1
               inc = 0.

               do while ( stack_ptr .GE. 1 )

                   myx =  stackX(stack_ptr)
                   myz =  stackZ(stack_ptr)

                   stack_ptr = stack_ptr - 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                   if ( qc_data(myx,myz) .LT. 1.e-5 ) cycle !for cloud!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   label_data(myx,myz) = label
                   inc = inc + dx * depth(myz)

                   ! override the current ptr
                   do cu_move = 1 , CONNECT
                       if ( C_DIM .EQ. 2 ) then
                           nex = myx + move_x_2d( cu_move )
                           nez = myz + move_z_2d( cu_move )
                       else
                           nex = myx + move_x_3d( cu_move )
                           nez = myz + move_z_3d( cu_move )
                       endif

                       if ( nex == NX+1 ) nex = 1  ! for periodic boundary
                       if ( nex == 0 ) nex = NX

                       if ( nex .GT. NX  .OR.  nez .GT. NZ ) cycle
                       if ( nex .LT. 1   .OR.  nez .LT. 1  ) cycle
                       if ( search_flag(nex,nez) .EQ. 1 ) cycle
                       stackX(stack_ptr + 1) = nex
                       stackZ(stack_ptr + 1) = nez
                       search_flag(nex,nez) = 1
                       stack_ptr = stack_ptr + 1
                   end do
               end do

               if ( inc .NE. 0 ) then
                   cld_inc_size( label ) = inc
                   label = label + 1
                   inc   = 0
               endif
           enddo
        enddo

    end subroutine

subroutine cloud_type( label_data, qc3d_data, qr3d_data, height, depth, hsfc, cltype, cld_one_type )
!! Algorithm of cloud classification from Cloudsat
!! Ref.: http://irina.eas.gatech.edu/EAS_Fall2008/CloudSat_ATBD_L2_cloud_clas.pdf
integer, dimension(NX,NZ), intent(in) :: label_data
real(kind=4), dimension(NX,NZ), intent(in) :: qc3d_data, qr3d_data
real(kind=4), dimension(NZ), intent(in) :: height, depth
real(kind=4), intent(in) :: hsfc
integer(kind=1), dimension(NX,NZ), intent(out) :: cltype
character(2), dimension(128), intent(out) :: cld_one_type
real(kind=4), dimension(NX,128) :: thick, Ze, base, top, rain, Zeheight
real(kind=4), dimension(128) :: Max10dbH, maxtop, minbase, maxdz, maxZe, &
                                meantop, meanbase, meandz, meanZe, meanHeight, &
                                devtop, devbase, devZe
integer, dimension(128) :: typenum, indprec, intenseprec, veryintenseprec, length
logical, dimension(128) :: deep_flag, conv_flag
character(2), dimension(8), parameter :: typename =  &
                                    (/"Hc","As","Ac","St","Sc","Cu","Ns","Dc"/)
integer :: i, k, n, maxn, k1, k2, m, countn
real :: cloudF, Ztemp
logical, parameter :: debug = .false.
real, dimension(NZ) :: cloudFlev
logical, dimension(NZ) :: mask

maxn = maxval(label_data(:,:))
thick = 0.
cltype = 0
Ze = -9.99e10
Max10dbH = -9.99e10
cloudF = 0.
cloudFlev = 0.
base = 9.99e10
top = -9.99e10
Zeheight = -9.99e10
rain = -9.99e10
maxtop = -9.99e10
minbase = 9.99e10
maxdz = -9.99e10
maxZe = -9.99e10

do i=1,NX
    if (any(label_data(i,:) > 0)) then
        where(label_data(i,:) > 0) cloudFlev(:) = cloudFlev(:) + 1./real(NX)
    else
        cycle
    endif
    do n=1,maxn
        if (.not. any(label_data(i,:)==n)) cycle
        k1=0
        k2=0
        do k=1,NZ
            if (label_data(i,k)==n) then
                base(i,n) = height(k) - hsfc
                k1 = k
                exit
            endif
        enddo
        if (all(label_data(i,1:k1-1)==0) .and. all(qr3d_data(i,1:k1) > 1.e-6)) then
            rain(i,n) = qr3d_data(i,1)
        endif
        do k=NZ,k1,-1
            if (label_data(i,k)==n) then
                top(i,n) = height(k) - hsfc
                k2 = k
                exit
            endif
        enddo
        if(k1==0 .and. k2==0) cycle
        do k=k1,k2
            thick(i,n) = thick(i,n) + depth(k)
            Ztemp = (log10(qr3d_data(i,k)+0.02*qc3d_data(i,k)) + 5.) * 20.
            if (Ze(i,n) < Ztemp) then
                Ze(i,n) = Ztemp
                Zeheight(i,n) = height(k)
            endif
            if (Ztemp > 10. .and. height(k) > Max10dbH(n)) Max10dbH(n) = height(k)
        enddo
   enddo
enddo

where(Ze < -40.) Ze = -40.
where(Ze > 50.) Ze = 50.
maxtop = maxval(top, dim=1)/1000.
minbase = minval(base, dim=1)/1000.
maxdz = maxval(thick, dim=1)/1000.
maxZe = maxval(Ze, dim=1)
Max10dbH = Max10dbH/1000.
indprec = count(rain > 1.e-5, dim=1)
intenseprec = count(rain > 1.e-4, dim=1)
veryintenseprec = count(rain > 1.e-3, dim=1)

meantop = 0.
meanbase = 0.
meandz = 0.
meanZe = 0.
meanHeight = 0.
length = 0
do n=1,maxn
    countn = 0
    do i=1,NX
        if(top(i,n) > 0.) then
            meantop(n) = meantop(n) + top(i,n)
            meanbase(n) = meanbase(n) + base(i,n)
            meandz(n) = meandz(n) + thick(i,n)
            meanZe(n) = meanZe(n) + Ze(i,n)
            meanHeight(n) = meanHeight(n) + Zeheight(i,n)
            countn = countn + 1
        endif
    enddo
    meantop(n) = meantop(n)/real(countn)/1000.
    meanbase(n) = meanbase(n)/real(countn)/1000.
    meandz(n) = meandz(n)/real(countn)/1000.
    meanZe(n) = meanZe(n)/real(countn)
    meanHeight(n) = meanHeight(n)/real(countn)/1000.
    length(n) = countn * 4
enddo

devtop = 0.
devbase = 0.
devZe = 0.
do n=1,maxn
    countn = 0
    do i=1,NX
        if(top(i,n) > 0.) then
            devtop(n) = devtop(n) + (top(i,n) - meantop(n))**2
            devbase(n) = devbase(n) + (base(i,n) - meanbase(n))**2
            devZe(n) = devZe(n) + (Ze(i,n) - meanZe(n))**2
            countn = countn + 1
        endif
    enddo
    devtop(n) = sqrt(devtop(n)/real(countn))/1000.
    devbase(n) = sqrt(devbase(n)/real(countn))/1000.
    devZe(n) = sqrt(devZe(n)/real(countn))
enddo

cld_one_type = "--"
do n=1,maxn
    mask = .false.
    do k=1,NZ
        if(any(label_data(:,k)==n)) mask = .true.
    enddo
    cloudF = maxval(cloudFlev, mask)
    if (indprec(n) > 0 .and. meantop(n) > 2.5) then  ! prec. cloud
        deep_flag(n) = .false.
        conv_flag(n) = .false.
        if ((maxtop(n)>12. .and. Max10dbH(n)>8.2) .or. &
            maxtop(n)>14. .and. meandz(n)>12. .or. &
            meantop(n)>8.5 .and. Max10dbH(n)>8.4) then
            deep_flag(n) = .true.
            conv_flag(n) = .true.
        endif
        if ((intenseprec(n) > 5 .or. maxZe(n) > 14. .or. meanZe(n) > 4.) .and. &
            (length(n)<80. .or. devtop(n)>0.5 .or. meanZe(n) > 4.) .and. &
            (meantop(n)-Max10dbH(n) < 0.34) .and. & 
            Max10dbH(n) > 3. .and. meandz(n) < 5.) then
            conv_flag(n) = .true.
        endif
        if (minbase(n) < 3.) then !? MeanbaseT > -6.
            if (meandz(n) < 2.6 .and. meanbase(n) < 2. .and. maxtop(n) < 4.5) then
                if (maxtop(n) < 4.5 .and. meandz(n) < 2.5) then
                    typenum(n) = 4  ! St
                    if(debug) print*,n,'prec->yes'
                elseif (meantop(n) < 3.5 .and. max10dbH(n) < 3. .and. &
                        (intenseprec(n) < 1 .or. length(n) > 100) .and. &
                        veryintenseprec(n) < 1) then
                    typenum(n) = 5  ! Sc
                    if(debug) print*,n,'prec->yes'
                else
                    typenum(n) = 6  ! Cu
                    if(debug) print*,n,'prec->yes->1'
                endif
            elseif (meandz(n) <= 6. .and. length(n) < 75 .and. &
                    meanZe(n)+devZe(n) >= 6. .and. maxtop(n) <=7 .and. &
                    devtop(n) >= 0.3) then
                typenum(n) = 6  ! Cu
                    if(debug) print*,n,'prec->yes->2'
            elseif ((indprec(n) >= 13 .or. (maxdz(n) < 8. .and. length(n) > 60) .or. &
                    (maxdz(n) > 7. .and. max10dbH(n) < 3.5 .and. length(n) > 30) .or. &
                    (meandz(n) < 10.2 .and. length(n) > 45 .and. meanZe(n) < 10. .and. &
                     max10dbH(n) < 4.2) .and. (meandz(n) < 8.2 .and. &
                     length(n) > 25. .and. meanZe(n) < 10. .and. &
                     max10dbH(n) < 4.2) .and. (meandz(n) < 8.2 .and. &
                     length(n) > 25. .and. meanZe(n) < 10. .and. max10dbH(n) < 2.)) & 
                    .and. (veryintenseprec(n) < 1 .or. (max10dbH(n) < 5. .and. &
                            veryintenseprec(n) > 0)) &
                    .and. (real(indprec(n)*4)/real(length(n)) > 0.3 .or. &
                           (meandz(n) > 5. .and. length(n) >= 240)) &
                    .and. (.not. conv_flag(n) .or. (length(n) > 100 .and. &
                           meanbase(n) < 1.8)) &
                    .and. (.not. deep_flag(n))) then
                typenum(n) = 7  ! Ns
                    if(debug) print*,n,'prec->yes'
            elseif (((meantop(n)-meanheight(n) < 2.1 .and. meanZe(n) < 5. .and. &
                      devtop(n) < 0.3) .or. (meantop(n) > 4. .and. maxZe(n) < 10. .and. &
                      meanZe(n) < -1 .and. maxtop(n) < 7.) .or. (meantop(n) > 4. .and. &
                      maxZe(n) < 7. .and. meanZe(n) < 0. .and. maxtop(n) < 7.) .or. &
                     (meandz(n) < 5. .and. meanZe(n) < 0.5 .and. &
                      devtop(n) < 0.45) .or. (meantop(n) > 3. .and. &
                      meandz(n) > 2.5 .and. meandz(n) < 4.5 .and. devtop(n) < 0.3)) &
                    .and. (.not. conv_flag(n))) then
                typenum(n) = 3  ! Ac
                    if(debug) print*,n,'prec->yes'
            elseif (meandz(n) < 5. .or. (meandz(n) < 6. .and. maxtop(n) < 6.5)) then
                typenum(n) = 6  ! Cu
                    if(debug) print*,n,'prec->yes->3',meandz(n)
            else
                typenum(n) = 8  ! Dc
                    if(debug) print*,n,'prec->yes'
            endif
        else
            if (maxtop(n) < 3.9 .and. meandz(n) < 2.5 .and. meanbase(n) < 1.5 &
                    .and. (length(n) > 50 .or. maxZe(n) < 10.)) then
                typenum(n) = 5  ! Sc
                    if(debug) print*,n,'prec->no'
            elseif (meandz(n) < 2.5 .and. meanbase(n) > 1.8) then
                typenum(n) = 3  ! Ac
                    if(debug) print*,n,'prec->no'
            elseif ((length(n) < 50  .and. meanZe(n)+devZe(n) > 7.) .or. &
                    (length(n) < 70 .and. meanZe(n)+devZe(n) > 12.) .or. &
                    (conv_flag(n) .and. meandz(n) < 3.)) then
                typenum(n) = 6  ! Cu
                    if(debug) print*,n,'prec->no', meanZe(n),devZe(n)
            else
                typenum(n) = 7  ! Ns
                    if(debug) print*,n,'prec->no'
            endif
        endif  ! end prec. cloud

    elseif ((meanbase(n)>5. .and. meanheight(n)>5 .and. meanZe(n)<-3.) .or. &
            meanbase(n)>10.) then  ! high cloud
        if ((meanZe(n) < 0.05 .and. meanheight(n) > 7.5 .and. minbase(n) > 5. .and. &
             meandz(n) < 6.1 .and. meanbase(n) > 5.5) .or. meanbase(n) > 10.) then
            typenum(n) = 1  ! Hc
                    if(debug) print*,n,'high'
        elseif (meanbase(n) > 2.) then
            typenum(n) = 2  ! As
                    if(debug) print*,n,'high',meanbase(n),meanheight(n),meanZe(n)
        elseif (meandz(n) < 6.) then
            typenum(n) = 6  ! Cu
                    if(debug) print*,n,'high'
        else
            typenum(n) = 8  ! Dc
                    if(debug) print*,n,'high'
        endif  ! end high cloud

    elseif (meanHeight(n)<2. .or. meanbase(n)<1.5) then  !low cloud
        if (cloudF < 0.25) then
            typenum(n) = 6  ! Cu
                    if(debug) print*,n,'low->1'
        elseif (maxtop(n) < 3. .and. meanbase(n) < 1.8 .and. intenseprec(n) < 1) then
            typenum(n) = 4  ! St
                    if(debug) print*,n,'low'
        elseif (devZe(n)/meanZe(n) > 0.3 .and. maxtop(n) > 3. .and. &
                maxtop(n) < 9.5 .and. meanZe(n) < 2. .and. meandz(n) < 8. .and. &
                (meanbase(n) > 1.8 .or. (meanbase(n) > 1. .and. maxtop(n) > 3.5) .or. &
                 (meanZe(n) < -5. .and. maxtop(n) > 3.5 .and. meandz(n) > 2.))) then
            typenum(n) = 3  ! Ac
                    if(debug) print*,n,'low'
        elseif (meandz(n) > 8. .and. meanZe(n) < 0.) then
            typenum(n) = 2  ! As
                    if(debug) print*,n,'low'
        elseif ((meandz(n) > 2. .or. intenseprec(n) < 1 .or. maxtop(n) >= 3.) .and. &
                meandz(n) < 7. .and. meanZe(n) > 0. .and. length(n) < 100) then
            typenum(n) = 6  ! Cu
                    if(debug) print*,n,'low'
        elseif ((meandz(n) > 2. .or. maxtop(n) >= 4.) .and. meandz(n) < 7. .and. &
                meanZe(n) > -5 .and. length(n) > 56) then
            typenum(n) = 7  ! Ns
                    if(debug) print*,n,'low'
        elseif (meandz(n) > 8.) then
            typenum(n) = 8  ! Dc
                    if(debug) print*,n,'low'
        else
            typenum(n) = 5  ! Sc
                    if(debug) print*,n,'low',meanbase(n),meanHeight(n)
        endif  ! end low cloud

    else  ! middle cloud
        if (meanbase(n) < 8. .and. ((meanheight(n) < 2.2 .and. maxtop(n) < 3.) .or. &
                    (meanbase(n) < 2. .and. maxtop(n) < 3.) .or. &
                    (meanbase(n) < 1.85 .and. maxtop(n) < 4. .and. &
                     meandz(n) < 2. .and. meantop(n) < 3.2))) then
            typenum(n) = 5  ! Sc
                    if(debug) print*,n,'middle'
        elseif (meanbase(n) < 2. .and. meandz(n) < 1. .and. &
                maxtop(n)-minbase(n) < 1.5) then !? mintop ?
            typenum(n) = 4  ! St
                    if(debug) print*,n,'middle'
        elseif (meanbase(n) < 1.2 .and. meandz(n) > 4.5 .and. length(n) > 52) then
            typenum(n) = 7  ! Ns
                    if(debug) print*,n,'middle'
        elseif (meanbase(n) < 2. .and. (cloudF < 0.25 .or. &
                    (cloudF < 0.5 .and. meanheight(n) < 5.))) then
            typenum(n) = 6  ! Cu
                    if(debug) print*,n,'middle'
        elseif (meanbase(n) < 8. .and. ((meanbase(n) > 1.8 .and. &
                 meanbase(n) < 3.5 .and. meandz(n) < 2.7 .and. meanZe(n) < 0.) &
                .or. (meanbase(n) > 1.8 .and. meanbase(n) < 7. .and. &
                    (meandz(n) < 1.25 .and. meanZe(n) < -10 .and. &
                     meanheight(n) < 7.5) .or. (meandz(n) < 1.75 .and. &
                     meanZe(n) < -10 .and. meanheight(n) < 8.3 .and. &
                     maxtop(n) > 9.5) .or. (meandz(n) < 2.5 .and. meanZe(n) < -8 .and. &
                     maxtop(n) < 7.5)) &
                .or. (meanZe(n)+devZe(n) > 0. .and. meanZe(n)-devZe(n) < -15 .and. &
                      devbase(n) > 0.6 .and. meanbase(n) < 5.5 .and. &
                      meandz(n) < 4. .and. maxtop(n) > 8.7 .and. devtop(n) < 0.6) &
                .or. (devbase(n)/meandz(n) > 0.5 .and. maxtop(n) < 5.) &
                .or. (cloudF < 0.65 .or. meandz(n) < 1. .or. maxtop(n) < 3.8) &
                .or. (((meanbase(n) < 1.8 .and. devbase(n) > 0.6 .and. &
                        meanZe(n) < -5.) .or. (meanbase(n) < 1.2 .and. &
                        devbase(n) > 0.48 .and. meanZe(n) < -1.5)) .and. &
                      meantop(n) > 3. .and. meantop(n) < 8. .and. &
                      (length(n) < 100 .or. meanZe(n) < -15.)) &
                .or. (meanbase(n) > 1.7 .and. meantop(n) < 8. .and. &
                      maxtop(n) < 9.5 .and. (length(n) < 100 .or. meanZe(n) < -3. .or. &
                      (devbase(n) > 1. .and. devtop(n) < 0.3)) .and. &
                      ((minbase(n) < 0.8 .and. devbase(n) > 0.48) .or. &
                       (minbase(n) < 1.2 .and. devbase(n) > 0.57) .or. &
                       (minbase(n) < 1.5 .and. devbase(n) > 0.9))))) then
            typenum(n) = 3  ! Ac
                    if(debug) print*,n,'middle'
        else
            typenum(n) = 2  ! As
                    if(debug) print*,n,'middle'
        endif  ! end middle cloud
    endif
    cld_one_type(n) = typename( typenum(n) )
    cltype = 0
    do k=1,NZ
        do i=1,NX
            if (label_data(i,k) < 1) cycle
            cltype(i,k) = typenum(label_data(i,k))
        enddo
    enddo
enddo
end subroutine

subroutine check_nf90(err)
integer, intent(in) :: err
if (err /= nf90_noerr) then
  print*, "ERROR: ", trim(nf90_strerror(err)), err
  stop
endif
end subroutine

subroutine surface_height( lon, lat, hsfc )
real(kind=4), intent(in) :: lon, lat
real(kind=4), intent(out) :: hsfc
real(kind=4), parameter :: g = 9.80616
character(*), parameter :: topofile = "/data/W.eddie/SPCAM/USGS-gtopo30_0.9x1.25_remap_c051027.nc"
integer, parameter :: nlon=288, nlat=192
integer :: toponcid, phisid, lonid, latid, xidx, yidx, i
real(kind=8) :: phis, lonall(nlon), latall(nlat)
call check_nf90( nf90_open(topofile, NF90_NOWRITE, toponcid) )
call check_nf90( nf90_inq_varid(toponcid, "lon", lonid) )
call check_nf90( nf90_get_var(toponcid, lonid, lonall) )
call check_nf90( nf90_inq_varid(toponcid, "lat", latid) )
call check_nf90( nf90_get_var(toponcid, latid, latall) )
xidx = 0
do i=1,nlon
  if(abs(real(lon,kind=8)-lonall(i)) < 1.e-4 ) then
    xidx = i
    exit
  endif
enddo
yidx = 0
do i=1,nlat
  if(abs(real(lat,kind=8)-latall(i)) < 1.e-4 ) then
    yidx = i
    exit
  endif
enddo
call check_nf90( nf90_inq_varid(toponcid, "PHIS", phisid) )
call check_nf90( nf90_get_var(toponcid, phisid, phis, (/xidx, yidx/)) )
call check_nf90( nf90_close(toponcid) )
hsfc = real(phis,kind=4)/g
end subroutine

end program
