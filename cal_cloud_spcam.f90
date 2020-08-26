program cldsize
    use netcdf
    implicit none
!    include "netcdf.inc"
    Integer             , parameter :: NX=64 , NY=1, NZ=24 , NT=1488
    real(kind=4), parameter         :: dx=4000. ![m]
    Integer             , parameter :: GRIDSIZE   = NX*NZ
    Character ( len=* ) , parameter :: FILE_NAME  = "Q_120-20_0001-01.nc"
    Character ( len=* ) , parameter :: SIZE_FILE  = "Q_120-20_0001-01.txt"
    Logical             , parameter :: WRITE_BACK = .TRUE.
    Integer             , parameter :: CONNECT    = 4
    Integer             , parameter :: C_DIM      = 2
    Integer                         :: Time, x, label, i, j, temp, k, t, num
    Integer                         :: label_data(NX,NZ) , cloudflag(NX,NZ)
    Real(kind=4)                    :: cld_inc_size( GRIDSIZE )
    real(kind=4), dimension(NX,NZ,NT) :: qcalldata, qralldata
    Real(kind=4), dimension(NX,NZ)    :: qc3d_data, qr3d_data, size_data
    real(kind=4), dimension(NZ)       :: height, depth
    integer(kind=1), dimension(NX,NZ) :: cltype
    integer :: ncid, timeid, zid, xid, qid, rid, tsize, labid, sizid, hid, typid
    integer :: latid, lonid, status
    real(kind=4) :: lon, lat, hsfc
    character(2), dimension(GRIDSIZE) :: cld_one_type

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
       call cloud_type( label_data, qr3d_data, height, depth, hsfc, cld_inc_size, cltype, cld_one_type )

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

subroutine cloud_type( label_data, qr3d_data, height, depth, hsfc, cld_inc_size, cltype, cld_one_type )
integer, dimension(NX,NZ), intent(in) :: label_data
real(kind=4), dimension(NX,NZ), intent(in) :: qr3d_data
real(kind=4), dimension(NZ), intent(in) :: height, depth
real(kind=4), intent(in) :: hsfc, cld_inc_size(GRIDSIZE)
integer(kind=1), dimension(NX,NZ), intent(out) :: cltype
character(2), dimension(GRIDSIZE), intent(out) :: cld_one_type
character(2), dimension(8), parameter :: typename =  &
                                    (/"Hc","As","Ac","St","Sc","Cu","Ns","Dc"/)
integer, dimension(8,GRIDSIZE) :: typecount
integer :: i, k, n, maxn, typenum, k1, k2, m
real :: base, top, clsize, rain, thick
cltype = 0
typecount = 0
do i=1,NX
   maxn = maxval(label_data(i,:))
   do n=1,maxn
      typenum = 10
      k1=0
      k2=0
      do k=1,NZ
         if (label_data(i,k)==n) then
            base = height(k) - hsfc
            rain = qr3d_data(i,k)
            k1 = k
            exit
         endif
      enddo
      do k=NZ,k1,-1
         if (label_data(i,k)==n) then
            top = height(k)
            k2 = k
            exit
         endif
      enddo
      if(k1==0 .and. k2 ==0) cycle
      clsize = cld_inc_size(n)
      thick = 0.
      do k=k1,k2
         thick = thick + depth(k)
      enddo

      if (base > 7000.) then
         typenum = 1
      elseif (base > 4000.) then
         if (rain > 1.e-6) then
            typenum = 3
         else
            typenum = 2
         endif
      elseif ((thick > 5000.) .and. (rain > 1.e-5)) then
         if (clsize*1000./thick > 100.) then
            typenum = 7
         else
            typenum = 8
         endif
      elseif (base <= 2000.) then
         if (clsize*1000./thick > 30.) then
            if (rain >= 1.e-5) then
               typenum = 5
            else
               typenum = 4
            endif
         else
            typenum = 6
         endif
      else
         if ((base > 3000.) .or. (clsize*1000./thick > 100.)) then
            if (rain > 1.e-6) then
               typenum = 3
            else
               typenum = 2
            endif
         else
            typenum = 6
         endif
      endif

      do k=1,NZ
         if (label_data(i,k) == n) cltype(i,k) = typenum
      enddo
      typecount(typenum,n) = typecount(typenum,n) + 1
   enddo
enddo
cld_one_type = "--"
maxn = maxval(label_data)
do n=1,maxn
   cld_one_type(n) = typename( maxloc(typecount(:,n), dim=1) )
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
character(*), parameter :: topofile = "/data/W.eddie/SPCAM/USGS-gtopo30_1.9x2.5_remap_c050602.nc"
integer :: toponcid, phisid, lonid, latid, xidx, yidx, i
real(kind=8) :: phis, lonall(144), latall(96)
call check_nf90( nf90_open(topofile, NF90_NOWRITE, toponcid) )
call check_nf90( nf90_inq_varid(toponcid, "lon", lonid) )
call check_nf90( nf90_get_var(toponcid, lonid, lonall) )
call check_nf90( nf90_inq_varid(toponcid, "lat", latid) )
call check_nf90( nf90_get_var(toponcid, latid, latall) )
xidx = 0
do i=1,144
  if(abs(real(lon,kind=8)-lonall(i)) < 1.e-4 ) then
    xidx = i
    exit
  endif
enddo
yidx = 0
do i=1,96
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
