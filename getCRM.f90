program getCRM
use netcdf
implicit none
character(20), parameter :: casename="CPL64"
character(99), parameter :: path="/data/W.eddie/SPCAM/"//trim(casename)//"/"
real(kind=4), parameter :: target_lon=120., target_lat=20.
integer, dimension(12), parameter :: dayend=(/31,28,31,30,31,30,31,31,30,31,30,31/)
integer :: i, xidx, yidx, year, month
integer :: ncid, lonid, latid, z3id
real(kind=4) :: lon(144), lat(96), hgt(24)
logical :: first
character(4) :: yyyy
character(2) :: mm, dd
real(kind=4), dimension(64,24,1488)

call check_nf90( nf90_open(trim(path)//trim(casename)//".cam.h0.0001-01-01-00000.nc", &
                           NF90_NOWRITE, ncid) )
call check_nf90( nf90_inq_varid(ncid, "lon", lonid) )
call check_nf90( nf90_inq_varid(ncid, "lat", latid) )
call check_nf90( nf90_get_var(ncid, lonid, lon ) )
call check_nf90( nf90_get_var(ncid, latid, lat ) )
do i=1,144
   if( abs(lon(i)-target_lon) < abs(lon(i+1)-target_lon) ) then
     xidx = i
     exit
   endif
enddo
do i=1,96
   if( abs(lat(i)-target_lat) < abs(lat(i+1)-target_lat) ) then
     yidx = i
     exit
   endif
enddo

call check_nf90( nf90_inq_varid(ncid, "Z3", z3id) )
call check_nf90( nf90_get_var(ncid, z3id, hgt, (/xidx, yidx, 3/), (/1, 1, 24/)) )
hgt = hgt(24:1:-1)
call check_nf90( nf90_close(ncid) )

write(outlon,'(I3)') nint(target_lon)
write(outlat,'(I3)') nint(target_lat)

first = .True.
do year=1,10
   write(yyyy,'(I0.4)') year
   do month=1,12
      write(mm,'(I0.2)') month
      print*, yyyy,mm
      do day=1,dayend(month)
         write(dd,'(I0.2)') day
         crmt1 = (day-1)*48 + 1
         filename = &
            trim(path)//trim(casename)//".cam.h1."//year//"-"//month//"-"//dd//"-00000.nc"
         call check_nf90( nf90_open(filename, NF90_NOWRITE, ncid) )
         if (first) then
            call check_nf90( nf90_inq_varid(ncid, "LAT_15s_to_30n", latid) )
            call check_nf90( nf90_get_var(ncid, latid, lat(1:24)) )
            call check_nf90( nf90_inq_varid(ncid, "LON_60e_to_180e", lonid) )
            call check_nf90( nf90_get_var(ncid, lonid, lon(1:49)) )
            do i=1,49
               if( abs(lon(i)-target_lon) < abs(lon(i+1)-target_lon) ) then
                 xidx = i
                 exit
               endif
            enddo
            do i=1,24
               if( abs(lat(i)-target_lat) < abs(lat(i+1)-target_lat) ) then
                 yidx = i
                 exit
               endif
            enddo

            first = .False.
         endif
         call check_nf90( nf90_inq_varid(ncid, "CRM_QI_LON_60e_to_180e_LAT_15s_to_30n", &
                                         qtotid) )
         call check_nf90( nf90_get_var(ncid, qtotid, qtot(:,:,crmt1:day*48), &
                                       (/xidx, yidx, 1, 1, 1, 1/), (/1, 1, 64, 1, 24, 48/)) )
         call check_nf90( nf90_inq_varid(ncid, "CRM_QC_LON_60e_to_180e_LAT_15s_to_30n", &
                                         qciid) )
         call check_nf90( nf90_get_var(ncid, qciid, qci(:,:,crmt1:day*48), &
                                       (/xidx, yidx, 1, 1, 1, 1/), (/1, 1, 64, 1, 24, 48/)) )
         call check_nf90( nf90_close(ncid) )
      enddo
      outfile = &
        "/data/W.eddie/cloudsize/"//trim(casename)//"/Q_"//trim(outlon)//"-"//trim(outlat)//"_"//yyyy//"-"//mm//".nc"
      call check_nf90( nf90_create(outfile, NF90_NETCDF4, ncid) )
      call check_nf90( nf90_def_dim(ncid, "crm_nx", 64, crmxid) )
      call check_nf90( nf90_def_dim(ncid, "crm_nz", 24, crmzid) )
      call check_nf90( nf90_def_dim(ncid, "time", dayend(month)*48, timeid) )

   enddo
enddo

contains

subroutine check_nf90(err)
integer, intent(in) :: err
if (err /= nf90_noerr) then
  print*, "ERROR: ", nf90_strerror(err), err
  stop
endif
end subroutine

end program
