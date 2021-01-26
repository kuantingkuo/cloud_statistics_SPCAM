program getCRM
use netcdf
implicit none
character(20), parameter :: casename="tc20170601"
character(99), parameter :: path="/data/W.eddie/SPCAM/tc20170601/atm/hist/"
real(kind=4), parameter :: target_lon=95, target_lat=-8.52632 ! target grid
integer, parameter :: nlat=192, nlon=288, nlev=30, nx=64, nz=28, nt=720
character(99), parameter :: outpath="./"//trim(casename)//"/"
integer, dimension(12), parameter :: dom=(/31,28,31,30,31,30,31,31,30,31,30,31/)
integer :: i, xidx, yidx, year, month, day, crmt1, t, dayend, d
integer :: ncid, lonid, latid, z3id, timevid, qtotid, qciid, crmxid, crmzid
integer :: timeid, lat1id, lon1id, crmzvid, crmxvid, qprid, lastind
real(kind=4) :: lon(nlon), lat(nlat), hgt(nz)
real(kind=8) :: time(nt)
logical :: first, file_exist
character(3) :: outlon, outlat
character(4) :: yyyy
character(2) :: mm, dd
character(99) :: filename, outfile, timeunit
real(kind=4), dimension(nx,nz,nt) :: qci, qtot, qpr
integer(kind=8) :: counttime, crate
logical :: leap
call system_clock(count=counttime, count_rate=crate)

lastind = len_trim(casename)
yyyy = casename(lastind-7:lastind-4)
mm = casename(lastind-3:lastind-2)
dd = casename(lastind-1:lastind)
call check_nf90( nf90_open(trim(path)//trim(casename)//".cam.h0."// &
            yyyy//"-"//mm//"-"//dd//"-00000.nc", NF90_NOWRITE, ncid) )
call check_nf90( nf90_inq_varid(ncid, "lon", lonid) )
call check_nf90( nf90_inq_varid(ncid, "lat", latid) )
call check_nf90( nf90_get_var(ncid, lonid, lon ) )
call check_nf90( nf90_get_var(ncid, latid, lat ) )
call check_nf90( nf90_inq_varid(ncid, "time", timeid) )
call check_nf90( nf90_get_att(ncid, timeid, "units", timeunit) )
do i=1,nlon
   if( abs(lon(i)-target_lon) < abs(lon(i+1)-target_lon) ) then
     xidx = i
     exit
   endif
enddo
do i=1,nlat
   if( abs(lat(i)-target_lat) < abs(lat(i+1)-target_lat) ) then
     yidx = i
     exit
   endif
enddo
! read height coordinate for this grid
call check_nf90( nf90_inq_varid(ncid, "Z3", z3id) )
call check_nf90( nf90_get_var(ncid, z3id, hgt, (/xidx, yidx, 3/), (/1, 1, nz/)) )
hgt = hgt(nz:1:-1)
call check_nf90( nf90_close(ncid) )

call execute_command_line("mkdir -p "//trim(outpath), wait=.True.)
read(yyyy,*) year
read(mm,*) month
read(dd,*) day

first = .True.
do d=1,15
   write(yyyy,'(I0.4)') year
   if (mod(year,4)==0 .and. (mod(year,100)/=0 .or. mod(year,400)==0)) then
      leap = .true.
   else
      leap = .false.
   endif
   write(mm,'(I0.2)') month
   dayend = dom(month)
   if (month==2 .and. leap) dayend = 29
   call timer(counttime)
   print*,"open data..."
   write(dd,'(I0.2)') day
   print*, yyyy,mm,dd
   crmt1 = (d-1)*48 + 1
   filename = &
      trim(path)//trim(casename)//".cam.h1."//yyyy//"-"//mm//"-"//dd//"-00000.nc"
   call check_nf90( nf90_open(filename, NF90_NOWRITE, ncid) )
   if (first) then
      call check_nf90( nf90_inq_varid(ncid, "LAT_15s_to_30n", latid) )
      call check_nf90( nf90_get_var(ncid, latid, lat(1:24)) )
      call check_nf90( nf90_inq_varid(ncid, "LON_60e_to_180e", lonid) )
      call check_nf90( nf90_get_var(ncid, lonid, lon(1:49)) )
      do i=1,97
         if( abs(lon(i)-target_lon) < abs(lon(i+1)-target_lon) ) then
           xidx = i
           exit
         endif
      enddo
      do i=1,48
         if( abs(lat(i)-target_lat) < abs(lat(i+1)-target_lat) ) then
           yidx = i
           exit
         endif
      enddo
      write(outlon,'(I3)') nint(lon(xidx))
      write(outlat,'(I3)') nint(lat(yidx))
      outlon=adjustl(outlon)
      outlat=adjustl(outlat)
      first = .False.
   endif
   outfile = &
       trim(outpath)//"Q_"//trim(outlon)//"-"//trim(outlat)//".nc"
   INQUIRE(FILE=outfile, EXIST=file_exist)
   if (file_exist) then
       print*, trim(outfile), " exists."
!       call execute_command_line("rm -f "//outfile, wait=.True.)
       call system("./cal_cloud_spcam2.exe "//trim(outfile))
       stop
   endif

   call check_nf90( nf90_inq_varid(ncid, "time", timevid) )
   call check_nf90( nf90_inq_varid(ncid, "CRM_QI_LON_60e_to_180e_LAT_15s_to_30n", &
                                   qtotid) )
   call check_nf90( nf90_inq_varid(ncid, "CRM_QC_LON_60e_to_180e_LAT_15s_to_30n", &
                                   qciid) )
   do t=1,48
      call check_nf90( nf90_get_var(ncid, timevid, time(crmt1+t-1), (/t/)) )
      call check_nf90( nf90_get_var(ncid, qtotid, qtot(:,:,crmt1+t-1), &
                               (/xidx, yidx, 1, 1, 1, t/), (/1, 1, nx, 1, nz, 1/)) )
      call check_nf90( nf90_get_var(ncid, qciid, qci(:,:,crmt1+t-1), &
                               (/xidx, yidx, 1, 1, 1, t/), (/1, 1, nx, 1, nz, 1/)) )
   enddo
   call check_nf90( nf90_close(ncid) )
    day = day + 1
    if(day > dayend) then
        day = 1
        month = month + 1
    endif
    if(month > 12) then
        month = 1
        year = year + 1
    endif
enddo

print*,"processing dataset..."
qpr(:,:,:) = qtot(:,:,:) - qci(:,:,:)
call check_nf90( nf90_create(outfile, NF90_NETCDF4, ncid) )
call check_nf90( nf90_def_dim(ncid, "crm_nx", nx, crmxid) )
call check_nf90( nf90_def_dim(ncid, "crm_nz", nz, crmzid) )
call check_nf90( nf90_def_dim(ncid, "time", nt, timeid) )
call check_nf90( nf90_def_var(ncid, "lat", NF90_DOUBLE, lat1id) )
call check_nf90( nf90_def_var(ncid, "lon", NF90_DOUBLE, lon1id) )
call check_nf90( nf90_def_var(ncid, "time", NF90_DOUBLE, timeid, timevid) )
call check_nf90( nf90_def_var(ncid, "crm_nz", NF90_FLOAT, crmzid, crmzvid) )
call check_nf90( nf90_def_var(ncid, "crm_nx", NF90_INT64, crmxid, crmxvid) )
call check_nf90( nf90_def_var(ncid, "qci", NF90_FLOAT, &
                              (/crmxid, crmzid, timeid/), qciid) )
call check_nf90( nf90_def_var(ncid, "qpr", NF90_FLOAT, &
                              (/crmxid, crmzid, timeid/), qprid) )
call check_nf90( nf90_put_att(ncid, lat1id, "long_name", "latitude") )
call check_nf90( nf90_put_att(ncid, lat1id, "units", "degrees_north") )
call check_nf90( nf90_put_att(ncid, lon1id, "long_name", "longitude") )
call check_nf90( nf90_put_att(ncid, lon1id, "units", "degrees_east") )
call check_nf90( nf90_put_att(ncid, timevid, "long_name", "time") )
call check_nf90( nf90_put_att(ncid, timevid, "units", trim(timeunit)) )
call check_nf90( nf90_put_att(ncid, timevid, "calendar", "noleap") )
call check_nf90( nf90_put_att(ncid, crmzvid, "long_name", "Height") )
call check_nf90( nf90_put_att(ncid, crmzvid, "units", "m") )
call check_nf90( nf90_put_att(ncid, crmxvid, "units", "km") )
call check_nf90( nf90_put_att(ncid, qciid, "long_name", &
                              "CRM Cloud Water + Cloud Ice") )
call check_nf90( nf90_put_att(ncid, qciid, "units", "kg/kg") )
call check_nf90( nf90_put_att(ncid, qprid, "long_name", &
                              "CRM Precipitating Water (QPC+QPI)") )
call check_nf90( nf90_put_att(ncid, qprid, "units", "kg/kg") )

call check_nf90( nf90_enddef(ncid) )
print*,"output netcdf file..."
call check_nf90( nf90_put_var(ncid, lat1id, lat(yidx)) )
call check_nf90( nf90_put_var(ncid, lon1id, lon(xidx)) )
call check_nf90( nf90_put_var(ncid, timevid, time) )
call check_nf90( nf90_put_var(ncid, crmzvid, hgt) )
call check_nf90( nf90_put_var(ncid, crmxvid, (/ (i*4+2, i=0,63) /)) )
call check_nf90( nf90_put_var(ncid, qciid, qci) )
call check_nf90( nf90_put_var(ncid, qprid, qpr) )
call check_nf90( nf90_close(ncid) )
!call execute_command_line("./cal_cloud_spcam.exe "//trim(outfile), wait=.False.)
call system("./cal_cloud_spcam2.exe "//trim(outfile))

contains

subroutine check_nf90(err)
integer, intent(in) :: err
if (err /= nf90_noerr) then
  print*, "ERROR: ", nf90_strerror(err), err
  stop
endif
end subroutine

subroutine timer(starttime)
integer(kind=8), intent(inout) :: starttime
integer(kind=8) :: ic,crate 
call system_clock(count=ic, count_rate=crate)
print*, real(ic-starttime)/real(crate), "seconds"
starttime = ic
end subroutine

end program
