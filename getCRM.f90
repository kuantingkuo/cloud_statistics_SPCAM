program getCRM
use netcdf
implicit none
character(20), parameter :: casename="CPL64"
character(99), parameter :: path="/data/W.eddie/SPCAM/CPL64/"
real(kind=4), parameter :: target_lon=98, target_lat=-10 ! target grid
character(99), parameter :: outpath="./"//trim(casename)//"/"
integer, dimension(12), parameter :: dayend=(/31,28,31,30,31,30,31,31,30,31,30,31/)
integer :: i, xidx, yidx, year, month, day, crmt1, tsize, t
integer :: ncid, lonid, latid, z3id, timevid, qtotid, qciid, crmxid, crmzid
integer :: timeid, lat1id, lon1id, crmzvid, crmxvid, qprid
real(kind=4) :: lon(144), lat(96), time(1488), hgt(24)
logical :: first, file_exist
character(3) :: outlon, outlat
character(4) :: yyyy
character(2) :: mm, dd
character(99) :: filename, outfile
real(kind=4), dimension(64,24,1488) :: qci, qtot, qpr
integer(kind=8) :: counttime, crate
call system_clock(count=counttime, count_rate=crate)

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
! read height coordinate for this grid
call check_nf90( nf90_inq_varid(ncid, "Z3", z3id) )
call check_nf90( nf90_get_var(ncid, z3id, hgt, (/xidx, yidx, 3/), (/1, 1, 24/)) )
hgt = hgt(24:1:-1)
call check_nf90( nf90_close(ncid) )

call execute_command_line("mkdir -p "//trim(outpath), wait=.True.)

first = .True.
do year=1,10
   write(yyyy,'(I0.4)') year
   do month=1,12
      write(mm,'(I0.2)') month
      call timer(counttime)
      print*, yyyy,mm
      print*,"open data..."
      do day=1,dayend(month)
         write(dd,'(I0.2)') day
         crmt1 = (day-1)*48 + 1
         filename = &
            trim(path)//trim(casename)//".cam.h1."//yyyy//"-"//mm//"-"//dd//"-00000.nc"
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
            write(outlon,'(I3)') nint(lon(xidx))
            write(outlat,'(I3)') nint(lat(yidx))
            outlon=adjustl(outlon)
            outlat=adjustl(outlat)

            first = .False.
         endif
         call check_nf90( nf90_inq_varid(ncid, "time", timevid) )
         call check_nf90( nf90_inq_varid(ncid, "CRM_QI_LON_60e_to_180e_LAT_15s_to_30n", &
                                         qtotid) )
         call check_nf90( nf90_inq_varid(ncid, "CRM_QC_LON_60e_to_180e_LAT_15s_to_30n", &
                                         qciid) )
         do t=1,48
            call check_nf90( nf90_get_var(ncid, timevid, time(crmt1+t-1), (/t/)) )
            call check_nf90( nf90_get_var(ncid, qtotid, qtot(:,:,crmt1+t-1), &
                                     (/xidx, yidx, 1, 1, 1, t/), (/1, 1, 64, 1, 24, 1/)) )
            call check_nf90( nf90_get_var(ncid, qciid, qci(:,:,crmt1+t-1), &
                                     (/xidx, yidx, 1, 1, 1, t/), (/1, 1, 64, 1, 24, 1/)) )
         enddo
         call check_nf90( nf90_close(ncid) )
      enddo
      print*,"processing dataset..."
      tsize = dayend(month)*48
      qpr(:,:,1:tsize) = qtot(:,:,1:tsize) - qci(:,:,1:tsize)
      outfile = &
        trim(outpath)//"Q_"//trim(outlon)//"-"//trim(outlat)//"_"//yyyy//"-"//mm//".nc"
      INQUIRE(FILE=outfile, EXIST=file_exist)
      if (file_exist) then
         print*, trim(outfile), " exists. Delete it and create new one."
         call execute_command_line("rm -f "//outfile, wait=.True.)
      endif
      call check_nf90( nf90_create(outfile, NF90_NETCDF4, ncid) )
      call check_nf90( nf90_def_dim(ncid, "crm_nx", 64, crmxid) )
      call check_nf90( nf90_def_dim(ncid, "crm_nz", 24, crmzid) )
      call check_nf90( nf90_def_dim(ncid, "time", tsize, timeid) )
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
      call check_nf90( nf90_put_att(ncid, timevid, "units", &
                                    "days since 0001-01-01 00:00:00") )
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
      call check_nf90( nf90_put_var(ncid, timevid, time(1:tsize)) )
      call check_nf90( nf90_put_var(ncid, crmzvid, hgt) )
      call check_nf90( nf90_put_var(ncid, crmxvid, (/ (i*4+2, i=0,63) /)) )
      call check_nf90( nf90_put_var(ncid, qciid, qci(:,:,1:tsize)) )
      call check_nf90( nf90_put_var(ncid, qprid, qpr(:,:,1:tsize)) )
      call check_nf90( nf90_close(ncid) )
      call execute_command_line("./cal_cloud_spcam.exe "//trim(outfile), wait=.False.)
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

subroutine timer(starttime)
integer(kind=8), intent(inout) :: starttime
integer(kind=8) :: ic,crate 
call system_clock(count=ic, count_rate=crate)
print*, real(ic-starttime)/real(crate), "seconds"
starttime = ic
end subroutine

end program
