import xarray as xr
import dask
import numpy as np
import os
import time

starttime=time.time()
# Basic settings 
case="CPL64"
path="/data/W.eddie/SPCAM/"+case+"/"
target_lon=90.  # target longitude
target_lat=0.   # target latitude
outpath="/data/W.eddie/cloudtype/"+casename+"/"

# Read heiget coordinate in this grid
gcmf=xr.open_dataset(path+case+".cam.h0.0001-01-01-00000.nc")
hgt=gcmf.Z3.sel(lat=target_lat, lon=target_lon, method="nearest").values[0,:1:-1]
gcmf.close

for y in range(1, 11):
    year=str(y).zfill(4)
    for m in range(1, 13):
        month=str(m).zfill(2)
        endtime=time.time()
        print(endtime-starttime,"seconds")
        starttime=endtime
        print(year+month)

        print("open data...")
        crmf=xr.open_mfdataset(path+case+".cam.h1."+year+"-"+month+"-*-00000.nc", parallel=True, combine='by_coords')
        print("extract the target grid...")
        qtot=crmf.CRM_QI_LON_60e_to_180e_LAT_15s_to_30n.squeeze().sel(LAT_15s_to_30n=target_lat, LON_60e_to_180e=target_lon, method="nearest")
        qci=crmf.CRM_QC_LON_60e_to_180e_LAT_15s_to_30n.squeeze().sel(LAT_15s_to_30n=target_lat, LON_60e_to_180e=target_lon, method="nearest")
        crmf.close

        lon=str(int(np.round(qtot.LON_60e_to_180e.values)))
        lat=str(int(np.round(qtot.LAT_15s_to_30n.values)))
        outfile=outpath+"Q_"+lon+"-"+lat+"_"+year+"-"+month+".nc"

        if os.path.isfile(outfile):
           if os.stat(outfile).st_size in {18307429, 16536805, 17717221}:
              print("netcdf file allready exist")
              continue

        print("processing dataset...")
        del qci.attrs['mdims'],qci.attrs['basename']

        ds = xr.Dataset(
            {"qci": qci, "qpr": qtot-qci},
            coords={"crm_nz": hgt, "crm_nx": np.arange(0,256,4)},
        ).rename({'LON_60e_to_180e': 'lon', 'LAT_15s_to_30n': 'lat'}).reset_coords()

        ds['qpr'].attrs={'units':'kg/kg', 'long_name':'CRM Precipitating Water (QPC+QPI)'}
        ds['crm_nz'].attrs={'units':'m', 'long_name':'Height'}
        ds['crm_nx'].attrs={'units':'km'}

        print("output netcdf file...")
        ds.to_netcdf(outfile,
                     encoding={'time': {'_FillValue': None},
                               'lat': {'_FillValue': None},
                               'lon': {'_FillValue': None},
                               'crm_nz': {'_FillValue': None},
                               'qci': {'_FillValue': None},
                               'qpr': {'_FillValue': None}
                              })
        ds.close
        os.system("/data/W.eddie/cloudtype/cal_cloud_spcam.exe "+outfile+" &")
