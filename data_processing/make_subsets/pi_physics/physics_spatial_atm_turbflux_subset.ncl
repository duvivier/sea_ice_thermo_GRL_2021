;***************************************
; PROGRAM DESCRIPTION: Create timeseries of CESM ice volume
; INPUT DATA: CESM timeseries
; OUTPUT DATA: plot
; CREATOR: Alice DuVivier - Jan 2019
;***************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;***************************************

begin
;;;;;;;;;;
; Manual inputs
;;;;;;;;;;
;case = "b.e21.B1850.f09_g17.CMIP6-piControl.001b"
;case = "b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b"
varcode1 = "LHFLX"
varcode2 = "SHFLX"
varcode3 = "TURBFLX"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "./"
;data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case+"_spatial_timeseries/"
mask_dir = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/masks/"
grid_dir = "/glade/campaign/cgd/ppc/dbailey/"+case+"/atm/proc/tseries/month_1/"


;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get dimensions
b_file = addfile(grid_dir+case+".cam.h0."+varcode1+".087101-092012.nc","r")
lat1d = b_file->lat
nlat = dimsizes(lat1d)
lon1d = b_file->lon
nlon = dimsizes(lon1d)
delete(b_file)

; calculate lat2d/lon2d
lat2d = new((/nlat,nlon/),double)
lon2d = new((/nlat,nlon/),double)
n = 0
do n = 0,nlon-1
  lat2d(:,n) = lat1d
end do
n = 0
do n = 0,nlat-1
  lon2d(n,:) = lon1d
end do
lat2d!1 = lon2d!1
lon2d!0 = lat2d!0
delete([/nlat,nlon/])

; subset over Antarctic coast only
lat2d_sub = lat2d(0:40,:)
lon2d_sub = lon2d(0:40,:)

; get new dimensions
dims = dimsizes(lat2d_sub)
nlat = dims(0)
nlon = dims(1)

;;;;;;;;;;
; Load CESM Data
;;;;;;;;;;
nmons = 12
nyrs = 100
ntimes = nmons*nyrs
; list files
fin1 = addfile(data_dir+case+"_coastal_"+varcode1+".nc","r")
fin2 = addfile(data_dir+case+"_coastal_"+varcode2+".nc","r")
; load data
data_1 = fin1->data_coast
data_2 = fin2->data_coast
; load times
yy_b = fin1->yy_all
mm_b = fin1->mm_all

; calculate total turbluent flux
data_3 = data_1+data_2
copy_VarMeta(data_1,data_3)
data_3@long_name = "total turbulent heat flux"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Saving to netcdf file for plotting")
outdir = "./"
fout3 = outdir+case+"_coastal_"+varcode3
ncdf3 = addfile(fout3+".nc","c")
ncdf3->yy_all = yy_b
ncdf3->mm_all = mm_b
ncdf3->lat2d = lat2d_sub
ncdf3->lon2d = lon2d_sub
ncdf3->data_coast = data_3

print("Done with file")

;;;;;;;;;;;;;;;;;;;;;; END script
end
