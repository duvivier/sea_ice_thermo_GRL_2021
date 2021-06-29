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
;varcode = "U10"
; Choose LHFLX,PBLH,PSL,QFLX,QREFHT,SHFLX,SST,TREFHT,TS,U10

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "/glade/campaign/cgd/ppc/dbailey/"+case+"/atm/proc/tseries/month_1/"

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get dimensions
b_file = addfile(data_dir+case+".cam.h0."+varcode+".087101-092012.nc","r")
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
; Load CESM data
;;;;;;;;;;
nmons = 12
nyrs = 100
ntimes = nmons*nyrs

;;;; load the historical data ;;;;
; loop through each file to load
b_time = new((/ntimes/),"double")
data_coast = new((/ntimes,nlat,nlon/),"float")

files = systemfunc("ls "+data_dir+case+"*."+varcode+".*.nc")
nfiles = dimsizes(files)
; loop through the files
n = 0
f = 0
do f = 0, nfiles -1
   fin = addfile(files(f),"r")
   ; deal with time index
   dum_t = fin->time
   dim_t = dimsizes(dum_t)
   ; don't load year 1 files since they're not correct
   if (dum_t(0) .gt. 800)
      ; deal with indices
      if (n .eq. 0)
         ind_st = n
         ind_ed = dim_t -1
      else
         ind_st = ind_ed + 1
         ind_ed = ind_st+dim_t - 1
      end if
      b_time(ind_st:ind_ed) = dum_t
      ; put data in the arrays
      data_coast(ind_st:ind_ed,:,:) = fin->$varcode$(:,0:40,:)
      n = n+1
   end if
   delete([/fin,dum_t,dim_t/])
end do

;;;;;;;;;;
; Convert time arrays
;;;;;;;;;;
b_time_conv = cd_calendar(b_time,0)
; deal with both time periods
yy_b = b_time_conv(:,0)
mm_b = b_time_conv(:,1)

; need to convert months/years
mm_b = mm_b - 1
inds = ind(mm_b .eq. 0)
mm_b(inds) = 12
yy_b(inds) = yy_b(inds) - 1 
copy_VarCoords(data_coast(:,0,0),yy_b)
yy_b@long_name = "year"
copy_VarCoords(data_coast(:,0,0),mm_b)
mm_b@long_name = "month"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Saving to netcdf file for plotting")
outdir = "./"
fout = outdir+case+"_coastal_"+varcode
ncdf = addfile(fout+".nc","c")
ncdf->yy_all = yy_b
ncdf->mm_all = mm_b
ncdf->lat2d = lat2d_sub
ncdf->lon2d = lon2d_sub
ncdf->data_coast = data_coast

print("Done with file")

;;;;;;;;;;;;;;;;;;;;;; END script
end
