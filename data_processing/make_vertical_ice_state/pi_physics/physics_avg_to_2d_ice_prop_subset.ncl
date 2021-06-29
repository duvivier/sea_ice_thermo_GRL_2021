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
;varcode = "Tinz"
varcode0 = "vicen"
; Choose Sinz,Tinz

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "/glade/campaign/cgd/ppc/dbailey/"+case+"/ice/proc/tseries/month_1/"
mask_dir = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/masks/"
grid_dir = "/glade/p/cesm/omwg/grids/"

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
print("Loading CESM grid data")
; get dimensions
b_file = addfile(grid_dir+"gx1v7_grid.nc","r")
lat2d = b_file->TLAT
lon2d = b_file->TLONG
tarea = b_file->TAREA
kmt = b_file->KMT
z_t = b_file->z_t
; mask tarea where there is land
tarea = where(kmt .gt. 0.5, tarea, tarea@_FillValue)
; convert from cm^2 to m^2
tarea = tarea/100.0/100.0
;dims = dimsizes(lat2d)
;nlat = dims(0)
;nlon = dims(1)
delete(b_file)

; subset over Antarctic coast only
lat2d_sub = lat2d(0:36,:)
lon2d_sub = lon2d(0:36,:)
tarea_sub = tarea(0:36,:)

; get dimensions
dims = dimsizes(lat2d_sub)
nlat = dims(0)
nlon = dims(1)
nlevs = 8
ncat = 5

;;;;;;;;;;
; Get CESM masks for region
;;;;;;;;;;
print("get mask of region of interest")
;; get hemispheric mask
fin = addfile(mask_dir+"antarctic_ocean_masks.nc","r")
sh_mask    = fin->sh_mask
coast_mask = fin->coast_mask
Wed_mask   = fin->Wed_mask
Ross_mask  = fin->Ross_mask
Ind_mask   = fin->Ind_mask
Pac_mask   = fin->Pac_mask
BAm_mask   = fin->BAm_mask
delete(fin)

; now get areas for these regions
coast = coast_mask(0:36,:)
Wed = Wed_mask(0:36,:)      
Ross = Ross_mask(0:36,:)      
Ind = Ind_mask(0:36,:)      
Pac = Pac_mask(0:36,:)      
BAm = BAm_mask(0:36,:)      

;;;;;;;;;;
; Load CESM data
;;;;;;;;;;
nmons = 12
nyrs = 100
ntimes = nmons*nyrs

;;;; load the data ;;;;
; loop through each file to load
b_time = new((/ntimes/),"double")
data_raw = new((/ntimes,ncat,nlat,nlon/),"float")
data_aicen = new((/ntimes,ncat,nlat,nlon/),"float")

files = systemfunc("ls "+data_dir+case+"*."+varcode+".*.nc")
files0 = systemfunc("ls "+data_dir+case+"*."+varcode0+".*.nc")
nfiles = dimsizes(files)
; loop through the files
n = 0
f = 0
do f = 0, nfiles -1
   fin = addfile(files(f),"r")
   fin0 = addfile(files0(f),"r")
   ; deal with time index
   dum_t = fin->time
   dim_t = dimsizes(dum_t)
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
      dum = fin->$varcode$(:,:,:,0:36,:)
      ; average over vertical levels because they're equally spaced
      data_raw(ind_st:ind_ed,:,:,:) = dim_avg_n_Wrap(dum,2)
      delete(dum)
      data_aicen(ind_st:ind_ed,:,:,:) = fin0->$varcode0$(:,:,0:36,:)
      n = n+1
   end if
   delete([/fin,fin0,dum_t,dim_t/])
end do

;;;;;;;;;;
; Convert time arrays
;;;;;;;;;;
b_time_conv = cd_calendar(b_time,0)
yy_b = b_time_conv(:,0)
mm_b = b_time_conv(:,1)

; need to convert months/years
mm_b = mm_b - 1
inds = ind(mm_b .eq. 0)
mm_b(inds) = 12
yy_b(inds) = yy_b(inds) - 1
copy_VarCoords(data_raw(:,0,0,0),yy_b)
yy_b@long_name = "year"
copy_VarCoords(data_raw(:,0,0,0),mm_b)
mm_b@long_name = "month"

;;;;;;;;;;
; Calculate the bulk over the ice categories
;;;;;;;;;;
data_coast = dim_sum_n_Wrap((data_raw*data_aicen),1)/dim_sum_n_Wrap(data_aicen,1)
copy_VarCoords(data_raw(:,0,:,:),data_coast)
copy_VarMeta(data_raw,data_coast)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Saving to netcdf files for plotting")
outdir = "./"

fout = outdir+case+"_coastal_"+varcode
ncdf = addfile(fout+".nc","c")
ncdf->coast_mask = coast
ncdf->Wed_mask = Wed
ncdf->Ross_mask = Ross
ncdf->Ind_mask = Ind
ncdf->Pac_mask = Pac
ncdf->BAm_mask = BAm
ncdf->yy_all = yy_b
ncdf->mm_all = mm_b
ncdf->tarea = tarea_sub
ncdf->lat2d = lat2d_sub
ncdf->lon2d = lon2d_sub
ncdf->data_coast = data_coast
print("Done with file "+fout)
delete([/fout,ncdf/])

print("Done with all files")

;;;;;;;;;;;;;;;;;;;;;; END script
end
