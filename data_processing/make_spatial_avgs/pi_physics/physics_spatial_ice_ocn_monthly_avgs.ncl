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
;varcode = "hi" 

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case+"_spatial_timeseries/"
grid_dir = "/glade/p/cesm/omwg/grids/"

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
nmons = 12
mons = ispan(1,nmons,1)
nyrs = 51
ntimes = nmons*nyrs 

; how many regions
nreg = 6

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("reading in processed data")
fin1 = addfile(data_dir+case+"_coastal_"+varcode+".nc","r")

; load masks
coast_mask = fin1->coast_mask
yy_all = fin1->yy_all
mm_all = fin1->mm_all
lat2d = fin1->lat2d
lon2d = fin1->lon2d
tarea = fin1->tarea
; convert from cm^2 to m^2
tarea = tarea/100.0/100.0
dims = dimsizes(lat2d)
nlat = dims(0)
nlon = dims(1)

; convert time arrays
yy = onedtond(yy_all,(/nyrs,nmons/))
mm = onedtond(mm_all,(/nyrs,nmons/))

;;;;;;;;;;
; Get data
;;;;;;;;;;
; load data
var1 = fin1->data_coast(:,:,:)

; if it's a mass budget term, need to cm/day from kg/m2/s
if ((varcode.eq."sidmassgrowthwat").or.(varcode.eq."sidmassgrowthbot").or.(varcode.eq."sidmasssi").or.(varcode.eq."sidmassdyn").or.(varcode.eq."sidmassmelttop").or.(varcode.eq."sidmassmeltbot").or.(varcode.eq."sidmasslat"))
   print("converting units")
   rhoi = 917.0   ;kg/m3
   timeconv = 60.0*60.0*24.0
   var1 = var1*(1/rhoi)*100.0*timeconv
   var1@units = "cm/day"
end if

;;;;;;;;;;
; Get monthly averages
;;;;;;;;;;
print("Calculating monthly means")

; loop through decades
data_1_avg = new((/nmons,nlat,nlon/),typeof(var1))
data_1_std = new((/nmons,nlat,nlon/),typeof(var1))
data_1_var = new((/nmons,nlat,nlon/),typeof(var1))

   m = 0
   do m = 0, nmons -1
      inds = ind(mm_all.eq.(m+1))
      data_1_avg(m,:,:) = dim_avg_n(var1(inds,:,:),0)
      data_1_std(m,:,:) = dim_stddev_n(var1(inds,:,:),0)
      data_1_var(m,:,:) = dim_variance_n(var1(inds,:,:),0)
      delete(inds)
   end do

; mask all but coastal areas
data_1_avg = data_1_avg*conform(data_1_avg,coast_mask,(/1,2/))
data_1_std = data_1_std*conform(data_1_std,coast_mask,(/1,2/))
data_1_var = data_1_var*conform(data_1_var,coast_mask,(/1,2/))

; assign coordinates
data_1_avg@units = var1@units 
data_1_avg@long_name = var1@long_name
data_1_avg!0 = "months"
data_1_avg!1 = lat2d!0
data_1_avg!2 = lat2d!1
copy_VarMeta(data_1_avg, data_1_std)
copy_VarMeta(data_1_avg, data_1_var)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Saving to netcdf file for plotting")
outdir = "./"
fout = outdir+case+"_coastal_"+varcode+"_monthly_spatial_avgs"
ncdf = addfile(fout+".nc","c")
ncdf->yy_all = yy_all
ncdf->mm_all = mm_all
ncdf->lat2d = lat2d
ncdf->lon2d = lon2d
ncdf->data_avg = data_1_avg
ncdf->data_std = data_1_std
ncdf->data_var = data_1_var

print("Done with file")

;;;;;;;;;;;;;;;;;;;;;; END script
end
