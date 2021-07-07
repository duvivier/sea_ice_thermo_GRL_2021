;**************************************
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
;varcode = "meltb"

; set flag about if we need to get total over a region
if ((varcode.eq."frazil").or.(varcode.eq."congel").or.(varcode.eq."snoice").or.(varcode.eq."dvidtd").or.(varcode.eq."meltt").or.(varcode.eq."meltb").or.(varcode.eq."meltl").or.(varcode.eq."frazil_plus_congel").or.(varcode.eq."frazil_plus_congel_plus_snoice"))
  flag = "Y"
else
  flag = "N"
end if

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case+"_spatial_timeseries/"

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
nmons = 12
mons = ispan(1,nmons,1)
nyrs = 100
ntimes = nmons*nyrs

; how many regions
nreg = 6

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("get grid info and read in processed data")
fin1 = addfile(data_dir+case+"_coastal_"+varcode+".nc","r")
var1 = fin1->data_coast

; load masks
coast_mask = fin1->coast_mask
Wed_mask = fin1->Wed_mask
Ross_mask = fin1->Ross_mask
Ind_mask = fin1->Ind_mask
Pac_mask = fin1->Pac_mask
BAm_mask = fin1->BAm_mask
yy_all = fin1->yy_all
mm_all = fin1->mm_all
lat2d = fin1->lat2d
lon2d = fin1->lon2d
tarea = fin1->tarea
; convert from cm^2 to m^2
tarea = tarea/100.0/100.0

;;;;;;;;;;
; Deal with masks
;;;;;;;;;;
; choose non-coastal regions
ASL_mask = tofloat(where(ismissing(coast_mask), 1.0, tarea@_FillValue))
; choose the ASL area
ASL_mask = tofloat(where((lon2d.ge.190).and.(lon2d.le.270),ASL_mask,tarea@_FillValue))
; set fill value
ASL_mask@_FillValue = tofloat(tarea@_FillValue)

; now get total areas for these regions
ASL = tofloat(tarea*ASL_mask)
ASL_tot = dim_sum_n(dim_sum_n(ASL,1),0)

;;;;;;;;;;
; Get regional values
;;;;;;;;;;
print("getting regional data")
; now mask and average the data and put into a single array
data_all = new((/ntimes/),"float")

if (flag.eq."Y")
   ; if it's a mass budget term, get total mass growth over region and divide
   ; by total regional area to get "avg"
   data_all(:) = dim_sum_n(dim_sum_n(var1*conform(var1,ASL,(/1,2/)),2),1)/ASL_tot
else
   ; if it's any other variable, just average over the region of interest
   data_all(:) = dim_avg_n(dim_avg_n(var1*conform(var1,ASL_mask,(/1,2/)),2),1)
end if

;;;;;;;;;;
; Prep for output
;;;;;;;;;;
; add coordinates
data_all@units = var1@units
data_all@long_name = var1@long_name
data_all!0 = "time"
data_all&time = var1&time

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Saving to netcdf file for plotting")
outdir = "./"
fout = outdir+case+"_"+varcode+"_ASL_regional_avg"
ncdf = addfile(fout+".nc","c")
ncdf->yy_all = yy_all
ncdf->mm_all = mm_all
ncdf->data_all = data_all

print("Done with file")

;;;;;;;;;;;;;;;;;;;;;; END script
end
