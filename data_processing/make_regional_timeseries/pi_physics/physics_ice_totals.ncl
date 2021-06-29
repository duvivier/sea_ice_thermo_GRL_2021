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
;case = "b.e21.B1850.f09_g17.CMIP6-piControl.001_snow1"
;varcode = "hi"

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
; reset missing value info
coast_mask@_FillValue = -9999
coast_mask = where(ismissing(coast_mask), 0.0, coast_mask)
Wed_mask@_FillValue = -9999
Wed_mask = where(ismissing(Wed_mask), 0.0, Wed_mask)
Ross_mask@_FillValue = -9999
Ross_mask = where(ismissing(Ross_mask), 0.0, Ross_mask)
Ind_mask@_FillValue = -9999
Ind_mask = where(ismissing(Ind_mask), 0.0, Ind_mask)
Pac_mask@_FillValue = -9999
Pac_mask = where(ismissing(Pac_mask), 0.0, Pac_mask)
BAm_mask@_FillValue = -9999
BAm_mask = where(ismissing(BAm_mask), 0.0, BAm_mask)

; make mask of coastal regions
coast = tofloat(where((coast_mask .eq. 1.0),tarea,tarea@_FillValue))
Wed_coast = tofloat(where((coast_mask .eq. 1.0) .and. (Wed_mask .eq. 1.0),tarea,tarea@_FillValue))
Ross_coast = tofloat(where((coast_mask .eq. 1.0) .and. (Ross_mask .eq. 1.0),tarea,tarea@_FillValue))
Ind_coast = tofloat(where((coast_mask .eq. 1.0) .and. (Ind_mask .eq. 1.0),tarea,tarea@_FillValue))
Pac_coast = tofloat(where((coast_mask .eq. 1.0) .and. (Pac_mask .eq. 1.0),tarea,tarea@_FillValue))
BAm_coast = tofloat(where((coast_mask .eq. 1.0) .and. (BAm_mask .eq. 1.0),tarea,tarea@_FillValue))

;;;;;;;;;;
; Get regional values
;;;;;;;;;;
print("getting regional data")
; now mask and average the data and put into a single array
data_all = new((/6,ntimes/),"float")
; all coastal pts
data_all(0,:) = dim_sum_n(dim_sum_n(var1*conform(var1,coast,(/1,2/)),2),1)
; Ross sea
data_all(1,:) = dim_sum_n(dim_sum_n(var1*conform(var1,Ross_coast,(/1,2/)),2),1)
; BAm
data_all(2,:) = dim_sum_n(dim_sum_n(var1*conform(var1,BAm_coast,(/1,2/)),2),1)
; Weddell
data_all(3,:) = dim_sum_n(dim_sum_n(var1*conform(var1,Wed_coast,(/1,2/)),2),1)
; Indian
data_all(4,:) = dim_sum_n(dim_sum_n(var1*conform(var1,Ind_coast,(/1,2/)),2),1)
; Pacific
data_all(5,:) = dim_sum_n(dim_sum_n(var1*conform(var1,Pac_coast,(/1,2/)),2),1)

;;;;;;;;;;
; Prep for output
;;;;;;;;;;
; fix units if necessary
if (flag .eq. "Y")
   ; convert to millions of cubic meters per day growth
   data_all = data_all/(10.0^6)
   units = "10~S~6~N~ m~S~3~N~/day"
end if
if (flag .eq. "N")
   ;wind speeds
   units = var1@units
   if (varcode .eq. "aice")
      ; convert to km^2 ice area
      data_all = data_all/(10.0^6)
      units = "km~S~2~N~"
   end if
   if (varcode .eq. "hi")
      ; convert to km^3 ice area
      data_all = data_all/(10.0^9)
      units = "km~S~3~N~"
   end if
end if

; add coordinates
data_all@units = units
data_all@long_name = var1@long_name
data_all!0 = "region"
;regions = (/"all","Ross","BAm","Wed","Ind","Pac"/)
;data_all&region = regions
data_all!1 = "time"
data_all&time = var1&time

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Saving to netcdf file for plotting")
outdir = "./"
fout = outdir+case+"_coastal_"+varcode+"_regional_totals"
ncdf = addfile(fout+".nc","c")
ncdf->yy_all = yy_all
ncdf->mm_all = mm_all
ncdf->data_all = data_all

print("Done with file")

;;;;;;;;;;;;;;;;;;;;;; END script
end
