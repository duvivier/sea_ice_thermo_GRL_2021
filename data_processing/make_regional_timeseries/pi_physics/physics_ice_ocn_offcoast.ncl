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
;varcode = "fswthru"

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

; make mask of off-coastal regions
all_coast = tofloat(where((coast_mask .eq. 1.0),tarea@_FillValue,1.0))
all_coast@_FillValue = tofloat(tarea@_FillValue)
Wed_coast = tofloat(where((coast_mask .eq. 1.0) .and. (Wed_mask .eq. 1.0),tarea@_FillValue,1.0))
Wed_coast@_FillValue = tofloat(tarea@_FillValue)
Ross_coast = tofloat(where((coast_mask .eq. 1.0) .and. (Ross_mask .eq. 1.0),tarea@_FillValue,1.0))
Ross_coast@_FillValue = tofloat(tarea@_FillValue)
Ind_coast = tofloat(where((coast_mask .eq. 1.0) .and. (Ind_mask .eq. 1.0),tarea@_FillValue,1.0))
Ind_coast@_FillValue = tofloat(tarea@_FillValue)
Pac_coast = tofloat(where((coast_mask .eq. 1.0) .and. (Pac_mask .eq. 1.0),tarea@_FillValue,1.0))
Pac_coast@_FillValue = tofloat(tarea@_FillValue)
BAm_coast = tofloat(where((coast_mask .eq. 1.0) .and. (BAm_mask .eq. 1.0),tarea@_FillValue,1.0))
BAm_coast@_FillValue = tofloat(tarea@_FillValue)

; now get total areas for these regions
coast = tofloat(tarea*all_coast)
coast_tot = dim_sum_n(dim_sum_n(coast,1),0)
Wed = tofloat(tarea*Wed_coast)
Wed_tot = dim_sum_n(dim_sum_n(Wed,1),0)
Ross = tofloat(tarea*Ross_coast)
Ross_tot = dim_sum_n(dim_sum_n(Ross,1),0)
Ind = tofloat(tarea*Ind_coast)
Ind_tot = dim_sum_n(dim_sum_n(Ind,1),0)
Pac = tofloat(tarea*Pac_coast)
Pac_tot = dim_sum_n(dim_sum_n(Pac,1),0)
BAm = tofloat(tarea*BAm_coast)
BAm_tot = dim_sum_n(dim_sum_n(BAm,1),0)

;;;;;;;;;;
; Get regional values
;;;;;;;;;;
print("getting regional data")
; now mask and average the data and put into a single array
data_all = new((/6,ntimes/),"float")

if (flag.eq."Y")
   ; if it's a mass budget term, get total mass growth over region and divide
   ; by total regional area to get "avg"
   ; all coastal pts
   data_all(0,:) = dim_sum_n(dim_sum_n(var1*conform(var1,coast,(/1,2/)),2),1)/coast_tot
   ; Ross sea
   data_all(1,:) = dim_sum_n(dim_sum_n(var1*conform(var1,Ross,(/1,2/)),2),1)/Ross_tot
   ; BAm
   data_all(2,:) = dim_sum_n(dim_sum_n(var1*conform(var1,BAm,(/1,2/)),2),1)/BAm_tot
   ; Weddell
   data_all(3,:) = dim_sum_n(dim_sum_n(var1*conform(var1,Wed,(/1,2/)),2),1)/Wed_tot
   ; Indian
   data_all(4,:) = dim_sum_n(dim_sum_n(var1*conform(var1,Ind,(/1,2/)),2),1)/Ind_tot
   ; Pacific
   data_all(5,:) = dim_sum_n(dim_sum_n(var1*conform(var1,Pac,(/1,2/)),2),1)/Pac_tot
else
   ; if it's any other variable, just average over the region of interest
   ; all coastal pts
   data_all(0,:) = dim_avg_n(dim_avg_n(var1*conform(var1,all_coast,(/1,2/)),2),1)
   ; Ross sea
   data_all(1,:) = dim_avg_n(dim_avg_n(var1*conform(var1,Ross_coast,(/1,2/)),2),1)
   ; BAm
   data_all(2,:) = dim_avg_n(dim_avg_n(var1*conform(var1,BAm_coast,(/1,2/)),2),1)
   ; Weddell
   data_all(3,:) = dim_avg_n(dim_avg_n(var1*conform(var1,Wed_coast,(/1,2/)),2),1)
   ; Indian
   data_all(4,:) = dim_avg_n(dim_avg_n(var1*conform(var1,Ind_coast,(/1,2/)),2),1)
   ; Pacific
   data_all(5,:) = dim_avg_n(dim_avg_n(var1*conform(var1,Pac_coast,(/1,2/)),2),1)
end if

;;;;;;;;;;
; Prep for output
;;;;;;;;;;
; add coordinates
data_all@units = var1@units
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
fout = outdir+case+"_offcoastal_"+varcode+"_regional_avgs"
ncdf = addfile(fout+".nc","c")
ncdf->yy_all = yy_all
ncdf->mm_all = mm_all
ncdf->data_all = data_all

print("Done with file")

;;;;;;;;;;;;;;;;;;;;;; END script
end
