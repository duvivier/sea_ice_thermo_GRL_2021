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
case1 = "b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b"
case2 = "b.e21.B1850.f09_g17.CMIP6-piControl.001b"
;varcode1 = "CLDLOW"
;varcode1 = "PBLH"
varcode1 = "TURBFLX"

; set years we want to average over
yr_min = 871
yr_max = 970
nyrs = (yr_max-yr_min)+1

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir1 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case1+"_spatial_timeseries/"
data_dir2 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case2+"_spatial_timeseries/"
grid_dir = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/masks/"

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("get grid info")
fin0 = addfile(grid_dir+case2+".cam.h0.LANDFRAC.087101-092012.nc","r")
; get land only over Antarctic, 0 is ocean point, 1 is all land, fractions depend.
landfrac = fin0->LANDFRAC(0,0:40,:)

print("reading in processed data")
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode1+".nc","r")
data_1 = fin1->data_coast
lat2d = fin1->lat2d
lon2d = fin1->lon2d
yy_1 = fin1->yy_all
mm_1 = fin1->mm_all
dims = dimsizes(lat2d)
nlat = dims(0)
nlon = dims(1)
delete([/fin1,dims/])

fin2 = addfile(data_dir2+case2+"_coastal_"+varcode1+".nc","r")
data_2 = fin2->data_coast
yy_2 = fin2->yy_all
mm_2 = fin2->mm_all

; set units
units = data_1@units

;;;;;;;;;;
; Get winter means (AMJJAS)
;;;;;;;;;;
inds = ind(((mm_1.ge.4).and.(mm_1.le.9)).and.((yy_1.ge.yr_min).and.(yy_1.le.yr_max)))
data_1_avg = dim_avg_n(data_1(inds,:,:),0)
data_1_var = dim_variance_n(data_1(inds,:,:),0)
delete(inds)
inds = ind(((mm_2.ge.4).and.(mm_2.le.9)).and.((yy_2.ge.yr_min).and.(yy_2.le.yr_max)))
data_2_avg = dim_avg_n(data_2(inds,:,:),0)
data_2_var = dim_variance_n(data_2(inds,:,:),0)
delete(inds)

;;;;;;;;;;
; Finding differences
;;;;;;;;;;
; calculate differences
diff = data_2_avg - data_1_avg

; calculate percent difference
pcnt_diff = ((data_2_avg - data_1_avg)/data_1_avg)*100.0

; calculate significance
prob = 100.*(1-ttest(data_1_avg,data_1_var,(nyrs*6),data_2_avg,data_2_var,(nyrs*6),True,False))

; mask where the differences are not significant
diff = where((prob.ge.95.0),diff,diff@_FillValue)
pcnt_diff = where((prob.ge.95.0),pcnt_diff,pcnt_diff@_FillValue)

;;;;;;;;;;
; Deal with masks
;;;;;;;;;;
; want masks over ocean/ice and south of 60S
coast_mask = lat2d
coast_mask = where((landfrac.le.0.5),1.0,lat2d@_FillValue)

; make mask of coastal regions
all_coast = tofloat(coast_mask)
Ross_coast = all_coast
Ross_coast = where(((lon2d.ge.160.0).and.(lon2d.lt.220.0).and.(lat2d.le.-66.0).and.(landfrac.le.0.5)),1.0,all_coast@_FillValue)
BAm_coast = all_coast
BAm_coast = where(((lon2d.ge.220.0).and.(lon2d.lt.300.0).and.(lat2d.le.-66.0).and.(landfrac.le.0.5)),1.0,all_coast@_FillValue)
Wed_coast = all_coast
Wed_coast = where((((lon2d.ge.300.0).or.(lon2d.lt.20.0)).and.(lat2d.le.-66.0).and.(landfrac.le.0.5)),1.0,all_coast@_FillValue)
Ind_coast = all_coast
Ind_coast = where(((lon2d.ge.20.0).and.(lon2d.lt.90.0).and.(lat2d.le.-60.0).and.(landfrac.le.0.5)),1.0,all_coast@_FillValue)
Pac_coast = all_coast
Pac_coast = where(((lon2d.ge.90.0).and.(lon2d.lt.160.0).and.(lat2d.le.-60.0).and.(landfrac.le.0.5)),1.0,all_coast@_FillValue)

;;;;;;;;;;
; Calculate totals by regions
;;;;;;;;;;
if (varcode1 .eq. "SHFLX")
;if (varcode1 .eq. "TURBFLX")
   ; if looking at fluxes, calculate the total energy difference over winter
   ; diff (W/m2)*tarea (m2)* seconds in winter (s) gives units of Joules
   ; THEN sum over significant points by region
   ; total days in winter Apr-Sept = 183
   nsec = 183*24*60*60



return

   print("Calculating regional total significant differences for "+varcode1+" "+units)
   print("All coast: "+tostring(sum(diff*all_coast)))
   print("Indian Ocean coast: "+tostring(sum(diff*Ind_coast)))
   print("Pacific Ocean coast: "+tostring(sum(diff*Pac_coast)))
   print("Ross Sea coast: "+tostring(sum(diff*Ross_coast)))
   print("BAM Sea coast: "+tostring(sum(diff*BAm_coast)))
   print("Weddell Sea coast: "+tostring(sum(diff*Wed_coast)))
else
   ; PBLH, needs average, not sum
   print("Calculating regional average significant differences for "+varcode1+" "+units)
   print("All coast: "+tostring(avg(diff*all_coast)))
   print("Indian Ocean coast: "+tostring(avg(diff*Ind_coast)))
   print("Pacific Ocean coast: "+tostring(avg(diff*Pac_coast)))
   print("Ross Sea coast: "+tostring(avg(diff*Ross_coast)))
   print("BAM Sea coast: "+tostring(avg(diff*BAm_coast)))
   print("Weddell Sea coast: "+tostring(avg(diff*Wed_coast)))
   print("All coast: "+tostring(avg(pcnt_diff*all_coast)))
   print("Indian Ocean coast: "+tostring(avg(pcnt_diff*Ind_coast)))
   print("Pacific Ocean coast: "+tostring(avg(pcnt_diff*Pac_coast)))
   print("Ross Sea coast: "+tostring(avg(pcnt_diff*Ross_coast)))
   print("BAM Sea coast: "+tostring(avg(pcnt_diff*BAm_coast)))
   print("Weddell Sea coast: "+tostring(avg(pcnt_diff*Wed_coast)))


end if

;;;;;;;;;;;;;;;;;;;;;; END script
end
