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

varcode1 = "frazil_plus_congel_plus_snoice"
;varcode1 = "congel"
;varcode2 = "frazil"
;varcode3 = "snoice"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir1 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case1+"_spatial_timeseries/"
data_dir2 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case2+"_spatial_timeseries/"

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
nmons = 12
ndperm = (/31,28,31,30,31,30,31,31,30,31,30,31/)
nyrs = 100
ntimes = nmons*nyrs 

; how many regions
reg_str = (/"all","Ross","BAm","Wed","Ind","Pac"/)
nreg = dimsizes(reg_str)

;;;;;;;;;;
; Deal with grids
;;;;;;;;;;
print("reading in grid info")
; code is that first # is the variable name and second # is case name

fin1 = addfile(data_dir1+case1+"_coastal_"+varcode1+".nc","r")
; load masks
coast_mask = fin1->coast_mask
lat2d = fin1->lat2d
lon2d = fin1->lon2d
tarea = fin1->tarea
delete(fin1)
;; convert from cm^2 to m^2
;tarea = tarea/100.0/100.0

; reset masks
;coast_mask@_FillValue = -9999
;coast_mask = where(ismissing(coast_mask), 0.0, coast_mask)
;coast = tofloat(where((coast_mask .eq. 1.0),tarea,tarea@_FillValue))

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("reading in processed data")
; code is that first # is the variable name and second # is case name

; case 1
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode1+".nc","r")
; read in data
var1 = fin1->data_coast
; convert time arrays
yy_all_1 = fin1->yy_all
mm_all_1 = fin1->mm_all
yy_1 = onedtond(yy_all_1,(/nyrs,nmons/))
mm_1 = onedtond(mm_all_1,(/nyrs,nmons/))

; case 2
fin2 = addfile(data_dir2+case2+"_coastal_"+varcode1+".nc","r")
; read in data
var2 = fin2->data_coast
; convert time arrays
yy_all_2 = fin2->yy_all
mm_all_2 = fin2->mm_all
yy_2 = onedtond(yy_all_2,(/nyrs,nmons/))
mm_2 = onedtond(mm_all_2,(/nyrs,nmons/))

;;;;;;;;;;
; Get total volume by month
;;;;;;;;;;
print("Calculating volume per cell")
; cm/day*cm2/cell*days/month = cm3/month per cell

data_1_vol = new(dimsizes(var1),typeof(var1))
data_2_vol = new(dimsizes(var2),typeof(var2))
data_1_vol_coast = new(dimsizes(var1),typeof(var1))
data_2_vol_coast = new(dimsizes(var2),typeof(var2))

m = 0
do m = 0, nmons -1
   inds = ind(mm_all_1 .eq. (m+1))
   tmp1 = var1(inds,:,:)
   data_1_vol(inds,:,:) = ndperm(m)*tmp1*conform(tmp1,tofloat(tarea),(/1,2/))
   data_1_vol_coast(inds,:,:) = ndperm(m)*tmp1*conform(tmp1,tofloat(tarea),(/1,2/))*conform(tmp1,coast_mask,(/1,2/))
   tmp2 = var2(inds,:,:)
   data_2_vol(inds,:,:) = ndperm(m)*tmp2*conform(tmp2,tofloat(tarea),(/1,2/))
   data_2_vol_coast(inds,:,:) = ndperm(m)*tmp2*conform(tmp2,tofloat(tarea),(/1,2/))*conform(tmp2,coast_mask,(/1,2/))
   delete([/inds,tmp1,tmp2/])
end do

;; convert to m3 from cm3
;data_1_vol = data_1_vol/(100.*100.*100.)
;data_2_vol = data_2_vol/(100.*100.*100.)
;data_1_vol_coast = data_1_vol_coast/(100.*100.*100.)
;data_2_vol_coast = data_2_vol_coast/(100.*100.*100.)

;;;;;;;;;;
; Sum over all points
;;;;;;;;;;
var1_all_prod = dim_sum_n(data_1_vol,(/1,2/))
var2_all_prod = dim_sum_n(data_2_vol,(/1,2/))
var1_all_prod_coast = dim_sum_n(data_1_vol_coast,(/1,2/))
var2_all_prod_coast = dim_sum_n(data_2_vol_coast,(/1,2/))

;;;;;;;;;;
; Get fraction of total ice production
;;;;;;;;;;
var_1_frac = 100.*(var1_all_prod_coast/var1_all_prod)
var_2_frac = 100.*(var2_all_prod_coast/var2_all_prod)

;;;;;;;;;;
; Get seasonal average fractions
;;;;;;;;;;
seas_str = (/"AMJJAS"/)
inds = ind((mm_all_1.ge.4).and.(mm_all_1.le.9))

var_1_avg = avg(var_1_frac(inds))
var_2_avg = avg(var_2_frac(inds))

print("Fraction of total ice production in coastal region")
print(case1+": "+tostring(var_1_avg))
print(case2+": "+tostring(var_2_avg))

;;;;;;;;;;;;;;;;;;;;;; END script
end
