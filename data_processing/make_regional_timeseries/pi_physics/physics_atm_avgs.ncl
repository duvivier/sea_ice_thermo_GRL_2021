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
;varcode = "U10"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case+"_spatial_timeseries/"
grid_dir = "/glade/campaign/cgd/ppc/dbailey/"+case+"/atm/proc/tseries/month_1/"

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
print("get grid info")
fin0 = addfile(grid_dir+case+".cam.h0.LANDFRAC.087101-092012.nc","r")
; get land only over Antarctic, 0 is ocean point, 1 is all land, fractions depend.
landfrac = fin0->LANDFRAC(0,0:40,:)

; load the main file
fin1 = addfile(data_dir+case+"_coastal_"+varcode+".nc","r")
var1 = fin1->data_coast

; load grid data
yy_all = fin1->yy_all
mm_all = fin1->mm_all
lat2d = fin1->lat2d
lon2d = fin1->lon2d

;;;;;;;;;;
; Deal with masks
;;;;;;;;;;
coast_mask = lat2d
coast_mask = where(((landfrac.le.0.5).and.(lat2d.le.-60.0)),1.0,lat2d@_FillValue)

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
; Get regional values
;;;;;;;;;;
print("getting regional data")
; now mask and average the data and put into a single array
data_all = new((/6,ntimes/),"float")

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
fout = outdir+case+"_coastal_"+varcode+"_regional_avgs"
ncdf = addfile(fout+".nc","c")
ncdf->yy_all = yy_all
ncdf->mm_all = mm_all
ncdf->data_all = data_all

print("Done with file")

;;;;;;;;;;;;;;;;;;;;;; END script
end
