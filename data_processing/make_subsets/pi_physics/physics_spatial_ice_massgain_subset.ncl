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
varcode1 = "frazil"
varcode2 = "congel"
varcode3 = "snoice"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir = "./"
;data_dir = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case+"_spatial_timeseries/"
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
; Load CESM Data
;;;;;;;;;;
nmons = 12
nyrs = 100
ntimes = nmons*nyrs

; list files
fin1 = addfile(data_dir+case+"_coastal_"+varcode1+".nc","r")
fin2 = addfile(data_dir+case+"_coastal_"+varcode2+".nc","r")
fin3 = addfile(data_dir+case+"_coastal_"+varcode3+".nc","r")
; load data
data_1 = fin1->data_coast
data_2 = fin2->data_coast
data_3 = fin3->data_coast
; load times
yy_b = fin1->yy_all
mm_b = fin1->mm_all

; calculate growth
data_frazil_congel = data_1+data_2
copy_VarMeta(data_1,data_frazil_congel)
data_frazil_congel@long_name = "Sea-Ice growth frazil+congelation"

data_all = data_1+data_2+data_3
copy_VarMeta(data_1,data_all)
data_all@long_name = "Sea-Ice growth frazil+congelation+snow"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Output to a netcdf so I can plot together
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Saving to netcdf file for plotting")
outdir = "./"
fout1 = outdir+case+"_coastal_frazil_plus_congel"
ncdf1 = addfile(fout1+".nc","c")
ncdf1->coast_mask = coast
ncdf1->Wed_mask = Wed
ncdf1->Ross_mask = Ross
ncdf1->Ind_mask = Ind
ncdf1->Pac_mask = Pac
ncdf1->BAm_mask = BAm
ncdf1->yy_all = yy_b
ncdf1->mm_all = mm_b
ncdf1->tarea = tarea_sub
ncdf1->lat2d = lat2d_sub
ncdf1->lon2d = lon2d_sub
ncdf1->data_coast = data_frazil_congel

print("file 1 done")

fout2 = outdir+case+"_coastal_frazil_plus_congel_plus_snoice"
ncdf2 = addfile(fout2+".nc","c")
ncdf2->coast_mask = coast
ncdf2->Wed_mask = Wed
ncdf2->Ross_mask = Ross
ncdf2->Ind_mask = Ind
ncdf2->Pac_mask = Pac
ncdf2->BAm_mask = BAm
ncdf2->yy_all = yy_b
ncdf2->mm_all = mm_b
ncdf2->tarea = tarea_sub
ncdf2->lat2d = lat2d_sub
ncdf2->lon2d = lon2d_sub
ncdf2->data_coast = data_all

print("Done with file")

;;;;;;;;;;;;;;;;;;;;;; END script
end
