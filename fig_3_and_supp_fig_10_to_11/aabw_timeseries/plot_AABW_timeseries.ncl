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
case_1 = "b.e21.B1850.f09_g17.CMIP6-piControl.001_bl99b"
case_2 = "b.e21.B1850.f09_g17.CMIP6-piControl.001b"

; set years we want to average over
yr_min = 871
yr_max = 970
nyrs = (yr_max-yr_min)+1

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory

data_dir_1 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case_1+"_aabw/"
data_dir_2 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case_2+"_aabw/"
;data_dir_1 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case_1+"_aabw/"
;data_dir_2 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case_2+"_aabw/"

;;;;;;;;;;
; Load CESM data
;;;;;;;;;;
nmons = 12
ntimes = nmons*nyrs
; other dimensions
nsig = 86
nlat_aux_grid = 395

;;;; load the data ;;;;
; case 1
print("Loading "+case_1)
time_1 = new((/ntimes/),"double")
data_1 = new((/ntimes,nsig,nlat_aux_grid/),"float")

files_1 = systemfunc("ls "+data_dir_1+case_1+"*.nc")
nfiles = dimsizes(files_1)
; loop through the files
f = 0
do f = 0, nfiles -1
   fin = addfile(files_1(f),"r")
   ; deal with time index
   dum_t = fin->time
   ; deal with indices
   dim_t = dimsizes(dum_t)
   if (f .eq. 0)
      ind_st = f
      ind_ed = dim_t -1
      ; load some grid info
      lat = fin->lat_aux_grid
      sig2 = fin->sigma
   else
      ind_st = ind_ed + 1
      ind_ed = ind_st+dim_t - 1
   end if
   time_1(ind_st:ind_ed) = dum_t
   ; put data in the arrays
   data_1(ind_st:ind_ed,:,:) = fin->MOC(:,0,:,:)
   delete([/fin,dum_t,dim_t/])
end do
delete(nfiles)

; case 2
print("Loading "+case_2)
time_2 = new((/ntimes/),"double")
data_2 = new((/ntimes,nsig,nlat_aux_grid/),"float")

files_2 = systemfunc("ls "+data_dir_2+case_2+"*.nc")
nfiles = dimsizes(files_1)
; loop through the files
f = 0
do f = 0, nfiles -1
   fin = addfile(files_2(f),"r")
   ; deal with time index
   dum_t = fin->time
   ; deal with indices
   dim_t = dimsizes(dum_t)
   if (f .eq. 0)
      ind_st = f
      ind_ed = dim_t -1
      ; load some grid info
      ;lat = fin->lat_aux_grid
      ;sig2 = fin->sigma
   else
      ind_st = ind_ed + 1
      ind_ed = ind_st+dim_t - 1
   end if
   time_2(ind_st:ind_ed) = dum_t
   ; put data in the arrays
   data_2(ind_st:ind_ed,:,:) = fin->MOC(:,0,:,:)
   delete([/fin,dum_t,dim_t/])
end do
delete(nfiles)

;;;;;;;;;;
; Convert time arrays
;;;;;;;;;;
; time arrays are identical
time_conv = cd_calendar(time_1,0)
yy_all_1 = time_conv(:,0)
mm_all_1 = time_conv(:,1)

; need to convert months/years
mm_all_1 = mm_all_1 - 1
inds = ind(mm_all_1 .eq. 0)
mm_all_1(inds) = 12
yy_all_1(inds) = yy_all_1(inds) - 1
yy_all_1@long_name = "year"
mm_all_1@long_name = "month"
delete(inds)

yy_1 = onedtond(yy_all_1,(/nyrs,nmons/))
mm_1 = onedtond(mm_all_1,(/nyrs,nmons/))

;;;;;;;;;;
; Calculate AABW
;;;;;;;;;;
print("Calculating AABW")
; set cutoff latitude
ilat_60S = min(ind(lat.ge.-60.0))

; min MOC S of 60S
; basically find the minimum value across latitudes,then across signma levels
; then take absolute value to get a timeseries of AABW by month
AABW_1 = abs(dim_min(dim_min(data_1(:,:,:ilat_60S))))
AABW_2 = abs(dim_min(dim_min(data_2(:,:,:ilat_60S))))

;;;;;;;;;;
; Calculate annual mean
;;;;;;;;;;
print("Calculating annual means")
data_annual = new((/2,nyrs/),typeof(AABW_1))

y = 0
do y = 0, nyrs-1
   ; annual means
   inds = ind(yy_all_1 .eq. yy_1(y,0))
   ; case 1
   data_annual(0,y) = dim_avg_n(AABW_1(inds),0)
   ; case 2
   data_annual(1,y) = dim_avg_n(AABW_2(inds),0)
   delete(inds)
end do

;;;;;;;;;;
; Find if mean differences are significant
;;;;;;;;;;
print("Calculating significance")
; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level

; averages over the years of interest
inds = ind((yy_1(:,0).ge.yr_min).and.(yy_1(:,0).le.yr_max))

; annual means
prob = 100.*(1-ttest(avg(data_annual(0,inds)),variance(data_annual(0,inds)),dimsizes(data_annual(0,inds)),avg(data_annual(1,inds)),variance(data_annual(1,inds)),dimsizes(data_annual(1,inds)),True,False))

prob_var = 100.*(1-ftest(variance(data_annual(0,inds)),dimsizes(data_annual(0,inds)),variance(data_annual(1,inds)),dimsizes(data_annual(1,inds)),0))

   print("Probability annual average is different: "+tostring(prob))
   print(case_1+" avg: "+tostring(avg(data_annual(0,inds))))
   print(case_2+" avg: "+tostring(avg(data_annual(1,inds))))
   print("Probability variability is different: "+tostring(prob_var))
   print(case_1+" stddev: "+tostring(stddev(data_annual(0,inds))))
   print(case_2+" stddev: "+tostring(stddev(data_annual(1,inds))))

;;;;;;;;;;
; Get running means and combine to plot
;;;;;;;;;;
nyrsavg = 10
data_plot = new((/4,nyrs/),float)

; annual average
data_plot(0:1,:) = data_annual
data_plot(2:3,:) = runave_n_Wrap(data_annual,nyrsavg,0,1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting timeseries")
   outdir = "./"
;;;;;;;;;;;;;
; XY Cross Section Plot Resources
;;;;;;;;;;;;;
; set resources used for all plots
  opts                             = True              ; plot mods desired
  opts@gsnFrame                    = False             ; Do not draw plot
  opts@gsnDraw                     = False             ; Do not advance frame
  opts@gsnMaximize                 = False             ; maximize plot
  opts@vpWidthF                    = 1.5               ; set aspect ratio
  opts@vpHeightF                   = 0.5

; set title resources
  opts@tiMainOn                    = False
  opts@tiMainFontHeightF           = 0.02 ;0.03
  opts@tiMainFontThicknessF        = 1.5
  ;opts@trXMinF                     = 0
  ;opts@trXMaxF                     = 50
  opts@gsnYRefLine                 = 0.0
  opts@tmXBOn                      = True
  opts@tiXAxisFontHeightF          = 0.02 ;0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.02 ;0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.02 ;0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.02 ;0.03              ; set tick mark label size for y axis
  opts@tiYAxisString               = data_1@units
  opts@tiXAxisString               = "Years"

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025
  opts@xyMarkLineMode              = "Lines"
  opts@xyLineColors                = (/"grey","lightgoldenrod","black","goldenrod"/)
  opts@xyLineThicknesses           = (/2,2,4.0,4.0/)
  opts@xyDashPattern               = 0

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 1.5                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 1.5                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True             ; fill the page

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make actual plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 fout = outdir+"/pi_physics_coastal_AABW"
 pref = "e) Annual mean AABW production"
 wks = gsn_open_wks("ps" ,fout)
 plot = new(1, graphic)
 plot(0) = gsn_csm_xy(wks,yy_1(:,0),data_plot,opts)

  ; make panel plot
  pnlres@txString  = pref
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")
  delete([/fout,plot/])

print("Done with seasonal plots")

;;;;;;;;;;;;;;;;;;;;;; END script
end
