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
; Calculate annual mean
;;;;;;;;;;
print("Calculating annual means")
; case 1
inds = ind((yy_all_1.ge.yr_min).and.(yy_all_1.le.yr_max))
data_1_annual_avg = dim_avg_n_Wrap(data_1(inds,:,:),0)
data_1_annual_var = dim_variance_n_Wrap(data_1(inds,:,:),0)
; case 2
data_2_annual_avg = dim_avg_n_Wrap(data_2(inds,:,:),0)
data_2_annual_var = dim_variance_n_Wrap(data_2(inds,:,:),0)

; get difference
diff = data_2_annual_avg - data_1_annual_avg

; is difference significant
prob = 100.*(1-ttest(data_1_annual_avg,data_1_annual_var,ntimes,data_2_annual_avg,data_2_annual_var,ntimes,True,False))

copy_VarMeta(data_1_annual_avg,diff)
copy_VarMeta(data_1_annual_avg,prob)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting transect in sigma space")
   outdir = "./"

main_cntrs = fspan(-6,6,13)
diff_cntrs = fspan(-0.8,0.8,17)
colors = "temp_diff_18lev"

;;;;;;;;;;;;;
; XY Cross Section Plot Resources
;;;;;;;;;;;;;
; set resources used for all plots
  opts                             = True              ; plot mods desired
  opts@gsnFrame                    = False             ; Do not draw plot
  opts@gsnDraw                     = False             ; Do not advance frame
  opts@gsnMaximize                 = False             ; maximize plot
  opts@gsnLeftString               = ""
  opts@gsnRightString              = ""
  opts@gsnCenterString             = ""
  ;opts@vpWidthF                    = 1  ;1.25               ; set aspect ratio
  ;opts@vpHeightF                   = 1 ;0.75

; set title resources
  opts@tiMainOn                    = False
  opts@tiMainFontHeightF           = 0.02 ;0.03
  opts@tiMainFontThicknessF        = 1.5
  opts@trXMinF                     = -79
  opts@trXMaxF                     = -40 ;-60
  ;opts@trYMinF                     = 32
  ;opts@trYMaxF                     = 38
  opts@trYReverse                  = True
  opts@tmXBOn                      = True
  opts@tiXAxisFontHeightF          = 0.02 ;0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.02 ;0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.02 ;0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.02 ;0.03              ; set tick mark label size for y axis
  opts@tiYAxisString               = "Sigma level (kg m~S~-3~N~)"
  opts@tiXAxisString               = "Latitude"

  ; contour resources for everything
  opts@cnFillOn                    = True          ; turn on color
  opts@cnFillMode                  = "AreaFill"
  opts@cnLinesOn                   = False
  opts@lbLabelBarOn                = False         ; turn off label bars
  opts@cnInfoLabelOn               = False         ; turn off contour info labels
  opts@cnLineLabelsOn              = True
  opts@cnMissingValFillColor       = "grey85"
  opts@cnLevelSelectionMode        = "ExplicitLevels"
  ;opts@cnLevels                    = main_cntrs

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.025 ;0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 1.5                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 1.5                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelLabelBar           = True
  pnlres@lbOrientation              = "horizontal"
  pnlres@lbTitleOn                  = True
  pnlres@lbTitleString              = data_1@units
  pnlres@lbTitlePosition            = "Right"
  pnlres@lbTitleDirection           = "Across"
  pnlres@lbTitleFontHeightF         = 0.01 ;0.013
  pnlres@lbLabelStride              = 2
  pnlres@lbLabelFontHeightF         = 0.01 ;0.013

; set resources for 95% statistical significance overlay
  res_95                            = True               ; reset options
  res_95@gsnLeftString              = ""
  res_95@gsnRightString             = ""
  res_95@gsnFrame                   = False             ; Do not draw plot
  res_95@gsnDraw                    = False             ; Do not advance fr
  res_95@cnInfoLabelOn              = False             ; Turn off info on contours at bottom
  res_95@lbLabelBarOn               = False             ; Turn off label bar
  res_95@cnFillOn                   = True              ; turn on filled contours
  res_95@cnMonoFillColor            = True
  res_95@cnMonoFillPattern          = False
  res_95@cnLevelSelectionMode       = "ManualLevels"    ; set manual contour
  res_95@cnMinLevelValF             = 0.0              ; set min contour level specified above
  res_95@cnMaxLevelValF             = 95.              ; set max contour level specified above
  res_95@cnLevelSpacingF            = 95.0              ; set contour spacing specified above
  res_95@cnFillPatterns             = (/-1,17,-1/)      ; stipling where significant and above threshold

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make actual plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;
; Make annual  mean figure
;;;;;;;;;;;;;
fout = outdir+"pi_physics_MOC_transect_MUSHY"

; open wks
  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,colors)

  nplots = 1
  plot = new(nplots, graphic)
  n = 0
  opts@gsnCenterString = "a) Annual Mean MOC - MUSHY"
  opts@cnLevels                    = main_cntrs
  cplot = gsn_csm_contour(wks,data_2_annual_avg(:,:),opts)
  plot(n) = cplot

; panel plot
print("Paneling")
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,plot/])

;;;;;;;;;;;;;
; Make annual  mean figure
;;;;;;;;;;;;;
fout = outdir+"pi_physics_MOC_transect_BL99"

; open wks
  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,colors)

  nplots = 1
  plot = new(nplots, graphic)
  n = 0
  opts@gsnCenterString = "BL99 MOC"
  opts@cnLevels                    = main_cntrs
  cplot = gsn_csm_contour(wks,data_1_annual_avg(:,:),opts)
  plot(n) = cplot

; panel plot
print("Paneling")
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,plot/])


;;;;;;;;;;;;;
; Make winter mean figure
;;;;;;;;;;;;;
fout = outdir+"pi_physics_MOC_transect_diff"

; open wks
  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,colors)

  nplots = 1
  plot = new(nplots, graphic)
  n = 0
  opts@gsnCenterString = "b) MOC annual mean difference (MUSHY-BL99)"
  delete(opts@cnLevels)
  opts@cnLevels                    = diff_cntrs
  opts@cnLineLabelsOn              = False
  cplot = gsn_csm_contour(wks,diff(:,:),opts)
  sig = gsn_csm_contour(wks,prob(:,:),res_95)
  overlay(cplot,sig)
  plot(n) = cplot

; panel plot
print("Paneling")
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")


;;;;;;;;;;;;;;;;;;;;;; END script
end
