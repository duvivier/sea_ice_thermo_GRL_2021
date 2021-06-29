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
;varcode1 = "TEMP"
;TEMP,SALT,PD

; set years we want to average over
yr_min = 871
yr_max = 970
nyrs = (yr_max-yr_min)+1

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
;data_dir1 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case1+"_spatial_timeseries/"
;data_dir2 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case2+"_spatial_timeseries/"
data_dir1 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case1+"_spatial_timeseries/"
data_dir2 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case2+"_spatial_timeseries/"

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("reading in processed data")
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode1+".nc","r")
data_1 = fin1->data_coast
z_t = fin1->z_t
; change vertical grid to m
z_t = z_t/100.0
z_t@units = "m"
ndepth = dimsizes(z_t)
depth_orig = -1.0*z_t
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

if (varcode1.eq."PD")
   data_1 = data_1*1000.0
   data_2 = data_2*1000.0
   data_1@units = "kg/m3"
end if
if (varcode1.eq."SALT")
   scale = 1000.0
   data_1 = data_1*scale
   data_2 = data_2*scale
   data_1@units = "g/kg"
end if

; set units
units = data_1@units

;;;;;;;;;;
; Caculate zonal means
;;;;;;;;;;
print("calculating zonal means")
data_1_zonal = zonalAve(data_1)
data_2_zonal = zonalAve(data_2)
; assign lat
lat1d = lat2d(:,0)
data_1_zonal&nlat = lat1d
data_2_zonal&nlat = lat1d

;;;;;;;;;;;;
; Vertically interpolate
;;;;;;;;;;;;
print("interpolating to regular depth")
; create vertical spacing in meters
  depth_max = -2000
  n_space = 5
  depth_grid = ispan(-1,depth_max,n_space)
  depth_grid!0 = "Vertical"
  depth_grid@units = "m"
  n_depth = dimsizes(depth_grid)

data_1_interp = linint1_n_Wrap(depth_orig,data_1_zonal,False,depth_grid,0,1)
data_2_interp = linint1_n_Wrap(depth_orig,data_2_zonal,False,depth_grid,0,1)

;;;;;;;;;;
; Get winter means (AMJJAS)
;;;;;;;;;;
print("calculating winter mean")
inds = ind(((mm_1.ge.4).and.(mm_1.le.9)).and.((yy_1.ge.yr_min).and.(yy_1.le.yr_max)))
data_1_avg = dim_avg_n_Wrap(data_1_interp(inds,:,:),0)
data_1_var = dim_variance_n_Wrap(data_1_interp(inds,:,:),0)
delete(inds)
inds = ind(((mm_2.ge.4).and.(mm_2.le.9)).and.((yy_2.ge.yr_min).and.(yy_2.le.yr_max)))
data_2_avg = dim_avg_n_Wrap(data_2_interp(inds,:,:),0)
data_2_var = dim_variance_n_Wrap(data_2_interp(inds,:,:),0)
delete(inds)

;;;;;;;;;;
; Finding differences
;;;;;;;;;;
; calculate differences
diff = data_2_avg - data_1_avg
copy_VarMeta(data_1_avg,diff)

; calculate significance
prob = 100.*(1-ttest(data_1_avg,data_1_var,(nyrs*6),data_2_avg,data_2_var,(nyrs*6),True,False))
copy_VarMeta(data_1_avg,prob)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting differences "+varcode1)
   outdir = "./"

; set contour levels based on varcode
if (varcode1.eq."TEMP")
   pref = "f) Temperature" 
   ;main_cntrs = fspan(-1.0,1.0,11)
   main_cntrs = fspan(-0.5,0.5,11)
end if
if (varcode1.eq."SALT")
   pref = "d) Salinity"
   ;main_cntrs = fspan(-0.1,0.1,11)
   main_cntrs = fspan(-0.07,0.07,15)
end if
if (varcode1.eq."PD")
   pref = "b) Potential density" 
   ;main_cntrs = fspan(-0.1,0.1,11)
   main_cntrs = fspan(-0.07,0.07,15)
end if

colors = "temp_diff_18lev"

; set general resources
  res                             = True          ; plot mods desired
  res@gsnDraw                     = False         ; don't draw the plot yet
  res@gsnFrame                    = False         ; or the frame
  res@cnLinesOn                   = False         ; turn off contour lines
  res@txFontHeightF               = 0.02          ; set subtitle text font
  res@gsnLeftString               = ""
  res@gsnRightString              = ""
  res@gsnCenterString             = ""
  ;res@vpWidthF                    = 1 ;1.25
  ;res@vpHeightF                   = 1 ;0.75
  res@trYMaxF                     = -1
  res@trYMinF                     = -1000
  res@tiYAxisString               = "Depth (m)"
  res@tiYAxisFontHeightF          = 0.02
  res@tmYLLabelFontHeightF        = 0.02
  ;res@trXMaxF                     = -60
  ;res@trXMinF                     = -80
  res@tiXAxisString               = "Latitude"
  res@tiXAxisFontHeightF          = 0.02
  res@tmXBLabelFontHeightF        = 0.02
  res@tmYLLabelFontHeightF        = 0.02

  ; contour resources for everything
  res@cnFillOn                    = True          ; turn on color
  res@cnFillMode                  = "AreaFill"
  res@cnLinesOn                   = False
  res@lbLabelBarOn                = False         ; turn off label bars
  res@cnInfoLabelOn               = False         ; turn off contour info labels
  res@cnLineLabelsOn              = True ;False
  res@cnMissingValFillColor       = "grey85"
  res@cnLevelSelectionMode        = "ExplicitLevels"
  res@cnLevels                    = main_cntrs

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 2                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelLabelBar           = True
  pnlres@lbOrientation              = "horizontal"
  pnlres@lbTitleOn                  = True
  pnlres@lbTitleString              = units
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
; Make winter mean figure
;;;;;;;;;;;;;
fout = outdir+"pi_physics_"+varcode1+"_zonal_mean_transect_winter_diff"

; open wks
  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,colors)

  nplots = 1
  plot = new(nplots, graphic)
  
  n = 0   
  res@gsnCenterString = pref+" zonal mean"
  cplot = gsn_csm_contour(wks,diff(:,:),res)
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

delete([/fout,plot/])

print("Done with plots")

;;;;;;;;;;;;;;;;;;;;;; END script
end
