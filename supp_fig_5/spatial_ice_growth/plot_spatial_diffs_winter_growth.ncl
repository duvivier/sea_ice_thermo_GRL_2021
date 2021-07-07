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
varcode1 = "congel"
varcode2 = "frazil"
varcode3 = "snoice"
varcode4 = "aice"

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

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("reading in processed data")
; order of data label is data_N1N2
; N1 = variable
; N2 = case

; case 1
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode1+".nc","r")
data_11 = fin1->data_coast
coast_mask = fin1->coast_mask
lat2d = fin1->lat2d
lon2d = fin1->lon2d
yy_1 = fin1->yy_all
mm_1 = fin1->mm_all
dims = dimsizes(lat2d)
nlat = dims(0)
nlon = dims(1)
delete([/fin1,dims/])
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode2+".nc","r")
data_21 = fin1->data_coast
delete(fin1)
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode3+".nc","r")
data_31 = fin1->data_coast
delete(fin1)
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode4+".nc","r")
data_41 = fin1->data_coast
delete(fin1)

; case 2
fin2 = addfile(data_dir2+case2+"_coastal_"+varcode1+".nc","r")
data_12 = fin2->data_coast
yy_2 = fin2->yy_all
mm_2 = fin2->mm_all
delete(fin2)
fin2 = addfile(data_dir2+case2+"_coastal_"+varcode2+".nc","r")
data_22 = fin2->data_coast
delete(fin2)
fin2 = addfile(data_dir2+case2+"_coastal_"+varcode3+".nc","r")
data_32 = fin2->data_coast
delete(fin2)
fin2 = addfile(data_dir2+case2+"_coastal_"+varcode4+".nc","r")
data_42 = fin2->data_coast
delete(fin2)

; set units
units = data_11@units

;;;;;;;;;;
; Get winter means (AMJJAS)
;;;;;;;;;;
inds = ind(((mm_1.ge.4).and.(mm_1.le.9)).and.((yy_1.ge.yr_min).and.(yy_1.le.yr_max)))
data_11_avg = dim_avg_n(data_11(inds,:,:),0)
data_11_var = dim_variance_n(data_11(inds,:,:),0)
data_21_avg = dim_avg_n(data_21(inds,:,:),0)
data_21_var = dim_variance_n(data_21(inds,:,:),0)
data_31_avg = dim_avg_n(data_31(inds,:,:),0)
data_31_var = dim_variance_n(data_31(inds,:,:),0)
data_41_avg = dim_avg_n(data_41(inds,:,:),0)
data_41_var = dim_variance_n(data_41(inds,:,:),0)
delete(inds)
inds = ind(((mm_2.ge.4).and.(mm_2.le.9)).and.((yy_2.ge.yr_min).and.(yy_2.le.yr_max)))
data_12_avg = dim_avg_n(data_12(inds,:,:),0)
data_12_var = dim_variance_n(data_12(inds,:,:),0)
data_22_avg = dim_avg_n(data_22(inds,:,:),0)
data_22_var = dim_variance_n(data_22(inds,:,:),0)
data_32_avg = dim_avg_n(data_32(inds,:,:),0)
data_32_var = dim_variance_n(data_32(inds,:,:),0)
data_42_avg = dim_avg_n(data_42(inds,:,:),0)
data_42_var = dim_variance_n(data_42(inds,:,:),0)
delete(inds)

print("masking where aice lt 15%")
data_11_avg = where(data_41_avg.le.0.15,data_11_avg@_FillValue,data_11_avg)
data_21_avg = where(data_41_avg.le.0.15,data_21_avg@_FillValue,data_21_avg)
data_31_avg = where(data_41_avg.le.0.15,data_31_avg@_FillValue,data_31_avg)
data_12_avg = where(data_42_avg.le.0.15,data_12_avg@_FillValue,data_12_avg)
data_22_avg = where(data_42_avg.le.0.15,data_22_avg@_FillValue,data_22_avg)
data_32_avg = where(data_42_avg.le.0.15,data_32_avg@_FillValue,data_32_avg)

;;;;;;;;;;
; Finding differences
;;;;;;;;;;
; calculate differences
diff1 = data_12_avg - data_11_avg
diff2 = data_22_avg - data_21_avg
diff3 = data_32_avg - data_31_avg

; calculate significance
prob1 = 100.*(1-ttest(data_11_avg,data_11_var,(nyrs*6),data_12_avg,data_12_var,(nyrs*6),True,False))
prob2 = 100.*(1-ttest(data_21_avg,data_21_var,(nyrs*6),data_22_avg,data_22_var,(nyrs*6),True,False))
prob3 = 100.*(1-ttest(data_31_avg,data_31_var,(nyrs*6),data_32_avg,data_32_var,(nyrs*6),True,False))

; assign lat and lon
diff1@lat2d = lat2d
diff1@lon2d = lon2d
prob1@lat2d = lat2d
prob1@lon2d = lon2d
diff2@lat2d = lat2d
diff2@lon2d = lon2d
prob2@lat2d = lat2d
prob2@lon2d = lon2d
diff3@lat2d = lat2d
diff3@lon2d = lon2d
prob3@lat2d = lat2d
prob3@lon2d = lon2d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting differences in growth terms")
   outdir = "./"

; set contour levels based on varcode
main_cntrs = fspan(-1.6,1.6,17)
;main_cntrs = fspan(-1.5,1.5,31)
colors = "BlRe"

; set general resources
  res                             = True          ; plot mods desired
  res@gsnDraw                     = False         ; don't draw the plot yet
  res@gsnFrame                    = False         ; or the frame
  res@cnLinesOn                   = False         ; turn off contour lines
  res@txFontHeightF               = 0.030          ; set subtitle text font
  res@gsnLeftString               = ""
  res@gsnRightString              = ""
  res@gsnCenterString             = ""
  res@gsnAddCyclic                = True          ; global data, wrap around

; Map resources for main dataset
  res@gsnPolar                    = "SH"
  res@mpMaxLatF                   = -62.0  ; -60.
  res@mpGridAndLimbOn             = False            ; Turn on lat/lon grid
  res@mpFillOn                    = True             ; Enables map fill
  res@mpLandFillColor             = "grey"
  res@mpGeophysicalLineColor      = "black"          ; Color of continents
  res@mpGeophysicalLineThicknessF = 1.5
  res@gsnPolarLabelSpacing        = 45
  res@gsnPolarLabelFontHeightF    = 0.02

  ; contour resources for everything
  res@cnFillOn                    = True          ; turn on color
  res@cnFillMode                  = "AreaFill"
  res@cnLinesOn                   = False
  res@lbLabelBarOn                = False         ; turn off label bars
  res@cnInfoLabelOn               = False         ; turn off contour info labels
  res@cnLineLabelsOn              = False
  res@cnMissingValFillColor       = "grey"
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
  pnlres@lbLabelStride              = 2
  pnlres@lbTitleFontHeightF         = 0.013
 
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

; line resources
  lnres                           = True
  lnres@gsLineThicknessF          = 1.0
  lnres@gsLineColor               = "black"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make actual plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;
; Make winter mean figure
;;;;;;;;;;;;;
  fout = outdir+"pi_physics_ice_growth_winter_diff"

; open wks
  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,colors)

  nplots = 3
  plot = new(nplots, graphic)
  dum0 = new(nplots, graphic)
  dum1 = new(nplots, graphic)
  dum2 = new(nplots, graphic)
  dum3 = new(nplots, graphic)
  dum4 = new(nplots, graphic)

  n = 0   
  res@gsnCenterString = "a) Congelation ice growth"
  cplot = gsn_csm_contour_map(wks,diff1,res)
  sig = gsn_csm_contour(wks,prob1,res_95)
  overlay(cplot,sig)
  plot(n) = cplot
  dum0(n) = gsn_add_polyline(wks,plot(n),(/20.0,20.0/),(/-60.0,-90.0/),lnres)
  dum1(n) = gsn_add_polyline(wks,plot(n),(/90.0,90.0/),(/-60.0,-90.0/),lnres)
  dum2(n) = gsn_add_polyline(wks,plot(n),(/160.0,160.0/),(/-60.0,-90.0/),lnres)
  dum3(n) = gsn_add_polyline(wks,plot(n),(/220.0,220.0/),(/-60.0,-90.0/),lnres)
  dum4(n) = gsn_add_polyline(wks,plot(n),(/300.0,300.0/),(/-60.0,-90.0/),lnres)
  n = 1   
  res@gsnCenterString = "b) Frazil ice growth"
  cplot = gsn_csm_contour_map(wks,diff2,res)
  sig = gsn_csm_contour(wks,prob2,res_95)
  overlay(cplot,sig)
  plot(n) = cplot
  dum0(n) = gsn_add_polyline(wks,plot(n),(/20.0,20.0/),(/-60.0,-90.0/),lnres)
  dum1(n) = gsn_add_polyline(wks,plot(n),(/90.0,90.0/),(/-60.0,-90.0/),lnres)
  dum2(n) = gsn_add_polyline(wks,plot(n),(/160.0,160.0/),(/-60.0,-90.0/),lnres)
  dum3(n) = gsn_add_polyline(wks,plot(n),(/220.0,220.0/),(/-60.0,-90.0/),lnres)
  dum4(n) = gsn_add_polyline(wks,plot(n),(/300.0,300.0/),(/-60.0,-90.0/),lnres)
  n = 2   
  res@gsnCenterString = "c) Snow-ice growth"
  cplot = gsn_csm_contour_map(wks,diff3,res)
  sig = gsn_csm_contour(wks,prob3,res_95)
  overlay(cplot,sig)
  plot(n) = cplot
  dum0(n) = gsn_add_polyline(wks,plot(n),(/20.0,20.0/),(/-60.0,-90.0/),lnres)
  dum1(n) = gsn_add_polyline(wks,plot(n),(/90.0,90.0/),(/-60.0,-90.0/),lnres)
  dum2(n) = gsn_add_polyline(wks,plot(n),(/160.0,160.0/),(/-60.0,-90.0/),lnres)
  dum3(n) = gsn_add_polyline(wks,plot(n),(/220.0,220.0/),(/-60.0,-90.0/),lnres)
  dum4(n) = gsn_add_polyline(wks,plot(n),(/300.0,300.0/),(/-60.0,-90.0/),lnres)

; panel plot
print("Paneling")
  ; make panel plot
  gsn_panel(wks,(/plot/),(/1,3/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,plot/])

print("Done with plots")

;;;;;;;;;;;;;;;;;;;;;; END script
end
