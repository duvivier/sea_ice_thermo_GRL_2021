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
varcode1 = "Tinz"
; Tinz and Sinz

if ((varcode1.eq."aice").or.(varcode1.eq."hi").or.(varcode1.eq."congel").or.(varcode1.eq."frazil").or.(varcode1.eq."snoice").or.(varcode1.eq."frazil_plus_congel_plus_snoice").or.(varcode1.eq."fresh").or.(varcode1.eq."fhocn").or.(varcode1.eq."fswthru").or.(varcode1.eq."totChl"))
   flag = "Y"
else
   flag = "N"
end if

nyrs = 100

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
data_dir1 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/analysis/pi_physics/version_5/vertical_ice_state/"
data_dir2 = data_dir1

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("reading in processed data")
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode1+".nc","r")
data_1 = fin1->data_coast
if (flag.eq."Y")
   coast_mask = fin1->coast_mask
end if
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

; convert units to Kelvin
;data_1 = data_1+273.15
;data_2 = data_2+273.15
;data_1@units = "K"

; set units
units = data_1@units

;;;;;;;;;;
; Get winter means (AMJJAS)
;;;;;;;;;;
inds = ind((mm_1.ge.4).and.(mm_1.le.9))
data_1_avg = dim_avg_n(data_1(inds,:,:),0)
data_1_var = dim_variance_n(data_1(inds,:,:),0)
delete(inds)
inds = ind((mm_2.ge.4).and.(mm_2.le.9))
data_2_avg = dim_avg_n(data_2(inds,:,:),0)
data_2_var = dim_variance_n(data_2(inds,:,:),0)
delete(inds)

if (flag.eq."Y")
   print("masking regionally")
   data_1_avg = data_1_avg*coast_mask
   data_2_avg = data_2_avg*coast_mask
end if

;;;;;;;;;;
; Finding differences
;;;;;;;;;;
; calculate differences
diff = data_2_avg - data_1_avg

; calculate significance
prob = 100.*(1-ttest(data_1_avg,data_1_var,(nyrs*6),data_2_avg,data_2_var,(nyrs*6),True,False))

; assign lat and lon
data_1_avg@lat2d = lat2d
data_1_avg@lon2d = lon2d
data_2_avg@lat2d = lat2d
data_2_avg@lon2d = lon2d
diff@lat2d = lat2d
diff@lon2d = lon2d
prob@lat2d = lat2d
prob@lon2d = lon2d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting differences "+varcode1)
   outdir = "./"

; set contour levels based on varcode
if (varcode1.eq."Tinz")
   pref = "Sea ice internal temp ("+units+")"
   main_cntrs = fspan(-4,4,15)
end if
if (varcode1.eq."Sinz")
   pref = "Sea ice internal salinity ("+units+")"
   main_cntrs = fspan(-8,8,15)
end if

colors = "temp_diff_18lev"

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
  res@mpMaxLatF                   = -60.0
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
  fout = outdir+"pi_physics_"+varcode1+"_winter_diff"

; open wks
  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,colors)

  nplots = 1
  plot = new(nplots, graphic)
  dum0 = new(nplots, graphic)
  dum1 = new(nplots, graphic)
  dum2 = new(nplots, graphic)
  dum3 = new(nplots, graphic)
  dum4 = new(nplots, graphic)

  n = 0   
  res@gsnCenterString = pref
  cplot = gsn_csm_contour_map(wks,diff,res)
  sig = gsn_csm_contour(wks,prob,res_95)
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
  gsn_panel(wks,(/plot/),(/1,1/),pnlres)                   ; create final panel plot
  delete(wks)
print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

delete([/fout,plot/])

print("Done with plots")

;;;;;;;;;;;;;;;;;;;;;; END script
end
