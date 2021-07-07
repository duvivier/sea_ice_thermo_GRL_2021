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
;varcode1 = "aice"
varcode0 = "aice"

; set years we want to average over
yr_min = 871
yr_max = 970
nyrs = (yr_max-yr_min)+1

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory

data_dir1 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case1+"_spatial_timeseries/"
data_dir2 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case2+"_spatial_timeseries/"
;data_dir1 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case1+"_spatial_timeseries/"
;data_dir2 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case2+"_spatial_timeseries/"

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("reading in processed data")
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode0+".nc","r")
data_0a = fin1->data_coast
delete(fin1)
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode1+".nc","r")
data_1 = fin1->data_coast
coast_mask = fin1->coast_mask
lat2d = fin1->lat2d
lon2d = fin1->lon2d
yy_1 = fin1->yy_all
mm_1 = fin1->mm_all
dims = dimsizes(lat2d)
nlat = dims(0)
nlon = dims(1)
delete([/fin1,dims/])

fin2 = addfile(data_dir2+case2+"_coastal_"+varcode0+".nc","r")
data_0b = fin2->data_coast
delete(fin2)
fin2 = addfile(data_dir2+case2+"_coastal_"+varcode1+".nc","r")
data_2 = fin2->data_coast
yy_2 = fin2->yy_all
mm_2 = fin2->mm_all
delete(fin2)

if (varcode1 .eq. "aice")
   data_1 = data_1*100.0
   data_2 = data_2*100.0
   data_1@units = "%"
end if

; set units
units = data_1@units

;;;;;;;;;;
; Get winter means (AMJJAS)
;;;;;;;;;;
inds = ind(((mm_1.ge.4).and.(mm_1.le.9)).and.((yy_1.ge.yr_min).and.(yy_1.le.yr_max)))
data_1_avg = dim_avg_n(data_1(inds,:,:),0)
data_1_var = dim_variance_n(data_1(inds,:,:),0)
data_0a_avg = dim_avg_n(data_0a(inds,:,:),0)
delete(inds)
inds = ind(((mm_2.ge.4).and.(mm_2.le.9)).and.((yy_2.ge.yr_min).and.(yy_2.le.yr_max)))
data_2_avg = dim_avg_n(data_2(inds,:,:),0)
data_2_var = dim_variance_n(data_2(inds,:,:),0)
data_0b_avg = dim_avg_n(data_0b(inds,:,:),0)
delete(inds)

; mask where ice concentration is less than 15%
data_1_avg = where(data_0a_avg.ge.0.15,data_1_avg,data_1_avg@_FillValue)
data_2_avg = where(data_0b_avg.ge.0.15,data_2_avg,data_2_avg@_FillValue)

; mask where not sector of interest
data_1_avg = where((lon2d.ge.190).and.(lon2d.le.270),data_1_avg,data_1_avg@_FillValue)
data_2_avg = where((lon2d.ge.190).and.(lon2d.le.270),data_2_avg, data_2_avg@_FillValue)

; mask all points that are coastal
data_1_avg = where(ismissing(coast_mask), data_1_avg,data_1_avg@_FillValue)
data_2_avg = where(ismissing(coast_mask), data_2_avg,data_2_avg@_FillValue)

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
if (varcode1.eq."aice")
   pref = "Difference in sea ice concentration ("+units+")"
   main_cntrs = fspan(-5,5,11)
end if
if (varcode1.eq."hi")
   pref = "Difference in sea ice thickness ("+units+")"
   main_cntrs = fspan(-0.2,0.2,11)
end if
if (varcode1.eq."frazil_plus_congel_plus_snoice")
   pref = "Total thermodynamic ice growth ("+units+")"
   main_cntrs = fspan(-2,2,13)
end if
if (varcode1.eq."snoice")
   pref = "Snow ice growth ("+units+")"
   main_cntrs = fspan(-0.4,0.4,9)
end if
if (varcode1.eq."frazil")
   pref = "Frazil ice growth ("+units+")"
   main_cntrs = fspan(-2,2,13)
end if
if (varcode1.eq."congel")
   pref = "Congelation ice growth ("+units+")"
   main_cntrs = fspan(-1,1,9)
end if
if (varcode1.eq."U10")
   pref = "10m wind speed ("+units+")"
   main_cntrs = fspan(-0.3,0.3,9)
end if
if (varcode1.eq."PSL")
   pref = "Sea level pressure ("+units+")"
   main_cntrs = fspan(-4.0,4.0,17)
end if
if (varcode1.eq."PBLH")
   pref = "Planetary boundary layer height ("+units+")"
   main_cntrs = fspan(-30.0,30.0,13)
end if
if (varcode1.eq."PRECT")
   pref = "Total precipitation ("+units+")"
   main_cntrs = fspan(-0.03,0.03,11)
end if
if (varcode1.eq."TURBFLX")
   pref = "Total turbulent heat flux ("+units+")"
   main_cntrs = fspan(-8,8,9)
end if
if (varcode1.eq."LHFLX")
   pref = "Total latent heat flux ("+units+")"
   main_cntrs = fspan(-3,3,13)
end if
if (varcode1.eq."SHFLX")
   pref = "Total sensible heat flux ("+units+")"
   main_cntrs = fspan(-5,5,11)
end if
if (varcode1.eq."TREFHT")
   pref = "2m temperature ("+units+")"
   main_cntrs = fspan(-5,5,11)
end if
if (varcode1.eq."QREFHT")
   pref = "2m humidity ("+units+")"
   main_cntrs = fspan(-5,5,11)
end if
if (varcode1.eq."TEMP_100m")
   pref = "Temperature at 100m depth ("+units+")"
   main_cntrs = fspan(-0.5,0.5,11)
end if
if (varcode1.eq."SALT_100m")
   pref = "Salinity at 100m depth ("+units+")"
   main_cntrs = fspan(-0.1,0.1,11)
end if
if (varcode1.eq."PD_100m")
   pref = "Potential density at 100m depth ("+units+")"
   main_cntrs = fspan(-0.1,0.1,11)
end if
if (varcode1.eq."HMXL")
   pref = "Ocean mixed layer depth ("+units+")"
   main_cntrs = fspan(-30.0, 30.0,15)
end if
if (varcode1.eq."fhocn")
   pref = "Heat flux ice to ocean ("+units+")"
   main_cntrs = fspan(-10.0,10.0,11)
end if
if (varcode1.eq."fresh")
   pref = "Freshwater flux ice to ocean ("+units+")"
   main_cntrs = fspan(-1.0,1.0,11)
end if
if (varcode1.eq."fswthru")
   pref = "Shortwave flux through ice ("+units+")"
   main_cntrs = fspan(-2,2,15)
end if
if (varcode1.eq."totChl_150m")
   pref = "Total chlorophyll integrated over top 150m ("+units+")"
   main_cntrs = fspan(-0.02,0.02,15)
end if
if (varcode1.eq."totChl_100m")
   pref = "Total chlorophyll integrated over top 100m ("+units+")"
   main_cntrs = fspan(-0.02,0.02,15)
end if
if (varcode1.eq."totChl_50m")
   pref = "Total chlorophyll integrated over top 50m ("+units+")"
   main_cntrs = fspan(-0.02,0.02,15)
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
  fout = outdir+"pi_physics_"+varcode1+"_winter_diff_ASL"

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
