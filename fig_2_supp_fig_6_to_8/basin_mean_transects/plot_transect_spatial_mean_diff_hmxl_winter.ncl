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
varcode0 = "HMXL"
;varcode1 = "SALT"
;TEMP,SALT,PD

; set years we want to average over
yr_min = 871
yr_max = 970
nyrs = (yr_max-yr_min)+1

ntime = 12*nyrs

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
coast_mask = fin1->coast_mask
Wed_mask = fin1->Wed_mask
Ross_mask = fin1->Ross_mask
Ind_mask = fin1->Ind_mask
Pac_mask = fin1->Pac_mask
BAm_mask = fin1->BAm_mask
yy_1 = fin1->yy_all
mm_1 = fin1->mm_all
dims = dimsizes(lat2d)
nlat = dims(0)
nlon = dims(1)
delete([/fin1,dims/])
; load HMXL
fin1 = addfile(data_dir1+case1+"_coastal_"+varcode0+".nc","r")
data_1a = fin1->data_coast
delete(fin1)

fin2 = addfile(data_dir2+case2+"_coastal_"+varcode1+".nc","r")
data_2 = fin2->data_coast
yy_2 = fin2->yy_all
mm_2 = fin2->mm_all
delete(fin2)
fin2 = addfile(data_dir2+case2+"_coastal_"+varcode0+".nc","r")
data_2a = fin2->data_coast
delete(fin2)

if (varcode1.eq."PD")
   data_1 = data_1*1000.0
   data_2 = data_2*1000.0
   ; now convert to sigma
   data_1 = data_1-1000.0
   data_2 = data_2-1000.0
   data_1@units = "kg/m3"
end if
if (varcode1.eq."SALT")
   scale = 1000.0
   data_1 = data_1*scale
   data_2 = data_2*scale
   data_1@units = "g/kg"
end if

; convert HMXL to m and flip to negative
data_1a = data_1a/-100.0
data_2a = data_2a/-100.0
data_1a@units = "m"
data_2a@units = "m"

; set units
units = data_1@units
units_hmxl = data_1a@units

;;;;;;;;;;
; Create transects
;;;;;;;;;;
print("Create transects")
; we want to mask by basin and then average over latitudes
reg_str = (/"Ross","Amundsen-Bellingshausen","Weddell","Indian","Pacific"/)
nreg = dimsizes(reg_str)

; get mixed layer depths by region
print("Calculating HMXL for transect")
data_1a_trans = new((/nreg,ntime,nlat/),float)
data_2a_trans = new((/nreg,ntime,nlat/),float)
; Ross sea
n = 0
tmp = data_1a*conform(data_1a,Ross_mask,(/1,2/))
data_1a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp)
tmp = data_2a*conform(data_2a,Ross_mask,(/1,2/))
data_2a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp)
; AmB
n = 1
tmp = data_1a*conform(data_1a,BAm_mask,(/1,2/))
data_1a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp)
tmp = data_2a*conform(data_2a,BAm_mask,(/1,2/))
data_2a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp)
; Weddell sea
n = 2
tmp = data_1a*conform(data_1a,Wed_mask,(/1,2/))
data_1a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp)
tmp = data_2a*conform(data_2a,Wed_mask,(/1,2/))
data_2a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp); Indian Ocean
n = 3
tmp = data_1a*conform(data_1a,Ind_mask,(/1,2/))
data_1a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp)
tmp = data_2a*conform(data_2a,Ind_mask,(/1,2/))
data_2a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp)
; Pacific Ocean
n = 4
tmp = data_1a*conform(data_1a,Pac_mask,(/1,2/))
data_1a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp)
tmp = data_2a*conform(data_2a,Pac_mask,(/1,2/))
data_2a_trans(n,:,:) = dim_avg_n(tmp,2)
delete(tmp)

; assign meta data
data_1a_trans!0 = "regions"
data_1a_trans!1 = "time"
data_1a_trans&time = data_1&time
data_1a_trans!2 = "nlat"
data_1a_trans&nlat = lat2d(:,0)
copy_VarMeta(data_1a_trans,data_2a_trans)


; get depth transects now
print("Calculating "+varcode1+" for transect")
data_1_trans = new((/nreg,ntime,ndepth,nlat/),float)
data_2_trans = new((/nreg,ntime,ndepth,nlat/),float)
; Ross sea
n = 0
tmp = data_1*conform(data_1,Ross_mask,(/2,3/))
data_1_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)
tmp = data_2*conform(data_2,Ross_mask,(/2,3/))
data_2_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)
; AmB
n = 1
tmp = data_1*conform(data_1,BAm_mask,(/2,3/))
data_1_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)
tmp = data_2*conform(data_2,BAm_mask,(/2,3/))
data_2_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)
; Weddell sea
n = 2
tmp = data_1*conform(data_1,Wed_mask,(/2,3/))
data_1_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)
tmp = data_2*conform(data_2,Wed_mask,(/2,3/))
data_2_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)
; Indian Ocean
n = 3
tmp = data_1*conform(data_1,Ind_mask,(/2,3/))
data_1_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)
tmp = data_2*conform(data_2,Ind_mask,(/2,3/))
data_2_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)
; Pacific Ocean
n = 4
tmp = data_1*conform(data_1,Pac_mask,(/2,3/))
data_1_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)
tmp = data_2*conform(data_2,Pac_mask,(/2,3/))
data_2_trans(n,:,:,:) = dim_avg_n(tmp,3)
delete(tmp)

; assign meta data
data_1_trans!0 = "regions"
data_1_trans!1 = "time"
data_1_trans&time = data_1&time
data_1_trans!2 = "z_t"
data_1_trans&z_t = data_1&z_t
data_1_trans!3 = "nlat"
data_1_trans&nlat = lat2d(:,0)
copy_VarMeta(data_1_trans,data_2_trans)

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

data_1_interp = linint1_n_Wrap(depth_orig,data_1_trans,False,depth_grid,0,2)
data_2_interp = linint1_n_Wrap(depth_orig,data_2_trans,False,depth_grid,0,2)

;;;;;;;;;;
; Get winter means (AMJJAS)
;;;;;;;;;;
print("Getting winter averages")
inds = ind(((mm_1.ge.4).and.(mm_1.le.9)).and.((yy_1.ge.yr_min).and.(yy_1.le.yr_max)))
data_1_avg = dim_avg_n_Wrap(data_1_interp(:,inds,:,:),1)
data_1_var = dim_variance_n_Wrap(data_1_interp(:,inds,:,:),1)
data_1a_avg = dim_avg_n_Wrap(data_1a_trans(:,inds,:),1)
delete(inds)
inds = ind(((mm_2.ge.4).and.(mm_2.le.9)).and.((yy_2.ge.yr_min).and.(yy_2.le.yr_max)))
data_2_avg = dim_avg_n_Wrap(data_2_interp(:,inds,:,:),1)
data_2_var = dim_variance_n_Wrap(data_2_interp(:,inds,:,:),1)
data_2a_avg = dim_avg_n_Wrap(data_2a_trans(:,inds,:),1)
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
   pref1 = "c) Temperature"
   pref2 = "f) Temperature"
   colors = "cmocean_thermal"
   main_cntrs = fspan(-2.0,2.5,19)
   diff_cntrs = fspan(-0.8,0.8,17)
end if
if (varcode1.eq."SALT")
   pref1 = "b) Salinity"
   pref2 = "e) Salinity"
   colors = "cmocean_haline"
   main_cntrs = fspan(33.2,34.8,17)
   diff_cntrs = fspan(-0.14,0.14,15)
end if
if (varcode1.eq."PD")
   pref1 = "a) Sigma"
   pref2 = "d) Sigma"
   colors = "cmocean_dense"
   main_cntrs = fspan(26.8,28.2,15)
   ;main_cntrs = fspan(1026.8,1028.2,15)
   diff_cntrs = fspan(-0.14,0.14,15)
end if
colors_diff = "temp_diff_18lev"

; set general resources
  res                             = True          ; plot mods desired
  res@gsnDraw                     = False         ; don't draw the plot yet
  res@gsnFrame                    = False         ; or the frame
  res@cnLinesOn                   = False         ; turn off contour lines
  res@txFontHeightF               = 0.013          ; set subtitle text font
  res@gsnLeftString               = ""
  res@gsnRightString              = ""
  res@gsnCenterString             = ""
  ;res@vpWidthF                    = 1 ;1.25
  ;res@vpHeightF                   = 1 ;0.75
  res@trYMaxF                     = -1
  res@trYMinF                     = -600 ;-1000
  res@tiYAxisString               = "Depth (m)"
  res@tiYAxisFontHeightF          = 0.013
  res@tmYLLabelFontHeightF        = 0.013
  ;res@trXMaxF                     = -60
  ;res@trXMinF                     = -80
  res@tiXAxisString               = "Latitude"
  res@tiXAxisFontHeightF          = 0.013
  res@tmXBLabelFontHeightF        = 0.013

  ; contour resources for everything
  res@cnFillOn                    = True          ; turn on color
  res@cnFillMode                  = "AreaFill"
  res@cnLinesOn                   = False
  res@lbLabelBarOn                = False         ; turn off label bars
  res@cnInfoLabelOn               = False         ; turn off contour info labels
  ;if (varcode1 .ne. "PD")
     res@cnLineLabelsOn              = True
  ;else
  ;   res@cnLineLabelsOn              = False
  ;end if
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
  pnlres@lbOrientation              = "horizontal" ;"vertical"
  pnlres@lbTitleOn                  = True
  pnlres@lbTitleString              = units
  pnlres@lbTitlePosition            = "Right" ;"Top"
  pnlres@lbTitleDirection           = "Across"
  pnlres@lbTitleFontHeightF         = 0.01
  pnlres@lbLabelStride              = 2
  pnlres@lbLabelFontHeightF         = 0.01

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

; set polyres to show mld from case1
  polyres1                         = True
  polyres1@gsMarkerIndex           = 4               ; polymarker style
  polyres1@gsMarkerSizeF           = 3               ; polymarker size
  polyres1@gsMarkerColor           = "black"

; set polyres to show mld from case 2
  polyres2                         = True
  polyres2@gsMarkerIndex           = 2               ; polymarker style
  polyres2@gsMarkerSizeF           = 3                ; polymarker size
  polyres2@gsMarkerColor           = "black"          ; polymarker color

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make actual plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;
; Make winter mean figures
;;;;;;;;;;;;;
r = 0
do r = 0, nreg -1
  fout = outdir+"pi_physics_"+varcode1+"_"+reg_str(r)+"_regional_MUSHY_mean_transect"
  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,colors)

  nplots = 1
  plot = new(nplots, graphic)
  n = 0   
  res@gsnCenterString = pref1+" "+reg_str(r)+" sector"
  ;res@gsnCenterString = pref1+" "+reg_str(r)+" sector winter mean - MUSHY"
  cplot = gsn_csm_contour(wks,data_2_avg(r,:,:),res)
  dum1 = gsn_add_polymarker(wks,cplot,lat2d(:,0),data_2a_avg(r,:),polyres2)
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
end do
print("Finished MUSHY mean plots")

;;;;;;;;;;;;;
; Make winter mean figures
;;;;;;;;;;;;;
r = 0
do r = 0, nreg -1
  fout = outdir+"pi_physics_"+varcode1+"_"+reg_str(r)+"_regional_BL99_mean_transect"
  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,colors)

  nplots = 1
  plot = new(nplots, graphic)
  n = 0
  res@gsnCenterString = pref1+" "+reg_str(r)+" sector"
  ;res@gsnCenterString = pref1+" "+reg_str(r)+" sector winter mean - BL99"
  cplot = gsn_csm_contour(wks,data_2_avg(r,:,:),res)
  dum1 = gsn_add_polymarker(wks,cplot,lat2d(:,0),data_1a_avg(r,:),polyres1)
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
end do
print("Finished BL99 mean plots")

;;;;;;;;;;;;;
; Make winter diff figures
;;;;;;;;;;;;;
r = 0
do r = 0, nreg -1
   fout = outdir+"pi_physics_"+varcode1+"_"+reg_str(r)+"_regional_diff_transect"
  wks_type = "ps"
  wks = gsn_open_wks(wks_type,fout)
  gsn_define_colormap(wks,colors_diff)

  delete(res@cnLevels)
  res@cnLevels                    = diff_cntrs

  nplots = 1
  plot = new(nplots, graphic)
  n = 0
  res@gsnCenterString = pref2+" "+reg_str(r)+" sector"
  ;res@gsnCenterString = pref2+" "+reg_str(r)+" sector winter diff (MUSHY-BL99)"
  cplot = gsn_csm_contour(wks,diff(r,:,:),res)
  sig = gsn_csm_contour(wks,prob(r,:,:),res_95)
  overlay(cplot,sig)
  dum1 = gsn_add_polymarker(wks,cplot,lat2d(:,0),data_1a_avg(r,:),polyres1)
  dum2 = gsn_add_polymarker(wks,cplot,lat2d(:,0),data_2a_avg(r,:),polyres2)
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
end do
print("Finished difference plots")

print("Done with plots")

;;;;;;;;;;;;;;;;;;;;;; END script
end
