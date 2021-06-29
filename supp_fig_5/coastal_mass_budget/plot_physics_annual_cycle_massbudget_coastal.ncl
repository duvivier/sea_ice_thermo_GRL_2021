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
varcode4 = "meltt"
varcode5 = "meltb"
varcode6 = "meltl"
varcode7 = "dvidtd"
; one more than above because of net
nvars = 7

; set years we want to average over
yr_min = 871 
yr_max = 970

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
;data_dir1 = "/glade/p/cgd/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case1+"_regional_timeseries/"
;data_dir2 = "/glade/p/cgd/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case2+"_regional_timeseries/"
data_dir1 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case1+"_regional_timeseries/"
data_dir2 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case2+"_regional_timeseries/"

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
print("reading in processed data")
; case 1
fin11 = addfile(data_dir1+case1+"_coastal_"+varcode1+"_regional_avgs.nc","r")
fin21 = addfile(data_dir1+case1+"_coastal_"+varcode2+"_regional_avgs.nc","r")
fin31 = addfile(data_dir1+case1+"_coastal_"+varcode3+"_regional_avgs.nc","r")
fin41 = addfile(data_dir1+case1+"_coastal_"+varcode4+"_regional_avgs.nc","r")
fin51 = addfile(data_dir1+case1+"_coastal_"+varcode5+"_regional_avgs.nc","r")
fin61 = addfile(data_dir1+case1+"_coastal_"+varcode6+"_regional_avgs.nc","r")
fin71 = addfile(data_dir1+case1+"_coastal_"+varcode7+"_regional_avgs.nc","r")

; load data
var11 = fin11->data_all
var21 = fin21->data_all
var31 = fin31->data_all
var41 = -1.0*fin41->data_all
var51 = -1.0*fin51->data_all
var61 = -1.0*fin61->data_all
var71 = fin71->data_all
; get total mass budget
var81 = var11+var21+var31+var41+var51+var61+var71
; get total thermodynamic change
var91 = var11+var21+var31+var41+var51+var61
; convert time arrays
yy_all_1 = fin11->yy_all
mm_all_1 = fin11->mm_all
yy_1 = onedtond(yy_all_1,(/nyrs,nmons/))
mm_1 = onedtond(mm_all_1,(/nyrs,nmons/))

; case 2
fin12 = addfile(data_dir2+case2+"_coastal_"+varcode1+"_regional_avgs.nc","r")
fin22 = addfile(data_dir2+case2+"_coastal_"+varcode2+"_regional_avgs.nc","r")
fin32 = addfile(data_dir2+case2+"_coastal_"+varcode3+"_regional_avgs.nc","r")
fin42 = addfile(data_dir2+case2+"_coastal_"+varcode4+"_regional_avgs.nc","r")
fin52 = addfile(data_dir2+case2+"_coastal_"+varcode5+"_regional_avgs.nc","r")
fin62 = addfile(data_dir2+case2+"_coastal_"+varcode6+"_regional_avgs.nc","r")
fin72 = addfile(data_dir2+case2+"_coastal_"+varcode7+"_regional_avgs.nc","r")

; load data
var12 = fin12->data_all
var22 = fin22->data_all
var32 = fin32->data_all
var42 = -1.0*fin42->data_all
var52 = -1.0*fin52->data_all
var62 = -1.0*fin62->data_all
var72 = fin72->data_all
; get total mass budget
var82 = var12+var22+var32+var42+var52+var62+var72
; get total thermodynamic change
var92 = var12+var22+var32+var42+var52+var62
; convert time arrays
yy_all_2 = fin12->yy_all
mm_all_2 = fin12->mm_all
yy_2 = onedtond(yy_all_2,(/nyrs,nmons/))
mm_2 = onedtond(mm_all_2,(/nyrs,nmons/))

; set units
units = var11@units

;;;;;;;;;;
; Get monthly averages and check if significantly different
;;;;;;;;;;
print("Calculating monthly means and checking significance")
; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level

; loop through months, decades, regions
data_1_avg = new((/nreg,nvars,nmons/),typeof(var11))
data_2_avg = new((/nreg,nvars,nmons/),typeof(var11))
data_1_std = new((/nreg,nvars,nmons/),typeof(var11))
data_2_std = new((/nreg,nvars,nmons/),typeof(var11))
sig = new((/nreg,nvars,nmons/),typeof(var11))

m = 0
do m = 0, nmons -1
   r = 0
   do r = 0, nreg -1
      ; case1
      inds1 = ind((mm_all_1.eq.mons(m)).and.((yy_all_1.ge.yr_min).and.(yy_all_1.le.yr_max)))
      inds2 = ind((mm_all_2.eq.mons(m)).and.((yy_all_2.ge.yr_min).and.(yy_all_2.le.yr_max)))
      ; total mass budget
      v = 0
      tmp1 = ndtooned(var81(r,inds1))
      tmp2 = ndtooned(var82(r,inds2))
      data_1_avg(r,v,m) = avg(tmp1)
      data_2_avg(r,v,m) = avg(tmp2)
      data_1_std(r,v,m) = stddev(tmp1)
      data_2_std(r,v,m) = stddev(tmp2)
      sig(r,v,m) = 100.*(1.-ttest(avg(tmp1),variance(tmp1),dimsizes(tmp1),avg(tmp2),variance(tmp2),dimsizes(tmp2),True,False))
      delete([/tmp1,tmp2/])
      ; total thermo
      v = 1
      tmp1 = ndtooned(var91(r,inds1))
      tmp2 = ndtooned(var92(r,inds2))
      data_1_avg(r,v,m) = avg(tmp1)
      data_2_avg(r,v,m) = avg(tmp2)
      data_1_std(r,v,m) = stddev(tmp1)
      data_2_std(r,v,m) = stddev(tmp2)
      sig(r,v,m) = 100.*(1.-ttest(avg(tmp1),variance(tmp1),dimsizes(tmp1),avg(tmp2),variance(tmp2),dimsizes(tmp2),True,False))
      delete([/tmp1,tmp2/])
      ; congel
      v = 2
      tmp1 = ndtooned(var11(r,inds1))
      tmp2 = ndtooned(var12(r,inds2))
      data_1_avg(r,v,m) = avg(tmp1)
      data_2_avg(r,v,m) = avg(tmp2)
      data_1_std(r,v,m) = stddev(tmp1)
      data_2_std(r,v,m) = stddev(tmp2)
      sig(r,v,m) = 100.*(1.-ttest(avg(tmp1),variance(tmp1),dimsizes(tmp1),avg(tmp2),variance(tmp2),dimsizes(tmp2),True,False))
      delete([/tmp1,tmp2/])
      ; frazil
      v = 3
      tmp1 = ndtooned(var21(r,inds1))
      tmp2 = ndtooned(var22(r,inds2))
      data_1_avg(r,v,m) = avg(tmp1)
      data_2_avg(r,v,m) = avg(tmp2)
      data_1_std(r,v,m) = stddev(tmp1)
      data_2_std(r,v,m) = stddev(tmp2)
      sig(r,v,m) = 100.*(1.-ttest(avg(tmp1),variance(tmp1),dimsizes(tmp1),avg(tmp2),variance(tmp2),dimsizes(tmp2),True,False))
      delete([/tmp1,tmp2/])
      ; snowice
      v = 4
      tmp1 = ndtooned(var31(r,inds1))
      tmp2 = ndtooned(var32(r,inds2))
      data_1_avg(r,v,m) = avg(tmp1)
      data_2_avg(r,v,m) = avg(tmp2)
      data_1_std(r,v,m) = stddev(tmp1)
      data_2_std(r,v,m) = stddev(tmp2)
      sig(r,v,m) = 100.*(1.-ttest(avg(tmp1),variance(tmp1),dimsizes(tmp1),avg(tmp2),variance(tmp2),dimsizes(tmp2),True,False))
      delete([/tmp1,tmp2/])
      ; bottom melt
      v = 5
      tmp1 = ndtooned(var51(r,inds1))
      tmp2 = ndtooned(var52(r,inds2))
      data_1_avg(r,v,m) = avg(tmp1)
      data_2_avg(r,v,m) = avg(tmp2)
      data_1_std(r,v,m) = stddev(tmp1)
      data_2_std(r,v,m) = stddev(tmp2)
      sig(r,v,m) = 100.*(1.-ttest(avg(tmp1),variance(tmp1),dimsizes(tmp1),avg(tmp2),variance(tmp2),dimsizes(tmp2),True,False))
      delete([/tmp1,tmp2/])
      ; dynamics
      v = 6
      tmp1 = ndtooned(var71(r,inds1))
      tmp2 = ndtooned(var72(r,inds2))
      data_1_avg(r,v,m) = avg(tmp1)
      data_2_avg(r,v,m) = avg(tmp2)
      data_1_std(r,v,m) = stddev(tmp1)
      data_2_std(r,v,m) = stddev(tmp2)
      sig(r,v,m) = 100.*(1.-ttest(avg(tmp1),variance(tmp1),dimsizes(tmp1),avg(tmp2),variance(tmp2),dimsizes(tmp2),True,False))
      delete([/tmp1,tmp2/])
   end do
end do

;;;;;;;;;;
; Get a +/1 2sigma array for plotting
;;;;;;;;;;
data_2_2sig = new((/nreg,nvars*2,nmons/),typeof(var11))
data_2_2sig(:,0:(nvars-1),:) = data_2_avg + (2*(data_2_std))
data_2_2sig(:,nvars:(nvars*2-1),:) = data_2_avg - (2*(data_2_std))

;;;;;;;;;;
; Get difference array for plotting
;;;;;;;;;;
diff = new((/nreg,nvars*2,nmons/),typeof(var11))
; mean differences (Mushy - BL99)
diff(:,0:(nvars-1),:) = data_2_avg - data_1_avg
; now deal with significance
diff(:,nvars:(nvars*2-1),:) = where((sig.ge.95.),diff(:,0:(nvars-1),:),diff@_FillValue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; Plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print("Plotting mass budget")
   outdir = "./"

label_inds = mons
label_str = tostring(mons)

;;;;;;;;;;;;;
; XY Cross Section Plot Resources
;;;;;;;;;;;;;
; set resources used for all plots
  opts                             = True              ; plot mods desired
  opts@gsnFrame                    = False             ; Do not draw plot 
  opts@gsnDraw                     = False             ; Do not advance frame
  opts@gsnMaximize                 = False             ; maximize plot
  opts@vpWidthF                    = 1.25  ;1.75               ; set aspect ratio
  opts@vpHeightF                   = 0.75 ;0.55

; set title resources
  opts@tiMainOn                    = True
  opts@tiMainFontHeightF           = 0.03
  opts@tiMainFontThicknessF        = 1.5
  opts@tmXBMode                    = "Explicit"
  opts@tmXBValues                  = label_inds
  opts@tmXBLabels                  = label_str
  opts@trXMinF                     = 1
  opts@trXMaxF                     = 12
  opts@tiXAxisString               = "Month"
  opts@tmXBOn                      = True
  opts@trYMinF                     = -1.5
  opts@trYMaxF                     = 1.5
  opts@gsnYRefLine                 = 0.0
  opts@tiYAxisString               = units
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis

; line resources for all variables
  opts@xyMarkLineMode              = "MarkLines"
  opts@xyDashPattern               = 0
  opts@xyLineColors                = (/"black","grey50","steelblue1","blue2","darkturquoise","firebrick","goldenrod","black","grey50","steelblue1","blue2","darkturquoise","firebrick","goldenrod"/)
  opts@xyLineThicknessF            = 2.0
  opts@xyDashPattern               = 0
  opts@xyMarkerSizeF               = 0.025
  opts@xyMarkers                   = (/1,1,1,1,1,1,1,9,9,9,9,9,9,9/)
  opts@xyMarkerColors              = opts@xyLineColors

; line resources for two sigma
  opts2sig                         = opts
  opts2sig@xyMarkLineMode          = "Markers"
  delete(opts2sig@xyMarkers)
  opts2sig@xyMarker                = 1
  opts2sig@xyLineColors            = (/"black","grey50","steelblue1","blue2","darkturquoise","firebrick","goldenrod","black","grey50","steelblue1","blue2","darkturquoise","firebrick","goldenrod"/)

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 3                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 3                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True	       ; fill the page

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make actual plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; region labels
reg_str = (/"Ross", "Amundsen Bellingshausen", "Weddell", "Indian Ocean", "Pacific Ocean"/)
nplots = dimsizes(reg_str)

;;;;;;;;;;;;;
; Make all regions mean state plot
;;;;;;;;;;;;;
; set outdir and fout names
fout = outdir+case2+"_all_regions_coastal_mass_budget_avg_annual_cycle"
wks = gsn_open_wks("ps" ,fout)
plot = new(nplots, graphic)
   opts@trYMinF                     = -1.5
   opts@trYMaxF                     = 1.5
   r = 0
   do r = 0, nplots -1
      opts@tiMainString = "MUSHY monthly "+reg_str(r)+" coastal mean"
      plot(r) = gsn_csm_xy(wks,mons,data_2_avg(r+1,:,:),opts)
   end do
; make panel plot
gsn_panel(wks,(/plot/),(/nplots,1/),pnlres)                   ; create final panel plot
delete(wks)
print("converting to png")
system("convert -trim -border 5 -bordercolor white -density 300 "+  \
       "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Done with mean stateplots")
delete([/plot,fout/])

;;;;;;;;;;;;;
; Make all regions difference plot
;;;;;;;;;;;;;
; set outdir and fout names
fout = outdir+case1+"_"+case2+"_diff_all_regions_coastal_mass_budget_avg_annual_cycle"
wks = gsn_open_wks("ps" ,fout)
plot = new(nplots, graphic)
   opts@trYMinF                     = -0.4
   opts@trYMaxF                     = 0.6
   r = 0
   do r = 0, nplots -1
      opts@tiMainString = "Monthly "+reg_str(r)+" coastal mean difference"
      plot(r) = gsn_csm_xy(wks,mons,diff(r+1,:,:),opts)
   end do
; make panel plot
gsn_panel(wks,(/plot/),(/nplots,1/),pnlres)                   ; create final panel plot
delete(wks)
print("converting to png")
system("convert -trim -border 5 -bordercolor white -density 300 "+  \
       "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")

print("Done with diff plots")
delete([/plot,fout/])

;;;;;;;;;;;;;;;;;;;;;; END script
end
