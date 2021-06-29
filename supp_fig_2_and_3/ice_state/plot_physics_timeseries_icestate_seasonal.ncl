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

varcode1 = "aice"
varcode2 = "hi"

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
;data_dir1 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case1+"_regional_timeseries/"
;data_dir2 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case2+"_regional_timeseries/"

data_dir1 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case1+"_regional_timeseries/"
data_dir2 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case2+"_regional_timeseries/"


;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
nmons = 12
nyrs = 100
ntimes = nmons*nyrs 

; how many regions
reg_str = (/"all","Ross","BAm","Wed","Ind","Pac"/)
nreg = dimsizes(reg_str)

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("reading in processed data")
; code is that first # is the variable name and second # is case name

; case 1
fin11 = addfile(data_dir1+case1+"_coastal_"+varcode1+"_regional_totals.nc","r")
fin21 = addfile(data_dir1+case1+"_coastal_"+varcode2+"_regional_totals.nc","r")
; read in data
var11 = fin11->data_all
var21 = fin21->data_all
; convert time arrays
yy_all_1 = fin11->yy_all
mm_all_1 = fin11->mm_all
yy_1 = onedtond(yy_all_1,(/nyrs,nmons/))
mm_1 = onedtond(mm_all_1,(/nyrs,nmons/))

; case 2
fin12 = addfile(data_dir2+case2+"_coastal_"+varcode1+"_regional_totals.nc","r")
fin22 = addfile(data_dir2+case2+"_coastal_"+varcode2+"_regional_totals.nc","r")
; read in data
var12 = fin12->data_all
var22 = fin22->data_all
; convert time arrays
yy_all_2 = fin12->yy_all
mm_all_2 = fin12->mm_all
yy_2 = onedtond(yy_all_2,(/nyrs,nmons/))
mm_2 = onedtond(mm_all_2,(/nyrs,nmons/))

; set units
units1 = var11@units
units2 = var21@units

;;;;;;;;;;
; Get seasonal averages
;;;;;;;;;;
seas_str = (/"AMJ","JAS","AMJJAS"/)
nseas = dimsizes(seas_str)
print("Calculating seasonal means")
; seasonal arrays
data_1 = new((/nseas,nreg,2,nyrs/),typeof(var11))
data_2 = new((/nseas,nreg,2,nyrs/),typeof(var11))

; loop through years and calculate seasonal means
y = 0
do y = 0, nyrs-2
   ; case1
   c = 0
   ; AMJ
   s = 0
   inds = ind(((mm_all_1.ge.4).and.(mm_all_1.le.6)).and.(yy_all_1.eq.yy_1(y,0)))
   data_1(s,:,c,y) = dim_avg_n(var11(:,inds),1)
   data_2(s,:,c,y) = dim_avg_n(var21(:,inds),1)
   delete(inds)
   ; JAS
   s = 1
   inds = ind(((mm_all_1.ge.7).and.(mm_all_1.le.9)).and.(yy_all_1.eq.yy_1(y,0)))
   data_1(s,:,c,y) = dim_avg_n(var11(:,inds),1)
   data_2(s,:,c,y) = dim_avg_n(var21(:,inds),1)
   delete(inds)   
   ; AMJJAS
   s = 2
   inds = ind(((mm_all_1.ge.4).and.(mm_all_1.le.9)).and.(yy_all_1.eq.yy_1(y,0)))
   data_1(s,:,c,y) = dim_avg_n(var11(:,inds),1)
   data_2(s,:,c,y) = dim_avg_n(var21(:,inds),1)
   delete(inds)

   ; case2
   c = 1
   ; AMJ
   s = 0
   inds = ind(((mm_all_2.ge.4).and.(mm_all_2.le.6)).and.(yy_all_2.eq.yy_2(y,0)))
   data_1(s,:,c,y) = dim_avg_n(var12(:,inds),1)
   data_2(s,:,c,y) = dim_avg_n(var22(:,inds),1)
   delete(inds)
   ; JAS
   s = 1
   inds = ind(((mm_all_2.ge.7).and.(mm_all_2.le.9)).and.(yy_all_2.eq.yy_2(y,0)))
   data_1(s,:,c,y) = dim_avg_n(var12(:,inds),1)
   data_2(s,:,c,y) = dim_avg_n(var22(:,inds),1)
   delete(inds)  
   ; AMJJAS
   s = 2
   inds = ind(((mm_all_2.ge.4).and.(mm_all_2.le.9)).and.(yy_all_2.eq.yy_2(y,0)))
   data_1(s,:,c,y) = dim_avg_n(var12(:,inds),1)
   data_2(s,:,c,y) = dim_avg_n(var22(:,inds),1)
   delete(inds)
end do

;;;;;;;;;;
; Find if mean differences are significant
;;;;;;;;;;
print("Calculating significance")
; Use student's t-test to see if ensemble means differ
; if probability is less than 0.05 then we know at 95% confidence level the means are statistically different
; with equation below, if values are >= 95, we know they're different at 95% level

s = 0
do s = 0, nseas -1
  r = 0
  do r = 0, nreg-1
   c = 0
   ; varcode1
   prob = 100.*(1-ttest(avg(data_1(s,r,c,:)),variance(data_1(s,r,c,:)),dimsizes(data_1(s,r,c,:)),avg(data_1(s,r,c+1,:)),variance(data_1(s,r,c+1,:)),dimsizes(data_1(s,r,c+1,:)),True,False))
   print(varcode1+" "+seas_str(s)+" "+reg_str(r)+" region, probability average is different: "+tostring(prob))
   print(case1+" avg: "+tostring(avg(data_1(s,r,c,:))))
   print(case2+" avg: "+tostring(avg(data_1(s,r,c+1,:))))
   delete(prob)
   prob = 100.*(1-ftest(variance(data_1(s,r,c,:)),dimsizes(data_1(s,r,c,:)),variance(data_1(s,r,c+1,:)),dimsizes(data_1(s,r,c+1,:)),0))
   print(varcode1+" "+seas_str(s)+" "+reg_str(r)+" region, probability variance is different: "+tostring(prob))
   print(case1+" std dev: "+tostring(stddev(data_1(s,r,c,:))))
   print(case2+" std dev: "+tostring(stddev(data_1(s,r,c+1,:))))
   delete(prob)

   ; varcode2
   prob = 100.*(1-ttest(avg(data_2(s,r,c,:)),variance(data_2(s,r,c,:)),dimsizes(data_2(s,r,c,:)),avg(data_2(s,r,c+1,:)),variance(data_2(s,r,c+1,:)),dimsizes(data_2(s,r,c+1,:)),True,False))
   print(varcode2+" "+seas_str(s)+" "+reg_str(r)+" region, probability average is different: "+tostring(prob))
   print(case1+" avg: "+tostring(avg(data_2(s,r,c,:))))
   print(case2+" avg: "+tostring(avg(data_2(s,r,c+1,:))))
   delete(prob)
   prob = 100.*(1-ftest(variance(data_2(s,r,c,:)),dimsizes(data_2(s,r,c,:)),variance(data_2(s,r,c+1,:)),dimsizes(data_2(s,r,c+1,:)),0))
   print(varcode2+" "+seas_str(s)+" "+reg_str(r)+" region, probability variance is different: "+tostring(prob))
   print(case1+" std dev: "+tostring(stddev(data_2(s,r,c,:))))
   print(case2+" std dev: "+tostring(stddev(data_2(s,r,c+1,:))))
  delete(prob)

   end do
end do

;;;;;;;;;;
; Get running means
;;;;;;;;;;
nyrsavg = 10
data_1_ave = runave_n_Wrap(data_1,nyrsavg,0,3)
data_2_ave = runave_n_Wrap(data_2,nyrsavg,0,3)

;;;;;;;;;;
; combine to plot
;;;;;;;;;;
data_1_plot = new((/nseas,nreg,4,nyrs/),typeof(var11))
data_2_plot = new((/nseas,nreg,4,nyrs/),typeof(var11))

data_1_plot(:,:,0:1,:) = data_1
data_1_plot(:,:,2:3,:) = data_1_ave
data_2_plot(:,:,0:1,:) = data_2
data_2_plot(:,:,2:3,:) = data_2_ave

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
  opts@vpWidthF                    = 1.25  ;1.75               ; set aspect ratio
  opts@vpHeightF                   = 0.75 ;0.55

; set title resources
  opts@tiMainOn                    = True
  opts@tiMainFontHeightF           = 0.03
  opts@tiMainFontThicknessF        = 1.5
  ;opts@trXMinF                     = 0
  ;opts@trXMaxF                     = 50
  opts@gsnYRefLine                 = 0.0
  opts@tmXBOn                      = True
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025
  opts@xyMarkLineMode              = "Lines"
  opts@xyLineColors                = (/"grey","lightgoldenrod","black","goldenrod"/)
  opts@xyLineThicknesses           = (/1.5,1.5,4.0,4.0/)
  opts@xyDashPattern               = 0

; for the ice state plots
  opts1 = opts
  opts1@tiYAxisString               = "Total ice area "+units1
  opts2 = opts
  opts2@tiYAxisString               = "Total ice volume "+units2

; Set panel resources
  pnlres                            = True
  pnlres@txFontHeightF              = 0.013
  pnlres@gsnPaperOrientation        = "landscape"      ; landscape mode
  pnlres@gsnMaximize                = True             ; fill the page
  pnlres@gsnPanelYWhiteSpacePercent = 1.5                ; Add white space b/w plots.
  pnlres@gsnPanelXWhiteSpacePercent = 1.5                ; Add white space b/w plots.
  pnlres@gsnPanelScalePlotIndex     = 1
  pnlres@gsnMaximize                = True	       ; fill the page

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make actual plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;loop through seasons
s = 0
do s = 0, nseas -1
print("Plotting "+seas_str(s))
pref = "Regional transient total "+varcode1+" "+seas_str(s)
fout = outdir+"/pi_physics_coastal_"+varcode1+"_"+seas_str(s)
 wks = gsn_open_wks("ps" ,fout)
 plot = new(6, graphic)
   ; Ice state
   r = 0
   do r = 0, nreg -1
      opts1@tiMainString = reg_str(r)+" Coast"
      plot(r) = gsn_csm_xy(wks,yy_1(:,0),data_1_plot(s,r,:,:),opts1)
   end do
  ; make panel plot
  pnlres@txString  = pref
  gsn_panel(wks,(/plot/),(/3,2/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")
  delete([/fout,plot/])

end do

;loop through seasons
s = 0
do s = 0, nseas -1
print("Plotting "+seas_str(s))
pref = "Regional transient total "+varcode2+" "+seas_str(s)
fout = outdir+"/pi_physics_coastal_"+varcode2+"_"+seas_str(s)
 wks = gsn_open_wks("ps" ,fout)
 plot = new(6, graphic)
   r = 0
   do r = 0, nreg -1
      opts2@tiMainString = reg_str(r)+" Coast"
      plot(r) = gsn_csm_xy(wks,yy_1(:,0),data_2_plot(s,r,:,:),opts2)
   end do
  ; make panel plot
  pnlres@txString  = pref
  gsn_panel(wks,(/plot/),(/3,2/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")
  delete([/fout,plot/])

end do

print("Done with seasonal plots")

;;;;;;;;;;;;;;;;;;;;;; END script
end
