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
;varcode1 = "HMXL"

; set years we want to average over
yr_min = 871
yr_max = 970
;nyrs = (yr_max-yr_min)+1
nyrs = 100

;;;;;;;;;;
; Get paths to files
;;;;;;;;;;
print("Establishing data paths")
; cesm directory
;data_dir1 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case1+"_regional_timeseries/"
;data_dir2 = "/glade/p/cgd/ppc/duvivier/cesm2_antarctic_polynya/DATA/"+case2+"_regional_timeseries/"
data_dir1 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case1+"_regional_timeseries/"
data_dir2 = "/project/ppc/duvivier/projects/cesm2_antarctic_polynya/DATA/"+case2+"_regional_timeseries/"

;;;;;;;;;;
; Get CESM dimensions, etc
;;;;;;;;;;
nmons = 12
mons = ispan(1,nmons,1)
ntimes = nmons*nyrs

; how many regions
nreg = 6

;;;;;;;;;;
; Get data
;;;;;;;;;;
print("reading in processed data")
; case 1
fin11 = addfile(data_dir1+case1+"_all_"+varcode1+"_regional_avgs.nc","r")
; load data
var11 = fin11->data_all
; convert time arrays
yy_all_1 = fin11->yy_all
mm_all_1 = fin11->mm_all
yy_1 = onedtond(yy_all_1,(/nyrs,nmons/))
mm_1 = onedtond(mm_all_1,(/nyrs,nmons/))

; case 2
fin12 = addfile(data_dir2+case2+"_all_"+varcode1+"_regional_avgs.nc","r")
; load data
var12 = fin12->data_all
; convert time arrays
yy_all_2 = fin12->yy_all
mm_all_2 = fin12->mm_all
yy_2 = onedtond(yy_all_2,(/nyrs,nmons/))
mm_2 = onedtond(mm_all_2,(/nyrs,nmons/))

; convert to m from cm
var11 = var11/100.0
var11@units = "m"
var12 = var12/100.0
var12@units = "m"

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
data_1_avg = new((/nreg,nmons/),typeof(var11))
data_2_avg = new((/nreg,nmons/),typeof(var11))
sig = new((/nreg,nmons/),typeof(var11))

m = 0
do m = 0, nmons -1
   r = 0
   do r = 0, nreg -1
      inds1 = ind((mm_all_1.eq.mons(m)).and.((yy_all_1.ge.yr_min).and.(yy_all_1.le.yr_max)))
      inds2 = ind((mm_all_2.eq.mons(m)).and.((yy_all_2.ge.yr_min).and.(yy_all_2.le.yr_max)))
      tmp1 = ndtooned(var11(r,inds1))
      tmp2 = ndtooned(var12(r,inds2))
      data_1_avg(r,m) = avg(tmp1)
      data_2_avg(r,m) = avg(tmp2)
      sig(r,m) = 100.*(1.-ttest(avg(tmp1),variance(tmp1),dimsizes(tmp1),avg(tmp2),variance(tmp2),dimsizes(tmp2),True,False))
      delete([/tmp1,tmp2/])
   end do
end do

;;;;;;;;;;
; Get difference array for plotting
;;;;;;;;;;
diff = new((/nreg,2,nmons/),typeof(var11))
; mean differences (Mushy - BL99)
diff(:,0,:) = data_2_avg - data_1_avg
; now deal with significance
diff(:,1,:) = where((sig.ge.95.),diff(:,0,:),diff@_FillValue)

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
  ;opts@trYMinF                     = -1.5 ;-0.7
  ;opts@trYMaxF                     = 1.5 ;0.5
  opts@gsnYRefLine                 = 0.0
  opts@tiYAxisString               = units
  opts@tiXAxisFontHeightF          = 0.03              ; set font size for x axis
  opts@tmXBLabelFontHeightF        = 0.03              ; set tick mark label size for x axis
  opts@tiYAxisFontHeightF          = 0.03              ; set font size for y axis
  opts@tmYLLabelFontHeightF        = 0.03              ; set tick mark label size for y axis

; line resources for all ensembles
  opts@xyDashPattern               = 0
  opts@xyMarker                    = 1
  opts@xyMarkerSizeF               = 0.025
  opts@xyMarkLineMode              = "MarkLines"
  opts@xyLineColors                = (/"black","black"/)
  opts@xyLineThicknessF            = 1.5 
  opts@xyDashPattern               = 0
  opts@xyMarkers                   = (/1,9/)
  ;opts@xyMarkers                   = (/1,16/)
  opts@xyMarkerColors              = opts@xyLineColors

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
;;;;;;;;;;;;;
; Make case 1 annual mean figure
;;;;;;;;;;;;;
; set outdir and fout names
fout = outdir+case1+"_all_"+varcode1+"_avg_annual_cycle"
wks = gsn_open_wks("ps" ,fout)
plot = new(5, graphic)

   ; Ross coastal
   r = 0
   opts@tiMainString = "Ross Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_1_avg(r+1,:),opts)
   ; Admunsen Bellingshausen coastal
   r = 1
   opts@tiMainString = "Admundsen Bellingshausen Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_1_avg(r+1,:),opts)
   ; Weddell coastal
   r = 2
   opts@tiMainString = "Weddell Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_1_avg(r+1,:),opts)
   ; Indian coastal
   r = 3
   opts@tiMainString = "Indian Ocean Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_1_avg(r+1,:),opts)
   ; Pacific coastal
   r = 4
   opts@tiMainString = "Pacific Ocean Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_1_avg(r+1,:),opts)

  ; make panel plot
  gsn_panel(wks,(/plot/),(/2,3/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")
  print("Done with "+case1)
;;;;;;;;;;;;;
; Make case 2 annual mean figure
;;;;;;;;;;;;;
; set outdir and fout names
fout = outdir+case2+"_all_"+varcode1+"_avg_annual_cycle"
wks = gsn_open_wks("ps" ,fout)
plot = new(5, graphic)

   ; Ross coastal
   r = 0
   opts@tiMainString = "Ross Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_2_avg(r+1,:),opts)
   ; Admunsen Bellingshausen coastal
   r = 1
   opts@tiMainString = "Admundsen Bellingshausen Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_2_avg(r+1,:),opts)
   ; Weddell coastal
   r = 2
   opts@tiMainString = "Weddell Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_2_avg(r+1,:),opts)
   ; Indian coastal
   r = 3
   opts@tiMainString = "Indian Ocean Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_2_avg(r+1,:),opts)
   ; Pacific coastal
   r = 4
   opts@tiMainString = "Pacific Ocean Coast"
   plot(r) = gsn_csm_xy(wks,mons,data_2_avg(r+1,:),opts)

  ; make panel plot
  gsn_panel(wks,(/plot/),(/2,3/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")
  print("Done with "+case2)

return
;;;;;;;;;;;;;
; Make annual mean differences
;;;;;;;;;;;;;
; set outdir and fout names
fout = outdir+case1+"_"+case2+"_diffs_all_"+varcode1+"_avg_annual_cycle"
wks = gsn_open_wks("ps" ,fout)
plot = new(5, graphic)
  
; change a few resources
  ;opts@trYMinF                     = -0.4
  ;opts@trYMaxF                     = 0.6

   ; Ross coastal
   r = 0
   opts@tiMainString = "Ross Coast"
   plot(r) = gsn_csm_xy(wks,mons,diff(r+1,:,:),opts)
   ; Admunsen Bellingshausen coastal
   r = 1
   opts@tiMainString = "Admundsen Bellingshausen Coast"
   plot(r) = gsn_csm_xy(wks,mons,diff(r+1,:,:),opts)
   ; Weddell coastal
   r = 2
   opts@tiMainString = "Weddell Coast"
   plot(r) = gsn_csm_xy(wks,mons,diff(r+1,:,:),opts)
   ; Indian coastal
   r = 3
   opts@tiMainString = "Indian Ocean Coast"
   plot(r) = gsn_csm_xy(wks,mons,diff(r+1,:,:),opts)
   ; Pacific coastal
   r = 4
   opts@tiMainString = "Pacific Ocean Coast"
   plot(r) = gsn_csm_xy(wks,mons,diff(r+1,:,:),opts)

  ; make panel plot
  gsn_panel(wks,(/plot/),(/2,3/),pnlres)                   ; create final panel plot
  delete(wks)
  print("converting to png")
  system("convert -trim -border 5 -bordercolor white -density 300 "+  \
         "-rotate -90 -trim -border 5 "+fout+".ps "+fout+".png")
  print("Done with diffs")

print("Done with plots")

;;;;;;;;;;;;;;;;;;;;;; END script
end
