; script from SYeager to calculate MOC in sigma2 coordinates from
; single variable time series. Using for CESM2 CMIP6 runs...
; LLandrum
; June 2019
; ==============================================================
; NCL script to compute POP MOC field in sigma2 coordinates
; offline from POP netcdf history files.
;
; NOTE:  designed for a collection of POP monthly timeseries files 
;	containing single variables.
; ==============================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/glade/u/home/yeager/ncl/yeager_util.ncl"

;external MOC "/glade/u/home/yeager/POP_tools/computeMOCoffline/lib/MOCoffline.so"
external MOC "/glade/u/home/yeager/POP_tools/computeMOCoffline/lib/MOCoffline.new.so"

begin
  TimeDate = systemfunc("date")

  indir  = getenv("indir")
  outdir  = getenv("outdir")
  casename  = getenv("casename")
  timestamp  = getenv("timestamp")
  year  = tointeger(getenv("year"))

; casename = "b.e21.B1850.f09_g17.CMIP6-piControl.001"
; indir = "/glade/collections/cdg/timeseries-cmip6/"+casename+"/ocn/proc/tseries/month_1/"
; outdir = "/glade/scratch/yeager/g.DPLE.GECOIAF.T62_g16.009.chey/"
; timestamp = "070001-079912"
; year = 700

  fin_temp = casename+".pop.h.XXXX."+timestamp+".nc"
  yrstr = sprinti("%0.4i", year)
  fout = outdir+casename+".pop.h."+yrstr+"01-"+yrstr+"12.MOCsig2.nc"
  print("generating "+fout+"...")

; ==============================================================
; Get MOC template from POP file
; ==============================================================
  fin = indir+str_sub_str(fin_temp,"XXXX","MOC")
  f        = addfile (fin, "r")
  time = f->time
  simyear = time/365.
  t0 = min(ind(simyear.gt.year))
  t1 = t0+11
  delete(time)
  delete(simyear)
  MOCz   = f->MOC(t0:t1,:,0,:,:)            ; (nt,2,3,61,395)
  MOCz   = 0.0         ; keep meta data
  transport_regions    = f->transport_regions
  lat_aux_grid         = todouble(f->lat_aux_grid)
  dims  = dimsizes(transport_regions)  ; (0)=2 , (1)=256
  ntr = dims(0)
  delete(dims)
; lat_aux_grid = todouble(fspan(-89.5,89.5,359))
  nyaux = dimsizes(lat_aux_grid)
; ==============================================================
; Define sigma coordinates  
; ==============================================================
    ref_depth = 2000.   ; in m

      tmp1 = array_append_record( 28.+ispan(0,34,1)*0.2, 35.+ispan(0,9,1)*0.1, 0)
      tmp2 = array_append_record( tmp1, 36.+ispan(0,40,1)*0.05, 0)
    SIGMA2 = tmp2      ; [kg/m^3]
      delete(tmp1)
      delete(tmp2)

    SIGMA2@units = "kg/m^3"
    SIGMA2@reference_depth = sprintf("%4.0f",ref_depth)+"m"
    nsig = dimsizes(SIGMA2)

    SIGMA2_top=SIGMA2
    tmp = 0.5*(SIGMA2(1:nsig-1)-SIGMA2(0:nsig-2))
    SIGMA2_top(0) = SIGMA2(0)-tmp(0)
    SIGMA2_top(1:nsig-1) = SIGMA2(0:nsig-2)+tmp

    SIGMA2_bot=SIGMA2
    SIGMA2_bot(nsig-1) = SIGMA2(nsig-1)+tmp(nsig-2)
    SIGMA2_bot(0:nsig-2)= SIGMA2(0:nsig-2)+tmp
    delete(tmp)

; ==============================================================
; Read in other needed variables from POP output netcdf
; ==============================================================
  f        = addfile (indir+str_sub_str(fin_temp,"XXXX","PD"), "r")
  pd     = f->PD(t0:t1,:,:,:)
  f        = addfile (indir+str_sub_str(fin_temp,"XXXX","TEMP"), "r")
  temp   = f->TEMP(t0:t1,:,:,:)
  f        = addfile (indir+str_sub_str(fin_temp,"XXXX","SALT"), "r")
  salt   = f->SALT(t0:t1,:,:,:)
  f        = addfile (indir+str_sub_str(fin_temp,"XXXX","UVEL"), "r")
  u_e   = f->UVEL(t0:t1,:,:,:)
  f        = addfile (indir+str_sub_str(fin_temp,"XXXX","VVEL"), "r")
  v_e   = f->VVEL(t0:t1,:,:,:)
  f        = addfile (indir+str_sub_str(fin_temp,"XXXX","WVEL"), "r")
  w_e   = f->WVEL(t0:t1,:,:,:)
  ulat   = f->ULAT
  ulon   = f->ULONG
  tlat   = f->TLAT
  tlon   = f->TLONG
  angle  = f->ANGLE
  kmt  = f->KMT
  kmu  = f->KMU
  dxu    = f->DXU
  dyu    = f->DYU
  rmask  = f->REGION_MASK
  tarea  = f->TAREA
  uarea  = f->UAREA
  time   = f->time(t0:t1)
  simyear = time/365.
  z_t   = f->z_t
  z_w   = f->z_w
  dz   = f->dz
  dzw   = f->dzw
  hu   = f->HU
  ht   = f->HT
  delete(f)
  if (isatt(salt,"scale_factor")) then
    salt = salt*salt@scale_factor
  end if

; ==============================================================
; U-grid mask
; ==============================================================
  umask = ismissing(u_e)

; ==============================================================
; dimensions
; ==============================================================
  dims = dimsizes(temp)
  nt = dims(0)
  nz = dims(1)
  ny = dims(2)
  nx = dims(3)
  delete(dims)
  km     = max(kmt)

; ==============================================================
; convert to mks units
; ==============================================================
  dz = dz/100.
  dz@units = "meters"
  dzw = dzw/100.
  dzw@units = "meters"
  z_t = z_t/100.
  z_t@units = "meters"
  z_w = z_w/100.
  z_w@units = "meters"
  u_e = u_e/100.
  u_e@units = "m/s"
  v_e = v_e/100.
  v_e@units = "m/s"
  w_e = w_e/100.
  w_e@units = "m/s"
  dxu = dxu/100.
  dxu@units = "m"
  dyu = dyu/100.
  dyu@units = "m"
  hu = hu/100.
  hu@units = "m"
  ht = ht/100.
  ht@units = "m"

  uarea = uarea/100./100.
  uarea@units = "m^2"
  tarea = tarea/100./100.
  tarea@units = "m^2"
  z_bot=z_w
  z_bot=z_w+dz
  z_top=z_w

; ==============================================================
; Compute PD referenced to ref_depth on T grid
; ==============================================================
;   TimeDate = systemfunc("date")
    pdt = new((/nt,nz,ny,nx/),typeof(pd),pd@_FillValue)
    do it=0,nt-1
    do iz=0,nz-1
;      tmp = rho_mwjf(todouble(temp(it,iz,:,:)),todouble(salt(it,iz,:,:)),ref_depth)
       tmp = rho_mwjf(temp(it,iz,:,:),salt(it,iz,:,:),ref_depth)
       pdt(it,iz,:,:) = tmp
       delete(tmp)
    end do
    end do
    pdt@coordinates = temp@coordinates
    pdt@long_name = "Potential Density Ref to "+ref_depth+"m"
    pdt = (pdt-1.)*1000.
    pdt@units = "kg/m^3"

; ==============================================================
; Regrid PDT onto U,V grid by simply averaging four surrounding points
; ==============================================================
    pdu = pdt
    pdu@coordinates = u_e@coordinates
    tmppdt = new((/4,nt,nz,ny,nx/),typeof(pdt),pdt@_FillValue)
    tmppdt(0,:,:,:,:) = (/ pdt /)                         ; PDT(i,j)
    tmppdt(1,:,:,:,:) = (/ eoshift(pdt,2,1) /)            ; PDT(i,j+1)
    tmpip1 = eoshift(pdt,1,1)
    tmppdt(2,:,:,:,:) = (/ tmpip1 /)                      ; PDT(i+1,j)
    tmppdt(3,:,:,:,:) = (/ eoshift(tmpip1,2,1) /)         ; PDT(i+1,j+1)
    pdu = (/ dim_avg_n_Wrap(tmppdt,0) /)
    delete(tmppdt)
    pdu = where(umask,pdu@_FillValue,pdu)

; ==============================================================
; Remove any inversions in PDU field & ensure monotonicity
; ==============================================================
    do it=0,nt-1
    do izz=1,nz-1
       iz = nz-izz
       tmp1 = pdu(it,iz,:,:)
       tmp2 = pdu(it,iz-1,:,:)
       tmp2 = where((.not.ismissing(tmp1)).and.(.not.ismissing(tmp2)).and.(tmp2.gt.tmp1),tmp1,tmp2)
       pdu(it,iz-1,:,:) = (/ tmp2-1.e-5 /)
       delete(tmp1)
       delete(tmp2)
    end do
    end do
    pdu = where(umask,pdu@_FillValue,pdu)
;   wallClockElapseTime(TimeDate, "finished computing PDU",0)

; ==============================================================
; Find isopycnal layer depths on U-grid
; ==============================================================
;   TimeDate = systemfunc("date")
    zsigu_top = new((/nt,nsig,ny,nx/),typeof(z_t),1.e30)
    zsigu_top!1 = "sigma2"
    zsigu_top&sigma2=SIGMA2
    zsigu_top@long_name="Depth of top of sigma layer"
    zsigu_top@units="m"
    zsigu_bot = new((/nt,nsig,ny,nx/),typeof(z_t),1.e30)
    zsigu_bot!1 = "sigma2"
    zsigu_bot&sigma2=SIGMA2
    zsigu_bot@long_name="Depth of bottom of sigma layer"
    zsigu_bot@units="m"

;; offline routine to compute depths:
    MOC::sig2z(nt,nz,ny,nx,nsig,pdu,pdu@_FillValue,z_t,z_bot,SIGMA2_top, \
               SIGMA2_bot,zsigu_top,zsigu_bot,zsigu_top@_FillValue)
;   wallClockElapseTime(TimeDate, "computed zsigu_top,zsigu_bot",0)

; ==============================================================
; Calculate volume fluxes: UDYDZ & VDXDZ (m^3/s) for each layer in 
; z_coord on U-grid
; ==============================================================
; TimeDate = systemfunc("date")
  uedydz = u_e*conform(u_e,dyu,(/2,3/))*conform(u_e,dz,(/1/))  ;m^3/s
  vedxdz = v_e*conform(v_e,dxu,(/2,3/))*conform(v_e,dz,(/1/))  ;m^3/s
; wallClockElapseTime(TimeDate, "computed uedydz, etc",0)

; ==============================================================
; Bin volume fluxes into sigma coordinates
; ==============================================================
  wedxdy_sig= new((/nt,nsig,ny,nx/),double,1.d30)       ; defined on T-grid
  vedxdz_sig= new((/nt,nsig,ny,nx/),double,1.d30)       ; defined on U-grid
  uedydz_sig= new((/nt,nsig,ny,nx/),double,1.d30)       ; defined on U-grid

; TimeDate = systemfunc("date")
  MOC::sig2fluxconv(nt,nz,ny,nx,nsig,kmt,z_top,z_bot,dz,zsigu_top,zsigu_bot, \
               uedydz,vedxdz,uedydz@_FillValue,uedydz_sig,vedxdz_sig,wedxdy_sig,vedxdz_sig@_FillValue)
; wallClockElapseTime(TimeDate, "finished sig2fluxconv",0)

  wedxdy_sig_sumz = dim_sum_n_Wrap(wedxdy_sig,1)

; ==============================================================
; Compute WDXDY (m^3/s) in sigma-coord by cumulative integral
; of dWDXDY obtained from sig2fluxconv
; ==============================================================
; TimeDate = systemfunc("date")
;; need dim_cumsum with option "1", because otherwise can get 
;; non-zero W even where sigma surface is undefined
  workw = dim_cumsum_n(wedxdy_sig(:,(nsig-1):0,:,:),1,1)
  wedxdy_sig = (/ workw(:,(nsig-1):0,:,:) /)
; wallClockElapseTime(TimeDate, "finished computed wedxdy_sig, etc",0)

; ==============================================================
; compute MOC using w-method: by integrating zonally, then meridionally
; ==============================================================
; TimeDateMOC = systemfunc("date")

  MOCnew = new((/nt,ntr,nsig,nyaux/),typeof(MOCz),MOCz@_FillValue)
  MOCnew = 0.
  MOCnew!0 = MOCz!0             ; time
  MOCnew!1 = MOCz!1             ; transport_reg
  MOCnew!2 = "sigma"            ; sigma-coord
  MOCnew!3 = MOCz!3             ; latitude
  MOCnew&time = MOCz&time
  MOCnew&lat_aux_grid = lat_aux_grid
  MOCnew&sigma = SIGMA2
  MOCnew@units = "Sv"

  REGION_MASK_LAT_AUX        = new((/ntr,ny,nx/),integer)
  REGION_MASK_LAT_AUX(0,:,:) = rmask
  REGION_MASK_LAT_AUX(1,:,:) = rmask
  REGION_MASK_LAT_AUX(0,:,:) = where(rmask.gt.0, 1, 0)
  REGION_MASK_LAT_AUX(1,:,:) = where(rmask.ge.6.and.rmask.le.11, 1, 0)

  TMP1  = new((/nt,ntr,nsig,nyaux/), double, 1d30)
  TMP2  = new((/nt,ntr,nsig,nyaux/), double, 1d30)
  TMP3  = new((/nt,ntr,nsig,nyaux/), double, 1d30)

; ==============================================================
; Compute zonal integral for all T-cells within a latitude band
; ==============================================================
; Use below for MOCoffline.so:
; MOC::moczonalint(nyaux,nx,ny,nsig,nt,ntr,tlat,lat_aux_grid,REGION_MASK_LAT_AUX, \
;       wedxdy_sig,wedxdy_sig,wedxdy_sig,wedxdy_sig@_FillValue,TMP1,TMP2,TMP3,TMP1@_FillValue)
; Use below for MOCoffline.new.so:
  MOC::moczonalint(nyaux,nx,ny,nsig,nt,ntr,tlat,lat_aux_grid,REGION_MASK_LAT_AUX, \
        wedxdy_sig,wedxdy_sig@_FillValue,TMP1,TMP1@_FillValue)
; wallClockElapseTime(TimeDate, "finished moczonalint",0)

; ==============================================================
; NOw, the meridional integral...
; =============================================================
  TMP4 = dim_cumsum(TMP1,1)
; TMP5 = dim_cumsum(TMP2,1)
; TMP6 = dim_cumsum(TMP3,1)

  MOCnew(:,0,:,:) =  (/ tofloat(TMP4(:,0,:,:)) /)
  MOCnew(:,1,:,:) =  (/ tofloat(TMP4(:,1,:,:)) /)
  delete(TMP1)
  delete(TMP2)
  delete(TMP3)
  delete(TMP4)

; ==============================================================
; compute MOC addition at Atlantic southern boundary
; ==============================================================
; TimeDateMOC = systemfunc("date")
  ; Find starting j-index for Atlantic region
  lat_aux_atl_start = nyaux
  do n=1,nyaux-1
    section = (tlat.ge.lat_aux_grid(n-1)).and.(tlat.lt.lat_aux_grid(n)).and.(REGION_MASK_LAT_AUX(1,:,:).eq.1)
    if (any(section).and.n.lt.lat_aux_atl_start) then
       lat_aux_atl_start = n-1
    end if
  end do

  ; compute VDXDZ in sigma-coord on TLONG,ULAT grid
  vedxdz_sig = where(ismissing(vedxdz_sig),0.,vedxdz_sig)

  WORK1 = vedxdz_sig
  WORK1 = (/ 0.5*(vedxdz_sig + eoshift(vedxdz_sig,1,-1)) /)

  atlmask = conform(WORK1, REGION_MASK_LAT_AUX(1,:,:), (/2,3/))
  j = lat_aux_atl_start
  TMP1 = where(atlmask(:,:,j,:).eq.1,WORK1(:,:,j-1,:),0.0)

  MOCs = new((/nt,nsig/),typeof(MOCz),MOCz@_FillValue)
  MOCs(:,:) = doubletofloat(dim_sum(TMP1))    ; zonal integral

  delete(WORK1)
  delete(TMP1)

  WORK4 = -MOCs
  MOCs(:,0:nsig-1:-1) = dim_cumsum(WORK4(:,0:nsig-1:-1),1)

  WORK5 = conform(MOCnew(:,1,0:nsig-1,:),MOCs,(/0,1/))
  MOCnew(:,1,0:nsig-1,:) = MOCnew(:,1,0:nsig-1,:)+WORK5
; wallClockElapseTime(TimeDate, "finished Atlantic southern boundary",0)

; ==============================================================
; Convert to Sverdrups
; ==============================================================
  MOCnew = MOCnew*1.0e-6

; ==============================================================
; Save to netcdf
; ==============================================================
  system("/bin/rm -f "+fout)
  outf = addfile(fout,"c")
  filedimdef( outf, "time", 0, True )
  outf->MOC=MOCnew
  outf->transport_regions=transport_regions

  wallClockElapseTime(TimeDate, "DONE: created "+fout,0)

end
