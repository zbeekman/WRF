load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

begin
  DOMAIN="ARCTAS_CARB_alt"
  STARTH="00"
  STARTD="27"
  STARTM="06"
  STARTY="2008"
  EMISDIR="/data/wrf/DATA/emissions/NEI2005/wrfemiss/"+DOMAIN+"/"+STARTY+STARTM+"/"
  inname=EMISDIR+"wrfchemi_d01_"+STARTY+"-"+STARTM+"-"+STARTD+"_"+STARTH+":00:00"
  INFILE  = addfile(inname+".nc", "r")

  att_names = getvaratts(INFILE)           ; get INFILE global attributes
  names = getfilevarnames(INFILE)       ; get INFILE variable names
  vartypes = getfilevartypes(INFILE,names)       ; get INFILE variable types


  do ihour = 0,23
     outname="wrfchemi_d01_" +STARTY+"-"+STARTM+"-"+STARTD+"_"+sprinti("%0.2i",ihour)+":00:00"
     system("rm "+outname+".nc")
     OUTFILE = addfile(outname+".nc","c")

     if(.not.all(ismissing(att_names))) then
       do i = 0,dimsizes(att_names)-1
   ;      print("copy_fileatts: global attributes->" + att_names(i))
         OUTFILE@$att_names(i)$ = INFILE@$att_names(i)$     ; copy INFILE global attributes
       end do
     end if
;
; loop over variables and copy to single hour file
;
     do i=0, dimsizes(names)-1
;       print("writing ncdf: i,name="+ i+ "   "+ names(i)+" @ hour = "+ihour)
       var = INFILE->$names(i)$
       dimvar = dimsizes(var)
       rank = dimsizes(dimvar)
       vartype = vartypes(i)
       if(rank .eq. 2) then
          var_plus=new((/1,dimvar(1)/),vartype)
          var_plus(0,:) = var(ihour,:)
          OUTFILE->$names(i)$ = var_plus
       else if(rank .eq. 3) then
          var_plus=new((/1,dimvar(1),dimvar(2)/),vartype)
          var_plus(0,:,:) = var(ihour,:,:)
          OUTFILE->$names(i)$ = var_plus
       else if(rank .eq. 4) then
          var_plus=new((/1,dimvar(1),dimvar(2),dimvar(3)/),vartype)
          var_plus(0,:,:,:) = var(ihour,:,:,:)
          OUTFILE->$names(i)$ = var_plus
       end if
       end if
       end if
       delete(var_plus)
       delete(var)
       delete(dimvar)
       delete(rank)
     end do
     delete(OUTFILE)
     print("Emissions dispersed to "+outname)
     system("mv "+outname+".nc "+outname)
  end do
end
