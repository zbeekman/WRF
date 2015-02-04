load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

begin
 nx0=79
 ny0=79
 nz0=30
 nx=79
 ny=79
 nz=30

 SCALE = 0.5; needs to be modified depending on sensitivity window (0.5 - 3hr, 1.0 - 6hr, etc..)

 adj_points0=(/28,64,3/)

; adj_points0=(/(/28,64,1/), \
;               (/35,53,1/)/)
;               (/39,37,5/)/)

; adj_points0=(/(/28,64,3/), \
;               (/35,53,4/)/)

; adj_points0=(/(/28,64,3/), \
;               (/35,53,4/), \
;               (/39,37,5/)/)

 adj_dim=dimsizes(adj_points0)
 ndims=dimsizes(adj_dim)
 adj_points      = adj_points0

 if(ndims .gt. 1) then
    adj_points(:,0) = adj_points0(:,0) * nx/nx0
    adj_points(:,1) = adj_points0(:,1) * ny/ny0
    adj_points(:,2) = adj_points0(:,2) * nz/nz0
    ninverse=adj_dim(0)
 else
    adj_points(0) = adj_points0(0) * nx/nx0
    adj_points(1) = adj_points0(1) * ny/ny0
    adj_points(2) = adj_points0(2) * nz/nz0
    ninverse=1
 end if

; TRANSFORM0= (/ 0,   0, -2/)

; TRANSFORM0= (/ -5,   2,  -1/)

; TRANSFORM0= (/(/ 0,   0,  -2/), \
;               (/ 0,   0,  -1/)/)

; TRANSFORM0= (/(/ -5,   2,  -1/), \
;               (/ -5,   0,  1/)/)

; TRANSFORM0= (/(/ -8,   0,  0/), \
;               (/  0,  -2,  0/), \
;               (/ -5,  -2,  0/), \
;               (/ -5,   2,  0/), \
;               (/ -5,   0,  1/)/)

 TRANSFORM0= (/(/ -8,   0,  -3/), \
               (/  0,   5,  -5/), \
               (/ -5,   5,  -4/)/)

; TRANSFORM0= (/(/ -8,   0,  -2/), \
;               (/ -5,   2,  -1/), \
;               (/ -5,   0,  1/)/)

; TRANSFORM0= (/(/ -8,   0,  -2/), \
;               (/ -5,   2,  -1/), \
;               (/ 0,   0,   -2/)/)

; TRANSFORM0= (/(/ -8,   0,  -2/), \
;               (/ -5,   2,  -1/), \
;               (/ -5,   0,   1/), \
;               (/  0,   0,  -2/), \
;               (/  0,   0,  -1/)/)


; TRANSFORM0= (/(/ -8,   0,  -2/), \
;               (/  0,  -2,  -2/), \
;               (/ -5,  -2,  -2/), \
;               (/ -5,   2,  -1/), \
;               (/ -5,   0,   1/)/)

; TRANSFORM0= (/(/ -8,   0,  -2/), \
;               (/  0,  -2,  -2/), \
;               (/ -5,  -2,  -2/), \
;               (/ -5,   2,  -1/), \
;               (/  0,  -2,   3/), \
;               (/ -5,  -2,   3/), \
;               (/ -5,   0,   1/)/)
 
 TRANSFORM=TRANSFORM0
 TR_dim=dimsizes(TRANSFORM0)

 if(dimsizes(TR_dim) .gt. 1)then
    TRANSFORM(:,0) = doubletoint(TRANSFORM0(:,0) * nx/nx0*SCALE)
    TRANSFORM(:,1) = doubletoint(TRANSFORM0(:,1) * ny/ny0*SCALE)
;    TRANSFORM(:,2) = doubletoint(TRANSFORM0(:,2) * nz/nz0*SCALE)
    firatio=TR_dim(0)
 else
    delete(TRANSFORM)
    TRANSFORM=new((/2,3/),integer)
    TRANSFORM=0
    TRANSFORM(0,0) = doubletoint(TRANSFORM0(0) * nx/nx0*SCALE)
    TRANSFORM(0,1) = doubletoint(TRANSFORM0(1) * ny/ny0*SCALE)
    TRANSFORM(0,2) = TRANSFORM0(2)
;    TRANSFORM(2) = doubletoint(TRANSFORM0(:,2) * nz/nz0*SCALE)
    firatio=1
 end if

 tl_points=new((/firatio*ninverse,3/),integer)
 do na = 0,ninverse-1
   do nt = 0,firatio-1
      if(ndims .gt. 1) then
         tl_points(nt+firatio*na,:)=adj_points(na,:)+TRANSFORM(nt,:)
      else
         tl_points(nt+firatio*na,:)=adj_points(:)+TRANSFORM(nt,:)
      end if
      do index=0,2
         if(tl_points(nt+firatio*na,index) .le. 0 ) then
            tl_points(nt+firatio*na,index) = 1
         end if
      end do
   end do
 end do

 if(ndims .gt. 1) then
    asciiwrite("locations_f",sprinti("%3.0i",adj_points(:,0)) +"    " +\
                              sprinti("%3.0i",adj_points(:,1)) +"    " +\
                              sprinti("%3.0i",adj_points(:,2))) 
 else
    asciiwrite("locations_f",sprinti("%3.0i",adj_points(0)) +"    " +\
                              sprinti("%3.0i",adj_points(1)) +"    " +\
                              sprinti("%3.0i",adj_points(2))) 

 end if
 asciiwrite("locations_i",sprinti("%3.0i",tl_points(:,0)) +"    " +\
                           sprinti("%3.0i",tl_points(:,1)) +"    " +\
                           sprinti("%3.0i",tl_points(:,2))) 

end
