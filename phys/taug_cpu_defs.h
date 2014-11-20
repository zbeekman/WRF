#ifndef _CUDA
! changed to arguments for thread safety (could reduce this list a bit)
      real  :: pavel(:,:)
      real  :: wx1(:,:)
      real  :: wx2(:,:)
      real  :: wx3(:,:)
      real  :: wx4(:,:)
      real  :: coldry(:,:)
      integer  :: laytrop(:)
      integer  :: jp(:,:)
      integer  :: jt(:,:)
      integer  :: jt1(:,:)
      real  :: colh2o(:,:)
      real  :: colco2(:,:)
      real  :: colo3(:,:)
      real  :: coln2o(:,:)
      real  :: colco(:,:)
      real  :: colch4(:,:)
      real  :: colo2(:,:)
      real  :: colbrd(:,:)
      integer  :: indself(:,:)
      integer  :: indfor(:,:)
      real  :: selffac(:,:)
      real  :: selffrac(:,:)
      real  :: forfac(:,:)
      real  :: forfrac(:,:)
      integer  :: indminor(:,:)
      real  :: minorfrac(:,:)
      real  :: scaleminor(:,:)
      real  :: scaleminorn2(:,:)
      real  :: fac00(:,:), fac01(:,:), fac10(:,:), fac11(:,:)
      real  :: rat_h2oco2(:,:),rat_h2oco2_1(:,:), &
               rat_h2oo3(:,:),rat_h2oo3_1(:,:), &
               rat_h2on2o(:,:),rat_h2on2o_1(:,:), &
               rat_h2och4(:,:),rat_h2och4_1(:,:), &
               rat_n2oco2(:,:),rat_n2oco2_1(:,:), &
               rat_o3co2(:,:),rat_o3co2_1(:,:)
                                                      !    Dimensions: (ncol,nlayers)
      real  :: tauaa(:,:,:)
                                                      !    Dimensions: (ncol,nlayers,ngptlw)
     
      integer  :: nspad(:)
      integer  :: nspbd(:)
      real  :: oneminusd 
#endif
