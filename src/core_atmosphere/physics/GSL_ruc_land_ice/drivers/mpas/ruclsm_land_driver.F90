! Author: Tanya Smirnova (Tanya.Smirnova@noaa.gov)
! 15 march 2024 - RUC LSM is updated to the version used in RRFS version 1.
! Most changes are in the snow model.
!
#define lsmruc_dbg_lvl 3000
!wrf:model_layer:physics
!
module ruclsm_land_driver

! notes for perturbations of soil properties (judith berner)
! perturbations are applied in subroutine soilprob to array hydro;
! soilprop is called from subroutine sfctmp which is called from subroutine lsmruc;
! subroutine lsmruc had two new 3d fields: pattern_spp_lsm (in) and field_sf(inout);
!    their vertical dimension is number of atmospheric levels (kms:kme) - (suboptimal, but easiest hack)
!    field_sf is used to pass perturbed fields of hydrop up to model (and output) driver;
! in argument list to sfctmp the arrays are passed as pattern_spp_lsm(i,1:nzs,j), and exist henceforth as
! column arrays;
! in the subroutines below sfctmp (snow and snowsoil) the fields are called rstochcol,fieldcol_sf
! to reflect their dimension rstochcol (1:nzs)

#if defined(mpas)
use mpas_atmphys_constants, g => gravity,rhowater=>rho_w, piconst => pii
use mpas_atmphys_utilities, only: physics_error_fatal,physics_message
use mpas_log, only: mpas_log_write
#define em_core 1
#define fatal_error(m) call physics_error_fatal( m )
#define wrf_at_debug_level(lsmruc_dbg_lvl) .false.
#else
  use module_model_constants
  use module_wrf_error
#define fatal_error(m) call wrf_error_fatal( m )
#endif


contains

!-----------------------------------------------------------------
    subroutine ruc_land_driver(spp_lsm,                          &
#if (em_core==1)
                   pattern_spp_lsm,field_sf,                     &
#endif
                   dt,ktau,nsl,                                  &
#if (em_core==1)
                   lakemodel,lakemask,                           &
                   graupelncv,snowncv,rainncv,                   &
#endif
                   zs,rainbl,snow,snowh,snowc,frzfrac,frpcpn,    &
                   rhosnf,precipfr,                              & ! pass it out to module_diagnostics
                   z3d,p8w,t3d,qv3d,qc3d,rho3d,                  & !p8w in [pa]
                   glw,gsw,emiss,chklowq, chs,                   &
                   flqc,flhc,mavail,canwat,vegfra,alb,znt,       &
                   z0,snoalb,albbck,lai,                         &  !new
                   mminlu, landusef, nlcat, mosaic_lu,           &
                   mosaic_soil, soilctop, nscat,                 &  !new
                   qsfc,qsg,qvg,qcg,dew,soilt1,tsnav,            &
                   tbot,ivgtyp,isltyp,xland,                     &
                   iswater,isice,xice,xice_threshold,            &
                   cp,rovcp,g0,lv,stbolt,                        &
                   soilmois,sh2o,smavail,smmax,                  &
                   tso,soilt,hfx,qfx,lh,                         &
                   sfcrunoff,udrunoff,acrunoff,sfcexc,           &
                   sfcevp,grdflx,snowfallac,acsnow,snom,         &
                   smfr3d,keepfr3dflag,                          &
                   myjpbl,shdmin,shdmax,rdlai2d,                 &
                   globalcells,                                  &
                   ids,ide, jds,jde, kds,kde,                    &
                   ims,ime, jms,jme, kms,kme,                    &
                   its,ite, jts,jte, kts,kte                     )
!-----------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------
!
! the RUC LSM model is described in:
!  Smirnova, T.G., J.M. Brown, and S.G. Benjamin, 1997:
!     performance of different soil model configurations in simulating
!     ground surface temperature and surface fluxes.
!     mon. wea. rev. 125, 1870-1884.
!  Smirnova, T.G., J.M. Brown, and D. Kim, 2000: parameterization of
!     cold-season processes in the maps land-surface scheme.
!     j. geophys. res. 105, 4077-4086.
!-----------------------------------------------------------------
!-- dt            time step (second)
!   ktau - number of time step
!   nsl  - number of soil layers
!   nzs  - number of levels in soil
!   zs   - depth of soil levels (m)
!-- rainbl    - accumulated rain in [mm] between the pbl calls
!-- rainncv         one time step grid scale precipitation (mm/step)
!-- snow - snow water equivalent [mm]
!-- frazfrac - fraction of frozen precipitation
!-- precipfr (mm) - time step frozen precipitation
!-- snowc       flag indicating snow coverage (1 for snow cover)
!-- z3d         heights (m)
!-- p8w         3d pressure (pa)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- qc3d - 3d cloud water mixing ratio (kg/kg)
!-- rho3d - 3d air density (kg/m^3)
!-- glw         downward long wave flux at ground surface (w/m^2)
!-- gsw         absorbed short wave flux at ground surface (w/m^2)
!-- emiss       surface emissivity (between 0 and 1)
!   flqc - surface exchange coefficient for moisture (kg/m^2/s)
!   flhc - surface exchange coefficient for heat [w/m^2/s/degreek]
!   sfcexc - surface exchange coefficient for heat [m/s]
!   canwat - canopy moisture content (mm)
!   vegfra - vegetation fraction (between 0 and 100)
!   alb - surface albedo (between 0 and 1)
!   snoalb - maximum snow albedo (between 0 and 1)
!   albbck - snow-free albedo (between 0 and 1)
!   znt - roughness length [m]
!-- tbot        soil temperature at lower boundary (k)
!   ivgtyp - usgs vegetation type (24 classes)
!   isltyp - stasgo soil type (16 classes)
!-- xland       land mask (1 for land, 2 for water)
!-- cp          heat capacity at constant pressure for dry air (j/kg/k)
!-- g0          acceleration due to gravity (m/s^2)
!-- lv          latent heat of melting (j/kg)
!-- stbolt      stefan-boltzmann constant (w/m^2/k^4)
!   soilmois - soil moisture content (volumetric fraction)
!   tso - soil temp (k)
!-- soilt       surface temperature (k)
!-- hfx         upward heat flux at the surface (w/m^2)
!-- qfx         upward moisture flux at the surface (kg/m^2/s)
!-- lh          upward latent heat flux (w/m^2)
!   sfcrunoff - ground surface runoff [mm]
!   udrunoff - underground runoff [mm]
!   acrunoff - run-total surface runoff [mm]
!   sfcevp - total evaporation in [kg/m^2]
!   grdflx - soil heat flux (w/m^2: negative, if downward from surface)
!   snowfallac - run-total snowfall accumulation [m]
!   acsnow - run-toral swe of snowfall [mm]
!-- chklowq - is either 0 or 1 (so far set equal to 1).
!--           used only in myjpbl.
!-- tice - sea ice temperture (c)
!-- rhosice - sea ice density (kg m^-3)
!-- capice - sea ice volumetric heat capacity (j/m^3/k)
!-- thdifice - sea ice thermal diffusivity (m^2/s)
!--
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-------------------------------------------------------------------------
!   integer,     parameter            ::     nzss=5
!   integer,     parameter            ::     nddzs=2*(nzss-2)

   integer,     parameter            ::     nvegclas=24+3

   real,       intent(in   )    ::     dt
   logical,    intent(in   )    ::     myjpbl,frpcpn
   integer,    intent(in   )    ::     spp_lsm
   integer,    intent(in   )    ::     nlcat, nscat, mosaic_lu, mosaic_soil
   integer,    intent(in   )    ::     ktau, nsl, isice, iswater, &
                                       ims,ime, jms,jme, kms,kme, &
                                       ids,ide, jds,jde, kds,kde, &
                                       its,ite, jts,jte, kts,kte
   integer,   intent(in   )     ::     globalcells(ims:ime)

#if (em_core==1)
   real,    dimension( ims:ime, kms:kme, jms:jme ),optional::    pattern_spp_lsm
   real,    dimension( ims:ime, kms:kme, jms:jme ),optional::    field_sf
#endif
   real,    dimension( ims:ime, 1  :nsl, jms:jme )         ::    field_sf_loc

   real,    dimension( ims:ime, kms:kme, jms:jme )            , &
            intent(in   )    ::                           qv3d, &
                                                          qc3d, &
                                                           p8w, &
                                                         rho3d, &
                                                           t3d, &
                                                           z3d

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                       rainbl, &
                                                            glw, &
                                                            gsw, &
                                                         albbck, &
                                                           chs , &
                                                           xice, &
                                                          xland, &
                                                         vegfra, &
                                                           tbot

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(inout   )    ::         flqc,flhc


#if (em_core==1)
   real,       optional, dimension( ims:ime , jms:jme ),         &
               intent(in   )    ::                   graupelncv, &
                                                        snowncv, &
                                                        rainncv
   real,       dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                     lakemask
   integer,    intent(in   )    ::                    lakemodel
#endif

   real, dimension( ims:ime , jms:jme ), intent(in )::   shdmax
   real, dimension( ims:ime , jms:jme ), intent(in )::   shdmin
   logical, intent(in) :: rdlai2d

   real,       dimension( 1:nsl), intent(in   )      ::      zs

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(inout)    ::                               &
                                                           snow, &
                                                          snowh, &
                                                          snowc, &
                                                         canwat, & ! new
                                                         snoalb, &
                                                            alb, &
                                                          emiss, &
                                                            lai, &
                                                         mavail, &
                                                         sfcexc, &
                                                            z0 , &
                                                            znt

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                               &
                                                        frzfrac

   integer,    dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                       ivgtyp, &
                                                         isltyp
   character(len=*), intent(in   )    ::                 mminlu
   real,     dimension( ims:ime , 1:nlcat, jms:jme ), intent(in):: landusef
   real,     dimension( ims:ime , 1:nscat, jms:jme ), intent(in):: soilctop

   real, intent(in   )          ::         cp,rovcp,g0,lv,stbolt,xice_threshold

   real,       dimension( ims:ime , 1:nsl, jms:jme )           , &
               intent(inout)    ::                 soilmois,sh2o,tso

   real,       dimension( ims:ime, jms:jme )                   , &
               intent(inout)    ::                        soilt, &
                                                            hfx, &
                                                            qfx, &
                                                             lh, &
                                                         sfcevp, &
                                                      sfcrunoff, &
                                                       udrunoff, &
                                                       acrunoff, &
                                                         grdflx, &
                                                         acsnow, &
                                                           snom, &
                                                            qvg, &
                                                            qcg, &
                                                            dew, &
                                                           qsfc, &
                                                            qsg, &
                                                        chklowq, &
                                                         soilt1, &
                                                          tsnav

   real,       dimension( ims:ime, jms:jme )                   , &
               intent(inout)    ::                      smavail, &
                                                          smmax

   real,       dimension( its:ite, jts:jte )    ::               &
                                                             pc, &
                                                        runoff1, &
                                                        runoff2, &
                                                         emissl, &
                                                           zntl, &
                                                        lmavail, &
                                                          smelt, &
                                                           snoh, &
                                                          snflx, &
                                                           edir, &
                                                             ec, &
                                                            ett, &
                                                         sublim, &
                                                           sflx, &
                                                            smf, &
                                                          evapl, &
                                                          prcpl, &
                                                         seaice, &
                                                        infiltr
!-- energy and water budget variables:
   real,       dimension( its:ite, jts:jte )    ::               &
                                                         budget, &
                                                       acbudget, &
                                                    waterbudget, &
                                                  acwaterbudget, &
                                                       smtotold, &
                                                        snowold, &
                                                      canwatold


   real,       dimension( ims:ime, 1:nsl, jms:jme)               &
                                             ::    keepfr3dflag, &
                                                         smfr3d

   real,       dimension( ims:ime, jms:jme ), intent(out)     :: &
                                                         rhosnf, & ! rho of snowfall
                                                       precipfr, & ! time-step frozen precip
                                                     snowfallac
!-- soil/snow properties
   real                                                          &
                             ::                           rhocs, &
                                                       rhonewsn, &
                                                          rhosn, &
                                                      rhosnfall, &
                                                           bclh, &
                                                            dqm, &
                                                           ksat, &
                                                           psis, &
                                                           qmin, &
                                                          qwrtz, &
                                                            ref, &
                                                           wilt, &
                                                        canwatr, &
                                                       snowfrac, &
                                                          snhei, &
                                                           snwe

   real                                      ::              cn, &
                                                         sat,cw, &
                                                           c1sn, &
                                                           c2sn, &
                                                         kqwrtz, &
                                                           kice, &
                                                            kwt


   real,     dimension(1:nsl)                ::          zsmain, &
                                                         zshalf, &
                                                         dtdzs2

   real,     dimension(1:2*(nsl-2))          ::           dtdzs

   real,     dimension(1:5001)               ::             tbq


   real,     dimension( 1:nsl )              ::         soilm1d, &
                                                          tso1d, &
                                                        soilice, &
                                                        soiliqw, &
                                                       smfrkeep

   real,     dimension( 1:nsl )              ::          keepfr

   real,     dimension( 1:nlcat )            ::          lufrac
   real,     dimension( 1:nscat )            ::          soilfrac

   real                           ::                        rsm, &
                                                      snweprint, &
                                                     snheiprint

   real                           ::                     prcpms, &
                                                        newsnms, &
                                                      prcpncliq, &
                                                       prcpncfr, &
                                                      prcpculiq, &
                                                       prcpcufr, &
                                                           patm, &
                                                          patmb, &
                                                           tabs, &
                                                          qvatm, &
                                                          qcatm, &
                                                          q2sat, &
                                                         conflx, &
                                                            rho, &
                                                           qkms, &
                                                           tkms, &
                                                        snowrat, &
                                                       grauprat, &
                                                       graupamt, &
                                                         icerat, &
                                                          curat, &
                                                       infiltrp
   real      ::  cq,r61,r273,arp,brp,x,evs,eis
   real      ::  cropfr, cropsm, newsm, factor

   real      ::  meltfactor, ac,as, wb
   integer   ::  nroot
   integer   ::  iland,isoil,iforest

   integer   ::  i,j,k,nzs,nzs1,nddzs
   integer   ::  k1,l,k2,kp,km
   integer   ::  globalcellid
   character (len=132) :: message
   logical   :: print_flag = .false.
   real,dimension(ims:ime,1:nsl,jms:jme) :: rstoch
   real,dimension(1:nsl)::rstoch_temp,field_sf_temp
   real,dimension(ims:ime,jms:jme)::emisso,vegfrao,albo,snoalbo
   real,dimension(its:ite,jts:jte)::emisslo

!-----------------------------------------------------------------

       call ruc_land( spp_lsm, lakemodel, lakemask,  &
                    rhosnf, precipfr, mosaic_lu, &
                    qsg, qvg, dew,       &
                    soilt1, tsnav, acrunoff,  &
                    snowfallac, keepfr3dflag, stbolt,      &
                    graupelncv, snowncv, rainncv,   &
                    qcg, flqc, flhc,      &
                    landusef, nlcat, soilctop,     &
                    nscat, iswater, chklowq,   &
                    sfcexc, sfcevp, myjpbl,      &
                    dt, ktau, nsl,   &
                    rainbl, snow, snowh,     &
                    snowc, frzfrac, frpcpn,      &
                    z3d, p8w, t3d,         &
                    qv3d, qc3d, rho3d,       &
                    glw, gsw, emiss, &
                    chs, mavail, canwat,    &
                    vegfra, alb, znt,       &
                    z0, snoalb, albbck,&
                    lai, mminlu, qsfc,      &
                    tbot, ivgtyp, isltyp,    &
                    xland, isice, xice,      &
                    xice_threshold, cp, rovcp,         &
                    g0, lv, soilmois,     &
                    sh2o, smavail, smmax,    &
                    tso, soilt, hfx,       &
                    qfx, lh, sfcrunoff, &
                    udrunoff, grdflx, acsnow,    &
                    snom, shdmin, rdlai2d,     &
                    mosaic_soil, zs, smfr3d,    &
                    shdmax, globalCells,                           &
                    ids, ide, jds, jde, kds, kde,  &
                    ims, ime, jms, jme, kms, kme,  &
                    its, ite, jts, jte, kts, kte&
                          )

   end subroutine ruc_land_driver

end module ruclsm_land_driver
