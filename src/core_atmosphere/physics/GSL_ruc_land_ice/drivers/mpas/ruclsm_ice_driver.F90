! Author: Tanya Smirnova (Tanya.Smirnova@noaa.gov)
! 15 march 2024 - RUC LSM is updated to the version used in RRFS version 1.
! Most changes are in the snow model.
!
#define iceruc_dbg_lvl 3000
!wrf:model_layer:physics
!
module ruclsm_ice_driver

#if defined(mpas)
use mpas_atmphys_constants, g => gravity,rhowater=>rho_w, piconst => pii
use mpas_atmphys_utilities, only: physics_error_fatal,physics_message
use mpas_log, only: mpas_log_write
#define em_core 1
#define fatal_error(m) call physics_error_fatal( m )
#define wrf_at_debug_level(iceruc_dbg_lvl) .false.
#else
  use module_model_constants
  use module_wrf_error
#define fatal_error(m) call wrf_error_fatal( m )
#endif


contains

!-----------------------------------------------------------------
   subroutine ruc_ice_driver(spp_lsm,                       &
!-- in
#if (em_core==1)
              pattern_spp_lsm,field_sf,                     &
#endif
              dt,ktau,nsl,                                  &
#if (em_core==1)
              graupelncv,snowncv,rainncv,                   &
#endif
              zs,rainbl,frzfrac,frpcpn,rhosnf,precipfr,     &
              z3d,p8w,t3d,qv3d,qc3d,rho3d,                  & 
              glw,gsw,emiss,chklowq, chs,                   &
              flqc,flhc,alb,znt,z0,snoalb,albbck, mminlu,   &
              tbot,ivgtyp,isltyp,xland,myjpbl,              &
              iswater,isice,xice,xice_threshold,            &
!-- constants                   
              cp,rovcp,g0,lv,stbolt,                        &
!-- in/out                   
              tso,soilt,soilt1,tsnav,snow,snowh,snowc,      &
              qsfc,qsg,qvg,qcg,                             &
!-- out                   
              hfx,qfx,lh,dew,sfcrunoff,acrunoff,            &
              sfcevp,grdflx,snowfallac,acsnow,snom,         &
              globalcells,                                  &
              ids,ide, jds,jde, kds,kde,                    &
              ims,ime, jms,jme, kms,kme,                    &
              its,ite, jts,jte, kts,kte                     )
!-----------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------
!
! the RUC ice model is a part of RUC LSM described in:
!  Smirnova, T.G., J.M. Brown, and S.G. Benjamin, 1997:
!     performance of different soil model configurations in simulating
!     ground surface temperature and surface fluxes.
!     mon. wea. rev. 125, 1870-1884, 
!     https://doi.org/10.1175/1520-0493(1997)125%3C1870:PODSMC%3E2.0.CO;2 
!  Smirnova, T.G., J.M. Brown, and D. Kim, 2000: parameterization of
!     cold-season processes in the maps land-surface scheme.
!     j. geophys. res. 105, 4077-4086, https://doi.org/10.1029/1999JD901047
!  Smirnova, T. G., J. M. Brown, S. G. Benjamin, and J. S. Kenyon, 2016: 
!     Modifications to the Rapid Update Cycle land surface model (RUC LSM) 
!     available in the weather Research and forecasting model. 
!     Mon. Wea. Rev., 144, 1851–1865, https://doi.org/10.1175/MWR-D-15-0198.1  
!-----------------------------------------------------------------
!-- dt         time step (second)
!   ktau       number of time step
!   nsl        number of soil layers
!   nzs        number of levels in soil
!   zs         depth of soil levels (m)
!-- rainbl     accumulated rain in [mm] between the pbl calls
!-- rainncv    one time step grid scale precipitation (mm/step)
!-- snow       snow water equivalent [mm]
!-- frazfrac   fraction of frozen precipitation
!-- precipfr (mm) time step frozen precipitation
!-- snowc      flag indicating snow coverage (1 for snow cover)
!-- z3d        heights (m)
!-- p8w        3d pressure (pa)
!-- t3d        temperature (k)
!-- qv3d       3d water vapor mixing ratio (kg/kg)
!-- qc3d       3d cloud water mixing ratio (kg/kg)
!-- rho3d      3d air density (kg/m^3)
!-- glw        downward long wave flux at ground surface (w/m^2)
!-- gsw        absorbed short wave flux at ground surface (w/m^2)
!-- emiss      surface emissivity (between 0 and 1)
!   flqc       surface exchange coefficient for moisture (kg/m^2/s)
!   flhc       surface exchange coefficient for heat [w/m^2/s/degreek]
!   alb        surface albedo (between 0 and 1)
!   snoalb     maximum snow albedo (between 0 and 1)
!   albbck     snow-free albedo (between 0 and 1)
!   znt        roughness length [m]
!-- tbot       soil temperature at lower boundary (k)
!   ivgtyp     vegetation type
!   isltyp     soil type (16 classes)
!-- xland      land mask (1 for land, 2 for water)
!-- cp         heat capacity at constant pressure for dry air (j/kg/k)
!-- g0         acceleration due to gravity (m/s^2)
!-- lv         latent heat of melting (j/kg)
!-- stbolt     stefan-boltzmann constant (w/m^2/k^4)
!   tso        3d ice temp (k)
!-- soilt      ice surface temperature (k)
!-- hfx        upward heat flux at the surface (w/m^2)
!-- qfx        upward moisture flux at the surface (kg/m^2/s)
!-- lh         upward latent heat flux (w/m^2)
!   sfcrunoff  surface runoff [mm]
!   acrunoff   run-total surface runoff [mm]
!   sfcevp     total evaporation in [kg/m^2]
!   grdflx     soil heat flux (w/m^2: negative, if downward from surface)
!   snowfallac run-total snowfall accumulation [m]
!   acsnow     run-toral swe of snowfall [mm]
!-- chklowq    is either 0 or 1 (so far set equal to 1).
!--            used only in myjpbl.
!-- tice       sea ice temperture (c)
!-- rhosice    sea ice density (kg m^-3)
!-- capice     sea ice volumetric heat capacity (j/m^3/k)
!-- thdifice   sea ice thermal diffusivity (m^2/s)
!--
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-------------------------------------------------------------------------

   integer,     parameter       ::     nvegclas=24+3

   real,       intent(in   )    ::     dt
   logical,    intent(in   )    ::     myjpbl,frpcpn
   integer,    intent(in   )    ::     spp_lsm
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

   real,    dimension( ims:ime, kms:kme, jms:jme )             , &
            intent(in   )    ::                            qv3d, &
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
                                                           tbot

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                    flqc,flhc


#if (em_core==1)
   real,       optional, dimension( ims:ime , jms:jme ),         &
               intent(in   )    ::                   graupelncv, &
                                                        snowncv, &
                                                        rainncv
#endif


   real,       dimension( 1:nsl), intent(in   )      ::      zs

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(inout)    ::                      chklowq, &
                                                         snoalb, &
                                                           snow, &
                                                          snowh, &
                                                          snowc, &
                                                            alb, &
                                                            znt, &
                                                             z0, &
                                                          emiss

   real,       dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                      frzfrac

   integer,    dimension( ims:ime , jms:jme ),                   &
               intent(in   )    ::                       ivgtyp, &
                                                         isltyp
   character(len=*), intent(in   )    ::                 mminlu

   real, intent(in   )          ::         cp,rovcp,g0,lv,stbolt,xice_threshold

   real,       dimension( ims:ime , 1:nsl, jms:jme )           , &
               intent(inout)    ::                       tso

   real,       dimension( ims:ime, jms:jme )                   , &
               intent(inout)    ::                        soilt, &
                                                            qvg, &
                                                            qcg, &
                                                            dew, &
                                                           qsfc, &
                                                            qsg, &
                                                         soilt1, &
                                                          tsnav
   real,       dimension( ims:ime, jms:jme )                   , &
               intent(inout)    ::                          hfx, &
                                                            qfx, &
                                                             lh, &
                                                         sfcevp, &
                                                      sfcrunoff, &
                                                       acrunoff, &
                                                         grdflx, &
                                                         acsnow, &
                                                           snom

   real,       dimension( ims:ime, jms:jme ), intent(out)     :: &
                                                         rhosnf, & ! rho of snowfall
                                                       precipfr, & ! time-step frozen precip
                                                     snowfallac

!-- local variables
   real,       dimension( its:ite, jts:jte )    ::               &
                                                        runoff1, &
                                                        runoff2, &
                                                         emissl, &
                                                         budget, &
                                                           zntl, &
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
                                                       snowfrac, &
                                                          snhei, &
                                                           snwe

   real                                      ::              cn, &
                                                         sat,cw, &
                                                           c1sn, &
                                                           c2sn


   real,     dimension(1:nsl)                ::          zsmain, &
                                                         zshalf, &
                                                         dtdzs2

   real,     dimension(1:nsl)                ::          soilice,&
                                                         soiliqw           

   real,     dimension(1:2*(nsl-2))          ::           dtdzs

   real,     dimension(1:5001)               ::             tbq


   real,     dimension( 1:nsl )              ::           tso1d

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
   real      ::  qsn
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

          call ruc_ice( spp_lsm, iswater  , chklowq   , &
                    rhosnf     , precipfr , chs       , &
                    qsg        , qvg      , dew       , &
                    soilt1     , tsnav    , acrunoff  , &
                    snowfallac , stbolt   , qsfc      , &
                    graupelncv , snowncv  , rainncv   , &
                    qcg        , flqc     , flhc      , &
                    dt         , ktau     , nsl       , &
                    rainbl     , snow     , snowh     , &
                    snowc      , frzfrac  , frpcpn    , &
                    z3d        , p8w      , t3d       , &
                    qv3d       , qc3d     , rho3d     , &
                    glw        , gsw      , emiss     , &
                    alb        , znt      , zs        , &
                    z0         , snoalb   , albbck    , &
                    tbot       , ivgtyp   , isltyp    , &
                    xland      , isice    , xice      , &
                    xice_threshold, cp    , rovcp     , &
                    g0         , lv       , snom      , &
                    tso        , soilt    , hfx       , &
                    qfx        , lh       , sfcrunoff , &
                    grdflx     , acsnow   , sfcevp    , &
                    myjpbl     , mminlu   , globalCells,&
                    ids, ide, jds, jde, kds, kde,       &
                    ims, ime, jms, jme, kms, kme,       &
                    its, ite, jts, jte, kts, kte        &
                          )
       call sfcdiags_ruclsm( &
                    hfx    , qfx  , tsk , qsfc, chs ,   & 
                    chs2   , cqs2 , t2  , th2 , q2  ,   &
                    psfc2d , t3d  , qv3d, dz  , cp  , R_d, &
                    rovcp  , rho3d, p3d , snow, cqs ,   &
                    ids, ide , jds , jde , kds , kde ,  &
                    ims, ime , jms , jme , kms , kme ,  &
                    )

!-----------------------------------------------------------------
   end subroutine ruc_ice_driver
!------------------------------------------------------------------------

end module ruclsm_ice_driver
