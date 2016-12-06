module mod_nsi_element_assembly

  !------------------------------------------------------------------------
  !> @addtogroup NastinMatrixAssembly
  !> @{
  !> @file    nsi_elmmat.f90
  !> @author  Guillaume Houzeaux
  !> @brief   Navier-Stokes system element assembly and other element 
  !>          calculations
  !> @details Compute element matrices
  !>
  !>          \verbatim
  !> 
  !>          Without enrichement
  !>          -------------------
  !>            +-        +  +- -+     +-  -+
  !>            | Auu Aup |  | u |     | bu |
  !>            |         |  |   |  =  |    |
  !>            | Apu App |  | p |     | bp |
  !>            +-       -+  +- -+     +-  -+
  !>
  !>          With enrichement
  !>          ----------------
  !>            +-            +  +- -+     +-  -+
  !>            | Auu Aup Auq |  | u |     | bu |
  !>            |             |  |   |     |    |
  !>            | Apu App Apq |  | p |  =  | bp |
  !>            |             |  |   |     |    |
  !>            | Aqu Aqp Aqq |  | q |     | bq | 
  !>            +-           -+  +- -+     +-  -+
  !>
  !>            q = Aqq^-1 ( bq - Aqu u - Aqp p )
  !>
  !>            Auu <= Auu - Auq Aqq^-1 Aqu
  !>            Aup <= Aup - Auq Aqq^-1 Aqp
  !>            bu  <= bu  - Auq Aqq^-1 bq
  !>
  !>            Apu <= Apu - Apq Aqq^-1 Aqu 
  !>            App <= App - Apq Aqq^-1 Aqp
  !>            bp  <= bp  - Apq Aqq^-1 bq
  !>
  !>          \endverbatim
  !>
  !> @{
  !------------------------------------------------------------------------
   
  use def_parameters
  implicit none

  private 

  public :: nsi_element_assembly_asgs_oss     
  public :: nsi_element_assembly_split_oss     

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   ASGS and OSS
  !> @details Assembly of Navier Stokes equations using ASGS or OSS
  !>          Variational Multiscale Model.
  !>          Includes subgriscale tracking in time and convection
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_element_assembly_asgs_oss(             &
       pnode,pgaus,gpden,gpvis,gppor,gpgvi,gpsp1,       &
       gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
       gpadv,gpvep,gprhs,gprhc,rmomu,rcont,elpre,elbub, &
       elauu,elaup,elapp,elapu,elrbu,elrbp,rmom2,gpst1, &
       gpgve,gprh2,gprhs_sgs,elvel,ellum,dtinv_loc,     &
       pbubl,gpsha_bub,gpcar_bub,gppre,elauq,elapq,     &
       elaqu,elaqp,elaqq,elrbq,&
       ! Original global variables
       kfl_lumped,&
       mnode,ntens,&
       kfl_duatss,&
       fact_duatss,&
       kfl_stabi_nsi,&
       fvins_nsi,fcons_nsi,bemol_nsi,kfl_regim_nsi,&
       fvela_nsi,kfl_rmom2_nsi,kfl_press_nsi,&
       kfl_p1ve2_nsi,kfl_linea_nsi,pabdf_nsi,&
       kfl_confi_nsi,nbdfp_nsi,kfl_sgsti_nsi,&
       kfl_nota1_nsi,kfl_limit_nsi,kfl_penal_nsi,&
       penal_nsi,&
       kfl_bubbl_nsi)

    integer(ip), intent(in)    :: kfl_lumped
    integer(ip), intent(in)    :: mnode
    integer(ip), intent(in)    :: ntens
    integer(ip), intent(in)    :: kfl_duatss
    integer(ip), intent(in)    :: fact_duatss
    integer(ip), intent(in)    :: kfl_stabi_nsi
    real(rp),    intent(in)    :: fvins_nsi
    real(rp),    intent(in)    :: fcons_nsi
    real(rp),    intent(in)    :: bemol_nsi
    integer(ip), intent(in)    :: kfl_regim_nsi
    real(rp),    intent(in)    :: fvela_nsi(3)
    integer(ip), intent(in)    :: kfl_rmom2_nsi
    integer(ip), intent(in)    :: kfl_press_nsi
    integer(ip), intent(in)    :: kfl_p1ve2_nsi
    integer(ip), intent(in)    :: kfl_linea_nsi
    real(rp),    intent(in)    :: pabdf_nsi(*)
    integer(ip), intent(in)    :: kfl_confi_nsi
    integer(ip), intent(in)    :: nbdfp_nsi
    integer(ip), intent(in)    :: kfl_sgsti_nsi
    integer(ip), intent(in)    :: kfl_nota1_nsi
    integer(ip), intent(in)    :: kfl_limit_nsi
    integer(ip), intent(in)    :: kfl_penal_nsi
    real(rp),    intent(in)    :: penal_nsi
    integer(ip), intent(in)    :: kfl_bubbl_nsi

    integer(ip), intent(in)    :: pnode,pgaus
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpgvi(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gptt1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gptt2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gplap(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)
    real(rp),    intent(in)    :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhs_sgs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gprhc(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gprh2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: rmomu(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: rcont(VECTOR_SIZE,ndime,pnode,pgaus)
    real(rp),    intent(in)    :: elpre(VECTOR_SIZE,pnode,*)
    real(rp),    intent(in)    :: elbub(VECTOR_SIZE)
    ! Element matrices
    real(rp),    intent(out)   :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(out)   :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(out)   :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(out)   :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp),    intent(in)    :: gpst1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(out)   :: ellum(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    integer(ip), intent(in)    :: pbubl(VECTOR_SIZE)
    real(rp),    intent(in)    :: gpsha_bub(VECTOR_SIZE,pgaus)                    ! Ne
    real(rp),    intent(in)    :: gpcar_bub(VECTOR_SIZE,ndime,pgaus)              ! dNe/dxi
    real(rp),    intent(in)    :: gppre(VECTOR_SIZE,pgaus) 
    ! Enrichement Element matrices
    real(rp),    intent(out)   :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp),    intent(out)   :: elapq(VECTOR_SIZE,pnode,1)
    real(rp),    intent(out)   :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp),    intent(out)   :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp),    intent(out)   :: elaqq(VECTOR_SIZE,1,1)
    real(rp),    intent(out)   :: elrbq(VECTOR_SIZE,1)
    ! Indices and local variables
    integer(ip)                :: inode,jnode,idofv
    integer(ip)                :: idof1,jdime,ivect
    integer(ip)                :: igaus,idime,jdofv,kdime
    real(rp)                   :: xvis2,xvisc,one_rp
    ! Vectorized local arrays
    real(rp)                   :: p1ve2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp)                   :: p1vec(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: p2vec(VECTOR_SIZE,ndime,pnode,pgaus)
    real(rp)                   :: p2sca(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: wgrgr(VECTOR_SIZE,pnode,pnode,pgaus)
    real(rp)                   :: wgrvi(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: fact0(VECTOR_SIZE),fact1(VECTOR_SIZE)
    real(rp)                   :: fact2(VECTOR_SIZE),fact3(VECTOR_SIZE)
    real(rp)                   :: fact4(VECTOR_SIZE),fact5(VECTOR_SIZE)
    real(rp)                   :: fact6(VECTOR_SIZE),fact7(VECTOR_SIZE)
    real(rp)                   :: fact8(VECTOR_SIZE),factz(VECTOR_SIZE)
    real(rp)                   :: factx(VECTOR_SIZE),facty(VECTOR_SIZE)
    real(rp)                   :: facx1(VECTOR_SIZE),facy1(VECTOR_SIZE)
    real(rp)                   :: ugraN(VECTOR_SIZE)
    real(rp)                   :: gramugraN(VECTOR_SIZE)
    real(rp)                   :: penal(VECTOR_SIZE)
    real(rp)                   :: gprhh(VECTOR_SIZE,ndime,pgaus)
    real(rp)                   :: taupr(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpveo(VECTOR_SIZE,ndime)
    real(rp)                   :: gpcar1ji(VECTOR_SIZE)
    real(rp)                   :: gpcar2ji(VECTOR_SIZE)
    real(rp)                   :: gpcar3ji(VECTOR_SIZE)
    real(rp)                   :: p2sca_bub(VECTOR_SIZE,pgaus)
    real(rp)                   :: p2vec_bub(VECTOR_SIZE,ndime,pgaus)

    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu11,elauu21,elauu31
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu12,elauu22,elauu32
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu13,elauu23,elauu33
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4))       :: factvec1,factvec2,factvec3,factvec4

    !yor! the comments marked with 'yor' explain the recent changes for optimization purpose, within EoCoE project. To be removed later.
    !yor! optimization test : for memory padding purpose : 4*((pnode+3)/4) will give the closest number >= pnode which is a multiple of 4.
    !yor! We have to think of a clean way to write it. The idea could be expanded in the whole code, with the padding size as parameter.

#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    !----------------------------------------------------------------------
    !
    ! Local variables
    !
    !----------------------------------------------------------------------
    !
    ! Viscous term factor
    !
    xvis2 = 0.0_rp
    if( fvins_nsi > 0.9_rp ) then
       xvisc = 1.0_rp
       if( fvins_nsi > 1.9_rp ) xvis2 = 2.0_rp/3.0_rp
    else
       xvisc = 0.0_rp
    end if
    !
    ! Do not stabilize convection, reaction, coriolis
    !
    if( kfl_nota1_nsi == 1 ) then
       one_rp = 0.0_rp 
    else
       one_rp = 1.0_rp
    end if

#ifdef OPENACC
    do ivect = 1,VECTOR_SIZE
#endif

    !----------------------------------------------------------------------
    !
    ! Initialization
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elaup(DEF_VECT,:,:)   = 0.0_rp
    elapu(DEF_VECT,:,:)   = 0.0_rp
    elapp(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    elrbp(DEF_VECT,:)     = 0.0_rp
 
    elauu11(DEF_VECT,:,:) = 0.0_rp
    elauu21(DEF_VECT,:,:) = 0.0_rp
    elauu12(DEF_VECT,:,:) = 0.0_rp
    elauu22(DEF_VECT,:,:) = 0.0_rp
    elauu31(DEF_VECT,:,:) = 0.0_rp
    elauu32(DEF_VECT,:,:) = 0.0_rp
    elauu13(DEF_VECT,:,:) = 0.0_rp
    elauu23(DEF_VECT,:,:) = 0.0_rp
    elauu33(DEF_VECT,:,:) = 0.0_rp
    !
    ! Incompressible vs low Mach
    !
    if( kfl_regim_nsi == 3 ) then
       taupr(DEF_VECT,1:pgaus) = 1.0_rp
    else
       taupr(DEF_VECT,1:pgaus) = 1.0_rp / max(gpst1(DEF_VECT,1:pgaus),zeror)
    end if
    ! Nest1
    !----------------------------------------------------------------------
    !
    ! Test functions
    !
    !----------------------------------------------------------------------
    !
    ! P1VEC =  v * [ (tau1'/tau1) - tau1' * sig ] + tau1' * rho * (uc.grad)v
    ! P2SCA =  ( tau2^{-1} * tau2' ) * q   
    ! P2VEC =  tau1' * grad(q)  ( and * rho if in Low-Mach, but not yet)
    ! WGRVI =  grad(Ni) . grad(mu) + mu * Delta(Ni)
    ! WGRGR =  grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus

       fact1(DEF_VECT) = gpsp1(DEF_VECT,igaus) * one_rp * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
       fact2(DEF_VECT) = (gptt1(DEF_VECT,igaus)-gpsp1(DEF_VECT,igaus) * one_rp * gppor(DEF_VECT,igaus)) * gpvol(DEF_VECT,igaus)
       fact3(DEF_VECT) = gpsp1(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       fact4(DEF_VECT) = gptt2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode

          ugraN(DEF_VECT)     = 0.0_rp
          gramugraN(DEF_VECT) = 0.0_rp

          do idime = 1,ndime
             ugraN(DEF_VECT)                   = ugraN(DEF_VECT)     + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             gramugraN(DEF_VECT)               = gramugraN(DEF_VECT) + gpgvi(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             p2vec(DEF_VECT,idime,inode,igaus) = fact3(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
          end do

          p1vec(DEF_VECT,inode,igaus)          = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) + fact1(DEF_VECT) * ugraN(DEF_VECT)
          p2sca(DEF_VECT,inode,igaus)          = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
          wgrvi(DEF_VECT,inode,igaus)          = gramugraN(DEF_VECT) + gpvis(DEF_VECT,igaus) * gplap(DEF_VECT,inode,igaus)

          do jnode = 1,pnode           
             wgrgr(DEF_VECT,inode,jnode,igaus) = 0.0_rp
             do idime = 1,ndime
                wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
             end do
          end do

       end do
    end do
    ! Nest2
    !
    ! P1VE2: Off-diagonal test function
    !
    if( kfl_p1ve2_nsi /= 0 ) then
       !
       ! tau1'*2*rho*(w x v)
       ! x-equation: v = (v,0,0) => w x v = (    0 , wz*v , -wy*v)     
       ! y-equation: v = (0,v,0) => w x v = (-wz*v ,    0 ,  wx*v)     
       ! z-equation: v = (0,0,v) => w x v = ( wy*v ,-wx*v ,     0)     
       ! 
       p1ve2 = 0.0_rp
       if( ndime == 2 ) then
          do igaus = 1,pgaus           
             fact3(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus) * fvela_nsi(3)
             factz(DEF_VECT) = gpsp1(DEF_VECT,igaus)* fact3(DEF_VECT) * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode 
                fact1(DEF_VECT)                 =   factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                p1ve2(DEF_VECT,1,2,inode,igaus) =   fact1(DEF_VECT)
                p1ve2(DEF_VECT,2,1,inode,igaus) = - fact1(DEF_VECT)
             end do
          end do
       else
          do igaus = 1,pgaus   
             fact8(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus) * gpsp1(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
             factx(DEF_VECT) = fact8(DEF_VECT)  * fvela_nsi(1)
             facty(DEF_VECT) = fact8(DEF_VECT)  * fvela_nsi(2)
             factz(DEF_VECT) = fact8(DEF_VECT)  * fvela_nsi(3)
             do inode = 1,pnode
                p1ve2(DEF_VECT,1,2,inode,igaus) =  factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! x-equation
                p1ve2(DEF_VECT,1,3,inode,igaus) = -facty(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! x-equation
                p1ve2(DEF_VECT,2,1,inode,igaus) = -factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! y-equation
                p1ve2(DEF_VECT,2,3,inode,igaus) =  factx(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! y-equation
                p1ve2(DEF_VECT,3,1,inode,igaus) =  facty(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! z-equation
                p1ve2(DEF_VECT,3,2,inode,igaus) = -factx(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! z-equation    
             end do
          end do
       end if
    end if

    ! Nest3
    !----------------------------------------------------------------------
    !
    ! Auu
    !
    !----------------------------------------------------------------------
    !
    !  Laplacian  form: A=0, B=0, eps(u) = 1/2 grad(u)
    !  Divergence form: A=1, B=0, eps(u) = 1/2 ( grad(u) + grad(u)^t )
    !  Complete   form: A=1, B=1, eps(u) = 1/2 ( grad(u) + grad(u)^t ) - 1/3 div(u) I
    !
    !    ( div(u) , tau2' div(v) )             (1)   <=   tau2' * duj/dxj * dvi/dxi
    !  + ( rmomu(DEF_VECT,u) , p1vec(DEF_VECT,v) )               (2)       
    !  + ( mu dui/dxj , dvi/dxj )              (3)   <=   ( 2 mu eps(u) , grad(v) ) = mu * ( dui/dxj + duj/dxi - 2/3 div(u) delta_ij ) dvi/dxj
    !  + A * ( mu duj/dxi , dvi/dxj )          (4)        divergence
    !  - 2/3 * B * ( mu (div u) , dvi/dxi )    (5)        complete
    !  + ( mu d^2ui/dxj^2 , vi )               (6)   <= - ( div[-2 mu eps(u)] , v ) = - d/dxj ( -2*mu* ( dui/dxj + duj/dxi - 2/3*mu div(u) delta_ij ) * vi
    !  + ( dmu/dxj dui/dxj , vi )              (7)        laplacian
    !  + A * ( mu  d^2uj/dxidxj , vi )         (8)        divergence                           
    !  + A * ( dmu/dxj duj/dxi , vi )          (9)        divergence
    !  - 2/3 * B * ( dmu/dxi (div u) , vi )   (10)        complete
    !  - 2/3 * B * ( mu d(div u)/dxi , vi )   (11)        complete         
    !  
    !  1. The terms ( div(u) , tau2' div(v) ) and - 2/3 * B * ( mu div u, dvi/dxj ) are assembled together as
    !     ( ( tau2' - 2/3 B mu ) div u, dvi/dxi ). Therefore, term (5) is assembled always although it is zero
    !  2. Terms (3), (6) and (7) are assembled just below as they are the same for all momentum equations 
    !  3. Terms (4), (8), (9), (10) and (11) are assmbled later on
    !
    !  (3)        x:     mu * ( dux/dx * dv/dx  +  dux/dy * dv/dy  +  dux/dz * dv/dz ) 
    !             y:     mu * ( duy/dx * dv/dx  +  duy/dy * dv/dy  +  duy/dz * dv/dz ) 
    !             z:     mu * ( duz/dx * dv/dx  +  duz/dy * dv/dy  +  duz/dz * dv/dz )   
    !
    !  (4)        x: A * mu * ( dux/dx * dv/dx  +  duy/dx * dv/dy  +  duz/dx * dv/dz )
    !             y: A * mu * ( dux/dy * dv/dx  +  duy/dy * dv/dy  +  duz/dy * dv/dz )
    !             z: A * mu * ( dux/dz * dv/dx  +  duy/dz * dv/dy  +  duz/dz * dv/dz )
    !
    !  (6)        x:  mu * ( d^2ux/dx^2  +  d^2ux/dy^2  +  d^2ux/dz^2 ) v 
    !             y:  mu * ( d^2uy/dx^2  +  d^2uy/dy^2  +  d^2uy/dz^2 ) v 
    !             z:  mu * ( d^2uz/dx^2  +  d^2uz/dy^2  +  d^2uz/dz^2 ) v 
    !
    !  (7)        x: ( dmu/dx * dux/dx  +  dmu/dy * dux/dy  +  dmu/dz * dux/dz ) v 
    !             y: ( dmu/dx * duy/dx  +  dmu/dy * duy/dy  +  dmu/dz * duy/dz ) v 
    !             y: ( dmu/dx * duz/dx  +  dmu/dy * duz/dy  +  dmu/dz * duz/dz ) v 
    !
    !  (8)        x: A * mu * ( d^2ux/dxdx  +  d^2uy/dxdy  +  d^2uz/dxdz ) v
    !             y: A * mu * ( d^2ux/dxdy  +  d^2uy/dydy  +  d^2uz/dydz ) v
    !             z: A * mu * ( d^2ux/dxdz  +  d^2uy/dzdy  +  d^2uz/dzdz ) v
    !
    !  (9)        x: A * ( dmu/dx * dux/dx  +  dmu/dy * duy/dx  +  dmu/dz * duz/dx ) v
    !             y: A * ( dmu/dx * dux/dy  +  dmu/dy * duy/dy  +  dmu/dz * duz/dy ) v
    !             z: A * ( dmu/dx * dux/dz  +  dmu/dy * duy/dz  +  dmu/dz * duz/dz ) v
    !
    !  (10)       x: -2/3 B * dmu/dx ( dux/dx  +  duy/dy  +  duz/dz ) * v
    !             y: -2/3 B * dmu/dy ( dux/dx  +  duy/dy  +  duz/dz ) * v
    !             z: -2/3 B * dmu/dz ( dux/dx  +  duy/dy  +  duz/dz ) * v
    !
    !  (11)       x: -2/3 B * mu * ( d^2ux/dxdx + d^2uy/dydx + d^2uz/dzdx ) v
    !             y: -2/3 B * mu * ( d^2ux/dxdy + d^2uy/dydy + d^2uz/dzdy ) v
    !             z: -2/3 B * mu * ( d^2ux/dxdz + d^2uy/dydz + d^2uz/dzdz ) v
    !
    !  (1) + (5)  x: ( tau2' - 2/3 B mu ) ( dux/dx  +  duy/dy  +  duz/dz ) * dv/dx
    !             y: ( tau2' - 2/3 B mu ) ( dux/dx  +  duy/dy  +  duz/dz ) * dv/dy
    !             z: ( tau2' - 2/3 B mu ) ( dux/dx  +  duy/dy  +  duz/dz ) * dv/dz
    !
    if( ndime == 2 ) then

       do igaus = 1,pgaus

          fact0(DEF_VECT) = ( gpsp2(DEF_VECT,igaus) - xvis2 * gpvis(DEF_VECT,igaus) ) * gpvol(DEF_VECT,igaus)                                      ! (1) + (5)
          fact6(DEF_VECT) = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

          do inode = 1,pnode
             fact1(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) 
             fact2(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus)
             fact5(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)

             do jnode = 1,pnode              
                fact4(DEF_VECT)               =   p1vec(DEF_VECT,inode,igaus) * rmomu(DEF_VECT,jnode,igaus) &                                          ! (2):       ( rmomu(DEF_VECT,u) , p1vec(DEF_VECT,v) ) 
                     &                          + fact5(DEF_VECT) * wgrvi(DEF_VECT,jnode,igaus)              &                                         ! (3) + (6): ( mu d^2ui/dxj^2 + dmu/dxj dui/dxj, vi )
                     &                          + fact6(DEF_VECT) * wgrgr(DEF_VECT,inode,jnode,igaus)                                                  ! (7):       ( mu dui/dxj , dvi/dxj ) 

                elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * rcont(DEF_VECT,1,jnode,igaus) + fact4(DEF_VECT)      ! Auu_xx
                elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * rcont(DEF_VECT,1,jnode,igaus)                        ! Auu_yx
                elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * rcont(DEF_VECT,2,jnode,igaus)                        ! Auu_xy
                elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * rcont(DEF_VECT,2,jnode,igaus) + fact4(DEF_VECT)      ! Auu_yy
             end do

          end do

       end do

    else

       do igaus = 1,pgaus    

          fact0(DEF_VECT) = ( gpsp2(DEF_VECT,igaus) - xvis2 * gpvis(DEF_VECT,igaus) ) * gpvol(DEF_VECT,igaus)                                          ! (1) + (5)
          fact6(DEF_VECT) = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

          factvec1(DEF_VECT,1:pnode) = gpcar(DEF_VECT,1,1:pnode,igaus)
          factvec2(DEF_VECT,1:pnode) = gpcar(DEF_VECT,2,1:pnode,igaus)
          factvec3(DEF_VECT,1:pnode) = gpcar(DEF_VECT,3,1:pnode,igaus)

          do jnode = 1,pnode              
             fact1(DEF_VECT) = fact0(DEF_VECT)*rcont(DEF_VECT,1,jnode,igaus) 
             fact2(DEF_VECT) = fact0(DEF_VECT)*rcont(DEF_VECT,2,jnode,igaus)
             fact3(DEF_VECT) = fact0(DEF_VECT)*rcont(DEF_VECT,3,jnode,igaus)
             fact5(DEF_VECT) = wgrvi(DEF_VECT,jnode,igaus)*gpvol(DEF_VECT,igaus)

             do inode = 1,pnode
                factvec4(DEF_VECT,inode) = p1vec(DEF_VECT,inode,igaus) * rmomu(DEF_VECT,jnode,igaus) + fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) &
                     + fact6(DEF_VECT) * wgrgr(DEF_VECT,inode,jnode,igaus)
             end do

             do inode = 1,pnode
                elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * factvec1(DEF_VECT,inode) + factvec4(DEF_VECT,inode)  ! Auu_xx
                elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * factvec1(DEF_VECT,inode)                             ! Auu_xy
                elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + fact3(DEF_VECT) * factvec1(DEF_VECT,inode)                             ! Auu_xz

                elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * factvec2(DEF_VECT,inode)                             ! Auu_yx
                elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * factvec2(DEF_VECT,inode) + factvec4(DEF_VECT,inode)  ! Auu_yy
                elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + fact3(DEF_VECT) * factvec2(DEF_VECT,inode)                             ! Auu_yz

                elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * factvec3(DEF_VECT,inode)                             ! Auu_zx
                elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * factvec3(DEF_VECT,inode)                             ! Auu_zy
                elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + fact3(DEF_VECT) * factvec3(DEF_VECT,inode) + factvec4(DEF_VECT,inode)  ! Auu_zz
             end do

          end do

       end do

    end if
    ! Nest4
    if( fvins_nsi > 0.9_rp ) then

       if( ndime == 2 ) then

          do igaus = 1,pgaus

             fact4(DEF_VECT) = gpgvi(DEF_VECT,1,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = gpgvi(DEF_VECT,2,igaus) * gpvol(DEF_VECT,igaus)
             fact7(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

             do inode = 1,pnode
                fact0(DEF_VECT) = fact7(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) *(1.0_rp-xvis2)
                fact1(DEF_VECT) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                fact2(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                facx1(DEF_VECT) = fact1(DEF_VECT) * xvis2
                facy1(DEF_VECT) = fact2(DEF_VECT) * xvis2    
                do jnode = 1,pnode              
                   !
                   !                                         (4)  ( mu duj/dxi , dvi/dxj )
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                   ! Auu_xx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                   ! Auu_xy
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                   ! Auu_yx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                   ! Auu_yy
                   !
                   !                                         (10) - 2/3 * B * ( dmu/dxi (div u) , vi )
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - facx1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                   ! Auu_xx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) - facx1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                   ! Auu_xy
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) - facy1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                   ! Auu_yx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - facy1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                   ! Auu_yy
                   !
                   !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
                   !                                               + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi )
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,1,jnode,igaus) ! Auu_xx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,jnode,igaus) ! Auu_xy
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,jnode,igaus) ! Auu_yx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,2,jnode,igaus) ! Auu_yy
                   !
                   !                                        

                end do

             end do
          end do

       else

          do igaus = 1,pgaus

             fact4(DEF_VECT) = gpgvi(DEF_VECT,1,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = gpgvi(DEF_VECT,2,igaus) * gpvol(DEF_VECT,igaus)
             fact6(DEF_VECT) = gpgvi(DEF_VECT,3,igaus) * gpvol(DEF_VECT,igaus)
             fact7(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)
             !          factx(DEF_VECT) = xvis2 * fact4(DEF_VECT)
             !          facty(DEF_VECT) = xvis2 * fact5(DEF_VECT)
             !          factz(DEF_VECT) = xvis2 * fact6(DEF_VECT)

             !yor! [Originally line 451] reverted loops order, refactored the factors, with an effort to vectorize each instruction 
             !
             !                                (4)  ( mu duj/dxi , dvi/dxj )
             factvec1(DEF_VECT,1:pnode) = gpcar(DEF_VECT,1,1:pnode,igaus)
             factvec2(DEF_VECT,1:pnode) = gpcar(DEF_VECT,2,1:pnode,igaus)
             factvec3(DEF_VECT,1:pnode) = gpcar(DEF_VECT,3,1:pnode,igaus)

             do jnode = 1,pnode
                gpcar1ji(DEF_VECT) = gpcar(DEF_VECT,1,jnode,igaus)
                gpcar2ji(DEF_VECT) = gpcar(DEF_VECT,2,jnode,igaus)
                gpcar3ji(DEF_VECT) = gpcar(DEF_VECT,3,jnode,igaus)
                !
                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec1(DEF_VECT,inode) * gpcar1ji(DEF_VECT)         ! Auu_xx
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec1(DEF_VECT,inode) * gpcar2ji(DEF_VECT)         ! Auu_yx
                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec1(DEF_VECT,inode) * gpcar3ji(DEF_VECT)         ! Auu_zx

                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec2(DEF_VECT,inode) * gpcar1ji(DEF_VECT)         ! Auu_xy
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec2(DEF_VECT,inode) * gpcar2ji(DEF_VECT)         ! Auu_yy
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec2(DEF_VECT,inode) * gpcar3ji(DEF_VECT)         ! Auu_zy

                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec3(DEF_VECT,inode) * gpcar1ji(DEF_VECT)         ! Auu_xz
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec3(DEF_VECT,inode) * gpcar2ji(DEF_VECT)         ! Auu_yz
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec3(DEF_VECT,inode) * gpcar3ji(DEF_VECT)         ! Auu_zz
                end do
             end do
             !
             !                     (10) - 2/3 * B * ( dmu/dxi (div u) , vi )

             do inode = 1,pnode
                factvec1(DEF_VECT,inode) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * xvis2 !yor! (ex facx1) this copy is for pure esthaetical purpose
                factvec2(DEF_VECT,inode) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * xvis2 !yor! (ex facy1)
                factvec3(DEF_VECT,inode) = fact6(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * xvis2 !yor! (ex facz1)  
             end do

             do jnode = 1,pnode
                gpcar1ji(DEF_VECT) = gpcar(DEF_VECT,1,jnode,igaus)
                gpcar2ji(DEF_VECT) = gpcar(DEF_VECT,2,jnode,igaus)
                gpcar3ji(DEF_VECT) = gpcar(DEF_VECT,3,jnode,igaus)

                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - factvec1(DEF_VECT,inode) * gpcar1ji(DEF_VECT)                                ! Auu_xx 
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) - factvec1(DEF_VECT,inode) * gpcar2ji(DEF_VECT)                                ! Auu_xy
                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) - factvec1(DEF_VECT,inode) * gpcar3ji(DEF_VECT)                                ! Auu_xz

                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) - factvec2(DEF_VECT,inode) * gpcar1ji(DEF_VECT)                                ! Auu_yx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - factvec2(DEF_VECT,inode) * gpcar2ji(DEF_VECT)                                ! Auu_yy
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) - factvec2(DEF_VECT,inode) * gpcar3ji(DEF_VECT)                                ! Auu_yz

                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) - factvec3(DEF_VECT,inode) * gpcar1ji(DEF_VECT)                                ! Auu_zx
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) - factvec3(DEF_VECT,inode) * gpcar2ji(DEF_VECT)                                ! Auu_zy
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) - factvec3(DEF_VECT,inode) * gpcar3ji(DEF_VECT)                                ! Auu_zz
                end do
             end do

             !                                                                   
             !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
             !                                                + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi ) 
             do inode = 1,pnode
                factvec4(DEF_VECT,inode) = fact7(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)*(1.0_rp-xvis2) !yor!(fact0(DEF_VECT))
                factvec1(DEF_VECT,inode) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) !yor!(ex fact1(DEF_VECT))
                factvec2(DEF_VECT,inode) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) !yor!(ex fact2(DEF_VECT))
                factvec3(DEF_VECT,inode) = fact6(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) !yor!(ex fact3(DEF_VECT))
             end do

             do jnode = 1,pnode
                gpcar1ji(DEF_VECT) = gpcar(DEF_VECT,1,jnode,igaus)
                gpcar2ji(DEF_VECT) = gpcar(DEF_VECT,2,jnode,igaus)
                gpcar3ji(DEF_VECT) = gpcar(DEF_VECT,3,jnode,igaus)

                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode) * gpcar1ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,1,jnode,igaus) ! Auu_xx
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode) * gpcar2ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,4,jnode,igaus) ! Auu_yx
                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode) * gpcar3ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,5,jnode,igaus) ! Auu_zx

                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + factvec2(DEF_VECT,inode) * gpcar1ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,4,jnode,igaus) ! Auu_xy
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + factvec2(DEF_VECT,inode) * gpcar2ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,2,jnode,igaus) ! Auu_yy
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + factvec2(DEF_VECT,inode) * gpcar3ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,6,jnode,igaus) ! Auu_zy

                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + factvec3(DEF_VECT,inode) * gpcar1ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,5,jnode,igaus) ! Auu_xz
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + factvec3(DEF_VECT,inode) * gpcar2ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,6,jnode,igaus) ! Auu_yz
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + factvec3(DEF_VECT,inode) * gpcar3ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,3,jnode,igaus) ! Auu_zz
                end do
             end do
          end do
       end if
    end if

    ! Nest5
    if( fcons_nsi > 0.1_rp ) then   
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode
                do jnode = 1,pnode     
                   fact1(DEF_VECT) = 0.0_rp
                   fact2(DEF_VECT) = 0.0_rp
                   do idime = 1,2
                      fact1(DEF_VECT) = fact1(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                      fact2(DEF_VECT) = fact2(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) 
                   end do
                   fact1(DEF_VECT) = fact1(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
                   fact2(DEF_VECT) = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT) 

                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_xx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_yy
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode
                do jnode = 1,pnode
                   fact1(DEF_VECT) = 0.0_rp
                   fact2(DEF_VECT) = 0.0_rp
                   do idime = 1,3
                      fact1(DEF_VECT) = fact1(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                      fact2(DEF_VECT) = fact2(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) 
                   end do
                   fact1(DEF_VECT) = fact1(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
                   fact2(DEF_VECT) = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT) 

                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_xx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_yy
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_zz
                end do
             end do
          end do
       end if
    end if

    ! Nest6
    !----------------------------------------------------------------------
    !
    ! Auu: Off-diagonal terms
    !
    !----------------------------------------------------------------------

    if( kfl_rmom2_nsi /= 0 ) then
       !
       ! p1vec * rmom2(DEF_VECT,
       !
       !yor! [Originally loop on line 577] Reordered loops to solve indirections problems 
       if (ndime == 2) then
          do igaus = 1,pgaus
             do jnode = 1,pnode
                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,1,1,jnode,igaus)
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,2,1,jnode,igaus)
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,1,2,jnode,igaus)
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,2,2,jnode,igaus)            
                end do
             end do
          end do
       else !ndime==3
          do igaus = 1,pgaus
             factvec1(DEF_VECT,1:pnode) = p1vec(DEF_VECT,1:pnode,igaus)
             do jnode = 1,pnode
                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,1,jnode,igaus)
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,1,jnode,igaus)
                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,1,jnode,igaus)
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,2,jnode,igaus)
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,2,jnode,igaus)
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,2,jnode,igaus)
                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,3,jnode,igaus)
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,3,jnode,igaus)
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,3,jnode,igaus)
                end do
             end do
          end do
       end if

       if( kfl_p1ve2_nsi /= 0 ) then
          !
          ! p1ve2 * rmomu
          ! p1ve2 * rmom2(DEF_VECT,
          !
          !yor! optim : [originally line 597] reverted loops, refactored 
          if(ndime==2)then   
             do igaus=1,pgaus
                do jnode=1,pnode
                   do inode=1,pnode
                      elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      do kdime = 1,ndime
                         elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                         elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                      end do
                   end do
                end do
             end do
          else ! ndime==3
             do igaus=1,pgaus
                do jnode=1,pnode
                   do inode=1,pnode
                      elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu31(DEF_VECT,inode,jnode)=elauu31(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu32(DEF_VECT,inode,jnode)=elauu32(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu13(DEF_VECT,inode,jnode)=elauu13(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,3,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu23(DEF_VECT,inode,jnode)=elauu23(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,3,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu33(DEF_VECT,inode,jnode)=elauu33(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,3,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      do kdime = 1,ndime
                         elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu31(DEF_VECT,inode,jnode)=elauu31(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                         elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                         elauu32(DEF_VECT,inode,jnode)=elauu32(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                         elauu13(DEF_VECT,inode,jnode)=elauu13(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,3,jnode,igaus)
                         elauu23(DEF_VECT,inode,jnode)=elauu23(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,3,jnode,igaus)
                         elauu33(DEF_VECT,inode,jnode)=elauu33(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,3,jnode,igaus)
                      end do
                   end do
                end do
             end do
          end if
       end if
    end if
    ! Nest7
    !
    ! Bemol:
    ! - b * ( v , rho * (a.grad) u ) - b * ( u , rho * (a.grad) v ) - b ( rho*div(a) u , v )
    !
    if( abs(bemol_nsi).gt.1.0e-9_rp ) then
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             fact0(DEF_VECT) = bemol_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = fact0(DEF_VECT) * ( gpgve(DEF_VECT,1,1,igaus) + gpgve(DEF_VECT,2,2,igaus) )
             do inode = 1,pnode
                do jnode = 1,pnode              
                   !
                   ! fact2(DEF_VECT) = - b * ( v , rho * (a.grad) u )
                   ! fact3(DEF_VECT) = - b * ( u , rho * (a.grad) v ) 
                   ! fact4(DEF_VECT) = - b * ( rho * div(a) u , v ) 
                   !
                   fact2(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus) ) * gpsha(DEF_VECT,jnode,igaus)
                   fact3(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,jnode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,jnode,igaus) ) * gpsha(DEF_VECT,inode,igaus)
                   fact4(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus) 
                   fact2(DEF_VECT) = fact0(DEF_VECT) * fact2(DEF_VECT)
                   fact3(DEF_VECT) = fact0(DEF_VECT) * fact3(DEF_VECT)

                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             fact0(DEF_VECT) = bemol_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = fact0(DEF_VECT) * ( gpgve(DEF_VECT,1,1,igaus) + gpgve(DEF_VECT,2,2,igaus) + gpgve(DEF_VECT,3,3,igaus))
             do inode = 1,pnode
                do jnode = 1,pnode
                   !
                   ! fact2(DEF_VECT) = - b * ( v , rho * (a.grad) u )
                   ! fact3(DEF_VECT) = - b * ( u , rho * (a.grad) v ) 
                   ! fact4(DEF_VECT) = - b * ( rho * div(a) u , v ) 
                   !
                   fact2(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus) &
                        + gpadv(DEF_VECT,3,igaus) * gpcar(DEF_VECT,3,inode,igaus)) * gpsha(DEF_VECT,jnode,igaus)
                   fact3(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,jnode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,jnode,igaus) &
                        + gpadv(DEF_VECT,3,igaus) * gpcar(DEF_VECT,3,jnode,igaus)) * gpsha(DEF_VECT,inode,igaus)
                   fact4(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                   fact2(DEF_VECT) = fact0(DEF_VECT) * fact2(DEF_VECT)
                   fact3(DEF_VECT) = fact0(DEF_VECT) * fact3(DEF_VECT)

                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                end do
             end do
          end do
       end if

    end if
    ! Nest8
    !
    ! Newton-Raphson
    !
    if( kfl_linea_nsi == 2 ) then
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             do inode = 1,pnode
                fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                do jnode = 1,pnode
                   fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   ! Auu_1i
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus) ! rho * ux * dux^i-1/dx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus) ! rho * uy * dux^i-1/dy
                   ! Auu_2i
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus) ! rho * ux * duy^i-1/dx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus) ! rho * uy * duy^i-1/dy
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             do inode = 1,pnode
                fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                do jnode = 1,pnode
                   fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   ! Auu_1i
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus) ! rho * ux * dux^i-1/dx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus) ! rho * uy * dux^i-1/dy
                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,1,igaus) ! rho * uz * dux^i-1/dz
                   ! Auu_2i
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus) ! rho * ux * duy^i-1/dx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus) ! rho * uy * duy^i-1/dy
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,2,igaus) ! rho * uz * duy^i-1/dz
                   ! Auu_3i
                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,3,igaus) ! rho * ux * duz^i-1/dx
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,3,igaus) ! rho * uy * duz^i-1/dy
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,3,igaus) ! rho * uz * duz^i-1/dz

                end do
             end do
          end do
       end if
    end if
    ! Nest9
    !
    ! Lumped evolution matrix (only backward euler)
    !
    if( kfl_lumped == 1 ) then 
       !
       ! Remove Galerkin term and add lumped term 
       ! 
       if( ndime == 2 ) then
          stop
       else
          do igaus =1, pgaus
             gpveo(DEF_VECT,1:3) = 0.0_rp
             do inode = 1,pnode
                gpveo(DEF_VECT,1) = gpveo(DEF_VECT,1) + elvel(DEF_VECT,1, inode, 2)* gpsha(DEF_VECT,inode, igaus)
                gpveo(DEF_VECT,2) = gpveo(DEF_VECT,2) + elvel(DEF_VECT,2, inode, 2)* gpsha(DEF_VECT,inode, igaus)
                gpveo(DEF_VECT,3) = gpveo(DEF_VECT,3) + elvel(DEF_VECT,3, inode, 2)* gpsha(DEF_VECT,inode, igaus)
             end do
             do inode =1, pnode
                fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)*dtinv_loc(DEF_VECT)
                elauu11(DEF_VECT,inode,inode) =  elauu11(DEF_VECT,inode,inode) + fact0(DEF_VECT)
                elauu22(DEF_VECT,inode,inode) =  elauu22(DEF_VECT,inode,inode) + fact0(DEF_VECT)
                elauu33(DEF_VECT,inode,inode) =  elauu33(DEF_VECT,inode,inode) + fact0(DEF_VECT)
                do idime = 1,ndime
                   elrbu(DEF_VECT,idime,inode)     =  elrbu(DEF_VECT,idime,inode)   - fact0(DEF_VECT)*gpveo(DEF_VECT,idime)
                   elrbu(DEF_VECT,idime,inode)     =  elrbu(DEF_VECT,idime,inode)   + fact0(DEF_VECT)*elvel(DEF_VECT,idime, inode, 2)
                end do
                do jnode =1, pnode !yor! for this case, better write a second loop and revert inode and jnode
                   elauu11(DEF_VECT,inode,jnode) =  elauu11(DEF_VECT,inode,jnode) - fact0(DEF_VECT)*gpsha(DEF_VECT,jnode, igaus) 
                   elauu22(DEF_VECT,inode,jnode) =  elauu22(DEF_VECT,inode,jnode) - fact0(DEF_VECT)*gpsha(DEF_VECT,jnode, igaus) 
                   elauu33(DEF_VECT,inode,jnode) =  elauu33(DEF_VECT,inode,jnode) - fact0(DEF_VECT)*gpsha(DEF_VECT,jnode, igaus) 
                end do
             end do
          end do
       end if
    else if( kfl_lumped == 2 ) then 
       !
       ! No time term have been added up to now: add Galerkin term
       !
       do igaus = 1,pgaus
          fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * dtinv_loc(DEF_VECT)
          if( ndime == 2 ) then
             do inode = 1, pnode
                elauu11(DEF_VECT,inode,inode) = elauu11(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elauu22(DEF_VECT,inode,inode) = elauu22(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             end do
          else ! ndime==3
             do inode = 1, pnode
                elauu11(DEF_VECT,inode,inode) = elauu11(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elauu22(DEF_VECT,inode,inode) = elauu22(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elauu33(DEF_VECT,inode,inode) = elauu33(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             end do
          end if
          do inode = 1, pnode
             do idime = 1,ndime   
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * elvel(DEF_VECT,idime,inode,2)
             end do
          end do
       end do
    end if
    ! Nest10
    !
    ! Dual time step preconditioner
    !
    if( kfl_duatss == 1 ) then
       ellum = 0.0_rp 
       do igaus =1, pgaus
          do inode =1, pnode
             fact0(DEF_VECT)       = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * dtinv_loc(DEF_VECT) * real(fact_duatss-1,rp)
             ellum(DEF_VECT,inode) = ellum(DEF_VECT,inode) + fact0(DEF_VECT) 
          end do
       end do
    end if

    ! Nest11
    !yor! elauu final matrix assembly
    if( ndime == 2 ) then  
       do jnode = 1,pnode
          do inode = 1,pnode
             elauu(DEF_VECT,2*inode-1,2*jnode-1)=elauu11(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,2*inode  ,2*jnode-1)=elauu21(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,2*inode-1,2*jnode  )=elauu12(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,2*inode  ,2*jnode  )=elauu22(DEF_VECT,inode,jnode)
          end do
       end do
    else ! ndime==3
       do jnode = 1,pnode
          do inode = 1,pnode
             elauu(DEF_VECT,3*inode-2,3*jnode-2)=elauu11(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-1,3*jnode-2)=elauu21(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode  ,3*jnode-2)=elauu31(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-2,3*jnode-1)=elauu12(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-1,3*jnode-1)=elauu22(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode  ,3*jnode-1)=elauu32(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-2,3*jnode  )=elauu13(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-1,3*jnode  )=elauu23(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode  ,3*jnode  )=elauu33(DEF_VECT,inode,jnode)
          end do
       end do
    end if
    ! Nest12
    !----------------------------------------------------------------------
    !
    ! Aup
    !
    !----------------------------------------------------------------------
    !
    ! Pressure: - ( p , div(v)  ) + ( grad(p) , p1vec(DEF_VECT,v)-v ) ( if kfl_press_nsi = 1 )
    ! Pressure: + ( grad(p) , p1vec(DEF_VECT,v) )                     ( if kfl_press_nsi = 0 )
    !  
    if( kfl_press_nsi == 1 ) then
       do igaus = 1,pgaus
          do jnode = 1,pnode
             fact2(DEF_VECT) = gpvol(DEF_VECT,igaus)*gpsha(DEF_VECT,jnode,igaus)
             do inode = 1,pnode
                fact1(DEF_VECT) = -gpvol(DEF_VECT,igaus)*gpsha(DEF_VECT,inode,igaus)+p1vec(DEF_VECT,inode,igaus)
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) - fact2(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) &
                        &                                  + fact1(DEF_VECT) * gpcar(DEF_VECT,idime,jnode,igaus)
                end do
             end do
          end do
       end do
    else
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) + p1vec(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                end do
             end do
          end do
       end do
    end if
    
    ! Nest13
    if( kfl_p1ve2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   do jdime = 1,ndime
                      elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) + &
                           gpcar(DEF_VECT,jdime,jnode,igaus) * p1ve2(DEF_VECT,idime,jdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end do
    end if

    ! Nest14
    !----------------------------------------------------------------------
    !
    ! Apu
    !
    !----------------------------------------------------------------------
    !
    ! ( div(u) , (tau2^{-1}*tau2')*q )   --- or ( div (rho u), q) if Low-Mach 
    ! + ( rho*(uc.grad)u + 2*rho*(w x u) + sig*u -div[2*mu*eps(u)], tau1' grad(q) )  (only laplacian form of div[2 mu eps(u)] )
    !
    if (kfl_regim_nsi == 3 .and. kfl_confi_nsi == 1) then ! Penalization of pressure, never imposse pressure
       do igaus =1, pgaus
          penal(DEF_VECT) = 1.0e-4_rp*gpden(DEF_VECT,igaus)/gpvis(DEF_VECT,igaus)
          do inode =1, pnode
             fact0(DEF_VECT) = penal(DEF_VECT) * p2sca(DEF_VECT,inode,igaus)
             do jnode =inode +1, pnode
                fact1(DEF_VECT) = fact0(DEF_VECT)*gpsha(DEF_VECT,jnode,igaus)
                elapp(DEF_VECT,inode,jnode) = elapp(DEF_VECT,inode,jnode) + fact1(DEF_VECT)
                elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + fact1(DEF_VECT)
             end do
             elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode, igaus)
             elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact0(DEF_VECT) * gppre(DEF_VECT,igaus)
          end do
       end do
    end if

    ! Nest15
#ifdef matiaslma

    if (ndime==2) then
       do igaus =1, pgaus
          do inode =1,pnode
             p2vec(DEF_VECT,1,inode,igaus)= p2vec(DEF_VECT,1,inode,igaus)*gpden(DEF_VECT,igaus)
             p2vec(DEF_VECT,2,inode,igaus)= p2vec(DEF_VECT,2,inode,igaus)*gpden(DEF_VECT,igaus)
             p2sca(DEF_VECT,inode, igaus) = p2sca(DEF_VECT,inode, igaus)*gpden(DEF_VECT,igaus)
          end do
          fact0(DEF_VECT) = gpden(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
          do inode =1, pnode           
             idof1 = 2*inode-1
             idof2 = idof1+1
             fact1(DEF_VECT) = gpsha(DEF_VECT,inode,igaus)* fact0(DEF_VECT)
             do jnode =1, pnode
                elapu(DEF_VECT,jnode, idof1) = elapu(DEF_VECT,jnode, idof1) -fact1(DEF_VECT)             * gpcar(DEF_VECT,1,jnode,igaus)  &
                     &                                  + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,1,jnode,igaus)
                elapu(DEF_VECT,jnode, idof2) = elapu(DEF_VECT,jnode, idof2) -fact1(DEF_VECT)             * gpcar(DEF_VECT,2,jnode,igaus)  &
                     &                                  + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,2,jnode,igaus)
             end do
          end do
       end do
    else !ndime =3
       do igaus =1, pgaus
          do inode =1,pnode
             p2vec(DEF_VECT,1:ndime,inode,igaus) = p2vec(DEF_VECT,1:ndime,inode,igaus)*gpden(DEF_VECT,igaus)
             p2sca(DEF_VECT,inode, igaus)        = p2sca(DEF_VECT,inode, igaus) *gpden(DEF_VECT,igaus)
          end do
          fact0(DEF_VECT) = gpden(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
          do inode =1, pnode
             idof1 = 3*inode - 2 
             idof2 = idof1+1
             idof3 = idof2+1
             fact1(DEF_VECT) = gpsha(DEF_VECT,inode,igaus)* fact0(DEF_VECT)
             do jnode =1, pnode 
                elapu(DEF_VECT,jnode, idof1) = elapu(DEF_VECT,jnode, idof1) - fact1(DEF_VECT)              * gpcar(DEF_VECT,1,jnode,igaus) &
                     &                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,1,jnode,igaus)
                elapu(DEF_VECT,jnode, idof2) = elapu(DEF_VECT,jnode, idof2) - fact1(DEF_VECT)              * gpcar(DEF_VECT,2,jnode,igaus) &
                     &                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,2,jnode,igaus)
                elapu(DEF_VECT,jnode, idof3) = elapu(DEF_VECT,jnode, idof3) - fact1(DEF_VECT)              * gpcar(DEF_VECT,3,jnode,igaus) &
                     &                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,3,jnode,igaus)
             end do
          end do
       end do
    end if

#else
    do igaus = 1,pgaus
       do inode = 1,pnode
          do jnode = 1,pnode
             do idime = 1,ndime
                idof1 = (inode-1)*ndime+idime
                elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) + rcont(DEF_VECT,idime,inode,igaus) * p2sca(DEF_VECT,jnode,igaus)&
                     &                                  + rmomu(DEF_VECT,inode,igaus)       * p2vec(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do
    end do

#endif
    ! Nest16
    !yor! optim [originally loop on line 985] reordered inner loops 
    if( kfl_rmom2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do jdime = 1,ndime
                jdofv = (jnode-1)*ndime+jdime
                do inode = 1,pnode
                   do kdime = 1,ndime
                      elapu(DEF_VECT,inode,jdofv) = elapu(DEF_VECT,inode,jdofv) + rmom2(DEF_VECT,kdime,jdime,jnode,igaus) * p2vec(DEF_VECT,kdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end do
    end if

    ! Nest17
    !----------------------------------------------------------------------
    !
    ! App
    !
    !----------------------------------------------------------------------
    !
    ! Pressure: ( grad(p) , tau1' grad(q) )
    ! 
    do igaus = 1,pgaus
       do inode = 1,pnode
          do jnode = 1,pnode
             do idime = 1,ndime
                elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + p2vec(DEF_VECT,idime,jnode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             end do
          end do
       end do
    end do
    !call nsi_element_system_output(&
    !     pnode,elauu(1,:,:),elaup(1,:,:),elapp(1,:,:),elapu(1,:,:),elrbu(1,:,:),elrbp(1,:),&
    !     elauq(1,:,:),elapq(1,:,:),elaqu(1,:,:),elaqp(1,:,:),elaqq(1,:,:),elrbq(1,:))
    !stop

    ! Nest18
    !
    ! Penalization
    !
    do igaus = 1,pgaus
       fact1(DEF_VECT) = penal_nsi * gpvol(DEF_VECT,igaus)
       do inode = 1,pnode
          elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
          elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * elpre(DEF_VECT,inode,1) 
       end do
    end do

    ! Nest19
    !----------------------------------------------------------------------
    !
    ! bu and bp
    !
    !----------------------------------------------------------------------
    !
    ! bu = ( f , p1vec(DEF_VECT,v) ) 
    ! bp = ( f , tau1' grad(q) )
    !
    if( ndime == 2 ) then
       do igaus = 1,pgaus
          gprhs_sgs(DEF_VECT,1,igaus) = gprhs_sgs(DEF_VECT,1,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,1,igaus)
          gprhs_sgs(DEF_VECT,2,igaus) = gprhs_sgs(DEF_VECT,2,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,2,igaus)
          gprhh(DEF_VECT,1,igaus)     = gprhs(DEF_VECT,1,igaus)     + gprhs_sgs(DEF_VECT,1,igaus)
          gprhh(DEF_VECT,2,igaus)     = gprhs(DEF_VECT,2,igaus)     + gprhs_sgs(DEF_VECT,2,igaus)
          do inode = 1,pnode
             elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,1,igaus)
             elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,2,igaus)
             elrbp(DEF_VECT,inode)   = elrbp(DEF_VECT,inode)   + p2vec(DEF_VECT,1,inode,igaus) * gprhh(DEF_VECT,1,igaus) &
                  &                                            + p2vec(DEF_VECT,2,inode,igaus) * gprhh(DEF_VECT,2,igaus) 
          end do
       end do
    else
       do igaus = 1,pgaus
          gprhs_sgs(DEF_VECT,1,igaus) = gprhs_sgs(DEF_VECT,1,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,1,igaus)
          gprhs_sgs(DEF_VECT,2,igaus) = gprhs_sgs(DEF_VECT,2,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,2,igaus)
          gprhs_sgs(DEF_VECT,3,igaus) = gprhs_sgs(DEF_VECT,3,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,3,igaus)
          gprhh(DEF_VECT,1, igaus)    = gprhs(DEF_VECT,1,igaus)     + gprhs_sgs(DEF_VECT,1,igaus)
          gprhh(DEF_VECT,2, igaus)    = gprhs(DEF_VECT,2,igaus)     + gprhs_sgs(DEF_VECT,2,igaus)
          gprhh(DEF_VECT,3, igaus)    = gprhs(DEF_VECT,3,igaus)     + gprhs_sgs(DEF_VECT,3,igaus)
          do inode = 1,pnode
             elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,1,igaus)
             elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,2,igaus)
             elrbu(DEF_VECT,3,inode) = elrbu(DEF_VECT,3,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,3,igaus)
             elrbp(DEF_VECT,inode)   = elrbp(DEF_VECT,inode)   + p2vec(DEF_VECT,1,inode,igaus) * gprhh(DEF_VECT,1,igaus) &
                  &                                            + p2vec(DEF_VECT,2,inode,igaus) * gprhh(DEF_VECT,2,igaus) &
                  &                                            + p2vec(DEF_VECT,3,inode,igaus) * gprhh(DEF_VECT,3,igaus)
          end do
       end do
    end if

    ! Nest20
    if( kfl_p1ve2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                do kdime = 1,ndime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                        + p1ve2(DEF_VECT,idime,kdime,inode,igaus) * gprhh(DEF_VECT,kdime,igaus)
                end do
             end do
          end do
       end do
    end if

    ! Nest21
    !
    ! low Mach regime: we add to the right hand side the residue of the continuity eq
    !
    do igaus = 1,pgaus
       fact0(DEF_VECT) = gprhc(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
       fact1(DEF_VECT) = gpsp2(DEF_VECT,igaus) * gprh2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
       do inode = 1,pnode
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact1(DEF_VECT)        * gpcar(DEF_VECT,idime,inode,igaus) ! ( rhs, tau2' * div(v) )
          end do
          elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)                             ! ( rhs , q)
       end do
    end do

    ! Nest22
    !
    ! Newton-Raphson
    !
    if( kfl_linea_nsi == 2 ) then
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             fact0(DEF_VECT) =  gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
             do inode = 1,pnode
                fact1(DEF_VECT)         = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,1,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,1,igaus) )
                elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,2,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,2,igaus) )
             end do
          end do
       else
          do igaus = 1,pgaus
             fact0(DEF_VECT) =  gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
             do inode = 1,pnode
                fact1(DEF_VECT)         = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) & 
                     + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,1,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,1,igaus) &
                     + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,1,igaus) )
                elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) &
                     + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,2,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,2,igaus) &
                     + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,2,igaus) )
                elrbu(DEF_VECT,3,inode) = elrbu(DEF_VECT,3,inode) &
                     + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,3,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,3,igaus) &
                     + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,3,igaus) )
             end do
          end do
       end if
    end if

    ! Nest23
    !
    ! Projection: ( L*(v) , P )
    !
    if( kfl_stabi_nsi == 1 ) then

       do igaus = 1, pgaus
          do inode = 1,pnode
             fact0(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
             do idime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) - gprhs_sgs(DEF_VECT,idime,igaus) * fact0(DEF_VECT)
             end do
          end do
       end do

    end if

#ifdef OPENACC
    end do
#endif

    !--------------------------------------------------------------------
    !
    ! Pressure bubble
    !
    !--------------------------------------------------------------------

    if( maxval(pbubl) == 1 ) then
       !
       ! Initialization
       !
       elauq = 0.0_rp
       elapq = 0.0_rp
       elaqu = 0.0_rp
       elaqp = 0.0_rp
       elaqq = 0.0_rp
       elrbq = 0.0_rp
       !
       ! Test functions
       !
       do igaus = 1,pgaus
          p2sca_bub(DEF_VECT,igaus) = gptt2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
          do idime = 1,ndime
             p2vec_bub(DEF_VECT,idime,igaus) = gpsp1(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
          end do
       end do
       ! 
       ! Momentum equations
       !
       if( kfl_press_nsi == 1 ) then
          do igaus = 1,pgaus
             fact2(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
             do inode = 1,pnode
                fact1(DEF_VECT) = -gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) + p1vec(DEF_VECT,inode,igaus)
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) - fact2(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) &
                        &                                            + fact1(DEF_VECT) * gpcar_bub(DEF_VECT,idime,igaus)
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + p1vec(DEF_VECT,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                end do
             end do
          end do
       end if

       if( kfl_p1ve2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   do jdime = 1,ndime
                      elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + &
                           gpcar_bub(DEF_VECT,jdime,igaus) * p1ve2(DEF_VECT,idime,jdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end if
       !
       ! Bubble equation
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                idof1 = (inode-1)*ndime+idime
                elaqu(DEF_VECT,1,idof1) = elaqu(DEF_VECT,1,idof1) + rcont(DEF_VECT,idime,inode,igaus) * p2sca_bub(DEF_VECT,igaus) &
                     &                                            + rmomu(DEF_VECT,inode,igaus)       * p2vec_bub(DEF_VECT,idime,igaus)
             end do
          end do
       end do
       if( kfl_rmom2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do jnode = 1,pnode
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   do kdime = 1,ndime
                      elaqu(DEF_VECT,1,jdofv) = elaqu(DEF_VECT,1,jdofv) + rmom2(DEF_VECT,kdime,jdime,jnode,igaus) * p2vec_bub(DEF_VECT,kdime,igaus)
                   end do
                end do
             end do
          end do
       end if
       !
       ! Also pressure equation...
       !
       do igaus = 1,pgaus 
          do inode = 1,pnode
             do idime = 1,ndime
                elapq(DEF_VECT,inode,1) = elapq(DEF_VECT,inode,1) + p2vec(DEF_VECT,idime,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                elaqp(DEF_VECT,1,inode) = elaqp(DEF_VECT,1,inode) + p2vec_bub(DEF_VECT,idime,igaus)   * gpcar(DEF_VECT,idime,inode,igaus)
             end do
          end do
          do idime = 1,ndime
             elaqq(DEF_VECT,1,1)             = elaqq(DEF_VECT,1,1)             + p2vec_bub(DEF_VECT,idime,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
             gprhs_sgs(DEF_VECT,idime,igaus) = gprhs_sgs(DEF_VECT,idime,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,idime,igaus)
             gprhh(DEF_VECT,idime,igaus)     = gprhs(DEF_VECT,idime,igaus)     + gprhs_sgs(DEF_VECT,idime,igaus)
             elrbq(DEF_VECT,1)               = elrbq(DEF_VECT,1)               + p2vec_bub(DEF_VECT,idime,igaus) * gprhh(DEF_VECT,idime,igaus)
          end do
       end do
       do igaus = 1,pgaus
          elrbq(DEF_VECT,1) = elrbq(DEF_VECT,1) + gprhc(DEF_VECT,igaus) * p2sca_bub(DEF_VECT,igaus) ! ( rhs , q_bub)
       end do
       !
       ! Penalization
       !
       do igaus = 1,pgaus
          elaqq(DEF_VECT,1,1) = elaqq(DEF_VECT,1,1) + penal_nsi * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) 
          elrbq(DEF_VECT,1)   = elrbq(DEF_VECT,1)   + penal_nsi * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * elbub(DEF_VECT)
       end do

    end if

  end subroutine nsi_element_assembly_asgs_oss

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   ASGS and OSS
  !> @details Assembly of Navier Stokes equations using the split OSS
  !>          Variational Multiscale Model.
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_element_assembly_split_oss(&
       pnode,pgaus,gpden,gpvis,gppor,gpsp1,gpsp2,gpvol,   &
       gpsha,gpcar,gpadv,gpvep,gpgrp,gprhs,gprhc,gpvel,   &
       gpsgs,elvel,elpre,elbub,elauu,elaup,elapp,elapu,   &
       elrbu,elrbp,dtinv_loc,dtsgs,pbubl,gpsha_bub,       &
       gpcar_bub,elauq,elapq,elaqu,elaqp,elaqq,elrbq,&
       ! Original global variables
       kfl_lumped,&
       mnode,ntens,&
       kfl_duatss,&
       fact_duatss,&
       kfl_stabi_nsi,&
       fvins_nsi,fcons_nsi,bemol_nsi,kfl_regim_nsi,&
       fvela_nsi,kfl_rmom2_nsi,kfl_press_nsi,&
       kfl_p1ve2_nsi,kfl_linea_nsi,pabdf_nsi,&
       kfl_confi_nsi,nbdfp_nsi,kfl_sgsti_nsi,&
       kfl_nota1_nsi,kfl_limit_nsi,kfl_penal_nsi,&
       penal_nsi,&
       kfl_bubbl_nsi)

    integer(ip), intent(in)    :: kfl_lumped
    integer(ip), intent(in)    :: mnode
    integer(ip), intent(in)    :: ntens
    integer(ip), intent(in)    :: kfl_duatss
    integer(ip), intent(in)    :: fact_duatss
    integer(ip), intent(in)    :: kfl_stabi_nsi
    real(rp),    intent(in)    :: fvins_nsi
    real(rp),    intent(in)    :: fcons_nsi
    real(rp),    intent(in)    :: bemol_nsi
    integer(ip), intent(in)    :: kfl_regim_nsi
    real(rp),    intent(in)    :: fvela_nsi(3)
    integer(ip), intent(in)    :: kfl_rmom2_nsi
    integer(ip), intent(in)    :: kfl_press_nsi
    integer(ip), intent(in)    :: kfl_p1ve2_nsi
    integer(ip), intent(in)    :: kfl_linea_nsi
    real(rp),    intent(in)    :: pabdf_nsi(*)
    integer(ip), intent(in)    :: kfl_confi_nsi
    integer(ip), intent(in)    :: nbdfp_nsi
    integer(ip), intent(in)    :: kfl_sgsti_nsi
    integer(ip), intent(in)    :: kfl_nota1_nsi
    integer(ip), intent(in)    :: kfl_limit_nsi
    integer(ip), intent(in)    :: kfl_penal_nsi
    real(rp),    intent(in)    :: penal_nsi
    integer(ip), intent(in)    :: kfl_bubbl_nsi

    integer(ip), intent(in)    :: pnode,pgaus
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gpgrp(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhc(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvel(VECTOR_SIZE,ndime,pgaus,*)
    real(rp),    intent(in)    :: gpsgs(VECTOR_SIZE,ndime,pgaus,*)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(in)    :: elpre(VECTOR_SIZE,pnode,*)
    real(rp),    intent(in)    :: elbub(VECTOR_SIZE)
    ! Matrices
    real(rp),    intent(out)   :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(out)   :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(out)   :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(out)   :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elrbp(VECTOR_SIZE,pnode)
    ! Others
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    real(rp),    intent(in)    :: dtsgs(VECTOR_SIZE)
    integer(ip), intent(in)    :: pbubl(VECTOR_SIZE)
    real(rp),    intent(in)    :: gpsha_bub(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpcar_bub(VECTOR_SIZE,ndime,pgaus)
    ! Enrichement Element matrices
    real(rp),    intent(out)   :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp),    intent(out)   :: elapq(VECTOR_SIZE,pnode,1)
    real(rp),    intent(out)   :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp),    intent(out)   :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp),    intent(out)   :: elaqq(VECTOR_SIZE,1,1)
    real(rp),    intent(out)   :: elrbq(VECTOR_SIZE,1)
    ! Local arrays
    real(rp)                   :: wgrgr(VECTOR_SIZE,pnode,pnode,pgaus)
    real(rp)                   :: agrau(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: gpsp1_p(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpsp1_v(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpsp2_v(VECTOR_SIZE,pgaus)
    real(rp)                   :: c1(VECTOR_SIZE)
    real(rp)                   :: c2(VECTOR_SIZE)
    real(rp)                   :: c3(VECTOR_SIZE)
    real(rp)                   :: c4(VECTOR_SIZE)
    real(rp)                   :: alpha(VECTOR_SIZE)
    real(rp)                   :: beta(VECTOR_SIZE)
    real(rp)                   :: fact0(VECTOR_SIZE)
    real(rp)                   :: fact1(VECTOR_SIZE)
    real(rp)                   :: fact2(VECTOR_SIZE)
    real(rp)                   :: fact3(VECTOR_SIZE)
    real(rp)                   :: fact4(VECTOR_SIZE)
    real(rp)                   :: fact5(VECTOR_SIZE)
    real(rp)                   :: fact6(VECTOR_SIZE)
    real(rp)                   :: fact7(VECTOR_SIZE)
    real(rp)                   :: fact8(VECTOR_SIZE)
    real(rp)                   :: gpveo(VECTOR_SIZE,3)
    real(rp)                   :: fact1_p(VECTOR_SIZE)
    real(rp)                   :: dtinv_mod(VECTOR_SIZE)
    integer(ip)                :: inode,jnode,jdime
    integer(ip)                :: idofv,jdof2,jdof3,ivect
    integer(ip)                :: idof1,idof3,idof2,igaus
    integer(ip)                :: idime,jdof1,jdofv,itime


#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
    
    dtinv_mod = dtinv_loc

    !----------------------------------------------------------------------
    !
    ! possibility of using only pressure stabilization - not ready with limiter - nor with shock capturing
    !
    !----------------------------------------------------------------------

    gpsp1_p = gpsp1
    gpsp1_v = gpsp1
    gpsp2_v = gpsp2

    if( kfl_nota1_nsi == 1 ) gpsp1_v = 0.0_rp 

    !----------------------------------------------------------------------
    !
    ! Initialization
    !
    !----------------------------------------------------------------------

    elrbp = 0.0_rp
    elrbu = 0.0_rp
    elapp = 0.0_rp
    elauu = 0.0_rp
    elaup = 0.0_rp
    elapu = 0.0_rp

    !----------------------------------------------------------------------
    !
    ! Test functions
    !
    !----------------------------------------------------------------------

    !
    ! AGRAU = rho * (a.grad) Ni
    ! WGRGR = grad(Ni) . grad(Nj)
    !
    agrau(DEF_VECT,:,:)   = 0.0_rp
    wgrgr(DEF_VECT,:,:,:) = 0.0_rp 
    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             agrau(DEF_VECT,inode,igaus) =  agrau(DEF_VECT,inode,igaus) + &
                  &                         gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
          end do
          agrau(DEF_VECT,inode,igaus) =  gpden(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus) 
          do jnode = 1,pnode
             do idime = 1,ndime
                wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + &
                     &                              gpcar(DEF_VECT,idime,inode,igaus)*gpcar(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do
    end do

    !----------------------------------------------------------------------
    !
    ! Auu
    !
    !----------------------------------------------------------------------
    !
    ! Galerkin + ( tau2 * div(u) , div(v) ) + ( tau1 * rho*a.grad(u), rho*a.grad(v) )
    !
    do igaus = 1,pgaus

       fact0(DEF_VECT) = gpsp2_v(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       fact6(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)
       fact7(DEF_VECT) = gpsp1_v(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       fact8(DEF_VECT) = pabdf_nsi(1) * gpden(DEF_VECT,igaus) * dtinv_mod(DEF_VECT) + gppor(DEF_VECT,igaus)

       do inode = 1,pnode
          do idime = 1,ndime

             idofv           = (inode-1)*ndime+idime
             fact1(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)      

             do jnode = 1,pnode    
                !
                ! div(u) * tau2' * div(v)
                !
                do jdime = 1,ndime                   
                   jdofv                       = (jnode-1)*ndime+jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,jdime,jnode,igaus)                   
                end do
                !
                ! ( rho/dt N_j + s Nj + rho*(a.grad)Nj ) Ni 
                ! + mu * grad(Ni) . grad(Nj)
                ! + t1 * rho*(a.grad)Nj * rho*(a.grad)Ni
                !
                jdofv           = (jnode-1)*ndime+idime
                fact4(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                fact5(DEF_VECT) = fact4(DEF_VECT) * ( agrau(DEF_VECT,jnode,igaus) + fact8(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) ) & 
                     &           + fact6(DEF_VECT) *   wgrgr(DEF_VECT,inode,jnode,igaus) &                                               
                     &           + fact7(DEF_VECT) *   agrau(DEF_VECT,jnode,igaus) * agrau(DEF_VECT,inode,igaus)   
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact5(DEF_VECT)

             end do
          end do
       end do
    end do
    !
    ! ( mu*duj/dxi , dv/dxj ) (only div form)
    !
    if( fvins_nsi > 0.9_rp ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                idofv = (inode-1)*ndime + idime
                do jnode = 1,pnode
                   fact1(DEF_VECT) = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)     
                   do jdime = 1,ndime
                      jdofv                       = (jnode-1)*ndime + jdime
                      elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,jdime,inode,igaus)
                   end do
                end do
                if( fvins_nsi == 2.0_rp ) then
                   fact1(DEF_VECT) = -2.0_rp / 3.0_rp * gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                   do jnode = 1,pnode
                      do jdime = 1,ndime
                         jdofv                       = (jnode-1)*ndime + jdime
                         elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,jdime,jnode,igaus)
                      end do
                   end do
                end if
             end do
          end do
       end do
    end if

    !call nsi_element_system_output(&
    !     pnode,elauu(1,:,:),elaup(1,:,:),elapp(1,:,:),elapu(1,:,:),elrbu(1,:,:),elrbp(1,:),&
    !     elauq(1,:,:),elapq(1,:,:),elaqu(1,:,:),elaqp(1,:,:),elaqq(1,:,:),elrbq(1,:))
    !stop
    !
    ! Lumped evolution matrix (only backward euler)
    !
    if( kfl_lumped == 1 ) then 
       !
       ! Remove Galerkin term and add lumped term 
       ! 
       if( ndime == 2 ) then
          stop
       else
          do igaus = 1,pgaus
             gpveo(DEF_VECT,1:3) = 0.0_rp
             do inode = 1,pnode
                do idime = 1,ndime
                   gpveo(DEF_VECT,idime) = gpveo(DEF_VECT,idime) + elvel(DEF_VECT,idime,inode,2) * gpsha(DEF_VECT,inode,igaus)
                end do
             end do
             do inode = 1,pnode
                idof1                       = 3*inode-2
                idof2                       = 3*inode-1
                idof3                       = 3*inode
                fact0(DEF_VECT)             = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * dtinv_mod(DEF_VECT)
                elauu(DEF_VECT,idof1,idof1) = elauu(DEF_VECT,idof1,idof1) + fact0(DEF_VECT)
                elauu(DEF_VECT,idof2,idof2) = elauu(DEF_VECT,idof2,idof2) + fact0(DEF_VECT)
                elauu(DEF_VECT,idof3,idof3) = elauu(DEF_VECT,idof3,idof3) + fact0(DEF_VECT)
                do idime = 1,ndime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) - fact0(DEF_VECT) * gpveo(DEF_VECT,idime)
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact0(DEF_VECT) * elvel(DEF_VECT,idime,inode,2)
                end do
                do jnode = 1,pnode 
                   jdof1                       = 3*jnode-2
                   jdof2                       = 3*jnode-1
                   jdof3                       = 3*jnode
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) 
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) - fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) 
                   elauu(DEF_VECT,idof3,jdof3) = elauu(DEF_VECT,idof3,jdof3) - fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) 
                end do
             end do
          end do
       end if

    else if( kfl_lumped == 2 ) then 
       !
       ! No time term have been added up to now: add Galerkin term
       !
       do igaus = 1,pgaus
          fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * dtinv_mod(DEF_VECT)
          do inode = 1, pnode
             fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             do idime = 1,ndime
                idof1                       = (inode-1) * ndime + idime
                elauu(DEF_VECT,idof1,idof1) = elauu(DEF_VECT,idof1,idof1) + fact1(DEF_VECT)
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact1(DEF_VECT) * elvel(DEF_VECT,idime,inode,2)
             end do
          end do
       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Apu and Aup
    !
    !----------------------------------------------------------------------
    !
    ! ( div(u) , q ) and - ( p , div(v) ) 
    !
    if( ndime == 2 ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             idof1 = 2*inode-1
             idof2 = 2*inode
             do jnode = 1,pnode
                fact0(DEF_VECT)             = gpvol(DEF_VECT,igaus)       * gpsha(DEF_VECT,jnode,igaus) 
                fact1(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,1,inode,igaus)
                fact2(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,2,inode,igaus)
                elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) + fact1(DEF_VECT)
                elapu(DEF_VECT,jnode,idof2) = elapu(DEF_VECT,jnode,idof2) + fact2(DEF_VECT)
                elaup(DEF_VECT,idof1,jnode) = elaup(DEF_VECT,idof1,jnode) - fact1(DEF_VECT)
                elaup(DEF_VECT,idof2,jnode) = elaup(DEF_VECT,idof2,jnode) - fact2(DEF_VECT)
             end do
          end do
       end do
    else
       do igaus = 1,pgaus
          do inode = 1,pnode
             idof1 = 3*inode-2
             idof2 = 3*inode-1
             idof3 = 3*inode
             do jnode = 1,pnode
                fact0(DEF_VECT)             = gpvol(DEF_VECT,igaus)       * gpsha(DEF_VECT,jnode,igaus) 
                fact1(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,1,inode,igaus)
                fact2(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,2,inode,igaus)
                fact3(DEF_VECT)             = fact0(DEF_VECT)             * gpcar(DEF_VECT,3,inode,igaus)
                elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) + fact1(DEF_VECT)
                elapu(DEF_VECT,jnode,idof2) = elapu(DEF_VECT,jnode,idof2) + fact2(DEF_VECT)
                elapu(DEF_VECT,jnode,idof3) = elapu(DEF_VECT,jnode,idof3) + fact3(DEF_VECT) 
                elaup(DEF_VECT,idof1,jnode) = elaup(DEF_VECT,idof1,jnode) - fact1(DEF_VECT)
                elaup(DEF_VECT,idof2,jnode) = elaup(DEF_VECT,idof2,jnode) - fact2(DEF_VECT)
                elaup(DEF_VECT,idof3,jnode) = elaup(DEF_VECT,idof3,jnode) - fact3(DEF_VECT)
             end do
          end do
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! App
    !
    !----------------------------------------------------------------------
    !
    ! Pressure: ( tau1' * grad(p) , grad(q) )
    ! 
    if( kfl_stabi_nsi /= -1 ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             do jnode = inode+1,pnode
                fact1(DEF_VECT)             = gpsp1_p(DEF_VECT,igaus) * wgrgr(DEF_VECT,jnode,inode,igaus) * gpvol(DEF_VECT,igaus)
                elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + fact1(DEF_VECT)
                elapp(DEF_VECT,inode,jnode) = elapp(DEF_VECT,inode,jnode) + fact1(DEF_VECT)
             end do
             fact1(DEF_VECT)             = gpsp1_p(DEF_VECT,igaus) * wgrgr(DEF_VECT,inode,inode,igaus) * gpvol(DEF_VECT,igaus)
             elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact1(DEF_VECT)
          end do
       end do
    end if
    !
    ! Penalization
    !
    do igaus = 1,pgaus
       fact1(DEF_VECT) = penal_nsi * gpvol(DEF_VECT,igaus)
       do inode = 1,pnode
          elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
          elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * elpre(DEF_VECT,inode,1) 
       end do
    end do

    !----------------------------------------------------------------------
    !
    ! bu and bp
    !
    ! P1  = P [ tau1' * rho * a.grad(u) ]
    ! P1' = P1 + tau1' * rho * u'n / dt
    !
    ! P2  = P [ tau1' * ( grad(p) - f ) ]
    ! P2' = P2 + tau1' * rho * u'n / dt + tau1' * f 
    !
    !----------------------------------------------------------------------
    !
    ! Limiter
    !
    if( kfl_stabi_nsi == -1 ) then

       gpvep(DEF_VECT,:,:) = 0.0_rp 

    else if( kfl_limit_nsi == -1 ) then

       gpvep(DEF_VECT,:,:) = 0.0_rp

    else if( kfl_limit_nsi > 0 ) then

       do igaus = 1,pgaus
          c1(DEF_VECT) = 0.0_rp
          c2(DEF_VECT) = 0.0_rp
          c3(DEF_VECT) = 0.0_rp
          do idime = 1,ndime
             c4(DEF_VECT) = 0.0_rp
             do inode = 1,pnode
                c4(DEF_VECT) = c4(DEF_VECT) + agrau(DEF_VECT,inode,igaus) * elvel(DEF_VECT,idime,inode,1)
             end do
             c4(DEF_VECT) = gpsp1(DEF_VECT,igaus) * c4(DEF_VECT)
             c1(DEF_VECT) = c1(DEF_VECT) + ( gpvep(DEF_VECT,idime,igaus) - c4(DEF_VECT) )**2
             c3(DEF_VECT) = c3(DEF_VECT) + gpvep(DEF_VECT,idime,igaus) * gpvep(DEF_VECT,idime,igaus)
             c2(DEF_VECT) = c2(DEF_VECT) + c4(DEF_VECT) * c4(DEF_VECT)
          end do
          c3(DEF_VECT)   = sqrt( c2(DEF_VECT) ) + sqrt( c3(DEF_VECT) )
          c1(DEF_VECT)   = sqrt( c1(DEF_VECT) )
          beta(DEF_VECT) = c1(DEF_VECT) / ( c3(DEF_VECT) + epsilon(1.0_rp) )
          if( kfl_limit_nsi == 1 ) then
             alpha(DEF_VECT) = min(1.0_rp,2.0_rp*(1.0_rp-beta(DEF_VECT)))
          else if( kfl_limit_nsi == 2 ) then
             alpha(DEF_VECT) = 0.5_rp*(tanh(20.0_rp*(beta(DEF_VECT)-0.8_rp))+1.0_rp)
          end if
          do idime = 1,ndime
             gpvep(DEF_VECT,idime,igaus) = alpha(DEF_VECT) * gpvep(DEF_VECT,idime,igaus)
          end do
       end do

    end if
    !
    ! P2 <= P2 + tau1' * f
    !
    if( kfl_stabi_nsi == -1 ) then
       gpgrp(DEF_VECT,:,:) = 0.0_rp
    else
       do igaus = 1,pgaus
          do idime = 1,ndime
             gpgrp(DEF_VECT,idime,igaus) = gpgrp(DEF_VECT,idime,igaus) + gpsp1_p(DEF_VECT,igaus) * gprhs(DEF_VECT,idime,igaus)
          end do
       end do
       !
       ! P1 <= P1 + tau1' * rho * u'n / dt
       ! P2 <= P2 + tau1' * rho * u'n / dt
       !
       if( kfl_sgsti_nsi == 1 ) then
          do igaus = 1,pgaus 
             fact1(DEF_VECT)    = gpden(DEF_VECT,igaus) * dtsgs(DEF_VECT) * gpsp1_v(DEF_VECT,igaus)
             fact1_p (DEF_VECT) = gpden(DEF_VECT,igaus) * dtsgs(DEF_VECT) * gpsp1_p(DEF_VECT,igaus)
             do idime = 1,ndime
                gpvep(DEF_VECT,idime,igaus) = gpvep(DEF_VECT,idime,igaus) + fact1(DEF_VECT)   * gpsgs(DEF_VECT,idime,igaus,2)
                gpgrp(DEF_VECT,idime,igaus) = gpgrp(DEF_VECT,idime,igaus) + fact1_p(DEF_VECT) * gpsgs(DEF_VECT,idime,igaus,2)
             end do
          end do
       end if
    end if
    !
    ! bu = ( f + rho*u^n/dt , v ) + ( rho * a.grad(v) , tau1' * rho u'^n/dt + P1 ) 
    !    = ( f + rho*u^n/dt , v ) + ( rho * a.grad(v) , P1' ) 
    !
    ! bp = ( f + rho*u'^n/dt , tau1' grad(q) ) + ( P2 , grad(q) )
    !    = ( P2' , grad(q) ) 
    ! 
    do igaus = 1,pgaus
       fact4(DEF_VECT) = gpden(DEF_VECT,igaus) * dtinv_mod(DEF_VECT)
       do itime = 2,nbdfp_nsi
          do idime = 1,ndime
             gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - pabdf_nsi(itime) * fact4(DEF_VECT) * gpvel(DEF_VECT,idime,igaus,itime)
          end do
       end do
       do inode = 1,pnode
          fact1(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f + rho*u^n/dt , v )
          fact3(DEF_VECT) = gpvol(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus)  ! ( rho * a.grad(v) , P1' ) 
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact1(DEF_VECT) * gprhs(DEF_VECT,idime,igaus) &
                  &                                                    + fact3(DEF_VECT) * gpvep(DEF_VECT,idime,igaus)              
          end do
          elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) + gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * gprhc(DEF_VECT,igaus)                ! ( rhs, q )
          do idime = 1,ndime
             elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) + gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,inode,igaus) * gpgrp(DEF_VECT,idime,igaus) ! ( P2' , grad(q) ) 
          end do
       end do
    end do

    !--------------------------------------------------------------------
    !
    ! Pressure bubble
    !
    !--------------------------------------------------------------------

    if( maxval(pbubl) == 1 ) then
       if( kfl_stabi_nsi /= -1 ) then
          write(6,*) 'BUBBLE NOT CODED FOR SPLIT OSS'
          stop
       end if
       !
       ! Initialization
       !
       elauq = 0.0_rp
       elapq = 0.0_rp
       elaqu = 0.0_rp
       elaqp = 0.0_rp
       elaqq = 0.0_rp
       elrbq = 0.0_rp
       !
       ! Auq and Aqu
       !
       if( kfl_press_nsi == 1 ) then
          do igaus = 1,pgaus
             fact1(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
             do inode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) - fact1(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
                   elaqu(DEF_VECT,1,idofv) = elaqu(DEF_VECT,1,idofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) 
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             fact1(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
             do inode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                   elaqu(DEF_VECT,1,idofv) = elaqu(DEF_VECT,1,idofv) + fact1(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) 
                end do
             end do
          end do
       end if
       !
       ! Penalization and others
       !
       do igaus = 1,pgaus
          elaqq(DEF_VECT,1,1) = elaqq(DEF_VECT,1,1) + gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * penal_nsi
          elrbq(DEF_VECT,1)   = elrbq(DEF_VECT,1)   + gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * penal_nsi * elbub(DEF_VECT) 
          elrbq(DEF_VECT,1)   = elrbq(DEF_VECT,1)   + gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * gprhc(DEF_VECT,igaus) 
       end do

    end if

  end subroutine nsi_element_assembly_split_oss

end module mod_nsi_element_assembly
!> @}
