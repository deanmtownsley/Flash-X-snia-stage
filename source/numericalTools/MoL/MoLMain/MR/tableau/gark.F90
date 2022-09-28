!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief GARK tableau utilities

!> @ingroup MoLMR
!! @brief Utilities for setting up a specified GARK tableau.
!!
!! @details
!!    Available methods currently include (list by runtime parameter values for slowMethod):
!!       - `mr-gark3` : Third-order IMEX-MRI-GARK3b scheme in [1]
!!       - `mr-gark4` : Fourth-order IMEX-MRI-GARK4 scheme in [1]
!!
!!    The tableau are all given in the form:
!!
!!    @f[
!!       \begin{array}{c|ccc}
!!          c_1    & A_{11} & \dots & A_{1n}\\
!!          \vdots & \vdots & \ddots& \vdots\\
!!          c_n    & A_{n1} & \dots & A_{nn}\\
!!           1     & b_{1}  & \dots & b_{n} \\
!!          \hline
!!                 & b_1    & \dots & b_n\\
!!       \end{array}
!!    @f]
!!
!!    - For the implicit methods, the matrices a_ij are lower-triangular
!!    - For the explicit methods, the matrices a_ij are strictly lower-triangular
!!
!!    All tableau follow a "solve decoupled" formulation that allows for alternating
!!    slow- and fast-stages
!!
!!    @par References
!!    @parblock
!!    [1] Rujeko Chinomona and Daniel R. Reynolds,
!!        Implicit-Explicit Multirate Infinitesimal GARK Methods,
!!        SIAM Journal on Scientific Computing 2021 43:5, A3082-A3113,
!!        https://doi.org/10.1137/20M1354349
!!    @endparblock
module gark

   implicit none

contains

   !> Implicit-tableau time interpolant
   !!
   !! @param i    Stage index
   !! @param j    Intermediate state index
   !! @param gamK Implicit tableaus
   !! @param kmax Number of tableaus
   function gamTau(i, j, tau, gamK, kmax)
      implicit none

      real :: gamTau
      integer, intent(in) :: i, j
      real, intent(in) :: tau
      real, dimension(:, :, :), intent(in) :: gamK
      integer, intent(in) :: kmax

      integer :: k

      gamTau = 0d0
      do k = 1, kmax
         gamTau = gamTau + gamK(i, j, k)*tau**(k - 1)
      end do
   end function gamTau

   !> Explicit-tableau time interpolant
   !!
   !! @param i    Stage index
   !! @param j    Intermediate state index
   !! @param wK   Explicit tableaus
   !! @param kmax Number of tableaus
   function wTau(i, j, tau, wK, kmax)

      implicit none

      real :: wTau
      integer, intent(in) :: i, j
      real, intent(in) :: tau
      real, dimension(:, :, :), intent(in) :: wK
      integer, intent(in) :: kmax

      integer :: k

      wTau = 0d0
      do k = 1, kmax
         wTau = wTau + wK(i, j, k)*tau**(k - 1)
      end do
   end function wTau

   !> Third-order IMEX-MRI-GARK3b tableau
   !!
   !! @param kmax   Number of tableaus (for each type)
   !! @param gamK   Implicit tableaus
   !! @param wK     Explicit tableaus
   !! @param cS     Timing coefficients
   !! @param order  Order of the method
   !! @param stages Number of stages
   subroutine gark3_init(kmax, gamK, wK, cS, order, stages)
      implicit none

      integer, intent(out) :: kmax
      real, allocatable, dimension(:, :, :), intent(out) :: gamK, wK
      real, allocatable, dimension(:), intent(out) :: cS
      integer, intent(out) :: order, stages

      integer :: i, j, l, ii, jj

      stages = 8
      order = 3

      kmax = 1

      allocate (cS(8))
      cS(1) = 0d0
      cS(2) = 0.4358665215084589994160194511935568425d0
      cS(3) = 0.4358665215084589994160194511935568425d0
      cS(4) = 0.7179332607542294997080097255967784213d0
      cS(5) = 0.7179332607542294997080097255967784213d0
      cS(6) = 1d0
      cS(7) = 1d0
      cS(8) = 1d0

      allocate (gamK(8, 8, 1))
      allocate (wK(8, 8, 1))

      gamK = 0d0
      gamK(2, 1, 1) = 0.4358665215084589994160194511935568425d0

      gamK(3, 1, 1) = -0.4358665215084589994160194511935568425d0
      gamK(3, 3, 1) = 0.4358665215084589994160194511935568425d0

      gamK(4, 1, 1) = 0.0414273753564414837153799230278275639d0
      gamK(4, 3, 1) = 0.2406393638893290165766103513753940148d0

      gamK(5, 1, 1) = -0.0414273753564414837153799230278275639d0
      gamK(5, 3, 1) = -0.3944391461520175157006395281657292786d0
      gamK(5, 5, 1) = 0.4358665215084589994160194511935568425d0

      gamK(6, 1, 1) = 0.1123373143006047802633543416889605123d0
      gamK(6, 3, 1) = 0.1051807513648115027700693049638099167d1
      gamK(6, 5, 1) = -0.8820780887029493076720571169238381009d0

      gamK(7, 1, 1) = -0.1123373143006047802633543416889605123d0
      gamK(7, 3, 1) = -0.1253776037178754576562056399779976346d0
      gamK(7, 5, 1) = -0.1981516034899787614964594695265986957d0
      gamK(7, 7, 1) = 0.4358665215084589994160194511935568425d0

      wK = 0d0
      wK(2, 1, 1) = 0.4358665215084589994160194511935568425d0

      wK(4, 1, 1) = -0.1750145285570467590610670000018749059d0
      wK(4, 3, 1) = 0.4570812678028172593530572744050964846d0

      wK(5, 1, 1) = 0.6042689307721552209333459437020635774d-01
      wK(5, 3, 1) = -0.6042689307721552209333459437020635774d-01

      wK(6, 1, 1) = 0.1195213959425454440038786034027936869d0
      wK(6, 3, 1) = -0.1843725226689661917898533950296297650d1
      wK(6, 5, 1) = 0.2006270569992886974186645621296725542d1

      wK(7, 1, 1) = -0.5466585780430528451745431084418669343d0
      wK(7, 3, 1) = 0.2000000000000000000000000000000000000d1
      wK(7, 5, 1) = -0.1453341421956947154825456891558133066d1

      wK(8, 1, 1) = 0.1058582960718796387223774594771849530d0
      wK(8, 3, 1) = 0.6555675011400702509752889543247306350d0
      wK(8, 5, 1) = -0.1197292318720408889113685864995472431d1
      wK(8, 7, 1) = 0.4358665215084589994160194511935568425d0
   end subroutine gark3_init

   !> Fourth-order IMEX-MRI-GARK4 tableau
   !! @copydetails gark3_init
   subroutine gark4_init(kmax, gamK, wK, cS, order, stages)
      implicit none

      integer, intent(out) :: kmax
      real, allocatable, dimension(:, :, :), intent(out) :: gamK, wK
      real, allocatable, dimension(:), intent(out) :: cS
      integer, intent(out) :: order, stages

      integer :: i, j, l, ii, jj

      stages = 12
      order = 4

      kmax = 2

      allocate (cS(12))
      cS(1) = 0d0
      cS(2:3) = 1d0/2d0
      cS(4:5) = 5d0/8d0
      cS(6:7) = 3d0/4d0
      cS(8:9) = 7d0/8d0
      cS(10:12) = 1d0

      allocate (gamK(12, 12, 2))
      allocate (wK(12, 12, 2))

      gamK = 0d0
      gamK(2, 1, 1) = 0.5d0

      gamK(3, 1, 1) = -0.25d0
      gamK(3, 3, 1) = 0.25d0

      gamK(4, 1, 1) = -3.97728124810848818306703385146227889d0
      gamK(4, 3, 1) = 4.10228124810848818306703385146227889d0

      gamK(5, 1, 1) = -0.0690538874140169123272414708480937406d0
      gamK(5, 3, 1) = -0.180946112585983087672758529151906259d0
      gamK(5, 5, 1) = 0.25d0

      gamK(6, 1, 1) = -1.76176766375792052886337896482241241d0
      gamK(6, 3, 1) = 2.69452469837729861015533815079146138d0
      gamK(6, 5, 1) = -0.807757034619378081291959185969048978d0

      gamK(7, 1, 1) = 0.555872179155396948730508100958808496d0
      gamK(7, 3, 1) = -0.679914050157999501395850152788348695d0
      gamK(7, 5, 1) = -0.125958128997397447334657948170459801d0
      gamK(7, 7, 1) = 0.25d0

      gamk(8, 1, 1) = -5.84017602872495595444642665754106511d0
      gamk(8, 3, 1) = 8.17445668429191508919127080571071637d0
      gamk(8, 5, 1) = 0.125958128997397447334657948170459801d0
      gamk(8, 7, 1) = -2.33523878456435658207950209634011106d0

      gamK(9, 1, 1) = -1.9067926451678118080947593050360523d0
      gamK(9, 3, 1) = -1.54705781138512393363298457924938844d0
      gamK(9, 5, 1) = 4.12988801314935030595449173802031322d0
      gamK(9, 7, 1) = -0.926037556596414564226747853734872477d0
      gamK(9, 9, 1) = 0.25d0

      gamK(10, 1, 1) = 3.33702815168872605455765278252966252d0
      gamK(10, 3, 1) = 1.54705781138512393363298457924938844d0
      gamK(10, 5, 1) = -4.12988801314935030595449173802031322d0
      gamK(10, 7, 1) = 0.926037556596414564226747853734872477d0
      gamK(10, 9, 1) = -1.55523550652091424646289347749361021d0

      gamK(11, 1, 1) = -0.821293629221007618720524112312446752d0
      gamK(11, 3, 1) = 0.328610356068599988551677264268969646d0
      gamK(11, 5, 1) = 0.678001812102026694142641232421139516d0
      gamK(11, 7, 1) = -0.342779287862800022896645471462060708d0
      gamK(11, 9, 1) = -0.0925392510868190410771489129156017025d0
      gamK(11, 11, 1) = 0.25d0

      gamK(4, 1, 2) = 8.70456249621697636613406770292455778d0
      gamK(4, 3, 2) = -8.70456249621697636613406770292455778d0

      gamK(6, 1, 2) = 3.91164310234387488238124087134101229d0
      gamK(6, 3, 2) = -5.02715717158263104496515924327911025d0
      gamK(6, 5, 2) = 1.11551406923875616258391837193809796d0

      gamK(8, 1, 2) = 10.8186076991391180114318371131645132d0
      gamK(8, 3, 2) = -14.9890852682678311755908413058447354d0
      gamK(8, 7, 2) = 4.17047756912871316415900419268022213d0

      gamK(10, 1, 2) = -2.61047101304182849292578695498722043d0
      gamK(10, 9, 2) = 2.61047101304182849292578695498722043d0

      wK = 0d0
      wK(2, 1, 1) = 0.5d0

      wK(4, 1, 1) = -1.91716534363662868878172216064946905d0
      wK(4, 3, 1) = 2.04216534363662868878172216064946905d0

      wK(5, 1, 1) = -0.404751031801105942697915907046990469d0
      wK(5, 3, 1) = 0.404751031801105942697915907046990469d0

      wK(6, 1, 1) = 11.4514660224922163666569802860263173d0
      wK(6, 3, 1) = -30.2107574752650427144064781557395061d0
      wK(6, 5, 1) = 18.8842914527728263477494978697131888d0

      wK(7, 1, 1) = -0.709033564760261450684711672946330144d0
      wK(7, 3, 1) = 1.03030720858751876652616190884004718d0
      wK(7, 5, 1) = -0.321273643827257315841450235893717036d0

      wK(8, 1, 1) = -29.9954871645582843984091068494419927d0
      wK(8, 3, 1) = 37.605982774991801805364896856243857d0
      wK(8, 5, 1) = 0.321273643827257315841450235893717036d0
      wK(8, 7, 1) = -7.80676925426077472279724024269558129d0

      wK(9, 1, 1) = 3.10466505427296211633876939184912422d0
      wK(9, 3, 1) = -2.43032501975716229713206592741556636d0
      wK(9, 5, 1) = -1.90547930115152463521920165948384213d0
      wK(9, 7, 1) = 1.23113926663572481601249819505028427d0

      wK(10, 1, 1) = -2.42442954775204786987587591435551401d0
      wK(10, 3, 1) = 2.43032501975716229713206592741556636d0
      wK(10, 5, 1) = 1.90547930115152463521920165948384213d0
      wK(10, 7, 1) = -1.23113926663572481601249819505028427d0
      wK(10, 9, 1) = -0.555235506520914246462893477493610215d0

      wK(11, 1, 1) = -0.010441350444797485902945189451653542d0
      wK(11, 3, 1) = 0.0726030361465507450515210450548814161d0
      wK(11, 5, 1) = -0.128827595167726095223945409857642431d0
      wK(11, 7, 1) = 0.112935535009382356613944010712215408d0
      wK(11, 9, 1) = -0.0462696255434095205385744564578008512d0

      wK(12, 1, 1) = -0.81085227877621013281757892286079321d0
      wK(12, 3, 1) = 0.25600731992204924350015621921408823d0
      wK(12, 5, 1) = 0.806829407269752789366586642278781947d0
      wK(12, 7, 1) = -0.455714822872182379510589482174276116d0
      wK(12, 9, 1) = -0.0462696255434095205385744564578008512d0
      wK(12, 11, 1) = 0.25d0

      wK(4, 1, 2) = 4.0843306872732573775634443212989381d0
      wK(4, 3, 2) = -4.0843306872732573775634443212989381d0

      wK(6, 1, 2) = -21.8434299813822208479181287579586536d0
      wK(6, 3, 2) = 59.6120128869278735434171244973850312d0
      wK(6, 5, 2) = -37.7685829055456526954989957394263776d0

      wK(8, 1, 2) = 61.6590414586370916981876370447766458d0
      wK(8, 3, 2) = -77.2725799671586411437821175301678084d0
      wK(8, 7, 2) = 15.6135385085215494455944804853911626d0

      wK(10, 1, 2) = -1.11047101304182849292578695498722043d0
      wK(10, 9, 2) = 1.11047101304182849292578695498722043d0
   end subroutine gark4_init

end module gark
