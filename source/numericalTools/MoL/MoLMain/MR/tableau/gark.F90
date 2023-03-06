!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
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

      gamTau = 0.0
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

      wTau = 0.0
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
      cS(1) = 0.0
      cS(2) = 0.4358665215084589994160194511935568425
      cS(3) = 0.4358665215084589994160194511935568425
      cS(4) = 0.7179332607542294997080097255967784213
      cS(5) = 0.7179332607542294997080097255967784213
      cS(6) = 1.0
      cS(7) = 1.0
      cS(8) = 1.0

      allocate (gamK(8, 8, 1))
      allocate (wK(8, 8, 1))

      gamK = 0.0
      gamK(2, 1, 1) = 0.4358665215084589994160194511935568425

      gamK(3, 1, 1) = -0.4358665215084589994160194511935568425
      gamK(3, 3, 1) = 0.4358665215084589994160194511935568425

      gamK(4, 1, 1) = 0.0414273753564414837153799230278275639
      gamK(4, 3, 1) = 0.2406393638893290165766103513753940148

      gamK(5, 1, 1) = -0.0414273753564414837153799230278275639
      gamK(5, 3, 1) = -0.3944391461520175157006395281657292786
      gamK(5, 5, 1) = 0.4358665215084589994160194511935568425

      gamK(6, 1, 1) = 0.1123373143006047802633543416889605123
      gamK(6, 3, 1) = 0.1051807513648115027700693049638099167e1
      gamK(6, 5, 1) = -0.8820780887029493076720571169238381009

      gamK(7, 1, 1) = -0.1123373143006047802633543416889605123
      gamK(7, 3, 1) = -0.1253776037178754576562056399779976346
      gamK(7, 5, 1) = -0.1981516034899787614964594695265986957
      gamK(7, 7, 1) = 0.4358665215084589994160194511935568425

      wK = 0.0
      wK(2, 1, 1) = 0.4358665215084589994160194511935568425

      wK(4, 1, 1) = -0.1750145285570467590610670000018749059
      wK(4, 3, 1) = 0.4570812678028172593530572744050964846

      wK(5, 1, 1) = 0.6042689307721552209333459437020635774e-01
      wK(5, 3, 1) = -0.6042689307721552209333459437020635774e-01

      wK(6, 1, 1) = 0.1195213959425454440038786034027936869
      wK(6, 3, 1) = -0.1843725226689661917898533950296297650e1
      wK(6, 5, 1) = 0.2006270569992886974186645621296725542e1

      wK(7, 1, 1) = -0.5466585780430528451745431084418669343
      wK(7, 3, 1) = 0.2000000000000000000000000000000000000e1
      wK(7, 5, 1) = -0.1453341421956947154825456891558133066e1

      wK(8, 1, 1) = 0.1058582960718796387223774594771849530
      wK(8, 3, 1) = 0.6555675011400702509752889543247306350
      wK(8, 5, 1) = -0.1197292318720408889113685864995472431e1
      wK(8, 7, 1) = 0.4358665215084589994160194511935568425
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
      cS(1) = 0.0
      cS(2:3) = 1.0/2.0
      cS(4:5) = 5.0/8.0
      cS(6:7) = 3.0/4.0
      cS(8:9) = 7.0/8.0
      cS(10:12) = 1.0

      allocate (gamK(12, 12, 2))
      allocate (wK(12, 12, 2))

      gamK = 0.0
      gamK(2, 1, 1) = 0.5

      gamK(3, 1, 1) = -0.25
      gamK(3, 3, 1) = 0.25

      gamK(4, 1, 1) = -3.97728124810848818306703385146227889
      gamK(4, 3, 1) = 4.10228124810848818306703385146227889

      gamK(5, 1, 1) = -0.0690538874140169123272414708480937406
      gamK(5, 3, 1) = -0.180946112585983087672758529151906259
      gamK(5, 5, 1) = 0.25

      gamK(6, 1, 1) = -1.76176766375792052886337896482241241
      gamK(6, 3, 1) = 2.69452469837729861015533815079146138
      gamK(6, 5, 1) = -0.807757034619378081291959185969048978

      gamK(7, 1, 1) = 0.555872179155396948730508100958808496
      gamK(7, 3, 1) = -0.679914050157999501395850152788348695
      gamK(7, 5, 1) = -0.125958128997397447334657948170459801
      gamK(7, 7, 1) = 0.25

      gamk(8, 1, 1) = -5.84017602872495595444642665754106511
      gamk(8, 3, 1) = 8.17445668429191508919127080571071637
      gamk(8, 5, 1) = 0.125958128997397447334657948170459801
      gamk(8, 7, 1) = -2.33523878456435658207950209634011106

      gamK(9, 1, 1) = -1.9067926451678118080947593050360523
      gamK(9, 3, 1) = -1.54705781138512393363298457924938844
      gamK(9, 5, 1) = 4.12988801314935030595449173802031322
      gamK(9, 7, 1) = -0.926037556596414564226747853734872477
      gamK(9, 9, 1) = 0.25

      gamK(10, 1, 1) = 3.33702815168872605455765278252966252
      gamK(10, 3, 1) = 1.54705781138512393363298457924938844
      gamK(10, 5, 1) = -4.12988801314935030595449173802031322
      gamK(10, 7, 1) = 0.926037556596414564226747853734872477
      gamK(10, 9, 1) = -1.55523550652091424646289347749361021

      gamK(11, 1, 1) = -0.821293629221007618720524112312446752
      gamK(11, 3, 1) = 0.328610356068599988551677264268969646
      gamK(11, 5, 1) = 0.678001812102026694142641232421139516
      gamK(11, 7, 1) = -0.342779287862800022896645471462060708
      gamK(11, 9, 1) = -0.0925392510868190410771489129156017025
      gamK(11, 11, 1) = 0.25

      gamK(4, 1, 2) = 8.70456249621697636613406770292455778
      gamK(4, 3, 2) = -8.70456249621697636613406770292455778

      gamK(6, 1, 2) = 3.91164310234387488238124087134101229
      gamK(6, 3, 2) = -5.02715717158263104496515924327911025
      gamK(6, 5, 2) = 1.11551406923875616258391837193809796

      gamK(8, 1, 2) = 10.8186076991391180114318371131645132
      gamK(8, 3, 2) = -14.9890852682678311755908413058447354
      gamK(8, 7, 2) = 4.17047756912871316415900419268022213

      gamK(10, 1, 2) = -2.61047101304182849292578695498722043
      gamK(10, 9, 2) = 2.61047101304182849292578695498722043

      wK = 0.0
      wK(2, 1, 1) = 0.5

      wK(4, 1, 1) = -1.91716534363662868878172216064946905
      wK(4, 3, 1) = 2.04216534363662868878172216064946905

      wK(5, 1, 1) = -0.404751031801105942697915907046990469
      wK(5, 3, 1) = 0.404751031801105942697915907046990469

      wK(6, 1, 1) = 11.4514660224922163666569802860263173
      wK(6, 3, 1) = -30.2107574752650427144064781557395061
      wK(6, 5, 1) = 18.8842914527728263477494978697131888

      wK(7, 1, 1) = -0.709033564760261450684711672946330144
      wK(7, 3, 1) = 1.03030720858751876652616190884004718
      wK(7, 5, 1) = -0.321273643827257315841450235893717036

      wK(8, 1, 1) = -29.9954871645582843984091068494419927
      wK(8, 3, 1) = 37.605982774991801805364896856243857
      wK(8, 5, 1) = 0.321273643827257315841450235893717036
      wK(8, 7, 1) = -7.80676925426077472279724024269558129

      wK(9, 1, 1) = 3.10466505427296211633876939184912422
      wK(9, 3, 1) = -2.43032501975716229713206592741556636
      wK(9, 5, 1) = -1.90547930115152463521920165948384213
      wK(9, 7, 1) = 1.23113926663572481601249819505028427

      wK(10, 1, 1) = -2.42442954775204786987587591435551401
      wK(10, 3, 1) = 2.43032501975716229713206592741556636
      wK(10, 5, 1) = 1.90547930115152463521920165948384213
      wK(10, 7, 1) = -1.23113926663572481601249819505028427
      wK(10, 9, 1) = -0.555235506520914246462893477493610215

      wK(11, 1, 1) = -0.010441350444797485902945189451653542
      wK(11, 3, 1) = 0.0726030361465507450515210450548814161
      wK(11, 5, 1) = -0.128827595167726095223945409857642431
      wK(11, 7, 1) = 0.112935535009382356613944010712215408
      wK(11, 9, 1) = -0.0462696255434095205385744564578008512

      wK(12, 1, 1) = -0.81085227877621013281757892286079321
      wK(12, 3, 1) = 0.25600731992204924350015621921408823
      wK(12, 5, 1) = 0.806829407269752789366586642278781947
      wK(12, 7, 1) = -0.455714822872182379510589482174276116
      wK(12, 9, 1) = -0.0462696255434095205385744564578008512
      wK(12, 11, 1) = 0.25

      wK(4, 1, 2) = 4.0843306872732573775634443212989381
      wK(4, 3, 2) = -4.0843306872732573775634443212989381

      wK(6, 1, 2) = -21.8434299813822208479181287579586536
      wK(6, 3, 2) = 59.6120128869278735434171244973850312
      wK(6, 5, 2) = -37.7685829055456526954989957394263776

      wK(8, 1, 2) = 61.6590414586370916981876370447766458
      wK(8, 3, 2) = -77.2725799671586411437821175301678084
      wK(8, 7, 2) = 15.6135385085215494455944804853911626

      wK(10, 1, 2) = -1.11047101304182849292578695498722043
      wK(10, 9, 2) = 1.11047101304182849292578695498722043
   end subroutine gark4_init

end module gark
