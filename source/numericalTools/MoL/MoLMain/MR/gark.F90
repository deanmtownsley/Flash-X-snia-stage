!!****if* source/numericalTools/MoL/MoLMain/MR/gark
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!  NAME
!!    gark
!!
!!  SYNOPSIS
!!    use gark
!!
!!  DESCRIPTION
!!    Utilities for setting up a specified GARK tableau.
!!
!!    Available methods currently include (list by runtime parameter values for mr_slowMethod):
!!       mr-gark3 : Third-order IMEX-MRI-GARK3b scheme in [1]
!!       mr-gark4 : Fourth-order IMEX-MRI-GARK4 scheme in [1]
!!
!!    The tableau are all given in the form:
!!
!!       c_1 | a_11 ... a_1n
!!        .  |  .   .    .
!!        .  |  .    .   .
!!        .  |  .     .  .
!!       c_n | a_n1 ... a_nn
!!        1  | b_1  ... b_n
!!       -------------------
!!           | b_1  ... b_n
!!
!!    For the implicit methods, the matrices a_ij are lower-triangular
!!    For the explicit methods, the matrices a_ij are strictly lower-triangular
!!
!!    All tableau follow a "solve decoupled" formulation that allows for alternating
!!    slow- and fast-stages
!!
!! REFERENCES
!!    [1] Implicit-Explicit Multirate Infinitesimal GARK Methods
!!        Rujeko Chinomona and Daniel R. Reynolds
!!        SIAM Journal on Scientific Computing 2021 43:5, A3082-A3113
!!        https://doi.org/10.1137/20M1354349
!!
!!  NOTES
!!
!!***
module gark

   implicit none

contains

   function mr_gamTau(i, j, tau)
      use mr_data, only: mr_gamK, mr_kmax
      implicit none

      real :: mr_gamTau
      integer, intent(in) :: i, j
      real, intent(in) :: tau

      integer :: k

      mr_gamTau = 0d0
      do k = 1, mr_kmax
         mr_gamTau = mr_gamTau + mr_gamK(i, j, k)*tau**(k - 1)
      end do
   end function mr_gamTau

   function mr_wTau(i, j, tau)
      use mr_data, only: mr_wK, mr_kmax

      implicit none

      real :: mr_wTau
      integer, intent(in) :: i, j
      real, intent(in) :: tau

      integer :: k

      mr_wTau = 0d0
      do k = 1, mr_kmax
         mr_wTau = mr_wTau + mr_wK(i, j, k)*tau**(k - 1)
      end do
   end function mr_wTau

   subroutine gark3_init()
      use mr_data, only: mr_nstages_slow, mr_kmax, mr_cS, mr_gamK, mr_wK

      implicit none

      integer :: i, j, l, ii, jj

      mr_nstages_slow = 8

      mr_kmax = 1

      allocate (mr_cS(8))
      mr_cS(1) = 0d0
      mr_cS(2) = 0.4358665215084589994160194511935568425d0
      mr_cS(3) = 0.4358665215084589994160194511935568425d0
      mr_cS(4) = 0.7179332607542294997080097255967784213d0
      mr_cS(5) = 0.7179332607542294997080097255967784213d0
      mr_cS(6) = 1d0
      mr_cS(7) = 1d0
      mr_cS(8) = 1d0

      allocate (mr_gamK(8, 8, 1))
      allocate (mr_wK(8, 8, 1))

      mr_gamK = 0d0
      mr_gamK(2, 1, 1) = 0.4358665215084589994160194511935568425d0

      mr_gamK(3, 1, 1) = -0.4358665215084589994160194511935568425d0
      mr_gamK(3, 3, 1) = 0.4358665215084589994160194511935568425d0

      mr_gamK(4, 1, 1) = 0.0414273753564414837153799230278275639d0
      mr_gamK(4, 3, 1) = 0.2406393638893290165766103513753940148d0

      mr_gamK(5, 1, 1) = -0.0414273753564414837153799230278275639d0
      mr_gamK(5, 3, 1) = -0.3944391461520175157006395281657292786d0
      mr_gamK(5, 5, 1) = 0.4358665215084589994160194511935568425d0

      mr_gamK(6, 1, 1) = 0.1123373143006047802633543416889605123d0
      mr_gamK(6, 3, 1) = 0.1051807513648115027700693049638099167d1
      mr_gamK(6, 5, 1) = -0.8820780887029493076720571169238381009d0

      mr_gamK(7, 1, 1) = -0.1123373143006047802633543416889605123d0
      mr_gamK(7, 3, 1) = -0.1253776037178754576562056399779976346d0
      mr_gamK(7, 5, 1) = -0.1981516034899787614964594695265986957d0
      mr_gamK(7, 7, 1) = 0.4358665215084589994160194511935568425d0

      mr_wK = 0d0
      mr_wK(2, 1, 1) = 0.4358665215084589994160194511935568425d0

      mr_wK(4, 1, 1) = -0.1750145285570467590610670000018749059d0
      mr_wK(4, 3, 1) = 0.4570812678028172593530572744050964846d0

      mr_wK(5, 1, 1) = 0.6042689307721552209333459437020635774d-01
      mr_wK(5, 3, 1) = -0.6042689307721552209333459437020635774d-01

      mr_wK(6, 1, 1) = 0.1195213959425454440038786034027936869d0
      mr_wK(6, 3, 1) = -0.1843725226689661917898533950296297650d1
      mr_wK(6, 5, 1) = 0.2006270569992886974186645621296725542d1

      mr_wK(7, 1, 1) = -0.5466585780430528451745431084418669343d0
      mr_wK(7, 3, 1) = 0.2000000000000000000000000000000000000d1
      mr_wK(7, 5, 1) = -0.1453341421956947154825456891558133066d1

      mr_wK(8, 1, 1) = 0.1058582960718796387223774594771849530d0
      mr_wK(8, 3, 1) = 0.6555675011400702509752889543247306350d0
      mr_wK(8, 5, 1) = -0.1197292318720408889113685864995472431d1
      mr_wK(8, 7, 1) = 0.4358665215084589994160194511935568425d0
   end subroutine gark3_init

   subroutine gark4_init()
      use mr_data, only: mr_nstages_slow, mr_kmax, mr_cS, mr_gamK, mr_wK

      implicit none

      integer :: i, j, l, ii, jj

      mr_nstages_slow = 12

      mr_kmax = 2

      allocate (mr_cS(12))
      mr_cS(1) = 0d0
      mr_cS(2:3) = 1d0/2d0
      mr_cS(4:5) = 5d0/8d0
      mr_cS(6:7) = 3d0/4d0
      mr_cS(8:9) = 7d0/8d0
      mr_cS(10:12) = 1d0

      allocate (mr_gamK(12, 12, 2))
      allocate (mr_wK(12, 12, 2))

      mr_gamK = 0d0
      mr_gamK(2, 1, 1) = 0.5d0

      mr_gamK(3, 1, 1) = -0.25d0
      mr_gamK(3, 3, 1) = 0.25d0

      mr_gamK(4, 1, 1) = -3.97728124810848818306703385146227889d0
      mr_gamK(4, 3, 1) = 4.10228124810848818306703385146227889d0

      mr_gamK(5, 1, 1) = -0.0690538874140169123272414708480937406d0
      mr_gamK(5, 3, 1) = -0.180946112585983087672758529151906259d0
      mr_gamK(5, 5, 1) = 0.25d0

      mr_gamK(6, 1, 1) = -1.76176766375792052886337896482241241d0
      mr_gamK(6, 3, 1) = 2.69452469837729861015533815079146138d0
      mr_gamK(6, 5, 1) = -0.807757034619378081291959185969048978d0

      mr_gamK(7, 1, 1) = 0.555872179155396948730508100958808496d0
      mr_gamK(7, 3, 1) = -0.679914050157999501395850152788348695d0
      mr_gamK(7, 5, 1) = -0.125958128997397447334657948170459801d0
      mr_gamK(7, 7, 1) = 0.25d0

      mr_gamk(8, 1, 1) = -5.84017602872495595444642665754106511d0
      mr_gamk(8, 3, 1) = 8.17445668429191508919127080571071637d0
      mr_gamk(8, 5, 1) = 0.125958128997397447334657948170459801d0
      mr_gamk(8, 7, 1) = -2.33523878456435658207950209634011106d0

      mr_gamK(9, 1, 1) = -1.9067926451678118080947593050360523d0
      mr_gamK(9, 3, 1) = -1.54705781138512393363298457924938844d0
      mr_gamK(9, 5, 1) = 4.12988801314935030595449173802031322d0
      mr_gamK(9, 7, 1) = -0.926037556596414564226747853734872477d0
      mr_gamK(9, 9, 1) = 0.25d0

      mr_gamK(10, 1, 1) = 3.33702815168872605455765278252966252d0
      mr_gamK(10, 3, 1) = 1.54705781138512393363298457924938844d0
      mr_gamK(10, 5, 1) = -4.12988801314935030595449173802031322d0
      mr_gamK(10, 7, 1) = 0.926037556596414564226747853734872477d0
      mr_gamK(10, 9, 1) = -1.55523550652091424646289347749361021d0

      mr_gamK(11, 1, 1) = -0.821293629221007618720524112312446752d0
      mr_gamK(11, 3, 1) = 0.328610356068599988551677264268969646d0
      mr_gamK(11, 5, 1) = 0.678001812102026694142641232421139516d0
      mr_gamK(11, 7, 1) = -0.342779287862800022896645471462060708d0
      mr_gamK(11, 9, 1) = -0.0925392510868190410771489129156017025d0
      mr_gamK(11, 11, 1) = 0.25d0

      mr_gamK(4, 1, 2) = 8.70456249621697636613406770292455778d0
      mr_gamK(4, 3, 2) = -8.70456249621697636613406770292455778d0

      mr_gamK(6, 1, 2) = 3.91164310234387488238124087134101229d0
      mr_gamK(6, 3, 2) = -5.02715717158263104496515924327911025d0
      mr_gamK(6, 5, 2) = 1.11551406923875616258391837193809796d0

      mr_gamK(8, 1, 2) = 10.8186076991391180114318371131645132d0
      mr_gamK(8, 3, 2) = -14.9890852682678311755908413058447354d0
      mr_gamK(8, 7, 2) = 4.17047756912871316415900419268022213d0

      mr_gamK(10, 1, 2) = -2.61047101304182849292578695498722043d0
      mr_gamK(10, 9, 2) = 2.61047101304182849292578695498722043d0

      mr_wK = 0d0
      mr_wK(2, 1, 1) = 0.5d0

      mr_wK(4, 1, 1) = -1.91716534363662868878172216064946905d0
      mr_wK(4, 3, 1) = 2.04216534363662868878172216064946905d0

      mr_wK(5, 1, 1) = -0.404751031801105942697915907046990469d0
      mr_wK(5, 3, 1) = 0.404751031801105942697915907046990469d0

      mr_wK(6, 1, 1) = 11.4514660224922163666569802860263173d0
      mr_wK(6, 3, 1) = -30.2107574752650427144064781557395061d0
      mr_wK(6, 5, 1) = 18.8842914527728263477494978697131888d0

      mr_wK(7, 1, 1) = -0.709033564760261450684711672946330144d0
      mr_wK(7, 3, 1) = 1.03030720858751876652616190884004718d0
      mr_wK(7, 5, 1) = -0.321273643827257315841450235893717036d0

      mr_wK(8, 1, 1) = -29.9954871645582843984091068494419927d0
      mr_wK(8, 3, 1) = 37.605982774991801805364896856243857d0
      mr_wK(8, 5, 1) = 0.321273643827257315841450235893717036d0
      mr_wK(8, 7, 1) = -7.80676925426077472279724024269558129d0

      mr_wK(9, 1, 1) = 3.10466505427296211633876939184912422d0
      mr_wK(9, 3, 1) = -2.43032501975716229713206592741556636d0
      mr_wK(9, 5, 1) = -1.90547930115152463521920165948384213d0
      mr_wK(9, 7, 1) = 1.23113926663572481601249819505028427d0

      mr_wK(10, 1, 1) = -2.42442954775204786987587591435551401d0
      mr_wK(10, 3, 1) = 2.43032501975716229713206592741556636d0
      mr_wK(10, 5, 1) = 1.90547930115152463521920165948384213d0
      mr_wK(10, 7, 1) = -1.23113926663572481601249819505028427d0
      mr_wK(10, 9, 1) = -0.555235506520914246462893477493610215d0

      mr_wK(11, 1, 1) = -0.010441350444797485902945189451653542d0
      mr_wK(11, 3, 1) = 0.0726030361465507450515210450548814161d0
      mr_wK(11, 5, 1) = -0.128827595167726095223945409857642431d0
      mr_wK(11, 7, 1) = 0.112935535009382356613944010712215408d0
      mr_wK(11, 9, 1) = -0.0462696255434095205385744564578008512d0

      mr_wK(12, 1, 1) = -0.81085227877621013281757892286079321d0
      mr_wK(12, 3, 1) = 0.25600731992204924350015621921408823d0
      mr_wK(12, 5, 1) = 0.806829407269752789366586642278781947d0
      mr_wK(12, 7, 1) = -0.455714822872182379510589482174276116d0
      mr_wK(12, 9, 1) = -0.0462696255434095205385744564578008512d0
      mr_wK(12, 11, 1) = 0.25d0

      mr_wK(4, 1, 2) = 4.0843306872732573775634443212989381d0
      mr_wK(4, 3, 2) = -4.0843306872732573775634443212989381d0

      mr_wK(6, 1, 2) = -21.8434299813822208479181287579586536d0
      mr_wK(6, 3, 2) = 59.6120128869278735434171244973850312d0
      mr_wK(6, 5, 2) = -37.7685829055456526954989957394263776d0

      mr_wK(8, 1, 2) = 61.6590414586370916981876370447766458d0
      mr_wK(8, 3, 2) = -77.2725799671586411437821175301678084d0
      mr_wK(8, 7, 2) = 15.6135385085215494455944804853911626d0

      mr_wK(10, 1, 2) = -1.11047101304182849292578695498722043d0
      mr_wK(10, 9, 2) = 1.11047101304182849292578695498722043d0
   end subroutine gark4_init

   subroutine gark_init()
      use mr_data, only: mr_slowMethod, mr_kmax, mr_gamBar, mr_wBar, &
                         mr_gamK, mr_wK, mr_nstages_fast, mr_nstages_slow

      implicit none

      integer :: k

      select case (mr_slowMethod)
      case ("mr-gark3")
         call gark3_init

      case ("mr-gark4")
         call gark4_init
      end select

      allocate (mr_gamBar(mr_nstages_slow, mr_nstages_slow))
      allocate (mr_wBar(mr_nstages_slow, mr_nstages_slow))

      mr_gamBar = 0d0
      mr_wBar = 0d0

      do k = 1, mr_kmax
         mr_gamBar = mr_gamBar + mr_gamK(:, :, k)/k
         mr_wBar = mr_wBar + mr_wK(:, :, k)/k
      end do ! k
   end subroutine gark_init

end module gark
