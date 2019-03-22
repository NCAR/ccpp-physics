module mo_rrtmgp_lw_cloud_optics
  use machine,          only: kind_phys
  use physparam,        only: ilwcliq, ilwcice, iovrlw
  use mersenne_twister, only: random_setseed, random_number, random_stat

  implicit none

  ! Parameter used for RRTMG cloud-optics
  integer,parameter :: &
       nBandsLW_RRTMG = 16
  ! ipat is bands index for ebert & curry ice cloud (for iflagice=1)
  integer,dimension(nBandsLW_RRTMG),parameter :: &
       ipat = (/ 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5 /)
  real(kind_phys), parameter :: &
       absrain  = 0.33e-3, & ! Rain drop absorption coefficient \f$(m^{2}/g)\f$ .
       abssnow0 = 1.5,     & ! Snow flake absorption coefficient (micron), fu coeff
       abssnow1 = 2.34e-3    ! Snow flake absorption coefficient \f$(m^{2}/g)\f$, ncar coef
  real(kind_phys), parameter :: &
       cldmin  = 1e-20_kind_phys

  ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
  ! and 1.80) as a function of total column water vapor.  the function
  ! has been defined to minimize flux and cooling rate errors in these bands
  ! over a wide range of precipitable water values.
  real (kind_phys), dimension(nbandsLW_RRTMG) :: &
       a0 = (/ 1.66,  1.55,  1.58,  1.66,  1.54, 1.454,  1.89,  1.33,      &
              1.668,  1.66,  1.66,  1.66,  1.66,  1.66,  1.66,  1.66 /), &
       a1 = (/ 0.00,  0.25,  0.22,  0.00,  0.13, 0.446, -0.10,  0.40,     &
             -0.006,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /), &
       a2 = (/ 0.00, -12.0, -11.7,  0.00, -0.72,-0.243,  0.19,-0.062,  &
               0.414,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /)
  real(kind_phys),parameter :: &
       diffusivityLow   = 1.50, & ! Minimum diffusivity angle for bands 2-3 and 5-9
       diffusivityHigh  = 1.80, & ! Maximum diffusivity angle for bands 2-3 and 5-9
       diffusivityB1410 = 1.66    ! Diffusivity for bands 1, 4, and 10
 
  ! RRTMG LW cloud property coefficients
  real(kind_phys) , dimension(58,nBandsLW_RRTMG),parameter :: &
       absliq1 = reshape(source=(/ &
       1.64047e-03_kind_phys, 6.90533e-02_kind_phys, 7.72017e-02_kind_phys, 7.78054e-02_kind_phys, 7.69523e-02_kind_phys, & !1
       7.58058e-02_kind_phys, 7.46400e-02_kind_phys, 7.35123e-02_kind_phys, 7.24162e-02_kind_phys, 7.13225e-02_kind_phys, & !1
       6.99145e-02_kind_phys, 6.66409e-02_kind_phys, 6.36582e-02_kind_phys, 6.09425e-02_kind_phys, 5.84593e-02_kind_phys, & !1
       5.61743e-02_kind_phys, 5.40571e-02_kind_phys, 5.20812e-02_kind_phys, 5.02245e-02_kind_phys, 4.84680e-02_kind_phys, & !1
       4.67959e-02_kind_phys, 4.51944e-02_kind_phys, 4.36516e-02_kind_phys, 4.21570e-02_kind_phys, 4.07015e-02_kind_phys, & !1
       3.92766e-02_kind_phys, 3.78747e-02_kind_phys, 3.64886e-02_kind_phys, 3.53632e-02_kind_phys, 3.41992e-02_kind_phys, & !1
       3.31016e-02_kind_phys, 3.20643e-02_kind_phys, 3.10817e-02_kind_phys, 3.01490e-02_kind_phys, 2.92620e-02_kind_phys, & !1
       2.84171e-02_kind_phys, 2.76108e-02_kind_phys, 2.68404e-02_kind_phys, 2.61031e-02_kind_phys, 2.53966e-02_kind_phys, & !1
       2.47189e-02_kind_phys, 2.40678e-02_kind_phys, 2.34418e-02_kind_phys, 2.28392e-02_kind_phys, 2.22586e-02_kind_phys, & !1
       2.16986e-02_kind_phys, 2.11580e-02_kind_phys, 2.06356e-02_kind_phys, 2.01305e-02_kind_phys, 1.96417e-02_kind_phys, & !1
       1.91682e-02_kind_phys, 1.87094e-02_kind_phys, 1.82643e-02_kind_phys, 1.78324e-02_kind_phys, 1.74129e-02_kind_phys, & !1
       1.70052e-02_kind_phys, 1.66088e-02_kind_phys, 1.62231e-02_kind_phys,                                               & !1
       2.19486e-01_kind_phys, 1.80687e-01_kind_phys, 1.59150e-01_kind_phys, 1.44731e-01_kind_phys, 1.33703e-01_kind_phys, & !2
       1.24355e-01_kind_phys, 1.15756e-01_kind_phys, 1.07318e-01_kind_phys, 9.86119e-02_kind_phys, 8.92739e-02_kind_phys, & !2
       8.34911e-02_kind_phys, 7.70773e-02_kind_phys, 7.15240e-02_kind_phys, 6.66615e-02_kind_phys, 6.23641e-02_kind_phys, & !2
       5.85359e-02_kind_phys, 5.51020e-02_kind_phys, 5.20032e-02_kind_phys, 4.91916e-02_kind_phys, 4.66283e-02_kind_phys, & !2
       4.42813e-02_kind_phys, 4.21236e-02_kind_phys, 4.01330e-02_kind_phys, 3.82905e-02_kind_phys, 3.65797e-02_kind_phys, & !2
       3.49869e-02_kind_phys, 3.35002e-02_kind_phys, 3.21090e-02_kind_phys, 3.08957e-02_kind_phys, 2.97601e-02_kind_phys, & !2
       2.86966e-02_kind_phys, 2.76984e-02_kind_phys, 2.67599e-02_kind_phys, 2.58758e-02_kind_phys, 2.50416e-02_kind_phys, & !2
       2.42532e-02_kind_phys, 2.35070e-02_kind_phys, 2.27997e-02_kind_phys, 2.21284e-02_kind_phys, 2.14904e-02_kind_phys, & !2
       2.08834e-02_kind_phys, 2.03051e-02_kind_phys, 1.97536e-02_kind_phys, 1.92271e-02_kind_phys, 1.87239e-02_kind_phys, & !2
       1.82425e-02_kind_phys, 1.77816e-02_kind_phys, 1.73399e-02_kind_phys, 1.69162e-02_kind_phys, 1.65094e-02_kind_phys, & !2
       1.61187e-02_kind_phys, 1.57430e-02_kind_phys, 1.53815e-02_kind_phys, 1.50334e-02_kind_phys, 1.46981e-02_kind_phys, & !2
       1.43748e-02_kind_phys, 1.40628e-02_kind_phys, 1.37617e-02_kind_phys,                                               & !2
       2.95174e-01_kind_phys, 2.34765e-01_kind_phys, 1.98038e-01_kind_phys, 1.72114e-01_kind_phys, 1.52083e-01_kind_phys, & !3
       1.35654e-01_kind_phys, 1.21613e-01_kind_phys, 1.09252e-01_kind_phys, 9.81263e-02_kind_phys, 8.79448e-02_kind_phys, & !3
       8.12566e-02_kind_phys, 7.44563e-02_kind_phys, 6.86374e-02_kind_phys, 6.36042e-02_kind_phys, 5.92094e-02_kind_phys, & !3
       5.53402e-02_kind_phys, 5.19087e-02_kind_phys, 4.88455e-02_kind_phys, 4.60951e-02_kind_phys, 4.36124e-02_kind_phys, & !3
       4.13607e-02_kind_phys, 3.93096e-02_kind_phys, 3.74338e-02_kind_phys, 3.57119e-02_kind_phys, 3.41261e-02_kind_phys, & !3
       3.26610e-02_kind_phys, 3.13036e-02_kind_phys, 3.00425e-02_kind_phys, 2.88497e-02_kind_phys, 2.78077e-02_kind_phys, & !3
       2.68317e-02_kind_phys, 2.59158e-02_kind_phys, 2.50545e-02_kind_phys, 2.42430e-02_kind_phys, 2.34772e-02_kind_phys, & !3
       2.27533e-02_kind_phys, 2.20679e-02_kind_phys, 2.14181e-02_kind_phys, 2.08011e-02_kind_phys, 2.02145e-02_kind_phys, & !3
       1.96561e-02_kind_phys, 1.91239e-02_kind_phys, 1.86161e-02_kind_phys, 1.81311e-02_kind_phys, 1.76673e-02_kind_phys, & !3
       1.72234e-02_kind_phys, 1.67981e-02_kind_phys, 1.63903e-02_kind_phys, 1.59989e-02_kind_phys, 1.56230e-02_kind_phys, & !3
       1.52615e-02_kind_phys, 1.49138e-02_kind_phys, 1.45791e-02_kind_phys, 1.42565e-02_kind_phys, 1.39455e-02_kind_phys, & !3
       1.36455e-02_kind_phys, 1.33559e-02_kind_phys, 1.30761e-02_kind_phys,                                               & !3
       3.00925e-01_kind_phys, 2.36949e-01_kind_phys, 1.96947e-01_kind_phys, 1.68692e-01_kind_phys, 1.47190e-01_kind_phys, & !4
       1.29986e-01_kind_phys, 1.15719e-01_kind_phys, 1.03568e-01_kind_phys, 9.30028e-02_kind_phys, 8.36658e-02_kind_phys, & !4
       7.71075e-02_kind_phys, 7.07002e-02_kind_phys, 6.52284e-02_kind_phys, 6.05024e-02_kind_phys, 5.63801e-02_kind_phys, & !4
       5.27534e-02_kind_phys, 4.95384e-02_kind_phys, 4.66690e-02_kind_phys, 4.40925e-02_kind_phys, 4.17664e-02_kind_phys, & !4
       3.96559e-02_kind_phys, 3.77326e-02_kind_phys, 3.59727e-02_kind_phys, 3.43561e-02_kind_phys, 3.28662e-02_kind_phys, & !4
       3.14885e-02_kind_phys, 3.02110e-02_kind_phys, 2.90231e-02_kind_phys, 2.78948e-02_kind_phys, 2.69109e-02_kind_phys, & !4
       2.59884e-02_kind_phys, 2.51217e-02_kind_phys, 2.43058e-02_kind_phys, 2.35364e-02_kind_phys, 2.28096e-02_kind_phys, & !4
       2.21218e-02_kind_phys, 2.14700e-02_kind_phys, 2.08515e-02_kind_phys, 2.02636e-02_kind_phys, 1.97041e-02_kind_phys, & !4
       1.91711e-02_kind_phys, 1.86625e-02_kind_phys, 1.81769e-02_kind_phys, 1.77126e-02_kind_phys, 1.72683e-02_kind_phys, & !4
       1.68426e-02_kind_phys, 1.64344e-02_kind_phys, 1.60427e-02_kind_phys, 1.56664e-02_kind_phys, 1.53046e-02_kind_phys, & !4
       1.49565e-02_kind_phys, 1.46214e-02_kind_phys, 1.42985e-02_kind_phys, 1.39871e-02_kind_phys, 1.36866e-02_kind_phys, & !4
       1.33965e-02_kind_phys, 1.31162e-02_kind_phys, 1.28453e-02_kind_phys,                                               & !4
       2.64691e-01_kind_phys, 2.12018e-01_kind_phys, 1.78009e-01_kind_phys, 1.53539e-01_kind_phys, 1.34721e-01_kind_phys, & !5
       1.19580e-01_kind_phys, 1.06996e-01_kind_phys, 9.62772e-02_kind_phys, 8.69710e-02_kind_phys, 7.87670e-02_kind_phys, & !5
       7.29272e-02_kind_phys, 6.70920e-02_kind_phys, 6.20977e-02_kind_phys, 5.77732e-02_kind_phys, 5.39910e-02_kind_phys, & !5
       5.06538e-02_kind_phys, 4.76866e-02_kind_phys, 4.50301e-02_kind_phys, 4.26374e-02_kind_phys, 4.04704e-02_kind_phys, & !5
       3.84981e-02_kind_phys, 3.66948e-02_kind_phys, 3.50394e-02_kind_phys, 3.35141e-02_kind_phys, 3.21038e-02_kind_phys, & !5
       3.07957e-02_kind_phys, 2.95788e-02_kind_phys, 2.84438e-02_kind_phys, 2.73790e-02_kind_phys, 2.64390e-02_kind_phys, & !5
       2.55565e-02_kind_phys, 2.47263e-02_kind_phys, 2.39437e-02_kind_phys, 2.32047e-02_kind_phys, 2.25056e-02_kind_phys, & !5
       2.18433e-02_kind_phys, 2.12149e-02_kind_phys, 2.06177e-02_kind_phys, 2.00495e-02_kind_phys, 1.95081e-02_kind_phys, & !5
       1.89917e-02_kind_phys, 1.84984e-02_kind_phys, 1.80269e-02_kind_phys, 1.75755e-02_kind_phys, 1.71431e-02_kind_phys, & !5
       1.67283e-02_kind_phys, 1.63303e-02_kind_phys, 1.59478e-02_kind_phys, 1.55801e-02_kind_phys, 1.52262e-02_kind_phys, & !5
       1.48853e-02_kind_phys, 1.45568e-02_kind_phys, 1.42400e-02_kind_phys, 1.39342e-02_kind_phys, 1.36388e-02_kind_phys, & !5
       1.33533e-02_kind_phys, 1.30773e-02_kind_phys, 1.28102e-02_kind_phys,                                               & !5
       8.81182e-02_kind_phys, 1.06745e-01_kind_phys, 9.79753e-02_kind_phys, 8.99625e-02_kind_phys, 8.35200e-02_kind_phys, & !6
       7.81899e-02_kind_phys, 7.35939e-02_kind_phys, 6.94696e-02_kind_phys, 6.56266e-02_kind_phys, 6.19148e-02_kind_phys, & !6
       5.83355e-02_kind_phys, 5.49306e-02_kind_phys, 5.19642e-02_kind_phys, 4.93325e-02_kind_phys, 4.69659e-02_kind_phys, & !6
       4.48148e-02_kind_phys, 4.28431e-02_kind_phys, 4.10231e-02_kind_phys, 3.93332e-02_kind_phys, 3.77563e-02_kind_phys, & !6
       3.62785e-02_kind_phys, 3.48882e-02_kind_phys, 3.35758e-02_kind_phys, 3.23333e-02_kind_phys, 3.11536e-02_kind_phys, & !6
       3.00310e-02_kind_phys, 2.89601e-02_kind_phys, 2.79365e-02_kind_phys, 2.70502e-02_kind_phys, 2.62618e-02_kind_phys, & !6
       2.55025e-02_kind_phys, 2.47728e-02_kind_phys, 2.40726e-02_kind_phys, 2.34013e-02_kind_phys, 2.27583e-02_kind_phys, & !6
       2.21422e-02_kind_phys, 2.15522e-02_kind_phys, 2.09869e-02_kind_phys, 2.04453e-02_kind_phys, 1.99260e-02_kind_phys, & !6
       1.94280e-02_kind_phys, 1.89501e-02_kind_phys, 1.84913e-02_kind_phys, 1.80506e-02_kind_phys, 1.76270e-02_kind_phys, & !6
       1.72196e-02_kind_phys, 1.68276e-02_kind_phys, 1.64500e-02_kind_phys, 1.60863e-02_kind_phys, 1.57357e-02_kind_phys, & !6
       1.53975e-02_kind_phys, 1.50710e-02_kind_phys, 1.47558e-02_kind_phys, 1.44511e-02_kind_phys, 1.41566e-02_kind_phys, & !6
       1.38717e-02_kind_phys, 1.35960e-02_kind_phys, 1.33290e-02_kind_phys,                                               & !6
       4.32174e-02_kind_phys, 7.36078e-02_kind_phys, 6.98340e-02_kind_phys, 6.65231e-02_kind_phys, 6.41948e-02_kind_phys, & !7
       6.23551e-02_kind_phys, 6.06638e-02_kind_phys, 5.88680e-02_kind_phys, 5.67124e-02_kind_phys, 5.38629e-02_kind_phys, & !7
       4.99579e-02_kind_phys, 4.86289e-02_kind_phys, 4.70120e-02_kind_phys, 4.52854e-02_kind_phys, 4.35466e-02_kind_phys, & !7
       4.18480e-02_kind_phys, 4.02169e-02_kind_phys, 3.86658e-02_kind_phys, 3.71992e-02_kind_phys, 3.58168e-02_kind_phys, & !7
       3.45155e-02_kind_phys, 3.32912e-02_kind_phys, 3.21390e-02_kind_phys, 3.10538e-02_kind_phys, 3.00307e-02_kind_phys, & !7
       2.90651e-02_kind_phys, 2.81524e-02_kind_phys, 2.72885e-02_kind_phys, 2.62821e-02_kind_phys, 2.55744e-02_kind_phys, & !7
       2.48799e-02_kind_phys, 2.42029e-02_kind_phys, 2.35460e-02_kind_phys, 2.29108e-02_kind_phys, 2.22981e-02_kind_phys, & !7
       2.17079e-02_kind_phys, 2.11402e-02_kind_phys, 2.05945e-02_kind_phys, 2.00701e-02_kind_phys, 1.95663e-02_kind_phys, & !7
       1.90824e-02_kind_phys, 1.86174e-02_kind_phys, 1.81706e-02_kind_phys, 1.77411e-02_kind_phys, 1.73281e-02_kind_phys, & !7
       1.69307e-02_kind_phys, 1.65483e-02_kind_phys, 1.61801e-02_kind_phys, 1.58254e-02_kind_phys, 1.54835e-02_kind_phys, & !7
       1.51538e-02_kind_phys, 1.48358e-02_kind_phys, 1.45288e-02_kind_phys, 1.42322e-02_kind_phys, 1.39457e-02_kind_phys, & !7
       1.36687e-02_kind_phys, 1.34008e-02_kind_phys, 1.31416e-02_kind_phys,                                               & !7
       1.41881e-01_kind_phys, 7.15419e-02_kind_phys, 6.30335e-02_kind_phys, 6.11132e-02_kind_phys, 6.01931e-02_kind_phys, & !8
       5.92420e-02_kind_phys, 5.78968e-02_kind_phys, 5.58876e-02_kind_phys, 5.28923e-02_kind_phys, 4.84462e-02_kind_phys, & !8
       4.60839e-02_kind_phys, 4.56013e-02_kind_phys, 4.45410e-02_kind_phys, 4.31866e-02_kind_phys, 4.17026e-02_kind_phys, & !8
       4.01850e-02_kind_phys, 3.86892e-02_kind_phys, 3.72461e-02_kind_phys, 3.58722e-02_kind_phys, 3.45749e-02_kind_phys, & !8
       3.33564e-02_kind_phys, 3.22155e-02_kind_phys, 3.11494e-02_kind_phys, 3.01541e-02_kind_phys, 2.92253e-02_kind_phys, & !8
       2.83584e-02_kind_phys, 2.75488e-02_kind_phys, 2.67925e-02_kind_phys, 2.57692e-02_kind_phys, 2.50704e-02_kind_phys, & !8
       2.43918e-02_kind_phys, 2.37350e-02_kind_phys, 2.31005e-02_kind_phys, 2.24888e-02_kind_phys, 2.18996e-02_kind_phys, & !8
       2.13325e-02_kind_phys, 2.07870e-02_kind_phys, 2.02623e-02_kind_phys, 1.97577e-02_kind_phys, 1.92724e-02_kind_phys, & !8
       1.88056e-02_kind_phys, 1.83564e-02_kind_phys, 1.79241e-02_kind_phys, 1.75079e-02_kind_phys, 1.71070e-02_kind_phys, & !8
       1.67207e-02_kind_phys, 1.63482e-02_kind_phys, 1.59890e-02_kind_phys, 1.56424e-02_kind_phys, 1.53077e-02_kind_phys, & !8
       1.49845e-02_kind_phys, 1.46722e-02_kind_phys, 1.43702e-02_kind_phys, 1.40782e-02_kind_phys, 1.37955e-02_kind_phys, & !8
       1.35219e-02_kind_phys, 1.32569e-02_kind_phys, 1.30000e-02_kind_phys,                                               & !8
       6.72726e-02_kind_phys, 6.61013e-02_kind_phys, 6.47866e-02_kind_phys, 6.33780e-02_kind_phys, 6.18985e-02_kind_phys, & !9
       6.03335e-02_kind_phys, 5.86136e-02_kind_phys, 5.65876e-02_kind_phys, 5.39839e-02_kind_phys, 5.03536e-02_kind_phys, & !9
       4.71608e-02_kind_phys, 4.63630e-02_kind_phys, 4.50313e-02_kind_phys, 4.34526e-02_kind_phys, 4.17876e-02_kind_phys, & !9
       4.01261e-02_kind_phys, 3.85171e-02_kind_phys, 3.69860e-02_kind_phys, 3.55442e-02_kind_phys, 3.41954e-02_kind_phys, & !9
       3.29384e-02_kind_phys, 3.17693e-02_kind_phys, 3.06832e-02_kind_phys, 2.96745e-02_kind_phys, 2.87374e-02_kind_phys, & !9
       2.78662e-02_kind_phys, 2.70557e-02_kind_phys, 2.63008e-02_kind_phys, 2.52450e-02_kind_phys, 2.45424e-02_kind_phys, & !9
       2.38656e-02_kind_phys, 2.32144e-02_kind_phys, 2.25885e-02_kind_phys, 2.19873e-02_kind_phys, 2.14099e-02_kind_phys, & !9
       2.08554e-02_kind_phys, 2.03230e-02_kind_phys, 1.98116e-02_kind_phys, 1.93203e-02_kind_phys, 1.88482e-02_kind_phys, & !9
       1.83944e-02_kind_phys, 1.79578e-02_kind_phys, 1.75378e-02_kind_phys, 1.71335e-02_kind_phys, 1.67440e-02_kind_phys, & !9
       1.63687e-02_kind_phys, 1.60069e-02_kind_phys, 1.56579e-02_kind_phys, 1.53210e-02_kind_phys, 1.49958e-02_kind_phys, & !9
       1.46815e-02_kind_phys, 1.43778e-02_kind_phys, 1.40841e-02_kind_phys, 1.37999e-02_kind_phys, 1.35249e-02_kind_phys, & !9
       1.32585e-02_kind_phys, 1.30004e-02_kind_phys, 1.27502e-02_kind_phys,                                               & !9
       7.97040e-02_kind_phys, 7.63844e-02_kind_phys, 7.36499e-02_kind_phys, 7.13525e-02_kind_phys, 6.93043e-02_kind_phys, & !10
       6.72807e-02_kind_phys, 6.50227e-02_kind_phys, 6.22395e-02_kind_phys, 5.86093e-02_kind_phys, 5.37815e-02_kind_phys, & !10
       5.14682e-02_kind_phys, 4.97214e-02_kind_phys, 4.77392e-02_kind_phys, 4.56961e-02_kind_phys, 4.36858e-02_kind_phys, & !10
       4.17569e-02_kind_phys, 3.99328e-02_kind_phys, 3.82224e-02_kind_phys, 3.66265e-02_kind_phys, 3.51416e-02_kind_phys, & !10
       3.37617e-02_kind_phys, 3.24798e-02_kind_phys, 3.12887e-02_kind_phys, 3.01812e-02_kind_phys, 2.91505e-02_kind_phys, & !10
       2.81900e-02_kind_phys, 2.72939e-02_kind_phys, 2.64568e-02_kind_phys, 2.54165e-02_kind_phys, 2.46832e-02_kind_phys, & !10
       2.39783e-02_kind_phys, 2.33017e-02_kind_phys, 2.26531e-02_kind_phys, 2.20314e-02_kind_phys, 2.14359e-02_kind_phys, & !10
       2.08653e-02_kind_phys, 2.03187e-02_kind_phys, 1.97947e-02_kind_phys, 1.92924e-02_kind_phys, 1.88106e-02_kind_phys, & !10
       1.83483e-02_kind_phys, 1.79043e-02_kind_phys, 1.74778e-02_kind_phys, 1.70678e-02_kind_phys, 1.66735e-02_kind_phys, & !10
       1.62941e-02_kind_phys, 1.59286e-02_kind_phys, 1.55766e-02_kind_phys, 1.52371e-02_kind_phys, 1.49097e-02_kind_phys, & !10
       1.45937e-02_kind_phys, 1.42885e-02_kind_phys, 1.39936e-02_kind_phys, 1.37085e-02_kind_phys, 1.34327e-02_kind_phys, & !10
       1.31659e-02_kind_phys, 1.29075e-02_kind_phys, 1.26571e-02_kind_phys,                                               & !10
       1.49438e-01_kind_phys, 1.33535e-01_kind_phys, 1.21542e-01_kind_phys, 1.11743e-01_kind_phys, 1.03263e-01_kind_phys, & !11
       9.55774e-02_kind_phys, 8.83382e-02_kind_phys, 8.12943e-02_kind_phys, 7.42533e-02_kind_phys, 6.70609e-02_kind_phys, & !11
       6.38761e-02_kind_phys, 5.97788e-02_kind_phys, 5.59841e-02_kind_phys, 5.25318e-02_kind_phys, 4.94132e-02_kind_phys, & !11
       4.66014e-02_kind_phys, 4.40644e-02_kind_phys, 4.17706e-02_kind_phys, 3.96910e-02_kind_phys, 3.77998e-02_kind_phys, & !11
       3.60742e-02_kind_phys, 3.44947e-02_kind_phys, 3.30442e-02_kind_phys, 3.17079e-02_kind_phys, 3.04730e-02_kind_phys, & !11
       2.93283e-02_kind_phys, 2.82642e-02_kind_phys, 2.72720e-02_kind_phys, 2.61789e-02_kind_phys, 2.53277e-02_kind_phys, & !11
       2.45237e-02_kind_phys, 2.37635e-02_kind_phys, 2.30438e-02_kind_phys, 2.23615e-02_kind_phys, 2.17140e-02_kind_phys, & !11
       2.10987e-02_kind_phys, 2.05133e-02_kind_phys, 1.99557e-02_kind_phys, 1.94241e-02_kind_phys, 1.89166e-02_kind_phys, & !11
       1.84317e-02_kind_phys, 1.79679e-02_kind_phys, 1.75238e-02_kind_phys, 1.70983e-02_kind_phys, 1.66901e-02_kind_phys, & !11
       1.62983e-02_kind_phys, 1.59219e-02_kind_phys, 1.55599e-02_kind_phys, 1.52115e-02_kind_phys, 1.48761e-02_kind_phys, & !11
       1.45528e-02_kind_phys, 1.42411e-02_kind_phys, 1.39402e-02_kind_phys, 1.36497e-02_kind_phys, 1.33690e-02_kind_phys, & !11
       1.30976e-02_kind_phys, 1.28351e-02_kind_phys, 1.25810e-02_kind_phys,                                               & !11
       3.71985e-02_kind_phys, 3.88586e-02_kind_phys, 3.99070e-02_kind_phys, 4.04351e-02_kind_phys, 4.04610e-02_kind_phys, & !12
       3.99834e-02_kind_phys, 3.89953e-02_kind_phys, 3.74886e-02_kind_phys, 3.54551e-02_kind_phys, 3.28870e-02_kind_phys, & !12
       3.32576e-02_kind_phys, 3.22444e-02_kind_phys, 3.12384e-02_kind_phys, 3.02584e-02_kind_phys, 2.93146e-02_kind_phys, & !12
       2.84120e-02_kind_phys, 2.75525e-02_kind_phys, 2.67361e-02_kind_phys, 2.59618e-02_kind_phys, 2.52280e-02_kind_phys, & !12
       2.45327e-02_kind_phys, 2.38736e-02_kind_phys, 2.32487e-02_kind_phys, 2.26558e-02_kind_phys, 2.20929e-02_kind_phys, & !12
       2.15579e-02_kind_phys, 2.10491e-02_kind_phys, 2.05648e-02_kind_phys, 1.99749e-02_kind_phys, 1.95704e-02_kind_phys, & !12
       1.91731e-02_kind_phys, 1.87839e-02_kind_phys, 1.84032e-02_kind_phys, 1.80315e-02_kind_phys, 1.76689e-02_kind_phys, & !12
       1.73155e-02_kind_phys, 1.69712e-02_kind_phys, 1.66362e-02_kind_phys, 1.63101e-02_kind_phys, 1.59928e-02_kind_phys, & !12
       1.56842e-02_kind_phys, 1.53840e-02_kind_phys, 1.50920e-02_kind_phys, 1.48080e-02_kind_phys, 1.45318e-02_kind_phys, & !12
       1.42631e-02_kind_phys, 1.40016e-02_kind_phys, 1.37472e-02_kind_phys, 1.34996e-02_kind_phys, 1.32586e-02_kind_phys, & !12
       1.30239e-02_kind_phys, 1.27954e-02_kind_phys, 1.25728e-02_kind_phys, 1.23559e-02_kind_phys, 1.21445e-02_kind_phys, & !12
       1.19385e-02_kind_phys, 1.17376e-02_kind_phys, 1.15417e-02_kind_phys,                                               & !12
       3.11868e-02_kind_phys, 4.48357e-02_kind_phys, 4.90224e-02_kind_phys, 4.96406e-02_kind_phys, 4.86806e-02_kind_phys, & !13
       4.69610e-02_kind_phys, 4.48630e-02_kind_phys, 4.25795e-02_kind_phys, 4.02138e-02_kind_phys, 3.78236e-02_kind_phys, & !13
       3.74266e-02_kind_phys, 3.60384e-02_kind_phys, 3.47074e-02_kind_phys, 3.34434e-02_kind_phys, 3.22499e-02_kind_phys, & !13
       3.11264e-02_kind_phys, 3.00704e-02_kind_phys, 2.90784e-02_kind_phys, 2.81463e-02_kind_phys, 2.72702e-02_kind_phys, & !13
       2.64460e-02_kind_phys, 2.56698e-02_kind_phys, 2.49381e-02_kind_phys, 2.42475e-02_kind_phys, 2.35948e-02_kind_phys, & !13
       2.29774e-02_kind_phys, 2.23925e-02_kind_phys, 2.18379e-02_kind_phys, 2.11793e-02_kind_phys, 2.07076e-02_kind_phys, & !13
       2.02470e-02_kind_phys, 1.97981e-02_kind_phys, 1.93613e-02_kind_phys, 1.89367e-02_kind_phys, 1.85243e-02_kind_phys, & !13
       1.81240e-02_kind_phys, 1.77356e-02_kind_phys, 1.73588e-02_kind_phys, 1.69935e-02_kind_phys, 1.66392e-02_kind_phys, & !13
       1.62956e-02_kind_phys, 1.59624e-02_kind_phys, 1.56393e-02_kind_phys, 1.53259e-02_kind_phys, 1.50219e-02_kind_phys, & !13
       1.47268e-02_kind_phys, 1.44404e-02_kind_phys, 1.41624e-02_kind_phys, 1.38925e-02_kind_phys, 1.36302e-02_kind_phys, & !13
       1.33755e-02_kind_phys, 1.31278e-02_kind_phys, 1.28871e-02_kind_phys, 1.26530e-02_kind_phys, 1.24253e-02_kind_phys, & !13
       1.22038e-02_kind_phys, 1.19881e-02_kind_phys, 1.17782e-02_kind_phys,                                               & !13
       1.58988e-02_kind_phys, 3.50652e-02_kind_phys, 4.00851e-02_kind_phys, 4.07270e-02_kind_phys, 3.98101e-02_kind_phys, & !14
       3.83306e-02_kind_phys, 3.66829e-02_kind_phys, 3.50327e-02_kind_phys, 3.34497e-02_kind_phys, 3.19609e-02_kind_phys, & !14
       3.13712e-02_kind_phys, 3.03348e-02_kind_phys, 2.93415e-02_kind_phys, 2.83973e-02_kind_phys, 2.75037e-02_kind_phys, & !14
       2.66604e-02_kind_phys, 2.58654e-02_kind_phys, 2.51161e-02_kind_phys, 2.44100e-02_kind_phys, 2.37440e-02_kind_phys, & !14
       2.31154e-02_kind_phys, 2.25215e-02_kind_phys, 2.19599e-02_kind_phys, 2.14282e-02_kind_phys, 2.09242e-02_kind_phys, & !14
       2.04459e-02_kind_phys, 1.99915e-02_kind_phys, 1.95594e-02_kind_phys, 1.90254e-02_kind_phys, 1.86598e-02_kind_phys, & !14
       1.82996e-02_kind_phys, 1.79455e-02_kind_phys, 1.75983e-02_kind_phys, 1.72584e-02_kind_phys, 1.69260e-02_kind_phys, & !14
       1.66013e-02_kind_phys, 1.62843e-02_kind_phys, 1.59752e-02_kind_phys, 1.56737e-02_kind_phys, 1.53799e-02_kind_phys, & !14
       1.50936e-02_kind_phys, 1.48146e-02_kind_phys, 1.45429e-02_kind_phys, 1.42782e-02_kind_phys, 1.40203e-02_kind_phys, & !14
       1.37691e-02_kind_phys, 1.35243e-02_kind_phys, 1.32858e-02_kind_phys, 1.30534e-02_kind_phys, 1.28270e-02_kind_phys, & !14
       1.26062e-02_kind_phys, 1.23909e-02_kind_phys, 1.21810e-02_kind_phys, 1.19763e-02_kind_phys, 1.17766e-02_kind_phys, & !14
       1.15817e-02_kind_phys, 1.13915e-02_kind_phys, 1.12058e-02_kind_phys,                                               & !14
       5.02079e-03_kind_phys, 2.17615e-02_kind_phys, 2.55449e-02_kind_phys, 2.59484e-02_kind_phys, 2.53650e-02_kind_phys, & !15
       2.45281e-02_kind_phys, 2.36843e-02_kind_phys, 2.29159e-02_kind_phys, 2.22451e-02_kind_phys, 2.16716e-02_kind_phys, & !15
       2.11451e-02_kind_phys, 2.05817e-02_kind_phys, 2.00454e-02_kind_phys, 1.95372e-02_kind_phys, 1.90567e-02_kind_phys, & !15
       1.86028e-02_kind_phys, 1.81742e-02_kind_phys, 1.77693e-02_kind_phys, 1.73866e-02_kind_phys, 1.70244e-02_kind_phys, & !15
       1.66815e-02_kind_phys, 1.63563e-02_kind_phys, 1.60477e-02_kind_phys, 1.57544e-02_kind_phys, 1.54755e-02_kind_phys, & !15
       1.52097e-02_kind_phys, 1.49564e-02_kind_phys, 1.47146e-02_kind_phys, 1.43684e-02_kind_phys, 1.41728e-02_kind_phys, & !15
       1.39762e-02_kind_phys, 1.37797e-02_kind_phys, 1.35838e-02_kind_phys, 1.33891e-02_kind_phys, 1.31961e-02_kind_phys, & !15
       1.30051e-02_kind_phys, 1.28164e-02_kind_phys, 1.26302e-02_kind_phys, 1.24466e-02_kind_phys, 1.22659e-02_kind_phys, & !15
       1.20881e-02_kind_phys, 1.19131e-02_kind_phys, 1.17412e-02_kind_phys, 1.15723e-02_kind_phys, 1.14063e-02_kind_phys, & !15
       1.12434e-02_kind_phys, 1.10834e-02_kind_phys, 1.09264e-02_kind_phys, 1.07722e-02_kind_phys, 1.06210e-02_kind_phys, & !15
       1.04725e-02_kind_phys, 1.03269e-02_kind_phys, 1.01839e-02_kind_phys, 1.00436e-02_kind_phys, 9.90593e-03_kind_phys, & !15
       9.77080e-03_kind_phys, 9.63818e-03_kind_phys, 9.50800e-03_kind_phys,                                               & !15
       5.64971e-02_kind_phys, 9.04736e-02_kind_phys, 8.11726e-02_kind_phys, 7.05450e-02_kind_phys, 6.20052e-02_kind_phys, & !16
       5.54286e-02_kind_phys, 5.03503e-02_kind_phys, 4.63791e-02_kind_phys, 4.32290e-02_kind_phys, 4.06959e-02_kind_phys, & !16
       3.74690e-02_kind_phys, 3.52964e-02_kind_phys, 3.33799e-02_kind_phys, 3.16774e-02_kind_phys, 3.01550e-02_kind_phys, & !16
       2.87856e-02_kind_phys, 2.75474e-02_kind_phys, 2.64223e-02_kind_phys, 2.53953e-02_kind_phys, 2.44542e-02_kind_phys, & !16
       2.35885e-02_kind_phys, 2.27894e-02_kind_phys, 2.20494e-02_kind_phys, 2.13622e-02_kind_phys, 2.07222e-02_kind_phys, & !16
       2.01246e-02_kind_phys, 1.95654e-02_kind_phys, 1.90408e-02_kind_phys, 1.84398e-02_kind_phys, 1.80021e-02_kind_phys, & !16
       1.75816e-02_kind_phys, 1.71775e-02_kind_phys, 1.67889e-02_kind_phys, 1.64152e-02_kind_phys, 1.60554e-02_kind_phys, & !16
       1.57089e-02_kind_phys, 1.53751e-02_kind_phys, 1.50531e-02_kind_phys, 1.47426e-02_kind_phys, 1.44428e-02_kind_phys, & !16
       1.41532e-02_kind_phys, 1.38734e-02_kind_phys, 1.36028e-02_kind_phys, 1.33410e-02_kind_phys, 1.30875e-02_kind_phys, & !16
       1.28420e-02_kind_phys, 1.26041e-02_kind_phys, 1.23735e-02_kind_phys, 1.21497e-02_kind_phys, 1.19325e-02_kind_phys, & !16
       1.17216e-02_kind_phys, 1.15168e-02_kind_phys, 1.13177e-02_kind_phys, 1.11241e-02_kind_phys, 1.09358e-02_kind_phys, & !16
       1.07525e-02_kind_phys, 1.05741e-02_kind_phys, 1.04003e-02_kind_phys/),                                             & !16
       shape=(/58,nBandsLW_RRTMG/))
  
  real(kind_phys), dimension(2),parameter :: &
       absice0 = (/0.005,1.0/)

  real(kind_phys), dimension(2,5),parameter :: &
       absice1 = reshape(source=(/  &
       0.0036, 1.136,  0.0068, 0.600, 0.0003, 1.338, 0.0016, 1.166,  0.0020, 1.118  /),&
       shape=(/2,5/))

  real(kind_phys), dimension(43, nBandsLW_RRTMG),parameter :: &
       absice2 =  reshape(source=(/ &
       7.798999e-02_kind_phys, 6.340479e-02_kind_phys, 5.417973e-02_kind_phys, 4.766245e-02_kind_phys, 4.272663e-02_kind_phys,  & !1
       3.880939e-02_kind_phys, 3.559544e-02_kind_phys, 3.289241e-02_kind_phys, 3.057511e-02_kind_phys, 2.855800e-02_kind_phys,  & !1
       2.678022e-02_kind_phys, 2.519712e-02_kind_phys, 2.377505e-02_kind_phys, 2.248806e-02_kind_phys, 2.131578e-02_kind_phys,  & !1
       2.024194e-02_kind_phys, 1.925337e-02_kind_phys, 1.833926e-02_kind_phys, 1.749067e-02_kind_phys, 1.670007e-02_kind_phys,  & !1
       1.596113e-02_kind_phys, 1.526845e-02_kind_phys, 1.461739e-02_kind_phys, 1.400394e-02_kind_phys, 1.342462e-02_kind_phys,  & !1
       1.287639e-02_kind_phys, 1.235656e-02_kind_phys, 1.186279e-02_kind_phys, 1.139297e-02_kind_phys, 1.094524e-02_kind_phys,  & !1
       1.051794e-02_kind_phys, 1.010956e-02_kind_phys, 9.718755e-03_kind_phys, 9.344316e-03_kind_phys, 8.985139e-03_kind_phys,  & !1
       8.640223e-03_kind_phys, 8.308656e-03_kind_phys, 7.989606e-03_kind_phys, 7.682312e-03_kind_phys, 7.386076e-03_kind_phys,  & !1
       7.100255e-03_kind_phys, 6.824258e-03_kind_phys, 6.557540e-03_kind_phys,                                                  & !1
       2.784879e-02_kind_phys, 2.709863e-02_kind_phys, 2.619165e-02_kind_phys, 2.529230e-02_kind_phys, 2.443225e-02_kind_phys,  & !2
       2.361575e-02_kind_phys, 2.284021e-02_kind_phys, 2.210150e-02_kind_phys, 2.139548e-02_kind_phys, 2.071840e-02_kind_phys,  & !2
       2.006702e-02_kind_phys, 1.943856e-02_kind_phys, 1.883064e-02_kind_phys, 1.824120e-02_kind_phys, 1.766849e-02_kind_phys,  & !2
       1.711099e-02_kind_phys, 1.656737e-02_kind_phys, 1.603647e-02_kind_phys, 1.551727e-02_kind_phys, 1.500886e-02_kind_phys,  & !2
       1.451045e-02_kind_phys, 1.402132e-02_kind_phys, 1.354084e-02_kind_phys, 1.306842e-02_kind_phys, 1.260355e-02_kind_phys,  & !2
       1.214575e-02_kind_phys, 1.169460e-02_kind_phys, 1.124971e-02_kind_phys, 1.081072e-02_kind_phys, 1.037731e-02_kind_phys,  & !2
       9.949167e-03_kind_phys, 9.526021e-03_kind_phys, 9.107615e-03_kind_phys, 8.693714e-03_kind_phys, 8.284096e-03_kind_phys,  & !2
       7.878558e-03_kind_phys, 7.476910e-03_kind_phys, 7.078974e-03_kind_phys, 6.684586e-03_kind_phys, 6.293589e-03_kind_phys,  & !2
       5.905839e-03_kind_phys, 5.521200e-03_kind_phys, 5.139543e-03_kind_phys,                                                  & !2
       1.065397e-01_kind_phys, 8.005726e-02_kind_phys, 6.546428e-02_kind_phys, 5.589131e-02_kind_phys, 4.898681e-02_kind_phys,  & !3
       4.369932e-02_kind_phys, 3.947901e-02_kind_phys, 3.600676e-02_kind_phys, 3.308299e-02_kind_phys, 3.057561e-02_kind_phys,  & !3
       2.839325e-02_kind_phys, 2.647040e-02_kind_phys, 2.475872e-02_kind_phys, 2.322164e-02_kind_phys, 2.183091e-02_kind_phys,  & !3
       2.056430e-02_kind_phys, 1.940407e-02_kind_phys, 1.833586e-02_kind_phys, 1.734787e-02_kind_phys, 1.643034e-02_kind_phys,  & !3
       1.557512e-02_kind_phys, 1.477530e-02_kind_phys, 1.402501e-02_kind_phys, 1.331924e-02_kind_phys, 1.265364e-02_kind_phys,  & !3
       1.202445e-02_kind_phys, 1.142838e-02_kind_phys, 1.086257e-02_kind_phys, 1.032445e-02_kind_phys, 9.811791e-03_kind_phys,  & !3
       9.322587e-03_kind_phys, 8.855053e-03_kind_phys, 8.407591e-03_kind_phys, 7.978763e-03_kind_phys, 7.567273e-03_kind_phys,  & !3
       7.171949e-03_kind_phys, 6.791728e-03_kind_phys, 6.425642e-03_kind_phys, 6.072809e-03_kind_phys, 5.732424e-03_kind_phys,  & !3
       5.403748e-03_kind_phys, 5.086103e-03_kind_phys, 4.778865e-03_kind_phys,                                                  & !3
       1.804566e-01_kind_phys, 1.168987e-01_kind_phys, 8.680442e-02_kind_phys, 6.910060e-02_kind_phys, 5.738174e-02_kind_phys,  & !4
       4.902332e-02_kind_phys, 4.274585e-02_kind_phys, 3.784923e-02_kind_phys, 3.391734e-02_kind_phys, 3.068690e-02_kind_phys,  & !4
       2.798301e-02_kind_phys, 2.568480e-02_kind_phys, 2.370600e-02_kind_phys, 2.198337e-02_kind_phys, 2.046940e-02_kind_phys,  & !4
       1.912777e-02_kind_phys, 1.793016e-02_kind_phys, 1.685420e-02_kind_phys, 1.588193e-02_kind_phys, 1.499882e-02_kind_phys,  & !4
       1.419293e-02_kind_phys, 1.345440e-02_kind_phys, 1.277496e-02_kind_phys, 1.214769e-02_kind_phys, 1.156669e-02_kind_phys,  & !4
       1.102694e-02_kind_phys, 1.052412e-02_kind_phys, 1.005451e-02_kind_phys, 9.614854e-03_kind_phys, 9.202335e-03_kind_phys,  & !4
       8.814470e-03_kind_phys, 8.449077e-03_kind_phys, 8.104223e-03_kind_phys, 7.778195e-03_kind_phys, 7.469466e-03_kind_phys,  & !4
       7.176671e-03_kind_phys, 6.898588e-03_kind_phys, 6.634117e-03_kind_phys, 6.382264e-03_kind_phys, 6.142134e-03_kind_phys,  & !4
       5.912913e-03_kind_phys, 5.693862e-03_kind_phys, 5.484308e-03_kind_phys,                                                  & !4
       2.131806e-01_kind_phys, 1.311372e-01_kind_phys, 9.407171e-02_kind_phys, 7.299442e-02_kind_phys, 5.941273e-02_kind_phys,  & !5
       4.994043e-02_kind_phys, 4.296242e-02_kind_phys, 3.761113e-02_kind_phys, 3.337910e-02_kind_phys, 2.994978e-02_kind_phys,  & !5
       2.711556e-02_kind_phys, 2.473461e-02_kind_phys, 2.270681e-02_kind_phys, 2.095943e-02_kind_phys, 1.943839e-02_kind_phys,  & !5
       1.810267e-02_kind_phys, 1.692057e-02_kind_phys, 1.586719e-02_kind_phys, 1.492275e-02_kind_phys, 1.407132e-02_kind_phys,  & !5
       1.329989e-02_kind_phys, 1.259780e-02_kind_phys, 1.195618e-02_kind_phys, 1.136761e-02_kind_phys, 1.082583e-02_kind_phys,  & !5
       1.032552e-02_kind_phys, 9.862158e-03_kind_phys, 9.431827e-03_kind_phys, 9.031157e-03_kind_phys, 8.657217e-03_kind_phys,  & !5
       8.307449e-03_kind_phys, 7.979609e-03_kind_phys, 7.671724e-03_kind_phys, 7.382048e-03_kind_phys, 7.109032e-03_kind_phys,  & !5
       6.851298e-03_kind_phys, 6.607615e-03_kind_phys, 6.376881e-03_kind_phys, 6.158105e-03_kind_phys, 5.950394e-03_kind_phys,  & !5
       5.752942e-03_kind_phys, 5.565019e-03_kind_phys, 5.385963e-03_kind_phys,                                                  & !5
       1.546177e-01_kind_phys, 1.039251e-01_kind_phys, 7.910347e-02_kind_phys, 6.412429e-02_kind_phys, 5.399997e-02_kind_phys,  & !6
       4.664937e-02_kind_phys, 4.104237e-02_kind_phys, 3.660781e-02_kind_phys, 3.300218e-02_kind_phys, 3.000586e-02_kind_phys,  & !6
       2.747148e-02_kind_phys, 2.529633e-02_kind_phys, 2.340647e-02_kind_phys, 2.174723e-02_kind_phys, 2.027731e-02_kind_phys,  & !6
       1.896487e-02_kind_phys, 1.778492e-02_kind_phys, 1.671761e-02_kind_phys, 1.574692e-02_kind_phys, 1.485978e-02_kind_phys,  & !6
       1.404543e-02_kind_phys, 1.329489e-02_kind_phys, 1.260066e-02_kind_phys, 1.195636e-02_kind_phys, 1.135657e-02_kind_phys,  & !6
       1.079664e-02_kind_phys, 1.027257e-02_kind_phys, 9.780871e-03_kind_phys, 9.318505e-03_kind_phys, 8.882815e-03_kind_phys,  & !6
       8.471458e-03_kind_phys, 8.082364e-03_kind_phys, 7.713696e-03_kind_phys, 7.363817e-03_kind_phys, 7.031264e-03_kind_phys,  & !6
       6.714725e-03_kind_phys, 6.413021e-03_kind_phys, 6.125086e-03_kind_phys, 5.849958e-03_kind_phys, 5.586764e-03_kind_phys,  & !6
       5.334707e-03_kind_phys, 5.093066e-03_kind_phys, 4.861179e-03_kind_phys,                                                  & !6
       7.583404e-02_kind_phys, 6.181558e-02_kind_phys, 5.312027e-02_kind_phys, 4.696039e-02_kind_phys, 4.225986e-02_kind_phys,  & !7
       3.849735e-02_kind_phys, 3.538340e-02_kind_phys, 3.274182e-02_kind_phys, 3.045798e-02_kind_phys, 2.845343e-02_kind_phys,  & !7
       2.667231e-02_kind_phys, 2.507353e-02_kind_phys, 2.362606e-02_kind_phys, 2.230595e-02_kind_phys, 2.109435e-02_kind_phys,  & !7
       1.997617e-02_kind_phys, 1.893916e-02_kind_phys, 1.797328e-02_kind_phys, 1.707016e-02_kind_phys, 1.622279e-02_kind_phys,  & !7
       1.542523e-02_kind_phys, 1.467241e-02_kind_phys, 1.395997e-02_kind_phys, 1.328414e-02_kind_phys, 1.264164e-02_kind_phys,  & !7
       1.202958e-02_kind_phys, 1.144544e-02_kind_phys, 1.088697e-02_kind_phys, 1.035218e-02_kind_phys, 9.839297e-03_kind_phys,  & !7
       9.346733e-03_kind_phys, 8.873057e-03_kind_phys, 8.416980e-03_kind_phys, 7.977335e-03_kind_phys, 7.553066e-03_kind_phys,  & !7
       7.143210e-03_kind_phys, 6.746888e-03_kind_phys, 6.363297e-03_kind_phys, 5.991700e-03_kind_phys, 5.631422e-03_kind_phys,  & !7
       5.281840e-03_kind_phys, 4.942378e-03_kind_phys, 4.612505e-03_kind_phys,                                                  & !7
       9.022185e-02_kind_phys, 6.922700e-02_kind_phys, 5.710674e-02_kind_phys, 4.898377e-02_kind_phys, 4.305946e-02_kind_phys,  & !8
       3.849553e-02_kind_phys, 3.484183e-02_kind_phys, 3.183220e-02_kind_phys, 2.929794e-02_kind_phys, 2.712627e-02_kind_phys,  & !8
       2.523856e-02_kind_phys, 2.357810e-02_kind_phys, 2.210286e-02_kind_phys, 2.078089e-02_kind_phys, 1.958747e-02_kind_phys,  & !8
       1.850310e-02_kind_phys, 1.751218e-02_kind_phys, 1.660205e-02_kind_phys, 1.576232e-02_kind_phys, 1.498440e-02_kind_phys,  & !8
       1.426107e-02_kind_phys, 1.358624e-02_kind_phys, 1.295474e-02_kind_phys, 1.236212e-02_kind_phys, 1.180456e-02_kind_phys,  & !8
       1.127874e-02_kind_phys, 1.078175e-02_kind_phys, 1.031106e-02_kind_phys, 9.864433e-03_kind_phys, 9.439878e-03_kind_phys,  & !8
       9.035637e-03_kind_phys, 8.650140e-03_kind_phys, 8.281981e-03_kind_phys, 7.929895e-03_kind_phys, 7.592746e-03_kind_phys,  & !8
       7.269505e-03_kind_phys, 6.959238e-03_kind_phys, 6.661100e-03_kind_phys, 6.374317e-03_kind_phys, 6.098185e-03_kind_phys,  & !8
       5.832059e-03_kind_phys, 5.575347e-03_kind_phys, 5.327504e-03_kind_phys,                                                  & !8
       1.294087e-01_kind_phys, 8.788217e-02_kind_phys, 6.728288e-02_kind_phys, 5.479720e-02_kind_phys, 4.635049e-02_kind_phys,  & !9
       4.022253e-02_kind_phys, 3.555576e-02_kind_phys, 3.187259e-02_kind_phys, 2.888498e-02_kind_phys, 2.640843e-02_kind_phys,  & !9
       2.431904e-02_kind_phys, 2.253038e-02_kind_phys, 2.098024e-02_kind_phys, 1.962267e-02_kind_phys, 1.842293e-02_kind_phys,  & !9
       1.735426e-02_kind_phys, 1.639571e-02_kind_phys, 1.553060e-02_kind_phys, 1.474552e-02_kind_phys, 1.402953e-02_kind_phys,  & !9
       1.337363e-02_kind_phys, 1.277033e-02_kind_phys, 1.221336e-02_kind_phys, 1.169741e-02_kind_phys, 1.121797e-02_kind_phys,  & !9
       1.077117e-02_kind_phys, 1.035369e-02_kind_phys, 9.962643e-03_kind_phys, 9.595509e-03_kind_phys, 9.250088e-03_kind_phys,  & !9
       8.924447e-03_kind_phys, 8.616876e-03_kind_phys, 8.325862e-03_kind_phys, 8.050057e-03_kind_phys, 7.788258e-03_kind_phys,  & !9
       7.539388e-03_kind_phys, 7.302478e-03_kind_phys, 7.076656e-03_kind_phys, 6.861134e-03_kind_phys, 6.655197e-03_kind_phys,  & !9
       6.458197e-03_kind_phys, 6.269543e-03_kind_phys, 6.088697e-03_kind_phys,                                                  & !9
       1.593628e-01_kind_phys, 1.014552e-01_kind_phys, 7.458955e-02_kind_phys, 5.903571e-02_kind_phys, 4.887582e-02_kind_phys,  & !10
       4.171159e-02_kind_phys, 3.638480e-02_kind_phys, 3.226692e-02_kind_phys, 2.898717e-02_kind_phys, 2.631256e-02_kind_phys,  & !10
       2.408925e-02_kind_phys, 2.221156e-02_kind_phys, 2.060448e-02_kind_phys, 1.921325e-02_kind_phys, 1.799699e-02_kind_phys,  & !10
       1.692456e-02_kind_phys, 1.597177e-02_kind_phys, 1.511961e-02_kind_phys, 1.435289e-02_kind_phys, 1.365933e-02_kind_phys,  & !10
       1.302890e-02_kind_phys, 1.245334e-02_kind_phys, 1.192576e-02_kind_phys, 1.144037e-02_kind_phys, 1.099230e-02_kind_phys,  & !10
       1.057739e-02_kind_phys, 1.019208e-02_kind_phys, 9.833302e-03_kind_phys, 9.498395e-03_kind_phys, 9.185047e-03_kind_phys,  & !10
       8.891237e-03_kind_phys, 8.615185e-03_kind_phys, 8.355325e-03_kind_phys, 8.110267e-03_kind_phys, 7.878778e-03_kind_phys,  & !10
       7.659759e-03_kind_phys, 7.452224e-03_kind_phys, 7.255291e-03_kind_phys, 7.068166e-03_kind_phys, 6.890130e-03_kind_phys,  & !10
       6.720536e-03_kind_phys, 6.558794e-03_kind_phys, 6.404371e-03_kind_phys,                                                  & !10
       1.656227e-01_kind_phys, 1.032129e-01_kind_phys, 7.487359e-02_kind_phys, 5.871431e-02_kind_phys, 4.828355e-02_kind_phys,  & !11
       4.099989e-02_kind_phys, 3.562924e-02_kind_phys, 3.150755e-02_kind_phys, 2.824593e-02_kind_phys, 2.560156e-02_kind_phys,  & !11
       2.341503e-02_kind_phys, 2.157740e-02_kind_phys, 2.001169e-02_kind_phys, 1.866199e-02_kind_phys, 1.748669e-02_kind_phys,  & !11
       1.645421e-02_kind_phys, 1.554015e-02_kind_phys, 1.472535e-02_kind_phys, 1.399457e-02_kind_phys, 1.333553e-02_kind_phys,  & !11
       1.273821e-02_kind_phys, 1.219440e-02_kind_phys, 1.169725e-02_kind_phys, 1.124104e-02_kind_phys, 1.082096e-02_kind_phys,  & !11
       1.043290e-02_kind_phys, 1.007336e-02_kind_phys, 9.739338e-03_kind_phys, 9.428223e-03_kind_phys, 9.137756e-03_kind_phys,  & !11
       8.865964e-03_kind_phys, 8.611115e-03_kind_phys, 8.371686e-03_kind_phys, 8.146330e-03_kind_phys, 7.933852e-03_kind_phys,  & !11
       7.733187e-03_kind_phys, 7.543386e-03_kind_phys, 7.363597e-03_kind_phys, 7.193056e-03_kind_phys, 7.031072e-03_kind_phys,  & !11
       6.877024e-03_kind_phys, 6.730348e-03_kind_phys, 6.590531e-03_kind_phys,                                                  & !11
       9.194591e-02_kind_phys, 6.446867e-02_kind_phys, 4.962034e-02_kind_phys, 4.042061e-02_kind_phys, 3.418456e-02_kind_phys,  & !12
       2.968856e-02_kind_phys, 2.629900e-02_kind_phys, 2.365572e-02_kind_phys, 2.153915e-02_kind_phys, 1.980791e-02_kind_phys,  & !12
       1.836689e-02_kind_phys, 1.714979e-02_kind_phys, 1.610900e-02_kind_phys, 1.520946e-02_kind_phys, 1.442476e-02_kind_phys,  & !12
       1.373468e-02_kind_phys, 1.312345e-02_kind_phys, 1.257858e-02_kind_phys, 1.209010e-02_kind_phys, 1.164990e-02_kind_phys,  & !12
       1.125136e-02_kind_phys, 1.088901e-02_kind_phys, 1.055827e-02_kind_phys, 1.025531e-02_kind_phys, 9.976896e-03_kind_phys,  & !12
       9.720255e-03_kind_phys, 9.483022e-03_kind_phys, 9.263160e-03_kind_phys, 9.058902e-03_kind_phys, 8.868710e-03_kind_phys,  & !12
       8.691240e-03_kind_phys, 8.525312e-03_kind_phys, 8.369886e-03_kind_phys, 8.224042e-03_kind_phys, 8.086961e-03_kind_phys,  & !12
       7.957917e-03_kind_phys, 7.836258e-03_kind_phys, 7.721400e-03_kind_phys, 7.612821e-03_kind_phys, 7.510045e-03_kind_phys,  & !12
       7.412648e-03_kind_phys, 7.320242e-03_kind_phys, 7.232476e-03_kind_phys,                                                  & !12
       1.437021e-01_kind_phys, 8.872535e-02_kind_phys, 6.392420e-02_kind_phys, 4.991833e-02_kind_phys, 4.096790e-02_kind_phys,  & !13
       3.477881e-02_kind_phys, 3.025782e-02_kind_phys, 2.681909e-02_kind_phys, 2.412102e-02_kind_phys, 2.195132e-02_kind_phys,  & !13
       2.017124e-02_kind_phys, 1.868641e-02_kind_phys, 1.743044e-02_kind_phys, 1.635529e-02_kind_phys, 1.542540e-02_kind_phys,  & !13
       1.461388e-02_kind_phys, 1.390003e-02_kind_phys, 1.326766e-02_kind_phys, 1.270395e-02_kind_phys, 1.219860e-02_kind_phys,  & !13
       1.174326e-02_kind_phys, 1.133107e-02_kind_phys, 1.095637e-02_kind_phys, 1.061442e-02_kind_phys, 1.030126e-02_kind_phys,  & !13
       1.001352e-02_kind_phys, 9.748340e-03_kind_phys, 9.503256e-03_kind_phys, 9.276155e-03_kind_phys, 9.065205e-03_kind_phys,  & !13
       8.868808e-03_kind_phys, 8.685571e-03_kind_phys, 8.514268e-03_kind_phys, 8.353820e-03_kind_phys, 8.203272e-03_kind_phys,  & !13
       8.061776e-03_kind_phys, 7.928578e-03_kind_phys, 7.803001e-03_kind_phys, 7.684443e-03_kind_phys, 7.572358e-03_kind_phys,  & !13
       7.466258e-03_kind_phys, 7.365701e-03_kind_phys, 7.270286e-03_kind_phys,                                                  & !13
       1.288870e-01_kind_phys, 8.160295e-02_kind_phys, 5.964745e-02_kind_phys, 4.703790e-02_kind_phys, 3.888637e-02_kind_phys,  & !14
       3.320115e-02_kind_phys, 2.902017e-02_kind_phys, 2.582259e-02_kind_phys, 2.330224e-02_kind_phys, 2.126754e-02_kind_phys,  & !14
       1.959258e-02_kind_phys, 1.819130e-02_kind_phys, 1.700289e-02_kind_phys, 1.598320e-02_kind_phys, 1.509942e-02_kind_phys,  & !14
       1.432666e-02_kind_phys, 1.364572e-02_kind_phys, 1.304156e-02_kind_phys, 1.250220e-02_kind_phys, 1.201803e-02_kind_phys,  & !14
       1.158123e-02_kind_phys, 1.118537e-02_kind_phys, 1.082513e-02_kind_phys, 1.049605e-02_kind_phys, 1.019440e-02_kind_phys,  & !14
       9.916989e-03_kind_phys, 9.661116e-03_kind_phys, 9.424457e-03_kind_phys, 9.205005e-03_kind_phys, 9.001022e-03_kind_phys,  & !14
       8.810992e-03_kind_phys, 8.633588e-03_kind_phys, 8.467646e-03_kind_phys, 8.312137e-03_kind_phys, 8.166151e-03_kind_phys,  & !14
       8.028878e-03_kind_phys, 7.899597e-03_kind_phys, 7.777663e-03_kind_phys, 7.662498e-03_kind_phys, 7.553581e-03_kind_phys,  & !14
       7.450444e-03_kind_phys, 7.352662e-03_kind_phys, 7.259851e-03_kind_phys,                                                  & !14
       8.254229e-02_kind_phys, 5.808787e-02_kind_phys, 4.492166e-02_kind_phys, 3.675028e-02_kind_phys, 3.119623e-02_kind_phys,  & !15
       2.718045e-02_kind_phys, 2.414450e-02_kind_phys, 2.177073e-02_kind_phys, 1.986526e-02_kind_phys, 1.830306e-02_kind_phys,  & !15
       1.699991e-02_kind_phys, 1.589698e-02_kind_phys, 1.495199e-02_kind_phys, 1.413374e-02_kind_phys, 1.341870e-02_kind_phys,  & !15
       1.278883e-02_kind_phys, 1.223002e-02_kind_phys, 1.173114e-02_kind_phys, 1.128322e-02_kind_phys, 1.087900e-02_kind_phys,  & !15
       1.051254e-02_kind_phys, 1.017890e-02_kind_phys, 9.873991e-03_kind_phys, 9.594347e-03_kind_phys, 9.337044e-03_kind_phys,  & !15
       9.099589e-03_kind_phys, 8.879842e-03_kind_phys, 8.675960e-03_kind_phys, 8.486341e-03_kind_phys, 8.309594e-03_kind_phys,  & !15
       8.144500e-03_kind_phys, 7.989986e-03_kind_phys, 7.845109e-03_kind_phys, 7.709031e-03_kind_phys, 7.581007e-03_kind_phys,  & !15
       7.460376e-03_kind_phys, 7.346544e-03_kind_phys, 7.238978e-03_kind_phys, 7.137201e-03_kind_phys, 7.040780e-03_kind_phys,  & !15
       6.949325e-03_kind_phys, 6.862483e-03_kind_phys, 6.779931e-03_kind_phys,                                                  & !15
       1.382062e-01_kind_phys, 8.643227e-02_kind_phys, 6.282935e-02_kind_phys, 4.934783e-02_kind_phys, 4.063891e-02_kind_phys,  & !16
       3.455591e-02_kind_phys, 3.007059e-02_kind_phys, 2.662897e-02_kind_phys, 2.390631e-02_kind_phys, 2.169972e-02_kind_phys,  & !16
       1.987596e-02_kind_phys, 1.834393e-02_kind_phys, 1.703924e-02_kind_phys, 1.591513e-02_kind_phys, 1.493679e-02_kind_phys,  & !16
       1.407780e-02_kind_phys, 1.331775e-02_kind_phys, 1.264061e-02_kind_phys, 1.203364e-02_kind_phys, 1.148655e-02_kind_phys,  & !16
       1.099099e-02_kind_phys, 1.054006e-02_kind_phys, 1.012807e-02_kind_phys, 9.750215e-03_kind_phys, 9.402477e-03_kind_phys,  & !16
       9.081428e-03_kind_phys, 8.784143e-03_kind_phys, 8.508107e-03_kind_phys, 8.251146e-03_kind_phys, 8.011373e-03_kind_phys,  & !16
       7.787140e-03_kind_phys, 7.577002e-03_kind_phys, 7.379687e-03_kind_phys, 7.194071e-03_kind_phys, 7.019158e-03_kind_phys,  & !16
       6.854061e-03_kind_phys, 6.697986e-03_kind_phys, 6.550224e-03_kind_phys, 6.410138e-03_kind_phys, 6.277153e-03_kind_phys,  & !16
       6.150751e-03_kind_phys, 6.030462e-03_kind_phys, 5.915860e-03_kind_phys/), &                                                !16
       shape=(/43,nBandsLW_RRTMG/))
  
  real(kind_phys) , dimension(46,nBandsLW_RRTMG),parameter :: &
       absice3 =  reshape(source=(/ &
       3.110649e-03_kind_phys, 4.666352e-02_kind_phys, 6.606447e-02_kind_phys, 6.531678e-02_kind_phys, 6.012598e-02_kind_phys,  & !1
       5.437494e-02_kind_phys, 4.906411e-02_kind_phys, 4.441146e-02_kind_phys, 4.040585e-02_kind_phys, 3.697334e-02_kind_phys,  & !1
       3.403027e-02_kind_phys, 3.149979e-02_kind_phys, 2.931596e-02_kind_phys, 2.742365e-02_kind_phys, 2.577721e-02_kind_phys,  & !1
       2.433888e-02_kind_phys, 2.307732e-02_kind_phys, 2.196644e-02_kind_phys, 2.098437e-02_kind_phys, 2.011264e-02_kind_phys,  & !1
       1.933561e-02_kind_phys, 1.863992e-02_kind_phys, 1.801407e-02_kind_phys, 1.744812e-02_kind_phys, 1.693346e-02_kind_phys,  & !1
       1.646252e-02_kind_phys, 1.602866e-02_kind_phys, 1.562600e-02_kind_phys, 1.524933e-02_kind_phys, 1.489399e-02_kind_phys,  & !1
       1.455580e-02_kind_phys, 1.423098e-02_kind_phys, 1.391612e-02_kind_phys, 1.360812e-02_kind_phys, 1.330413e-02_kind_phys,  & !1
       1.300156e-02_kind_phys, 1.269801e-02_kind_phys, 1.239127e-02_kind_phys, 1.207928e-02_kind_phys, 1.176014e-02_kind_phys,  & !1
       1.143204e-02_kind_phys, 1.109334e-02_kind_phys, 1.074243e-02_kind_phys, 1.037786e-02_kind_phys, 9.998198e-03_kind_phys,  & !1
       9.602126e-03_kind_phys,                                                                                                  & !1
       3.984966e-04_kind_phys, 1.681097e-02_kind_phys, 2.627680e-02_kind_phys, 2.767465e-02_kind_phys, 2.700722e-02_kind_phys,  & !2
       2.579180e-02_kind_phys, 2.448677e-02_kind_phys, 2.323890e-02_kind_phys, 2.209096e-02_kind_phys, 2.104882e-02_kind_phys,  & !2
       2.010547e-02_kind_phys, 1.925003e-02_kind_phys, 1.847128e-02_kind_phys, 1.775883e-02_kind_phys, 1.710358e-02_kind_phys,  & !2
       1.649769e-02_kind_phys, 1.593449e-02_kind_phys, 1.540829e-02_kind_phys, 1.491429e-02_kind_phys, 1.444837e-02_kind_phys,  & !2
       1.400704e-02_kind_phys, 1.358729e-02_kind_phys, 1.318654e-02_kind_phys, 1.280258e-02_kind_phys, 1.243346e-02_kind_phys,  & !2
       1.207750e-02_kind_phys, 1.173325e-02_kind_phys, 1.139941e-02_kind_phys, 1.107487e-02_kind_phys, 1.075861e-02_kind_phys,  & !2
       1.044975e-02_kind_phys, 1.014753e-02_kind_phys, 9.851229e-03_kind_phys, 9.560240e-03_kind_phys, 9.274003e-03_kind_phys,  & !2
       8.992020e-03_kind_phys, 8.713845e-03_kind_phys, 8.439074e-03_kind_phys, 8.167346e-03_kind_phys, 7.898331e-03_kind_phys,  & !2
       7.631734e-03_kind_phys, 7.367286e-03_kind_phys, 7.104742e-03_kind_phys, 6.843882e-03_kind_phys, 6.584504e-03_kind_phys,  & !2
       6.326424e-03_kind_phys,                                                                                                  & !2
       6.933163e-02_kind_phys, 8.540475e-02_kind_phys, 7.701816e-02_kind_phys, 6.771158e-02_kind_phys, 5.986953e-02_kind_phys,  & !3
       5.348120e-02_kind_phys, 4.824962e-02_kind_phys, 4.390563e-02_kind_phys, 4.024411e-02_kind_phys, 3.711404e-02_kind_phys,  & !3
       3.440426e-02_kind_phys, 3.203200e-02_kind_phys, 2.993478e-02_kind_phys, 2.806474e-02_kind_phys, 2.638464e-02_kind_phys,  & !3
       2.486516e-02_kind_phys, 2.348288e-02_kind_phys, 2.221890e-02_kind_phys, 2.105780e-02_kind_phys, 1.998687e-02_kind_phys,  & !3
       1.899552e-02_kind_phys, 1.807490e-02_kind_phys, 1.721750e-02_kind_phys, 1.641693e-02_kind_phys, 1.566773e-02_kind_phys,  & !3
       1.496515e-02_kind_phys, 1.430509e-02_kind_phys, 1.368398e-02_kind_phys, 1.309865e-02_kind_phys, 1.254634e-02_kind_phys,  & !3
       1.202456e-02_kind_phys, 1.153114e-02_kind_phys, 1.106409e-02_kind_phys, 1.062166e-02_kind_phys, 1.020224e-02_kind_phys,  & !3
       9.804381e-03_kind_phys, 9.426771e-03_kind_phys, 9.068205e-03_kind_phys, 8.727578e-03_kind_phys, 8.403876e-03_kind_phys,  & !3
       8.096160e-03_kind_phys, 7.803564e-03_kind_phys, 7.525281e-03_kind_phys, 7.260560e-03_kind_phys, 7.008697e-03_kind_phys,  & !3
       6.769036e-03_kind_phys,                                                                                                  & !3
       1.765735e-01_kind_phys, 1.382700e-01_kind_phys, 1.095129e-01_kind_phys, 8.987475e-02_kind_phys, 7.591185e-02_kind_phys,  & !4
       6.554169e-02_kind_phys, 5.755500e-02_kind_phys, 5.122083e-02_kind_phys, 4.607610e-02_kind_phys, 4.181475e-02_kind_phys,  & !4
       3.822697e-02_kind_phys, 3.516432e-02_kind_phys, 3.251897e-02_kind_phys, 3.021073e-02_kind_phys, 2.817876e-02_kind_phys,  & !4
       2.637607e-02_kind_phys, 2.476582e-02_kind_phys, 2.331871e-02_kind_phys, 2.201113e-02_kind_phys, 2.082388e-02_kind_phys,  & !4
       1.974115e-02_kind_phys, 1.874983e-02_kind_phys, 1.783894e-02_kind_phys, 1.699922e-02_kind_phys, 1.622280e-02_kind_phys,  & !4
       1.550296e-02_kind_phys, 1.483390e-02_kind_phys, 1.421064e-02_kind_phys, 1.362880e-02_kind_phys, 1.308460e-02_kind_phys,  & !4
       1.257468e-02_kind_phys, 1.209611e-02_kind_phys, 1.164628e-02_kind_phys, 1.122287e-02_kind_phys, 1.082381e-02_kind_phys,  & !4
       1.044725e-02_kind_phys, 1.009154e-02_kind_phys, 9.755166e-03_kind_phys, 9.436783e-03_kind_phys, 9.135163e-03_kind_phys,  & !4
       8.849193e-03_kind_phys, 8.577856e-03_kind_phys, 8.320225e-03_kind_phys, 8.075451e-03_kind_phys, 7.842755e-03_kind_phys,  & !4
       7.621418e-03_kind_phys,                                                                                                  & !4
       2.339673e-01_kind_phys, 1.692124e-01_kind_phys, 1.291656e-01_kind_phys, 1.033837e-01_kind_phys, 8.562949e-02_kind_phys,  & !5
       7.273526e-02_kind_phys, 6.298262e-02_kind_phys, 5.537015e-02_kind_phys, 4.927787e-02_kind_phys, 4.430246e-02_kind_phys,  & !5
       4.017061e-02_kind_phys, 3.669072e-02_kind_phys, 3.372455e-02_kind_phys, 3.116995e-02_kind_phys, 2.894977e-02_kind_phys,  & !5
       2.700471e-02_kind_phys, 2.528842e-02_kind_phys, 2.376420e-02_kind_phys, 2.240256e-02_kind_phys, 2.117959e-02_kind_phys,  & !5
       2.007567e-02_kind_phys, 1.907456e-02_kind_phys, 1.816271e-02_kind_phys, 1.732874e-02_kind_phys, 1.656300e-02_kind_phys,  & !5
       1.585725e-02_kind_phys, 1.520445e-02_kind_phys, 1.459852e-02_kind_phys, 1.403419e-02_kind_phys, 1.350689e-02_kind_phys,  & !5
       1.301260e-02_kind_phys, 1.254781e-02_kind_phys, 1.210941e-02_kind_phys, 1.169468e-02_kind_phys, 1.130118e-02_kind_phys,  & !5
       1.092675e-02_kind_phys, 1.056945e-02_kind_phys, 1.022757e-02_kind_phys, 9.899560e-03_kind_phys, 9.584021e-03_kind_phys,  & !5
       9.279705e-03_kind_phys, 8.985479e-03_kind_phys, 8.700322e-03_kind_phys, 8.423306e-03_kind_phys, 8.153590e-03_kind_phys,  & !5
       7.890412e-03_kind_phys,                                                                                                  & !5
       1.145369e-01_kind_phys, 1.174566e-01_kind_phys, 9.917866e-02_kind_phys, 8.332990e-02_kind_phys, 7.104263e-02_kind_phys,  & !6
       6.153370e-02_kind_phys, 5.405472e-02_kind_phys, 4.806281e-02_kind_phys, 4.317918e-02_kind_phys, 3.913795e-02_kind_phys,  & !6
       3.574916e-02_kind_phys, 3.287437e-02_kind_phys, 3.041067e-02_kind_phys, 2.828017e-02_kind_phys, 2.642292e-02_kind_phys,  & !6
       2.479206e-02_kind_phys, 2.335051e-02_kind_phys, 2.206851e-02_kind_phys, 2.092195e-02_kind_phys, 1.989108e-02_kind_phys,  & !6
       1.895958e-02_kind_phys, 1.811385e-02_kind_phys, 1.734245e-02_kind_phys, 1.663573e-02_kind_phys, 1.598545e-02_kind_phys,  & !6
       1.538456e-02_kind_phys, 1.482700e-02_kind_phys, 1.430750e-02_kind_phys, 1.382150e-02_kind_phys, 1.336499e-02_kind_phys,  & !6
       1.293447e-02_kind_phys, 1.252685e-02_kind_phys, 1.213939e-02_kind_phys, 1.176968e-02_kind_phys, 1.141555e-02_kind_phys,  & !6
       1.107508e-02_kind_phys, 1.074655e-02_kind_phys, 1.042839e-02_kind_phys, 1.011923e-02_kind_phys, 9.817799e-03_kind_phys,  & !6
       9.522962e-03_kind_phys, 9.233688e-03_kind_phys, 8.949041e-03_kind_phys, 8.668171e-03_kind_phys, 8.390301e-03_kind_phys,  & !6
       8.114723e-03_kind_phys,                                                                                                  & !6
       1.222345e-02_kind_phys, 5.344230e-02_kind_phys, 5.523465e-02_kind_phys, 5.128759e-02_kind_phys, 4.676925e-02_kind_phys,  & !7
       4.266150e-02_kind_phys, 3.910561e-02_kind_phys, 3.605479e-02_kind_phys, 3.342843e-02_kind_phys, 3.115052e-02_kind_phys,  & !7
       2.915776e-02_kind_phys, 2.739935e-02_kind_phys, 2.583499e-02_kind_phys, 2.443266e-02_kind_phys, 2.316681e-02_kind_phys,  & !7
       2.201687e-02_kind_phys, 2.096619e-02_kind_phys, 2.000112e-02_kind_phys, 1.911044e-02_kind_phys, 1.828481e-02_kind_phys,  & !7
       1.751641e-02_kind_phys, 1.679866e-02_kind_phys, 1.612598e-02_kind_phys, 1.549360e-02_kind_phys, 1.489742e-02_kind_phys,  & !7
       1.433392e-02_kind_phys, 1.380002e-02_kind_phys, 1.329305e-02_kind_phys, 1.281068e-02_kind_phys, 1.235084e-02_kind_phys,  & !7
       1.191172e-02_kind_phys, 1.149171e-02_kind_phys, 1.108936e-02_kind_phys, 1.070341e-02_kind_phys, 1.033271e-02_kind_phys,  & !7
       9.976220e-03_kind_phys, 9.633021e-03_kind_phys, 9.302273e-03_kind_phys, 8.983216e-03_kind_phys, 8.675161e-03_kind_phys,  & !7
       8.377478e-03_kind_phys, 8.089595e-03_kind_phys, 7.810986e-03_kind_phys, 7.541170e-03_kind_phys, 7.279706e-03_kind_phys,  & !7
       7.026186e-03_kind_phys,                                                                                                  & !7
       6.711058e-02_kind_phys, 6.918198e-02_kind_phys, 6.127484e-02_kind_phys, 5.411944e-02_kind_phys, 4.836902e-02_kind_phys,  & !8
       4.375293e-02_kind_phys, 3.998077e-02_kind_phys, 3.683587e-02_kind_phys, 3.416508e-02_kind_phys, 3.186003e-02_kind_phys,  & !8
       2.984290e-02_kind_phys, 2.805671e-02_kind_phys, 2.645895e-02_kind_phys, 2.501733e-02_kind_phys, 2.370689e-02_kind_phys,  & !8
       2.250808e-02_kind_phys, 2.140532e-02_kind_phys, 2.038609e-02_kind_phys, 1.944018e-02_kind_phys, 1.855918e-02_kind_phys,  & !8
       1.773609e-02_kind_phys, 1.696504e-02_kind_phys, 1.624106e-02_kind_phys, 1.555990e-02_kind_phys, 1.491793e-02_kind_phys,  & !8
       1.431197e-02_kind_phys, 1.373928e-02_kind_phys, 1.319743e-02_kind_phys, 1.268430e-02_kind_phys, 1.219799e-02_kind_phys,  & !8
       1.173682e-02_kind_phys, 1.129925e-02_kind_phys, 1.088393e-02_kind_phys, 1.048961e-02_kind_phys, 1.011516e-02_kind_phys,  & !8
       9.759543e-03_kind_phys, 9.421813e-03_kind_phys, 9.101089e-03_kind_phys, 8.796559e-03_kind_phys, 8.507464e-03_kind_phys,  & !8
       8.233098e-03_kind_phys, 7.972798e-03_kind_phys, 7.725942e-03_kind_phys, 7.491940e-03_kind_phys, 7.270238e-03_kind_phys,  & !8
       7.060305e-03_kind_phys,                                                                                                  & !8
       1.236780e-01_kind_phys, 9.222386e-02_kind_phys, 7.383997e-02_kind_phys, 6.204072e-02_kind_phys, 5.381029e-02_kind_phys,  & !9
       4.770678e-02_kind_phys, 4.296928e-02_kind_phys, 3.916131e-02_kind_phys, 3.601540e-02_kind_phys, 3.335878e-02_kind_phys,  & !9
       3.107493e-02_kind_phys, 2.908247e-02_kind_phys, 2.732282e-02_kind_phys, 2.575276e-02_kind_phys, 2.433968e-02_kind_phys,  & !9
       2.305852e-02_kind_phys, 2.188966e-02_kind_phys, 2.081757e-02_kind_phys, 1.982974e-02_kind_phys, 1.891599e-02_kind_phys,  & !9
       1.806794e-02_kind_phys, 1.727865e-02_kind_phys, 1.654227e-02_kind_phys, 1.585387e-02_kind_phys, 1.520924e-02_kind_phys,  & !9
       1.460476e-02_kind_phys, 1.403730e-02_kind_phys, 1.350416e-02_kind_phys, 1.300293e-02_kind_phys, 1.253153e-02_kind_phys,  & !9
       1.208808e-02_kind_phys, 1.167094e-02_kind_phys, 1.127862e-02_kind_phys, 1.090979e-02_kind_phys, 1.056323e-02_kind_phys,  & !9
       1.023786e-02_kind_phys, 9.932665e-03_kind_phys, 9.646744e-03_kind_phys, 9.379250e-03_kind_phys, 9.129409e-03_kind_phys,  & !9
       8.896500e-03_kind_phys, 8.679856e-03_kind_phys, 8.478852e-03_kind_phys, 8.292904e-03_kind_phys, 8.121463e-03_kind_phys,  & !9
       7.964013e-03_kind_phys,                                                                                                  & !9
       1.655966e-01_kind_phys, 1.134205e-01_kind_phys, 8.714344e-02_kind_phys, 7.129241e-02_kind_phys, 6.063739e-02_kind_phys,  & !10
       5.294203e-02_kind_phys, 4.709309e-02_kind_phys, 4.247476e-02_kind_phys, 3.871892e-02_kind_phys, 3.559206e-02_kind_phys,  & !10
       3.293893e-02_kind_phys, 3.065226e-02_kind_phys, 2.865558e-02_kind_phys, 2.689288e-02_kind_phys, 2.532221e-02_kind_phys,  & !10
       2.391150e-02_kind_phys, 2.263582e-02_kind_phys, 2.147549e-02_kind_phys, 2.041476e-02_kind_phys, 1.944089e-02_kind_phys,  & !10
       1.854342e-02_kind_phys, 1.771371e-02_kind_phys, 1.694456e-02_kind_phys, 1.622989e-02_kind_phys, 1.556456e-02_kind_phys,  & !10
       1.494415e-02_kind_phys, 1.436491e-02_kind_phys, 1.382354e-02_kind_phys, 1.331719e-02_kind_phys, 1.284339e-02_kind_phys,  & !10
       1.239992e-02_kind_phys, 1.198486e-02_kind_phys, 1.159647e-02_kind_phys, 1.123323e-02_kind_phys, 1.089375e-02_kind_phys,  & !10
       1.057679e-02_kind_phys, 1.028124e-02_kind_phys, 1.000607e-02_kind_phys, 9.750376e-03_kind_phys, 9.513303e-03_kind_phys,  & !10
       9.294082e-03_kind_phys, 9.092003e-03_kind_phys, 8.906412e-03_kind_phys, 8.736702e-03_kind_phys, 8.582314e-03_kind_phys,  & !10
       8.442725e-03_kind_phys,                                                                                                  & !10
       1.775615e-01_kind_phys, 1.180046e-01_kind_phys, 8.929607e-02_kind_phys, 7.233500e-02_kind_phys, 6.108333e-02_kind_phys,  & !11
       5.303642e-02_kind_phys, 4.696927e-02_kind_phys, 4.221206e-02_kind_phys, 3.836768e-02_kind_phys, 3.518576e-02_kind_phys,  & !11
       3.250063e-02_kind_phys, 3.019825e-02_kind_phys, 2.819758e-02_kind_phys, 2.643943e-02_kind_phys, 2.487953e-02_kind_phys,  & !11
       2.348414e-02_kind_phys, 2.222705e-02_kind_phys, 2.108762e-02_kind_phys, 2.004936e-02_kind_phys, 1.909892e-02_kind_phys,  & !11
       1.822539e-02_kind_phys, 1.741975e-02_kind_phys, 1.667449e-02_kind_phys, 1.598330e-02_kind_phys, 1.534084e-02_kind_phys,  & !11
       1.474253e-02_kind_phys, 1.418446e-02_kind_phys, 1.366325e-02_kind_phys, 1.317597e-02_kind_phys, 1.272004e-02_kind_phys,  & !11
       1.229321e-02_kind_phys, 1.189350e-02_kind_phys, 1.151915e-02_kind_phys, 1.116859e-02_kind_phys, 1.084042e-02_kind_phys,  & !11
       1.053338e-02_kind_phys, 1.024636e-02_kind_phys, 9.978326e-03_kind_phys, 9.728357e-03_kind_phys, 9.495613e-03_kind_phys,  & !11
       9.279327e-03_kind_phys, 9.078798e-03_kind_phys, 8.893383e-03_kind_phys, 8.722488e-03_kind_phys, 8.565568e-03_kind_phys,  & !11
       8.422115e-03_kind_phys,                                                                                                  & !11
       9.465447e-02_kind_phys, 6.432047e-02_kind_phys, 5.060973e-02_kind_phys, 4.267283e-02_kind_phys, 3.741843e-02_kind_phys,  & !12
       3.363096e-02_kind_phys, 3.073531e-02_kind_phys, 2.842405e-02_kind_phys, 2.651789e-02_kind_phys, 2.490518e-02_kind_phys,  & !12
       2.351273e-02_kind_phys, 2.229056e-02_kind_phys, 2.120335e-02_kind_phys, 2.022541e-02_kind_phys, 1.933763e-02_kind_phys,  & !12
       1.852546e-02_kind_phys, 1.777763e-02_kind_phys, 1.708528e-02_kind_phys, 1.644134e-02_kind_phys, 1.584009e-02_kind_phys,  & !12
       1.527684e-02_kind_phys, 1.474774e-02_kind_phys, 1.424955e-02_kind_phys, 1.377957e-02_kind_phys, 1.333549e-02_kind_phys,  & !12
       1.291534e-02_kind_phys, 1.251743e-02_kind_phys, 1.214029e-02_kind_phys, 1.178265e-02_kind_phys, 1.144337e-02_kind_phys,  & !12
       1.112148e-02_kind_phys, 1.081609e-02_kind_phys, 1.052642e-02_kind_phys, 1.025178e-02_kind_phys, 9.991540e-03_kind_phys,  & !12
       9.745130e-03_kind_phys, 9.512038e-03_kind_phys, 9.291797e-03_kind_phys, 9.083980e-03_kind_phys, 8.888195e-03_kind_phys,  & !12
       8.704081e-03_kind_phys, 8.531306e-03_kind_phys, 8.369560e-03_kind_phys, 8.218558e-03_kind_phys, 8.078032e-03_kind_phys,  & !12
       7.947730e-03_kind_phys,                                                                                                  & !12
       1.560311e-01_kind_phys, 9.961097e-02_kind_phys, 7.502949e-02_kind_phys, 6.115022e-02_kind_phys, 5.214952e-02_kind_phys,  & !13
       4.578149e-02_kind_phys, 4.099731e-02_kind_phys, 3.724174e-02_kind_phys, 3.419343e-02_kind_phys, 3.165356e-02_kind_phys,  & !13
       2.949251e-02_kind_phys, 2.762222e-02_kind_phys, 2.598073e-02_kind_phys, 2.452322e-02_kind_phys, 2.321642e-02_kind_phys,  & !13
       2.203516e-02_kind_phys, 2.096002e-02_kind_phys, 1.997579e-02_kind_phys, 1.907036e-02_kind_phys, 1.823401e-02_kind_phys,  & !13
       1.745879e-02_kind_phys, 1.673819e-02_kind_phys, 1.606678e-02_kind_phys, 1.544003e-02_kind_phys, 1.485411e-02_kind_phys,  & !13
       1.430574e-02_kind_phys, 1.379215e-02_kind_phys, 1.331092e-02_kind_phys, 1.285996e-02_kind_phys, 1.243746e-02_kind_phys,  & !13
       1.204183e-02_kind_phys, 1.167164e-02_kind_phys, 1.132567e-02_kind_phys, 1.100281e-02_kind_phys, 1.070207e-02_kind_phys,  & !13
       1.042258e-02_kind_phys, 1.016352e-02_kind_phys, 9.924197e-03_kind_phys, 9.703953e-03_kind_phys, 9.502199e-03_kind_phys,  & !13
       9.318400e-03_kind_phys, 9.152066e-03_kind_phys, 9.002749e-03_kind_phys, 8.870038e-03_kind_phys, 8.753555e-03_kind_phys,  & !13
       8.652951e-03_kind_phys,                                                                                                  & !13
       1.559547e-01_kind_phys, 9.896700e-02_kind_phys, 7.441231e-02_kind_phys, 6.061469e-02_kind_phys, 5.168730e-02_kind_phys,  & !14
       4.537821e-02_kind_phys, 4.064106e-02_kind_phys, 3.692367e-02_kind_phys, 3.390714e-02_kind_phys, 3.139438e-02_kind_phys,  & !14
       2.925702e-02_kind_phys, 2.740783e-02_kind_phys, 2.578547e-02_kind_phys, 2.434552e-02_kind_phys, 2.305506e-02_kind_phys,  & !14
       2.188910e-02_kind_phys, 2.082842e-02_kind_phys, 1.985789e-02_kind_phys, 1.896553e-02_kind_phys, 1.814165e-02_kind_phys,  & !14
       1.737839e-02_kind_phys, 1.666927e-02_kind_phys, 1.600891e-02_kind_phys, 1.539279e-02_kind_phys, 1.481712e-02_kind_phys,  & !14
       1.427865e-02_kind_phys, 1.377463e-02_kind_phys, 1.330266e-02_kind_phys, 1.286068e-02_kind_phys, 1.244689e-02_kind_phys,  & !14
       1.205973e-02_kind_phys, 1.169780e-02_kind_phys, 1.135989e-02_kind_phys, 1.104492e-02_kind_phys, 1.075192e-02_kind_phys,  & !14
       1.048004e-02_kind_phys, 1.022850e-02_kind_phys, 9.996611e-03_kind_phys, 9.783753e-03_kind_phys, 9.589361e-03_kind_phys,  & !14
       9.412924e-03_kind_phys, 9.253977e-03_kind_phys, 9.112098e-03_kind_phys, 8.986903e-03_kind_phys, 8.878039e-03_kind_phys,  & !14
       8.785184e-03_kind_phys,                                                                                                  & !14
       1.102926e-01_kind_phys, 7.176622e-02_kind_phys, 5.530316e-02_kind_phys, 4.606056e-02_kind_phys, 4.006116e-02_kind_phys,  & !15
       3.579628e-02_kind_phys, 3.256909e-02_kind_phys, 3.001360e-02_kind_phys, 2.791920e-02_kind_phys, 2.615617e-02_kind_phys,  & !15
       2.464023e-02_kind_phys, 2.331426e-02_kind_phys, 2.213817e-02_kind_phys, 2.108301e-02_kind_phys, 2.012733e-02_kind_phys,  & !15
       1.925493e-02_kind_phys, 1.845331e-02_kind_phys, 1.771269e-02_kind_phys, 1.702531e-02_kind_phys, 1.638493e-02_kind_phys,  & !15
       1.578648e-02_kind_phys, 1.522579e-02_kind_phys, 1.469940e-02_kind_phys, 1.420442e-02_kind_phys, 1.373841e-02_kind_phys,  & !15
       1.329931e-02_kind_phys, 1.288535e-02_kind_phys, 1.249502e-02_kind_phys, 1.212700e-02_kind_phys, 1.178015e-02_kind_phys,  & !15
       1.145348e-02_kind_phys, 1.114612e-02_kind_phys, 1.085730e-02_kind_phys, 1.058633e-02_kind_phys, 1.033263e-02_kind_phys,  & !15
       1.009564e-02_kind_phys, 9.874895e-03_kind_phys, 9.669960e-03_kind_phys, 9.480449e-03_kind_phys, 9.306014e-03_kind_phys,  & !15
       9.146339e-03_kind_phys, 9.001138e-03_kind_phys, 8.870154e-03_kind_phys, 8.753148e-03_kind_phys, 8.649907e-03_kind_phys,  & !15
       8.560232e-03_kind_phys,                                                                                                  & !15
       1.688344e-01_kind_phys, 1.077072e-01_kind_phys, 7.994467e-02_kind_phys, 6.403862e-02_kind_phys, 5.369850e-02_kind_phys,  & !16
       4.641582e-02_kind_phys, 4.099331e-02_kind_phys, 3.678724e-02_kind_phys, 3.342069e-02_kind_phys, 3.065831e-02_kind_phys,  & !16
       2.834557e-02_kind_phys, 2.637680e-02_kind_phys, 2.467733e-02_kind_phys, 2.319286e-02_kind_phys, 2.188299e-02_kind_phys,  & !16
       2.071701e-02_kind_phys, 1.967121e-02_kind_phys, 1.872692e-02_kind_phys, 1.786931e-02_kind_phys, 1.708641e-02_kind_phys,  & !16
       1.636846e-02_kind_phys, 1.570743e-02_kind_phys, 1.509665e-02_kind_phys, 1.453052e-02_kind_phys, 1.400433e-02_kind_phys,  & !16
       1.351407e-02_kind_phys, 1.305631e-02_kind_phys, 1.262810e-02_kind_phys, 1.222688e-02_kind_phys, 1.185044e-02_kind_phys,  & !16
       1.149683e-02_kind_phys, 1.116436e-02_kind_phys, 1.085153e-02_kind_phys, 1.055701e-02_kind_phys, 1.027961e-02_kind_phys,  & !16
       1.001831e-02_kind_phys, 9.772141e-03_kind_phys, 9.540280e-03_kind_phys, 9.321966e-03_kind_phys, 9.116517e-03_kind_phys,  & !16
       8.923315e-03_kind_phys, 8.741803e-03_kind_phys, 8.571472e-03_kind_phys, 8.411860e-03_kind_phys, 8.262543e-03_kind_phys,  & !16
       8.123136e-03_kind_phys/),                                                                                                & !16
       shape=(/46,nBandsLW_RRTMG/))      
contains
  ! #######################################################################################
  ! subroutine rrtmgp_lw_cloud_optics
  ! #######################################################################################
  subroutine rrtmgp_lw_cloud_optics(ncol, nlay, nBandsLW, cld_lwp, cld_ref_liq, cld_iwp,  &
       cld_ref_ice, cld_rwp, cld_ref_rain, cld_swp, cld_ref_snow, cld_frac, tau_cld)
    ! Inputs
    integer,intent(in) :: &
         nBandsLW,     & ! Number of spectral bands
         ncol,         & ! Number of horizontal gridpoints
         nlay            ! Number of vertical layers
    real(kind_phys), dimension(ncol,nlay), intent(in) :: &
         cld_frac,     & ! Cloud-fraction                         (1)
         cld_lwp,      & ! Cloud liquid water path                (g/m2)
         cld_ref_liq,  & ! Effective radius (liquid)              (micron)
         cld_iwp,      & ! Cloud ice water path                   (g/m2)
         cld_ref_ice,  & ! Effective radius (ice)                 (micron)
         cld_rwp,      & ! Cloud rain water path                  (g/m2)
         cld_ref_rain, & ! Effective radius (rain-drop)           (micron)
         cld_swp,      & ! Cloud snow-water path                  (g/m2)
         cld_ref_snow    ! Effective radius (snow-flake)          (micron)

     ! Outputs
     real(kind_phys),dimension(nBandsLW,ncol,nlay),intent(out) :: &
          tau_cld       

     ! Local variables
     integer :: ij,ik,ib,index,ia
     real(kind_phys) :: factor,fint,cld_ref_iceTemp
     real(kind_phys),dimension(ncol,nlay) :: tau_snow, tau_rain
     real(kind_phys),dimension(nBandsLW,ncol,nlay) :: tau_liq, tau_ice

     if (ilwcliq .gt. 0) then 
        do ij=1,ncol
           do ik=1,nlay
              if (cld_frac(ij,ik) .gt. 0._kind_phys) then
                 ! Rain optical-depth (No band dependence)
                 tau_rain(ij,ik) = absrain*cld_rwp(ij,ik)

                 ! Snow optical-depth (No band dependence)
                 if (cld_swp(ij,ik) .gt. 0._kind_phys .and. cld_ref_snow(ij,ik) .gt. 10._kind_phys) then
                    tau_snow(ij,ik) = abssnow0*1.05756*cld_swp(ij,ik)/cld_ref_snow(ij,ik)
                 else
                    tau_snow(ij,ik) = 0._kind_phys
                 endif

                 ! Liquid water opitcal-depth
                 if (cld_lwp(ij,ik) .le. 0._kind_phys) then
                    tau_liq(:,ij,ik) = 0._kind_phys
                 else
                    if (ilwcliq .eq. 1) then
                       factor = cld_ref_liq(ij,ik) - 1.5_kind_phys
                       index  = max( 1, min( 57, int( factor ) ))
                       fint   = factor - float(index)
                       do ib=1,nBandsLW
                          tau_liq(ib,ij,ik) = max(0._kind_phys, cld_lwp(ij,ik)*(absliq1(index,ib) + &
                               fint*(absliq1(index+1,ib)-absliq1(index,ib)) ))
                       enddo
                    endif
                 endif

                 ! Ice water optical-depth
                 if (cld_iwp(ij,ik) .le. 0._kind_phys) then
                    tau_ice(:,ij,ik) = 0._kind_phys
                 else
                    ! 1) Ebert and curry approach for all particle sizes. (bound between 13-130microns)
                    if (ilwcice .eq. 1) then
                       cld_ref_iceTemp = min(130._kind_phys, max(13._kind_phys,real(cld_ref_ice(ij,ik))))
                       do ib=1,nBandsLW
                          ia = ipat(ib)             ! eb_&_c band index for ice cloud coeff
                          tau_ice(ib,ij,ik) = max(0._kind_phys, cld_iwp(ij,ik)*(absice1(1,ia) + absice1(2,ia)/cld_ref_iceTemp) )
                       enddo

                    ! 2) Streamer approach for ice effective radius between 5.0 and 131.0 microns
                    !    and ebert and curry approach for ice eff radius greater than 131.0 microns.
                    !    no smoothing between the transition of the two methods
                    elseif (ilwcice .eq. 2) then
                       factor = (cld_ref_ice(ij,ik) - 2._kind_phys) / 3._kind_phys
                       index  = max( 1, min( 42, int( factor ) ))
                       fint   = factor - float(index)
                       do ib = 1, nBandsLW
                          tau_ice(ib,ij,ik) = max(0._kind_phys, cld_iwp(ij,ik)*(absice2(index,ib) + &
                               fint*(absice2(index+1,ib) - absice2(index,ib)) ))
                       enddo
                    ! 3) Fu's approach for ice effective radius between 4.8 and 135 microns
                    !    (generalized effective size from 5 to 140 microns)
                    elseif (ilwcice .eq. 3) then
                       cld_ref_iceTemp = max(5._kind_phys, 1.0315_kind_phys*cld_ref_ice(ij,ik))              ! v4.71 value
                       factor = (cld_ref_iceTemp - 2._kind_phys) / 3._kind_phys
                       index  = max( 1, min( 45, int( factor ) ))
                       fint   = factor - float(index)
                       do ib = 1, nBandsLW
                          tau_ice(ib,ij,ik) = max(0._kind_phys, cld_iwp(ij,ik)*(absice3(index,ib) + &
                               fint*(absice3(index+1,ib) - absice3(index,ib)) ))
                       enddo
                    endif
                 endif 
              else
                 tau_rain(ij,ik)  = 0._kind_phys
                 tau_snow(ij,ik)  = 0._kind_phys
                 tau_liq(:,ij,ik) = 0._kind_phys
                 tau_ice(:,ij,ik) = 0._kind_phys
              endif
              ! Cloud optical depth
              do ib = 1, nBandsLW
                 tau_cld(ib,ij,ik) = tau_ice(ib,ij,ik) + tau_liq(ib,ij,ik) + tau_rain(ij,ik) + tau_snow(ij,ik)
              enddo
           end do
        end do
     endif
   end subroutine rrtmgp_lw_cloud_optics
   ! #######################################################################################
   ! SUBROUTINE mcica_subcol_lw
   ! #######################################################################################
   subroutine mcica_subcol_lw(ncol, nlay, ngpts, cld_frac, icseed,  dzlyr, de_lgth, cld_frac_mcica)
     ! Inputs
    integer,intent(in) :: &
         ncol,         & ! Number of horizontal gridpoints
         nlay,         & ! Number of vertical layers
         ngpts           ! Number of spectral g-points
    integer,dimension(ncol),intent(in) :: &
         icseed          ! Permutation seed for each column.
    real(kind_phys), dimension(ncol), intent(in) :: &
         de_lgth         ! Cloud decorrelation length (km)
    real(kind_phys), dimension(ncol,nlay), intent(in) :: &
         cld_frac,     & ! Cloud-fraction  
         dzlyr           ! Layer thinkness (km)
    ! Outputs
    real(kind_phys),dimension(ngpts,ncol,nlay),intent(out) :: &
         cld_frac_mcica
    ! Local variables
    type(random_stat) :: stat
    integer :: icol,n,k,k1
    real(kind_phys) :: tem1
    real(kind_phys),dimension(ngpts) :: rand1D
    real(kind_phys),dimension(nlay*ngpts) :: rand2D
    real(kind_phys),dimension(ngpts,nlay) :: cdfunc,cdfun2
    real(kind_phys),dimension(nlay) :: fac_lcf
    logical,dimension(ngpts,nlay) :: lcloudy

    ! Loop over all columns
    do icol=1,ncol
       ! Call random_setseed() to advance random number generator by "icseed" values.
       call random_setseed(icseed(icol),stat)

       ! ###################################################################################
       ! Sub-column set up according to overlapping assumption:
       !  - For random overlap, pick a random value at every level 
       !  - For max-random overlap, pick a random value at every level
       !  - For maximum overlap, pick same random numebr at every level
       ! ###################################################################################
       select case ( iovrlw )
       ! ###################################################################################
       ! 0) Random overlap
       ! ###################################################################################
       case( 0 )
          call random_number(rand2D,stat)
          k1 = 0
          do n = 1, ngpts
             do k = 1, nlay
                k1 = k1 + 1
                cdfunc(n,k) = rand2d(k1)
             enddo
          enddo

       ! ###################################################################################
       ! 1) Maximum-random overlap
       ! ###################################################################################
       case(1)
          call random_number(rand2D,stat)
          k1 = 0
          do n = 1, ngpts
             do k = 1, nlay
                k1 = k1 + 1
                cdfunc(n,k) = rand2d(k1)
             enddo
          enddo          

          ! First pick a random number for bottom (or top) layer.
          ! then walk up the column: (aer's code)
          ! if layer below is cloudy, use the same rand num in the layer below
          ! if layer below is clear,  use a new random number
          do k = 2, nlay
             k1 = k - 1
             tem1 = 1._kind_phys - cld_frac(icol,k1)
             do n = 1, ngpts
                if ( cdfunc(n,k1) > tem1 ) then
                   cdfunc(n,k) = cdfunc(n,k1)
                else
                   cdfunc(n,k) = cdfunc(n,k) * tem1
                endif
             enddo
          enddo

       ! ###################################################################################
       ! 2) Maximum overlap
       ! ###################################################################################
       case(2)
          call random_number(rand1d,stat)
          do n = 1, ngpts
             tem1 = rand1d(n)
             do k = 1, nlay
                cdfunc(n,k) = tem1
             enddo
          enddo

       ! ###################################################################################
       ! 3) Decorrelation length
       ! ###################################################################################
       case(3)
          ! Compute overlapping factors based on layer midpoint distances and decorrelation 
          ! depths
          do k = nlay, 2, -1
             fac_lcf(k) = exp( -0.5 * (dzlyr(iCol,k)+dzlyr(iCol,k-1)) / de_lgth(iCol) )
          enddo

          ! Setup 2 sets of random numbers
          call random_number ( rand2d, stat )
          k1 = 0
          do k = 1, nlay
             do n = 1, ngpts
                k1 = k1 + 1
                cdfunc(n,k) = rand2d(k1)
             enddo
          enddo
          !
          call random_number ( rand2d, stat )
          k1 = 0
          do k = 1, nlay
             do n = 1, ngpts
                k1 = k1 + 1
                cdfun2(n,k) = rand2d(k1)
             enddo
          enddo

          ! Then working from the top down:
          !   if a random number (from an independent set -cdfun2) is smaller then the
          !   scale factor: use the upper layer's number,  otherwise use a new random
          !   number (keep the original assigned one).
          do k = nlay-1, 1, -1
             k1 = k + 1
             do n = 1, ngpts
                if ( cdfun2(n,k) <= fac_lcf(k1) ) then
                   cdfunc(n,k) = cdfunc(n,k1)
                endif
             enddo
          enddo

       end select
       
       ! ###################################################################################
       ! Generate subcolumn cloud mask (0/1 for clear/cloudy)
       ! ###################################################################################
       do k = 1, nlay
          tem1 = 1._kind_phys - cld_frac(icol,k)
          do n = 1, ngpts
             lcloudy(n,k) = cdfunc(n,k) >= tem1
             if (lcloudy(n,k)) then
                cld_frac_mcica(n,icol,k) = 1._kind_phys
             else
                cld_frac_mcica(n,icol,k) = 0._kind_phys
             endif
          enddo
       enddo
    enddo ! END LOOP OVER COLUMNS
  end subroutine mcica_subcol_lw

end module mo_rrtmgp_lw_cloud_optics
