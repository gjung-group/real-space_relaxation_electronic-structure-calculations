program GenMoire

   implicit none

   !double precision :: phideg = 21.051724d0 ! Rotation angle in degrees
   !double precision :: phideg = 6.01d0 ! Rotation angle in degrees
   double precision :: phideg = 1.12d0 ! Rotation angle in degrees
   double precision, parameter :: degtol=0.01d0, lattol = 0.01d0 ! Tolerance
   double precision, parameter :: h = 3.3d0, dz = 35.0d0
   !double precision, parameter :: aG = 2.46d0, aBN0 = 6.5796819680864944d0!2.49d0 ! Lattice parameters
   double precision, parameter :: aG = 2.46019d0, aBN0 = 2.46019d0 ! Lattice parameters
   character(*), parameter :: select = 'delta'
   logical, parameter :: fdfbuild = .true.

   double precision, parameter :: basis(2,2) = reshape((/1.0d0/3.0d0,1.0d0/3.0d0, 2.0d0/3.0d0,2.0d0/3.0d0/),(/2,2/))
   double precision :: Rot(3,3)
   double precision, parameter :: pi = acos(-1.0d0)
   character(*), parameter :: grSp(2) = (/'C','D'/), bnSp(2) =(/'E','F'/)
   integer, parameter :: grN(2) = (/6,6/), bnN(2) = (/5,7/)

   double precision :: delta, deltol, phi, g, lambda
   integer :: n1, n2, nMoire1, nMoire2, n(4)
   double precision, pointer :: Xg(:,:), Xbn(:,:)
   double precision, pointer :: Xg1(:,:), Xbn1(:,:)
   integer :: ncell, i

   double precision :: yShift
   double precision :: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, a1(3), a1Rot(3), a2(3), a2Rot(3), a3(3), a3Rot(3)

   integer :: firstN


   n1 = 0
   n2 = 100
   nMoire1 = 1
   nMoire2 = 100
   delta = aBN0/aG - 1.0d0
   deltol = (aBN0 + lattol)/aG - 1.0d0 - delta
   call MoireFind(n1,n2,nMoire1,nMoire2,phideg,delta,degtol,deltol,select,n)
   !firstN = XXX
   !n(1) = firstN
   !n(2) = firstN - 1
   !n(3) = firstN - 1
   !n(4) = firstN
   if (fdfbuild) then
      g = n(1)**2 + n(2)**2 + n(1)*n(2)
      delta = sqrt(real(n(3)**2 + n(4)**2 + n(3)*n(4))/g)
      phi = acos((2.0d0*n(1)*n(3)+2.0d0*n(2)*n(4) + n(1)*n(4) + n(2)*n(3))/(2.0d0*delta*g))
      if (delta<1.0d0) then
         delta = 1.0d0/delta
      else
         i = n(1)
         n(1) = n(3)
         n(3) = i
         i = n(2)
         n(2) = n(4)
         n(4) = i
      end if
      write(*,*)
      write(*,'(a,f13.6)') 'Final angle:', phi*180.0d0/pi
      write(*,'(a,f13.6)') 'Final BN lattice parameter:', aG*delta
      ncell = n(1)**2 + n(1)*n(2) + n(2)**2
      allocate(Xg(2,ncell*2))
      allocate(Xg1(2,ncell*2))
      write(*,'(a,i8)') 'Atoms for graphene:', ncell*2
      ncell = n(3)**2 + n(3)*n(4) + n(4)**2
      allocate(Xbn(2,ncell*2))
      allocate(Xbn1(2,ncell*2))
      write(*,'(a,i8)') 'Atoms for BN:', ncell*2
      write(*,*)
      call MoireConstruct(Xg,n(1:2),basis,0.0d0,1.0d0)
      call MoireConstruct(Xbn,n(3:4),basis,-phi*180.0d0/pi,delta)
      call MoireWrite(Xg,Xbn,n,aG,-h,dz,'moire.fdf',grSp,bnSp,grN,bnN)
      !yShift = (aG/3.0)/(n(1)*aG)
      !print*, yShift
      !print*, basis 
      !call MoireConstruct(Xg,n(1:2),basis,0.0d0,1.0d0)
      !call MoireConstruct(Xbn,n(3:4),basis,-phi*180.0d0/pi,delta)
      !call MoireConstructShift(Xg1,n(1:2),basis,0.0d0,1.0d0,yShift)
      !call MoireConstructShift(Xbn1,n(3:4),basis,-phi*180.0d0/pi,delta,yShift)
      !call MoireWrite4L(Xg,Xbn,Xg1,Xbn1,n,aG,-h,dz,'moire.fdf',grSp,bnSp,grN,bnN)
      delta = delta - 1.0d0
      lambda = (1.0d0+delta)*aG/sqrt(2.0d0*(1.0d0+delta)*(1.0d0-cos(phi))+delta**2)
      write(*,*) 'lambda1:', lambda
      lambda = aG*sqrt(1.0d0*n(1)**2+n(2)**2+n(1)*n(2))
      !lambda = aG*sqrt(1.0d0*n(1)**2+n(2)**2+n(1)*n(2))
      write(*,*) 'lambda2:', lambda
      !lambda = 3.0d0*aG**2/2.0d0
      !write(*,*) 'lambda3:', lambda
   end if

contains

subroutine MoireFind(n1,n2,nM1,nM2,phiin,delta,degtol,lattol,select,mout)

   implicit none

   double precision :: pi = acos(-1.0d0)

   integer, intent(in) :: n1, n2, nM1, nM2
   double precision, intent(in) :: phiin, delta, degtol, lattol
   character(*), intent(in) :: select
   integer, intent(out) :: mout(4)

   integer :: a, b, ap, bp, m(2,2), mp(2,2), mm(2,2)
   integer :: nn, N, nn2
   double precision :: p, alpha, g, angletol, pcomp, d, phi
   double precision :: ep, ea,e, ep0, ea0, e0
   double precision :: m0p(4), m0a(4), m0(4)

   angletol = degtol*pi/180.0d0
   phi = phiin*pi/180.0d0
   ep0 = huge(1.0d0)
   ea0 = huge(1.0d0)
   e0 = huge(1.0d0)
   write(*,'(a,f16.9)') 'Angle:', phiin
   write(*,'(a,f16.9)') 'Delta:', delta
   do a=n1,n2; do b=n1,a; do ap=n1,a; do bp=n1,a
      m(:,1) = (/a,b/)
      m(:,2) = (/-b,a+b/)
      mp(:,1) = (/ap,bp/)
      mp(:,2) = (/-bp,ap+bp/)
      mm = m-mp
      !N = a**2 + b**2 + a*b      
      !if (mod(N,3)==0) write(*,*) N, a, b, ap, bp
      N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
      if (N>=nM1 .and. N<=nM2) then
         g = a**2 + b**2 + a*b
         p = sqrt(real(ap**2+bp**2+ap*bp)/g)
         nn = 2*a*ap + 2*b*bp + b*ap + a*bp
         alpha = acos(real(nn,8)/(2.0d0*p*g))
         if (p < 1.0d0) then
            pcomp = 1.0d0/p
         else
            pcomp = p
         end if
         d = pcomp - 1.0d0
         ep = abs(d-delta)
         ea = abs(alpha-phi)
         if ((ep < lattol) .and. (ea < angletol)) then
         !if (abs(alpha-phi) < angletol) then
            if (ea < ea0) then
               ea0 = ea
               m0a = (/a,b,ap,bp/)
            end if
            if (ep < ep0) then
               ep0 = ep
               m0p = (/a,b,ap,bp/)
            end if
            !e = sqrt(ep**2/(delta**2+epsilon(1.0d0)) + ea**2/(phi**2+epsilon(1.0d0)))
            e = sqrt(ep**2/delta**2 + ea**2/phi**2)
            if (e < e0) then
               e0 = e
               m0 = (/a,b,ap,bp/)
            end if
            !write(*,'(4i6,2f16.9,i6)') a,b,ap,bp,pcomp, alpha*180.0d0/pi, N
         end if
      end if
   end do; end do; end do; end do
   write(*,*)
   write(*,'(a)') "                a     b     a'    b'       Delta            Angle       N   #1   #2"
   a = m0a(1)
   b = m0a(2)
   ap = m0a(3)
   bp = m0a(4)
   m(:,1) = (/a,b/)
   m(:,2) = (/-b,a+b/)
   mp(:,1) = (/ap,bp/)
   mp(:,2) = (/-bp,ap+bp/)
   mm = m-mp
   N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
   g = a**2 + b**2 + a*b
   p = sqrt(real(ap**2+bp**2+ap*bp)/g)
   nn = 2*a*ap + 2*b*bp + b*ap + a*bp
   alpha = acos(real(nn,8)/(2.0d0*p*g))
   if (p < 1.0d0) then
      p = 1.0d0/p
      a = m0a(3)
      b = m0a(4)
      ap = m0a(1)
      bp = m0a(2)
   end if
   d = p - 1.0d0
   nn = a**2 + a*b + b**2
   nn2 = ap**2 + ap*bp + bp**2
   write(*,'(a,4i6,2f16.9,i5,2i6)') 'For angle:  ',a,b,ap,bp,d,alpha*180.0d0/pi,N,nn, nn2
   a = m0p(1)
   b = m0p(2)
   ap = m0p(3)
   bp = m0p(4)
   m(:,1) = (/a,b/)
   m(:,2) = (/-b,a+b/)
   mp(:,1) = (/ap,bp/)
   mp(:,2) = (/-bp,ap+bp/)
   mm = m-mp
   N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
   g = a**2 + b**2 + a*b
   p = sqrt(real(ap**2+bp**2+ap*bp)/g)
   nn = 2*a*ap + 2*b*bp + b*ap + a*bp
   alpha = acos(real(nn,8)/(2.0d0*p*g))
   if (p < 1.0d0) then
      p = 1.0d0/p
      a = m0p(3)
      b = m0p(4)
      ap = m0p(1)
      bp = m0p(2)
   end if
   d = p - 1.0d0
   nn = a**2 + a*b + b**2
   nn2 = ap**2 + ap*bp + bp**2
   write(*,'(a,4i6,2f16.9,i5,2i6)') 'For delta:  ',a,b,ap,bp,d,alpha*180.0d0/pi,N, nn, nn2
   a = m0(1)
   b = m0(2)
   ap = m0(3)
   bp = m0(4)
   m(:,1) = (/a,b/)
   m(:,2) = (/-b,a+b/)
   mp(:,1) = (/ap,bp/)
   mp(:,2) = (/-bp,ap+bp/)
   mm = m-mp
   N = mm(1,1)*mm(2,2) - mm(1,2)*mm(2,1)
   g = a**2 + b**2 + a*b
   p = sqrt(real(ap**2+bp**2+ap*bp)/g)
   nn = 2*a*ap + 2*b*bp + b*ap + a*bp
   alpha = acos(real(nn,8)/(2.0d0*p*g))
   if (p < 1.0d0) then
      p = 1.0d0/p
      a = m0(3)
      b = m0(4)
      ap = m0(1)
      bp = m0(2)
   end if
   d = p - 1.0d0
   nn = a**2 + a*b + b**2
   nn2 = ap**2 + ap*bp + bp**2
   write(*,'(a,4i6,2f16.9,i5,2i6)') 'For both:   ',a,b,ap,bp,d,alpha*180.0d0/pi,N,nn,nn2
   if (select=='angle' .or. select=='ANGLE') then
      mout = m0a
   else if (select=='delta' .or. select=='DELTA') then
      mout = m0p
   else if (select=='both' .or. select=='BOTH') then
      mout = m0
   end if
   if (all(mout==0)) stop 'All integers are 0'

end subroutine MoireFind

subroutine MoireConstruct(X,mm,bss,angle,p)

   implicit none

   double precision, parameter :: pi = acos(-1.0d0)

   double precision, intent(out) :: X(:,:)
   integer, intent(in) :: mm(2)
   double precision, intent(in) :: bss(:,:), angle, p

   integer :: ncell, nbasis, i1, i2, j, k, n1, n2, mint(2,2)
   double precision :: ucell(2,2), Rcell(2,2), a, m(2,2)
   double precision :: f(2), sq
   double precision, pointer :: basis(:,:)
   integer, save :: nlayer=0

   ncell = mm(1)**2 + mm(1)*mm(2) + mm(2)**2
   nbasis = size(bss,2)
   allocate(basis(2,nbasis))
   if (size(X) /= 2*ncell*nbasis) then
      stop 'Wrong size of X array'
   end if
   a = angle*pi/180.0d0
   ucell(:,1) = (/cos(a),sin(a)/)
   ucell(:,2) = (/cos(a+pi/3.0d0),sin(a+pi/3)/)
   ucell = ucell*p
   m(:,1) = (/mm(1),mm(2)/)
   m(:,2) = (/-mm(2),mm(1)+mm(2)/)
   Rcell = matmul(ucell,m)
   basis = matmul(ucell,bss)
   sq = sqrt(real(ncell))
   mint = m
   n1 = gcd(mint(1,1),mint(2,1))
   n2 = ncell/n1
   nlayer = nlayer  + 1
   write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'Size of layer ',nlayer,':   ',n1,' x ',n2, '   (',mm(1),' x ',mm(2),')'
   !write(*,*) ucell
   k = 0
   !write(*,*) sqrt(Rcell(1,1)**2+Rcell(2,1)**2)/sqrt(ucell(1,1)**2+ucell(2,1)**2)
   do i1 =1,n1; do i2=1,n2
      do j=1,nbasis
         k = k + 1
         X(:,k) = (i1-1)*ucell(:,1) + (i2-1)*ucell(:,2) + basis(:,j)
      end do
   end do; end do
   do i1=1,ncell*nbasis
      !print*, Rcell(1,1), X(2,i1), X(1,i1)
      f(2) = Rcell(1,1)*(X(2,i1)-X(1,i1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      f(1) = (X(1,i1)-f(2)*Rcell(1,2))/Rcell(1,1)
      X(:,i1) = mod(f,1.0d0)
      if (X(1,i1)<0.0d0) X(1,i1) = X(1,i1) + 1.0d0
      if (X(2,i1)<0.0d0) X(2,i1) = X(2,i1) + 1.0d0
   end do
   deallocate(basis)

end subroutine MoireConstruct

subroutine MoireConstructShift(X,mm,bss,angle,p,yShift)

   implicit none

   double precision, parameter :: pi = acos(-1.0d0)

   double precision, intent(out) :: X(:,:)
   integer, intent(in) :: mm(2)
   double precision, intent(in) :: bss(:,:), angle, p

   integer :: ncell, nbasis, i1, i2, j, k, n1, n2, mint(2,2)
   double precision :: ucell(2,2), Rcell(2,2), a, m(2,2)
   double precision :: f(2), sq
   double precision, pointer :: basis(:,:)
   integer, save :: nlayer=0

   double precision :: yShift

   ncell = mm(1)**2 + mm(1)*mm(2) + mm(2)**2
   nbasis = size(bss,2)
   allocate(basis(2,nbasis))
   if (size(X) /= 2*ncell*nbasis) then
      stop 'Wrong size of X array'
   end if
   a = angle*pi/180.0d0
   ucell(:,1) = (/cos(a),sin(a)/)
   ucell(:,2) = (/cos(a+pi/3.0d0),sin(a+pi/3)/)
   ucell = ucell*p
   m(:,1) = (/mm(1),mm(2)/)
   m(:,2) = (/-mm(2),mm(1)+mm(2)/)
   Rcell = matmul(ucell,m)
   basis = matmul(ucell,bss)
   sq = sqrt(real(ncell))
   mint = m
   n1 = gcd(mint(1,1),mint(2,1))
   n2 = ncell/n1
   nlayer = nlayer  + 1
   write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'Size of layer ',nlayer,':   ',n1,' x ',n2, '   (',mm(1),' x ',mm(2),')'
   !write(*,*) ucell
   k = 0
   !write(*,*) sqrt(Rcell(1,1)**2+Rcell(2,1)**2)/sqrt(ucell(1,1)**2+ucell(2,1)**2)
   do i1 =1,n1; do i2=1,n2
      do j=1,nbasis
         k = k + 1
         X(:,k) = (i1-1)*ucell(:,1) + (i2-1)*ucell(:,2) + basis(:,j)
      end do
   end do; end do
   print*, yShift
   do i1=1,ncell*nbasis
      !print*, Rcell(1,1), X(2,i1), X(1,i1)
      f(2) = Rcell(1,1)*(X(2,i1)-X(1,i1)*Rcell(2,1)/Rcell(1,1))/(Rcell(2,2)*Rcell(1,1)-Rcell(1,2)*Rcell(2,1))
      f(1) = (X(1,i1)-f(2)*Rcell(1,2))/Rcell(1,1)
      !X(:,i1) = mod(f,1.0d0)
      !if (X(1,i1)<0.0d0) X(1,i1) = X(1,i1) + 1.0d0
      !if (X(2,i1)<0.0d0) X(2,i1) = X(2,i1) + 1.0d0

      X(:,i1) = [mod(f(1),1.0), mod(f(2),1.0)+yShift]!, 0.5_dp+d]
      if (X(1,i1)<0.0) X(1,i1) = mod(X(1,i1) + 1.0,1.0)
      if (X(2,i1)<0.0) X(2,i1) = mod(X(2,i1) + 1.0,1.0)
      if (X(1,i1)>1.0) X(1,i1) = mod(X(1,i1) - 1.0,1.0)
      if (X(2,i1)>1.0) X(2,i1) = mod(X(2,i1) - 1.0,1.0)

   end do
   deallocate(basis)

end subroutine MoireConstructShift

subroutine MoireWrite(Xs,Xo,mm,a,h,z,name,sp1,sp2,N1,N2)

   implicit none

   integer, parameter :: u=10, u2=11, u3=12, u4=13, u5=14, u6=15

   double precision, intent(in) :: Xs(:,:), Xo(:,:), a, h, z
   integer, intent(in) :: mm(4), N1(:), N2(:)
   character(*), intent(in) :: name, sp1(:), sp2(:)

   double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2)
   double precision :: rz
   integer :: sz1, sz2, i, j, k, nb1, nb2
   character(60), pointer :: Sp(:)
   integer, pointer :: Nel(:), num1(:), num2(:)
   double precision :: phirad, e1(3)

   open(u,FILE=name,STATUS='replace')
   open(u2,FILE="BLBL.frac",STATUS='replace')
   open(u3,FILE="BLBL.basis1",STATUS='replace')
   open(u4,FILE="BLBL.basis2",STATUS='replace')
   open(u5,FILE="BLBL.mol", STATUS='replace')
   open(u6,FILE="sublattices.dat", STATUS='replace')
   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)
   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
   Rcell = 0.0d0
   Rcell(1:2,1:2) = matmul(ucell,m1)
   Rcell(3,3) = z
   sz1 = size(Xs,2)
   sz2 = size(Xo,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
   allocate(Sp(nb1+nb2))
   allocate(Nel(nb1+nb2))
   allocate(num1(nb1))
   allocate(num2(nb2))
   Sp(1) = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = 1
out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
   end do out2
   write(u,'(a,i0)') 'NumberOfSpecies', k
   write(u,'(a)') '%block ChemicalSpeciesLabel'
   do i=1,k
      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
   end do
   !a1 = (/Rcell(1,1), Rcell(1,2), Rcell(1,3)/)
   !a2 = (/Rcell(2,1), Rcell(2,2), Rcell(2,3)/)
   !a3 = (/Rcell(3,1), Rcell(3,2), Rcell(3,3)/)
   a1 = Rcell(:,1)
   a2 = Rcell(:,2)
   a3 = Rcell(:,3)
   e1 = (/1,0,0/)
   phirad = -acos(dot_product(e1, a1)/(norm2(e1)*norm2(a1)))
   Rot = reshape((/cos(phirad),sin(phirad),0.0d0,-sin(phirad),cos(phirad),0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))
   a1Rot = matmul(Rot,a1)
   a2Rot = matmul(Rot,a2)
   a3Rot = matmul(Rot,a3)
   !print*, a1, a2, a3
   !print*, a1Rot, a2Rot, a3Rot
   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
   write(u,'(a,i0)') 'NumberOfAtoms  ', sz1+sz2
   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz2)
   write(u,'(a)') 'LatticeConstant  1.00 Ang'
   write(u,'(a)') '%block LatticeVectors'
   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
   write(u,'(a)') '%endblock LatticeVectors'
   write(u,*)
   write(u2,*)
   xlo = 0.0d0
   xhi = a1Rot(1)
   ylo = 0.0d0
   yhi = a2Rot(2)
   zlo = 0.0d0
   zhi = 35.0d0
   xy = a2Rot(1)
   xz = 0.0d0
   yz = 0.0d0
   write(u2,'(a,i0)') 'Lattice vectors:'
   write(u2,'(a,3f16.9)') 'a:', xhi, 0.0d0, 0.0d0
   write(u2,'(a,3f16.9)') 'b:', xy, yhi, 0.0d0
   write(u2,'(a,3f16.9)') 'c:', 0.0d0, 0.0d0, 35.0d0
   !write(u2,'(a,i0)') 'Lattice vectors:'
   !do i=1,3
   !   write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   !end do
   write(u2,*)
   write(u2,'(a)') 'Fractional coordinates:'
   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   write(u5,*)
   write(u5,'(i5,a)') sz1+sz2, ' atoms'
   write(u5,'(i1,a)') 2, ' atom types'
   write(u5,*)
   write(u5,'(2f16.9,a)') xlo, xhi, ' xlo xhi'
   write(u5,'(2f16.9,a)') ylo, yhi, ' ylo yhi'
   write(u5,'(2f16.9,a)') zlo, zhi, ' zlo zhi'
   write(u5,'(3f16.9,a)') xy, xz, yz, ' xy xz yz'
   write(u5,*)
   write(u5,'(a)') ' Masses'
   write(u5,*)
   write(u5,'(i1, 3f16.9)') 1, 12.0107
   write(u5,'(i1, 3f16.9)') 2, 12.0107
   write(u5,*)
   write(u5,'(a)') ' Atoms'
   write(u5,*)
   rz = (z/2.0d0 - h/2.0d0)/z
   j = 0
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xs(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 2', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i6,a,i6,a,3f16.9,3i3)') j, ' 2 ', 2, ' 0.0 ', Xs(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   rz = (z/2.0d0 + h/2.0d0)/z
   do i=1,sz2
      j = j+1
      write(u,'(3f16.9,i8)') Xo(:,i), rz, num2(mod(i+1,nb2)+1)
      write(u2,'(A8,3f16.9)') 'C', Xo(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xo(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
      !write(u5,'(i5,a,3f16.9,3i3)') j, ' 1', Xs(:,i), rz, 0, 0, 0 
      write(u5,'(i6,a,i6,a,3f16.9,3i3)') j, ' 1 ', 1, ' 0.0 ', Xo(:,i), rz, 0, 0, 0 
      write(u6,'(i5,i8)') j, num1(mod(i+1,nb1)+1)
   end do
   write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   close(u)

end subroutine MoireWrite

subroutine MoireWrite4L(Xs,Xo,Xs2,Xo2, mm,a,h,z,name,sp1,sp2,N1,N2)

   implicit none

   integer, parameter :: u=10, u2=11, u3=12, u4=13

   double precision, intent(in) :: Xs(:,:), Xo(:,:), a, h, z
   double precision, intent(in) :: Xs2(:,:), Xo2(:,:)
   integer, intent(in) :: mm(4), N1(:), N2(:)
   character(*), intent(in) :: name, sp1(:), sp2(:)

   double precision :: ucell(2,2), Rcell(3,3), m1(2,2), m2(2,2)
   double precision :: rz
   integer :: sz1, sz2, i, j, k, nb1, nb2
   character(60), pointer :: Sp(:)
   integer, pointer :: Nel(:), num1(:), num2(:)

   open(u,FILE=name,STATUS='replace')
   open(u2,FILE="BLBL.frac",STATUS='replace')
   open(u3,FILE="BLBL.basis1",STATUS='replace')
   open(u4,FILE="BLBL.basis2",STATUS='replace')
   ucell(:,1) = (/a,0.0d0/)
   ucell(:,2) = (/a/2.0d0,sqrt(3.0d0)*a/2.0d0/)
   m1(:,1) = (/mm(1),mm(2)/)
   m1(:,2) = (/-mm(2),mm(1)+mm(2)/)
   m2(:,1) = (/mm(3),mm(4)/)
   m2(:,2) = (/-mm(4),mm(3)+mm(4)/)
   Rcell = 0.0d0
   Rcell(1:2,1:2) = matmul(ucell,m1)
   Rcell(3,3) = z
   sz1 = size(Xs,2)
   sz2 = size(Xo,2)
   nb1 = size(sp1)
   nb2 = size(sp2)
   allocate(Sp(nb1+nb2))
   allocate(Nel(nb1+nb2))
   allocate(num1(nb1))
   allocate(num2(nb2))
   Sp(1) = sp1(1)
   Nel(1) = N1(1)
   k = 1
   num1(1) = 1
out1:do i=2,nb1
      do j=1,k
         if (sp1(i)==Sp(j)) then
            num1(i) = j
            cycle out1
         end if
      end do
      k = k+1
      num1(i) = k
      Sp(k) = sp1(i)
      Nel(k) = N1(i)
   end do out1
out2:do i=1,nb2
      do j=1,k
         if (sp2(i)==Sp(j)) then
            num2(i) = j
            cycle out2
         end if
      end do
      k = k+1
      num2(i) = k
      Sp(k) = sp2(i)
      Nel(k) = N2(i)
   end do out2
   write(u,'(a,i0)') 'NumberOfSpecies', k
   write(u,'(a)') '%block ChemicalSpeciesLabel'
   do i=1,k
      write(u,'(2i8,4x,a)') i, Nel(i), Sp(i)
   end do
   write(u,'(a)') '%endblock ChemicalSpeciesLabel'
   write(u,'(a,i0)') 'NumberOfAtoms  ', sz1+sz2
   write(u2,'(a,i0)') 'Number of atoms: ', (sz1+sz2)*2
   write(u,'(a)') 'LatticeConstant  1.00 Ang'
   write(u,'(a)') '%block LatticeVectors'
   write(u,'(3f16.9)') ((Rcell(j,i),j=1,3),i=1,3)
   write(u,'(a)') '%endblock LatticeVectors'
   write(u,*)
   write(u2,*)
   write(u2,'(a,i0)') 'Lattice vectors:'
   do i=1,3
      write(u2,'(a,3f16.9)') 'a:', Rcell(:,i)
   end do
   write(u2,*)
   write(u2,'(a)') 'Fractional coordinates:'
   write(u,'(a)') 'AtomicCoordinatesFormat     Fractional'
   write(u,'(a)') '%block AtomicCoordinatesAndAtomicSpecies'
   rz = (z/2.0d0 - 3/2*h)/z ! L1
   j = 0
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs2(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xs2(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs2(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
   end do
   rz = (z/2.0d0 - h/2.0d0)/z ! L2
   do i=1,sz1
      j = j+1
      write(u,'(3f16.9,i8)') Xs(:,i), rz, num1(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xs(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xs(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
   end do
   rz = (z/2.0d0 + h/2.0d0)/z ! L3
   do i=1,sz2
      j = j+1
      write(u,'(3f16.9,i8)') Xo(:,i), rz, num2(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xo(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xo(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 1', ' &'
   end do
   rz = (z/2.0d0 + 3/2*h)/z ! L4
   do i=1,sz2
      j = j+1
      write(u,'(3f16.9,i8)') Xo2(:,i), rz, num2(mod(i+1,nb1)+1)
      write(u2,'(A8,3f16.9)') 'C', Xo2(:,i), rz
      write(u3,'(a,3f16.9,a)') 'basis ', Xo2(:,i), rz, ' &'
      write(u4,'(a,i5,a,a)') 'basis ', j, ' 2', ' &'
   end do
   write(u,'(a)') '%endblock AtomicCoordinatesAndAtomicSpecies'

   close(u)

end subroutine MoireWrite4L

function gcd(i1,i2)

    integer :: gcd
    integer, intent(in) :: i1,i2
    integer :: a,b

    a = i1
    b = i2
    if (a > b) then
        gcd = a
        a = b
        b = gcd
    end if

    do
      gcd = mod(a, b)
      if (gcd == 0) exit
      a = b
      b = gcd
    end do

    gcd = b

end function gcd

end program GenMoire
