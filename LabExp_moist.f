c     Mani Rajagopal
c     Oct 1 2023
c     I copied this program from Kamal's version of mrbc 
c     and made the following changes
c     1) Make bottom plate as 0 and top plate as  1 as opposed to 
c        Kamal's version that treats bottom plate as 1
c        and top plate as 0. It didn't affect results except that 
c        plots need to be flipped
c        1 would be bottom and 0 would be top  
c     2) uses temperature in Kelvin to compute the density 
c         and RBC parameters as opposed to Kamal's version that 
c         uses Celsius. It impact is small on the profiles of 
c         scalar's mean and variance
c     3) added comments
c     4) Run multiple realizations 
c     5) structure input and output into different directories
c     6) output evolution of temperature and covariance of T and qv

      program LabExp
       implicit none
       integer M, L, idum, iy, Length
       integer N, Lo, Lp, Lm,j, fhdl_eddy, fhdl_etime
       integer NL(100000), ir(97), num_args, realz_int
       character(len=32), dimension(:), allocatable :: args
       character(len=:), allocatable :: input_dir,cmd_str,exp_name
       character(len=:), allocatable :: realization
       double precision Nt, Np, Na, Nd, Co, Cm, uK, vK, wK, PE
       double precision pmax, p, prob, pa, Sc
       double precision Tdif, T_o, p_atm, s_limit, s_m, pv_o, pvdif
       double precision dt, time, td, tmax, te, ts, random
       double precision u(100000), ua(100000), ur(100000), eu(100000)
       double precision v(100000), va(100000), vr(100000), ev(100000)
       double precision w(100000), wa(100000), wr(100000), ew(100000)
       double precision T(100000), Ta(100000), Tr(100000), eT(100000)
       double precision pv(100000), pva(100000), pvr(100000),
     +                  epv(100000)
       double precision s(100000), sa(100000), sr(100000),
     +                  esat(100000)
       double precision PL(100000), a(14), pdf(5,1000),
     +                  pdf_pv(5,1000), pdf_s(5,1000)
       double precision Tvar_evol(100000), time_evol(100000)
       double precision T_pv_covar(100000)
        

c      Command line arguments      
       num_args = command_argument_count()
       allocate(args(num_args))
       do j = 1, num_args
         call get_command_argument(j,args(j))
       end do
       
c      get input_dir, experiment name and realization value     
       if (num_args   < 3)then
           write(*,*) "Incorrect number of arguments"
           return
       else
           input_dir      = trim(args(1))
           exp_name       = trim(args(2)) 
           realization    = trim(args(3)) 
       endif     
       read(realization,*) realz_int
       write(*,'(6A)') "input:", input_dir," exp:",exp_name,
     +     " realization:", realization
       
c      change to input directory 
       cmd_str           = "input/"//input_dir
       call chdir(cmd_str)
       call system("pwd")

 10    format(5g16.6)
 20    format(g20.1,3g12.3,2i8)
       pmax = 0.1d0
       Nd   = 0.d0

c      read model params and initialization values      
       call readpar(idum,N,Lo,Lp,Lm,dt,td,tmax,a,Sc,p_atm,
     +               T_o,Tdif,s_limit,pv_o,pvdif)
       idum   = idum + realz_int;
       write(*,"(A,I8)") "Random seed Idum:", idum
       
       call init(N,u,v,w,T,pv,s)
       call zeroparam(NL,Nt,Np,Na,pdf,pdf_pv,pdf_s,time,te,ts,s_m)
       call zerovars(N,ua,ur,eu,va,vr,ev,wa,wr,ew,Ta,
     +               pva,Tr,pvr,eT,epv,sa,sr,esat,T_pv_covar)
       call LenProb(Lo,Lm,PL,a(14),Co,Cm)
       
c      change to parent directory       
       call chdir("../../") 
       write(*,*) "now at parent dir:"
       call system("pwd");
      
       cmd_str   = "mkdir -p output/"//exp_name//"/"//realization
       write(*,*) "cmd: ",cmd_str
       call system(cmd_str)
       cmd_str  = "cp input/"//input_dir//"/LabExppar.dat"
       cmd_str  =cmd_str//" output/"//exp_name//"/"//realization//"/"
       write(*,*) "cmd:",cmd_str
       call system(cmd_str)
       cmd_str   = "output/"//exp_name//"/"//realization
c      excecuting cd did not work. so using chdir       
       call chdir(cmd_str)
       call system("pwd")
       
       open(101, file="EddyDiagram.dat",status="unknown")
       open(102, file="Warnings.dat", status="unknown")
       fhdl_eddy = 103
       open(fhdl_eddy,file = "eddy.txt")
       fhdl_etime = 104
       open(fhdl_etime,file = "eddy_time.txt")
       write(6,*) "tmax,t", tmax, dt

c      Mani:Use random.dat to find random number generator's probability distribution
c       open(1001,file="random.dat");
c       do j = 1, 100000
c        write(1001,*) random(idum,iy,ir)
c       end do
c       close(1001);

c      begin ODT simulation      
       do while (time .le. tmax)
        time = time + dt
        Nt = Nt + 1.d0
c        write(*,*) "Nt:", Nt, " time:", time
c       track time and time step
        if (mod(Nt,1000000.0) == 0) then
           write(*,'(A24, F12.0, G12.2)') "Nt, time:",        Nt, time
           write(*,'(A24, 4F12.0)') "Nt, Np, Nd, Na:",  Nt, Np, Nd, Na
           write(*,'(A24, 2G12.2)') "tmax, dt:",        tmax, dt
           write(*,*)
           
c          exit
        endif

c       call diffusion  and compute stats
        if ((time-te) .ge. td) then
         call vis(N,u,v,w,time-te,a(9),a(10))
         call difT(N,T,pv,time-te,a(4),Sc)
         call compsup(N,T,pv,s,Tdif,T_o,p_atm)

         call statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,
     +                 time-te)
         call statsdens(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr,eT
     +                  ,epv,esat,time-te)
         call statscovar(N,time-te,T,pv,T_pv_covar); 
         call statsevol(N,Nd+1,Na,time,T,Tvar_evol,time_evol);
         call comppdf(N,T,pdf,time-te)
         call comppdf(N,pv,pdf_pv,time-te)
         call comppdfsup(N,s,pdf_s,time-te,s_limit,T_o,Tdif,s_m)
         
         Nd = Nd + 1.d0
         te = time
        endif

c       get random eddy size, L
        L = 3*Length(PL,a(14),Co,Cm,idum,iy,ir)
c       get random position
        M = 1 + int(random(idum,iy,ir)*(N-L))
c       probability for eddy of given size and is it energetically possible (buoyancy and kolmogorov scales)
        p = dt*prob(N,M,L,u,v,w,T,pv,a,uK,vK,wK,PE,T_o,Tdif,pv_o,pvdif)
c       if real i.e gt 0
        if (p .gt. 0.d0) then
c       if probability gt pmax then change the dt by same factor
         if (p .gt. pmax) then
          write(102,20) "Time step warning: ", time, dt, p, M, L
          dt = dt*pmax/p
          p = pmax
         endif
         pa = pa + p
         Np = Np + 1.d0
        endif

c       eddy is not called at every time step; so some iterations are wasted
c       w/o call to eddy or diffusion
        if (random(idum,iy,ir) .lt. p) then
c         write(6,*) "random value < p", p
c         call energy(N,M,M+L,a(2),u,v,w,T)
         call eddy(N,M,L,u,v,w,T,pv,uK,vK,wK,PE,a(12))
         write(fhdl_eddy,'(10000E16.4)') T(1:N)
c         call energy(N,M,M+L,a(2),u,v,w,T)
         call vis(N,u,v,w,time-te,a(9),a(10))
         call difT(N,T,pv,time-te,a(4),Sc)
         write(fhdl_eddy,'(10000E16.4)') T(1:N)
         write(fhdl_etime,'(2E16.4)') time, time-te
         
         call compsup(N,T,pv,s,Tdif,T_o,p_atm)
         call statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time-te)
         call statsdens(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr,eT
     +                  ,epv,esat,time-te)
         call statscovar(N,time-te,T,pv,T_pv_covar); 
         call statsevol(N,Nd,Na+1,time,T,Tvar_evol,time_evol);
         call comppdf(N,T,pdf,time-te)
         call comppdf(N,pv,pdf_pv,time-te)
         call comppdfsup(N,s,pdf_s,time-te,s_limit,T_o,Tdif,s_m)

         write(101,10) time, 1.d0*(M+(L/2))/(1.d0*N),
     +                (1.d0*L)/(2.d0*N)
         Na = Na + 1.d0
         NL(L/3) = NL(L/3) + 1
         te = time
        endif

        if (Np > 1.d4) then
          call raisedt(Np,dt,pa)
        endif
         call flush()
       enddo
       
       write (*,*) "outside loop"
       call vis(N,u,v,w,time-te,a(9),a(10))
       call difT(N,T,pv,time-te,a(4),Sc)
       call compsup(N,T,pv,s,Tdif,T_o,p_atm)

       call statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time-te)
       call statsdens(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr,eT
     +               ,epv,esat,time-te)
       call statscovar(N,time-te,T,pv,T_pv_covar);
       call statsevol(N,Nd,Na,time,T,Tvar_evol,time_evol);
       call comppdf(N,T,pdf,time-te)
       call comppdf(N,pv,pdf_pv,time-te)
       call comppdfsup(N,s,pdf_s,time-te,s_limit,T_o,Tdif,s_m)

       close(101)
       close(102)

       call outputstats(N,Na,Nt,Lo,Lm,NL,PL,pdf,pdf_pv,pdf_s,
     +                  time,dt,s_limit,s_m)
       call saveveldata(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,time)
       call savedensdata(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr
     +                  ,eT,epv,esat,time)
       call savecovar(N,time,T_pv_covar,Ta,Tr,pva,pvr);
       call saveevoldata(Nd,Na,Tvar_evol,time_evol)
       write (*,*) "Save data completed"
       stop
       end

       function prob(N,M,L,u,v,w,T,pv,a,uK,vK,wK,PE,T_o,Tdif,pv_o,pvdif)
       implicit none
       integer N, M, L, j
       double precision prob, u(100000), v(100000), w(100000)
       double precision T(100000), pv(100000), a(14), uK, vK,
     +                  wK, PE, KE
       double precision Tv(100000)
       double precision TK, p, x, psiK
       double precision Tdif, T_o, pv_o, pvdif
c       bot as 0 top as 1 change
        double precision T_b,pv_b,T_t,pv_t,T_cell,pv_cell
        double precision Tv_b,Tv_t,Tvdif
        integer Tvoob

       x       = (1.d0*L)/(1.d0*N)

       T_b     = T_o  + Tdif/2.0d0
       T_t     = T_o  - Tdif/2.0d0
       pv_b    = pv_o + pvdif/2.0d0
       pv_t    = pv_o - pvdif/2.0d0
       Tv_b    = (T_b+273.15) * (1.0d0 + pv_b/0.622d0)/(1.0d0+pv_b)
       Tv_t    = (T_t+273.15) * (1.0d0 + pv_t/0.622d0)/(1.0d0+pv_t)
       Tvdif   = (Tv_b - Tv_t) 
       Tvoob = 0

        do j = 1, N
c         bot as 0 top as 1 change
c          T_cell  = T_b  - T(j) * Tdif 
c          pv_cell = pv_b - pv(j) * pvdif
c          compute virtual temperature
c          Tv(j)   = (T_cell+ 273.15)*(1.0d0 + pv_cell/0.622d0)/
c     +              (1.0d0+pv_cell)
c         convert to value between 0 and 1
c          Tv(j)   = (Tv_b - Tv(j))/Tvdif
c         Find if Tv can be negative if Tv > Tv_b

c          Consider bot as 1 and top as 0
c          Tv(j)=((((T(j)-0.5d0)*Tdif+T_o)*
c     +          (1.d0+((pv(j)-0.5d0)*pvdif+pv_o)/0.622d0)
c     +          /(1.d0+((pv(j)-0.5d0)*pvdif+pv_o)))-T_o+Tdif*0.5d0)/Tdif

c          Consider bot as 0 and top as 1
          Tv(j) =(T_b - ((T_b - T(j)*Tdif)*
     +           (1.d0+ (pv_b - pv(j)*pvdif)/0.622d0)
     +          /(1.d0+ (pv_b - pv(j)*pvdif))))/Tdif

          if (Tv(j) < 0 .or. Tv(j) > 1) then
            !write(*,'(A,F12.8)') "Note: Tv < 0: ", Tv(j)
            Tvoob  = Tvoob + 1
          endif  

        enddo


       uK = psiK(N,M,L,u)
       vK = psiK(N,M,L,v)
       wK = psiK(N,M,L,w)
c      bot 0 and top 1 change ; remove "-" sign
       TK = -psiK(N,M,L,Tv)
       PE = a(2)*TK*x
       KE = ((1.d0 - a(12))*wK*wK) + (a(12)*((uK*uK)
     +      + (vK*vK))/2.d0)
       p = ((KE + PE)*x*x) - a(11)

       if (p .gt. 0.d0) then
        prob = a(1)*(1.d0-x)*dsqrt(p)*dexp(3.d0*a(14)/(x*N))/(x*x)
c        write(*, '(A,I6)') "Tv out of bounds:", Tvoob
c        write(*, '(A,F12.8)') "Tv min:", minval(Tv)
c        write(*, '(A,F12.8)') "Tv max:", maxval(Tv)
       else
        prob = 0.d0
       endif
       return
       end

       function psiK(N,M,L,psi)
       implicit none
       integer N, M, L, j
       double precision psiK, sum, z
       double precision psi(100000)
       sum = psi(M)*L/2.d0
       do j=1, L-1
        z = L -( 2*j)
        sum = sum + (psi(j+M)*z)
       enddo
       sum = sum - (psi(M+L)*L/2.d0)
       psiK = 4.d0*sum/(9.d0*L*L)
       return
       end

       subroutine eddy(N,M,L,u,v,w,T,pv,uK,vK,wK,PE,c)
       implicit none
       integer N, M, L
       double precision u(100000), v(100000), w(100000)
       double precision T(100000), pv(100000), uK, vK, wK, PE, c
       double precision qu, qv, qw, cu, cv, cw
       qu = dsqrt((c*((vK*vK) + (wK*wK))/2.d0) + ((1.d0-c)*uK*uK))
       if (uK .gt. 0.d0) then
        cu = 6.75d0*(qu - uK)/(1.d0*L)
       else
        cu = -6.75d0*(qu + uK)/(1.d0*L)
       endif
       qv = dsqrt((c*((uK*uK) + (wK*wK))/2.d0) + ((1.d0-c)*vK*vK))
       if (vK .gt. 0.d0) then
        cv = 6.75d0*(qv - vK)/(1.d0*L)
       else
        cv = -6.75d0*(qv + vK)/(1.d0*L)
       endif
       qw = dsqrt((c*((uK*uK) + (vK*vK))/2.d0) + ((1.d0-c)*wK*wK)
     +             + PE)
       if (wK .gt. 0.d0) then
        cw = 6.75d0*(qw - wK)/(1.d0*L)
       else
        cw = -6.75d0*(qw + wK)/(1.d0*L)
       endif
       call triplet(N,M,L,u)
       call triplet(N,M,L,v)
       call triplet(N,M,L,w)
       call triplet(N,M,L,T)
       call triplet(N,M,L,pv)
       call addK(N,M,L,u,cu)
       call addK(N,M,L,v,cv)
       call addK(N,M,L,w,cw)
       return
       end

       subroutine triplet(N,M,L,psi)
       implicit none
       integer N, M, L, Lo, j, k
       double precision psi(100000), x(100000)
       Lo = L/3
       do j = 1, Lo
        k = M + 3*(j-1)
        x(j) = psi(k)
       enddo
       do j=1, Lo
        k = M + L + 1 - (3*j)
        x(j+Lo) = psi(k)
       enddo
       do j = 1, Lo
        k = M + (3*j) - 1
        x(j+Lo+Lo) = psi(k)
       enddo
       do j=1, L
        k = M+j-1
        psi(k) = x(j)
       enddo
       return
       end

       subroutine addK(N,M,L,u,c)
       implicit none
       integer N, M, L, Lo, j, j1, j2, j3
       double precision c, y1, y2, y3
       double precision u(100000)
       Lo = L/3
       do j = 1, Lo
        y1 = -(2.d0*j)
        y2 = (4.d0*(j+Lo)) - (2.d0*L)
        y3 = (2.d0*L) - (2.d0*(j+Lo+Lo))
        j1 = M + j
        j2 = M + j + Lo
        j3 = M + j + Lo + Lo
        u(j1) = u(j1) + (c*y1)
        u(j2) = u(j2) + (c*y2)
        u(j3) = u(j3) + (c*y3)
       enddo
       return
       end

       subroutine difT(N,T,pv,dt,Pr,Sc2)
       implicit none
       integer N, j
       double precision T(100000), pv(100000), x(100000),
     +                   xpv(100000)
       double precision l(100000), lpv(100000), d(100000),
     +                   dpv(100000)
       double precision r(100000), rpv(100000)
       double precision dt, Pr, Sc2, De, Dv
       De = (dt*N*N)/(2.d0*Pr)
       Dv = (dt*N*N)/(2.d0*Sc2)
       l(1) = 0.d0
       d(1) = 1.d0 + (2.d0*De)
       r(1) = -De

       lpv(1) = 0.d0
       dpv(1) = 1.d0 + (2.d0*Dv)
       rpv(1) = -Dv

       do j = 2, N-1
         l(j) = -De
         d(j) = 1.d0 + (2.d0*De)
         r(j) = -De
         lpv(j) = -Dv
         dpv(j) = 1.d0 + (2.d0*Dv)
         rpv(j) = -Dv
       enddo
       d(N) = 1.d0
       l(N) = 0.d0
       r(N) = 0.d0
       dpv(N) = 1.d0
       lpv(N) = 0.d0
       rpv(N) = 0.d0

       x(1) = ((1.d0 - (2.d0*De))*T(1)) + (De*T(2))
       xpv(1) = ((1.d0 - (2.d0*Dv))*pv(1)) + (Dv*pv(2))

       do j = 2, N-1
         x(j) = ((1.d0-(2.d0*De))*T(j)) + (De*(T(j+1) + T(j-1)))
         xpv(j) = ((1.d0-(2.d0*Dv))*pv(j)) + (Dv*(pv(j+1)
     +            + pv(j-1)))
       enddo
       x(N) = T(N)
       xpv(N) = pv(N)
       call tridiagonal(N,l,d,r,x,T)
       call tridiagonal(N,lpv,dpv,rpv,xpv,pv)
       return
       end

       subroutine vis(N,u,v,w,dt,Uo,f)
       implicit none
       integer N, j
       double precision u(100000), v(100000), w(100000),
     +                  dt, Uo, f, De
       double precision l(100000), d(100000), r(100000)
       double precision xu(100000), xv(100000), xw(100000)
       double precision duf(100000), dvf(100000), cft, sft
       De = (dt*N*N)/(2.d0)
       l(1) = 0.d0
       d(1) = 1.d0 + (2.d0*De)
       r(1) = -De
       do j = 2, N-1
         l(j) = -De
         d(j) = 1.d0 + (2.d0*De)
         r(j) = -De
       enddo
       d(N) = 1.d0
       l(N) = 0.d0
       r(N) = 0.d0
       cft = dcos(f*dt) - 1.d0
       sft = dsin(f*dt)
       do j = 1, N-1
         duf(j) = ((u(j) - Uo)*cft) + (v(j)*sft)
         dvf(j) = (v(j)*cft) - ((u(j) - Uo)*sft)
       enddo
       xu(1) = ((1.d0-(2.d0*De))*u(1)) + (De*u(2))
       xv(1) = ((1.d0-(2.d0*De))*v(1)) + (De*v(2))
       xw(1) = ((1.d0-(2.d0*De))*w(1)) + (De*w(2))
       do j = 2, N-1
         xu(j) = ((1.d0-(2.d0*De))*u(j)) + (De*(u(j+1)+u(j-1)))
         xv(j) = ((1.d0-(2.d0*De))*v(j)) + (De*(v(j+1)+v(j-1)))
         xw(j) = ((1.d0-(2.d0*De))*w(j)) + (De*(w(j+1)+w(j-1)))
       enddo
       xu(N) = u(N)
       xv(N) = v(N)
       xw(N) = w(N)
       call tridiagonal(N,l,d,r,xu,u)
       call tridiagonal(N,l,d,r,xv,v)
       call tridiagonal(N,l,d,r,xw,w)
       do j = 1, N-1
         u(j) = u(j) + duf(j)
         v(j) = v(j) + dvf(j)
       enddo
       return
       end

       subroutine tridiagonal(N,l,d,r,x,y)
       implicit none
       integer N, j
       double precision l(100000), d(100000), r(100000)
       double precision x(100000), y(100000), b, g(100000)
       b = d(1)
       y(1) = x(1)/b
       do j = 2, N
         g(j) = r(j-1)/b
         b = d(j) - (l(j)*g(j))
         if (b .eq. 0.d0) pause 'Tridiagonal:  failure'
         y(j) = (x(j) - (l(j)*y(j-1)))/b
       enddo
       do j = N-1, 1, -1
         y(j) = y(j) - (g(j+1)*y(j+1))
       enddo
       return
       end

       subroutine raisedt(Np,dt,p)
       implicit none
       double precision Np, dt, p, pmin
       pmin = 1.d-3
       p = p/Np
       if (p .lt. (pmin/2.d0)) then
         dt = dt*2.d0;
       else
         dt = dt*pmin/p
       endif
       p = 0.d0
       Np = 0.d0
       return
       end

       function Length(PL,xp,Co,Cm,idum,iy,ir)
       implicit none
       integer Length, idum, iy, n
       integer ir(97)
       double precision xp, Co, Cm, r, random, x
       double precision PL(100000)
       r = random(idum,iy,ir)
       x = -xp/dlog((Co*r)+(Cm*(1.d0-r)))
       n = int(x)-1
       if (r .gt. PL(n)) then
 10     n = n + 1
        if (r .gt. PL(n)) goto 10
       endif
       if (r .lt. PL(n-1)) then
 20     n = n - 1
        if (r .lt. PL(n-1)) goto 20
       endif
       Length = n
       return
       end

       subroutine LenProb(Lo,Lm,P,xp,Co,Cm)
       implicit none
       integer Lo, Lm, L
       double precision P(100000), xp, Co, Cm, C, z
       Co = dexp(-xp/(1.d0*Lo))
       Cm = dexp(-xp/(1.d0*Lm))
       C = 0.d0
       do L = Lo, Lm
        z = dexp(-xp/(1.d0*L))*(dexp(xp/(L*(L+1.d0)))-1.d0)
        C = C + z
       enddo
       C = 1.d0/C
       do L = 1, Lo-1
        P(L) = 0.d0
       enddo
       do L = Lo, Lm
        z = dexp(-xp/(1.d0*L))*(dexp(xp/(L*(L+1.d0)))-1.d0)
        P(L) = P(L-1) + (C*z)
       enddo
       do L = Lm+1, 100000
        P(L) = 0.d0
       enddo
       return
       end

       subroutine comppdf(N,T,pdf,dt)
       implicit none
       integer N, j, k(5), i, w
       double precision dt, T(100000), pdf(5,1000), X
c removed '-' in X
       X = 1.d3
       w = 5
       k(1) = N/2
       k(2) = N/4
       k(3) = N/8
       k(4) = 5*w
       k(5) = 2*w
       do j = -w, w
        i = int(X*(T(j+k(1))))
        pdf(1,i) = pdf(1,i) + dt
        i = int(X*(T(j+k(2))))
        pdf(2,i) = pdf(2,i) + dt
        i = int(X*(T(j+k(3))))
        pdf(3,i) = pdf(3,i) + dt
        i = int(X*(T(j+k(4))))
        pdf(4,i) = pdf(4,i) + dt
        i = int(X*(T(j+k(5))))
        pdf(5,i) = pdf(5,i) + dt

       enddo
       return
       end

       subroutine comppdfsup(N,s,pdf,dt,s_limit,T_o,Tdif,s_m)
       implicit none
       integer N, j, k(5), i, w
       double precision dt, s(100000), pdf(5,1000), X
       double precision s_limit, T_o, Tdif, es_b, es_t
        double precision es_m, s_m

       es_b = 6.112*dexp(17.67*(T_o+Tdif/2)/(T_o+Tdif/2+243.5))
       es_t = 6.112*dexp(17.67*(T_o-Tdif/2)/(T_o-Tdif/2+243.5))
       es_m = 6.112*dexp(17.67*T_o/(T_o+243.5))
       s_m = ((es_b+es_t)/(2*es_m)-1)*100

       X = 1.d3
       w = 5
       k(1) = N/2
       k(2) = N/4
       k(3) = N/8
       k(4) = 5*w
       k(5) = 2*w

c      grid cells +/- 5 levels above and below kth cell
       do j = -w, w
c      don't know why the PDF is shifted by 500 position
c      I think this allows to sub-saturation values with negative s values
c      s_limit value is 30 and x is 1000, this allow for supersaturation to be
c      divided in to 33 bins
           i = int(X/2+X*(s(j+k(1))-s_m)/s_limit+1)
           pdf(1,i) = pdf(1,i) + dt

           i = int(X/2+X*(s(j+k(2))-s_m)/s_limit+1)
           pdf(2,i) = pdf(2,i) + dt

           i = int(X/2+X*(s(j+k(3))-s_m)/s_limit+1)
           pdf(3,i) = pdf(3,i) + dt

           i = int(X/2+X*(s(j+k(4))-s_m)/s_limit+1)
           pdf(4,i) = pdf(4,i) + dt

           i = int(X/2+X*(s(j+k(5))-s_m)/s_limit+1)
           pdf(5,i) = pdf(5,i) + dt

       enddo
       return
       end

       subroutine statsvel(N,u,ua,ur,eu,v,va,vr,ev,w,wa,wr,ew,dt)
       implicit none
       integer N, j
       double precision u(100000), ua(100000), ur(100000),
     +                  eu(100000)
       double precision v(100000), va(100000), vr(100000),
     +                  ev(100000)
       double precision w(100000), wa(100000), wr(100000),
     +                  ew(100000)
       double precision dt, dudz, dvdz, dwdz
       do j = 1, N
        ua(j) = ua(j)+(u(j)*dt)
        va(j) = va(j)+(v(j)*dt)
        wa(j) = wa(j)+(w(j)*dt)
        ur(j) = ur(j)+(u(j)*u(j)*dt)
        vr(j) = vr(j)+(v(j)*v(j)*dt)
        wr(j) = wr(j)+(w(j)*w(j)*dt)
       enddo
c      dz = 1/N
       dudz = (u(2) - u(1))*N
       dvdz = (v(2) - v(1))*N
       dwdz = (w(2) - w(1))*N
       eu(1) = eu(1)+(dudz*dudz*dt)
       ev(1) = ev(1)+(dvdz*dvdz*dt)
       ew(1) = ew(1)+(dwdz*dwdz*dt)
       do j = 2, N-1
         dudz = (u(j+1) - u(j-1))*N/2.d0
         dvdz = (v(j+1) - v(j-1))*N/2.d0
         dwdz = (w(j+1) - w(j-1))*N/2.d0
         eu(j) = eu(j)+(dudz*dudz*dt)
         ev(j) = ev(j)+(dvdz*dvdz*dt)
         ew(j) = ew(j)+(dwdz*dwdz*dt)
       enddo
       dudz = (u(N) - u(N-1))*N
       dvdz = (v(N) - v(N-1))*N
       dwdz = (w(N) - w(N-1))*N
       eu(N) = eu(N)+(dudz*dudz*dt)
       ev(N) = ev(N)+(dvdz*dvdz*dt)
       ew(N) = ew(N)+(dwdz*dwdz*dt)
       return
       end

       subroutine statsdens(N,T,pv,s,Ta,pva,sa,Tr,pvr
     +                     ,sr,eT,epv,esat,dt)
       implicit none
       integer N, j
       double precision T(100000), Ta(100000), Tr(100000),
     +                  eT(100000)
       double precision pv(100000), pva(100000), pvr(100000),
     +                  epv(100000)
       double precision s(100000), sa(100000), sr(100000),
     +                  esat(100000)


       double precision dt, delT, delpv, dels
       do j = 1, N
        Ta(j) = Ta(j)+(T(j)*dt)
        Tr(j) = Tr(j)+(T(j)*T(j)*dt)
        pva(j) = pva(j)+(pv(j)*dt)
        pvr(j) = pvr(j)+(pv(j)*pv(j)*dt)
        sa(j) = sa(j)+(s(j)*dt)
        sr(j) = sr(j)+(s(j)*s(j)*dt)
       enddo
       delT = T(2) - T(1)
       delpv = pv(2) - pv(1)
       dels  = s(2) - s(1)
       eT(1) = eT(1)+(delT*delT*N*N*dt)
       epv(1) = epv(1)+(delpv*delpv*N*N*dt)
       esat(1) = esat(1)+(dels*dels*N*N*dt)
       do j = 2, N-1
         delT = T(j+1) - T(j-1)
         eT(j) = eT(j)+(delT*delT*N*N*dt/4.d0)
         delpv = pv(j+1) - pv(j-1)
         epv(j) = epv(j)+(delpv*delpv*N*N*dt/4.d0)
         dels = s(j+1) - s(j-1)
         esat(j) = esat(j)+(dels*dels*N*N*dt/4.d0)

       enddo
       delT = T(N) - T(N-1)
       eT(N) = eT(N)+(delT*delT*N*N*dt)
       delpv = pv(N) - pv(N-1)
       epv(N) = epv(N)+(delpv*delpv*N*N*dt)
       dels = s(N) - s(N-1)
       esat(N) = esat(N)+(dels*dels*N*N*dt)

       return
       end
       
       subroutine statscovar(N,dt,T,pv,T_pv_covar)
          implicit none
          integer N,j
          double precision T(100000), pv(100000),T_pv_covar(100000)
          double precision dt
          
          do j = 1, N
              T_pv_covar(j)  = T_pv_covar(j) + T(j)*pv(j)*dt
          enddo
          
       return
       end
           
       subroutine statsevol(N,Nd,Na,time,T,Tvar_evol,time_evol)
       implicit none
       integer N, j, k, tstep_evol
       double precision T(100000),Tvar_evol(100000),time_evol(100000)
       double precision Tmean,Tvar, time, Na, Nd
         
         tstep_evol   = int(Na + Nd)
         Tmean        = sum(T(1:N))/N
         Tvar         = sum( (T(1:N)-Tmean)**2 )/N
         Tvar_evol(tstep_evol) = Tvar
         time_evol(tstep_evol)=time
c        write(*,*) "tstep_evol,time:", tstep_evol, time
c        write(*,*) "Tmean, Tvar", Tmean, Tvar
         
       return 
       end
       
       subroutine saveevoldata(Nd,Na,Tvar_evol,time_evol)
       implicit none
       double precision Tvar_evol(100000),time_evol(100000)
       double precision Nd,Na
       integer tstep_tot,j
       
       
       tstep_tot          = int(Nd+Na)
       
       open(100, file="evol_stats.dat", status="unknown")  
       
       do j = 1, tstep_tot
        write(100,'(2g12.4)') time_evol(j), Tvar_evol(j)
       enddo
       
       close(100)       
          
       
       return
       end
       
       subroutine saveveldata(N,u,ua,ur,eu,v,va,vr,ev,w
     +                        ,wa,wr,ew,tf)
       implicit none
       integer N, j, k
       double precision u(100000), ua(100000), ur(100000),
     +                  eu(100000)
       double precision v(100000), va(100000), vr(100000),
     +                  ev(100000)
       double precision w(100000), wa(100000), wr(100000),
     +                  ew(100000)
       double precision tf, z, small, dis, Re, etot
              
10     format(5g16.6)
20     format(g30.1,g15.3)
       small = 1.d-10
       dis = 0.d0
       do j = 1, N
        ua(j) = ua(j)/tf
        va(j) = va(j)/tf
        wa(j) = wa(j)/tf
        ur(j) = (ur(j)/tf)-(ua(j)*ua(j))
        vr(j) = (vr(j)/tf)-(va(j)*va(j))
        wr(j) = (wr(j)/tf)-(wa(j)*wa(j))
        eu(j) = eu(j)/tf
        ev(j) = ev(j)/tf
        ew(j) = ew(j)/tf
        dis = dis + eu(j) + ev(j) + ew(j)
       enddo
       dis = dis/(1.d0*N)
       Re = 0.d0
       k = 3*N/8
       do j = 1, N/4
        Re = Re + ur(j+k) + vr(j+k) + wr(j+k)
       enddo
       Re = dsqrt(4.d0*Re/(1.d0*N))
       write(6,20) "Core Reynolds Number, ReC: ", Re
       write(6,20) "Total KE dissipation:      ", dis
       write(6,20) "Wall stress, u:            ", ua(1)*N
       write(6,20) "                           ", (ua(N)-ua(N-1))*N
       write(6,20) "Wall stress, v:            ", va(1)*N
       write(6,20) "                           ", (va(N)-va(N-1))*N
       write(6,20) "Wall stress, w:            ", wa(1)*N
       write(6,20) "                           ", (wa(N)-wa(N-1))*N
       write(6,20)

       open(100, file="U.dat", status="unknown")
       open(200, file="V.dat", status="unknown")
       open(300, file="W.dat", status="unknown")
       open(400, file="eu.dat", status="unknown")           
       
       do j = 1, N
        z = (1.d0*j)/(1.d0*N)
        etot = eu(j) + ev(j) + ew(j)
        if (abs(ua(j)) .lt. small) ua(j) = 0.d0
        if (abs(va(j)) .lt. small) va(j) = 0.d0
        if (abs(wa(j)) .lt. small) wa(j) = 0.d0
        if (ur(j) .lt. 0.d0) ur(j) = 0.d0
        if (vr(j) .lt. 0.d0) vr(j) = 0.d0
        if (wr(j) .lt. 0.d0) wr(j) = 0.d0
        write(100,10) z, u(j), ua(j), dsqrt(ur(j))
        write(200,10) z, v(j), va(j), dsqrt(vr(j))
        write(300,10) z, w(j), wa(j), dsqrt(wr(j))
        write(400,10) z, eu(j), ev(j), ew(j), etot
       enddo
       
       close(100)
       close(200)
       close(300)
       close(400)
       return
       end

       subroutine savedensdata(N,T,pv,s,Ta,pva,sa,Tr,pvr,sr
     +                         ,eT,epv,esat,tf)
       implicit none
       integer N, j, k
       double precision T(100000), Ta(100000), Tr(100000),
     +                  eT(100000)
       double precision pv(100000), pva(100000), pvr(100000),
     +                  epv(100000)
       double precision s(100000), sa(100000), sr(100000),
     +                  esat(100000)

       double precision tf, z, small, dis, dis_pv, x, x_pv
      
10     format(5g16.6)
20     format(g30.1,g15.4)

       small = 1.d-10
       dis = 0.d0
       dis_pv = 0.d0
       do j = 1, N
        Ta(j) = Ta(j)/tf
        Tr(j) = (Tr(j)/tf)-(Ta(j)*Ta(j))
        eT(j) = eT(j)/tf
        dis = dis + eT(j)
        pva(j) = pva(j)/tf
        pvr(j) = (pvr(j)/tf)-(pva(j)*pva(j))
        epv(j) = epv(j)/tf
        dis_pv = dis_pv + epv(j)
        sa(j) = sa(j)/tf
        sr(j) = (sr(j)/tf)-(sa(j)*sa(j))
        esat(j) = esat(j)/tf

       enddo

       dis = dis/(1.d0*N)
       dis_pv = dis_pv/(1.d0*N)
       x = 0.d0
       x_pv = 0.d0
       k = 3*(N/8)
       do j = 1, N/4
        x = x + Tr(j+k)
        x_pv = x_pv + pvr(j+k)
       enddo
       x = dsqrt(4.d0*x/(1.d0*N))
       x_pv = dsqrt(4.d0*x_pv/(1.d0*N))

       open(100, file="T.dat", status="unknown")
       open(200, file="pv.dat", status="unknown")
       open(300, file="Nu.dat", status="unknown")
       open(400, file="s.dat", status="unknown")
       
       write(300,20) "Core Temperature and pv Fluct.:   ", x, x_pv
       write(300,20) "Total Temperature and pv dissip.: ", dis, dis_pv
       write(300,20) "Nusselt Number:            ", Ta(1)*N
       write(300,20) "                           ", (Ta(N)-Ta(N-1))*N
       write(300,20) "Sherwood Number:            ", pva(1)*N
       write(300,20) "                         ", (pva(N)-pva(N-1))*N

       
       do j = 1, N
        z = (1.d0*DBLE(j))/(1.d0*N)
        if (abs(Ta(j)) .lt. small) Ta(j) = 0.d0
        if (Tr(j) .lt. 0.d0) Tr(j) = 0.d0
        if (abs(pva(j)) .lt. small) pva(j) = 0.d0
        if (pvr(j) .lt. 0.d0) pvr(j) = 0.d0
        write(100,10) z, T(j), Ta(j), dsqrt(Tr(j))
        write(200,10) z, pv(j), pva(j), dsqrt(pvr(j))
        if (abs(sa(j)) .lt. small) sa(j) = 0.d0
        if (sr(j) .lt. 0.d0) sr(j) = 0.d0
        write(400,10) z, s(j), sa(j), dsqrt(sr(j))

       enddo


       close(100)
       close(200)
       close(300)
       close(400)

      return
       end

       subroutine savecovar(N,tf,T_pv_covar,Ta,Tr,pva,pvr)
        integer N,j
        double precision tf,Z,numer,denom
        double precision T_pv_covar(100000)
        double precision Ta(100000), Tr(100000)
        double precision  pva(100000),pvr(100000)
           
        do j = 1, N
          numer        = (T_pv_covar(j)/tf)-(Ta(j)*pva(j)) 
          denom        =  dsqrt(Tr(j))*dsqrt(pvr(j))
          T_pv_covar(j)= numer/denom;
        end do

        open(100, file="T_qv_covar.dat")      
        do j = 1, N
          z = (1.d0*DBLE(j))/(1.d0*N);
          write(100,'(2g16.6)') z,T_pv_covar(j)
        enddo     
        close(100)
           
       return
       end
       
       subroutine outputstats(N,Na,Nt,Lo,Lm,NL,PL,pdf,
     +                  pdf_pv,pdf_s,tf,dt,s_limit,s_m)

       implicit none
       integer N, NL(100000), Lo, Lm, j
       double precision Na, Nt, tf, dt, x, xo, x_sup
       double precision PL(100000), pdf(5,1000), pdf_pv(5,1000)
       double precision pdf_s(5,1000), s_limit, s_m
       
       
 10    format(6g16.8)
 20    format(g30.1,g15.4)
 
       open(100, file="PDF_T.dat", status="unknown")
       open(200, file="PDF_pv.dat", status="unknown")
       open(300, file="PDF_s.dat", status="unknown")
       open(500, file="PL.dat", status="unknown")
       open(600, file="out_ODT.dat", status="unknown")
       
       do j = Lo, Lm
        if (NL(j) .gt. 0) then
         x = (1.d0*NL(j))/Na
         xo = PL(j) - PL(j-1)
         write(500,10) dlog((3.d0*j)/(1.d0*N)), dlog(x), dlog(xo)
        endif
       enddo
       

       write(600,20) "Final Time Step:           ", dt
       write(600,20) "Average Eddy Time Step:    ", tf/Nt
       write(600,20) "Time Between Eddies:       ", tf/Na
       write(600,20) "Total Number of Eddies:    ", Na
       write(600,20) "Eddie Acceptance Rate:     ", Na/Nt
       write(600,20)


       do j = 1 , 1000
        x = j/1.d3
c       removed '-' above
        pdf(1,j) = pdf(1,j)/tf
        pdf(2,j) = pdf(2,j)/tf
        pdf(3,j) = pdf(3,j)/tf
        pdf(4,j) = pdf(4,j)/tf
        pdf(5,j) = pdf(5,j)/tf
        write(100,10) x, pdf(1,j), pdf(2,j), pdf(3,j), pdf(4,j),
     +              pdf(5,j)
        pdf_pv(1,j) = pdf_pv(1,j)/tf
        pdf_pv(2,j) = pdf_pv(2,j)/tf
        pdf_pv(3,j) = pdf_pv(3,j)/tf
        pdf_pv(4,j) = pdf_pv(4,j)/tf
        pdf_pv(5,j) = pdf_pv(5,j)/tf
        write(200,10) x, pdf_pv(1,j), pdf_pv(2,j), pdf_pv(3,j),
     +              pdf_pv(4,j), pdf_pv(5,j)

        pdf_s(1,j) = pdf_s(1,j)/tf
        pdf_s(2,j) = pdf_s(2,j)/tf
        pdf_s(3,j) = pdf_s(3,j)/tf
        pdf_s(4,j) = pdf_s(4,j)/tf
        pdf_s(5,j) = pdf_s(5,j)/tf

        x_sup = s_m-s_limit/2.d0+DBLE(j)*s_limit/1.d3
        write(300,10) x_sup, pdf_s(1,j), pdf_s(2,j), pdf_s(3,j),
     +              pdf_s(4,j), pdf_s(5,j)

       enddo
       close(100)
       close(200)
       close(300)
       close(500)
       close(600)

       return
       end

       subroutine init(N,u,v,w,T,pv,s)
       implicit none
       integer N, j
       double precision u(100000), v(100000), w(100000)
       double precision T(100000), pv(100000), s(100000), z
 10    format(5g16.8)
 
       open(100, file="U.dat", status="old")
       open(200, file="V.dat", status="old")
       open(300, file="W.dat", status="old")
       open(400, file="T.dat", status="old")
       open(500, file="pv.dat", status="old")
       open(600, file="s.dat", status="old")
      


       do j = 1, N
        read(100,10) z, u(j)
        read(200,10) z, v(j)
        read(300,10) z, w(j)
        read(400,10) z, T(j)
        read(500,10) z, pv(j)
        read(600,10) z, s(j)
       enddo
       
       close(100)
       close(200)
       close(300)
       close(400)
       close(500)
       close(600)
       
       return
       end

       subroutine zeroparam(NL,Nt,Np,Na,pdf,pdf_pv,pdf_s,time,te,
     +                       ts,s_m)

       implicit none
       integer NL(100000), j

       double precision Nt, Np, Na, pdf(5,1000), pdf_pv(5,1000)
     +                  , pdf_s(5,1000), time, te, ts, s_m
       Nt = 0.d0
       Np = 0.d0
       Na = 0.d0
       time = 0.d0
       te = 0.d0
       ts = 0.d0
       s_m = 0.d0
       do j=1, 1000
        pdf(1,j) = 0.d0
        pdf(2,j) = 0.d0
        pdf(3,j) = 0.d0
        pdf(4,j) = 0.d0
        pdf(5,j) = 0.d0
        pdf_pv(1,j) = 0.d0
        pdf_pv(2,j) = 0.d0
        pdf_pv(3,j) = 0.d0
        pdf_pv(4,j) = 0.d0
        pdf_pv(5,j) = 0.d0
        pdf_s(1,j) = 0.d0
        pdf_s(2,j) = 0.d0
        pdf_s(3,j) = 0.d0
        pdf_s(4,j) = 0.d0
        pdf_s(5,j) = 0.d0
       enddo
       do j=1, 100000
        NL(j) = 0
       enddo
       return
       end

       subroutine zerovars(N,ua,ur,eu,va,vr,ev,wa,wr,ew,
     +                  Ta,pva,Tr,pvr,eT,epv,sa,sr,esat,T_pv_covar)
       implicit none
       integer N, j, n_k
       double precision ua(100000), ur(100000), eu(100000)
       double precision va(100000), vr(100000), ev(100000)
       double precision wa(100000), wr(100000), ew(100000)
       double precision Ta(100000), Tr(100000), eT(100000)
       double precision pva(100000), pvr(100000), epv(100000)
       double precision sa(100000), sr(100000), esat(100000)
       double precision T_pv_covar(100000)

       do j=1, N
        ua(j) = 0.d0
        ur(j) = 0.d0
        eu(j) = 0.d0
        va(j) = 0.d0
        vr(j) = 0.d0
        ev(j) = 0.d0
        wa(j) = 0.d0
        wr(j) = 0.d0
        ew(j) = 0.d0
        Ta(j) = 0.d0
        Tr(j) = 0.d0
        eT(j) = 0.d0
        pva(j) = 0.d0
        pvr(j) = 0.d0
        epv(j) = 0.d0
        sa(j) = 0.d0
        sr(j) = 0.d0
        esat(j) = 0.d0
        T_pv_covar(j) = 0.d0
       enddo
       return
       end

       subroutine readpar(idum,N,Lo,Lp,Lm,dt,td,tmax,a,Sc,p_atm,T_o,Tdif
     +                    ,s_limit,pv_o,pvdif)
       implicit none
       integer idum, N, Lo, Lp, Lm
       double precision dt, tmax, td, a(14), Sc
       double precision p_atm, T_o, Tdif, s_limit, pv_o, pvdif
 10    format(3i12)
 20    format(3d15.4)
 30    format(g30.1,g15.4)
 40    format(g30.1,i15)
 

       open(100,file="LabExppar.dat", status="old")
       
       read(100,10) N
       read(100,20) a(2), a(4), Sc
       read(100,20) a(9), a(10)
       read(100,20) dt, tmax
       read(100,20) a(11), a(12)
       read(100,10) Lo, Lp, Lm
       read(100,10) idum
       read(100,20) p_atm, T_o, Tdif,s_limit, pv_o, pvdif
       close(100)
       a(14) = 2.d0*Lp
       a(1) = dexp(-a(14)/(1.d0*Lm))-dexp(-a(14)/(1.d0*Lo))
       a(1) = a(1)*N/(3.d0*a(14))
       td = 1.d1/(1.d0*N*N)
       if (a(10)*td .gt. 0.1d0) td = 0.1d0/a(10)
       write(6,30) "Dimensionless Parameters:  "
       write(6,30) "   Buoyancy, Temperature:  ", a(2)
       write(6,30) "   Prandtl Number and Sc:  ", a(4), Sc
       write(6,30) "   Wind Speed:             ", a(9)
       write(6,30) "   Coriolis Parameter:     ", a(10)
       write(6,30)
       write(6,30) "Model Parameters:          "
       write(6,30) "   ZC2 =                   ", a(11)
       write(6,30) "   KE mixing ratio =       ", a(12)
       write(6,30) "   a1 =                    ", a(1)
       write(6,30) "   Smallest eddy:      ", (3.d0*Lo)/(1.d0*N)
       write(6,30) "   Most likely eddy:   ", (3.d0*Lp)/(1.d0*N)
       write(6,30) "   Largest eddy:       ", (3.d0*Lm)/(1.d0*N)
       write(6,30)
       write(6,40) "Grid Points:               ", N
       write(6,30) "Total Simulation Time:     ", tmax
       write(6,30) "Diffusive Time Step:       ", td
       write(6,30) "Initial Eddy Time Step:    ", dt
       return
       end

       function random(idum,iy,ir)
       implicit none
       integer idum, iy, m, ia, ic, j
       double precision random, rm
       integer ir(97)
       m = 714025
       ia = 1366
       ic = 150889
       rm = 1.4005112d-6
       if (idum .lt. 0) then
        idum = mod(ic-idum,m)
        do j=1,97
         idum = mod(ia*idum+ic,m)
         ir(j) = idum
        enddo
        idum = mod(ia*idum+ic,m)
        iy = idum
       endif
       j = 1 + (97*iy)/m
       if (j .gt. 97 .or. j .lt. 1) then
        j = mod(j,97)
        if (j .ge. 0) j=-j
        j = j + 1
       endif
       iy = ir(j)
       idum = mod(ia*idum+ic,m)
       ir(j) = idum
       random = iy*rm
       return
       end

       subroutine energy(N,M1,M2,bT,u,v,w,T)
       implicit none
       integer N, M1,M2, j
       double precision bT, Eu, Ev, Ew, ET, EK
       double precision u(100000), v(100000), w(100000), T(100000)
 10    format(i12,3g16.6)
       Eu = 0.d0
       Ev = 0.d0
       Ew = 0.d0
       ET = 0.d0
       do j = M1, M2
        Eu = Eu + (u(j)*u(j))
        Ev = Ev + (v(j)*v(j))
        Ew = Ew + (w(j)*w(j))
        ET = ET + (T(j)*T(j))
       enddo
       Eu = Eu/(2.d0*N)
       Ev = Ev/(2.d0*N)
       Ew = Ew/(2.d0*N)
       EK = Eu + Ev + Ew
       ET = -ET*27.d0*bT/(8.d0*N*N)
       write(6, 10) M2-M1, EK+ET, EK, ET
       return

       end

       subroutine sumtoavg(Eth,Eth_a,Eqv,Eqv_a,Es,Es_a,dt)
       implicit none
       integer N, j, n_k
       double precision Eth(100000), Eth_a(100000)
       double precision Eqv(100000), Eqv_a(100000)
       double precision Es(100000), Es_a(100000)
       double precision fk(100000)
       double precision dt
       n_k =4096*4
       fk(1:100000) = 0.d0
21     format(5g16.6)
       open(755, file="E_spec_sum.dat", status="unknown")

       do j = 1, n_k
        Eth_a(j) = Eth_a(j)+(Eth(j)*dt)
        Eqv_a(j) = Eqv_a(j)+(Eqv(j)*dt)
        Es_a(j) = Es_a(j)+(Es(j)*dt)
        fk(j) = DBLE((j-1)*(2.d0/(n_k*N)))
        write(755,21) fk(j), Eth_a(j), Eqv_a(j), Es_a(j)

       enddo
       close(755)
       return
       end

       subroutine compsup(N,T,pv,s,Tdif,T_o,p_atm)
       implicit none
       integer N,j, xx
       double precision T(100000), pv(100000), s(100000)
       double precision Tdif, T_o, p_atm
       double precision pvT(2)
       double precision es_b, es_t, pv_t, pv_b, pv_sat, T_cell, T_b
 10    format(g16.8)
       pvT(1)= 0.d0
       pvT(2)= 0.d0
       es_b= 0.d0
       es_t= 0.d0
       pv_t= 0.d0
       pv_b= 0.d0
       pv_sat= 0.d0

       es_b = 6.112d0*dexp(17.67d0*(T_o+Tdif/2.d0)
     +        /(T_o+Tdif/2.d0+243.5))*100.d0
       es_t = 6.112d0*dexp(17.67d0*(T_o-Tdif/2.d0)
     +        /(T_o-Tdif/2.d0+243.5))*100.d0
       pv_t = 0.622d0*es_t/(p_atm-es_t)
       pv_b = 0.622d0*es_b/(p_atm-es_b)

       T_b  = T_o+Tdif/2.d0

       do j = 1, N
c      Seems like T and pv have values between 0 and 1;
c      multiplying by Tdiff will scale them back; then add to Ttop, which is 0

c      pv_sat = 6.112d0*dexp(17.67d0*(T(j)*Tdif+T_o-Tdif/2.d0) 
c     +               /(T(j)*Tdif+T_o-Tdif/2+243.5d0))*100.d0

c      bot 0 and top 1 change    
           T_cell =  T_b - T(j)*Tdif
           pv_sat = 6.112d0*dexp(17.67d0*T_cell
     +               /(T_cell+243.5d0))*100.d0

           pv_sat = 0.622d0*pv_sat/(p_atm-pv_sat)
c       any scalar can be obtained by scaling back for values between 0 and 1
       
c       s(j) = ((pv(j)*(pv_b-pv_t)+pv_t)/pv_sat-1.d0)*100.d0
c       bot 0 and top 1 change        
        s(j) = ((pv_b - pv(j)*(pv_b-pv_t))/pv_sat-1.d0)*100.d0

       enddo
c       write(*,*) " "
c       write(*,'(A,3F16.4)') "SS mean, min, max :", 
c     +      (sum(s(1:N))/N), minval(s(1:N)), maxval(s(1:N))
       
      return

      end
