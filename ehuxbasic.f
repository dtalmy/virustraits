      subroutine ehuxbasic(m,n,params,inits,hnt,vnt,htime,vtime,light)
      implicit none
      integer k,l,m,n
      logical light
      real*8 phi,muh,lambda,beta,loss,attach,growth
      real*8 S,I,V,params(4),inits(3)
      real*8 dt,t,dSdt,dIdt,dVdt
      real*8 hnt(m),htime(m),vnt(n),vtime(n)
      phi = exp(params(1))
      muh = exp(params(2))
      lambda = exp(params(3))
      beta = exp(params(4))
      S = inits(1)
      I = inits(2)
      V = inits(3)
      dt=900.0/86400.0
      loss = 0.1
      t = 0.0
      do l=1,1100
         do k=1,m
            if (abs(t+0.5*dt-htime(k)/24.0) < dt) then
               hnt(k)=S+I
            endif
         enddo
         do k=1,n
            if (abs(t+0.5*dt-vtime(k)/24.0) < dt) then
               vnt(k)=V
            endif
         enddo
         attach = phi
         growth = muh
         if (light .eqv. .true.) then
            if (modulo(t*24.0,24.0) .gt. 14.0) then
                attach = 0.0
            endif
         else
            if (t*24.0 .gt. 14.0) then
                attach = 0.0
            endif
            if (t*24.0 .gt. 30.0) then
                growth = 0.0
            endif
         endif
         t = t+dt
         dSdt = muh*S - phi*S*V
         dIdt = phi*S*V-lambda*I
         dVdt = beta*lambda*I - phi*S*V - loss*V
         S = max(S + dt * dSdt, 1e-30)
         I = max(I + dt * dIdt, 1e-30)
         V = max(V + dt * dVdt, 1e-30)    
      end do
      end
        
      
