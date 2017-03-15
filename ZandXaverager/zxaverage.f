c copy these averagers to your source code
c these two can be used to average in z and x direction
c usage 
c-----------------------------------------------------------------------
      subroutine my_avg()
      include 'SIZE'
      include 'TOTAL'
c  this example takes the average of vx,vy,vz in homogeneous x and z dir
      real vxa(lx1,ly1,lz1,lelt),vxa2(lx1,ly1,lz1,lelt)     
      real vya(lx1,ly1,lz1,lelt),vya2(lx1,ly1,lz1,lelt)     
      real vza(lx1,ly1,lz1,lelt),vza2(lx1,ly1,lz1,lelt)     
      integer i,j,n,nelxx,nelyy,nelzz,ifld,npts
      integer gs_avg_hndl
      save    gs_avg_hndl
      data    gs_avg_hndl / 0 /
      integer gs_avg_hndl2
      save    gs_avg_hndl2
      data    gs_avg_hndl2 / 0 /
      real uvw(ldim,lx1*ly1*lz1*lelt),xyz(ldim,lx1*ly1*lz1,lelt)

      nelxx = 10
      nelyy = 12
      nelzz = 8
      ifld = 1
cccccccc
c Z AVERAGE
cccccccc
      call z_avg(vxa,vx,gs_avg_hndl,nelxx,nelyy,nelzz,ifld)
      call z_avg(vya,vy,gs_avg_hndl,nelxx,nelyy,nelzz,ifld)
      call z_avg(vza,vz,gs_avg_hndl,nelxx,nelyy,nelzz,ifld)
c  the z averaged velocities are now in vxa,vya and vza arrays
cccccccc
c X AVERAGE
cccccccc
      call xx_avg(vxa2,vxa,gs_avg_hndl2,nelxx,nelyy,nelzz,ifld)
      call xx_avg(vya2,vya,gs_avg_hndl2,nelxx,nelyy,nelzz,ifld)
      call xx_avg(vza2,vza,gs_avg_hndl2,nelxx,nelyy,nelzz,ifld)
c  the x averaged velocities are now in vxa2,vya2 and vza2 arrays

c    now in order to get the values of these velocities along a line
c    we have to first specify the coordinates of the line
      npts = 1000
      do j=1,npts
         xyz(1,j) = 3
         xyz(2,j) = j*0.01
         xyz(3,j) = 1.
      enddo
c     xyz is points along the y-axis (from 0.01 to 10) at x=1,z=3

       call copy(vx,vxa2,lx1*ly1*lz1*lelt)
       call copy(vy,vya2,lx1*ly1*lz1*lelt)
       call copy(vz,vza2,lx1*ly1*lz1*lelt)
       call outpost(vx,vy,vz,pr,t,'   ') 
c      this outputs the x and z avg field
       call interp_v(uvw,xyz,1000)
c      uvw now contains the velocities along the line

      do i=1,npts
       write(6,*) i,uvw(1,i),uvw(2,i),uvw(3,i),'k10xzavg'
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine z_avg(ua,u,gs_avg_hndl,nelxx,nelyy,nelzz,ifld)
      include 'SIZE'
      include 'TOTAL'

c     Compute the z average of quantity u() - assumes global tens.prod.
      real u (lx1,ly1,lz1,lelt)
      real ua(lx1,ly1,lz1,lelt)

      integer gs_avg_hndl,e,ex,ey,ez,eg,nelxy,nelxx,nelyy

      nelxy = nelxx*nelyy

      if (gs_avg_hndl.eq.0) then
          call set_gs_zavg_hndl(gs_avg_hndl,nelxy,ifld)
      endif

      nel = nelfld(ifld)
      n   = nx1*ny1*nz1*nel

      call copy(ua,bm1,n)              ! Set the averaging weights
      call gs_op(gs_avg_hndl,ua,1,1,0) ! Sum weights over columns


      do i=1,n                          ! ua = (w_j*u_j)/( sum_i w_i)
         ua(i,1,1,1) = bm1(i,1,1,1)*u(i,1,1,1)/ua(i,1,1,1)
      enddo

      call gs_op(gs_avg_hndl,ua,1,1,0) ! Sum weighted values


      return
      end
c-----------------------------------------------------------------------
      subroutine set_gs_zavg_hndl(gs_avg_hndl,nelxy,ifld)

c     Set the z-average handle

      include 'SIZE'
      include 'TOTAL'

      integer gs_avg_hndl,e,ex,ey,ez,eg

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /c_is1/ glo_num(lx1,ly1,lz1,lelv)
      integer*8 glo_num,ex_g


      nel = nelfld(ifld)
      do e=1,nel
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelxy,1,1)

         ex_g = ex       ! Ensure int*8 promotion
         do k=1,nz1      ! Enumerate points in the x-y plane
            do j=1,ny1
            do i=1,nx1
               glo_num(i,j,k,e) = i+nx1*(j-1) + nx1*ny1*(ex_g-1)
            enddo
            enddo
         enddo

      enddo

      n = nel*nx1*ny1*nz1

      call gs_setup(gs_avg_hndl,glo_num,n,nekcomm,mp)

      return
      end
c-----------------------------------------------------------------------
      subroutine xx_avg(ua,u,gs_avg_hndl,nelxx,nelyy,nelzz,ifld)
      include 'SIZE'
      include 'TOTAL'

      real u (lx1,ly1,lz1,lelt)
      real ua(lx1,ly1,lz1,lelt)

      integer gs_avg_hndl,e,ex,ey,ez,eg,nelxx,nelyy,nelzz

      if (gs_avg_hndl.eq.0) then
          call set_gs_xxavg_hndl(gs_avg_hndl,nelxx,nelyy*nelzz,ifld)
      endif

      nel = nelfld(ifld)
      n   = nx1*ny1*nz1*nel

      call copy(ua,bm1,n)              ! Set the averaging weights
      call gs_op(gs_avg_hndl,ua,1,1,0) ! Sum weights over columns


      do i=1,n                          ! ua = (w_j*u_j)/( sum_i w_i)
         ua(i,1,1,1) = bm1(i,1,1,1)*u(i,1,1,1)/ua(i,1,1,1)
      enddo

      call gs_op(gs_avg_hndl,ua,1,1,0) ! Sum weighted values


      return
      end
c-----------------------------------------------------------------------
      subroutine set_gs_xxavg_hndl(gs_avg_hndl,nelxx,nelyz,ifld)

c     Set the z-average handle

      include 'SIZE'
      include 'TOTAL'

      integer gs_avg_hndl,e,ex,ey,ez,eg,nelyz,nelxx

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /c_is1/ glo_num(lx1,ly1,lz1,lelv)
      integer*8 glo_num,ex_g


      nel = nelfld(ifld)
      do e=1,nel
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelxx,nelyz,1)

         ex_g = ey       ! Ensure int*8 promotion
         do k=1,nz1      ! Enumerate points in the x-y plane
            do j=1,ny1
            do i=1,nx1
               glo_num(i,j,k,e) = j+ny1*(k-1) + ny1*nz1*(ex_g-1)
            enddo
            enddo
         enddo

      enddo

      n = nel*nx1*ny1*nz1

      call gs_setup(gs_avg_hndl,glo_num,n,nekcomm,mp)

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_v(uvw,xyz,n)
c
c     evaluate velocity for list of points xyz
c
c     Note:  -- modify
c     intpts to get rid off ' WARNING: point on boundary or ...'

      include 'SIZE'
      include 'TOTAL'

      real uvw(ldim,n),xyz(ldim,n)
      logical ifjac,ifpts

      parameter(nmax=lpart)
      common /rv_intp/ pts(ldim*nmax)
      common /iv_intp/ ihandle
      common /outtmp/ wrk(lx1*ly1*lz1*lelt,3)

      integer icalld,e
      save    icalld
      data    icalld /0/

      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt

      if (n.gt.nmax) call exitti ('ABORT: interp_v() n > nmax!$',n)

      if (nelgt.ne.nelgv) call exitti
     $   ('ABORT: interp_v() nelgt.ne.nelgv not yet supported!$',nelgv)

      do i=1,n                          ! ? not moving -> save?
         pts(i)     = xyz(1,i)
         pts(i + n) = xyz(2,i)
         if (if3d)  pts(i + n*2) = xyz(3,i)
      enddo

      if (icalld.eq.0) then             ! interpolation setup   !? intpts_done(ih_intp_v)?
        icalld = 1
        tolin  = 1.e-8
        call intpts_setup(tolin,ihandle)
      endif

      nflds  = ndim ! number of fields to interpolate

      ! pack working array
      call opcopy(wrk(1,1),wrk(1,2),wrk(1,3),vx,vy,vz)

      ! interpolate
      ifjac  = .true.           ! output transpose (of Jacobian)
      ifpts  = .true.            ! find points
      call intpts(wrk,nflds,pts,n,uvw,ifjac,ifpts,ihandle)      ! copy array instead?

      return
      end
c-----------------------------------------------------------------------
