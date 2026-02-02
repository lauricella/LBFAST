#include "defines.h"
#if !defined(_OPENACC)  
#error "To use this module the macros _OPENACC should be defined."
#endif

module lb_cuda_repulsive

   use vars
   use iso_c_binding
   use cudafor
   use mpi_template, only: coords,dostop,doerror,mydev,myrank,nprocs,nbuff,nbuffbvec
   use lb_cuda_vars

   implicit none
   

contains
#ifdef TWOCOMPONENT   
   attributes(global) subroutine thinfilm_scan_mark_kernel(flop,nx,ny,nz,coords,q_th,win,cosOppT,pwr,A_rep,isfluid, &
    rep_mask,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)
      implicit none
      
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      real(kind=db) :: q_th,win,cosOppT,pwr,A_rep	
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      integer(kind=isf), dimension(1:nx,1:ny,1:nz) :: rep_mask
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      integer :: i,j,k,gi,gj,gk,myblock,ii,jj,kk,iii,jjj,kkk
      
      integer :: di,dj,dk
	  integer :: diii,djjj,dkkk
	  real(kind=db) :: nix,niy,niz, dotn, qloc, qneig, face
	  real(kind=db) :: best_r2, r2, best_face
	  integer :: iii_best, jjj_best, kkk_best
	  logical :: found
	  real(kind=db), parameter :: eps = 1.0e-12_db
	  integer :: oii,ojj,okk
      integer :: oxblock,oyblock,ozblock,omyblock
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z
      
      rep_mask(i,j,k) = 0
   
      locauxfields_s(ii,jj,kk,3,myblock) = ZEROSTR
      locauxfields_s(ii,jj,kk,4,myblock) = ZEROSTR
      locauxfields_s(ii,jj,kk,5,myblock) = ZEROSTR

	  ! gate: interfacial cell (use clamped phi for q)
      qloc = real(phifields_s(ii,jj,kk,1,myblock),kind=db)

      qloc = min(max(qloc,0.0_db),1.0_db)
      qloc = qloc*(1.0_db - qloc)
      if (qloc < q_th) return
      nix = real(auxfields_s(ii,jj,kk,1,myblock),kind=db)
      niy = real(auxfields_s(ii,jj,kk,2,myblock),kind=db)
      niz = real(auxfields_s(ii,jj,kk,3,myblock),kind=db)

      best_r2   = HUGE(1.0_db)
      best_face = -1.0_db
      found     = .false.

      do di = -win, win
        do dj = -win, win
          do dk = -win, win
                  if (di==0 .and. dj==0 .and. dk==0) cycle

				  ! ---- 
				  iii = i + di
				  jjj = j + dj
				  kkk = k + dk
				  
				  if(abs(isfluid(iii,jjj,kkk)) .ne. 1)cycle

				  ! ---- minimum-image index differences
				  diii = iii - i
				 
				  djjj = jjj - j
				  
				  dkkk = kkk - k
				  
                  
				  r2 = real(diii,db)*real(diii,db) + real(djjj,db)*real(djjj,db) + real(dkkk,db)*real(dkkk,db)
				  if (r2 < eps) cycle

				  ! ---- neighbor interfacial gate (clamped) + similarity
				  
				  oxblock=(iii+2*TILE_DIMx-1)/TILE_DIMx   
                  oyblock=(jjj+2*TILE_DIMy-1)/TILE_DIMy     
                  ozblock=(kkk+2*TILE_DIMz-1)/TILE_DIMz 
                  omyblock=(oxblock-1)+(oyblock-1)*nxblock_d+(ozblock-1)*nxyblock_d+1
                  oii=iii-oxblock*TILE_DIMx+2*TILE_DIMx
                  ojj=jjj-oyblock*TILE_DIMy+2*TILE_DIMy
                  okk=kkk-ozblock*TILE_DIMz+2*TILE_DIMz
                  
				  qneig = real(phifields_s(oii,ojj,okk,1,omyblock),kind=db)
				  qneig = min(max(qneig,0.0_db),1.0_db)
				  qneig = qneig*(1.0_db - qneig)
				  if ( (qneig < q_th) .or. (abs(qneig - qloc) > 0.1_db*max(qloc,1.0e-12_db)) ) cycle

				  ! ---- facing condition (opposite normals): dotn <= cosOppT
				  dotn = nix*real(auxfields_s(oii,ojj,okk,1,omyblock),kind=db) &
				   + niy*real(auxfields_s(oii,ojj,okk,2,omyblock),kind=db) &
				   + niz*real(auxfields_s(oii,ojj,okk,3,omyblock),kind=db) 

				  if (dotn > cosOppT) cycle
				  face = 0.5_db*(1.0_db - dotn)   ! in [0,1]

				  ! ---- pick nearest; tie-break by larger 'face'
				  if (r2 < best_r2 - 1.0e-14_db) then
					best_r2 = r2; best_face = face
					iii_best = iii; jjj_best = jjj; kkk_best = kkk
					found   = .true.
				  else if (abs(r2 - best_r2) <= 1.0e-14_db) then
					if (face > best_face) then
					  best_face = face
					  iii_best = iii; jjj_best = jjj; kkk_best = kkk
					  found   = .true.
					end if
				  end if

				end do
			  end do
			end do

			if (found) then
              locauxfields_s(ii,jj,kk,3,myblock) = real(iii_best,kind=strdb)
			  locauxfields_s(ii,jj,kk,4,myblock) = real(jjj_best,kind=strdb)
			  locauxfields_s(ii,jj,kk,5,myblock) = real(kkk_best,kind=strdb)
			  
			  rep_mask(i,j,k) = 1
			end if
   
      return
      
   endsubroutine thinfilm_scan_mark_kernel
   
   attributes(global) subroutine repulsive_flux_normal_kernel(flop,nx,ny,nz,coords,width,q_th,win,cosOppT,pwr,A_rep,isfluid, &
    rep_mask,ntotphifields,ntotauxfields,ntotlocauxfields,phifields_s,auxfields_s,locauxfields_s)
      implicit none
      
      integer :: flop,nx,ny,nz
      integer, dimension(3) :: coords
      real(kind=db) :: width,q_th,win,cosOppT,pwr,A_rep	
      integer(kind=isf), dimension(1-nbuff:nx+nbuff,1-nbuff:ny+nbuff,1-nbuff:nz+nbuff) :: isfluid
      integer(kind=isf), dimension(1:nx,1:ny,1:nz) :: rep_mask
      integer :: ntotphifields,ntotauxfields,ntotlocauxfields
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nphifields,nblocks_d) :: phifields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nauxfields,nblocks_d) :: auxfields_s
      real(kind=strdb), dimension(TILE_DIMx,TILE_DIMy,TILE_DIMz,nlocauxfields,nblocks_d) :: locauxfields_s
      
      integer :: i,j,k,gi,gj,gk,myblock,ii,jj,kk,iii,jjj,kkk
      
      real(kind=db) :: q1,q2,qpair,qcl,loc_phi,loc_phi2
	  real(kind=db) :: nx1,ny1,nz1, nx2,ny2,nz2
	  real(kind=db) :: dx,dy,dz, r, rinv, face, arg_arcosh, ach, Wfilm, wdth
	  real(kind=db) :: nsx,nsy,nsz, nsmag,alpha,cap,scales
	  real(kind=db), parameter :: eps = 1.0e-9_db
      
	  integer :: oii,ojj,okk
      integer :: oxblock,oyblock,ozblock,omyblock
      
      i = (blockIdx%x-1) * TILE_DIMx + threadIdx%x
      j = (blockIdx%y-1) * TILE_DIMy + threadIdx%y
      k = (blockIdx%z-1) * TILE_DIMz + threadIdx%z
      
      if(abs(isfluid(i,j,k)) .ne. 1)return
      
      gi=nx*coords(1)+i
      gj=ny*coords(2)+j
      gk=nz*coords(3)+k
      
      myblock=blockIdx%x+blockIdx%y*nxblock_d+blockIdx%z*nxyblock_d+1
      
      ii=threadIdx%x
      jj=threadIdx%y
      kk=threadIdx%z

	  locauxfields_s(ii,jj,kk,6,myblock)=ZEROSTR
	  locauxfields_s(ii,jj,kk,7,myblock)=ZEROSTR
	  locauxfields_s(ii,jj,kk,8,myblock)=ZEROSTR
	
	  if (rep_mask(i,j,k) .ne. 1) return

      loc_phi=real(phifields_s(ii,jj,kk,1,myblock),kind=db)

	  q1 = loc_phi*(1.0_db - loc_phi)
	
	  if (q1 <= eps) return

	  iii = int(locauxfields_s(ii,jj,kk,3,myblock))
	  jjj = int(locauxfields_s(ii,jj,kk,4,myblock))
	  kkk = int(locauxfields_s(ii,jj,kk,5,myblock))

	  !line-of-centers
	  dx = real(iii - i,db)
	  dy = real(jjj - j,db)
	  dz = real(kkk - k,db)
	  r  = sqrt(dx*dx + dy*dy + dz*dz)
	  if (r <= eps) return
	  rinv = 1.0_db / r
	  dx = dx*rinv; dy = dy*rinv; dz = dz*rinv      ! u

	  !normals
	  nx1 = real(auxfields_s(ii,jj,kk,1,myblock),kind=db)
      ny1 = real(auxfields_s(ii,jj,kk,2,myblock),kind=db)
      nz1 = real(auxfields_s(ii,jj,kk,3,myblock),kind=db)
	  
	  oxblock=(iii+2*TILE_DIMx-1)/TILE_DIMx   
      oyblock=(jjj+2*TILE_DIMy-1)/TILE_DIMy     
      ozblock=(kkk+2*TILE_DIMz-1)/TILE_DIMz 
      omyblock=(oxblock-1)+(oyblock-1)*nxblock_d+(ozblock-1)*nxyblock_d+1
      oii=iii-oxblock*TILE_DIMx+2*TILE_DIMx
      ojj=jjj-oyblock*TILE_DIMy+2*TILE_DIMy
      okk=kkk-ozblock*TILE_DIMz+2*TILE_DIMz
      
	  nx2 = real(auxfields_s(oii,ojj,okk,1,omyblock),kind=db)
	  ny2 = real(auxfields_s(oii,ojj,okk,2,omyblock),kind=db)
	  nz2 = real(auxfields_s(oii,ojj,okk,3,omyblock),kind=db)

	  !facing factor in [0,1]
	  face = max( 0.0_db, -(nx1*nx2 + ny1*ny2 + nz1*nz2) )

	  if (face <= eps) return

	  !symmetric normal: bisector n1 - n2 (for facing sheets)
	  nsx = nx1 - nx2
	  nsy = ny1 - ny2
	  nsz = nz1 - nz2
	  nsmag = sqrt(nsx*nsx + nsy*nsy + nsz*nsz)
	  if (nsmag <= eps) then
	    !fallback to line-of-centers
	    nsx = dx; nsy = dy; nsz = dz
	  else
	    nsx = nsx / nsmag
	    nsy = nsy / nsmag
	    nsz = nsz / nsmag
	  end if

	  !orient so u·nsym >= 0  (partner will flip)
	  if (dx*nsx + dy*nsy + dz*nsz < 0.0_db) then
	    nsx = -nsx; nsy = -nsy; nsz = -nsz
	  end if

	  !symmetric magnitude from qpair
	  loc_phi2=real(phifields_s(oii,ojj,okk,1,omyblock),kind=db)

	  q2 = loc_phi2*(1.0_db - loc_phi2)
	  qpair = 0.5_db*(q1 + q2)
	  qcl   = min( max(qpair, eps), 0.25_db - eps )

	  arg_arcosh = 1.0_db / ( 2.0_db*sqrt(qcl) )
	  if (arg_arcosh <= 1.0_db) return

	  ach   = log( arg_arcosh + sqrt(arg_arcosh*arg_arcosh - 1.0_db) )
	  Wfilm = width * ach
	  if (Wfilm <= 0.0_db) return

	  wdth  = 1.0_db /(1.0 + wfilm**4.0) ! ( 1.0_db + (1.0_db/Wfilm)**4 )

	  !final purely-normal, symmetric repulsive flux
	  dx = A_rep * wdth * qcl * face * nsx
	  dy = A_rep * wdth * qcl * face * nsy
	  dz = A_rep * wdth * qcl * face * nsz
	
	  alpha = 1.5_db
	  cap   = alpha * (abs(dx)+abs(dy)+abs(dz)) 
	  scales = min(1.0_db, loc_phi / max(cap, 1.0e-9_db))

	  locauxfields_s(ii,jj,kk,6,myblock) = real(dx * scales,kind=strdb)
	  locauxfields_s(ii,jj,kk,7,myblock) = real(dy * scales,kind=strdb)
	  locauxfields_s(ii,jj,kk,8,myblock) = real(dz * scales,kind=strdb)
     
 end subroutine repulsive_flux_normal_kernel
#endif

endmodule lb_cuda_repulsive
