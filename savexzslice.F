c234567
      subroutine makexzslice(fftwmesh, nx,ny,nz,local_nz,slice)
	  integer nx,ny,nz,local_nz,myrank,nrank,stepcount
	  real fftwmesh(2*(nx/2+1),ny,local_nz)
	  real slice(nx,local_nz)
	  integer i,j,k

	  do k = 1, local_nz
	  do i = 1, nx
	     slice(i,k) = fftwmesh(i,1,k)
	  enddo
	  enddo
	  return
	  end
