!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!The amrvac_lib in a MATLAB class designed to load the MPI-AMRVAC
! xxxx.dat files. The library contains two user accessable functions: 
! get_data() - used to read out the data from file, 
! get_mesh() - udes to create the coordinat variables.
! 
! For specifications and usage of the functions consult the descriptors
! of the functions.
!
! version: 1 
! created by: Luka Poniatowski 2021 for matlab
! adapted by: Levin Hennicker 2021 for fortran
!
! Changes log:
! 04/05/21 - The output data comes now as a structure with keys
!            corresponding to the name s in heared.w_name.
!            To access the data instead of data(DIM,variableindex)
!            NOW data.varname(DIM) or data.('varname')(DIM)
!
module mod_amrvac_reader
!
use prog_type
use hdf5
use fund_const
use mod_sort
use mod_grid
use mod_interp1d
use mod_interp3d
use omp_lib
!
implicit none
!
!------------------------------datatypes--------------------------------
!
type header

   integer(i4b) :: version, treeoffset, blockoffset, nw, ndir, ndim, levmax, nleafs, nparents, it, levmax_usr
   integer(i4b), dimension(:), allocatable :: domain_nx, block_nx
   real(dp) :: global_time
   real(dp), dimension(:), allocatable :: xprobmin, xprobmax !are face values

   integer(i4b), dimension(:), allocatable :: iperiodic
   logical, dimension(:), allocatable :: periodic
   logical :: staggered
   character(len=16) :: geometry, physics_type
   character(len=16), dimension(:), allocatable :: w_names
    
!The physics parameters, such as gamma
   integer(i4b) :: n_params
   real(dp), dimension(:), allocatable :: parameters
   character(len=16), dimension(:), allocatable :: parameter_names

!Indexes for file output (for restarting)
   integer(i4b) :: snapshotnext, slicenext, collapsenext

   character(len=500) :: fname
!
end type header
!
!
!
type tree
!
   logical, dimension(:), allocatable :: leaf
   integer(i4b), dimension(:), allocatable :: ileaf, refinement_level
   integer(i4b), dimension(:), allocatable :: spatial_index1d   
   integer(i4b), dimension(:,:), allocatable :: spatial_index
   integer(i8b), dimension(:), allocatable :: offset_block

end type tree
!
!
!
type block
!
   integer(i4b) :: level
   integer(i4b), dimension(4) :: block_shape
   integer(i4b), dimension(:), allocatable :: indexes, ncells
   integer(i4b) :: nx, ny, nz
   real(dp), dimension(:), allocatable :: fblock1d
   real(dp), dimension(:,:,:,:), allocatable :: fblock
   real(dp), dimension(:), allocatable :: x_center, x_edgel, x_edger
   real(dp), dimension(:), allocatable :: x_coord, y_coord, z_coord
   !
   !level: level of the block
   !ncells: number of cells for all dimensions if complete domain had this refinement level
   !x_xenter: cell center coordinates of the block
   !x_edgel:  left face coordinates of the block
   !x_edger:  right face coordinates of the block
   !nxyz: number of cels within the block for each dimension
   !xyzcoord: x, y, z coordinates fo all cells within the block
   !indexes: index of the block for all dimensions
   !
!
end type block
!
!
!
type data
!
   integer(i4b) :: ndim, nw
   integer(i4b), dimension(4) :: data_shape   
   real(dp), dimension(:,:,:,:), allocatable :: data
   character(len=16), dimension(:), allocatable :: w_names
   type(grid) :: mesh
!
end type data
!
!
!
type alldata
   type(header) :: header
   type(tree) :: tree
   type(block), dimension(:), allocatable :: blocks
   type(data) :: data
end type alldata
!
!------------------------------subroutines------------------------------
!
contains  
!
  function get_data(fname, levmax_usr, stretching, interp_method, nd, grid_out)
    !
    !get all the data of a file with fname
    !
    !optional arguments:
    !levmax_usr:   maximum refinement level to read in each block
    !interp_method: interpolation method within the block when being refined (either constant cell-center values or linear interpolation)
    !stretching:   stretching factors of the grid (not implemented yet)
    !nd:           dimensions of the output grid
    !grid_out:     a specified grid onto which things shall be calculated
    !
    ! ... arguments
    type(alldata) :: get_data    
    character(len=*), intent(in) :: fname
    character(len=*), intent(in), optional :: interp_method
    integer(i4b), intent(in), optional :: levmax_usr
    real(dp), dimension(:), allocatable, intent(in), optional :: stretching
    integer(i4b), dimension(:), allocatable, intent(in), optional :: nd
    type(grid), intent(in), optional :: grid_out

    ! ... local scalars
    integer(i4b) :: err, alevmax_usr
    real(dp) :: tstart, tend
    !
    ! ... local characters
    character(len=10) :: ainterp_method

    ! ... local arrays
    integer(i4b), dimension(:), allocatable :: nx
    !
    ! ... local derived types
    type(header) :: fheader
    type(tree) :: ftree
    type(block), dimension(:), allocatable :: fblocks
    type(data) :: fdata
    type(grid) :: agrid_out, grid_in
    !
    !
    !
    if(.not.present(levmax_usr)) then
       alevmax_usr = 1000
    else
       alevmax_usr = levmax_usr
    endif
    !
    if(.not.present(interp_method)) then
       ainterp_method = 'lin'
    else
       ainterp_method = interp_method
    endif
    !    
    !
    fheader = get_header(fname, levmax_usr=alevmax_usr)
    ftree = get_tree(fheader)
    fblocks = get_blocks(fheader,ftree)
    grid_in = grid_internal(fheader, fblocks)

    !default resolution for the regridding: everything on maximum AMR level
    allocate(nx(fheader%ndim), stat=err)    
    if(.not.present(nd)) then
       nx = 2.**(fheader%levmax_usr-1)*fheader%domain_nx
    else
       nx = nd
    endif
    !
    if(.not.present(grid_out)) then
       !calculate the grid onto which things are remeshed       
       agrid_out = grid_external(grid_in, nx)
    else
       !take the specified grid to interpolate values
       agrid_out = grid_out
    endif

    !calculate the reconstruction of the data
    tstart = omp_get_wtime()
!    fdata = reconstruct2(fheader,fblocks,agrid_out)  !reconstruct onto grid_out mesh
    fdata = reconstruct(fheader,ftree,fblocks,interp_method=ainterp_method,grid_out=agrid_out)
    tend = omp_get_wtime()

    write(*,*) 'timing for reconstruction', tend-tstart
    write(*,*)

    get_data%header = fheader
    get_data%tree = ftree
    get_data%blocks = fblocks
    get_data%data = fdata
    !
  end function get_data
!
!-----------------------------------------------------------------------
!
  subroutine print_data(my_alldata)
    !
    ! ... arguments
    type(alldata) :: my_alldata

    call print_header(my_alldata%header)
    call print_tree(my_alldata%tree)
    call print_blocks(my_alldata%blocks, my_alldata%header)
    call print_reconstruct(my_alldata%data)
    !
  end subroutine print_data
!
!-----------------------------------------------------------------------
!
  subroutine fexist(fname)
    !
    !... arguments
    character(len=*), intent(in) :: fname
    !
    ! ... local logicals
    logical :: lcheck
    !
    inquire(file=trim(fname), exist=lcheck)
    if(.not.lcheck) then
       write(*,*) 'error in fexist: file "', trim(fname), '" does not exist'
       stop
    endif
    !
  end subroutine fexist
!
!-----------------------------------------------------------------------
!
  function get_header(fname, levmax_usr)
    !
    ! ... arguments
    character(len=*), intent(in) :: fname
    type(header) :: get_header
    integer(i4b), intent(in), optional :: levmax_usr
    !
    ! ... local scalars
    integer(i4b) :: err
    !
    call fexist(fname)
    !
    get_header%fname = trim(fname)
    !
    open(1, file=trim(fname), form='unformatted', access='stream')
       rewind(1)
       read(1) get_header%version
       read(1) get_header%treeoffset
       read(1) get_header%blockoffset
       read(1) get_header%nw
       read(1) get_header%ndir
       read(1) get_header%ndim
       read(1) get_header%levmax
       read(1) get_header%nleafs
       read(1) get_header%nparents
       read(1) get_header%it
       read(1) get_header%global_time
       allocate(get_header%xprobmin(get_header%ndim), stat=err)
       allocate(get_header%xprobmax(get_header%ndim), stat=err)
       allocate(get_header%domain_nx(get_header%ndim), stat=err)
       allocate(get_header%block_nx(get_header%ndim), stat=err)
       allocate(get_header%iperiodic(get_header%ndim), stat=err)       
       allocate(get_header%periodic(get_header%ndim), stat=err)
       read(1) get_header%xprobmin
       read(1) get_header%xprobmax
       read(1) get_header%domain_nx
       read(1) get_header%block_nx
       read(1) get_header%iperiodic
       get_header%periodic = get_header%iperiodic
       read(1) get_header%geometry
       read(1) get_header%staggered
       allocate(get_header%w_names(get_header%nw), stat=err)
       read(1) get_header%w_names
       read(1) get_header%physics_type
       read(1) get_header%n_params
       allocate(get_header%parameters(get_header%n_params), stat=err)
       allocate(get_header%parameter_names(get_header%n_params), stat=err)
       read(1) get_header%parameters
       read(1) get_header%parameter_names
       read(1) get_header%snapshotnext
       read(1) get_header%slicenext
       read(1) get_header%collapsenext
    !
    close(1)
    !
    if(present(levmax_usr)) then
       if(levmax_usr.le.0) then
          get_header%levmax_usr = 1    !lower refinement level not allowed
       elseif(levmax_usr.gt.get_header%levmax) then
          get_header%levmax_usr = get_header%levmax   !larger refinement level not allowed
       else
          get_header%levmax_usr = levmax_usr   !default
       endif
    else
       get_header%levmax_usr = get_header%levmax
    endif
    !
    !
  end function get_header
!
!-----------------------------------------------------------------------
!
  subroutine print_header(fheader)
    !
    ! ... arguments
    type(header), intent(in) :: fheader
    !
    ! ... local scalars
    !
    ! ... local logicals
    !
    write(*,*) '-------------HEADER-----------'
    write(*,*)
    write(*,*) 'fname', trim(fheader%fname)
    write(*,*) 'version', fheader%version
    write(*,*) 'treeoffset', fheader%treeoffset
    write(*,*) 'blockoffset', fheader%blockoffset
    write(*,*) 'nw', fheader%nw
    write(*,*) 'ndir', fheader%ndir
    write(*,*) 'ndim', fheader%ndim
    write(*,*) 'levmax', fheader%levmax
    write(*,*) 'nleafs', fheader%nleafs
    write(*,*) 'nparents', fheader%nparents
    write(*,*) 'it', fheader%it
    write(*,*) 'global_time', fheader%global_time
    write(*,*) 'xprobmin', fheader%xprobmin
    write(*,*) 'xprobmax', fheader%xprobmax
    write(*,*) 'domain_nx', fheader%domain_nx
    write(*,*) 'block_nx', fheader%block_nx
    write(*,*) 'periodic', fheader%periodic
    write(*,*) 'geometry', fheader%geometry
    write(*,*) 'staggered', fheader%staggered
    write(*,*) 'w_names', fheader%w_names
    write(*,*) 'physics_type', fheader%physics_type
    write(*,*) 'n_params', fheader%n_params
    write(*,*) 'parameters', fheader%parameters
    write(*,*) 'parameter_names', fheader%parameter_names
    write(*,*) 'snapshotnext', fheader%snapshotnext
    write(*,*) 'slicenext', fheader%slicenext
    write(*,*) 'collapsenext', fheader%collapsenext
    write(*,*)
    !
  end subroutine print_header
!
!-----------------------------------------------------------------------
!
  function get_tree(fheader)
    !
    ! ... arguments
    type(header), intent(in) :: fheader
    type(tree) :: get_tree
    !
    ! ... local scalars
    integer(i4b) :: err
    !
    call fexist(fheader%fname)
    !
    open(1, file=trim(fheader%fname), form='unformatted', access='stream')
    !
       call fseek(1,fheader%treeoffset, 0) !navigate to the begining of the tree info
       !
       ! secuentially read the data for the tree
       allocate(get_tree%ileaf(fheader%nleafs+fheader%nparents), stat=err)
       read(1) get_tree%ileaf
       get_tree%leaf=get_tree%ileaf
       allocate(get_tree%refinement_level(fheader%nleafs), stat=err)
       read(1) get_tree%refinement_level
       allocate(get_tree%spatial_index1d(fheader%ndim*fheader%nleafs), stat=err)
       allocate(get_tree%spatial_index(fheader%ndim,fheader%nleafs), stat=err)
       read(1) get_tree%spatial_index1d
       get_tree%spatial_index = reshape(get_tree%spatial_index1d, (/fheader%ndim, fheader%nleafs /))
       allocate(get_tree%offset_block(fheader%nleafs),stat=err)
       read(1) get_tree%offset_block
       get_tree%offset_block = get_tree%offset_block + (2*fheader%ndim*4)
    !                
    close(1)
    !   
  end function get_tree
!
!-----------------------------------------------------------------------
!
  subroutine print_tree(ftree)
    !
    ! ... arguments
    type(tree), intent(in) :: ftree
    !
    ! ... local scalars
    !
    ! ... local logicals
    !
    write(*,*) '--------------TREE------------'
    write(*,*)
    write(*,*) 'leaf:'
    write(*,*) ftree%leaf
    write(*,*) 'refinement_level:'
    write(*,*) ftree%refinement_level
    write(*,*) 'spatial_index:'       
    write(*,*) ftree%spatial_index
    write(*,*) 'offset_block:'       
    write(*,*) ftree%offset_block
    write(*,*)
    !
  end subroutine print_tree
!
!-----------------------------------------------------------------------
!
  function get_blocks(fheader,ftree)
    !
    ! ... arguments
    type(header), intent(in) :: fheader
    type(tree), intent(in) :: ftree
    type(block), dimension(:), allocatable :: get_blocks
    !
    ! ... local scalars
    integer(i4b) :: i, ind1, n_blocks, err
    !
    !... local arrays
    real(dp), dimension(:), allocatable :: xmin, xmax, deltax
    !
    ! ... local derived types
    type(block) :: fblock
    !
    call fexist(fheader%fname)
    !
    select case(fheader%ndim)
       case(1)
          fblock%nx = fheader%block_nx(1)
          fblock%ny = 1
          fblock%nz = 1
          allocate(fblock%x_coord(fblock%nx), stat=err)
          allocate(fblock%y_coord(fblock%ny), stat=err)
          allocate(fblock%z_coord(fblock%nz), stat=err)
       case(2)
          fblock%nx = fheader%block_nx(1)
          fblock%ny = fheader%block_nx(2)
          fblock%nz = 1
          allocate(fblock%x_coord(fblock%nx), stat=err)
          allocate(fblock%y_coord(fblock%ny), stat=err)
          allocate(fblock%z_coord(fblock%nz), stat=err)
       case(3)
          fblock%nx = fheader%block_nx(1)
          fblock%ny = fheader%block_nx(2)
          fblock%nz = fheader%block_nx(3)
          allocate(fblock%x_coord(fblock%nx), stat=err)
          allocate(fblock%y_coord(fblock%ny), stat=err)
          allocate(fblock%z_coord(fblock%nz), stat=err)
       case default
          stop 'error in get_blocks: fheader%ndim not properly set'
    end select
    !
    !define the length of the block-index array and the number of cells in each dimension,
    !and the center and edge coordinates of the block for each dimension
    allocate(fblock%indexes(fheader%ndim), stat=err)
    allocate(fblock%ncells(fheader%ndim), stat=err)
    allocate(fblock%x_center(fheader%ndim), stat=err)
    allocate(fblock%x_edgel(fheader%ndim), stat=err)    
    allocate(fblock%x_edger(fheader%ndim), stat=err)
    allocate(xmin(fheader%ndim), stat=err)
    allocate(xmax(fheader%ndim), stat=err)
    allocate(deltax(fheader%ndim), stat=err)
    xmin=fheader%xprobmin
    xmax=fheader%xprobmax
    !
    !define the block shape   
    fblock%block_shape = (/ fblock%nx, fblock%ny, fblock%nz, fheader%nw /)
    !
    allocate(fblock%fblock1d(fblock%block_shape(1)*fblock%block_shape(2)*fblock%block_shape(3)*fblock%block_shape(4)), stat=err)
    if(err.ne.0) stop 'error in get_blocks: allocation'
    allocate(fblock%fblock(fblock%block_shape(1), fblock%block_shape(2), fblock%block_shape(3), fblock%block_shape(4)), stat=err)
    if(err.ne.0) stop 'error in get_blocks: allocation'
    !
    !get the number of blocks
    n_blocks = fheader%nleafs
    !
    !create array to output blocks
    allocate(get_blocks(n_blocks), stat=err)
    !
    !
    open(1, file=trim(fheader%fname), form='unformatted', access='stream')
       !
       !for neach block in the dat file
       do ind1 = 1, n_blocks
          !
          !indices of the block
          fblock%indexes = ftree%spatial_index(:,ind1)
          !
          !store the level of the block
          fblock%level = ftree%refinement_level(ind1)
          !
          !store the number of cells that would be required for complete domain if this domain was refined with this level
          fblock%ncells = fheader%domain_nx/fheader%block_nx * 2**(fblock%level-1)
          !
          !the delta increments for each cell if complete domain was refined with this level
          deltax = (xmax-xmin)/fblock%ncells
          !
          !the center and edge coordinates of the cell
          fblock%x_center = xmin + ((2*fblock%indexes-1)*deltax)/2.d0
          fblock%x_edgel = fblock%x_center-deltax/2.d0
          fblock%x_edger = fblock%x_center+deltax/2.d0
          !
          !set the coordinates in the complete block
          if(fheader%ndim.ge.1) then
             do i=1, fblock%nx
                fblock%x_coord(i) = fblock%x_edgel(1) + (2*i-1)/2.d0 * (fblock%x_edger(1)-fblock%x_edgel(1))/fblock%nx
             enddo
          endif
          if(fheader%ndim.ge.2) then
             do i=1, fblock%ny
                fblock%y_coord(i) = fblock%x_edgel(2) + (2*i-1)/2.d0 * (fblock%x_edger(2)-fblock%x_edgel(2))/fblock%ny
             enddo
          endif
          if(fheader%ndim.ge.3) then
             do i=1, fblock%nz
                fblock%z_coord(i) = fblock%x_edgel(3) + (2*i-1)/2.d0 * (fblock%x_edger(3)-fblock%x_edgel(3))/fblock%nz
             enddo
          endif
          !
          !navigate to the beginning byte of the block (offset from the beginning of the file + n_ghostcell_ints)
          call fseek(1,ftree%offset_block(ind1),0)
          !
          !read the number of bytes needed to fiil the block
          read(1) fblock%fblock1d !are cell center values
          fblock%fblock = reshape(fblock%fblock1d, fblock%block_shape)
          !
          !store the block in the output
          get_blocks(ind1) = fblock
       end do
       !
    close(1)
    !
  end function get_blocks
!
!-----------------------------------------------------------------------
!
  subroutine print_blocks(fblocks, fheader)
    !
    ! ... arguments
    type(header) :: fheader
    type(block), dimension(:), allocatable :: fblocks
    !
    ! ... local scalars
    integer(i4b) :: ind1, n_blocks
    !
    ! ... local arrays
    type(block) :: fblock
    !
    !get the number of blocks
    n_blocks = fheader%nleafs
    !
    write(*,*) '--------------BLOCKS----------'
    write(*,*)
    
    do ind1 = 1, n_blocks
       fblock=fblocks(ind1)
       write(*,*) 'i', ind1
       write(*,*) fblock%fblock(1,1,1,:)
    end do
    write(*,*)
    !
  end subroutine print_blocks
!
!-----------------------------------------------------------------------
!
  function reconstruct(fheader,ftree,fblocks,interp_method, grid_out)
    !
    !reconstruct all the data onto a certain level
    !if grid_out is specified, data will be interpolated onto that grid
    !
    !typically two times faster than function reconstruct2 at the expense of memory
    !
    !... arguments
    type(header) :: fheader
    type(tree) :: ftree
    type(block), dimension(:), allocatable :: fblocks
    type(data) :: reconstruct
    character(len=*), intent(in), optional :: interp_method
    type(grid), intent(in), optional :: grid_out
    !
    ! ... local scalars
    integer(i4b) :: ndim, max_level, max_level_usr, max_level_diff, nw, err, n_blocks, ind1, ind2, level
    integer(i4b) :: i, j, k, istart, iend, jstart, jend, kstart, kend
    integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
    integer(i4b) :: nx1_start, nx1_end, nx2_start, nx2_end, nx3_start, nx3_end
    integer(i4b) :: nx1_dat, nx2_dat, nx3_dat, nx1_increment, nx2_increment, nx3_increment
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
    real(dp) :: fnorm
    !
    ! ... local characters
    character(len=10) :: ainterp_method
    !
    ! ... local arrays
    integer(i4b), dimension(:), allocatable :: nb, nx, nx_lev1, indexes
    integer(i4b), dimension(:), allocatable :: ix1_neighbour, ix2_neighbour, iy1_neighbour, iy2_neighbour, iz1_neighbour, iz2_neighbour
    integer(i4b), dimension(4) :: nblock
    integer(i4b), dimension(:), allocatable :: ix1_usr, ix2_usr, ix3_usr
    real(dp), dimension(:,:,:,:), allocatable :: data_tmp
    real(dp), dimension(:), allocatable :: x_coord_tmp, y_coord_tmp, z_coord_tmp
    real(dp), dimension(:), allocatable :: fsum
    !
    ! ... local derived types
    type(block) :: fblock
    !
    if(.not.present(interp_method)) then
       ainterp_method='lin'
    else
       ainterp_method=interp_method
    endif
    !
    !get info from the header
    ndim = fheader%ndim            !number of dimensions
    max_level = fheader%levmax     !maximum refinement level
    max_level_usr = fheader%levmax_usr !maximum refinement level that should be used
    max_level_diff = max_level-max_level_usr !difference of data maximum level and user specified one
    nw = fheader%nw                !number of variables
    allocate(nx_lev1(ndim), stat=err)
    allocate(nx(ndim), stat=err)
    allocate(nb(ndim), stat=err)
    nx_lev1 = fheader%domain_nx        !number of grid points on 1 level
    nx = 2**(max_level_usr-1)*nx_lev1  !number of grid points on last level for max_level_usr
    nb = fheader%block_nx              !number of grid points in each dimension in block
    !
    !
    write(*,*) '------------performing reconstruct-------------'
    write(*,*)
    write(*,'(a50,3i6)') 'maximum resolution on finest level', 2**(max_level-1)*fheader%domain_nx
    write(*,'(a50,3i6)') 'maximum resolution on finest user level', 2**(max_level_usr-1)*fheader%domain_nx    
    write(*,'(a50,3i6)') 'applied resolution', grid_out%xgrid%nd, grid_out%ygrid%nd, grid_out%zgrid%nd
    write(*,*)    
    !
    allocate(indexes(ndim), stat=err)
    allocate(fsum(nw), stat=err)
    fsum=0.d0
    !
    !depending on number of dimentions create 2,3 or 4 dimensional output matrix
    select case(ndim)
       case(1)
          allocate(data_tmp(1,1,nx(1),nw), stat=err)
          allocate(x_coord_tmp(nx(1)), stat=err)
          allocate(y_coord_tmp(1), stat=err)
          allocate(z_coord_tmp(1), stat=err)
       case(2)
          allocate(data_tmp(1,nx(1),nx(2), nw), stat=err)
          allocate(x_coord_tmp(nx(1)), stat=err)
          allocate(y_coord_tmp(nx(2)), stat=err)          
          allocate(z_coord_tmp(1), stat=err)
       case(3)
          allocate(data_tmp(nx(1),nx(2), nx(3), nw), stat=err)
          allocate(x_coord_tmp(nx(1)), stat=err)
          allocate(y_coord_tmp(nx(2)), stat=err)
          allocate(z_coord_tmp(nx(3)), stat=err)
       case default
          stop 'error in reconstruct: ndim not properly set'
    end select
    !
    !get the total number of blocks
    n_blocks = fheader%nleafs
    !
    !for each block
    !$omp parallel &
    !$omp private(ind1, fblock, level, ind2, indexes, nblock, nx1_start, nx1_end, nx1_dat, nx2_start, nx2_end, nx2_dat, nx3_start, nx3_end, nx3_dat, nx1_increment, nx2_increment, nx3_increment, i, j, k, ii, jj, kk, ix1_usr, ix2_usr, ix3_usr, fsum, fnorm, istart, iend, jstart, jend, kstart, kend)
    !$omp do schedule(static, 8)
    do ind1 = 1, n_blocks
       !
       fblock = fblocks(ind1)  !temp storage of the block data
       level = ftree%refinement_level(ind1) !refinement level of current block
       !
       !if refinement level is less then max refine
       if(level.lt.max_level) then
          do ind2 = 1, max_level - level
             fblock = refined_block(ndim,fblock,interp_method=ainterp_method)
          enddo
       endif
       !
       !get the spacial location of the block
       indexes = ftree%spatial_index(:,ind1)
       !
       !get the number of data points in the current block
       !this is needed as after the refinement this number is different from Nb
       nblock = fblock%block_shape
       !
       !depending on number of dimensions copy the data from block into corresponding pars of the output data
       select case(ndim)
          case(1)
             !based on the spetial indes, refinment level of the block and number of points compute the starting
             !index in the output data matrix, where the block has to be copied to
             if(max_level_usr.ne.max_level) stop 'error in reconstruct: max_level_usr ne max_level not implemented yet (see 3d)'
             nx1_start = 1 + (nb(1) * 2**(max_level - level))*(indexes(1) - 1)
             nx1_end = nx1_start + nblock(3) - 1 !the end index is then the strat + number of points on refined block
             data_tmp(1,1,nx1_start:nx1_end,:) = fblock%fblock(1,1,:,:)
          case(2)
             if(max_level_usr.ne.max_level) stop 'error in reconstruct: max_level_usr ne max_level not implemented yet (see 3d)'
             nx1_start = 1 + (nb(1) * 2**(max_level - level))*(indexes(1) - 1)
             nx1_end = nx1_start + nblock(2) - 1
             !
             nx2_start = 1 + (nb(2) * 2**(max_level - level))*(indexes(2) - 1)
             nx2_end = nx2_start + nblock(3) - 1
             !
             data_tmp(1,nx1_start:nx1_end,nx2_start:nx2_end,:) = fblock%fblock(1,:,:,:)
          case(3)
             nx1_increment = 2**max_level_diff
             nx2_increment = 2**max_level_diff
             nx3_increment = 2**max_level_diff
             !
             nx1_start = 1 + (nb(1)/nx1_increment * 2**(max_level - level))*(indexes(1) - 1)
             nx1_end = nx1_start + nblock(1)/nx1_increment - 1
             nx1_dat = nx1_end-nx1_start + 1
             !
             nx2_start = 1 + (nb(2)/nx2_increment * 2**(max_level - level))*(indexes(2) - 1)
             nx2_end = nx2_start + nblock(2)/nx2_increment - 1
             nx2_dat = nx2_end-nx2_start + 1
             !
             nx3_start = 1 + (nb(3)/nx3_increment * 2**(max_level - level))*(indexes(3) - 1)
             nx3_end = nx3_start + nblock(3)/nx3_increment - 1
             nx3_dat = nx3_end-nx3_start + 1
             !
             if(allocated(ix1_usr)) deallocate(ix1_usr)
             allocate(ix1_usr(nx1_dat))
             ix1_usr(1) = nx1_start
             do i=2, nx1_dat
                ix1_usr(i) = ix1_usr(i-1)+1
             enddo
             !
             if(allocated(ix2_usr)) deallocate(ix2_usr)
             allocate(ix2_usr(nx2_dat))
             ix2_usr(1) = nx2_start
             do i=2, nx2_dat
                ix2_usr(i) = ix2_usr(i-1)+1
             enddo
             !
             if(allocated(ix3_usr)) deallocate(ix3_usr)
             allocate(ix3_usr(nx3_dat))
             ix3_usr(1) = nx3_start
             do i=2, nx3_dat
                ix3_usr(i) = ix3_usr(i-1)+1
             enddo
             !
             x_coord_tmp(ix1_usr) = fblock%x_coord
             y_coord_tmp(ix2_usr) = fblock%y_coord
             z_coord_tmp(ix3_usr) = fblock%z_coord
!             write(*,*) fblock%z_coord
!             write(*,*)
             do i=1, nx1_dat
                istart = 1 + (i-1)*nx1_increment
                iend = istart + nx1_increment - 1
                do j=1, nx2_dat
                   jstart = 1 + (j-1)*nx2_increment
                   jend = jstart + nx2_increment - 1                   
                   do k=1, nx3_dat
                      kstart = 1+(k-1)*nx3_increment
                      kend = kstart + nx3_increment - 1
                      fsum = 0.d0
                      fnorm = 0.d0
                      do ii=istart, iend
                         do jj=jstart, jend
                            do kk=kstart, kend
                               fsum = fsum + fblock%fblock(ii,jj,kk,:)
                               fnorm = fnorm + 1.d0
                            enddo
                         enddo
                      enddo
                      data_tmp(ix1_usr(i),ix2_usr(j),ix3_usr(k),:) = fsum/fnorm
                   enddo
                enddo
             enddo
!             if(ind1.eq.10) stop 'nasflkf'
             !
          case default
             stop 'error in refine_blocks: ndim not properly specified'
       end select
       !
    enddo
    !$omp enddo
    !$omp end parallel
    
    !
    !create the output data structure for each name in header.w_names create corresponding section in data structure
    allocate(reconstruct%w_names(nw), stat=err)    
    reconstruct%nw = nw
    reconstruct%ndim = ndim
    reconstruct%w_names = fheader%w_names
    !
    if(.not.present(grid_out)) then
       !
       select case(ndim)
          case(1)
             allocate(reconstruct%data(1,1,nx(1),nw), stat=err)
             allocate(reconstruct%mesh%xgrid%coord(nx(1)), stat=err)
             allocate(reconstruct%mesh%ygrid%coord(1), stat=err)
             allocate(reconstruct%mesh%zgrid%coord(1), stat=err)
             reconstruct%data_shape=(/1,1,nx(1),nw/)
             reconstruct%mesh%xgrid%nd=nx(1)
             reconstruct%mesh%ygrid%nd=1
             reconstruct%mesh%zgrid%nd=1
          case(2)
             allocate(reconstruct%data(1,nx(1),nx(2), nw), stat=err)
             allocate(reconstruct%mesh%xgrid%coord(nx(1)), stat=err)
             allocate(reconstruct%mesh%ygrid%coord(nx(2)), stat=err)
             allocate(reconstruct%mesh%zgrid%coord(1), stat=err)
             reconstruct%data_shape=(/1,1,nx(1),nw/)
             reconstruct%mesh%xgrid%nd=nx(1)
             reconstruct%mesh%ygrid%nd=nx(2)
             reconstruct%mesh%zgrid%nd=1
          case(3)
             allocate(reconstruct%data(nx(1),nx(2), nx(3), nw), stat=err)
             allocate(reconstruct%mesh%xgrid%coord(nx(1)), stat=err)
             allocate(reconstruct%mesh%ygrid%coord(nx(2)), stat=err)
             allocate(reconstruct%mesh%zgrid%coord(nx(3)), stat=err)
             reconstruct%data_shape=(/nx(1),nx(2),nx(3),nw/)
             reconstruct%mesh%xgrid%nd=nx(1)
             reconstruct%mesh%ygrid%nd=nx(2)
             reconstruct%mesh%zgrid%nd=nx(3)
          case default
             stop 'error in reconstruct: ndim not properly set'
       end select
       !
       reconstruct%mesh%xgrid%coord=x_coord_tmp
       reconstruct%mesh%ygrid%coord=y_coord_tmp
       reconstruct%mesh%zgrid%coord=z_coord_tmp
       do ind1 = 1, nw
          reconstruct%data(:,:,:,ind1) = data_tmp(:,:,:,ind1)
       enddo
       !
    else
       !perform interpolation onto the output grid
       select case(ndim)
          case(1)
             stop 'error in reconstruct: to be implemented'
          case(2)             
             stop 'error in reconstruct: to be implemented'
          case(3)
             reconstruct%data_shape=(/ grid_out%xgrid%nd, grid_out%ygrid%nd, grid_out%zgrid%nd, nw/)
             allocate(reconstruct%data(grid_out%xgrid%nd, grid_out%ygrid%nd, grid_out%zgrid%nd, nw), stat=err)
             allocate(reconstruct%mesh%xgrid%coord(grid_out%xgrid%nd), stat=err)
             allocate(reconstruct%mesh%ygrid%coord(grid_out%ygrid%nd), stat=err)
             allocate(reconstruct%mesh%zgrid%coord(grid_out%zgrid%nd), stat=err)             
             reconstruct%mesh = grid_out

             !find the neighbour indices
             allocate(ix1_neighbour(grid_out%xgrid%nd), stat=err)
             allocate(ix2_neighbour(grid_out%xgrid%nd), stat=err)
             allocate(iy1_neighbour(grid_out%ygrid%nd), stat=err)
             allocate(iy2_neighbour(grid_out%ygrid%nd), stat=err)
             allocate(iz1_neighbour(grid_out%zgrid%nd), stat=err)
             allocate(iz2_neighbour(grid_out%zgrid%nd), stat=err)
             
             do i=1, grid_out%xgrid%nd
                call find_index(grid_out%xgrid%coord(i),x_coord_tmp, nx(1), iim2, iim1, ii, iip1)
                ix1_neighbour(i) = iim1
                ix2_neighbour(i) = ii
             enddo
             do i=1, grid_out%ygrid%nd
                call find_index(grid_out%ygrid%coord(i),y_coord_tmp, nx(2), iim2, iim1, ii, iip1)
                iy1_neighbour(i) = iim1
                iy2_neighbour(i) = ii
             enddo
             do i=1, grid_out%zgrid%nd
                call find_index(grid_out%zgrid%coord(i),z_coord_tmp, nx(3), iim2, iim1, ii, iip1)
                iz1_neighbour(i) = iim1
                iz2_neighbour(i) = ii
             enddo

             !
             !$omp parallel &
             !$omp private(j, k, acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff, iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1)
             !$omp do schedule(static)
             do i=1, grid_out%xgrid%nd
                iim1 = ix1_neighbour(i)
                ii = ix2_neighbour(i)
                do j=1, grid_out%ygrid%nd
                   jjm1 = iy1_neighbour(j)
                   jj = iy2_neighbour(j)
                   do k=1, grid_out%zgrid%nd
                      kkm1 = iz1_neighbour(k)
                      kk = iz2_neighbour(k)
                      !inverse distance weighting used to avoid extrapolation at the edges
                      call coeff3d_8p_idw(x_coord_tmp(iim1), y_coord_tmp(jjm1), z_coord_tmp(kkm1), &
                                          x_coord_tmp(ii),   y_coord_tmp(jjm1), z_coord_tmp(kkm1), &
                                          x_coord_tmp(iim1), y_coord_tmp(jj),   z_coord_tmp(kkm1), &
                                          x_coord_tmp(ii),   y_coord_tmp(jj),   z_coord_tmp(kkm1), &
                                          x_coord_tmp(iim1), y_coord_tmp(jjm1), z_coord_tmp(kk), &
                                          x_coord_tmp(ii),   y_coord_tmp(jjm1), z_coord_tmp(kk), &
                                          x_coord_tmp(iim1), y_coord_tmp(jj),   z_coord_tmp(kk), &
                                          x_coord_tmp(ii),   y_coord_tmp(jj),   z_coord_tmp(kk), &
                                          grid_out%xgrid%coord(i), grid_out%ygrid%coord(j), grid_out%zgrid%coord(k), &
                                          acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff, fac=2.d0)!
!
                      reconstruct%data(i,j,k,:) = acoeff*data_tmp(iim1,jjm1,kkm1,:) + bcoeff*data_tmp(ii,jjm1,kkm1,:) + &
                                                  ccoeff*data_tmp(iim1,jj  ,kkm1,:) + dcoeff*data_tmp(ii,jj  ,kkm1,:) + &
                                                  ecoeff*data_tmp(iim1,jjm1,kk  ,:) + fcoeff*data_tmp(ii,jjm1,kk  ,:) + &
                                                  gcoeff*data_tmp(iim1,jj  ,kk  ,:) + hcoeff*data_tmp(ii,jj  ,kk  ,:)
                   enddo
                enddo
             enddo
             !$omp enddo
             !$omp end parallel
                
          case default
             stop 'error in reconstruct: ndim not properly set'
       end select
       !
    endif
    !
    !
  end function reconstruct
!
!-----------------------------------------------------------------------
!
  function reconstruct2(fheader,fblocks,grid_out)
    !
    !reconstruct all the data directly from each block to an output grid grid_out
    !
    !slower than function reconstruct, but less memory consumption
    !
    !... arguments
    type(header) :: fheader
    type(block), dimension(:), allocatable :: fblocks
    type(grid), intent(in) :: grid_out    
    type(data) :: reconstruct2
    !
    ! ... local scalars
    integer(i4b) :: nw, n_blocks, ndim, max_level
    integer(i4b) :: i, j, k, ind1, err
    integer(i4b) :: iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1
    real(dp) :: xl, xr, yl, yr, zl, zr
    real(dp) :: acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff
    !
    ! ... local characters
    !
    ! ... local arrays
    integer(i4b), dimension(:), allocatable :: xindx_block, yindx_block, zindx_block
    !
    ! ... local logicals
    logical :: lfound
    !
    ! ... local derived types
    !
    !
    !get info from the header
    ndim = fheader%ndim            !number of dimensions
    nw = fheader%nw
    max_level=fheader%levmax
    allocate(xindx_block(max_level), stat=err)
    allocate(yindx_block(max_level), stat=err)
    allocate(zindx_block(max_level), stat=err)
    !
    write(*,*) '------------performing reconstruct2------------'
    write(*,*)
    write(*,'(a50,3i6)') 'maximum resolution on finest level', 2**(max_level-1)*fheader%domain_nx
    write(*,'(a50,3i6)') 'applied resolution', grid_out%xgrid%nd, grid_out%ygrid%nd, grid_out%zgrid%nd
    write(*,*)
    !
    reconstruct2%ndim = ndim
    reconstruct2%nw = nw
    allocate(reconstruct2%w_names(nw), stat=err)
    do i = 1, nw
       reconstruct2%w_names(i) = fheader%w_names(i)
    enddo    
    !
    !get the total number of blocks
    n_blocks = fheader%nleafs
    !
    !perform the required allocations
    allocate(reconstruct2%w_names(nw), stat=err)
    select case(ndim)
       case(1)
          stop 'not implemented yet'
       case(2)
          stop 'not implemented yet'
       case(3)          
          reconstruct2%data_shape = (/ grid_out%xgrid%nd, grid_out%ygrid%nd, grid_out%zgrid%nd,nw /)
!          write(*,*) reconstruct2%data_shape
          allocate(reconstruct2%data(grid_out%xgrid%nd,grid_out%ygrid%nd,grid_out%zgrid%nd,nw), stat=err)
          allocate(reconstruct2%mesh%xgrid%coord(grid_out%xgrid%nd), stat=err)
          allocate(reconstruct2%mesh%ygrid%coord(grid_out%ygrid%nd), stat=err)
          allocate(reconstruct2%mesh%zgrid%coord(grid_out%zgrid%nd), stat=err)
       case default
          stop 'error in reconstruct2: ndim not properly set'
    end select
    !
    !set the coordinate arrays
    reconstruct2%mesh = grid_out
    !
    !find the block index where each coordinate is located in
    !$omp parallel &
    !$omp private(j, k, ind1, xl, xr, yl, yr, zl, zr, lfound, acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff, iim2, iim1, ii, iip1, jjm2, jjm1, jj, jjp1, kkm2, kkm1, kk, kkp1)

    !$omp do schedule(static,8)
    do i=1, grid_out%xgrid%nd
       do j=1, grid_out%ygrid%nd
          do k=1, grid_out%zgrid%nd
             lfound = .false.
             do ind1=1, n_blocks
                xl = fblocks(ind1)%x_edgel(1)
                xr = fblocks(ind1)%x_edger(1)
                yl = fblocks(ind1)%x_edgel(2)
                yr = fblocks(ind1)%x_edger(2)
                zl = fblocks(ind1)%x_edgel(3)
                zr = fblocks(ind1)%x_edger(3)
                !
                if(grid_out%xgrid%coord(i).le.xr.and.grid_out%xgrid%coord(i).ge.xl.and.&
                   grid_out%ygrid%coord(j).le.yr.and.grid_out%ygrid%coord(j).ge.yl.and.&
                   grid_out%zgrid%coord(k).le.zr.and.grid_out%zgrid%coord(k).ge.zl) then
                   !
                   lfound=.true.
                   !perform interpolation
                   call find_index(grid_out%xgrid%coord(i), fblocks(ind1)%x_coord, fblocks(ind1)%nx, iim2, iim1, ii, iip1)
                   call find_index(grid_out%ygrid%coord(j), fblocks(ind1)%y_coord, fblocks(ind1)%ny, jjm2, jjm1, jj, jjp1)
                   call find_index(grid_out%zgrid%coord(k), fblocks(ind1)%z_coord, fblocks(ind1)%nz, kkm2, kkm1, kk, kkp1)
                   !
                   !                      call coeff3d_8p_lin(fblocks(ind1)%x_coord(iim1), fblocks(ind1)%x_coord(ii), &
                   !                                          fblocks(ind1)%y_coord(jjm1), fblocks(ind1)%y_coord(jj), &
                   !                                          fblocks(ind1)%z_coord(kkm1), fblocks(ind1)%z_coord(kk), &
                   !                                          grid_out%xgrid%coord(i), grid_out%ygrid%coord(j), grid_out%zgrid%coord(k), &
                   !                                          acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff)

                   !inverse distance weighting used to avoid extrapolation at the edges
                   call coeff3d_8p_idw(fblocks(ind1)%x_coord(iim1), fblocks(ind1)%y_coord(jjm1), fblocks(ind1)%z_coord(kkm1), &
                                       fblocks(ind1)%x_coord(ii),   fblocks(ind1)%y_coord(jjm1), fblocks(ind1)%z_coord(kkm1), &
                                       fblocks(ind1)%x_coord(iim1), fblocks(ind1)%y_coord(jj),   fblocks(ind1)%z_coord(kkm1), &
                                       fblocks(ind1)%x_coord(ii),   fblocks(ind1)%y_coord(jj),   fblocks(ind1)%z_coord(kkm1), &
                                       fblocks(ind1)%x_coord(iim1), fblocks(ind1)%y_coord(jjm1), fblocks(ind1)%z_coord(kk), &
                                       fblocks(ind1)%x_coord(ii),   fblocks(ind1)%y_coord(jjm1), fblocks(ind1)%z_coord(kk), &
                                       fblocks(ind1)%x_coord(iim1), fblocks(ind1)%y_coord(jj),   fblocks(ind1)%z_coord(kk), &
                                       fblocks(ind1)%x_coord(ii),   fblocks(ind1)%y_coord(jj),   fblocks(ind1)%z_coord(kk), &
                                       grid_out%xgrid%coord(i), grid_out%ygrid%coord(j), grid_out%zgrid%coord(k), &
                                       acoeff, bcoeff, ccoeff, dcoeff, ecoeff, fcoeff, gcoeff, hcoeff, fac=2.d0)
                   !
                   reconstruct2%data(i,j,k,:) = acoeff*fblocks(ind1)%fblock(iim1,jjm1,kkm1,:) + bcoeff*fblocks(ind1)%fblock(ii,jjm1,kkm1,:) + &
                                                ccoeff*fblocks(ind1)%fblock(iim1,jj  ,kkm1,:) + dcoeff*fblocks(ind1)%fblock(ii,jj  ,kkm1,:) + &
                                                ecoeff*fblocks(ind1)%fblock(iim1,jjm1,kk  ,:) + fcoeff*fblocks(ind1)%fblock(ii,jjm1,kk  ,:) + &
                                                gcoeff*fblocks(ind1)%fblock(iim1,jj  ,kk  ,:) + hcoeff*fblocks(ind1)%fblock(ii,jj  ,kk  ,:)
                   exit
                endif
             enddo
             if(.not.lfound) stop 'error in reconstruct2: no block found for given coordinates'
          enddo
       enddo
    enddo
    !$omp enddo
    !$omp end parallel
    !
  end function reconstruct2
!
!-----------------------------------------------------------------------
!
  subroutine print_reconstruct(fdata)
    !
    ! ... arguments
    type(data) :: fdata
    !
    ! ... local scalars
    integer(i4b) :: ind1, ind2, ind3
    !
    ! ... local arrays
    !
    ! ... local characters
    character(len=100) :: fstr1, fstr2
    !
    !get the number of blocks
    write(*,*) '--------------DATA------------'
    write(*,*)
    !
    !
    if(fdata%nw < 10) then
       write(fstr1,'(a1, i1, a4)') '(', fdata%nw+1, 'a16)'                
       write(fstr2,'(a1, i1, a7)') '(', fdata%nw+1, 'es16.4)'
    elseif(fdata%nw < 100) then
       write(fstr1,'(a1, i2, a4)') '(', fdata%nw+1, 'a16)'                
       write(fstr2,'(a1, i2, a7)') '(', fdata%nw+1, 'es16.4)'
    endif
    !
    !
    !
    select case(fdata%ndim)
       case(1)
          write(*,fstr1) 'x', fdata%w_names
          do ind1=1, fdata%data_shape(3)
!             write(*,fstr2) fmesh%x(ind1), fdata%data(1,1,ind1,:)
          enddo 
       case(2)
          write(*,*) 'print_reconstruct 2d to do'
       case(3)
          ind2 = int(fdata%data_shape(2)/2.)
          ind3 = int(fdata%data_shape(3)/2.)             
!          write(*,*) 'at (y,z)=', fdata%mesh%ygrid%coord(ind2), fdata%mesh%zgrid%coord(ind3)
          write(*,fstr1) 'x', fdata%w_names             
          do ind1 = 1, fdata%data_shape(1)
!             write(*,fstr2) fdata%mesh%xgrid%coord(ind1), fdata%data(ind1,ind2,ind3,:)
          enddo
    end select
    !
  end subroutine print_reconstruct
!
!-----------------------------------------------------------------------
!
  function refined_block(ndim,fblock,interp_method)
    !
    ! ... arguments
    integer(i4b), intent(in) :: ndim
    type(block), intent(in) :: fblock
    type(block) :: refined_block
    character(len=*), intent(in), optional :: interp_method
    !
    ! ... local scalars
    integer(i4b) :: nx1, nx2, nx3, nw, err, i
    integer(i4b) :: iinterp_method
    !
    !
    if(.not.present(interp_method)) then
       iinterp_method = 1   !use linear interpolations
    elseif(trim(interp_method).eq.'lin') then
       iinterp_method = 1
    elseif(trim(interp_method).eq.'const') then
       iinterp_method = 2       
    else
       stop 'error in refined_block: interp_method not properly specified'
    endif
    !
    !
    select case(ndim)
       case(1)
          nx1 = fblock%block_shape(3)
          nw = fblock%block_shape(4)
          !
          refined_block%block_shape = (/ 1, 1, 2*nx1, nw /)
          allocate(refined_block%fblock(1, 1, 2*nx1, nw), stat=err)
          !
          call block_interp1d_const(nx1, refined_block, fblock)
          !
       case(2)
          nx1 = fblock%block_shape(2)
          nx2 = fblock%block_shape(3)
          nw = fblock%block_shape(4)
          !
          refined_block%block_shape = (/ 1, 2*nx1, 2*nx2, nw /)
          allocate(refined_block%fblock(1, 2*nx1, 2*nx2, nw), stat=err)
          !
          call block_interp2d_const(nx1, nx2, refined_block, fblock)
          !
       case(3)
          nx1 = fblock%block_shape(1)
          nx2 = fblock%block_shape(2)
          nx3 = fblock%block_shape(3)
          nw = fblock%block_shape(4)
!          write(*,*) nx1, nx2, nx3
!          write(*,*) fblock%x_coord
!          write(*,*) fblock%x_edgel
!          write(*,*) fblock%indexes
!          write(*,*)


          !refine the mesh
          refined_block%level = fblock%level
          allocate(refined_block%ncells(3), stat=err)
          refined_block%ncells = 2*fblock%ncells
          refined_block%indexes = fblock%indexes
          refined_block%nx = 2*nx1
          refined_block%ny = 2*nx2
          refined_block%nz = 2*nx3
          allocate(refined_block%x_center(3), stat=err)
          allocate(refined_block%x_edgel(3), stat=err)
          allocate(refined_block%x_edger(3), stat=err)
          refined_block%x_center=fblock%x_center
          refined_block%x_edgel=fblock%x_edgel
          refined_block%x_edger=fblock%x_edger

          allocate(refined_block%x_coord(refined_block%nx), stat=err)
          allocate(refined_block%y_coord(refined_block%ny), stat=err)
          allocate(refined_block%z_coord(refined_block%nz), stat=err)          
          do i=1, refined_block%nx
             refined_block%x_coord(i) = refined_block%x_edgel(1) + (2*i-1)/2.d0 * (refined_block%x_edger(1)-refined_block%x_edgel(1))/refined_block%nx
          enddo
          do i=1, refined_block%ny
             refined_block%y_coord(i) = refined_block%x_edgel(2) + (2*i-1)/2.d0 * (refined_block%x_edger(2)-refined_block%x_edgel(2))/refined_block%ny
          enddo
          do i=1, refined_block%nz
             refined_block%z_coord(i) = refined_block%x_edgel(3) + (2*i-1)/2.d0 * (refined_block%x_edger(3)-refined_block%x_edgel(3))/refined_block%nz
          enddo
          !
          !
          !
          refined_block%block_shape = (/ 2*nx1, 2*nx2, 2*nx3, nw /)
          allocate(refined_block%fblock(2*nx1, 2*nx2, 2*nx3, nw), stat=err)
          !
          select case(iinterp_method)
             case(1)
                call block_interp3d_lin(nx1, nx2, nx3, refined_block, fblock)
             case(2)
                call block_interp3d_const(nx1, nx2, nx3, refined_block, fblock)
             case default
                stop 'error in refined_block: iinterp_method not properly specified'
          end select
          !
       case default
          stop 'error in refined_block: ndim not properly set'
    end select
    !
    !
  end function refined_block
!
!-----------------------------------------------------------------------
!
  subroutine block_interp1d_const(nx, refined_block, fblock)
    !
    !constant interpolation of refined block
    !
    ! ... arguments
    integer(i4b), intent(in) :: nx
    type(block), intent(in) :: fblock
    type(block), intent(inout) :: refined_block
    !
    ! ... local scalars
    integer(i4b) :: i, ii
    !
    ii = 1
    do i = 1, nx
       refined_block%fblock(1,1,ii,:) = fblock%fblock(1, 1, i,:)
       refined_block%fblock(1,1,ii+1,:) = fblock%fblock(1,1, i,:)
       ii=ii+2
    enddo
    !
  end subroutine block_interp1d_const
!
!-----------------------------------------------------------------------
!
  subroutine block_interp2d_const(nx, ny, refined_block, fblock)
    !
    !constant interpolation of refined block
    !
    ! ... arguments
    integer(i4b), intent(in) :: nx, ny
    type(block), intent(in) :: fblock
    type(block), intent(inout) :: refined_block
    !
    ! ... local scalars
    integer(i4b) :: i, j, ii, jj
    !
    ii = 1
    do i = 1, nx
       jj = 1
       do j = 1, ny
          refined_block%fblock(1,ii,jj,:) = fblock%fblock(1, i, j,:)
          refined_block%fblock(1,ii,jj+1,:) = fblock%fblock(1, i, j,:)
          refined_block%fblock(1,ii+1,jj,:) = fblock%fblock(1, i, j,:)
          refined_block%fblock(1,ii+1,jj+1,:) = fblock%fblock(1, i, j,:)
          jj=jj+2
       enddo
       ii=ii+2
    enddo
    !
  end subroutine block_interp2d_const
!
!-----------------------------------------------------------------------
!
  subroutine block_interp3d_const(nx, ny, nz, refined_block, fblock)
    !
    !constant interpolation of refined block
    !
    ! ... arguments
    integer(i4b), intent(in) :: nx, ny, nz
    type(block), intent(in) :: fblock
    type(block), intent(inout) :: refined_block
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, ii, jj, kk
    !
    ii = 1
    do i = 1, nx
       jj = 1
       do j = 1, ny
          kk = 1
          do k = 1, nz
             refined_block%fblock(ii,jj,kk,:) = fblock%fblock(i, j, k,:)
             refined_block%fblock(ii,jj,kk+1,:) = fblock%fblock(i, j, k,:)
             refined_block%fblock(ii,jj+1,kk,:) = fblock%fblock(i, j, k,:)
             refined_block%fblock(ii+1,jj,kk,:) = fblock%fblock(i, j, k,:)
             refined_block%fblock(ii+1,jj+1,kk,:) = fblock%fblock(i, j, k,:)
             refined_block%fblock(ii+1,jj,kk+1,:) = fblock%fblock(i, j, k,:)
             refined_block%fblock(ii,jj+1,kk+1,:) = fblock%fblock(i, j, k,:)
             refined_block%fblock(ii+1,jj+1,kk+1,:) = fblock%fblock(i, j, k,:)
             kk=kk+2
          enddo
          jj=jj+2
       enddo
       ii=ii+2
    enddo
    !
  end subroutine block_interp3d_const
!
!-----------------------------------------------------------------------
!
  subroutine block_interp3d_lin(nx, ny, nz, refined_block, fblock)
    !
    !linear interpolation of refined block, when possible
    !
    ! ... arguments
    integer(i4b), intent(in) :: nx, ny, nz
    type(block), intent(in) :: fblock
    type(block), intent(inout) :: refined_block
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, ii, jj, kk, im1, jm1, km1, ip1, jp1, kp1
    logical :: lcorner1, lcorner2, lcorner3, lcorner4, lcorner5, lcorner6, lcorner7, lcorner8, &
         ledge1, ledge2, ledge3, ledge4, ledge5, ledge6, ledge7, ledge8, ledge9, ledge10, ledge11, ledge12, &
         lplane1, lplane2, lplane3, lplane4, lplane5, lplane6
    real(dp), parameter :: w1=1.d0/4.d0, w2=3.d0/4.d0, w3=1.d0/16.d0, w4=3.d0/16.d0, w5=9.d0/16.d0, &
         w6=1.d0/64.d0, w7=3.d0/64.d0, w8=9.d0/64.d0, w9=27.d0/64.d0
    !
    !
    !
    ii = 1
    do i = 1, nx
       im1=i-1
       ip1=i+1
       jj = 1
       do j = 1, ny
          jm1=j-1
          jp1=j+1
          kk = 1
          do k = 1, nz
             km1=k-1
             kp1=k+1
             lcorner1 = i.eq.1.and.j.eq.1.and.k.eq.1
             lcorner2 = i.eq.nx.and.j.eq.1.and.k.eq.1
             lcorner3 = i.eq.nx.and.j.eq.ny.and.k.eq.1
             lcorner4 = i.eq.1.and.j.eq.ny.and.k.eq.1
             lcorner5 = i.eq.1.and.j.eq.1.and.k.eq.nz
             lcorner6 = i.eq.nx.and.j.eq.1.and.k.eq.nz
             lcorner7 = i.eq.nx.and.j.eq.ny.and.k.eq.nz
             lcorner8 = i.eq.1.and.j.eq.ny.and.k.eq.nz

             ledge1 = j.eq.1.and.k.eq.1
             ledge2 = i.eq.nx.and.k.eq.1
             ledge3 = j.eq.ny.and.k.eq.1
             ledge4 = i.eq.1.and.k.eq.1
             ledge5 = j.eq.1.and.k.eq.nz
             ledge6 = i.eq.nx.and.k.eq.nz
             ledge7 = j.eq.ny.and.k.eq.nz
             ledge8 = i.eq.1.and.k.eq.nz
             ledge9 = i.eq.1.and.j.eq.1
             ledge10= i.eq.nx.and.j.eq.1
             ledge11= i.eq.nx.and.j.eq.ny
             ledge12= i.eq.1.and.j.eq.ny

             lplane1 = k.eq.1
             lplane2 = k.eq.nz
             lplane3 = j.eq.1
             lplane4 = j.eq.ny
             lplane5 = i.eq.1
             lplane6 = i.eq.nx
             !corners of the block
             if(lcorner1) then
                refined_block%fblock(ii,jj,kk,:) = fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w1*fblock%fblock(ip1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w1*fblock%fblock(i, jp1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w3*fblock%fblock(ip1, jp1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w1*fblock%fblock(i, j, kp1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(ip1, j, kp1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(i, jp1, kp1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w6*fblock%fblock(ip1, jp1, kp1,:) + &
                                                         w9*fblock%fblock(i, j, k,:) + &
                                                         w7*fblock%fblock(i, jp1, kp1,:) + &
                                                         w7*fblock%fblock(ip1, j, kp1,:) + &
                                                         w7*fblock%fblock(ip1, jp1, k,:) + &
                                                         w8*fblock%fblock(ip1, j, k,:) + &
                                                         w8*fblock%fblock(i, jp1, k,:) + &
                                                         w8*fblock%fblock(i, j, kp1,:)
             elseif(lcorner2) then
                refined_block%fblock(ii,jj,kk,:) = w1*fblock%fblock(im1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w3*fblock%fblock(im1, jp1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w1*fblock%fblock(i, jp1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(im1, j, kp1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w1*fblock%fblock(i, j, kp1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w6*fblock%fblock(im1, jp1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(im1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jp1, kp1,:) + &
                                                       w7*fblock%fblock(im1, jp1, k,:) + &
                                                       w8*fblock%fblock(im1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w3*fblock%fblock(i, jp1, kp1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(lcorner3) then
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(im1, jm1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w1*fblock%fblock(i, jm1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w1*fblock%fblock(im1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w6*fblock%fblock(im1, jm1, kp1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, j, kp1,:) + &
                                                     w7*fblock%fblock(i, jm1, kp1,:) + &
                                                     w7*fblock%fblock(im1, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, kp1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(i, jm1, kp1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(im1, j, kp1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w1*fblock%fblock(i, j, kp1,:) + &
                     w2*fblock%fblock(i, j, k,:)
             elseif(lcorner4) then
                refined_block%fblock(ii,jj,kk,:) = w1*fblock%fblock(i, jm1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) =  w3*fblock%fblock(ip1, jm1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w1*fblock%fblock(ip1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) =  w3*fblock%fblock(i, jm1, kp1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w6*fblock%fblock(ip1, jm1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jm1, kp1,:) + &
                                                       w7*fblock%fblock(ip1, jm1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w1*fblock%fblock(i, j, kp1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) =  w3*fblock%fblock(ip1, j, kp1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(lcorner5) then
                refined_block%fblock(ii,jj,kk,:) = w1*fblock%fblock(i, j, km1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) =  w3*fblock%fblock(ip1, j, km1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) =  w3*fblock%fblock(i, jp1, km1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w6*fblock%fblock(ip1, jp1, km1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, jp1, k,:) + &
                                                       w7*fblock%fblock(ip1, j, km1,:) + &
                                                       w7*fblock%fblock(i, jp1, km1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii,jj,kk+1,:) = fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w1*fblock%fblock(ip1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w1*fblock%fblock(i, jp1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) =  w3*fblock%fblock(ip1, jp1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(lcorner6) then
                refined_block%fblock(ii,jj,kk,:) =  w3*fblock%fblock(im1, j, km1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w1*fblock%fblock(i, j, km1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w6*fblock%fblock(im1, jp1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, jp1, k,:) + &
                                                     w7*fblock%fblock(im1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jp1, km1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jp1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj+1,kk,:) =  w3*fblock%fblock(i, jp1, km1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w1*fblock%fblock(im1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) =  w3*fblock%fblock(im1, jp1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w1*fblock%fblock(i, jp1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
             elseif(lcorner7) then
                refined_block%fblock(ii,jj,kk,:) = w6*fblock%fblock(im1, jm1, km1,:) + &
                                                   w9*fblock%fblock(i, j, k,:) + &
                                                   w7*fblock%fblock(i, jm1, km1,:) + &
                                                   w7*fblock%fblock(im1, j, km1,:) + &
                                                   w7*fblock%fblock(im1, jm1, k,:) + &
                                                   w8*fblock%fblock(im1, j, k,:) + &
                                                   w8*fblock%fblock(i, jm1, k,:) + &
                                                   w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj,kk,:) =  w3*fblock%fblock(i, jm1, km1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) =  w3*fblock%fblock(im1, j, km1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w1*fblock%fblock(i, j, km1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) =  w3*fblock%fblock(im1, jm1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w1*fblock%fblock(i, jm1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w1*fblock%fblock(im1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = fblock%fblock(i, j, k,:)
             elseif(lcorner8) then
                refined_block%fblock(ii,jj,kk,:) =  w3*fblock%fblock(i, jm1, km1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w6*fblock%fblock(ip1, jm1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(ip1, jm1, k,:) + &
                                                     w7*fblock%fblock(ip1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jm1, km1,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:) + &
                                                     w8*fblock%fblock(ip1, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w1*fblock%fblock(i, j, km1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) =  w3*fblock%fblock(ip1, j, km1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w1*fblock%fblock(i, jm1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) =  w3*fblock%fblock(ip1, jm1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w1*fblock%fblock(ip1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
             elseif(ledge1) then
                refined_block%fblock(ii,jj,kk,:) = w1*fblock%fblock(im1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w1*fblock%fblock(ip1, j, k,:) + &
                                                     w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w3*fblock%fblock(im1, jp1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w3*fblock%fblock(ip1, jp1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(im1, j, kp1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(ip1, j, kp1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w6*fblock%fblock(im1, jp1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(im1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jp1, kp1,:) + &
                                                       w7*fblock%fblock(im1, jp1, k,:) + &
                                                       w8*fblock%fblock(im1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w6*fblock%fblock(ip1, jp1, kp1,:) + &
                                                         w9*fblock%fblock(i, j, k,:) + &
                                                         w7*fblock%fblock(i, jp1, kp1,:) + &
                                                         w7*fblock%fblock(ip1, j, kp1,:) + &
                                                         w7*fblock%fblock(ip1, jp1, k,:) + &
                                                         w8*fblock%fblock(ip1, j, k,:) + &
                                                         w8*fblock%fblock(i, jp1, k,:) + &
                                                         w8*fblock%fblock(i, j, kp1,:)
             elseif(ledge2) then
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(im1, jm1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w1*fblock%fblock(i, jm1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w3*fblock%fblock(im1, jp1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w1*fblock%fblock(i, jp1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w6*fblock%fblock(im1, jm1, kp1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, j, kp1,:) + &
                                                     w7*fblock%fblock(i, jm1, kp1,:) + &
                                                     w7*fblock%fblock(im1, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, kp1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(i, jm1, kp1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w6*fblock%fblock(im1, jp1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(im1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jp1, kp1,:) + &
                                                       w7*fblock%fblock(im1, jp1, k,:) + &
                                                       w8*fblock%fblock(im1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w3*fblock%fblock(i, jp1, kp1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(ledge3) then                    
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(im1, jm1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w3*fblock%fblock(ip1, jm1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w1*fblock%fblock(im1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w1*fblock%fblock(ip1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w6*fblock%fblock(im1, jm1, kp1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, j, kp1,:) + &
                                                     w7*fblock%fblock(i, jm1, kp1,:) + &
                                                     w7*fblock%fblock(im1, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, kp1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w6*fblock%fblock(ip1, jm1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jm1, kp1,:) + &
                                                       w7*fblock%fblock(ip1, jm1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(im1, j, kp1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w3*fblock%fblock(ip1, j, kp1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(ledge4) then
                refined_block%fblock(ii,jj,kk,:) = w1*fblock%fblock(i, jm1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w3*fblock%fblock(ip1, jm1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w1*fblock%fblock(i, jp1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w3*fblock%fblock(ip1, jp1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(i, jm1, kp1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w6*fblock%fblock(ip1, jm1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jm1, kp1,:) + &
                                                       w7*fblock%fblock(ip1, jm1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(i, jp1, kp1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w6*fblock%fblock(ip1, jp1, kp1,:) + &
                                                         w9*fblock%fblock(i, j, k,:) + &
                                                         w7*fblock%fblock(i, jp1, kp1,:) + &
                                                         w7*fblock%fblock(ip1, j, kp1,:) + &
                                                         w7*fblock%fblock(ip1, jp1, k,:) + &
                                                         w8*fblock%fblock(ip1, j, k,:) + &
                                                         w8*fblock%fblock(i, jp1, k,:) + &
                                                         w8*fblock%fblock(i, j, kp1,:)
             elseif(ledge5) then
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(im1, j, km1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w3*fblock%fblock(ip1, j, km1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w6*fblock%fblock(im1, jp1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, jp1, k,:) + &
                                                     w7*fblock%fblock(im1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jp1, km1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jp1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w6*fblock%fblock(ip1, jp1, km1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, jp1, k,:) + &
                                                       w7*fblock%fblock(ip1, j, km1,:) + &
                                                       w7*fblock%fblock(i, jp1, km1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii,jj,kk+1,:) = w1*fblock%fblock(im1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w1*fblock%fblock(ip1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(im1, jp1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w3*fblock%fblock(ip1, jp1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(ledge6) then
                refined_block%fblock(ii,jj,kk,:) = w6*fblock%fblock(im1, jm1, km1,:) + &
                                                   w9*fblock%fblock(i, j, k,:) + &
                                                   w7*fblock%fblock(i, jm1, km1,:) + &
                                                   w7*fblock%fblock(im1, j, km1,:) + &
                                                   w7*fblock%fblock(im1, jm1, k,:) + &
                                                   w8*fblock%fblock(im1, j, k,:) + &
                                                   w8*fblock%fblock(i, jm1, k,:) + &
                                                   w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj,kk,:) = w3*fblock%fblock(i, jm1, km1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w6*fblock%fblock(im1, jp1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, jp1, k,:) + &
                                                     w7*fblock%fblock(im1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jp1, km1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jp1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w3*fblock%fblock(i, jp1, km1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(im1, jm1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w1*fblock%fblock(i, jm1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(im1, jp1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w1*fblock%fblock(i, jp1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
             elseif(ledge7) then
                refined_block%fblock(ii,jj,kk,:) = w6*fblock%fblock(im1, jm1, km1,:) + &
                                                   w9*fblock%fblock(i, j, k,:) + &
                                                   w7*fblock%fblock(i, jm1, km1,:) + &
                                                   w7*fblock%fblock(im1, j, km1,:) + &
                                                   w7*fblock%fblock(im1, jm1, k,:) + &
                                                   w8*fblock%fblock(im1, j, k,:) + &
                                                   w8*fblock%fblock(i, jm1, k,:) + &
                                                   w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj,kk,:) = w6*fblock%fblock(ip1, jm1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(ip1, jm1, k,:) + &
                                                     w7*fblock%fblock(ip1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jm1, km1,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:) + &
                                                     w8*fblock%fblock(ip1, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w3*fblock%fblock(im1, j, km1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w3*fblock%fblock(ip1, j, km1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(im1, jm1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(ip1, jm1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w1*fblock%fblock(im1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w1*fblock%fblock(ip1, j, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
             elseif(ledge8) then
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(i, jm1, km1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w6*fblock%fblock(ip1, jm1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(ip1, jm1, k,:) + &
                                                     w7*fblock%fblock(ip1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jm1, km1,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:) + &
                                                     w8*fblock%fblock(ip1, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w3*fblock%fblock(i, jp1, km1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w6*fblock%fblock(ip1, jp1, km1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, jp1, k,:) + &
                                                       w7*fblock%fblock(ip1, j, km1,:) + &
                                                       w7*fblock%fblock(i, jp1, km1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii,jj,kk+1,:) = w1*fblock%fblock(i, jm1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(ip1, jm1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w1*fblock%fblock(i, jp1, k,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w3*fblock%fblock(ip1, jp1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(ledge9) then
                refined_block%fblock(ii,jj,kk,:) = w1*fblock%fblock(i, j, km1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w3*fblock%fblock(ip1, j, km1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w3*fblock%fblock(i, jp1, km1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w6*fblock%fblock(ip1, jp1, km1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, jp1, k,:) + &
                                                       w7*fblock%fblock(ip1, j, km1,:) + &
                                                       w7*fblock%fblock(i, jp1, km1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii,jj,kk+1,:) = w1*fblock%fblock(i, j, kp1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(ip1, j, kp1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(i, jp1, kp1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w6*fblock%fblock(ip1, jp1, kp1,:) + &
                                                         w9*fblock%fblock(i, j, k,:) + &
                                                         w7*fblock%fblock(i, jp1, kp1,:) + &
                                                         w7*fblock%fblock(ip1, j, kp1,:) + &
                                                         w7*fblock%fblock(ip1, jp1, k,:) + &
                                                         w8*fblock%fblock(ip1, j, k,:) + &
                                                         w8*fblock%fblock(i, jp1, k,:) + &
                                                         w8*fblock%fblock(i, j, kp1,:)
             elseif(ledge10) then
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(im1, j, km1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w1*fblock%fblock(i, j, km1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w6*fblock%fblock(im1, jp1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, jp1, k,:) + &
                                                     w7*fblock%fblock(im1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jp1, km1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jp1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w3*fblock%fblock(i, jp1, km1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(im1, j, kp1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w1*fblock%fblock(i, j, kp1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w6*fblock%fblock(im1, jp1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(im1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jp1, kp1,:) + &
                                                       w7*fblock%fblock(im1, jp1, k,:) + &
                                                       w8*fblock%fblock(im1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w3*fblock%fblock(i, jp1, kp1,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(ledge11) then
                refined_block%fblock(ii,jj,kk,:) = w6*fblock%fblock(im1, jm1, km1,:) + &
                                                   w9*fblock%fblock(i, j, k,:) + &
                                                   w7*fblock%fblock(i, jm1, km1,:) + &
                                                   w7*fblock%fblock(im1, j, km1,:) + &
                                                   w7*fblock%fblock(im1, jm1, k,:) + &
                                                   w8*fblock%fblock(im1, j, k,:) + &
                                                   w8*fblock%fblock(i, jm1, k,:) + &
                                                   w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj,kk,:) = w3*fblock%fblock(i, jm1, km1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w3*fblock%fblock(im1, j, km1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w1*fblock%fblock(i, j, km1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w6*fblock%fblock(im1, jm1, kp1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, j, kp1,:) + &
                                                     w7*fblock%fblock(i, jm1, kp1,:) + &
                                                     w7*fblock%fblock(im1, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, kp1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(i, jm1, kp1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(im1, j, kp1,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w1*fblock%fblock(i, j, kp1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
             elseif(ledge12) then
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(i, jm1, km1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w6*fblock%fblock(ip1, jm1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(ip1, jm1, k,:) + &
                                                     w7*fblock%fblock(ip1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jm1, km1,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:) + &
                                                     w8*fblock%fblock(ip1, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w1*fblock%fblock(i, j, km1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w3*fblock%fblock(ip1, j, km1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, km1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(i, jm1, kp1,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w6*fblock%fblock(ip1, jm1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jm1, kp1,:) + &
                                                       w7*fblock%fblock(ip1, jm1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w1*fblock%fblock(i, j, kp1,:) + &
                                                   w2*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w3*fblock%fblock(ip1, j, kp1,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, j, kp1,:) + &
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(lplane1) then
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(im1, jm1, k,:) + &
                                                   w4*fblock%fblock(im1, j, k,:) + &
                                                   w4*fblock%fblock(i, jm1, k,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w3*fblock%fblock(ip1, jm1, k,:) + &
                                                   w4*fblock%fblock(ip1, j, k,:) + &
                                                   w4*fblock%fblock(i, jm1, k,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w3*fblock%fblock(im1, jp1, k,:) + &
                                                   w4*fblock%fblock(im1, j, k,:) + &
                                                   w4*fblock%fblock(i, jp1, k,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w3*fblock%fblock(ip1, jp1, k,:) + &
                                                   w4*fblock%fblock(ip1, j, k,:) + &
                                                   w4*fblock%fblock(i, jp1, k,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w6*fblock%fblock(im1, jm1, kp1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, j, kp1,:) + &
                                                     w7*fblock%fblock(i, jm1, kp1,:) + &
                                                     w7*fblock%fblock(im1, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, kp1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w6*fblock%fblock(ip1, jm1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jm1, kp1,:) + &
                                                       w7*fblock%fblock(ip1, jm1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w6*fblock%fblock(im1, jp1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(im1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jp1, kp1,:) + &
                                                       w7*fblock%fblock(im1, jp1, k,:) + &
                                                       w8*fblock%fblock(im1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w6*fblock%fblock(ip1, jp1, kp1,:) + &
                                                         w9*fblock%fblock(i, j, k,:) + &
                                                         w7*fblock%fblock(i, jp1, kp1,:) + &
                                                         w7*fblock%fblock(ip1, j, kp1,:) + &
                                                         w7*fblock%fblock(ip1, jp1, k,:) + &
                                                         w8*fblock%fblock(ip1, j, k,:) + &
                                                         w8*fblock%fblock(i, jp1, k,:) + &
                                                         w8*fblock%fblock(i, j, kp1,:)
             elseif(lplane2) then
                refined_block%fblock(ii,jj,kk,:) = w6*fblock%fblock(im1, jm1, km1,:) + &
                                                   w9*fblock%fblock(i, j, k,:) + &
                                                   w7*fblock%fblock(i, jm1, km1,:) + &
                                                   w7*fblock%fblock(im1, j, km1,:) + &
                                                   w7*fblock%fblock(im1, jm1, k,:) + &
                                                   w8*fblock%fblock(im1, j, k,:) + &
                                                   w8*fblock%fblock(i, jm1, k,:) + &
                                                   w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj,kk,:) = w6*fblock%fblock(ip1, jm1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(ip1, jm1, k,:) + &
                                                     w7*fblock%fblock(ip1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jm1, km1,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:) + &
                                                     w8*fblock%fblock(ip1, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w6*fblock%fblock(im1, jp1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, jp1, k,:) + &
                                                     w7*fblock%fblock(im1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jp1, km1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jp1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w6*fblock%fblock(ip1, jp1, km1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, jp1, k,:) + &
                                                       w7*fblock%fblock(ip1, j, km1,:) + &
                                                       w7*fblock%fblock(i, jp1, km1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(im1, jm1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &                                                  
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(ip1, jm1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jm1, k,:) + &                                                  
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(im1, jp1, k,:) + &
                                                     w4*fblock%fblock(im1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &                                                  
                                                     w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w3*fblock%fblock(ip1, jp1, k,:) + &
                                                     w4*fblock%fblock(ip1, j, k,:) + &
                                                     w4*fblock%fblock(i, jp1, k,:) + &                                                  
                                                     w5*fblock%fblock(i, j, k,:)
             elseif(lplane3) then
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(im1, j, km1,:) + &
                                                   w4*fblock%fblock(im1, j, k,:) + &
                                                   w4*fblock%fblock(i, j, km1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w3*fblock%fblock(ip1, j, km1,:) + &
                                                   w4*fblock%fblock(ip1, j, k,:) + &
                                                   w4*fblock%fblock(i, j, km1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w6*fblock%fblock(im1, jp1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, jp1, k,:) + &
                                                     w7*fblock%fblock(im1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jp1, km1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jp1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w6*fblock%fblock(ip1, jp1, km1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, jp1, k,:) + &
                                                       w7*fblock%fblock(ip1, j, km1,:) + &
                                                       w7*fblock%fblock(i, jp1, km1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(im1, j, kp1,:) + &
                                                   w4*fblock%fblock(im1, j, k,:) + &
                                                   w4*fblock%fblock(i, j, kp1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(ip1, j, kp1,:) + &
                                                   w4*fblock%fblock(ip1, j, k,:) + &
                                                   w4*fblock%fblock(i, j, kp1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w6*fblock%fblock(im1, jp1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(im1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jp1, kp1,:) + &
                                                       w7*fblock%fblock(im1, jp1, k,:) + &
                                                       w8*fblock%fblock(im1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w6*fblock%fblock(ip1, jp1, kp1,:) + &
                                                         w9*fblock%fblock(i, j, k,:) + &
                                                         w7*fblock%fblock(i, jp1, kp1,:) + &
                                                         w7*fblock%fblock(ip1, j, kp1,:) + &
                                                         w7*fblock%fblock(ip1, jp1, k,:) + &
                                                         w8*fblock%fblock(ip1, j, k,:) + &
                                                         w8*fblock%fblock(i, jp1, k,:) + &
                                                         w8*fblock%fblock(i, j, kp1,:)
             elseif(lplane4) then
                refined_block%fblock(ii,jj,kk,:) = w6*fblock%fblock(im1, jm1, km1,:) + &
                                                   w9*fblock%fblock(i, j, k,:) + &
                                                   w7*fblock%fblock(i, jm1, km1,:) + &
                                                   w7*fblock%fblock(im1, j, km1,:) + &
                                                   w7*fblock%fblock(im1, jm1, k,:) + &
                                                   w8*fblock%fblock(im1, j, k,:) + &
                                                   w8*fblock%fblock(i, jm1, k,:) + &
                                                   w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj,kk,:) = w6*fblock%fblock(ip1, jm1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(ip1, jm1, k,:) + &
                                                     w7*fblock%fblock(ip1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jm1, km1,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:) + &
                                                     w8*fblock%fblock(ip1, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) =  w3*fblock%fblock(im1, j, km1,:) + &
                                                   w4*fblock%fblock(im1, j, k,:) + &
                                                   w4*fblock%fblock(i, j, km1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) =  w3*fblock%fblock(ip1, j, km1,:) + &
                                                   w4*fblock%fblock(ip1, j, k,:) + &
                                                   w4*fblock%fblock(i, j, km1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w6*fblock%fblock(im1, jm1, kp1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, j, kp1,:) + &
                                                     w7*fblock%fblock(i, jm1, kp1,:) + &
                                                     w7*fblock%fblock(im1, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, kp1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w6*fblock%fblock(ip1, jm1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jm1, kp1,:) + &
                                                       w7*fblock%fblock(ip1, jm1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) =  w3*fblock%fblock(im1, j, kp1,:) + &
                                                   w4*fblock%fblock(im1, j, k,:) + &
                                                   w4*fblock%fblock(i, j, kp1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) =  w3*fblock%fblock(ip1, j, kp1,:) + &
                                                   w4*fblock%fblock(ip1, j, k,:) + &
                                                   w4*fblock%fblock(i, j, kp1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
             elseif(lplane5) then
                refined_block%fblock(ii,jj,kk,:) = w3*fblock%fblock(i, jm1, km1,:) + &
                                                   w4*fblock%fblock(i, jm1, k,:) + &
                                                   w4*fblock%fblock(i, j, km1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk,:) = w6*fblock%fblock(ip1, jm1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(ip1, jm1, k,:) + &
                                                     w7*fblock%fblock(ip1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jm1, km1,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:) + &
                                                     w8*fblock%fblock(ip1, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w3*fblock%fblock(i, jp1, km1,:) + &
                                                   w4*fblock%fblock(i, jp1, k,:) + &
                                                   w4*fblock%fblock(i, j, km1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w6*fblock%fblock(ip1, jp1, km1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, jp1, k,:) + &
                                                       w7*fblock%fblock(ip1, j, km1,:) + &
                                                       w7*fblock%fblock(i, jp1, km1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii,jj,kk+1,:) = w3*fblock%fblock(i, jm1, kp1,:) + &
                                                   w4*fblock%fblock(i, jm1, k,:) + &
                                                   w4*fblock%fblock(i, j, kp1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w6*fblock%fblock(ip1, jm1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jm1, kp1,:) + &
                                                       w7*fblock%fblock(ip1, jm1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w3*fblock%fblock(i, jp1, kp1,:) + &
                                                   w4*fblock%fblock(i, jp1, k,:) + &
                                                   w4*fblock%fblock(i, j, kp1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w6*fblock%fblock(ip1, jp1, kp1,:) + &
                                                         w9*fblock%fblock(i, j, k,:) + &
                                                         w7*fblock%fblock(i, jp1, kp1,:) + &
                                                         w7*fblock%fblock(ip1, j, kp1,:) + &
                                                         w7*fblock%fblock(ip1, jp1, k,:) + &
                                                         w8*fblock%fblock(ip1, j, k,:) + &
                                                         w8*fblock%fblock(i, jp1, k,:) + &
                                                         w8*fblock%fblock(i, j, kp1,:)
             elseif(lplane6) then
                refined_block%fblock(ii,jj,kk,:) = w6*fblock%fblock(im1, jm1, km1,:) + &
                                                   w9*fblock%fblock(i, j, k,:) + &
                                                   w7*fblock%fblock(i, jm1, km1,:) + &
                                                   w7*fblock%fblock(im1, j, km1,:) + &
                                                   w7*fblock%fblock(im1, jm1, k,:) + &
                                                   w8*fblock%fblock(im1, j, k,:) + &
                                                   w8*fblock%fblock(i, jm1, k,:) + &
                                                   w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj,kk,:) = w3*fblock%fblock(i, jm1, km1,:) + &
                                                   w4*fblock%fblock(i, jm1, k,:) + &
                                                   w4*fblock%fblock(i, j, km1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w6*fblock%fblock(im1, jp1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, jp1, k,:) + &
                                                     w7*fblock%fblock(im1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jp1, km1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jp1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w3*fblock%fblock(i, jp1, km1,:) + &
                                                   w4*fblock%fblock(i, jp1, k,:) + &
                                                   w4*fblock%fblock(i, j, km1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj,kk+1,:) = w6*fblock%fblock(im1, jm1, kp1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, j, kp1,:) + &
                                                     w7*fblock%fblock(i, jm1, kp1,:) + &
                                                     w7*fblock%fblock(im1, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, kp1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w3*fblock%fblock(i, jm1, kp1,:) + &
                                                   w4*fblock%fblock(i, jm1, k,:) + &
                                                   w4*fblock%fblock(i, j, kp1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w6*fblock%fblock(im1, jp1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(im1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jp1, kp1,:) + &
                                                       w7*fblock%fblock(im1, jp1, k,:) + &
                                                       w8*fblock%fblock(im1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w3*fblock%fblock(i, jp1, kp1,:) + &
                                                   w4*fblock%fblock(i, jp1, k,:) + &
                                                   w4*fblock%fblock(i, j, kp1,:) + &
                                                   w5*fblock%fblock(i, j, k,:)
             else
                refined_block%fblock(ii,jj,kk,:) = w6*fblock%fblock(im1, jm1, km1,:) + &
                                                   w9*fblock%fblock(i, j, k,:) + &
                                                   w7*fblock%fblock(i, jm1, km1,:) + &
                                                   w7*fblock%fblock(im1, j, km1,:) + &
                                                   w7*fblock%fblock(im1, jm1, k,:) + &
                                                   w8*fblock%fblock(im1, j, k,:) + &
                                                   w8*fblock%fblock(i, jm1, k,:) + &
                                                   w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj,kk,:) = w6*fblock%fblock(ip1, jm1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(ip1, jm1, k,:) + &
                                                     w7*fblock%fblock(ip1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jm1, km1,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:) + &
                                                     w8*fblock%fblock(ip1, j, k,:)
                refined_block%fblock(ii,jj+1,kk,:) = w6*fblock%fblock(im1, jp1, km1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, jp1, k,:) + &
                                                     w7*fblock%fblock(im1, j, km1,:) + &
                                                     w7*fblock%fblock(i, jp1, km1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jp1, k,:) + &
                                                     w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii+1,jj+1,kk,:) = w6*fblock%fblock(ip1, jp1, km1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, jp1, k,:) + &
                                                       w7*fblock%fblock(ip1, j, km1,:) + &
                                                       w7*fblock%fblock(i, jp1, km1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, km1,:)
                refined_block%fblock(ii,jj,kk+1,:) = w6*fblock%fblock(im1, jm1, kp1,:) + &
                                                     w9*fblock%fblock(i, j, k,:) + &
                                                     w7*fblock%fblock(im1, j, kp1,:) + &
                                                     w7*fblock%fblock(i, jm1, kp1,:) + &
                                                     w7*fblock%fblock(im1, jm1, k,:) + &
                                                     w8*fblock%fblock(i, j, kp1,:) + &
                                                     w8*fblock%fblock(im1, j, k,:) + &
                                                     w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii+1,jj,kk+1,:) = w6*fblock%fblock(ip1, jm1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(ip1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jm1, kp1,:) + &
                                                       w7*fblock%fblock(ip1, jm1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:) + &
                                                       w8*fblock%fblock(ip1, j, k,:) + &
                                                       w8*fblock%fblock(i, jm1, k,:)
                refined_block%fblock(ii,jj+1,kk+1,:) = w6*fblock%fblock(im1, jp1, kp1,:) + &
                                                       w9*fblock%fblock(i, j, k,:) + &
                                                       w7*fblock%fblock(im1, j, kp1,:) + &
                                                       w7*fblock%fblock(i, jp1, kp1,:) + &
                                                       w7*fblock%fblock(im1, jp1, k,:) + &
                                                       w8*fblock%fblock(im1, j, k,:) + &
                                                       w8*fblock%fblock(i, jp1, k,:) + &
                                                       w8*fblock%fblock(i, j, kp1,:)
                refined_block%fblock(ii+1,jj+1,kk+1,:) = w6*fblock%fblock(ip1, jp1, kp1,:) + &
                                                         w9*fblock%fblock(i, j, k,:) + &
                                                         w7*fblock%fblock(i, jp1, kp1,:) + &
                                                         w7*fblock%fblock(ip1, j, kp1,:) + &
                                                         w7*fblock%fblock(ip1, jp1, k,:) + &
                                                         w8*fblock%fblock(ip1, j, k,:) + &
                                                         w8*fblock%fblock(i, jp1, k,:) + &
                                                         w8*fblock%fblock(i, j, kp1,:)
             endif
             kk=kk+2
          enddo
          jj=jj+2
       enddo
       ii=ii+2
    enddo
    !
  end subroutine block_interp3d_lin
!
!-----------------------------------------------------------------------
!
  subroutine output3d(fname, my_alldata, nstep)
    !
    !... arguments
    character(len=*) :: fname
    integer(i4b), intent(in), optional :: nstep
    type(alldata), intent(in) :: my_alldata
    !
    ! ... local scalars
    integer(i4b) :: i, j, k, ii, jj, kk, nx, ny, nz, nw, err, irho, ivelx, ively, ivelz, itgas, itrad, istep
    !
    ! ... output to hdf5
    integer(hid_t) :: file_id, dspace_id, dset_id, aspace_id, attr_id, group_id
    integer(hsize_t), dimension(1) :: dims_scalars = (/1/)
    integer(hsize_t), dimension(1) :: dims_x, dims_y, dims_z
    integer(hsize_t), dimension(3) :: dims
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: x, y, z
    real(dp), dimension(:,:,:), allocatable :: rho3d, velx3d, vely3d, velz3d, tgas3d, trad3d
    !
    ! ... local derived types
    type(data) :: fdata
    !
    !----------------------transform data to 3d arrays----------------------
    !
    !only every nstep data point
    if(.not.present(nstep)) then
       istep=1
    else
       istep=nstep
    endif
    !
    fdata = my_alldata%data
    !
    nx = fdata%data_shape(3)/istep
    ny = fdata%data_shape(2)/istep
    nz = fdata%data_shape(1)/istep
    nw = fdata%data_shape(4)

    !
    write(*,*) 'dimensions', nx, ny, nz
    !
    !
    allocate(x(nx),stat=err)
    allocate(y(ny),stat=err)
    allocate(z(nz),stat=err)
    allocate(rho3d(nx,ny,nz), stat=err)
    allocate(velx3d(nx,ny,nz), stat=err)
    allocate(vely3d(nx,ny,nz), stat=err)
    allocate(velz3d(nx,ny,nz), stat=err)
    allocate(tgas3d(nx,ny,nz), stat=err)
    allocate(trad3d(nx,ny,nz), stat=err)
    !
    irho=-1
    ivelx=-1
    ively=-1
    ivelz=-1
    itrad=-1
    itgas=-1
    
    
    do i=1, nw
       if(trim(fdata%w_names(i)).eq.'rho') irho=i
       if(trim(fdata%w_names(i)).eq.'v1') ivelz=i
       if(trim(fdata%w_names(i)).eq.'v2') ively=i
       if(trim(fdata%w_names(i)).eq.'v3') ivelx=i
       if(trim(fdata%w_names(i)).eq.'Trad') itrad=i
       if(trim(fdata%w_names(i)).eq.'Tgas') itgas=i
    enddo
    !write(*,*) irho, ivelz, itrad, itgas
    
    ii=1
    do i=1, nx
       x(i) = fdata%mesh%zgrid%coord(ii)
       ii=ii+istep
    enddo
    
    ii=1
    do i=1, ny
       y(i) = fdata%mesh%ygrid%coord(ii)
       ii=ii+istep
    enddo
    
    ii=1
    do i=1, nz
       z(i) = fdata%mesh%xgrid%coord(ii)
       ii=ii+istep
    enddo
    !
    !
    !
    ii=1
    do i=1, nx
       jj=1
       do j=1, ny
          kk=1
          do k=1, nz
             rho3d(i,j,k) = fdata%data(kk,jj,ii,irho)
             velx3d(i,j,k) = fdata%data(kk,jj,ii,ivelx)
             vely3d(i,j,k) = fdata%data(kk,jj,ii,ively)
             velz3d(i,j,k) = fdata%data(kk,jj,ii,ivelz)
             tgas3d(i,j,k) = fdata%data(kk,jj,ii,itgas)
             trad3d(i,j,k) = fdata%data(kk,jj,ii,itrad)
             kk=kk+istep
          enddo
          jj=jj+istep
       enddo
       ii=ii+istep
    enddo
    !
    !-------------------------output to hdf5-file---------------------------
    !
    dims_x = (/ nx /)
    dims_y = (/ ny /)
    dims_z = (/ nz /)
    dims = (/ nx, ny, nz /)
    !
    write(*,*) 'save model to file ', trim(fname)
    write(*,*)
     
    call h5open_f (err)
       call h5fcreate_f(trim(fname), h5f_acc_trunc_f, file_id, err)
          call h5gcreate_f(file_id, 'dimensions', group_id, err)
             call h5screate_simple_f(1, dims_scalars, aspace_id, err)
                call h5acreate_f(group_id, 'nx', h5t_native_integer, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, nx, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'ny', h5t_native_integer, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, ny, dims_scalars, err)
                call h5aclose_f(attr_id, err)
                call h5acreate_f(group_id, 'nz', h5t_native_integer, aspace_id, &
                                 attr_id, err)
                   call h5awrite_f(attr_id, h5t_native_integer, nz, dims_scalars, err)
                call h5aclose_f(attr_id, err)
             call h5sclose_f(aspace_id, err)
          call h5gclose_f(group_id, err)
    !
          call h5gcreate_f(file_id, 'coordinates', group_id, err)
             call h5screate_simple_f(1, dims_x, dspace_id, err)
                call h5dcreate_f(group_id, 'x', h5t_native_double, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, x, dims_x, err)
                call h5dclose_f(dset_id, err)
             call h5sclose_f(dspace_id, err)
             call h5screate_simple_f(1, dims_y, dspace_id, err)
                call h5dcreate_f(group_id, 'y', h5t_native_double, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, y, dims_y, err)
                call h5dclose_f(dset_id, err)
             call h5sclose_f(dspace_id, err)
             call h5screate_simple_f(1, dims_z, dspace_id, err)
                call h5dcreate_f(group_id, 'z', h5t_native_double, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, z, dims_z, err)
                call h5dclose_f(dset_id, err)
             call h5sclose_f(dspace_id, err)
          call h5gclose_f(group_id, err)
    !
          call h5gcreate_f(file_id, 'model', group_id, err)
             call h5screate_simple_f(3, dims, dspace_id, err)
                call h5dcreate_f(group_id, 'rho', h5t_native_double, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, rho3d, dims, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'velx', h5t_native_double, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, velx3d, dims, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'vely', h5t_native_double, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, vely3d, dims, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'velz', h5t_native_double, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, velz3d, dims, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'tgas', h5t_native_double, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, tgas3d, dims, err)
                call h5dclose_f(dset_id, err)
                call h5dcreate_f(group_id, 'trad', h5t_native_double, dspace_id, &
                                 dset_id, err)
                   call h5dwrite_f(dset_id, h5t_native_double, trad3d, dims, err)
                call h5dclose_f(dset_id, err)
             call h5sclose_f(dspace_id, err)
    
         call h5gclose_f(group_id, err)
       call h5fclose_f(file_id, err)
    call h5close_f(err)
    !
  end subroutine output3d
!
!-----------------------------------------------------------------------
!
  function grid_internal(fheader, fblocks)
    !
    !
    ! ... arguments  
    type(block), dimension(:), allocatable, intent(in) :: fblocks
    type(header), intent(in) :: fheader
    type(grid) :: grid_internal
    !
    ! ... local scalars
    integer(i4b) :: i, err, nx_internal, ny_internal, nz_internal
    integer(i4b) :: istart, iend, jstart, jend, kstart, kend
    integer(i4b) :: ndim, n_blocks
    integer(i4b), dimension(:), allocatable :: nb
    !
    ! ... local arrays
    real(dp), dimension(:), allocatable :: x_internal, y_internal, z_internal
    !
    ! ... local derived types
    type(block) :: fblock
    !
    !get info from the header
    ndim = fheader%ndim
    grid_internal%ndim = ndim
    n_blocks=fheader%nleafs
    !
    allocate(nb(ndim), stat=err)
    nb = fheader%block_nx
    !    
    !
    !
    !specify number of grid points for output grid and calculate the input grid distributions
    select case(ndim)
       case(1)
          stop 'not implemented yet'
       case(2)
          stop 'not implemented yet'
       case(3)
          nx_internal = nb(1)*n_blocks
          ny_internal = nb(2)*n_blocks
          nz_internal = nb(3)*n_blocks
          allocate(x_internal(nx_internal), stat=err)
          allocate(y_internal(ny_internal), stat=err)
          allocate(z_internal(nz_internal), stat=err)
          !
          !get all available coordinates
          istart = 1
          jstart = 1
          kstart = 1
          iend = istart + nb(1)-1
          jend = jstart + nb(2)-1
          kend = kstart + nb(3)-1
          do i=1, n_blocks
             fblock = fblocks(i)
             x_internal(istart:iend) = fblock%x_coord
             y_internal(jstart:jend) = fblock%y_coord
             z_internal(kstart:kend) = fblock%z_coord
             istart = istart + nb(1)
             jstart = jstart + nb(2)
             kstart = kstart + nb(3)
             iend = iend + nb(1)
             jend = jend + nb(2)
             kend = kend + nb(3)
          enddo
          !
          !sort the coordinates in ascending order
          call quicksort(x_internal, nx_internal, 1, nx_internal)
          call quicksort(y_internal, ny_internal, 1, ny_internal)
          call quicksort(z_internal, nz_internal, 1, nz_internal)
          !
          !find only the unique coordinates
          call uniq_elements(nx_internal, grid_internal%xgrid%nd, x_internal, grid_internal%xgrid%coord)
          call uniq_elements(ny_internal, grid_internal%ygrid%nd, y_internal, grid_internal%ygrid%coord)
          call uniq_elements(nz_internal, grid_internal%zgrid%nd, z_internal, grid_internal%zgrid%coord)
          !
          !replace inner and outer boundary grid points such that exact boundary is matched
          grid_internal%xgrid%coord(1) = fheader%xprobmin(1)
          grid_internal%ygrid%coord(1) = fheader%xprobmin(2)
          grid_internal%zgrid%coord(1) = fheader%xprobmin(3)          
          grid_internal%xgrid%coord(grid_internal%xgrid%nd) = fheader%xprobmax(1)          
          grid_internal%ygrid%coord(grid_internal%ygrid%nd) = fheader%xprobmax(2)          
          grid_internal%zgrid%coord(grid_internal%zgrid%nd) = fheader%xprobmax(3)
          !calculate the probability density function for the unique coordinates
          call calc_pdf2(grid_internal%xgrid)
          call calc_pdf2(grid_internal%ygrid)
          call calc_pdf2(grid_internal%zgrid)
          !
          !call print_grid1d(grid_internal%xgrid)
          !call print_grid1d(grid_internal%ygrid)
          !call print_grid1d(grid_internal%zgrid)
          !
       case default
          stop 'error in grid_external: ndim not specified'
    end select
    !
    !
    !
  end function grid_internal
!
!-----------------------------------------------------------------------
!
  function grid_external(grid_in, nd)
    !
    !
    ! ... arguments
    type(grid) :: grid_in
    type(grid) :: grid_external
    integer(i4b), dimension(:), allocatable, intent(in) :: nd
    !
    ! ... local scalars
    integer(i4b) :: ndim
    !
    !
    !get info from the header
    ndim = grid_in%ndim
    !
    !
    !specify number of grid points for output grid and calculate the input grid distributions
    select case(ndim)
       case(1)
          stop 'not implemented yet'
       case(2)
          stop 'not implemented yet'
       case(3)
          grid_external%xgrid%nd = nd(1)
          grid_external%ygrid%nd = nd(2)
          grid_external%zgrid%nd = nd(3)
          grid_external%xgrid = frecalc_grid1d(grid_in%xgrid, grid_external%xgrid%nd)
          grid_external%ygrid = frecalc_grid1d(grid_in%ygrid, grid_external%ygrid%nd)
          grid_external%zgrid = frecalc_grid1d(grid_in%zgrid, grid_external%zgrid%nd)
          !
          !call print_grid1d(grid_external%xgrid)
          !call print_grid1d(grid_external%ygrid)
          !call print_grid1d(grid_external%zgrid)
          !          
       case default
          stop 'error in grid_external: ndim not specified'
    end select
    !
    !
    !
  end function grid_external


end module mod_amrvac_reader
