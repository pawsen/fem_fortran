module processor

  ! This module contains input and output subroutines.

  ! **************************************************
  ! Note:  It is not necessary to alter in this module
  ! **************************************************

  use fedata
  implicit none

  private
!!$   public :: plot, input, output, plotdef, ploteval, &
!!$             plotevec, ploteig, numnode, &
!!$             stopwatch, plotmatlabdef, plotmatlabeval, &
!!$             plotmatlabevec, plotmatlabeig, rho_input

  public :: plot,  output, numnode, &
       stopwatch, plotmatlabdef, plotmatlabeval, &
       plotmatlabevec, plotmatlabeig, plotanim, plotdef

contains


  subroutine output(plate_type)
    !---------------------------------------------------------------------
    !  This subroutine outputs the results.
    !use fedata

    integer,optional, INTENT(IN) :: plate_type
    integer :: i,j, e, nen 

    write (*, *)
    write (*, '("Writing fem output file to ", a)') trim(filename)//'.out'
    open (10, file = trim(filename)//'.out')

    if (banded == 1) then ! Angiv kun bw, hvis banded storage bruges.
       write(10, '("Bandwith", g15.0)') bw
    end if
    write (10,'("banded type =" ,i2)') banded

    !$$$$$$    write(10, '("Samlet areal" ,g15.9)') arealet

    write(10, '("NG " ,i6)') ng

    if ((elem_type == 'PLATE') .or.  (elem_type == 'PLATE_GMSH')) then

       write (10,  '("plate_type: 1 = kirchoff, 2 = mindlin")')
       write (10,'("plate type=" ,i2)') plate_type
       write (10,*)
       write (10, '(" Node      Displacements              Slope            ")')
       write (10, '("number     w-direction     x-direction      y-direction")')
       do i = 1, nn
          write (10, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i, d(3*i-2), d(3*i-1), d(3*i)
       end do

       write (10, '(" Node                Applied/reaction forces")')
       write (10, '("number     1-direction     2-direction     3-direction")')
       do i = 1, nn
          write (10, '(i6, 1x, f15.9, 1x, f15.9, 1x, f15.9)') i, p(3*i-2), p(3*i-1), p(3*i)
       end do

       write (10, '("               ___________     ___________    ___________")')
       write (10, '("Summation ", f15.9, 1x, f15.9, 1x, f15.9)') &
            sum(p(1:neqn-2:3)), sum(p(2:neqn-1:3)),sum(p(3:neqn:3))

       if (antype == 'PIEZO') then ! print extra info
          write (10,*)
          write (10, '(" Node            Potential")')
          do i =1,nn !neqn+1,neqn+nn
             write (10, '(i6, 1x, f25.15)') i, d(neqn+i)
          end do

          write (10, '(" dof            lambda")')
          do i = neqn+nn+1, neqn+nn + nb
             write (10, '(i6, 1x, f25.15)') i, d(i)
          end do
          write (10,*)
       end if


    else ! plane42

       write (10, '(" Node            Displacements                Applied/reaction forces")')
       write (10, '("number     1-direction     2-direction      1-direction     2-direction")')
       do i = 1, nn
          write (10, '(i6, 1x, f25.15, 1x, f25.15, 2x, f15.9, 1x, f25.15)') i, d(2*i-1), d(2*i), p(2*i-1), p(2*i)
       end do
       write (10, '("                                            ___________     ___________")')
       write (10, '("                              Summation ", f15.9, 1x, f15.9)') &
            sum(p(1:neqn-1:2)), sum(p(2:neqn:2))

    end if

    write (10, *)
    write (10, '("Element                    Element strain                                   Element Stress")')
    write (10, '("number      1-direction     2-direction     3-direction      1-direction     2-direction     3-direction")')
    do e = 1, ne
       select case( element(e)%id )
       case( 1 )
          write (10, '(1x, i6, 1x, f25.15, 34x, f25.15)') e, strain(e, 1), stress(e, 1)
       case( 2 )
          write (10, '(1x, i6, 3(1x, f25.15), 1x, 3(1x, f25.15))') e, (strain(e, i), i = 1, 3) ,(stress(e, i), i = 1, 3)
       case( 3 )
          write (10, '(1x, i6, 3(1x, f15.9), 1x, 3(1x, f15.9))') e, (strain(e, i), i = 1, 3) ,(stress(e, i), i = 1, 3)
       end select
    end do


    close (10)

  end subroutine output


  subroutine plot(what, dest, color, title, data1, data2, data3)
    !---------------------------------------------------------------------
    !  This subroutine acts as an interface for calling
    !  different plot routines/tools.
    !
    !  INPUT:
    !     what CHARACTER: Specifies what to plot, can be 'un/deformed',
    !                     'deformed', 'vector', 'elements', or 'eigenmode'
    !     dest CHARACTER: Specifies where to plot, can be 'xwin',
    !                     'ps', 'gif', 'matlab' or '?'
    !     				  NOTE: matlab does not support: 'deformed'
    !     color CHARACTER: Sets the color mode, can be 'color' or 'gray'
    !     title CHARACTER: Sets a title for the plot, may be empty
    !
    !  OPTIONAL INPUT: (depending on the kind of plot)
    !     'un/deformed' no arguments needed
    !     'deformed' no arguments needed
    !     'vector' data1 = principal stress in direction 1
    !              data2 = principal stress in direction 2
    !              data3 = angel of principal stress direction 1
    !     'elements' data1 = vonmises stress in each element
    !     'eigenmode' data1(1) = lambda value
    !                 data2    = eigen vector
    !                 data3(1) = display time
    !                 data3(2) = time increment
    !
    !  NOTE that data1, data2, data3 are arrays of reals (8)
    !--------------------------------------------------------------------- 
    character(len=*), intent(in) :: what, dest, color, title
    real(8), dimension(:), intent(in), optional :: data1, data2, data3

    real(8) :: rvalue, rvalue2
    integer :: destination
    logical :: colors

    select case( dest )
    case( 'xwin' )
       destination = 0
    case( 'ps' )
       destination = 1
    case( 'gif' )
       destination = 2
    case( '?' )
       destination = 3 
    case( 'matlab' )
       ! All matlab routines are called from this case-select !
       select case( what )
       case( 'un/deformed' )
          call plotmatlabdef(title)
       case( 'elements' )
          if(.not.present(data1)) then
             write(*, *) 'Error: Data is missing for element value matlab plot'
             stop
          end if
          ! Call element plot
          call plotmatlabeval(title,data1)
       case( 'vector' )
          if(.not.present(data1) .or. .not.present(data2) .or. .not.present(data3)) then
             write(*, *) 'Error: Data is missing for vector matlab plot'
             stop
          end if
          ! call vector plot
          call plotmatlabevec(title,data1,data2,data3)
       case( 'eigenmode' )
          if(.not.present(data1) .or. .not.present(data2) .or. .not.present(data3)) then
             write(*, *) 'Error: Data is missing for eigen mode plot'
             stop
          end if
          !call plotmatlabeig(title,data1,data2,data3) ! linux_fejl
       case default
          write(*, *) 'Error: Unknown matlab plot type specified'
          stop

       end select

    case default
       write(*, *) 'Error: Unknown destination for plotting.'
       stop
    end select

    select case( color )
    case( 'color' )
       colors = .true.
    case( 'gray' )
       colors = .false.
    case default
       write(*, *) 'Error: Illegal option for color/gray'
       stop
    end select

    select case( what )
    case( 'un/deformed' )
       if ( dest == 'matlab' ) then
       else
          call plotdef
       endif
    case( 'deformed' )
       if ( dest == 'matlab' ) then
       else
          call plotanim(0, destination, 1, colors, .true., .false., .false., .true., &
               title, 1.d0, (/1.d0/), (/1.d0/), (/1.d0/), 1.d0, 1.d0, (/1.d0/))
       endif
    case( 'vector' )
       if(.not.present(data1) .or. .not.present(data2) .or. .not.present(data3)) then
          write(*, *) 'Error: Data is missing for vector plot'
          stop
       end if
       if ( dest == 'matlab' ) then
       else
          rvalue = maxval(data1)
          rvalue = max(rvalue, maxval(data2))
          call plotanim(0, destination, 1, colors, .false., .true., .false., .true., &
               title, rvalue, data1, data2, data3, 1.d0, 1.d0, (/1.d0/))
          !call plotevec(data1, data2, data3)
       endif
    case( 'elements' )
       if(.not.present(data1)) then
          write(*, *) 'Error: Data is missing for element value plot'
          stop
       end if
       if ( dest == 'matlab' ) then
       else
          rvalue = minval(data1)
          rvalue2 = maxval(data1)
          call plotanim(0, destination, 1, colors, .false., .false., .true., .true., &
               title, 1.d0, (/1.d0/), (/1.d0/), (/1.d0/), rvalue, rvalue2, data1)
       endif
    case( 'eigenmode' )
       if(.not.present(data1) .or. .not.present(data2) .or. .not.present(data3)) then
          write(*, *) 'Error: Data is missing for eigen mode plot'
          stop
       end if
       if ( dest == 'matlab' ) then
       else
          !call ploteig(data1(1), data2, data3(1), data3(2))
       endif
    case default
       write(*, *) 'Error: Unknown plot type specified'
       stop
    end select
  end subroutine plot



  subroutine plotdef
    !---------------------------------------------------------------------
    !  This subroutine plots the structure in its undeformed and
    !  deformed configurations using the pgplot graphics
    !  subroutine library (http://www.astro.caltech.edu/~tjp/pgplot/).
    !  Note that the library only accepts single precision real
    !  variables.

    real ::         xmin0, xmax0, ymin0, ymax0, lx0, ly0, lmax0
    real ::         xmin, xmax, ymin, ymax, lx, ly, lmax
    real ::         px(5), py(5), maxu
    integer ::      ier, pgopen

    integer :: e, i, nen, num

    if (plot2screen) then
       !      ier = pgopen(' ')
       ier = pgopen('/xwin ')
    else
       ier = pgopen(trim(filename)//'_01.ps/cps')
       write (*, *)
       write (*, '("Writing pgplot output file to ", a)') &
            trim(filename)//'_01.ps'
       !write (*, '("(To view, type:  ghostview ", a, " & <enter>)")') &
       !   trim(filename)//'_01.ps'
       write (*, '("(To view: open the file ", a, " with GSView)")') &
            trim(filename)//'_01.ps'
    endif
    if (ier .le. 0) then
       write (*, *) 'Error: pgplot pgopen.'
       stop
    end if

    call pgpap(3.,2.)

    call pgsubp(1,2)

    call pgpage
    call pgsvp(0.,.96,0.,.96)

    ! Find maximum size of UN-deformed structure 
    xmin0 = 1d10
    xmax0 = -1d10
    ymin0 = 1d10
    ymax0 = -1d10
    do e = 1,ne 
       nen = element(e)%numnode
       xmin0 = min(real(minval(x(element(e)%ix(1:nen),1))),xmin0)  
       xmax0 = max(real(maxval(x(element(e)%ix(1:nen),1))),xmax0)  
       ymin0 = min(real(minval(x(element(e)%ix(1:nen),2))),ymin0)
       ymax0 = max(real(maxval(x(element(e)%ix(1:nen),2))),ymax0)
    end do
    lx0 = xmax0 - xmin0
    ly0 = ymax0 - ymin0
    lmax0 = max(lx0, ly0)

    ! Find maximum size of deformed structure 
    xmin = 1d10
    xmax = -1d10
    ymin = 1d10
    ymax = -1d10
    do e = 1,ne 
       nen = element(e)%numnode
       xmin = min(real(minval(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))),xmin)  
       xmax = max(real(maxval(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))),xmax)  
       ymin = min(real(minval(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))),ymin)
       ymax = max(real(maxval(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))),ymax)
    end do
    lx = xmax - xmin
    ly = ymax - ymin
    lmax = max(lx, ly)

    call pgslct(1)

    call pgsclp(0)
    call pgbbuf

    ! Plotting area 1, panel 1:  structure 
    call pgslct(1)
    call pgpanl(1,1)

    call pgswin(-0.1*lmax0+xmin0, 1.1*lmax0+xmin0, -0.1*lmax0+ymin0, 1.1*lmax0+ymin0)

    call pgeras

    call pgsch(1.8)
    call pgsci(1)

    call pgptxt(.5*lmax0+xmin0, 1.07*lmax0+ymin0, 0.0, 0.5, 'Structure')

    do e = 1,ne
       call pgsci(7)
       nen = element(e)%numnode
       px(1:nen) = real(x(element(e)%ix(1:nen),1))
       py(1:nen) = real(x(element(e)%ix(1:nen),2))
       call pgpoly(nen,px,py)
       call pgsci(4)
       px(nen+1) = real(x(element(e)%ix(1),1))
       py(nen+1) = real(x(element(e)%ix(1),2))
       call pgline(nen+1,px,py)
    end do

    ! Plotting area 1, panel 2:  Structure with displacements
    call pgslct(1)
    call pgpanl(1,2)
    call pgeras

    call pgswin(-0.1*lmax+xmin, 1.1*lmax+xmin, -0.1*lmax+ymin, 1.1*lmax+ymin)

    call pgsch(1.8)
    call pgsci(1)
    call pgptxt(.5*lmax+xmin, 1.07*lmax+ymin, 0.0, 0.5, 'Deformed structure')

    maxu = 0.0d0
    do i = 1,neqn
       maxu = max(abs(real(d(i))),maxu)
    end do

    if (scale_def /= 1.) then
       write (*, *)
       write (*, '("Note: plot scale factor = ", f16.9)') scale_def
    endif

    do e = 1,ne
       call pgsci(7)
       nen = element(e)%numnode
       ! note:  the plot scale parameter is defined in the fedata module.
       px(1:nen) = real(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))
       py(1:nen) = real(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))
       call pgpoly(nen,px,py)
       call pgsci(4)
       px(nen+1) = real(x(element(e)%ix(1),1) + scale_def * d(2*element(e)%ix(1)-1))
       py(nen+1) = real(x(element(e)%ix(1),2) + scale_def * d(2*element(e)%ix(1)))
       call pgline(nen+1,px,py)
    end do

    call pgiden

    call pgebuf

    if (.not. plot2screen) then
       call pgclos
    endif

  end subroutine plotdef
!!$
!!$subroutine ploteval(eval)
!!$!---------------------------------------------------------------------
!!$!  This subroutine plots element values of the structure in its 
!!$!  undeformed configuration using the pgplot graphics 
!!$!  subroutine library (http://www.astro.caltech.edu/~tjp/pgplot/).
!!$!  Note that the library only accepts single precision real
!!$!  variables.
!!$
!!$   real(8), dimension(:), intent(in) :: eval
!!$
!!$   real :: xmin, xmax, ymin, ymax, lx, ly, lmax
!!$   real :: px(5), py(5), maxu
!!$   real :: mineval, maxeval  
!!$   integer :: ier, pgopen, cid, ncid
!!$
!!$   integer :: i, e, nen, num, mm, pp, nc
!!$   character(len = 9) :: strng
!!$
!!$   if (plot2screen) then
!!$!      ier = pgopen(' ')
!!$      ier = pgopen('/SS')
!!$   else
!!$      ier = pgopen(trim(filename)//'_02.ps/cps')
!!$      write (*, *)
!!$      write (*, '("Writing pgplot output file to ", a)') &
!!$      trim(filename)//'_02.ps'
!!$      !write (*, '("(To view, type:  ghostview ", a, " & <enter>)")') &
!!$      !trim(filename)//'_02.ps'
!!$      write (*, '("(To view: open the file ", a, " with GSView)")') &
!!$         trim(filename)//'_02.ps'      
!!$   endif
!!$   if (ier .le. 0) then
!!$      write (*, *) 'Error: pgplot pgopen.'
!!$      stop
!!$   end if
!!$
!!$   ! Set color definitions      
!!$!   call pgscr(27, 1., 0., 0.)   ! 2
!!$!   call pgscr(26, 1., 0.5, 0.)  ! 8
!!$!   call pgscr(25, 1., 1., 0.)   ! 7
!!$!   call pgscr(24, 0.5, 1., 0.)  ! 9
!!$!   call pgscr(23, 0., 1., 0.5)  ! 10
!!$!   call pgscr(22, 0., 1., 0.)   ! 3
!!$!   call pgscr(21, 0., 1., 1.)   ! 5
!!$!   call pgscr(20, 0., 0.5, 1.)  ! 11
!!$!   call pgscr(19, 0., 0., 1.)   ! 4
!!$!   call pgscr(18, 0.5, 0., 1.)  ! 12
!!$!   call pgscr(17, 1., 0., 1.)   ! 6
!!$!   call pgscr(16, 1., 0., 0.5)  ! 13
!!$   call pgscr(27, 1., 0., 0.)   ! 2
!!$   call pgscr(26, 1., .3, 0.)  ! 8
!!$   call pgscr(25, 1., .6, 0.)   ! 7
!!$   call pgscr(24, 1., 1., 0.)  ! 9
!!$   call pgscr(23, .5, 1., 0.)  ! 10
!!$   call pgscr(22, 0., 1., 0.)   ! 3
!!$   call pgscr(21, 0., 1., .0)   ! 5
!!$   call pgscr(20, 0., .8, .2)  ! 11
!!$   call pgscr(19, 0., .6, .4)   ! 4
!!$   call pgscr(18, 0., .4, .6)  ! 12
!!$   call pgscr(17, 0., .2, .8)   ! 6
!!$   call pgscr(16, 0 , .0, 1.)  ! 13
!!$   ncid = 27-16+1
!!$
!!$   call pgpap(5.,1.)
!!$
!!$   call pgsvp(0.,.96,0.,.96)
!!$
!!$   ! Find maximum size of structure 
!!$   xmin = 1d10
!!$   xmax = -1d10
!!$   ymin = 1d10
!!$   ymax = -1d10
!!$   do e = 1,ne 
!!$      nen = element(e)%numnode
!!$      xmin = min(real(minval(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))),xmin)  
!!$      xmax = max(real(maxval(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))),xmax)  
!!$      ymin = min(real(minval(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))),ymin)
!!$      ymax = max(real(maxval(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))),ymax)
!!$   end do
!!$   lx = xmax - xmin
!!$   ly = ymax - ymin
!!$   lmax = max(lx, ly)
!!$
!!$   call pgswin(-0.1*lmax+xmin, 2.0*lmax+xmin, -0.1*lmax+ymin, 2.0*lmax+ymin)
!!$
!!$   call pgsclp(0)
!!$   call pgbbuf
!!$
!!$   call pgeras
!!$
!!$   call pgsch(1.8)
!!$   call pgsci(1)
!!$
!!$   call pgptxt(.5*lmax+xmin, 1.8*lmax+xmin, 0.0, 0.5, 'Element values')
!!$
!!$   mineval = minval(eval)
!!$   maxeval = maxval(eval) 
!!$   do e = 1,ne
!!$!      if (maxeval == mineval) then
!!$      if (abs(maxeval-mineval)<1e-10) then ! Sigmund 2007 change
!!$         cid = 16
!!$      else
!!$         cid = 16 + int(ncid * (eval(e)-mineval)/(maxeval-mineval))
!!$      endif
!!$      if (cid < 16) then
!!$         cid = 16
!!$      elseif (cid > 27) then
!!$         cid = 27
!!$      endif
!!$      call pgsci(cid)
!!$      
!!$      nen = element(e)%numnode
!!$      px(1:nen) = real(x(element(e)%ix(1:nen),1))
!!$      py(1:nen) = real(x(element(e)%ix(1:nen),2))
!!$      call pgpoly(nen,px,py)
!!$      call pgsci(1)
!!$      px(nen+1) = real(x(element(e)%ix(1),1))
!!$      py(nen+1) = real(x(element(e)%ix(1),2))
!!$      call pgline(nen+1,px,py)
!!$   end do
!!$
!!$   ! Plot legend
!!$   do i = 1, ncid
!!$      px(1) = 1.80*lmax+xmin
!!$      py(1) = lmax/ncid*(i-1)+ymin
!!$      px(2) = 1.80*lmax+lmax/ncid+xmin
!!$      py(2) = lmax/ncid*(i-1)+ymin
!!$      px(3) = 1.80*lmax+lmax/ncid+xmin
!!$      py(3) = lmax/ncid*i+ymin
!!$      px(4) = 1.80*lmax+xmin
!!$      py(4) = lmax/ncid*i+ymin
!!$      px(5) = px(1)
!!$      py(5) = py(1)
!!$      call pgsci(16-1+i)
!!$      call pgpoly(4,px,py)
!!$      call pgsci(1)
!!$      call pgline(nen+1,px,py)
!!$      write (strng, '(g9.3)') mineval+(maxeval-mineval)/ncid*(i-1)
!!$      call pgptxt(1.75*lmax+xmin, lmax/12.*(i-1)-lmax/ncid/2.+ymin, 0., 1.0, strng)
!!$   end do
!!$   write (strng, '(g9.3)') maxeval
!!$   call pgptxt(1.75*lmax+xmin, lmax/12.*ncid-lmax/ncid/2.+ymin, 0., 1.0, strng)
!!$
!!$   call pgiden
!!$
!!$   call pgebuf
!!$
!!$   if (.not. plot2screen) then
!!$      call pgclos
!!$   endif
!!$
!!$end subroutine ploteval
!!$ 
!!$subroutine plotevec(p1, p2, ang)
!!$!---------------------------------------------------------------------
!!$! Plot vector field 
!!$
!!$   real(8), dimension(:), intent(in) :: p1, p2, ang
!!$
!!$   real :: xmin, xmax, ymin, ymax, lx, ly, lmax, clx, cly, clmax
!!$   real :: px(5), py(5), maxu
!!$   real :: minep1, maxep1, avepx, avepy, ca, sa, sf
!!$   integer :: ier, pgopen, cid, ncid
!!$   integer :: i, j, e, nen
!!$
!!$   if (plot2screen) then
!!$!      ier = pgopen(' ')
!!$      ier = pgopen('/SS')
!!$   else
!!$      ier = pgopen(trim(filename)//'_03.ps/cps')
!!$      write (*, *)
!!$      write (*, '("Writing pgplot output file to ", a)') &
!!$      trim(filename)//'_03.ps'
!!$!      write (*, '("(To view, type:  ghostview ", a, " & <enter>)")') &
!!$!      trim(filename)//'_03.ps'
!!$      write (*, '("(To view: open the file ", a, " with GSView)")') &
!!$         trim(filename)//'_03.ps'      
!!$   endif
!!$   if (ier .le. 0) then
!!$      write (*, *) 'Error: pgplot pgopen.'
!!$      stop
!!$   end if
!!$
!!$   call pgpap(5.,1.)
!!$
!!$   call pgsvp(0.,.96,0.,.96)
!!$
!!$   ! Find maximum size of structure 
!!$   xmin = 1d10
!!$   xmax = -1d10
!!$   ymin = 1d10
!!$   ymax = -1d10
!!$   clx = 0d0
!!$   cly = 0d0
!!$   do e = 1,ne 
!!$      nen = element(e)%numnode
!!$      xmin = min(real(minval(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))),xmin)  
!!$      xmax = max(real(maxval(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))),xmax)  
!!$      ymin = min(real(minval(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))),ymin)
!!$      ymax = max(real(maxval(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))),ymax)
!!$     
!!$      ! Define characteristic length
!!$      do i = 1, nen
!!$         do j = i+1, nen
!!$            clx = max(real(abs(x(element(e)%ix(i),1)-x(element(e)%ix(j),1))),clx)
!!$            cly = max(real(abs(x(element(e)%ix(i),2)-x(element(e)%ix(j),2))),cly)
!!$         end do
!!$      end do
!!$   end do
!!$   lx = xmax - xmin
!!$   ly = ymax - ymin
!!$   lmax = max(lx, ly)        
!!$   clmax = max(clx, cly)
!!$
!!$   call pgswin(-0.1*lmax+xmin, 1.1*lmax+xmin, -0.1*lmax+ymin, 1.1*lmax+ymin)
!!$
!!$   call pgsclp(0)
!!$   call pgbbuf
!!$
!!$   call pgeras
!!$
!!$   call pgsch(1.8)
!!$   call pgsci(1)
!!$
!!$   call pgptxt(.5*lmax+xmin, 1.1*lmax+ymin, 0.0, 0.5, 'Element vector plot')
!!$
!!$   if (scale_vec /= 1.) then
!!$      write (*, *)
!!$      write (*, '("Note: plot scale factor = ", f16.9)') scale_vec
!!$   endif
!!$
!!$   maxep1 = max(maxval(abs(p1)),maxval(abs(p2)))
!!$   do e = 1,ne
!!$      nen = element(e)%numnode
!!$      px(1:nen) = real(x(element(e)%ix(1:nen),1))
!!$      py(1:nen) = real(x(element(e)%ix(1:nen),2))
!!$
!!$      call pgsci(1)
!!$      px(nen+1) = real(x(element(e)%ix(1),1))
!!$      py(nen+1) = real(x(element(e)%ix(1),2))
!!$      call pgline(nen+1,px,py)
!!$
!!$      sf = maxep1*sqrt(2.)/clmax / scale_vec
!!$      avepx = sum(x(element(e)%ix(1:nen),1))/nen
!!$      avepy = sum(x(element(e)%ix(1:nen),2))/nen
!!$
!!$      ca = cos(-ang(e))
!!$      sa = sin(-ang(e))
!!$      px(1) = avepx - 0.5 * abs(p1(e)) * ca/sf
!!$      py(1) = avepy - 0.5 * abs(p1(e)) * sa/sf
!!$      px(2) = avepx + 0.5 * abs(p1(e)) * ca/sf
!!$      py(2) = avepy + 0.5 * abs(p1(e)) * sa/sf
!!$      if (p1(e) >= 0.) then
!!$         call pgsci(4)
!!$      else
!!$         call pgsci(2)
!!$      endif
!!$      call pgline(2,px,py)
!!$
!!$      ca = cos(-ang(e)+pi/2.)
!!$      sa = sin(-ang(e)+pi/2.)
!!$      px(1) = avepx - 0.5 * abs(p2(e)) * ca/sf
!!$      py(1) = avepy - 0.5 * abs(p2(e)) * sa/sf
!!$      px(2) = avepx + 0.5 * abs(p2(e)) * ca/sf
!!$      py(2) = avepy + 0.5 * abs(p2(e)) * sa/sf
!!$      if (p2(e) >= 0.) then
!!$         call pgsci(4)
!!$      else
!!$         call pgsci(2)
!!$      endif
!!$      call pgline(2,px,py)
!!$
!!$   end do
!!$
!!$   call pgiden
!!$
!!$   call pgebuf
!!$   
!!$   if (.not. plot2screen) then
!!$     call pgclos
!!$   endif
!!$ 
!!$end subroutine plotevec
!!$
!!$subroutine ploteig(eval, evec, timef, timei)
!!$!---------------------------------------------------------------------
!!$! This subroutine animates the deformation modes of the structure
!!$! using the pgplot graphics subroutine library
!!$! (http://www.astro.caltech.edu/~tjp/pgplot/).
!!$! Note that the library only accepts single precision real variables.
!!$
!!$   real(8), intent(in) :: eval, evec(:), timef, timei
!!$
!!$   real :: xmin, xmax, ymin, ymax, lx, ly, lmax
!!$   real :: px(5), py(5)
!!$   real :: mineval, maxeval
!!$   integer :: ier, pgopen
!!$
!!$   integer :: i, e, nen
!!$   character(len = 9) :: strng
!!$   real :: scale, time
!!$
!!$   if (plot2screen) then
!!$!      ier = pgopen(' ')
!!$      ier = pgopen('/SS')
!!$   else
!!$      ier = pgopen(trim(filename)//'_04.ps/cps')
!!$      write (*, *)
!!$      write (*, '("Writing pgplot output file to ", a)') &
!!$      trim(filename)//'_04.ps'
!!$      !write (*, '("(To view, type:  ghostview ", a, " & <enter>)")') &
!!$      !trim(filename)//'_04.ps'
!!$      write (*, '("(To view: open the file ", a, " with GSView)")') &
!!$         trim(filename)//'_04.ps'      
!!$   endif
!!$   if (ier .le. 0) then
!!$      write (*, *) 'Error: pgplot pgopen.'
!!$      stop
!!$   end if
!!$
!!$   ! set color definitions
!!$!   call pgscr(27, 1., 0., 0.)   ! 2
!!$!   call pgscr(26, 1., 0.5, 0.)  ! 8
!!$!   call pgscr(25, 1., 1., 0.)   ! 7
!!$!   call pgscr(24, 0.5, 1., 0.)  ! 9
!!$!   call pgscr(23, 0., 1., 0.5)  ! 10
!!$!   call pgscr(22, 0., 1., 0.)   ! 3
!!$!   call pgscr(21, 0., 1., 1.)   ! 5
!!$!   call pgscr(20, 0., 0.5, 1.)  ! 11
!!$!   call pgscr(19, 0., 0., 1.)   ! 4
!!$!   call pgscr(18, 0.5, 0., 1.)  ! 12
!!$!   call pgscr(17, 1., 0., 1.)   ! 6
!!$!   call pgscr(16, 1., 0., 0.5)  ! 13
!!$   call pgscr (1, 1.00, 1.00, 1.00)
!!$   call pgscr (0, 0.00, 0.00, 0.00)
!!$   call pgscr(27, 1., 0., 0.)   ! 2
!!$   call pgscr(26, 1., .3, 0.)  ! 8
!!$   call pgscr(25, 1., .6, 0.)   ! 7
!!$   call pgscr(24, 1., 1., 0.)  ! 9
!!$   call pgscr(23, .5, 1., 0.)  ! 10
!!$   call pgscr(22, 0., 1., 0.)   ! 3
!!$   call pgscr(21, 0., 1., .0)   ! 5
!!$   call pgscr(20, 0., .8, .2)  ! 11
!!$   call pgscr(19, 0., .6, .4)   ! 4
!!$   call pgscr(18, 0., .4, .6)  ! 12
!!$   call pgscr(17, 0., .2, .8)   ! 6
!!$   call pgscr(16, 0 , .0, 1.)  ! 13
!!$
!!$   call pgpap(5.,1.)
!!$
!!$   call pgsvp(0.,.96,0.,.96)
!!$
!!$   ! Find maximum size of structure
!!$   xmin = 1d10
!!$   xmax = -1d10
!!$   ymin = 1d10
!!$   ymax = -1d10
!!$   do e = 1,ne
!!$      nen = element(e)%numnode
!!$      xmin = min(real(minval(x(element(e)%ix(1:nen),1) + scale_def * evec(2*element(e)%ix(1:nen)-1))),xmin)
!!$      xmax = max(real(maxval(x(element(e)%ix(1:nen),1) + scale_def * evec(2*element(e)%ix(1:nen)-1))),xmax)
!!$      ymin = min(real(minval(x(element(e)%ix(1:nen),2) + scale_def * evec(2*element(e)%ix(1:nen)))),ymin)
!!$      ymax = max(real(maxval(x(element(e)%ix(1:nen),2) + scale_def * evec(2*element(e)%ix(1:nen)))),ymax)
!!$   end do
!!$   lx = xmax - xmin
!!$   ly = ymax - ymin
!!$   lmax = max(lx, ly)
!!$
!!$   call pgswin(-0.1*lmax+xmin, 1.1*lmax+xmin, -0.1*lmax+ymin, 1.1*lmax+ymin)
!!$
!!$   call pgsclp(0)
!!$
!!$   time = 0.
!!$   do
!!$      scale = sin(sqrt(eval)*time)
!!$
!!$      call pgbbuf
!!$      call pgeras
!!$
!!$      call pgsch(1.8)
!!$      call pgsci(1)
!!$      call pgptxt(.5*lmax+xmin, 1.07*lmax+ymin, 0.0, 0.5, 'Mode shape')
!!$
!!$      do e = 1,ne
!!$         call pgsci(7)
!!$         nen = element(e)%numnode
!!$         ! note:  the plot scale parameter is defined in the fedata module.
!!$         px(1:nen) = real(x(element(e)%ix(1:nen),1) + scale * scale_def * evec(2*element(e)%ix(1:nen)-1))
!!$         py(1:nen) = real(x(element(e)%ix(1:nen),2) + scale * scale_def * evec(2*element(e)%ix(1:nen)))
!!$         call pgpoly(nen,px,py)
!!$         call pgsci(4)
!!$         px(nen+1) = real(x(element(e)%ix(1),1) + scale * scale_def * evec(2*element(e)%ix(1)-1))
!!$         py(nen+1) = real(x(element(e)%ix(1),2) + scale * scale_def * evec(2*element(e)%ix(1)))
!!$         call pgline(nen+1,px,py)
!!$      end do
!!$
!!$      call pgiden
!!$
!!$      call pgebuf
!!$
!!$      time = time + timei
!!$      if (time > timef) exit
!!$
!!$   end do
!!$
!!$   if (.not. plot2screen) then
!!$      call pgclos
!!$   endif
!!$
!!$end subroutine ploteig


  subroutine plotanim(iter, device, window, color, deformed, vector, elements, &
       legend, title, maxp, p1, p2, ang, mineval, maxeval, eval)
    !-------------------------------------------------------------------------------
    ! This subroutine plots and animates  vector and/or element values in color 
    ! or monochrome with or without a legend.
    !  
    ! input: 
    ! iter     integer    iteration number: iter =  0    plot once to window
    !                                            >= 1    animate (after final plot, 
    !                                                    call plotanim with iter = -1)
    !                                            = -1    close animation window
    ! device   integer    device number: device = 0      plot to default device (/xwin)
    !                                           = 1      plot to postscript file
    !                                           = 2      plot to gif file
    !                                           = 3      plot to user device
    !                                           = 4      plot to animation
    ! window   integer    window number: 1 <= window <= 8
    ! color    logical    plot color: color = .true.     color    
    !                                       = .false.    grayscale
    !                     (for link1 elements, .true. plots values, .false. plots areas)
    ! deformed logical    plot deformed structure
    ! vector   logical    plot vector values (supply maxp, p1, p2, and ang)
    ! elements logical    plot element values (supply mineval, maxeval, eval)
    ! legend   logical    plot legend
    ! title    character  plot name
    ! maxp     real(8)  maximum p1 and p2 value for scaling
    ! p1(:)    real(8)  first principal values
    ! p2(:)    real(8)  second principal values
    ! ang(:)   real(8)  orientation
    ! mineval  real(8)  minimum eval value
    ! maxeval  real(8)  maximum eval value
    ! eval(:)  real(8)  element values   

    use fedata

    integer, intent(in) :: iter, device, window
    logical, intent(in) :: color, deformed, vector, elements, legend
    real(8), intent(in) :: maxp, mineval, maxeval
    real(8), dimension(:), intent(in) :: p1, p2, ang, eval
    character(len = *), intent(in) :: title

    real(4) :: xmin, xmax, ymin, ymax, lx, ly, lmax(8), clx, cly, clmax(8)
    real(4) :: px(5), py(5), maxu
    real(4) :: maxep1, avepx, avepy, ca, sa, sf
    integer :: ier, pgopen, cid, ncid, ncidmin, ncidmax
    integer :: i, j, e, nen
    integer :: panel(8)
    integer :: strnglen
    character(len = 12) :: strng
    character(len = 16) :: varp

    save panel, lmax, clmax, xmin, xmax, ymin, ymax

    if (window < 1 .or. window > 8) then
       write (*, *)
       write (*, '("Error: pgplot windows 1-8 available only.")')
       stop
    endif

    if (iter == 0 .or. iter == 1 .or. device == 4) then
       if (device == 0) then
          ! ændret af paw. Sat pgopen til altid at åbne xwin når device==0. Evt kan der skrives pgopen(''), så vil den åbne det der er specificeret med PGPLOT_DEV (ved compilation time)
          !se  http://www.astro.caltech.edu/~tjp/pgplot/subroutines.html#PGOPEN
          panel(window) = pgopen('/xwin ')
          !         panel(window) = pgopen('/SS')
       elseif (device == 1) then
          panel(window) = pgopen(trim(filename)//'.ps/cps')
       elseif (device == 2) then
          panel(window) = pgopen(trim(filename)//'.gif/gif')
       elseif (device == 3) then
          panel(window) = pgopen('?')
       elseif (device == 4) then
          varp(1:4) = 'film'
          write(varp(5:8),'(i4)') iter
          if (iter.lt.10) then 
             varp(5:7) = '000'
          endif
          if ((iter.ge.10).and.(iter.lt.100)) then 
             varp(5:6) = '00'
          endif
          if ((iter.ge.100).and.(iter.lt.1000)) then 
             varp(5:5) = '0'
          endif
          varp(9:16)='.gif/GIF' 	   
          panel(window) = pgopen(varp)
       else
          write (*, *) 
          write (*, '("Error: invalid pgplot device")')
       endif

       if (panel(window) .le. 0) then
          write (*, *) 'Error: pgplot pgopen.'
          stop
       end if

       ! Find maximum size of deformed structure
       xmin = 1e10
       xmax = -1e10
       ymin = 1e10
       ymax = -1e10
       clx = -1e10
       cly = -1e10
       do e = 1,ne
          nen = element(e)%numnode
          xmin = min(real(minval(x(element(e)%ix(1:nen),1))),xmin)
          xmax = max(real(maxval(x(element(e)%ix(1:nen),1))),xmax)
          ymin = min(real(minval(x(element(e)%ix(1:nen),2))),ymin)
          ymax = max(real(maxval(x(element(e)%ix(1:nen),2))),ymax)
          do i = 1, nen
             do j = i+1, nen
                clx = max(real(abs(x(element(e)%ix(i),1)-x(element(e)%ix(j),1))),clx)
                cly = max(real(abs(x(element(e)%ix(i),2)-x(element(e)%ix(j),2))),cly)
             end do
          end do
       end do
       ! Add space around structure
       ! xmin = xmin - 0.0*(xmax-xmin)
       ! xmax = xmax + 0.0*(xmax-xmin)
       ! ymin = ymin - 0.0*(ymax-ymin)
       ! ymax = ymax - 0.0*(ymax-ymin)

       xmin = xmin - 0.0*(xmax-xmin)
       xmax = xmax + 0.0*(xmax-xmin)
       ymin = ymin - 0.0*(ymax-ymin)
       ymax = ymax - 0.0*(ymax-ymin)


       lx = xmax - xmin
       ly = ymax - ymin
       lmax(window) = max(lx, ly)
       clmax(window) = max(clx, cly)

       call pgpap(7.,1.)
       call pgsvp(0.,.96,0.,.96)
       if (legend) then
          call pgswin(-0.1*lmax(window)+xmin, 1.9*lmax(window)+xmin, &
               -0.1*lmax(window)+ymin, 1.9*lmax(window)+ymin)
       else
          call pgswin(-0.1*lmax(window)+xmin, 1.1*lmax(window)+xmin, &
               -0.1*lmax(window)+ymin, 1.1*lmax(window)+ymin)
       endif
       call pgsclp(0)

    endif
    if (iter < 0) then
       call pgclos
       return
    endif

    call pgslct(panel(window))

    ! Set color definitions
    if (color) then
       !      call pgscr (1, 1.00, 1.00, 1.00)
       !      call pgscr (0, 0.00, 0.00, 0.00)
       !      call pgscr(27, 1., 0., 0.)   ! 2
       !      call pgscr(26, 1., 0.5, 0.)  ! 8
       !      call pgscr(25, 1., 1., 0.)   ! 7
       !      call pgscr(24, 0.5, 1., 0.)  ! 9
       !      call pgscr(23, 0., 1., 0.5)  ! 10
       !      call pgscr(22, 0., 1., 0.)   ! 3
       !      call pgscr(21, 0., 1., 1.)   ! 5
       !      call pgscr(20, 0., 0.5, 1.)  ! 11
       !      call pgscr(19, 0., 0., 1.)   ! 4
       !      call pgscr(18, 0.5, 0., 1.)  ! 12
       !      call pgscr(17, 1., 0., 1.)   ! 6
       !      call pgscr(16, 1., 0., 0.5)  ! 13

       call pgscr (1, 1.00, 1.00, 1.00)
       call pgscr (0, 0.00, 0.00, 0.00)
       call pgscr(27, 1., 0., 0.)   ! 2
       call pgscr(26, 1., .3, 0.)  ! 8
       call pgscr(25, 1., .6, 0.)   ! 7
       call pgscr(24, 1., 1., 0.)  ! 9
       call pgscr(23, .5, 1., 0.)  ! 10
       call pgscr(22, 0., 1., 0.)   ! 3
       call pgscr(21, 0., 1., .0)   ! 5
       call pgscr(20, 0., .8, .2)  ! 11
       call pgscr(19, 0., .6, .4)   ! 4
       call pgscr(18, 0., .4, .6)  ! 12
       call pgscr(17, 0., .2, .8)   ! 6
       call pgscr(16, 0 , .0, 1.)  ! 13
       ncidmin = 16
       ncidmax = 27
    else
       ! Set color definitions
       call pgscr (0, 1.00, 1.00, 1.00)
       call pgscr (1, 0.00, 0.00, 0.00)
       call pgscr (16, 1.00, 1.00, 1.00)
       call pgscr (17, 0.90, 0.90, 0.90)
       call pgscr (18, 0.80, 0.80, 0.80)
       call pgscr (19, 0.70, 0.70, 0.70)
       call pgscr (20, 0.60, 0.60, 0.60)
       call pgscr (21, 0.50, 0.50, 0.50)
       call pgscr (22, 0.40, 0.40, 0.40)
       call pgscr (23, 0.30, 0.30, 0.30)
       call pgscr (24, 0.20, 0.20, 0.20)  
       call pgscr (25, 0.10, 0.10, 0.10)
       call pgscr (26, 0.00, 0.00, 0.00)
       ncidmin = 16
       ncidmax = 26
    endif
    ncid = ncidmax-ncidmin+1

    call pgbbuf

    call pgeras

    call pgsch(1.8)
    call pgsci(1)

    ! Write title
    if (legend) then
       call pgsci(1)
       call pgptxt(.5*lmax(window)+xmin, 1.15*lmax(window)+ymin, 0.0, 0.5, title)
    end if

    ! Write iteration
    if (iter >= 1) then
       if (legend) then
          call pgsci(1)
          write (strng, '(i9)') iter
          strng = adjustl(strng)
          strnglen =  len_trim(strng)
          call pgptxt(.5*lmax(window)+xmin, 1.05*lmax(window)+ymin, 0.0, 0.5, 'Iteration #'//strng(1:strnglen))
       end if
    endif

    do e = 1,ne
       nen = element(e)%numnode
       if (nen == 2) then
          if (color) then
             if (abs(maxeval-mineval)<1e-10) then ! Sigmund 2007 change
                !            if (maxeval == mineval) then
                cid = ncidmin
             else
                cid = ncidmin + int(ncid * (eval(e)-mineval)/(maxeval-mineval))
             endif
             if (cid < ncidmin) then
                cid = ncidmin
             elseif (cid > ncidmax) then
                cid = ncidmax
             endif
             call pgsci(cid)
             if (deformed) then
                px(1:nen) = real(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))
                py(1:nen) = real(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))
             else
                px(1:nen) = real(x(element(e)%ix(1:nen),1))
                py(1:nen) = real(x(element(e)%ix(1:nen),2))
             endif
             call pgline(nen,px,py)
          else
             call pgsci(ncidmax)
             if (int(10*eval(e)/maxeval) > 0) then
                if (deformed) then
                   px(1:nen) = real(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))
                   py(1:nen) = real(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))
                else
                   px(1:nen) = real(x(element(e)%ix(1:nen),1))
                   py(1:nen) = real(x(element(e)%ix(1:nen),2))
                endif
                call pgslw(int(10*eval(e)/maxeval))
                call pgline(nen,px,py)
                call pgslw(1)
             endif
          endif
       else
          if (elements) then
             !            if (maxeval == mineval) then
             if (abs(maxeval-mineval)<1e-10) then ! Sigmund 2007 change
                cid = ncidmin
             else
                cid = ncidmin + int(ncid * (eval(e)-mineval)/(maxeval-mineval))
             endif
             if (cid < ncidmin) then
                cid = ncidmin
             elseif (cid > ncidmax) then
                cid = ncidmax
             endif
             call pgsci(cid)
             if (deformed) then
                px(1:nen) = real(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))
                py(1:nen) = real(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))
             else
                px(1:nen) = real(x(element(e)%ix(1:nen),1))
                py(1:nen) = real(x(element(e)%ix(1:nen),2))
             endif
             call pgpoly(nen,px,py)
          endif
          call pgsci(1)
          ! SIGMUNDS ADDITIONS!
          if (deformed) then
             px(1:nen) = real(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))
             py(1:nen) = real(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))
          else
             px(1:nen) = real(x(element(e)%ix(1:nen),1))
             py(1:nen) = real(x(element(e)%ix(1:nen),2))
          endif
          ! END SIGMUNDS ADDITIONS!
          if (deformed) then
             px(nen+1) = real(x(element(e)%ix(1),1) + scale_def * d(2*element(e)%ix(1)-1))
             py(nen+1) = real(x(element(e)%ix(1),2) + scale_def * d(2*element(e)%ix(1)))
          else
             px(nen+1) = real(x(element(e)%ix(1),1))
             py(nen+1) = real(x(element(e)%ix(1),2))
          endif
          call pgline(nen+1,px,py)
       endif

       if (vector) then
          sf = maxp*sqrt(2.)/clmax(window) / scale_vec
          if (deformed) then
             avepx = sum(x(element(e)%ix(1:nen),1) + scale_def * d(2*element(e)%ix(1:nen)-1))/nen
             avepy = sum(x(element(e)%ix(1:nen),2) + scale_def * d(2*element(e)%ix(1:nen)))/nen
          else
             avepx = sum(x(element(e)%ix(1:nen),1))/nen
             avepy = sum(x(element(e)%ix(1:nen),2))/nen
          endif

          ca = cos(-ang(e))
          sa = sin(-ang(e))
          px(1) = avepx - 0.5 * abs(p1(e)) * ca/sf
          py(1) = avepy - 0.5 * abs(p1(e)) * sa/sf
          px(2) = avepx + 0.5 * abs(p1(e)) * ca/sf
          py(2) = avepy + 0.5 * abs(p1(e)) * sa/sf
          if (p1(e) >= 0.) then
             call pgsci(4)
          else
             call pgsci(2)
          endif
          call pgline(2,px,py)

          ca = cos(-ang(e)+pi/2.)
          sa = sin(-ang(e)+pi/2.)
          px(1) = avepx - 0.5 * abs(p2(e)) * ca/sf
          py(1) = avepy - 0.5 * abs(p2(e)) * sa/sf
          px(2) = avepx + 0.5 * abs(p2(e)) * ca/sf
          py(2) = avepy + 0.5 * abs(p2(e)) * sa/sf
          if (p2(e) >= 0.) then
             call pgsci(4)
          else
             call pgsci(2)
          endif
          call pgline(2,px,py)
       endif

    end do

    ! Plot legend
    if (legend) then
       do i = 1, ncid
          px(1) = 1.80*lmax(window)+xmin
          py(1) = lmax(window)/ncid*(i-1)+ymin
          px(2) = 1.80*lmax(window)+lmax(window)/ncid+xmin
          py(2) = lmax(window)/ncid*(i-1)+ymin
          px(3) = 1.80*lmax(window)+lmax(window)/ncid+xmin
          py(3) = lmax(window)/ncid*i+ymin
          px(4) = 1.80*lmax(window)+xmin
          py(4) = lmax(window)/ncid*i+ymin
          px(5) = px(1)
          py(5) = py(1)
          call pgsci(ncidmin-1+i)
          call pgpoly(4,px,py)
          call pgsci(1)
          call pgline(5,px,py)
          write (strng, '(g9.3)') mineval+(maxeval-mineval)/ncid*(i-1)
          call pgptxt(1.75*lmax(window)+xmin, &
               lmax(window)/real(ncidmax-ncidmin+1)*(i-1)-lmax(window)/ncid/2.+ymin, 0., 1.0, strng)
       end do
       write (strng, '(g9.3)') maxeval
       call pgptxt(1.75*lmax(window)+xmin, &
            lmax(window)/real(ncidmax-ncidmin+1)*ncid-lmax(window)/ncid/2.+ymin, 0., 1.0, strng)
    endif

    call pgiden

    call pgebuf

    if (iter == 0 .or. device == 4) then
       call pgclos
    endif

  end subroutine plotanim

  integer function numnode(ename)
    !---------------------------------------------------------------------
    !  This function returns the number of nodes per element.
    !  Note that this function appears in the processor and
    !  fea modules.

    character(LEN = 12) :: ename

    if (ename == 'LINK1') THEN
       numnode = 2
    elseif (ename == 'PLANE42RECT') THEN
       numnode = 4
    elseif (ename == 'PLANE42') THEN
       numnode = 4
       ! SIGMUND CHANGE
    elseif (ename == 'PLATE41') THEN
       numnode = 4
    else
       write (*, *) 'ERROR: undefined element'
       stop
    endif

  end function numnode

  subroutine stopwatch(oper)
    !---------------------------------------------------------------------
    !   This subrotuine computes elapsed wallclock time after a second call.

    character(len = 4), intent(in) :: oper
    integer time_array_0(8), time_array_1(8)

    if (oper == 'STAR' .or. oper == 'star') then
       call date_and_time(values=time_array_0)
       start_time = time_array_0 (5) * 3600 + time_array_0 (6) * 60 &
            + time_array_0 (7) + 0.001 * time_array_0 (8)
    elseif (oper == 'STOP' .or. oper == 'stop') then
       call date_and_time(values=time_array_1)
       finish_time = time_array_1 (5) * 3600 + time_array_1 (6) * 60 &
            + time_array_1 (7) + 0.001 * time_array_1 (8)
       write (6, '(8x, 1a, 1f16.8)') 'elapsed wall clock time:', &
            finish_time - start_time 
    else
       write (*, *)
       write (*, '("Error: in Processor/stopwatch.")')
       stop
    endif
  end subroutine stopwatch

  subroutine plotmatlabdef(title)
    !---------------------------------------------------------------------
    ! Subroutine to plot the un/deformed structure using Matlab

    character(len=*), intent(in) :: title
    integer i, j, e

    ! write datafile        
    open (13, file = trim(filename)//'_plotdeformed_data.m')

    ! write nodal coordinates
    write(13,'("X = [")')
    do i = 1,size(x,1)
       write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
    enddo
    write(13,'("];")')
    write(13,'( )')

    ! write topology matrix
    write(13,'("IX = [")')
    do e = 1,ne
       write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
    enddo
    write(13,'("];")')

    ! write deformation vector
    write(13,'("D = [")')
    do i = 1,neqn
       write (13,'(f15.9,1x)') (d(i) )
    enddo
    write(13,'("];")')   
    close(13)

    ! Create matlab script
    open (13, file = trim(filename)//'_plotdeformed.m')
    write(13,*) '% Plotting Un-Deformed and Deformed Structure'
    write(13,*) 'close all'
    write(13,*) 'clear all'
    write(13,*) trim(filename), '_plotdeformed_data;'
    write(13,*) '% Make plot'
    write(13,*) 'figure'
    write(13,*) 'set(gcf,',"'",'Name',"','", trim(title) ,"'",')'

    ! Element dependent code
    ! NOTE: not possible to mix element types !!!
    if (element(1)%id .eq. 1) then
       write(13,*) 'subplot(2,1,2)'
       write(13,*) 'hold on'      
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '   edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];'
       write(13,*) '   xx = X(IX(e,1:2),1) + D(edof(1:2:4));'
       write(13,*) '   yy = X(IX(e,1:2),2) + D(edof(2:2:4));'
       write(13,*) '   plot(xx,yy,',"'",'b',"',","'",'LineWidth',"'",',1.5)'
       write(13,*) 'end'
       write(13,*) 'title(',"'",'Deformed',"'",')'            
       write(13,*) 'axis equal'
       write(13,*) 'xaxes = get(gca,',"'",'xlim',"'",');'
       write(13,*) 'yaxes = get(gca,',"'",'ylim',"'",');'
       write(13,*) 'axis off;'
       write(13,*) 'subplot(2,1,1)'
       write(13,*) 'hold on'
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '   xx = X(IX(e,1:2),1);'
       write(13,*) '   yy = X(IX(e,1:2),2);'
       write(13,*) '   plot(xx,yy,',"'",'b',"',","'",'LineWidth',"',",'1.5)'
       write(13,*) 'end'
       write(13,*) 'title(',"'",'Undeformed',"'",')'      
       write(13,*) 'axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...'
       write(13,*) ' min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);'
       write(13,*) 'axis off;'
    elseif (element(1)%id == 2) then
       write(13,*) 'subplot(2,1,2)'
       write(13,*) 'hold on'      
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)...'
       write(13,*) '        2*IX(e,3)-1 2*IX(e,3) 2*IX(e,4)-1 2*IX(e,4)];'
       write(13,*) '    xx = X(IX(e,1:4),1) + D(edof(1:2:8));'
       write(13,*) '    yy = X(IX(e,1:4),2) + D(edof(2:2:8));'
       write(13,*) '    patch(xx,yy,[1 1 0]);'
       write(13,*) 'end'
       write(13,*) 'title(',"'",'Deformed',"'",')'            
       write(13,*) 'axis equal;'
       write(13,*) 'xaxes = get(gca,',"'",'xlim',"'",');'
       write(13,*) 'yaxes = get(gca,',"'",'ylim',"'",');'
       write(13,*) 'axis off;'
       write(13,*) 'subplot(2,1,1)'
       write(13,*) 'hold on'      
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '    xx = X(IX(e,1:4),1);'
       write(13,*) '    yy = X(IX(e,1:4),2);'
       write(13,*) '    patch(xx,yy,[1 1 0]);'
       write(13,*) 'end'
       write(13,*) 'title(',"'",'Undeformed',"'",')'      
       write(13,*) 'axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...'
       write(13,*) ' min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);'
       write(13,*) 'axis equal;'
       write(13,*) 'axis off;'
       !ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ Tilfï¿½jet d 1.9-11 for shell41
    elseif (element(1)%id == 3) then
       write(13,*) 'subplot(2,1,2)'
       write(13,*) 'hold on'      
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '    edof = [3*IX(e,1)-2 3*IX(e,1)-1 3*IX(e,1)...'
       write(13,*) '            3*IX(e,2)-2 3*IX(e,2)-1 3*IX(e,2)...'
       write(13,*) '            3*IX(e,3)-2 3*IX(e,3)-1 3*IX(e,3)...'     
       write(13,*) '            3*IX(e,4)-2 3*IX(e,4)-1 3*IX(e,4)]; '
       write(13,*) '    xx = X(IX(e,1:4),1);'
       write(13,*) '    yy = X(IX(e,1:4),2);'
       write(13,*) '    zz = X(IX(e,1:4),3) + D(edof(1:3:10));'
       write(13,*) '    patch(xx,yy,zz,[1 1 0]);'
       write(13,*) 'end'
       write(13,*) 'title(',"'",'Deformed',"'",')'            
       write(13,*) 'axis equal;'
       write(13,*) 'xaxes = get(gca,',"'",'xlim',"'",');'
       write(13,*) 'yaxes = get(gca,',"'",'ylim',"'",');'
       write(13,*) 'axis off;'
       write(13,*) 'subplot(2,1,1)'
       write(13,*) 'hold on'      
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '    xx = X(IX(e,1:4),1);'
       write(13,*) '    yy = X(IX(e,1:4),2);'
       write(13,*) '    zz = X(IX(e,1:4),3);'
       write(13,*) '    patch(xx,yy,zz,[1 1 0]);'
       write(13,*) 'end'
       write(13,*) 'title(',"'",'Undeformed',"'",')'      
       write(13,*) 'axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...'
       write(13,*) ' min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);'
       write(13,*) 'axis equal;'
       write(13,*) 'axis off;'
    else
       write (*,'("Unsupported element type in Matlab routine:")')
       write (*,'("plotmatlab -> plotmatlabdef")')  
    endif
    ! End file and close
    write(13,*) 'hold off'
    write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
    close(13)

    ! Print to screen:
    print*,'Matlab plotting called: un/deformed'
    print*,'Generated data file: ', trim(filename)//'_plotdeformed_data.m' 
    print*,'and script file: ', trim(filename)//'_plotdeformed.m'
    print*,' '

  end subroutine plotmatlabdef

  subroutine plotmatlabeval(title,data1)
    ! Subroutine to plot the element values using Matlab

    character(len=*), intent(in) :: title
    real(8), dimension(:), intent(in) :: data1
    integer i, j, e

    ! write datafile     
    open (13, file = trim(filename)//'_plotelements_data.m')

    ! write nodal coordinates
    write(13,'("X = [")')
    do i = 1,size(x,1)
       write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
    enddo
    write(13,'("];")')
    write(13,'( )')

    ! write topology matrix
    write(13,'("IX = [")')
    do e = 1,ne
       write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
    enddo
    write(13,'("];")')

    ! write data1/element values
    write(13,'("plotval = [")')
    do i = 1,ne
       write (13,'(f32.15)') data1(i)
    enddo
    write(13,'("];")')   
    close(13) 

    ! Create matlab script
    open (13, file = trim(filename)//'_plotelements.m')
    write(13,*) '% Plotting Element Values'
    write(13,*) 'close all'
    write(13,*) 'clear all'
    write(13,*) trim(filename), '_plotelements_data;'
    write(13,*) '% Determine colorscale'
    write(13,*) 'cmap = colormap;'
    write(13,*) 'cinterp = linspace(min(plotval),max(plotval),size(cmap,1));'
    write(13,*) '% Make plot'
    write(13,*) 'title(',"'", trim(title) ,"'",')'
    write(13,*) 'hold on'
    ! Element dependent code
    ! NOTE: not possible to mix element types !!!
    if (element(1)%id .eq. 1) then
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '   [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
       write(13,*) '   xx = X(IX(e,1:2),1);'
       write(13,*) '   yy = X(IX(e,1:2),2);'
       write(13,*) '   plot(xx,yy,',"'",'Color',"'",',cmap(arr_pos,:),',"'",'Linewidth',"'",',1.5);'
       write(13,*) 'end'  
    elseif (element(1)%id == 2) then
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '   [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
       write(13,*) '   xx = X(IX(e,1:4),1);'
       write(13,*) '   yy = X(IX(e,1:4),2);'
       write(13,*) '   patch(xx,yy,cmap(arr_pos,:));'
       write(13,*) 'end'   
    else
       write (*,'("Unsupported element type in Matlab routine:")')
       write (*,'("plot -> plotmatlabelements")')  
    endif
    ! End file and close
    write(13,*) 'axis equal;'
    write(13,*) 'axis off;'
    write(13,*) 'caxis([min(plotval) max(plotval)]);'
    write(13,*) 'colorbar;'
    write(13,*) 'hold off'
    write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
    close(13) 

    ! Print to screen:
    print*,'Matlab plotting called: elements'
    print*,'Generated data file: ', trim(filename)//'_plotelements_data.m' 
    print*,'and script file: ', trim(filename)//'_plotelements.m'
    print*,' '

  end subroutine plotmatlabeval

  subroutine plotmatlabevec(title,ppp1,ppp2,pppang)
    !---------------------------------------------------------------------
    ! Vector field plot using Matlab

    real(8), dimension(:), intent(in) :: ppp1, ppp2, pppang
    character(len=*), intent(in) :: title   
    integer i, j, e

    ! write datafile     
    open (13, file = trim(filename)//'_plotvector_data.m')

    ! write nodal coordinates
    write(13,'("X = [")')
    do i = 1,size(x,1)
       write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
    enddo
    write(13,'("];")')
    write(13,'( )')

    ! write topology matrix
    write(13,'("IX = [")')
    do e = 1,ne
       write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
    enddo
    write(13,'("];")')

    ! write data1/element values
    write(13,'("vect = [")')
    do i = 1,ne
       write(13,*) ppp1(i), ppp2(i), pppang(i) 
    enddo
    write(13,'("];")')   
    close(13)

    ! Create matlab script
    open (13, file = trim(filename)//'_plotvector.m')
    write(13,*) '% Plotting Vector Field, i.e. principle stresses'
    write(13,*) 'close all'
    write(13,*) 'clear all'
    write(13,*) trim(filename), '_plotvector_data;'
    write(13,*) '% Make plot'
    ! Element dependent code
    ! NOTE: not possible to mix element types !!!
    if (element(1)%id .eq. 1) then
       print*,'LINK1 error: cannot do vector plot of truss structure'
       print*,'Files created are empty !!'
    elseif (element(1)%id == 2) then
       write(13,*) '% User scale parameter: See fedata - scale_vec'
       write(13,*) 'scale_vec = 1;'
       write(13,*) '% Define characteristic length'
       write(13,*) 'clx = 0;   cly = 0;'
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '    for i = 1:4'
       write(13,*) '        for j = i+1:4'
       write(13,*) '            clx = max(abs(  X(IX(e,i),1) - X(IX(e,j),1)   ),clx);'
       write(13,*) '            cly = max(abs(  X(IX(e,i),2) - X(IX(e,j),2)   ),cly);'
       write(13,*) '        end'
       write(13,*) '    end'
       write(13,*) 'end'
       write(13,*) 'clmax = max(clx, cly);'
       write(13,*) 'scal = max(max(abs(vect(:,1:2))))*sqrt(10)/clmax / scale_vec;'
       write(13,*) '% Make plot'
       write(13,*) 'figure'
       write(13,*) 'hold on'
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '    xx = X(IX(e,1:4),1);'
       write(13,*) '    yy = X(IX(e,1:4),2);'
       write(13,*) '    patch(xx,yy,[1 1 1]);'
       write(13,*) '    % Find approx. center for isoparametric element'
       write(13,*) '    xc = sum(X(IX(e,:),1))/size(IX,2);'
       write(13,*) '    yc = sum(X(IX(e,:),2))/size(IX,2);'
       write(13,*) '    % Directions for vect(:)'
       write(13,*) '    vec = [cos(-vect(e,3)) sin(-vect(e,3)) ...'
       write(13,*) '        cos(-vect(e,3)+pi/2) sin(-vect(e,3)+pi/2)];'
       write(13,*) '    % Plot magnitude and direction of vect_1'
       write(13,*) '    cc = ',"'", 'b',"'",';'
       write(13,*) '    if vect(e,1) < 0,    cc = ',"'",'r',"'",';     end'
       write(13,*) '    quiver(xc,yc,vec(1),vec(2),abs(vect(e,1))/scal,cc)'
       write(13,*) '    quiver(xc,yc,-vec(1),-vec(2),abs(vect(e,1))/scal,cc)'
       write(13,*) '    % Plot magnitude and direction of vect_2'
       write(13,*) '    cc = ',"'",'b',"'",';'
       write(13,*) '    if vect(e,2) < 0,    cc = ',"'",'r',"'",';     end'
       write(13,*) '    quiver(xc,yc,vec(3),vec(4),abs(vect(e,2))/scal,cc)'
       write(13,*) '    quiver(xc,yc,-vec(3),-vec(4),abs(vect(e,2))/scal,cc)'
       write(13,*) 'end   '
    else
       write (*,'("Unsupported element type in Matlab routine:")')
       write (*,'("plot -> plotmatlabevec")')  
    endif
    ! End file and close
    write(13,*) 'title( ',"'", trim(title),"'",')'
    write(13,*) 'axis equal;  axis off;  hold off'
    write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
    close(13)

    ! Print to screen:
    print*,'Matlab plotting called: Vector'
    print*,'Generated data file: ', trim(filename)//'_plotvector_data.m' 
    print*,'and script file: ', trim(filename)//'_plotvector.m'
    print*,' '

  end subroutine plotmatlabevec

  subroutine plotmatlabeig(title,lamb,evec,times)
    !---------------------------------------------------------------------
    ! Eigenmode plot using Matlab
    ! lamb  = eigenvalue
    ! evec  = eigenvector
    ! times = [totaltime, timeinterval]

    real(8), dimension(:), intent(in) :: evec,times
    real(8), intent(in) :: lamb
    character(len=*), intent(in) :: title   
    integer i, j, e

    ! write datafile        
    open (13, file = trim(filename)//'_ploteig_data.m')

    ! write nodal coordinates
    write(13,'("X = [")')
    do i = 1,size(x,1)
       write (13,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
    enddo
    write(13,'("];")')
    write(13,'( )')

    ! write topology matrix
    write(13,'("IX = [")')
    do e = 1,ne
       write (13,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
    enddo
    write(13,'("];")')

    ! write deformation vector
    write(13,'("D = [")')
    do i = 1,neqn
       write (13,'(f15.9,1x)') (evec(i) )
    enddo
    write(13,'("];")')   
    close(13)

    ! Plot script
    ! Create matlab script
    open (13, file = trim(filename)//'_ploteig.m')
    write(13,*) '% Plotting Eigenmodes'
    write(13,*) 'close all'
    write(13,*) 'clear all'
    write(13,*) trim(filename), '_ploteig_data;'
    write(13,*) 'lamb = ',lamb,';'
    write(13,*) 'timeint = ',times(2),';'
    write(13,*) 'timetot = ',times(1),';'
    write(13,*) '% Find max window size.'
    write(13,*) 'lxmin = min(X(:,1));        lxmax = max(X(:,1));'
    write(13,*) 'lymin = min(X(:,2));        lymax = max(X(:,2));'
    write(13,*) 'dxmin = min(D(1:2:end));    dxmax = max(D(1:2:end));'
    write(13,*) 'dymin = min(D(2:2:end));    dymax = max(D(2:2:end));'
    write(13,*) 'lxmin = lxmin - max(abs(dxmin),abs(dxmax))*1.05;'
    write(13,*) 'lxmax = lxmax + max(abs(dxmin),abs(dxmax))*1.05;'
    write(13,*) 'lymin = lymin - max(abs(dymin),abs(dymax))*1.05;'
    write(13,*) 'lymax = lymax + max(abs(dymin),abs(dymax))*1.05;'
    write(13,*) '% Make plot'
    write(13,*) 'figure'
    write(13,*) 'set(gcf,',"'",'color',"'",',[ 1  1 1]);'
    write(13,*) 'times = 0:timeint:timetot;'
    write(13,*) 'for i = 1:length(times)'
    write(13,*) '    tfact = sin(sqrt(lamb)*times(i));'
    write(13,*) '    clf;'
    write(13,*) '    hold on'
    write(13,*) '    for e = 1:size(IX,1)'
    ! Element dependent code
    ! NOTE: not possible to mix element types !!!
    if (element(1)%id .eq. 1) then
       write(13,*) '       edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];'
       write(13,*) '       xx = X(IX(e,1:2),1) + tfact*D(edof(1:2:4));'
       write(13,*) '       yy = X(IX(e,1:2),2) + tfact*D(edof(2:2:4));'
       write(13,*) '       plot(xx,yy,',"'",'b',"',","'",'LineWidth',"'",',1.5)'
    elseif (element(1)%id == 2) then
       write(13,*) '       edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)...'
       write(13,*) '          2*IX(e,3)-1 2*IX(e,3) 2*IX(e,4)-1 2*IX(e,4)];'
       write(13,*) '       xx = X(IX(e,1:4),1) + tfact*D(edof(1:2:8));'
       write(13,*) '       yy = X(IX(e,1:4),2) + tfact*D(edof(2:2:8));'
       write(13,*) '       patch(xx,yy,[1 1 0]);'
    else
       write (*,'("Unsupported element type in Matlab routine:")')
       write (*,'("plotmatlab -> plotmatlabdef")')  
    endif
    ! End file and close
    write(13,*) '    end'
    write(13,*) '    axis([lxmin lxmax lymin lymax])'
    write(13,*) '    axis off'
    write(13,*) '    title( ',"'", trim(title),"'",')'
    write(13,*) '    pause(0.01)'
    write(13,*) 'end'
    close(13)

    ! Print to screen:
    print*,'Matlab plotting called: eigenmode'
    print*,'Generated data file: ', trim(filename)//'_ploteig_data.m' 
    print*,'and script file: ', trim(filename)//'_ploteig.m'
    print*,' '

  end subroutine plotmatlabeig

end module processor
