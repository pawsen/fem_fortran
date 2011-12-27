MODULE plot_routiner


  use fedata
  use numeth
  IMPLICIT NONE

  PRIVATE

  !deformation, felt mv for struktur
  PUBLIC ::output_deformed,output_elements, output_term, plot_vector, write_structure
  PUBLIC :: WriteVTK
  PUBLIC ::output_forskydninger! skriver kun en vektor ud

  !matrix, vektor, mv
  PUBLIC ::output_vector, output_matrix,  output_fd_script

  ! animation mv
  PUBLIC :: output_anim, output_1DOF,output_data, output_center, plot_gauss

CONTAINS


  subroutine output_anim(Dcurr,n,plotval) !gemmer forskydninger og sp�ndinger for hvert kald. Til animation

    ! udskriver displ og stress/strain for hvert tidsskridt
	integer, INTENT(IN) :: n
    real(8),dimension(:), INTENT(IN) :: Dcurr
    real(8),optional, INTENT(IN) :: plotval(:)
    integer :: i,e,j
    Character(10000)  string
    ! m�ske character(len=*)


    if (n==1) then! write matlab datafile (topology) & plotting routine
       call make_dir

       call write_structure('_plotdeformed_topology.m.m')

       ! Create matlab script(plotting routine)
       open (13, file = trim(filename)//'dir/'//'plotelements.m')
       write(13,*) '% Plotting Element Values'
       write(13,*) 'close all'
       write(13,*) 'clear all'
       write(13,*) 'clc'
       write(13,*) '_plotdeformed_topology; % load topology'
       write(13,*) 'h = figure(1);'
       write(13,*) 'axis ([0 3.5 -.5 1.5]),axis equal, axis manual'
       write(13,*) 'hold on'
       write(13,*) '%plotval_upper = 1e-2;'
       write(13,*) '%plotval_lower = 0;'
       write(13,*) 'h = figure(1);'
       write(13,*) '%cmap = colormap;'
       write(13,*) '%cinterp = linspace(min(plotval_lower),max(plotval_upper),size(cmap,1));'
       write(13,*) 'for j = 1:1000' 
       write(13,*) '   clf'
       write(13,*) '   h = figure(1);'
       write(13,*) '   axis ([0 3.5 -.5 1.5]),axis equal, axis manual'
       write(13,*) '   importfile([',"'",trim(filename)//'def',"'",',num2str(j)])'
       write(13,*) '   D = data(:,2);'
       write(13,*) '   %DFFT(j)=data(end,2);'
       write(13,*) '   importfile([',"'",trim(filename)//'stress',"'",',num2str(j)])'
       write(13,*) '   plotval = data(:,2); %von Mises sp�nding'
       write(13,*) '   if j==1'
       write(13,*) '   	cmap = colormap;'
       write(13,*) '   	cinterp = linspace(min(plotval),max(plotval),size(cmap,1));'
       write(13,*) '   end'
       write(13,*) ''
       write(13,*) '   for e = 1:size(IX,1)'
       write(13,*) '      [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
       write(13,*) '      edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)...'
       write(13,*) '           2*IX(e,3)-1 2*IX(e,3) 2*IX(e,4)-1 2*IX(e,4)];'
       write(13,*) '      xx(edof(1:2:8)) = X(IX(e,1:4),1) + D(edof(1:2:8));'
       write(13,*) '      xx(edof(2:2:8)) = X(IX(e,1:4),2) + D(edof(2:2:8));'
       write(13,*) '      patch(xx(edof(1:2:8)),xx(edof(2:2:8)),cmap(arr_pos,:));%,[1 1 0]);'
       write(13,*) '   end'
       write(13,*) '   title(',"'",'Deformed',"'",')'
       write(13,*) '   M(j) = getframe(h);'
       write(13,*) 'end'
       write(13,*) ''
       write(13,*) '%movie(M,1,40) % play the movie'
       write(13,*) '%movie2avi(M,',"'",'vibration.avi',",'",'fps',",'",', 40',",'",'quality',"'",',10) %record the  movie'
       close(13)

       open (13, file = trim(filename)//'dir'//'importfile.m')
       write(13,*) 'function importfile(fileToRead1)'
       write(13,*) 'newData1 = importdata(fileToRead1);% Import the file'
       write(13,*) '% Create new variables in the base workspace from those fields.'
       write(13,*) 'vars = fieldnames(newData1);'
       write(13,*) 'for i = 1:length(vars)'
       write(13,*) '    assignin(',"'",'base',"'",',vars{i}, newData1.(vars{i}));'
       write(13,*) 'end'
       close(13)
    end if

    ! The following writes displacement & von mises stress to seperate files
    if (n<1e1) then ! giver filnavnet et nummer
       write(string,fmt=1) n
1      format (i1)
    elseif (n<1e2) then
       write(string,fmt=2) n
2      format (i2)
    elseif (n<1e3) then
       write(string,fmt=3) n
3      format (i3)
    elseif (n<1e4) then
       write(string,fmt=4) n
4      format (i4)
    else
       write(string,fmt=5) n
5      format (i5)
    endif

	! write displacement
    open (10, file = trim(filename)//'dir/'//trim(filename)//'def'//string)
    write (10, '("DOF       Displacement")')
    do i = 1, neqn
       write (10, '(i6, 1x, f32.15)') i, Dcurr(i)
    end do
    write(10,'("];")')
    close(10) 

    ! write stress
    if (present(plotval)) then
       open (10, file = trim(filename)//'dir/'//trim(filename)//'stress'//string)
       write(10,'("plotval = [")')
       do e = 1,ne
          write (10,'(i6,1x, f32.15)') e, plotval(e) 
       enddo
       write(10,'("];")')   
       close(10) 
    end if

  end subroutine output_anim

  subroutine output_vector(vec,title)

    ! Udskriver en vektor til matlab, fx compliance.

    real(8), dimension(:), INTENT(IN) :: vec
    character(len=*), intent(in) :: title
    integer :: i


    call make_dir

    open (14, file = trim(filename)//'dir/'//trim(title)//'.m')
    write(14,*) '% Konvergens af objektfunktionen'
    write(14,*) 'vec=['
	do i = 1,size(vec,1)
       write(14,*) vec(i)
 	end do
    write(14,*) '];'
    close(14)

  end subroutine output_vector

  subroutine output_forskydninger(matrix,title)
    ! Subroutine to plot the element values using Matlab

    real(8), dimension(:), INTENT(IN):: matrix
    character(len=*), intent(in) :: title
    integer i, j

    call make_dir

    ! write datafile     
    open (13, file = trim(filename)//'dir/'//trim(title)//'.m')

    ! write matrix
    write(13,'("mat = [")')
    do i=1,nn
       !print *,matrix(2*i-1),matrix(2*i)
       !      write (13,'(2(g32.25 2x))') matrix(2*i-1), matrix(2*i)
       write (13,*) matrix(2*i-1), matrix(2*i)

    end do
    write(13,'("];")')

    close(13)

    print*,'matrix printet'
  end subroutine output_forskydninger

  subroutine output_matrix(matrix,title)
    ! Subroutine to plot the element values using Matlab

    real(8), dimension(:,:), INTENT(IN):: matrix
    character(len=*), intent(in) :: title
    integer i, j

    call make_dir

    if (title == 'finite_diff_tjeck') then
       call output_fd_script
    end if

    ! write datafile     
    open (13, file = trim(filename)//'dir/'//trim(title)//'.m')

    ! write matrix
    write(13,'("mat = [")')
    do i=1,size(matrix,1)
       write (13,'(100(g32.15))') (matrix(i,j),j=1,size(matrix,2))
    end do
    write(13,'("];")')

    close(13)

    print*,'matrix printet'
  end subroutine output_matrix

  subroutine output_1DOF(D,n,nmax,deltaT,title)!

    ! Udskriver forskydningen for en enkelt frihedsgrad for hvert tidsskridt
	integer, INTENT(IN) :: n, nmax   
    real(8), INTENT(IN) :: D


    real(8), optional, INTENT(IN) :: deltaT
    character(len=*), intent(in) :: title

    if (n==1) then
       call make_dir


       open (14, file = trim(filename)//'dir/'//trim(title)//'.m')
       write(14,*) '% Plotting displacement for selected DOF'
       if (present(deltaT)) then; write(14,*) 'DeltaT=',deltaT,';' 
       endif
       write(14,*) 'D=['
    end if

    write(14,*) D

    if (n==nmax) then
       write(14,*) '];'
       close(14)
    end if

  end subroutine output_1DOF


  subroutine plot_vector(title,ppp1,ppp2,pppang)
    !---------------------------------------------------------------------
    ! Vector field plot using Matlab

    real(8), dimension(:), intent(in) :: ppp1, ppp2, pppang
    character(len=*), intent(in) :: title   
    integer i, j, e

    call make_dir
    ! write datafile   
    call write_structure('_plotvector_data.m',(/0d0/),(/DCMPLX(0d0,0d0)/),ppp1,ppp2,pppang)

    ! Create matlab script
    open (13, file = trim(filename)//'dir/'//trim(filename)//'_plotvector.m')

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
    elseif (element(1)%id == 2 .or. element(1)%id == 3) then
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

  end subroutine plot_vector

  SUBROUTINE output_fd_script

    call make_dir
    ! write datafile     
    open (13, file = trim(filename)//'dir/'//'fd_check_auto.m')

    write(13,*) 'clc ,clear all%, close all'
    write(13,*) 'finite_diff_tjeck;'
    write(13,*) 'fd = matrix;'
    write(13,*) ''
    write(13,*) '%fd(:,1) er fd_beregnede sensitiviter. fd(:,2) er analytiske sensitiviteter'
    write(13,*) ''
    write(13,*) 'figure(1),clf'
    write(13,*) 'plot(1:size(fd,1),fd(:,1),',"'",'xb',"'",')'
    write(13,*) 'hold on'
    write(13,*) 'plot(1:size(fd,1),fd(:,2),',"'",'sr',"'",')'
    write(13,*) 'hold off'
    write(13,*) 'legend(',"'",'fd\_beregnet',"'",",'",'eksakt',"'",",'",'location',"'",",'",'S',"'",')'
    write(13,*) ''
    write(13,*) 'figure(2),clf'
    write(13,*) 'fd_rel = abs(fd(:,2)-fd(:,1))./fd(:,2);'
    write(13,*) 'plot(1:size(fd,1),fd_rel,',"'",'xb',"'",')'
    write(13,*) 'axis([1 size(fd,1) -2 .01 ])%min(fd_rel) max(fd_rel)])'
    write(13,*) 'legend(',"'",'relativ forskel i sensitiviteter',"'",",'",'location',"'",",'",'NW',"'",')'
    write(13,*) ''
    write(13,*) 'min_rel_forskel = min(fd_rel)'
    write(13,*) 'max_rel_forskel = max(fd_rel)'
    write(13,*) ''
    write(13,*) '% find antallet af sensitiviter hvor den absolutte relative afvigelse er'
    write(13,*) '% over fx. 10%'
    write(13,*) 'find( abs(fd_rel) > 0.1)'
    write(13,*) ''
    write(13,*) '%absolut afvigelse'
    write(13,*) 'figure(3), clf'
    write(13,*) 'plot(1:size(fd,1),abs(fd(:,2)-fd(:,1)),',"'",'xb',"'",')'
    write(13,*) 'legend(',"'",'absolut forskel i sensitiviteter',"'",",'",'location',"'",",'",'NW',"'",')'
    close(13)
  END SUBROUTINE output_fd_script

 
  Subroutine write_structure(title,plotval,comp_plotval,ppp1,ppp2,pppang)
    use fedata

    character(len=*), intent(in) :: title
    real(8), optional,dimension(:), intent(in)::plotval, ppp1,ppp2,pppang
    COMPLEX(8), optional, INTENT(IN) :: comp_plotval(:)
    integer, parameter :: fid =13
    integer :: e, i, j

    open (fid, file = trim(filename)//'dir/'//trim(filename)//title)
    ! write nodal coordinates
    write(fid,'("X = [")')
    do i = 1,size(x,1)
       write (fid,'(3(f15.9,1x))') (x(i,j), j=1,size(x,2))
    enddo
    write(fid,'("];")')
    write(fid,'( )')

    ! write topology matrix
    write(fid,'("IX = [")')
    do e = 1,ne
       write (fid,'(20(i6,1x))') (element(e)%ix(i) ,i=1,size(element(e)%ix,1))
    enddo
    write(fid,'("];")')

    ! hvis title er 'plotdeformed_topology.m', skal KUN topologien udskrives
    if (trim(title) /= 'plotdeformed_topology.m') then
       ! write deformation vector
       !write(fid,'("i = sqrt(-1);")')
       if (harmonic) then
          write(fid,'("D_real = [")')
          do i = 1,neqn
             write (fid,10) REAL(dZ(i))            
          enddo
          write(fid,'("];")')

          write(fid,'("")')
          write(fid,'("D_imag = [")')
          do i = 1,neqn
             write (fid,10) DIMAG(dZ(i))
          end do
          write(fid,'("];")')
       else
          write(fid,'("D = [")')
          do i = 1,neqn
             write (fid,'(g32.15,1x)') (d(i) )
          enddo
          write(fid,'("];")')
       end if
    end if

!!$             if (DIMAG(dZ(i)) < 0) then
!!$                write (fid,10) REAL(dZ(i)),'-i*',-DIMAG(dZ(i))
!!$            else
!!$               write (fid,10)  REAL(dZ(i)),'+i*',DIMAG(dZ(i))
!!$             end if 

    if (present(plotval) .and. (trim(title) /= '_plotvector_data.m')) then
       ! write data1/element values

       if (harmonic) then
          write(fid,'("plotval_real = [")')
          do i = 1,size(comp_plotval,1)
             write (fid,10) REAL(comp_plotval(i))
          end do
          write(fid,'("];")')
          write(fid,'("")')

          write(fid,'("plotval_imag = [")')
          do i = 1,size(comp_plotval,1)
             write (fid,10) DIMAG(comp_plotval(i))
          end do
          write(fid,'("];")')  

       else
          write(fid,'("plotval = [")')
          do i = 1,size(plotval,1)
             write (fid,'(g32.15)') plotval(i)
          enddo
          write(fid,'("];")')  
       end if
    end if
10  format(g32.15)

    if (present(ppp1)) then
       !if (MAXVAL(DABS(ppp1))>1E-10) then
          ! write data1/element values
          write(fid,'("vect = [")')
          do i = 1,ne
             write(fid,*) ppp1(i), ppp2(i), pppang(i) 
          enddo
          write(fid,'("];")') 
       !end if
    end if

    close(fid) 

  end Subroutine write_structure


  SUBROUTINE output_deformed(title,felt, plotval,cmplx_plotval)

    ! udskriver den deformerede struktur til matlab-fil

    character(len=*), intent(in) :: title, felt
    real(8), dimension(:), intent(in) :: plotval
    complex(8),optional, dimension(:), intent(in) :: cmplx_plotval
    integer i, j, e

    call make_dir
    ! write datafile
    if (harmonic) then
       call write_structure('_plotdeformed_data.m',(/0d0/),cmplx_plotval)
    else
       call write_structure('_plotdeformed_data.m',plotval)
    end if

    ! Create matlab script
    open (13, file = trim(filename)//'dir/'//trim(filename)//'_plotdeformed.m')
    write(13,*) '% Plotting Un-Deformed and Deformed Structure'
    write(13,*) '%close all'
    write(13,*) 'clear all; clc'
    write(13,*) trim(filename), '_plotdeformed_data;'

    if (harmonic) then
       write(13,*) ''
       write(13,*) 'D = complex(D_real,D_imag);'
       write(13,*) 'plotval = complex(plotval_real,plotval_imag);'
       write(13,*) ''
       write(13,*) 'frekvens = 5e4;'
       write(13,*) 'T = 1/frekvens;'
       write(13,*) 't = T*0.5;'
       write(13,*) ''
       write(13,*) 'omega = frekvens*2*pi;'
       write(13,*) 'plotval= real(plotval*exp(i*omega*t));'
       write(13,*) 'D = real(D*exp(i*omega*t));'
       write(13,*) ''
    end if

    write(13,*) '% Make plot'
    write(13,*) 'figure(1), clf'
    write(13,*) 'set(gcf,',"'",'Name',"','", trim(title) ,"'",')'
    write(13,*) 'cmap = colormap;'
    write(13,*) 'cinterp = linspace(min(plotval),max(plotval),size(cmap,1));'

    ! Element dependent code
    write(13,*) 'subplot(2,1,1)'
    write(13,*) 'hold on'
    if (size(plotval,1)== ne) then
       if (element(1)%id == 2) then
          write(13,*) 'for e = 1:size(IX,1)'
          write(13,*) '   [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
          write(13,*) '    xx = X(IX(e,1:4),1);'
          write(13,*) '    yy = X(IX(e,1:4),2);'
          write(13,*) '   patch(xx,yy,cmap(arr_pos,:), ',"'",'EdgeColor',"'",',',"'",'none',"'",');'
          write(13,*) 'end'
       elseif (element(1)%id == 3 .or. elem_type == 'PLATE_GMSH') then
          write(13,*) 'for e = 1:size(IX,1)'
          write(13,*) '    [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
          write(13,*) '    edof = [3*IX(e,1)-2 3*IX(e,1)-1 3*IX(e,1)...'
          write(13,*) '         3*IX(e,2)-2 3*IX(e,2)-1 3*IX(e,2)...'
          write(13,*) '         3*IX(e,3)-2 3*IX(e,3)-1 3*IX(e,3)...'
          write(13,*) '         3*IX(e,4)-2 3*IX(e,4)-1 3*IX(e,4)]; '
          write(13,*) '    xx = X(IX(e,1:4),1);'
          write(13,*) '    yy = X(IX(e,1:4),2);'
          write(13,*) '    test = D(edof(1:3:12))*5e2;'
          write(13,*) '    zz = [0; 0; 0; 0]+test;'
          write(13,*) '    patch(xx,yy,zz,cmap(arr_pos,:), ',"'",'EdgeColor',"'",',',"'",'none',"'",');'
          write(13,*) 'end'
       END if
    else
       if (element(1)%id == 2) then
          write(13,*) 'for e = 1:size(IX,1)'
          write(13,*) '    xx = X(IX(e,1:4),1);'
          write(13,*) '    yy = X(IX(e,1:4),2);'
          write(13,*) '   patch(xx,yy,plotval(IX(e,1:4)), ',"'",'EdgeColor',"'",',',"'",'none',"'",');'
          write(13,*) 'end'
       elseif (element(1)%id == 3 .or. elem_type == 'PLATE_GMSH') then
          write(13,*) 'for e = 1:size(IX,1)'
          write(13,*) '    edof = [3*IX(e,1)-2 3*IX(e,1)-1 3*IX(e,1)...'
          write(13,*) '         3*IX(e,2)-2 3*IX(e,2)-1 3*IX(e,2)...'
          write(13,*) '         3*IX(e,3)-2 3*IX(e,3)-1 3*IX(e,3)...'
          write(13,*) '         3*IX(e,4)-2 3*IX(e,4)-1 3*IX(e,4)]; '
          write(13,*) '    xx = X(IX(e,1:4),1);'
          write(13,*) '    yy = X(IX(e,1:4),2);'
          write(13,*) '    test = D(edof(1:3:12))*5e2;'
          write(13,*) '    zz = [0; 0; 0; 0]+test;'
          write(13,*) '    patch(xx,yy,zz,plotval(IX(e,1:4)), ',"'",'EdgeColor',"'",',',"'",'none',"'",');'
          write(13,*) 'end'
       END if
    end if
    write(13,*) 'title(',"'",'Deformed + felt( ',trim(felt) ,')',"'",')' 
    write(13,*) 'axis equal;'
    !write(13,*) 'xaxes = get(gca,',"'",'xlim',"'",');'
    !write(13,*) 'yaxes = get(gca,',"'",'ylim',"'",');'
    write(13,*) '%axis off;'
    write(13,*) 'caxis([min(plotval) max(plotval)]);'
    write(13,*) 'colorbar;'
    write(13,*) ''

    ! PLOT2
    write(13,*) 'subplot(2,1,2)'
    write(13,*) 'hold on'
    if (element(1)%id == 2) then
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)...'
       write(13,*) '        2*IX(e,3)-1 2*IX(e,3) 2*IX(e,4)-1 2*IX(e,4)];'
       write(13,*) '    xx = X(IX(e,1:4),1) + D(edof(1:2:8));'
       write(13,*) '    yy = X(IX(e,1:4),2) + D(edof(2:2:8));'
       write(13,*) '    color = [norm(D(edof(1:2))); norm(D(edof(3:4))); norm(D(edof(5:6))); norm(D(edof(7:8)))];'
       write(13,*) '   patch(xx,yy,color, ',"'",'EdgeColor',"'",',',"'",'none',"'",');'
       write(13,*) 'end'
    elseif (element(1)%id == 3 .or. elem_type == 'PLATE_GMSH') then
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '    edof = [3*IX(e,1)-2 3*IX(e,1)-1 3*IX(e,1)...'
       write(13,*) '         3*IX(e,2)-2 3*IX(e,2)-1 3*IX(e,2)...'
       write(13,*) '         3*IX(e,3)-2 3*IX(e,3)-1 3*IX(e,3)...'
       write(13,*) '         3*IX(e,4)-2 3*IX(e,4)-1 3*IX(e,4)]; '
       write(13,*) '    xx = X(IX(e,1:4),1);'
       write(13,*) '    yy = X(IX(e,1:4),2);'
       write(13,*) '    test = D(edof(1:3:12))*5e2;'
       write(13,*) '    zz = [0; 0; 0; 0]+test;'
       write(13,*) '    patch(xx,yy,zz,D(edof(1:3:12)), ',"'",'EdgeColor',"'",',',"'",'none',"'",');'
       write(13,*) 'end'
    end if
    write(13,*) 'title(',"'",'Deformed + forskydninsfelt',"'",')'      
    !write(13,*) 'axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...'
    !write(13,*) ' min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);'
    write(13,*) 'axis equal;'
    write(13,*) '%axis off;'

    ! End file and close
    !write(13,*) '%caxis([min(plotval) max(plotval)]);'
    write(13,*) 'colorbar;'
    write(13,*) 'hold off'
    close(13)

    ! Print to screen:
    print*,'Matlab plotting called: un/deformed form plot_routiner'
    print*,' '

  end subroutine output_deformed

  SUBROUTINE output_elements(title, data1)

    ! Subroutine to plot the element values using Matlab

    character(len=*), intent(in) :: title
    real(8), dimension(:), intent(in) :: data1
    integer i, j, e

    call make_dir
    ! write datafile     
    call write_structure('_plotelements_data.m',data1)

    ! Create matlab script
    open (13, file = trim(filename)//'dir/'//trim(filename)//'_plotelements.m')
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


    if (size(data1,1)==ne) then
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '   [dummy,arr_pos] = min(abs(cinterp-plotval(e)));'
       write(13,*) '   xx = X(IX(e,1:4),1);'
       write(13,*) '   yy = X(IX(e,1:4),2);'
       write(13,*) '   patch(xx,yy,cmap(arr_pos,:));'
       write(13,*) 'end'  
    else ! there is a value for each node
       write(13,*) 'for e = 1:size(IX,1)'
       write(13,*) '   xx = X(IX(e,1:4),1);'
       write(13,*) '   yy = X(IX(e,1:4),2);'
       write(13,*) '   patch(xx,yy,plotval(IX(e,1:4)));'
       write(13,*) 'end'  
    end if

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
    print*,' '

  end subroutine output_elements



  ! transiente output-rutiner

  subroutine plot_gauss(nmax,deltaT, fc, T_0, delta)

    ! Subroutine plotter den benyttede gauss-b�lgepakke

    INTEGER, INTENT(IN):: nmax
    REAL(8), INTENT(IN):: deltaT, fc, T_0, delta

    ! Create matlab script(plotting routine)
    call make_dir
    open (13, file = trim(filename)//'dir/'//'gauss.m')

    write(13,*) '% GAUSS-PULS b�lgepakke'
    write(13,*) 'clear all'
    write(13,*) 'clc'
    write(13,*) 'fc = ',fc,';% center af b�lgepakke i frekvensdom�ne'
    write(13,*) 'deltaT = ',deltaT,';% tidsskridt'
    write(13,*) 'nmax = ',nmax,';'
    write(13,*) 'delta = ',delta,';% b�ndbredde'
    write(13,*) 't_0 = ',T_0,';%center af b�lgepakke i tidsdom�ne'
    write(13,*) 'n = linspace(1,nmax,nmax);'
    write(13,*) 't = deltaT*n;'
    write(13,*) ''
    write(13,*) 'f = cos(2*pi*fc.*(t-t_0)).*exp(-delta.*(t/t_0-1).^2); % excitationskraft;'
    write(13,*) ''
    write(13,*) ''
    write(13,*) 'fft_f   = fft(f) / length(f);'
    write(13,*) 'omega   = 1/(2*max(t)/length(t))* linspace(0, 1, round(length(t) / 2));'
    write(13,*) ''
    write(13,*) 'figure(2)'
    write(13,*) 'clf'
    write(13,*) 'subplot(2,1,1); plot(t,f);'
    write(13,*) 'subplot(2,1,2); plot(omega, abs( fft_f(1:round(length(f)/2)) ) );'

    close(13) 

  end subroutine plot_gauss


  subroutine output_data(title,rho,deltaT,dev)

    ! til b�lgehastighed

    character(len=*), intent(in) :: title
    real(8), INTENT(IN) :: deltaT, rho(:)
    integer,optional, INTENT(IN) :: dev
    real(8) :: dens, young, nu, lambda, mu, cL, cT

    ! linux_fejl
    !   dens = rho(element_end)*dens1+(1.0-rho(element_end))*dens2
    !   young = rho(element_end)*young1+(1.0-rho(element_end))*young2
    !   nu = rho(element_end)*nu1+(1.0-rho(element_end))*nu2

    ! Lame parameters  
    lambda = young*nu/ ((1.0+nu)*(1.0-2.0*nu)) 
    mu = young/(2.0*(1.0+nu))! shear modulus

    cL = dsqrt( (lambda+2.0*mu) / dens)
    cT = dsqrt( mu / dens)  

    call make_dir
    ! write datafile     
    open (14, file = trim(filename)//'dir/'//trim(title)//'.m')
    write(14,*) '% '
    write(14,*) 'deltaT= ',deltaT,';'
    write(14,*) 'cL= ',cL,';'
    write(14,*) 'cT= ',cT,';'
    write(14,*) 'dens= ',dens,';'
    write(14,*) 'young= ',young,';'
    write(14,*) 'nu= ',nu,';'
    write(14,*) 'lambda= ',lambda,';'
    write(14,*) 'mu= ',mu,';'
    if (present(dev)) then; write(14,*) 'dev= ',dev,';'; 
    end if
    close(14)


  end subroutine output_data

  subroutine output_center(D,n,nmax,deltaT,rho)

    ! udskriver D for de valgte frihedsgrader
    ! fx. forskydningen i x-retningen for alle centerknuder.
    ! bruges fx. til at animere b�lgens "vej" gennem bj�lken

    ! skriver alle forskydninger for et enkelt tidsskridt til samme fil
    ! dvs antallet af filer er lig antallet af tidsskridt

    integer, INTENT(IN) :: n, nmax
    real(8), DIMENSION(:), INTENT(IN) :: D, rho
    real(8), INTENT(IN) :: deltaT
    integer :: i,e,j
    Character(10000)  string

    if (n==1) then
       call output_data('data',rho,deltaT)
    end if

    ! The following writes displacement & von mises stress to seperate files

    if (n<1e1) then ! giver filnavnet et nummer
       write(string,fmt=1) n
1      format (i1)
    elseif (n<1e2) then
       write(string,fmt=2) n
2      format (i2)
    elseif (n<1e3) then
       write(string,fmt=3) n
3      format (i3)
    elseif (n<1e4) then
       write(string,fmt=4) n
4      format (i4)
    else
       write(string,fmt=5) n
5      format (i5)
    endif

    ! write displacement
    open (10, file = trim(filename)//'dir/'//trim(filename)//'wave'//string) ! skriv fil p�ny. Dvs slet indhold
    do i=1,SIZE(D,1)
       write(10,'(f32.15)') d(i)
    end do
    close(10)


  end subroutine output_center

  subroutine output_term(t_elem,tn)

    !---------------------------------------------------------------------
    !  This subroutine outputs the thermal results.

    integer :: i, e
    real(8), dimension(:), intent(IN) :: t_elem, tn

    write (*, *)
    write (*, '("Writing fem output file to ", a)') trim(filename)//'_term.out'
    open (10, file = trim(filename)//'_term.out')

    write (10, '(" Alpha     Kcond ")') 
    write (10, '(f8.3, 1x, f8.3)') mprop(element(1)%mat)%alpha, mprop(element(1)%mat)%kcond

    write (10, '(" Tstart      Thk ")')
    write (10, '(f8.3, 1x, f8.3)') mprop(element(1)%mat)%tstart, mprop(element(1)%mat)%thk

    write (10, *)
    write (10, '(" Node            Temperature                Element                      Temperature")')


    write (10, *)
    write (10, '(" Node      Temperature ")')
    do i = 1, nn
       write (10, '(i5, 6x, g25.15)') i, tn(i)
    end do

    write (10, *)
    write (10, '("Element    Temperature ")')
    do e = 1, ne
       write (10, '(i5, 6x, g25.15)') e, t_elem(e)
    end do

    close (10)

  end subroutine output_term

  !***********************************************************************
  subroutine WriteVTK(plotval,outfile,struct)
    !***********************************************************************
    ! outfile = name of file (.vtu)
    ! struct = 'deformed' or 'undeformed'
    character(len=*), intent(in) :: outfile, struct
    real(8), dimension(:), intent(in) :: plotval
    integer :: i

    open (unit=1,file=trim(outfile)//'.vtu')
    write(1,'(a70)')'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
    write(1,'(a20)')'  <UnstructuredGrid>'
    ! write(1,'(a27i5a17i5a2)')'    <Piece NumberOfPoints="',nn,'" NumberOfCells="',ne,'">' !linux_fejl
    write(1,'(a34)')'      <CellData Scalars="scalars">'
    write(1,'(a65)')'        <DataArray type="Float32" Name="VonMises" Format="ascii">'
    do i = 1,ne,1
       write(1,'(a10 e12.5)')'          ',plotval(i)
    end do
    write(1,'(a20)')'        </DataArray>'
    !------------
    !$$$$$$  write(1,'(a95)')'        <DataArray type="Float32" Name="PricipalStress1" NumberOfComponents="3" Format="ascii">'
    !$$$$$$        do i = 1,ne,1
    !$$$$$$         write(1,'(a10f8.3a1f8.3a1f8.3)')'          ',abs(Pstress1(i))*cos(-Psi(i)),' ',abs(Pstress1(i))*sin(-Psi(i)),' ',0.0
    !$$$$$$        end do
    !$$$$$$  write(1,'(a20)')'        </DataArray>'
    !$$$$$$ !------------
    !$$$$$$  write(1,'(a95)')'        <DataArray type="Float32" Name="PricipalStress2" NumberOfComponents="3" Format="ascii">'
    !$$$$$$        do i = 1,ne,1
    !$$$$$$         write(1,'(a10f8.3a1f8.3a1f8.3)')'          ', &
    !$$$$$$         abs(Pstress2(i))*cos(-Psi(i)+pi/2.),' ',abs(Pstress2(i))*sin(-Psi(i)+pi/2.),' ',0.0
    !$$$$$$        end do
    !$$$$$$  write(1,'(a20)')'        </DataArray>'
    !$$$$$$ !------------ 
    !$$$$$$  write(1,'(a86)')'        <DataArray type="Float32" Name="Stress" NumberOfComponents="3" Format="ascii">'
    !$$$$$$        do i = 1,ne,1
    !$$$$$$         write(1,'(a10f8.3a1f8.3a1f8.3)')'          ',stress(i,1),' ',stress(i,2),' ',0.0
    !$$$$$$        end do
    !$$$$$$  write(1,'(a20)')'        </DataArray>'
    !$$$$$$ !------------ 
    !$$$$$$  write(1,'(a86)')'        <DataArray type="Float32" Name="Strain" NumberOfComponents="3" Format="ascii">'
    !$$$$$$        do i = 1,ne,1
    !$$$$$$         write(1,'(a10f8.3a1f8.3a1f8.3)')'          ',strain(i,1),' ',strain(i,2),' ',0.0
    !$$$$$$        end do
    !$$$$$$  write(1,'(a20)')'        </DataArray>'
    !------------ 
    write(1,'(a17)')'      </CellData>'
    write(1,'(a14)')'      <Points>'
    write(1,'(a72)')'        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
    do i = 1,nn,1 
       select case( struct )
       case( 'undeformed' )
          write(1,'(a10f8.3a1f8.3a1f8.3)')'          ',x(i,1),' ',x(i,2),' ',0.0
       case( 'deformed' ) 
          write(1,'(a10f8.3a1f8.3a1f8.3)')'          ',x(i,1)+d(2*i-1),' ' &
               ,x(i,2)+d(2*i),' ',0.0
       case default
          write(*, *) 'Error: Illegal option for structure plot'
          stop
       end select
    end do
    write(1,'(a20)')'        </DataArray>'
    write(1,'(a15)')'      </Points>'
    write(1,'(a13)')'      <Cells>'
    write(1,'(a67)')'        <DataArray type="Int32" Name="connectivity" Format="ascii">'
    do i = 1,ne,1
       write(1,'(a10i5a1i5a1i5a1i5)')'          ',element(i)%ix(1)-1,' ',&
            element(i)%ix(2)-1,' ',element(i)%ix(3)-1,' ',element(i)%ix(4)-1
    end do
    write(1,'(a20)')'        </DataArray>'
    write(1,'(a62)')'        <DataArray type="Int32" Name="offsets" Format="ascii">'
    do i = 1,ne,1
       write(1,'(a10i5)')'          ',4*i
    end do
    write(1,'(a20)')'        </DataArray>'
    write(1,'(a60)')'        <DataArray type="Int32" Name="types" Format="ascii">'
    do i = 1,ne,1
       write(1,'(a10i1)')'          ',7
    end do
    write(1,'(a20)')'        </DataArray>'
    write(1,'(a14)')'      </Cells>'
    write(1,'(a12)')'    </Piece>'
    write(1,'(a21)')'  </UnstructuredGrid>'
    write(1,'(a10)')'</VTKFile>'
    close(1)
    return

  end subroutine WriteVTK


END MODULE plot_routiner
