module WriteEnsight

 ! This module contains subroutines specific to the PLANE42 element.

 IMPLICIT NONE
 
 PRIVATE
 PUBLIC :: WriteEnsightGeo,WriteEnsightVar,EnsightCase

CONTAINS

subroutine WriteEnsightGeo(outfile)!,WriteBinary)
!
!    ******************************************************************************
!                      subrout  WriteRectEnsightGeo
!    ******************************************************************************
!
!    writes mesh data in Ensight's ascii format for rectilinear geometry
      use fedata

	  character(len=*)		,INTENT(IN)::outfile
!      LOGICAL               ,INTENT(IN)::WriteBinary

      LOGICAL ::WriteBinary = .false.

!	character(LEN=80)::buffer
	character(LEN=80)::binary_form
	character(LEN=80)::file_description1,file_description2
	character(LEN=80)::node_id,element_id,element_type
	character(LEN=80)::part,description_part

      integer::FileUnit,i,j,npart
      integer::reclength

      FileUnit = 40

      binary_form      ='C Binary'
      file_description1='Ensight Model Geometry File Created by '
      file_description2='WriteEnsightGeo Routine'
      node_id          ='node id off'
      element_id       ='element id off'
      part             ='part'
      npart            =1
      element_type 	   ='quad4'
      reclength=80*8+4*(4*ne)+8*(3*ne)

      if (WriteBinary) then
        open (unit=FileUnit,file=trim(outfile)//'.geo', &
             form='UNFORMATTED')!,status='unknown')!access="direct",recl=reclength)
        write(unit=FileUnit) 	  binary_form                   &
                                 ,file_description1             &
                                 ,file_description2             &
                                 ,node_id                       &
                                 ,element_id                    &
                                 ,part,npart                    &
                                 ,description_part              &
                                 ,nn                            &
                                 ,(x(i,1),i=1,nn)       &
                                 ,(x(i,2),i=1,nn)       &
                                 ,(0d0,i=1,nn)                  &
                                 ,element_type                  &
                                 ,ne                            &
        ,(element(i)%ix(1),element(i)%ix(2),element(i)%ix(3),element(i)%ix(4),i=1,ne)

!$$$$$$         write(FileUnit) binary_form 
!$$$$$$         write(FileUnit)                         file_description1
!$$$$$$         write(FileUnit)                         file_description2
!$$$$$$         write(FileUnit)                         node_id
!$$$$$$        write(FileUnit)                          element_id
!$$$$$$         write(FileUnit)                         part
!$$$$$$         write(FileUnit)  npart
!$$$$$$         write(FileUnit)                         description_part
!$$$$$$         write(FileUnit)                         'coordinates'
!$$$$$$         write(FileUnit)                         nn
!$$$$$$         write(FileUnit)                         (x(i,1),i=1,nn)
!$$$$$$         write(FileUnit)                         (x(i,2),i=1,nn)
!$$$$$$         write(FileUnit)                         (0d0,i=1,nn)
!$$$$$$         write(FileUnit)                         element_type
!$$$$$$         write(FileUnit)                    ne
!$$$$$$         write(FileUnit) (element(i)%ix(1),element(i)%ix(2),element(i)%ix(3),element(i)%ix(4),i=1,ne)

      else
        open (unit=FileUnit,file=trim(outfile)//'.geo')
        write(FileUnit,'(A)') 'Ensight Model Geometry File Created by '
        write(FileUnit,'(A)') 'WriteRectEnsightGeo Routine'
        write(FileUnit,'(A)') 'node id off'!'node id given'
        write(FileUnit,'(A)') 'element id off'!'element id given'
		! undeformed size
        write(FileUnit,'(A)') 'part'
        write(FileUnit,'(i10)')npart
        write(FileUnit,'(A)') 'quad4' 					! beskrivelse
        write(FileUnit,'(A)') 'coordinates'
        write(FileUnit,'(i10)')nn						! number nodes
!        write(FileUnit,'(i10)') (i,i=1,nn)				! numbering of nodes
        write(FileUnit,'(E12.5)') (x(i,1),i=1,nn)		! x-coordinate
        write(FileUnit,'(E12.5)') (x(i,2),i=1,nn)		!
        write(FileUnit,'(E12.5)') (0d0,i=1,nn)			! z-coordinate = 0
        write(FileUnit,'(A)') 'quad4'					! element type
        write(FileUnit,'(i10)')ne						! number elements
!        write(FileUnit,'(i10)') (i,i=1,ne)				! numbering of elements

! fejl_linux
!		write(FileUnit,'(i10 i10 i10 i10)') &			! connectivity
!			(element(i)%ix(1),element(i)%ix(2),element(i)%ix(3),element(i)%ix(4),i=1,ne)

      end if
      close(FileUnit)
end subroutine WriteEnsightGeo


!    ******************************************************************************
!
!    WriteEnsightSca writes result data in Ensight's format
!
!    m1,m2,m3....: size of the variable in the x1,x2,x3 direction
!    ndv.........: number of dimension of the variable (1=>scalar   3=>vector) 
!    var.........: data to be written
!    Varname.....: word used to build filenames
!    imin,imax...: range of writting data in the x1 direction
!    jmin,jmax...: range of writting data in the x2 direction
!    kmin,kmax...: range of writting data in the x3 direction
!
!    ******************************************************************************

subroutine WriteEnsightVar(plotval,VarName)

	use fedata

	  character(len=*)		,INTENT(IN)::VarName
      REAL(8), dimension(:), intent(in) :: plotval
!      LOGICAL     ,INTENT(IN)::WriteBinary

      INTEGER::ndv = 1
	  LOGICAL:: WriteBinary = .false.
      character(len=80):: VarFileName,buffer
      character(len=80):: part,element_type
      integer::FileUnit,i,j,npart,m,reclength

      FileUnit = 40
      part ='part'
      npart=1
      element_type='quad4'
      reclength=80*3+4*(1)+8*((ne)*ndv)

      if (ndv.eq.1)VarFileName = trim(Varname)//'.scl'
      if (ndv.eq.3)VarFileName = trim(Varname)//'.vec'
      
!      write(*,'(5x,A)') VarFileName

	if(WriteBinary) then
        open (unit=FileUnit,file=VarFileName,&
        FORM     = 'UNFORMATTED')!, ACCESS='SEQUENTIAL')
       !      form='UNFORMATTED')!,access="direct",recl=reclength)
        write(FileUnit) VarFileName	
        write(FileUnit)	part
        write(FileUnit)	npart
        write(FileUnit) element_type
        write(FileUnit)	(plotval(i),i=1,ne)

      else
        open (unit=FileUnit,file=VarFileName)
        write(buffer,'(a,a)') Varname  
        write(FileUnit,'(A)')  buffer   
        write(FileUnit,'(A)') 'part'
        write(FileUnit,'(I10)')npart
        write(FileUnit,'(A)') 'quad4'					! element type
        write(FileUnit,'(e12.5)')(plotval(i),i=1,ne)
      endif
      close(FileUnit)

end  subroutine WriteEnsightVar



!
!    ******************************************************************************
!   EnsightCase helps to write a Ensight's case file
! 
!    VarName.....: Name of the variable
!    GeoName.....: Name of the geometrie
!    VarType.....: 1 => Scalar       3 => Vector
!    ntini.......: filename start number
!    nstop.......: filename end number
!    nprint......: filename increment
!
!    nfile.......: number of result files (time steps)
! 
subroutine EnsightCase(VarName,ntini,nstop,nprint)

     use fedata

      INTEGER,INTENT(IN)::nstop,ntini,nprint
      character(len=*),INTENT(IN)::Varname
      character(len=40):: plot_type 
      integer::FileUnit,i,nfile,VarType,element_val

!!$      VarType = 1
!!$      element_val = 1
!!$      
!!$
!!$      write(*,'(/2A)') ' Creating case file for Ensight and Paraview: ' &
!!$                       ,Varname
!!$
!!$      nfile=(nstop-ntini+1)/nprint		! number of files(=1 for non-transient)
!!$
!!$      FileUnit = 40
!!$      open(FileUnit,file=trim(filename)//'.case')
!!$
!!$      write(FileUnit,10) trim(filename)//'.geo'
!!$   10 format( &
!!$     'FORMAT'            ,/ ,		&
!!$     'type: ensight gold',//,		&
!!$     'GEOMETRY'          ,/ ,		&
!!$     'model:	',A         ,//,	&
!!$     'VARIABLE')
!!$
!!$ select case(element_val)
!!$   	case(0) ! knudev�rdi
!!$		plot_type = 'node'
!!$   	case(1)	! elementv�rdi
!!$		plot_type = 'element'
!!$ end select
!!$
!!$      if (nfile.eq.1) then
!!$        if(VarType.eq.1) &
!!$         write(FileUnit,15)trim(plot_type),trim(Varname),trim(filename)//'.scl'
!!$        if(VarType.eq.3) &
!!$         write(FileUnit,25)trim(plot_type),trim(Varname),trim(filename)//'.vec'
!!$      else
!!$        if(VarType.eq.1) &
!!$        	write(FileUnit,15)trim(plot_type),trim(Varname),trim(filename)//'******.scl'
!!$        if(VarType.eq.3) &
!!$        	write(FileUnit,25)trim(plot_type),trim(Varname),trim(filename)//'******.vec'
!!$        write(FileUnit,45) nfile,ntini,nprint
!!$        write(FileUnit,'(f15.3)') (ntini+float(i)*nprint,i=1,nfile)
!!$      endif
!!$
!!$      close(FileUnit)
!!$
!!$
!!$   15 format('scalar per ',A,': ',A,'   ', A)
!!$   25 format('vector per ',A,': ',A,'   ', A)
!!$
!!$   45 format(	&
!!$     /,'TIME            '      ,		&
!!$     /,'time set: 1     '      ,		&
!!$     /,'number of steps:'      ,i4 ,	&
!!$     /,'filename start number:',i10	&
!!$     /,'filename increment:'   ,i4		&
!!$     /,'time values: '					&
!!$     )

end subroutine EnsightCase

end module WriteEnsight
