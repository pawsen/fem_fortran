      SUBROUTINE WRITE_TO_ESG_RESULTS_FILE
!!
!! Copyright (c) by FMA Development, LLC, 28-SEP-2005
!!
!! Purpose: Write results to EnSight Gold C-binary files suitable for
!! EnSight or ParaView postprocessing. Note: EnSight's "UNFORMATTED"
!! Fortran file option (Fortran Binary) could not be used because 
!! ParaView does not support it. See the "FORM = 'BINARY'" directive
!! in the OPEN statements contained in this file.
!!
      USE shared_common_data
!!
!! The complete simulation data set.
!!
      USE indx_;           USE node_;           USE input_function_;
      USE beam_;           USE coord_;          USE sliding_interface_;
      USE value_;          USE force_;          USE nodal_constraints_;
      USE hexah_;          USE penta_;          USE nonreflecting_bc_;
      USE tetra_;          USE lsold_;          USE nodal_point_mass_;
      USE membq_;          USE membt_;          USE rigid_body_mass_;
      USE truss_;          USE platq_;          USE state_variables_;
      USE platt_;          USE motion_;         USE enumerated_sets_;
      USE spring_;         USE damper_;         USE displacement_bc_;
      USE stress_;         USE segment_;        USE contact_surface_;
      USE tied_bc_;        USE results_;        USE relink_scratch_;
      USE gauge1d_;        USE gauge2d_;        USE rigid_wall_bc_;
      USE gauge3d_;        USE massprop_;       USE include_file_;
      USE material_;       USE layering_;       USE sliding_node_;
      USE force_bc_;       USE node_set_;       USE contact_node_;
      USE nrbc_data_;      USE spring_bc_;      USE periodic_bc_;
      USE damper_bc_;      USE spot_weld_;      USE pressure_bc_;
      USE qa_record_;      USE plate_pair_;     USE segment_set_;
      USE body_force_;     USE section_2d_;     USE element_set_;
      USE section_1d_;     USE rigid_body_;     USE velocity_ic_;
      USE section_3d_;     USE extreme_value_;  USE centrifugal_force_;
      USE location_;       USE mean_stress_;    USE output_;         
      USE wedge_;          USE energy_flow_;    USE pyramid_;
      USE polyh_;

      USE precision_

      IMPLICIT    REAL(RTYPE) ( A-H, O-Z )
      IMPLICIT INTEGER(ITYPE) ( I-N      )
!!
!! Local variables.
      INTEGER(ITYPE), SAVE :: TStep = 0   ! Local count of time steps written
      INTEGER(ITYPE), SAVE :: MStep       ! Max number of time steps requested
      CHARACTER(8)         :: MEGO        ! Buffer for Material MatID output
      CHARACTER(5)         :: NEGO        ! Buffer for number-of-time-steps
      CHARACTER(80)        :: CBUFFER     ! EnSight Gold character const buffer
      CHARACTER(8),   SAVE :: ESG_ELEMENT_TYPE(11)
      INTEGER(ITYPE), SAVE :: MEL_COUNT(11)
      INTEGER(ITYPE), SAVE :: MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG
      INTEGER(ITYPE), SAVE :: NUMXL
      INTEGER(ITYPE)       :: NTUPLE(MXPHN) ! MXPHN comes from 'polyh_.f'
      INTEGER(ITYPE)       :: WedgeID
      LOGICAL              :: IOERROR
      LOGICAL,    EXTERNAL :: NEXT_NP_ID
      LOGICAL,    EXTERNAL :: NEXT_SEG_ID
      LOGICAL,        SAVE :: FIRST = .TRUE.
      LOGICAL,        SAVE :: INSERT_C_BINARY_QUE = .TRUE. 

      INTEGER(ITYPE), PARAMETER :: IBCOUNT = 1024

      REAL(RTYPE), DIMENSION(:), ALLOCATABLE, SAVE :: TIME_VALUES
!!
!! The following arrays are used to extract esg/vtk "parts" based on
!! material models.  For material "M," I = ELEMENTS_AND_NODES_USED(M) 
!! fills the following arrays with indecies.
!!
      INTEGER(ITYPE), DIMENSION(:), ALLOCATABLE, SAVE :: NELUSED ! Elems used by material M
      INTEGER(ITYPE), DIMENSION(:), ALLOCATABLE, SAVE :: NPTUSED ! Nodes used by material M
      INTEGER(ITYPE), DIMENSION(:), ALLOCATABLE, SAVE :: NPTNOWI ! Node's ESG output-index
!!
!! Contained functions
!!    INTEGER :: ELEMENTS_AND_NODES_USED ! Gens access arrays by material
!!    CHAR(8) :: ELEMENT_ESG_TYPE        ! Rtns EnSight Gold elem type
!!    INTEGER :: ELEMENT_N_TUPLE         ! Rtns connectivity count
!!    INTEGER :: MatID_DATA              ! Rtns material ID MatID
!!    INTEGER :: EleID_DATA              ! Rtns element ID  EleID
!!    REAL    :: STRESS_DATA             ! Rtns element stress
!!    REAL    :: BULK_STRAIN             ! Rtns element bulk strain
!!    REAL    :: STRAIN_ENERGY_DENSITY   ! Rtns element int energy
!!    REAL    :: PRESSURE                ! Rtns element pressure
!!    REAL    :: EFFECTIVE_STRESS        ! Rtns element eff stress
!!    REAL    :: MATERIAL_STATE          ! Rtns element material state
!!    INTEGER :: SEGMENTS_AND_NODES_USED ! Gens segment access arrays
!!    INTEGER :: SEGMENT_N_TUPLE         ! Rtns connectivity count
!!    INTEGER :: WEDGES_AND_NODES_USED   ! Gens wedge access arrays
!!    INTEGER :: WEDGE_N_TUPLE           ! Rtns connectivity count
!!    REAL    :: WEDGE_NP_COORDINATE     ! Rtns wedge np coordinates
!!    REAL    :: WEDGE_NP_VELOCITY       ! Rtns wedge np velocities
!!
!!
!!############################################################################
!! 1. CONSTANT FOR ALL-TIME DATA
!! "Create one-time geometry for subsequent ESG Data Files to reference."
!! This file will contain the undeformed configuration. It is generated the
!! first time this routine is called and contains the current velocities. 
!!
!! =========FIRST=============================================================
!! This file will contain the undeformed mesh with cell data for material 
!! index, initial physical volume, and initial critical time step, and with 
!! nodal point data for initial velocity conditions for examination.
!!
      IF (FIRST) THEN

#ifdef LANGUAGE_FORTRAN90
!!
!! ASSIGN is a PathScale pathf95 compiler unique procedure.
!! It causes the unformatted output to be Big_Endian (mips).
!! 
        CALL ASSIGN( "assign -N mips p:fmaego%", NERROR )
        IF (NERROR .NE. 0) THEN
          WRITE (MSG1,'(I8)') NERROR
          CALL USER_MESSAGE
     &      (
     &      MSGL//"FATAL"//
     &      MSGL//"WRITE_TO_ESG_RESULTS_FILE.999.99"//
     &      MSGL//"Call To pathf95 ASSIGN Failed."//
     &      MSGL//"Call PathScale. Error Number:"//MSG1
     &      )
        ENDIF
#endif
!!
!! Allocate storage for the number of time steps expected from the ESGFILE
!! input record entries. If the requested time increment is zero, limit the
!! number of time steps.
!!
      Tbgn = ESG_RESULTS_FILE%Begin
      Tinc = ESG_RESULTS_FILE%Delta
      Tend = MIN( TIME%Stop, ESG_RESULTS_FILE%End )

      IF (Tinc .GT. ZERO) THEN
        MStep = NINT( (Tend-Tbgn)/Tinc ) + 1
      ELSE
        MStep = 101
        WRITE (MSG1,'(I8)') MStep
        CALL USER_MESSAGE
     &    (
     &    MSGL//'INFORM'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.010.01'//
     &    MSGL//'Since "ESGFILE DELTA=0.0" Was Input And'//
     &    MSGL//'Due To The Nature Of The EnSight Gold Case'//
     &    MSGL//'File Format, The Number Of Time Steps Was'//
     &    MSGL//'Capped At'//MSG1//' Steps. Try A Non-Zero'//
     &    MSGL//'DELTA Value To Allow A Non-Zero Divide.'
     &    )
      ENDIF

      ALLOCATE ( TIME_VALUES(1:MStep) )
!!
!! Total finite element count. Note: Only one qudrilateral is produced for
!! each membrane (NUMM3, NUMM4) and shell (NUMP3, NUMP4) finite element.
!!
        NUMXL = NUMHX + NUMPX + NUMPY + NUMTX + NUMM3 + NUMP3 + NUMM4 + NUMP4 + NUMTR + NUMPH + NUMPG
!!
!! Count number of elements using each material model. Note: There may be more
!! material models specified then are actually referenced by elements in the
!! mesh. Hence, the effort to isolate materials that have a null usage count.
!!
        MXSNW = MAX ( NUMXL, NUMSG, NUMNP, NUMWX )
        ALLOCATE ( NELUSED(1:MXSNW), NPTUSED(1:NUMNP), NPTNOWI(1:NUMNP) )
!!
!! Define element and node counts for MATERIAL(*)%NElems and MATERIAL(*)%NNodes
!!
        DO M = 1,NUMMT
          I = ELEMENTS_AND_NODES_USED(M)
        ENDDO
!!
!! Initialize an indexable array used here for writing esg-element-types.
!!
        ESG_ELEMENT_TYPE( 1) = "hexa8"
        ESG_ELEMENT_TYPE( 2) = "penta6"
        ESG_ELEMENT_TYPE( 3) = "pyramid5"
        ESG_ELEMENT_TYPE( 4) = "tetra4"
        ESG_ELEMENT_TYPE( 5) = "tria3"
        ESG_ELEMENT_TYPE( 6) = "tria3"
        ESG_ELEMENT_TYPE( 7) = "quad4"
        ESG_ELEMENT_TYPE( 8) = "quad4"
        ESG_ELEMENT_TYPE( 9) = "bar2"
        ESG_ELEMENT_TYPE(10) = "nfaced"
        ESG_ELEMENT_TYPE(11) = "nsided"
!!
!! Open one-time geometry file fmaego.mesh.geom (EnSight Gold output).
!!
        IOERROR = .TRUE.
        OPEN
     &    (
     &    UNIT    =  IO_UNIT%LEGO,
     &    FILE    = 'fmaego.mesh.geom',
     &    STATUS  = 'UNKNOWN',
#ifdef _G95_
     &    FORM    = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &    FORM    = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &    FORM    = 'BINARY',                        !  Pathf95
#endif
     &    ERR     =  100
     &    )
        IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 100    IF (IOERROR) THEN
          CALL USER_MESSAGE
     &      (
     &      MSGL//'FATAL'//
     &      MSGL//'WRITE_TO_ESG_RESULTS_FILE.001.00'//
     &      MSGL//'Unable To Execute OPEN On: '//'fmaego.mesh.geom'
     &      )
        ELSE
!!
!! Initialize sequential Material part counter.
!!
          MPart = 0
!!
!! Initialize EnSight Gold C-binary data file (see FORM='BINARY' above).
!!
          CBUFFER = "C Binary"                         ! EnSight reader format que.
          WRITE (IO_UNIT%LEGO) CBUFFER
          CBUFFER = TRIM(JOB_ID_RECORD%CURRENT%TITLE)  ! User's job title record.
          WRITE (IO_UNIT%LEGO) CBUFFER
          CBUFFER = " "                                ! (unused subtitle)
          WRITE (IO_UNIT%LEGO) CBUFFER
          CBUFFER = "node id given"                    ! Expect to read nodal ID's
          WRITE (IO_UNIT%LEGO) CBUFFER
          CBUFFER = "element id given"                 ! Expect to read element ID's
          WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Turn individual material domains into esg/vtk "parts."
!!
          DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
            IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

              WRITE (MEGO,'(I8.8)') MATERIAL(M)%MatID

              NNodes = MATERIAL(M)%NNodes
              NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
              CBUFFER = "Part with MatID: "//MEGO
              WRITE (IO_UNIT%LEGO) CBUFFER
              CBUFFER = "coordinates"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) NNodes

              Nbgn = 1
              Nend = NREMAIN
              DO i = 1,NBLOCKS+1
                WRITE (IO_UNIT%LEGO) (NODE(NPTUSED(n))%ID, n = Nbgn,Nend)
                Nbgn = Nend + 1
                Nend = Nend + NBCOUNT
              ENDDO

              Nbgn = 1
              Nend = NREMAIN
              DO i = 1,NBLOCKS+1
                WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Px,KIND(0E0)), n = Nbgn,Nend)
                Nbgn = Nend + 1
                Nend = Nend + NBCOUNT
              ENDDO

              Nbgn = 1
              Nend = NREMAIN
              DO i = 1,NBLOCKS+1
                WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Py,KIND(0E0)), n = Nbgn,Nend)
                Nbgn = Nend + 1
                Nend = Nend + NBCOUNT
              ENDDO

              Nbgn = 1
              Nend = NREMAIN
              DO i = 1,NBLOCKS+1
                WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Pz,KIND(0E0)), n = Nbgn,Nend)
                Nbgn = Nend + 1
                Nend = Nend + NBCOUNT
              ENDDO

              Nend = 0
              LBCOUNT = IBCOUNT
              MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

              !! Predefined finite elements and their nodalizations.
              DO k = 1,9
                MXX = MEL_COUNT(k)
                IF (MXX .GT. 0) THEN
                  Nbgn = Nend + 1
                  Nend = Nend + MXX
                  CBUFFER = ESG_ELEMENT_TYPE(k)
                  WRITE (IO_UNIT%LEGO) CBUFFER
                  WRITE (IO_UNIT%LEGO) MXX

                  LBLOCKS = MXX / LBCOUNT
                  LREMAIN = MXX - LBCOUNT*LBLOCKS

                  Lbgn = Nbgn
                  Lend = Lbgn + LREMAIN - 1
                  DO i = 1,LBLOCKS+1
                    WRITE (IO_UNIT%LEGO) (EleID_DATA(NELUSED(n)), n = Lbgn,Lend)
                    Lbgn = Lend + 1
                    Lend = Lend + LBCOUNT
                  ENDDO

                  Lbgn = Nbgn
                  Lend = Lbgn + LREMAIN - 1
                  DO i = 1,LBLOCKS+1
                    WRITE (IO_UNIT%LEGO) (NTUPLE(1:ELEMENT_N_TUPLE(NELUSED(n))), n = Lbgn,Lend)
                    Lbgn = Lend + 1
                    Lend = Lend + LBCOUNT
                  ENDDO

                ENDIF
              ENDDO

              !!
              !! N-sided Polygons, if they use this material
              !!
              MXX = MEL_COUNT(10) !! Number of polyhedrons this material.
              IF (MXX .GT. 0) THEN
                !!
                !! Advance to loop on polyhedrons used by this material. 
                !!
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(10)
                WRITE (IO_UNIT%LEGO) CBUFFER
                WRITE (IO_UNIT%LEGO) MXX

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (EleID_DATA(NELUSED(n)), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO
                !!
                !! Polyhedron polygon-facet-counts (1:MPH)
                !!
                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (POLYHEDRON_FACET_COUNT(NELUSED(n)), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO
                !!
                !! Polygon nodal-point-counts NPGNP (1:NPHPG) per polyhedron (1:MPH) 
                !!
                DO n = Nbgn,Nend
                    WRITE (IO_UNIT%LEGO) NTUPLE(1:POLYGON_NP_COUNTERS(NELUSED(n)))
                ENDDO
                !!
                !! *ADVANCE* to a loop on polygons used by this material. 
                !!
                Nbgn = Nend + 1
                Nend = Nend + MEL_COUNT(11) !! Number of polygons this material, MPG.
                !!
                !! Polygon Nodal-point NPGNP-tuples per polygon (1:MPG)
                !!
                DO n = Nbgn,Nend
                  WRITE (IO_UNIT%LEGO) NTUPLE(1:ELEMENT_N_TUPLE(NELUSED(n))) 
                ENDDO
              ENDIF

            ENDIF
          ENDDO
        ENDIF
!!
!! =========SECOND============================================================
!! Write the individual node-sets segment-sets, and wedge-sets as separate 
!! parts. The idea is that EnSight and ParaView can load these parts one at 
!! a time "on top of" the above mesh file to confirm that the individual sets 
!! have been correctly specified. (The underlying mesh should first be colored 
!! "gray" so that the colors assigned by ParaView to each set standout and are 
!! more easily examined for correctness.)
!!
        DO M = 1,NUMNS
!!
!! Clear marker/index-sequence and translation arrays.
!!
          NELUSED = 0
          NPTUSED = 0
          NPTNOWI = 0
!!
!! Mark nodes used by this node set.
!!
          N = 0
          DO WHILE (NEXT_NP_ID(M,N))
            NPTUSED(N) = 1
          ENDDO
!!
!! Convert NPTUSED (and NELUSED) into a sequential index map.
!!
          K = 0
          DO N = 1,NUMNP
            IF (NPTUSED(N) .EQ. 1) THEN
              K = K + 1
              NPTNOWI(N) = K
              NPTUSED(K) = N
              NELUSED(K) = N
            ENDIF
          ENDDO
!!
!! For later use, record the number of nodes and elements (each node will be 
!! a "vertex element" for ParaView) that will be in the file for this NP_SET.
!!
          NNodes = K
          NElems = K

          IF (NNodes .GT. 0) THEN
            WRITE (MEGO,'(I8.8)') NODE_SET(M)%SetID
!!
!! Node Set: Nodal point ID's and coordinates.
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = Mpart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "Node Set ID: "//MEGO
            WRITE (IO_UNIT%LEGO) CBUFFER
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) NNodes
            WRITE (IO_UNIT%LEGO) (NODE(NPTUSED(n))%ID,                   n = 1,NNodes)
            WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Px,KIND(0E0)), n = 1,NNodes)
            WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Py,KIND(0E0)), n = 1,NNodes)
            WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Pz,KIND(0E0)), n = 1,NNodes)
!!
!! Node Set: Point-element inventory.
!!
            CBUFFER = "point"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) NElems  
            WRITE (IO_UNIT%LEGO) (NODE(NPTUSED(n))%ID, n = 1,NElems)
            WRITE (IO_UNIT%LEGO) (n,                   n = 1,NElems)
  
          ENDIF
        ENDDO
!!
!! Now, write segment set files.
!!
        DO M = 1,NUMSS
!!
!! The function SEGMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! segments-used array NELUSED and esg/esg output-index-translation array NPTNOWI for 
!! this segment set. If segment set M is not empty, the function returns 1, else 0. 
!!
          IF (SEGMENTS_AND_NODES_USED(M) .GT. 0) THEN
            WRITE (MEGO,'(I8.8)') SEGMENT_SET(M)%SetID
!!
!! Segment Set: Nodal point ID's and coordinates.
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = Mpart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "Segment Set ID: "//MEGO
            WRITE (IO_UNIT%LEGO) CBUFFER
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) NNodes
            WRITE (IO_UNIT%LEGO) (NODE(NPTUSED(n))%ID,                   n = 1,NNodes)
            WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Px,KIND(0E0)), n = 1,NNodes)
            WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Py,KIND(0E0)), n = 1,NNodes)
            WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Pz,KIND(0E0)), n = 1,NNodes)
!!
!! Segment Set: N-Sided polygon facet inventory for this segment set (SEGSET).
!!
            CBUFFER = "nsided"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) NElems  
            WRITE (IO_UNIT%LEGO) (SEGMENT(NELUSED(n))%PAR%SegID,         n = 1,NElems)
            WRITE (IO_UNIT%LEGO) (SEGMENT(NELUSED(n))%PAR%Knp,           n = 1,NElems)
            WRITE (IO_UNIT%LEGO) (NTUPLE(1:SEGMENT_N_TUPLE(NELUSED(n))), n = 1,NElems)
  
          ENDIF
        ENDDO
!!
!! Now, write wedge-set files. (Wedge-sets are generated internally based
!! on each specified solid-to-solid tied interface (Keyword: MPCON4).
!!
        WedgeID = 0
        DO M = 1,NUMC4
!!
!! The function WEDGES_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED,
!! wedges-used array NELUSED and esg output-index-translation array NPTNOWI for
!! this wedge set. If wedge set M is not empty, the function returns 1, else 0.
!!
          IF (WEDGES_AND_NODES_USED(M) .GT. 0) THEN
            WRITE (MEGO,'(I8.8)') SOLID_SOLID_INTERFACE(M)%MPCID
!!
!! Wedge Set: Nodal point ID's and coordinates.
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = Mpart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "Wedge Set ID: "//MEGO
            WRITE (IO_UNIT%LEGO) CBUFFER
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) NNodes
            WRITE (IO_UNIT%LEGO) (n,                                        n = 1,NNodes)
            WRITE (IO_UNIT%LEGO) (REAL(WEDGE_NP_COORDINATE(n,1),KIND(0E0)), n = 1,NNodes)  !  x-coordinate
            WRITE (IO_UNIT%LEGO) (REAL(WEDGE_NP_COORDINATE(n,2),KIND(0E0)), n = 1,NNodes)  !  y-coordinate
            WRITE (IO_UNIT%LEGO) (REAL(WEDGE_NP_COORDINATE(n,3),KIND(0E0)), n = 1,NNodes)  !  z-coordinate
!!
!! Wedge Set: Interface wedge set inventory for this tied-contact (MPCON4).
!!
            CBUFFER = "penta6"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) NElems  
            WRITE (IO_UNIT%LEGO) (WedgeID + n,                n = 1,NElems)
            WRITE (IO_UNIT%LEGO) (NTUPLE(1:WEDGE_N_TUPLE(n)), n = 1,NElems)

            WedgeID = WedgeID + NElems

          ENDIF
        ENDDO
!!
!! Now, write linked-pair-interface set. (All of the linked node-pairs 
!! will be written out as a single part.) (Keyword: MPCON5).
!!
!! The function LINKED_PAIR_NODES_USED(*) gen's nodal-points-used array NPTUSED,
!! linked-pairs-used array NELUSED and esg output-index-translation array NPTNOWI for
!! this linked-pair-interface set. If the linked-pair-interface set is not empty, 
!! the function returns 1, else 0.
!!
        M = 1  
        IF (LINKED_PAIR_NODES_USED(M) .GT. 0) THEN
!!
!! Use first linked-pair interface constraint ID as set ID.
!!
          WRITE (MEGO,'(I8.8)') LINKED_PAIR_INTERFACE(M)%MPCID 
!!
!!Linked-Pair Set: Nodal point ID's and coordinates.
!!
          CBUFFER = "part"
          WRITE (IO_UNIT%LEGO) CBUFFER
          MPart = Mpart + 1
          WRITE (IO_UNIT%LEGO) MPart
          CBUFFER = "Linked-Pair Set ID: "//MEGO
          WRITE (IO_UNIT%LEGO) CBUFFER
          CBUFFER = "coordinates"
          WRITE (IO_UNIT%LEGO) CBUFFER
          WRITE (IO_UNIT%LEGO) NNodes
          WRITE (IO_UNIT%LEGO) (NODE(NPTUSED(n))%ID,                   n = 1,NNodes)
          WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Px,KIND(0E0)), n = 1,NNodes)
          WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Py,KIND(0E0)), n = 1,NNodes)
          WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Pz,KIND(0E0)), n = 1,NNodes)
!!
!! Linked-Pair Set: Inventory for all linked-node-pair "dumb bells" (MPCON5).
!!
          CBUFFER = "bar2"
          WRITE (IO_UNIT%LEGO) CBUFFER
          WRITE (IO_UNIT%LEGO) NElems  
          WRITE (IO_UNIT%LEGO) (LINKED_PAIR_INTERFACE(NELUSED(n))%MPCID,   n = 1,NElems)
          WRITE (IO_UNIT%LEGO) (NTUPLE(1:LINKED_PAIR_N_TUPLE(NELUSED(n))), n = 1,NElems)

        ENDIF

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
!!
!! =========THIRD=============================================================
!! Write nodal point and cell (element) data.
!! 
!! Open an ESG nodal point velocity initial condition file.
!!
        IOERROR = .TRUE.
        OPEN
     &    (
     &    UNIT    =  IO_UNIT%LEGO,
     &    FILE    = 'fmaego.mesh.vics',
     &    STATUS  = 'UNKNOWN',
#ifdef _G95_
     &    FORM    = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &    FORM    = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &    FORM    = 'BINARY',                        !  Pathf95
#endif
     &    ERR     =  200
     &    )
        IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 200    IF (IOERROR) THEN
          CALL USER_MESSAGE
     &      (
     &      MSGL//'FATAL'//
     &      MSGL//'WRITE_TO_ESG_RESULTS_FILE.002.00'//
     &      MSGL//'Unable To Execute OPEN On: '//'fmaego.mesh.vics'
     &      )
        ELSE
!!
!! Initialize sequential Material part counter.
!!
          MPart = 0
!!
!! Initialize EnSight Gold static variable results file
!!
          CBUFFER = "Velocity Initial Conditions"
          WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
          DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
            IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

              NNodes = MATERIAL(M)%NNodes
              NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = MPart + 1
              WRITE (IO_UNIT%LEGO) MPart
              CBUFFER = "coordinates"
              WRITE (IO_UNIT%LEGO) CBUFFER

              Nbgn = 1
              Nend = NREMAIN
              DO i = 1,NBLOCKS+1
                WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vx,KIND(0E0)), n = Nbgn,Nend)
                Nbgn = Nend + 1
                Nend = Nend + NBCOUNT
              ENDDO

              Nbgn = 1
              Nend = NREMAIN
              DO i = 1,NBLOCKS+1
                WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vy,KIND(0E0)), n = Nbgn,Nend)
                Nbgn = Nend + 1
                Nend = Nend + NBCOUNT
              ENDDO

              Nbgn = 1
              Nend = NREMAIN
              DO i = 1,NBLOCKS+1
                WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vz,KIND(0E0)), n = Nbgn,Nend)
                Nbgn = Nend + 1
                Nend = Nend + NBCOUNT
              ENDDO

            ENDIF
          ENDDO
!!
!! Assign initial velocity values to node sets.
!!
          DO M = 1,NUMNS
!!
!! Clear marker/index-sequence and translation arrays.
!!
            NELUSED = 0
            NPTUSED = 0
            NPTNOWI = 0
!!
!! Mark nodes used by this node set.
!!
            N = 0
            DO WHILE (NEXT_NP_ID(M,N))
              NPTUSED(N) = 1
            ENDDO
!!
!! Convert NPTUSED (and NELUSED) into a sequential index map.
!!
            K = 0
            DO N = 1,NUMNP
              IF (NPTUSED(N) .EQ. 1) THEN
                K = K + 1
                NPTNOWI(N) = K
                NPTUSED(K) = N
                NELUSED(K) = N
              ENDIF
            ENDDO
!!
!! For later use, record the number of nodes and elements (each node will be 
!! a "vertex element" for ParaView) that will be in the file for this NP_SET.
!!
            NNodes = K
            NElems = K

            IF (NNodes .GT. 0) THEN
!!
!! Node Set: Velocity initial conditions.
!!
              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = MPart + 1
              WRITE (IO_UNIT%LEGO) MPart
              CBUFFER = "coordinates"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vx,KIND(0E0)), n = 1,NNodes)
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vy,KIND(0E0)), n = 1,NNodes)
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vz,KIND(0E0)), n = 1,NNodes)

            ENDIF
          ENDDO
!!
!! Assign initial velocity values to segment sets.
!!
          DO M = 1,NUMSS
!!
!! The function SEGMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! segments-used array NELUSED and esg/esg output-index-translation array NPTNOWI for 
!! this segment set. If segment set M is not empty, the function returns 1, else 0. 
!!
            IF (SEGMENTS_AND_NODES_USED(M) .GT. 0) THEN
!!
!! Segment Set: Velocity initial conditions.
!!
              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
              CBUFFER = "coordinates"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vx,KIND(0E0)), n = 1,NNodes)
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vy,KIND(0E0)), n = 1,NNodes)
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vz,KIND(0E0)), n = 1,NNodes)

            ENDIF
          ENDDO
!!
!! Assign initial velocity values to wedge sets.
!!
          DO M = 1,NUMC4
!!
!! The function WEDGES_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED,
!! wedges-used array NELUSED and esg output-index-translation array NPTNOWI for
!! this wedge set. If wedge set M is not empty, the function returns 1, else 0.
!!
            IF (WEDGES_AND_NODES_USED(M) .GT. 0) THEN
!!
!! Wedge Set: Velocity initial conditions.
!!
              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
              CBUFFER = "coordinates"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(WEDGE_NP_VELOCITY(n,1),KIND(0E0)), n = 1,NNodes)  !  x-coordinate
              WRITE (IO_UNIT%LEGO) (REAL(WEDGE_NP_VELOCITY(n,2),KIND(0E0)), n = 1,NNodes)  !  y-coordinate
              WRITE (IO_UNIT%LEGO) (REAL(WEDGE_NP_VELOCITY(n,3),KIND(0E0)), n = 1,NNodes)  !  z-coordinate

            ENDIF
          ENDDO
!!
!! Assign initial velocity values to linked-pair sets.
!!
!! The function LINKED_PAIR_NODES_USED(*) gen's nodal-points-used array 
!! NPTUSED, linked-pairs-used array NELUSED and esg output-index-translation 
!! array NPTNOWI for this linked-pair set. If linked-pair set M is not empty, 
!! the function returns 1, else 0.
!!
          M = 1
          IF (LINKED_PAIR_NODES_USED(M) .GT. 0) THEN
!!
!! Linked-Pair Set: Velocity initial conditions.
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = Mpart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vx,KIND(0E0)), n = 1,NNodes)
            WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vy,KIND(0E0)), n = 1,NNodes)
            WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vz,KIND(0E0)), n = 1,NNodes)

          ENDIF

          CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
        ENDIF
!! 
!! Open an ESG element material number file.
!!
        IOERROR = .TRUE.
        OPEN
     &    (
     &    UNIT    =  IO_UNIT%LEGO,
     &    FILE    = 'fmaego.mesh.mats',
     &    STATUS  = 'UNKNOWN',
#ifdef _G95_
     &    FORM    = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &    FORM    = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &    FORM    = 'BINARY',                        !  Pathf95
#endif
     &    ERR     =  301
     &    )
        IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 301    IF (IOERROR) THEN
          CALL USER_MESSAGE
     &      (
     &      MSGL//'FATAL'//
     &      MSGL//'WRITE_TO_ESG_RESULTS_FILE.003.00'//
     &      MSGL//'Unable To Execute OPEN On: '//'fmaego.mesh.mats'
     &      )
        ELSE
!!
!! Initialize sequential part counter.
!!
          MPart = 0
!!
!! Initialize EnSight Gold static variable results file.
!!
          CBUFFER = "Cell Material Color Index"
          WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
          DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
            IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

              NNodes = MATERIAL(M)%NNodes
              NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = MPart + 1
              WRITE (IO_UNIT%LEGO) MPart

              Nend = 0
              LBCOUNT = IBCOUNT
              MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

              DO k = 1,10
                MXX = MEL_COUNT(k)
                IF (MXX .GT. 0) THEN
                  Nbgn = Nend + 1
                  Nend = Nend + MXX
                  CBUFFER = ESG_ELEMENT_TYPE(k)         
                  WRITE (IO_UNIT%LEGO) CBUFFER

                  LBLOCKS = MXX / LBCOUNT
                  LREMAIN = MXX - LBCOUNT*LBLOCKS

                  Lbgn = Nbgn
                  Lend = Lbgn + LREMAIN - 1
                  DO i = 1,LBLOCKS+1
                    WRITE (IO_UNIT%LEGO) (REAL(M,KIND(0E0)), n = Lbgn,Lend)
                    Lbgn = Lend + 1
                    Lend = Lend + LBCOUNT
                  ENDDO

                ENDIF
              ENDDO

            ENDIF
          ENDDO
!!
!! Assign "MatID" values to node sets.
!!
          DO M = 1,NUMNS
!!
!! Clear marker/index-sequence and translation arrays.
!!
            NELUSED = 0
            NPTUSED = 0
            NPTNOWI = 0
!!
!! Mark nodes used by this node set.
!!
            N = 0
            DO WHILE (NEXT_NP_ID(M,N))
              NPTUSED(N) = 1
            ENDDO
!!
!! Convert NPTUSED (and NELUSED) into a sequential index map.
!!
            K = 0
            DO N = 1,NUMNP
              IF (NPTUSED(N) .EQ. 1) THEN
                K = K + 1
                NPTNOWI(N) = K
                NPTUSED(K) = N
                NELUSED(K) = N
              ENDIF
            ENDDO
!!
!! For later use, record the number of nodes and elements (each node will be 
!! a "vertex element" for ParaView) that will be in the file for this NP_SET.
!!
            NNodes = K
            NElems = K

            IF (NNodes .GT. 0) THEN

              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
!!
!! Node Set: Assign internal node set ID as "MatID."
!!
              CBUFFER = "point"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(NUMMT+M,KIND(0E0)), n = 1,NElems)

            ENDIF
          ENDDO
!!
!! Assign "MatID" values to segment sets.
!!
          DO M = 1,NUMSS
!!
!! The function SEGMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! segments-used array NELUSED and esg/esg output-index-translation array NPTNOWI for 
!! this segment set. If segment set M is not empty, the function returns 1, else 0. 
!!
            IF (SEGMENTS_AND_NODES_USED(M) .GT. 0) THEN

              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
!!
!! Segment Set: Assign internal segment set ID as "MatID."
!!
              CBUFFER = "nsided"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(NUMMT+NUMNS+M,KIND(0E0)), n = 1,NElems)

            ENDIF
          ENDDO
!!
!! Assign "MatID" values to wedge sets.
!!
          DO M = 1,NUMC4
!!
!! The function WEDGES_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED,
!! wedges-used array NELUSED and esg output-index-translation array NPTNOWI for
!! this wedge set. If wedge set M is not empty, the function returns 1, else 0.
!!
            IF (WEDGES_AND_NODES_USED(M) .GT. 0) THEN

              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
!!
!! Wedge Set: Assign internal wedge set ID as "MatID."
!!
              CBUFFER = "penta6"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(NUMMT+NUMNS+NUMSS+M,KIND(0E0)), n = 1,NElems)

            ENDIF
          ENDDO
!!
!! Assign "MatID" values to linked-pair sets.
!!
!! The function LINKED_PAIR_NODES_USED(*) gen's nodal-points-used array 
!! NPTUSED, linked-pairs-used array NELUSED and esg output-index-translation 
!! array NPTNOWI for this linked-pair set. If linked-pair set M is not empty, 
!! the function returns 1, else 0.
!!
          M = 1
          IF (LINKED_PAIR_NODES_USED(M) .GT. 0) THEN

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = Mpart + 1
            WRITE (IO_UNIT%LEGO) MPart
!!
!! Linked-Pair Set: Assign internal linked-pair set ID as "MatID."
!!
            CBUFFER = "bar2"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) (REAL(NUMMT+NUMNS+NUMSS+NUMC4+1,KIND(0E0)), n = 1,NElems)

          ENDIF

          CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
        ENDIF
!! 
!! Open an ESG element initial volume file.
!!
        IOERROR = .TRUE.
        OPEN
     &    (
     &    UNIT    =  IO_UNIT%LEGO,
     &    FILE    = 'fmaego.mesh.vols',
     &    STATUS  = 'UNKNOWN',
#ifdef _G95_
     &    FORM    = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &    FORM    = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &    FORM    = 'BINARY',                        !  Pathf95
#endif
     &    ERR     =  302
     &    )
        IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 302    IF (IOERROR) THEN
          CALL USER_MESSAGE
     &      (
     &      MSGL//'FATAL'//
     &      MSGL//'WRITE_TO_ESG_RESULTS_FILE.003.00'//
     &      MSGL//'Unable To Execute OPEN On: '//'fmaego.mesh.vols'
     &      )
        ELSE
!!
!! Initialize sequential part counter.
!!
          MPart = 0
!!
!! Initialize EnSight Gold static variable results file.
!!
          CBUFFER = "Cell Initial Volume Value"
          WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
          DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
            IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

              NNodes = MATERIAL(M)%NNodes
              NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = MPart + 1
              WRITE (IO_UNIT%LEGO) MPart

              Nend = 0
              LBCOUNT = IBCOUNT
              MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

              DO k = 1,10
                MXX = MEL_COUNT(k)
                IF (MXX .GT. 0) THEN
                  Nbgn = Nend + 1
                  Nend = Nend + MXX
                  CBUFFER = ESG_ELEMENT_TYPE(k)         
                  WRITE (IO_UNIT%LEGO) CBUFFER

                  LBLOCKS = MXX / LBCOUNT
                  LREMAIN = MXX - LBCOUNT*LBLOCKS

                  Lbgn = Nbgn
                  Lend = Lbgn + LREMAIN - 1
                  DO i = 1,LBLOCKS+1
                    WRITE (IO_UNIT%LEGO) (ELEMENT_VOLUME(NELUSED(n)), n = Lbgn,Lend)
                    Lbgn = Lend + 1
                    Lend = Lend + LBCOUNT
                  ENDDO

                ENDIF
              ENDDO

            ENDIF
          ENDDO
!!
!! Assign "Volume" values to node sets.
!!
          DO M = 1,NUMNS
!!
!! Clear marker/index-sequence and translation arrays.
!!
            NELUSED = 0
            NPTUSED = 0
            NPTNOWI = 0
!!
!! Mark nodes used by this node set.
!!
            N = 0
            DO WHILE (NEXT_NP_ID(M,N))
              NPTUSED(N) = 1
            ENDDO
!!
!! Convert NPTUSED (and NELUSED) into a sequential index map.
!!
            K = 0
            DO N = 1,NUMNP
              IF (NPTUSED(N) .EQ. 1) THEN
                K = K + 1
                NPTNOWI(N) = K
                NPTUSED(K) = N
                NELUSED(K) = N
              ENDIF
            ENDDO
!!
!! For later use, record the number of nodes and elements (each node will be 
!! a "vertex element" for ParaView) that will be in the file for this NP_SET.
!!
            NNodes = K
            NElems = K

            IF (NNodes .GT. 0) THEN

              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
!!
!! Node Set: Assign a null value for nodal point volume.
!!
              CBUFFER = "point"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(ZERO,KIND(0E0)), n = 1,NElems)

            ENDIF
          ENDDO
!!
!! Assign "Volume" values to segment sets.
!!
          DO M = 1,NUMSS
!!
!! The function SEGMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! segments-used array NELUSED and esg/esg output-index-translation array NPTNOWI for 
!! this segment set. If segment set M is not empty, the function returns 1, else 0. 
!!
            IF (SEGMENTS_AND_NODES_USED(M) .GT. 0) THEN

              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
!!
!! Segment Set: Assign "volume" value zero (0.0) to all segments in this set. 
!!
              CBUFFER = "nsided"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(ZERO,KIND(0E0)), n = 1,NElems)

            ENDIF
          ENDDO
!!
!! Assign "volume" values to wedge sets (can be negative, null, positive).
!!
          DO M = 1,NUMC4
!!
!! The function WEDGES_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED,
!! wedges-used array NELUSED and esg output-index-translation array NPTNOWI for
!! this wedge set. If wedge set M is not empty, the function returns 1, else 0.
!!
            IF (WEDGES_AND_NODES_USED(M) .GT. 0) THEN

              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
!!
!! Write this wedge's volume.
!!
              CBUFFER = "penta6"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (WEDGE_VOLUME(NELUSED(n)), n = 1,NElems)

            ENDIF
          ENDDO
!!
!! Assign "volume" values to linked-pair sets (distance apart).
!!
!! The function LINKED_PAIR_NODES_USED(*) gen's nodal-points-used array 
!! NPTUSED, linked-pairs-used array NELUSED and esg output-index-translation 
!! array NPTNOWI for this linked-pair set. If linked-pair set M is not empty, 
!! the function returns 1, else 0.
!!
          M = 1
          IF (LINKED_PAIR_NODES_USED(M) .GT. 0) THEN

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = Mpart + 1
            WRITE (IO_UNIT%LEGO) MPart
!!
!! Write this linked-pair set's volume (distance apart).
!!
            CBUFFER = "bar2"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) (LINKED_PAIR_VOLUME(NELUSED(n)), n = 1,NElems)

          ENDIF

          CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
        ENDIF
!! 
!! Open an ESG element initial critical time step file.
!!
        IOERROR = .TRUE.
        OPEN
     &    (
     &    UNIT    =  IO_UNIT%LEGO,
     &    FILE    = 'fmaego.mesh.cdts',
     &    STATUS  = 'UNKNOWN',
#ifdef _G95_
     &    FORM    = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &    FORM    = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &    FORM    = 'BINARY',                        !  Pathf95
#endif
     &    ERR     =  303
     &    )
        IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 303    IF (IOERROR) THEN
          CALL USER_MESSAGE
     &      (
     &      MSGL//'FATAL'//
     &      MSGL//'WRITE_TO_ESG_RESULTS_FILE.003.00'//
     &      MSGL//'Unable To Execute OPEN On: '//'fmaego.mesh.cdts'
     &      )
        ELSE
!!
!! Initialize sequential part counter.
!!
          MPart = 0
!!
!! Initialize EnSight Gold static variable results file.
!!
          CBUFFER = "Cell Initial Critical Dt"
          WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
          DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
            IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

              NNodes = MATERIAL(M)%NNodes
              NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = MPart + 1
              WRITE (IO_UNIT%LEGO) MPart

              Nend = 0
              LBCOUNT = IBCOUNT
              MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

              DO k = 1,10
                MXX = MEL_COUNT(k)
                IF (MXX .GT. 0) THEN
                  Nbgn = Nend + 1
                  Nend = Nend + MXX
                  CBUFFER = ESG_ELEMENT_TYPE(k)         
                  WRITE (IO_UNIT%LEGO) CBUFFER

                  LBLOCKS = MXX / LBCOUNT
                  LREMAIN = MXX - LBCOUNT*LBLOCKS

                  Lbgn = Nbgn
                  Lend = Lbgn + LREMAIN - 1
                  DO i = 1,LBLOCKS+1
                    WRITE (IO_UNIT%LEGO) (ELEMENT_CRITICAL_DT(NELUSED(n)), n = Lbgn,Lend)
                    Lbgn = Lend + 1
                    Lend = Lend + LBCOUNT
                  ENDDO

                ENDIF
              ENDDO

            ENDIF
          ENDDO
!!
!! Assign "Critical-Dt" values to node sets.
!!
          DO M = 1,NUMNS
!!
!! Clear marker/index-sequence and translation arrays.
!!
            NELUSED = 0
            NPTUSED = 0
            NPTNOWI = 0
!!
!! Mark nodes used by this node set.
!!
            N = 0
            DO WHILE (NEXT_NP_ID(M,N))
              NPTUSED(N) = 1
            ENDDO
!!
!! Convert NPTUSED (and NELUSED) into a sequential index map.
!!
            K = 0
            DO N = 1,NUMNP
              IF (NPTUSED(N) .EQ. 1) THEN
                K = K + 1
                NPTNOWI(N) = K
                NPTUSED(K) = N
                NELUSED(K) = N
              ENDIF
            ENDDO
!!
!! For later use, record the number of nodes and elements (each node will be 
!! a "vertex element" for ParaView) that will be in the file for this NP_SET.
!!
            NNodes = K
            NElems = K

            IF (NNodes .GT. 0) THEN

              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
!!
!! Node Set: Assign a null value for nodal point volume.
!!
              CBUFFER = "point"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(ZERO,KIND(0E0)), n = 1,NElems)

            ENDIF
          ENDDO
!!
!! Assign "Critical-Dt" values to segment sets.
!!
          DO M = 1,NUMSS
!!
!! The function SEGMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! segments-used array NELUSED and esg/esg output-index-translation array NPTNOWI for 
!! this segment set. If segment set M is not empty, the function returns 1, else 0. 
!!
            IF (SEGMENTS_AND_NODES_USED(M) .GT. 0) THEN

              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
!!
!! Segment Set: Assign internal segment set ID as "MatID."
!!
              CBUFFER = "nsided"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (REAL(ZERO,KIND(0E0)), n = 1,NElems)

            ENDIF
          ENDDO
!!
!! Assign "Critical-Dt" values to wedge sets; use value of parent element.
!!
          DO M = 1,NUMC4
!!
!! The function WEDGES_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED,
!! wedges-used array NELUSED and esg output-index-translation array NPTNOWI for
!! this wedge set. If wedge set M is not empty, the function returns 1, else 0.
!!
            IF (WEDGES_AND_NODES_USED(M) .GT. 0) THEN

              CBUFFER = "part"
              WRITE (IO_UNIT%LEGO) CBUFFER
              MPart = Mpart + 1
              WRITE (IO_UNIT%LEGO) MPart
!!
!! Write this wedge's volume.
!!
              CBUFFER = "penta6"
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) (WEDGE_CRITICAL_DT(NELUSED(n)), n = 1,NElems)

            ENDIF
          ENDDO
!!
!! Assign "Critical-Dt" values to linked-pair sets; use zero.
!!
!! The function LINKED_PAIR_NODES_USED(*) gen's nodal-points-used array 
!! NPTUSED, linked-pairs-used array NELUSED and esg output-index-translation 
!! array NPTNOWI for this linked-pair set. If linked-pair set M is not empty, 
!! the function returns 1, else 0.
!!
          M = 1
          IF (LINKED_PAIR_NODES_USED(M) .GT. 0) THEN

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = Mpart + 1
            WRITE (IO_UNIT%LEGO) MPart
!!
!! Write this linked-pairs's critical-dt.
!!
            CBUFFER = "bar2"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) (LINKED_PAIR_CRITICAL_DT(NELUSED(n)), n = 1,NElems)

          ENDIF

          CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
        ENDIF
!!
!! =========FOURTH============================================================
!! Write a *.case file to allow EnSight and ParaView to identify the mesh, 
!! node-set, side-set and wedge-set files as distinct "Parts."
!! 
!! Open an ESG case-file.
!!
        IOERROR = .TRUE.
        OPEN
     &    (
     &    UNIT    =  IO_UNIT%LEGO,
     &    FILE    = 'fmaego.mesh.case',
     &    STATUS  = 'UNKNOWN',
     &    FORM    = 'FORMATTED',
     &    ERR     =  400
     &    )
        IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 400    IF (IOERROR) THEN
          CALL USER_MESSAGE
     &      (
     &      MSGL//'FATAL'//
     &      MSGL//'WRITE_TO_ESG_RESULTS_FILE.004.00'//
     &      MSGL//'Unable To Execute OPEN On: '//'fmaego.mesh.case'
     &      )
        ELSE

          WRITE (IO_UNIT%LEGO,'(A)') "FORMAT"
          WRITE (IO_UNIT%LEGO,'(A)') "type: ensight gold"

          WRITE (IO_UNIT%LEGO,'(A)') "GEOMETRY"
          WRITE (IO_UNIT%LEGO,'(A)') "model: fmaego.mesh.geom"

          WRITE (IO_UNIT%LEGO,'(A)') "VARIABLE"
          WRITE (IO_UNIT%LEGO,'(A)') "vector per node:    Velocity-IC fmaego.mesh.vics"
          WRITE (IO_UNIT%LEGO,'(A)') "scalar per element: Material-No fmaego.mesh.mats"
          WRITE (IO_UNIT%LEGO,'(A)') "scalar per element: Initial-Vol fmaego.mesh.vols"
          WRITE (IO_UNIT%LEGO,'(A)') "scalar per element: Critical-Dt fmaego.mesh.cdts"

          CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
        ENDIF

        FIRST = .FALSE.
      ENDIF
!!
!!############################################################################
!! PART 2. RESULTS FOR THIS TIME STEP.
!! Initialize and append nodal point and cell (element) data.
!!
!! Increment EnSight Gold number-of-time-steps counter and store current time.
!!
      TStep = TStep + 1
      IF (TStep .LE. MStep) THEN
        WRITE (NEGO,'(I5)') TStep
        TIME_VALUES(TStep) = TIME%Total
      ELSE
        WRITE (MSG1,'(I8)') MStep
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.010.02'//
     &    MSGL//'TStep, The Current Write Counter, Has Exceeded MStep.'//
     &    MSGL//'TStep:'//MSG1//
     &    MSGL//'MStep:'//MSG2
     &    )
        RETURN    
      ENDIF
!!
!! =========FIRST=============================================================
!! Open geometry file fmaego.data.geom (EnSight Gold output) for nodal point 
!! and cell results to reference.
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.geom',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  500
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 500  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.005.00'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.geom'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Insert EnSight Gold C-binary data file que (see FORM='BINARY' above).
!!
        IF (INSERT_C_BINARY_QUE) THEN
          CBUFFER = "C Binary"                         ! EnSight reader format que.
          WRITE (IO_UNIT%LEGO) CBUFFER
          INSERT_C_BINARY_QUE = .FALSE.
        ENDIF
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"                  ! "Time serialized" que.
        WRITE (IO_UNIT%LEGO) CBUFFER
        CBUFFER = TRIM(JOB_ID_RECORD%CURRENT%TITLE)  ! User's job title record.
        WRITE (IO_UNIT%LEGO) CBUFFER
        CBUFFER = " "                                ! (unused subtitle)
        WRITE (IO_UNIT%LEGO) CBUFFER
        CBUFFER = "node id given"                    ! Expect to read nodal ID's
        WRITE (IO_UNIT%LEGO) CBUFFER
        CBUFFER = "element id given"                 ! Expect to read element ID's
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Turn individual material domains into esg/vtk "parts."
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            WRITE (MEGO,'(I8.8)') MATERIAL(M)%MatID

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = Mpart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "Part with MatID: "//MEGO
            WRITE (IO_UNIT%LEGO) CBUFFER
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER
            WRITE (IO_UNIT%LEGO) NNodes

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (NODE(NPTUSED(n))%ID, n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Px+MOTION(NPTUSED(n))%Ux,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Py+MOTION(NPTUSED(n))%Uy,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Pz+MOTION(NPTUSED(n))%Uz,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nend = 0
            LBCOUNT = IBCOUNT
            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

            !! Predefined finite elements and their nodalizations.
            DO k = 1,9
              MXX = MEL_COUNT(k)
              IF (MXX .GT. 0) THEN
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(k)         
                WRITE (IO_UNIT%LEGO) CBUFFER
                WRITE (IO_UNIT%LEGO) MXX

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (EleID_DATA(NELUSED(n)), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (NTUPLE(1:ELEMENT_N_TUPLE(NELUSED(n))), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

              ENDIF
            ENDDO

            !!
            !! N-sided Polygons, if they use this material
            !!
            MXX = MEL_COUNT(10) !! Number of polyhedrons this material.
            IF (MXX .GT. 0) THEN
              !!
              !! Advance to loop on polyhedrons used by this material. 
              !!
              Nbgn = Nend + 1
              Nend = Nend + MXX
              CBUFFER = ESG_ELEMENT_TYPE(10)
              WRITE (IO_UNIT%LEGO) CBUFFER
              WRITE (IO_UNIT%LEGO) MXX

              LBLOCKS = MXX / LBCOUNT
              LREMAIN = MXX - LBCOUNT*LBLOCKS

              Lbgn = Nbgn
              Lend = Lbgn + LREMAIN - 1
              DO i = 1,LBLOCKS+1
                WRITE (IO_UNIT%LEGO) (EleID_DATA(NELUSED(n)), n = Lbgn,Lend)
                Lbgn = Lend + 1
                Lend = Lend + LBCOUNT
              ENDDO
              !!
              !! Polyhedron polygon-facet-counts (1:MPH)
              !!
              Lbgn = Nbgn
              Lend = Lbgn + LREMAIN - 1
              DO i = 1,LBLOCKS+1
                WRITE (IO_UNIT%LEGO) (POLYHEDRON_FACET_COUNT(NELUSED(n)), n = Lbgn,Lend)
                Lbgn = Lend + 1
                Lend = Lend + LBCOUNT
              ENDDO
              !!
              !! Polygon nodal-point-counts NPGNP (1:NPHPG) per polyhedron (1:MPH) 
              !!
              DO n = Nbgn,Nend
                  WRITE (IO_UNIT%LEGO) NTUPLE(1:POLYGON_NP_COUNTERS(NELUSED(n)))
              ENDDO
              !!
              !! *ADVANCE* to a loop on polygons used by this material. 
              !!
              Nbgn = Nend + 1
              Nend = Nend + MEL_COUNT(11) !! Number of polygons this material, MPG.
              !!
              !! Polygon Nodal-point NPGNP-tuples per polygon (1:MPG)
              !!
              DO n = Nbgn,Nend
                WRITE (IO_UNIT%LEGO) NTUPLE(1:ELEMENT_N_TUPLE(NELUSED(n))) 
              ENDDO
            ENDIF

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
      ENDIF
!!
!! =========SECOND============================================================
!! Open and append ESG nodal point results to data files.
!!
!! Displacements...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.ndis',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  601
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 601  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.006.01'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.ndis'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Nodal Point Displacement Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Ux,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Uy,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Uz,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Velocities...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.nvel',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  602
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 602  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.006.02'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.nvel'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Nodal Point Velocity Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vx,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vy,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Vz,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Accelerations...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.nacc',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  603
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 603  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.006.03'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.nacc'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Nodal Point Acceleration Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Ax,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Ay,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(MOTION(NPTUSED(n))%Az,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Velocity Maximums (seen up to now)...
!!
      IF (NUMVX .GT. 0) THEN

      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.nvmx',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  604
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 604  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.006.03'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.vmx'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Nodal Point Velocity Maxima Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(VELOCITY_EXTREMA(NPTUSED(n))%Vx_max,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(VELOCITY_EXTREMA(NPTUSED(n))%Vy_max,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(VELOCITY_EXTREMA(NPTUSED(n))%Vz_max,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Velocity Minimums (seen up to now)...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.nvmn',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  605
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 605  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.006.03'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.nvmn'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Nodal Point Velocity Minima Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems
!!
!!***************************************************************************
!! To avoid any output buffer et cetera limitations, break output writes
!! into smaller batches of IBCOUNT each.  (Overflows have occurred already.)
!! Element(Cell) Blocking: L... counter items.
!! Nodal(Points) Blocking: N... counter items.
!!
      LBCOUNT=IBCOUNT;  LBLOCKS=NElems/LBCOUNT;  LREMAIN=NElems-LBCOUNT*LBLOCKS

      NBCOUNT=IBCOUNT;  NBLOCKS=NNodes/NBCOUNT;  NREMAIN=NNodes-NBCOUNT*NBLOCKS
!!
!!***************************************************************************
!!
            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart
            CBUFFER = "coordinates"
            WRITE (IO_UNIT%LEGO) CBUFFER

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(VELOCITY_EXTREMA(NPTUSED(n))%Vx_min,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(VELOCITY_EXTREMA(NPTUSED(n))%Vy_min,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

            Nbgn = 1
            Nend = NREMAIN
            DO i = 1,NBLOCKS+1
              WRITE (IO_UNIT%LEGO) (REAL(VELOCITY_EXTREMA(NPTUSED(n))%Vz_min,KIND(0E0)), n = Nbgn,Nend)
              Nbgn = Nend + 1
              Nend = Nend + NBCOUNT
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! End of if-block to see if velocity extrame data exists.
!!
      ENDIF
!!
!! =========THIRD=============================================================
!! Open and append ESG element (cell) results to data files.
!!
!! Internal Material Number...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.emat',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  701
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 701  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.01'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.emat'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Element Internal Material Index"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart

            Nend = 0
            LBCOUNT = IBCOUNT
            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

            DO k = 1,10
              MXX = MEL_COUNT(k)
              IF (MXX .GT. 0) THEN
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(k)         
                WRITE (IO_UNIT%LEGO) CBUFFER

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (REAL(M,KIND(0E0)), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

              ENDIF
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Internal Material Number...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.esta',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  711
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 711  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.01'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.esta'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Element Material Current State"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart

            Nend = 0
            LBCOUNT = IBCOUNT
            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

            DO k = 1,10
              MXX = MEL_COUNT(k)
              IF (MXX .GT. 0) THEN
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(k)         
                WRITE (IO_UNIT%LEGO) CBUFFER

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (MATERIAL_STATE(NELUSED(n)), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

              ENDIF
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Stress Flux (Stress*Velocity)...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.estv',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  712
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 712  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.02'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.estv'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Element Stress*Velocity Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart

            Nend = 0
            LBCOUNT = IBCOUNT
            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

            DO k = 1,10
              MXX = MEL_COUNT(k)
              IF (MXX .GT. 0) THEN
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(k)         
                WRITE (IO_UNIT%LEGO) CBUFFER

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRESS_FLUX(NELUSED(n),1), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRESS_FLUX(NELUSED(n),2), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRESS_FLUX(NELUSED(n),3), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

              ENDIF
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Stress...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.estr',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  702
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 702  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.02'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.estr'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Element Stress(xx,yy,zz,xy,xz,yz) Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart

            Nend = 0
            LBCOUNT = IBCOUNT
            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

            DO k = 1,10
              MXX = MEL_COUNT(k)
              IF (MXX .GT. 0) THEN
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(k)         
                WRITE (IO_UNIT%LEGO) CBUFFER

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRESS_DATA(NELUSED(n),1), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRESS_DATA(NELUSED(n),2), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRESS_DATA(NELUSED(n),3), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRESS_DATA(NELUSED(n),4), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRESS_DATA(NELUSED(n),5), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRESS_DATA(NELUSED(n),6), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

              ENDIF
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Bulk Strain...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.elnv',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  703
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 703  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.03'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.elnv'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Element Bulk Strain Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart

            Nend = 0
            LBCOUNT = IBCOUNT
            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

            DO k = 1,10
              MXX = MEL_COUNT(k)
              IF (MXX .GT. 0) THEN
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(k)         
                WRITE (IO_UNIT%LEGO) CBUFFER

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (BULK_STRAIN(NELUSED(n)), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

              ENDIF
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Strain Energy Density...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.esed',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  704
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 704  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.04'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.esed'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Element Strain Energy Density Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart

            Nend = 0
            LBCOUNT = IBCOUNT
            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

            DO k = 1,10
              MXX = MEL_COUNT(k)
              IF (MXX .GT. 0) THEN
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(k)         
                WRITE (IO_UNIT%LEGO) CBUFFER

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (STRAIN_ENERGY_DENSITY(NELUSED(n)), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

              ENDIF
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Pressure...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.eprs',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  705
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 705  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.05'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.eprs'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Element Pressure Density Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart

            Nend = 0
            LBCOUNT = IBCOUNT
            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

            DO k = 1,10
              MXX = MEL_COUNT(k)
              IF (MXX .GT. 0) THEN
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(k)         
                WRITE (IO_UNIT%LEGO) CBUFFER

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (PRESSURE(NELUSED(n)), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

              ENDIF
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!
!! Effective stress (deviatoric stress magnitude)...
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT     =  IO_UNIT%LEGO,
     &  FILE     = 'fmaego.data.edev',
     &  STATUS   = 'UNKNOWN',
#ifdef _G95_
     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
#endif
#ifdef _CVF_NT_
     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
#endif
#ifdef LANGUAGE_FORTRAN90
     &  FORM     = 'BINARY',                        !  Pathf95
#endif
     &  POSITION = 'APPEND',
     &  ERR      =  706
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 706  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.06'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.edev'
     &    )
      ELSE
!!
!! Initialize sequential Material part counter.
!!
        MPart = 0
!!
!! Start this-time-step block.
!!
        CBUFFER = "BEGIN TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CBUFFER = "Element Effective Stress Results"
        WRITE (IO_UNIT%LEGO) CBUFFER
!!
!! Loop on Material parts.
!!
        DO M = 1,NUMMT
!!
!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!! this material. If material M is used, the function returns 1, otherwise 0. 
!!
          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN

            NNodes = MATERIAL(M)%NNodes
            NElems = MATERIAL(M)%NElems

            CBUFFER = "part"
            WRITE (IO_UNIT%LEGO) CBUFFER
            MPart = MPart + 1
            WRITE (IO_UNIT%LEGO) MPart

            Nend = 0
            LBCOUNT = IBCOUNT
            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)

            DO k = 1,10
              MXX = MEL_COUNT(k)
              IF (MXX .GT. 0) THEN
                Nbgn = Nend + 1
                Nend = Nend + MXX
                CBUFFER = ESG_ELEMENT_TYPE(k)         
                WRITE (IO_UNIT%LEGO) CBUFFER

                LBLOCKS = MXX / LBCOUNT
                LREMAIN = MXX - LBCOUNT*LBLOCKS

                Lbgn = Nbgn
                Lend = Lbgn + LREMAIN - 1
                DO i = 1,LBLOCKS+1
                  WRITE (IO_UNIT%LEGO) (EFFECTIVE_STRESS(NELUSED(n)), n = Lbgn,Lend)
                  Lbgn = Lend + 1
                  Lend = Lend + LBCOUNT
                ENDDO

              ENDIF
            ENDDO

          ENDIF
        ENDDO
!!
!! Close out this-time-step block.
!!
        CBUFFER = "END TIME STEP"
        WRITE (IO_UNIT%LEGO) CBUFFER

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')

      ENDIF
!!!!!
!!!!! Maximum Principal Stress...
!!!!!
!!!      IOERROR = .TRUE.
!!!      OPEN
!!!     &  (
!!!     &  UNIT     =  IO_UNIT%LEGO,
!!!     &  FILE     = 'fmaego.data.emxs',
!!!     &  STATUS   = 'UNKNOWN',
!!!#ifdef _G95_
!!!     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
!!!#endif
!!!#ifdef _CVF_NT_
!!!     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
!!!#endif
!!!#ifdef LANGUAGE_FORTRAN90
!!!     &  FORM     = 'BINARY',                        !  Pathf95
!!!#endif
!!!     &  POSITION = 'APPEND',
!!!     &  ERR      =  707
!!!     &  )
!!!      IOERROR = .FALSE.
!!!!!
!!!!! Fatal error exit for failed OPEN operation.
!!!!!
!!! 707  IF (IOERROR) THEN
!!!        CALL USER_MESSAGE
!!!     &    (
!!!     &    MSGL//'FATAL'//
!!!     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.07'//
!!!     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.emxs'
!!!     &    )
!!!      ELSE
!!!!!
!!!!! Initialize sequential Material part counter.
!!!!!
!!!        MPart = 0
!!!!!
!!!!! Start this-time-step block.
!!!!!
!!!        CBUFFER = "BEGIN TIME STEP"
!!!        WRITE (IO_UNIT%LEGO) CBUFFER
!!!
!!!        CBUFFER = "Element Maximum Principal Stress Results"
!!!        WRITE (IO_UNIT%LEGO) CBUFFER
!!!!!
!!!!! Loop on Material parts.
!!!!!
!!!        DO M = 1,NUMMT
!!!!!
!!!!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!!!!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!!!!! this material. If material M is used, the function returns 1, otherwise 0. 
!!!!!
!!!          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN
!!!
!!!            NNodes = MATERIAL(M)%NNodes
!!!            NElems = MATERIAL(M)%NElems
!!!
!!!            CBUFFER = "part"
!!!            WRITE (IO_UNIT%LEGO) CBUFFER
!!!            MPart = MPart + 1
!!!            WRITE (IO_UNIT%LEGO) MPart
!!!
!!!            Nend = 0
!!!            LBCOUNT = IBCOUNT
!!!            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)
!!!
!!!            DO k = 1,10
!!!              MXX = MEL_COUNT(k)
!!!              IF (MXX .GT. 0) THEN
!!!                Nbgn = Nend + 1
!!!                Nend = Nend + MXX
!!!                CBUFFER = ESG_ELEMENT_TYPE(k)         
!!!                WRITE (IO_UNIT%LEGO) CBUFFER
!!!
!!!                LBLOCKS = MXX / LBCOUNT
!!!                LREMAIN = MXX - LBCOUNT*LBLOCKS
!!!
!!!                Lbgn = Nbgn
!!!                Lend = Lbgn + LREMAIN - 1
!!!                DO i = 1,LBLOCKS+1
!!!                  WRITE (IO_UNIT%LEGO) (MAX_PRINCIPAL_STRESS(NELUSED(n)), n = Lbgn,Lend)
!!!                  Lbgn = Lend + 1
!!!                  Lend = Lend + LBCOUNT
!!!                ENDDO
!!!
!!!              ENDIF
!!!            ENDDO
!!!
!!!          ENDIF
!!!        ENDDO
!!!!!
!!!!! Close out this-time-step block.
!!!!!
!!!        CBUFFER = "END TIME STEP"
!!!        WRITE (IO_UNIT%LEGO) CBUFFER
!!!
!!!        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
!!!
!!!      ENDIF
!!!!!
!!!!! Minimum Principal Stress...
!!!!!
!!!      IOERROR = .TRUE.
!!!      OPEN
!!!     &  (
!!!     &  UNIT     =  IO_UNIT%LEGO,
!!!     &  FILE     = 'fmaego.data.emns',
!!!     &  STATUS   = 'UNKNOWN',
!!!#ifdef _G95_
!!!     &  FORM     = 'UNFORMATTED', ACCESS='STREAM',  !  G95
!!!#endif
!!!#ifdef _CVF_NT_
!!!     &  FORM     = 'BINARY', CONVERT='BIG_ENDIAN',  !  CVF
!!!#endif
!!!#ifdef LANGUAGE_FORTRAN90
!!!     &  FORM     = 'BINARY',                        !  Pathf95
!!!#endif
!!!     &  POSITION = 'APPEND',
!!!     &  ERR      =  708
!!!     &  )
!!!      IOERROR = .FALSE.
!!!!!
!!!!! Fatal error exit for failed OPEN operation.
!!!!!
!!! 708  IF (IOERROR) THEN
!!!        CALL USER_MESSAGE
!!!     &    (
!!!     &    MSGL//'FATAL'//
!!!     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.007.08'//
!!!     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.emns'
!!!     &    )
!!!      ELSE
!!!!!
!!!!! Initialize sequential Material part counter.
!!!!!
!!!        MPart = 0
!!!!!
!!!!! Start this-time-step block.
!!!!!
!!!        CBUFFER = "BEGIN TIME STEP"
!!!        WRITE (IO_UNIT%LEGO) CBUFFER
!!!
!!!        CBUFFER = "Element Minimum Principal Stress Results"
!!!        WRITE (IO_UNIT%LEGO) CBUFFER
!!!!!
!!!!! Loop on Material parts.
!!!!!
!!!        DO M = 1,NUMMT
!!!!!
!!!!! The function ELEMENTS_AND_NODES_USED(*) gen's nodal-points-used array NPTUSED, 
!!!!! elements-used array NELUSED and vtk output-index-translation array NPTNOWI for 
!!!!! this material. If material M is used, the function returns 1, otherwise 0. 
!!!!!
!!!          IF (ELEMENTS_AND_NODES_USED(M) .GT. 0) THEN
!!!
!!!            NNodes = MATERIAL(M)%NNodes
!!!            NElems = MATERIAL(M)%NElems
!!!
!!!            CBUFFER = "part"
!!!            WRITE (IO_UNIT%LEGO) CBUFFER
!!!            MPart = MPart + 1
!!!            WRITE (IO_UNIT%LEGO) MPart
!!!
!!!            Nend = 0
!!!            LBCOUNT = IBCOUNT
!!!            MEL_COUNT = (/MHX,MPX,MPY,MTX,MM3,MP3,MM4,MP4,MTR,MPH,MPG/)
!!!
!!!            DO k = 1,10
!!!              MXX = MEL_COUNT(k)
!!!              IF (MXX .GT. 0) THEN
!!!                Nbgn = Nend + 1
!!!                Nend = Nend + MXX
!!!                CBUFFER = ESG_ELEMENT_TYPE(k)         
!!!                WRITE (IO_UNIT%LEGO) CBUFFER
!!!
!!!                LBLOCKS = MXX / LBCOUNT
!!!                LREMAIN = MXX - LBCOUNT*LBLOCKS
!!!
!!!                Lbgn = Nbgn
!!!                Lend = Lbgn + LREMAIN - 1
!!!                DO i = 1,LBLOCKS+1
!!!                  WRITE (IO_UNIT%LEGO) (MIN_PRINCIPAL_STRESS(NELUSED(n)), n = Lbgn,Lend)
!!!                  Lbgn = Lend + 1
!!!                  Lend = Lend + LBCOUNT
!!!                ENDDO
!!!
!!!              ENDIF
!!!            ENDDO
!!!
!!!          ENDIF
!!!        ENDDO
!!!!!
!!!!! Close out this-time-step block.
!!!!!
!!!        CBUFFER = "END TIME STEP"
!!!        WRITE (IO_UNIT%LEGO) CBUFFER
!!!
!!!        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
!!!
!!!      ENDIF
!!
!! =========FOURTH============================================================
!! Write a *.case file to allow EnSight and ParaView to identify the mesh, 
!! node-set, side-set and wedge-set files as distinct "Parts."
!! 
!! Open an ESG case-file.
!!
      IOERROR = .TRUE.
      OPEN
     &  (
     &  UNIT    =  IO_UNIT%LEGO,
     &  FILE    = 'fmaego.data.case',
     &  STATUS  = 'UNKNOWN',
     &  FORM    = 'FORMATTED',
     &  ERR     =  800
     &  )
      IOERROR = .FALSE.
!!
!! Fatal error exit for failed OPEN operation.
!!
 800  IF (IOERROR) THEN
        CALL USER_MESSAGE
     &    (
     &    MSGL//'FATAL'//
     &    MSGL//'WRITE_TO_ESG_RESULTS_FILE.001.08'//
     &    MSGL//'Unable To Execute OPEN On: '//'fmaego.data.case'
     &    )
      ELSE

        WRITE (IO_UNIT%LEGO,'(A)') "FORMAT"
        WRITE (IO_UNIT%LEGO,'(A)') "type: ensight gold"

        !                                                   ts fs
        WRITE (IO_UNIT%LEGO,'(A)') "GEOMETRY"
        WRITE (IO_UNIT%LEGO,'(A)') "model:                   1  1  fmaego.data.geom    change_coords_only"

        !                                                   ts fs
        WRITE (IO_UNIT%LEGO,'(A)') "VARIABLE"
        WRITE (IO_UNIT%LEGO,'(A)') "vector per node:         1  1  Displacement          fmaego.data.ndis"
        WRITE (IO_UNIT%LEGO,'(A)') "vector per node:         1  1  Velocity              fmaego.data.nvel"
        WRITE (IO_UNIT%LEGO,'(A)') "vector per node:         1  1  Acceleration          fmaego.data.nacc"

        IF (NUMVX .GT. 0) THEN
          WRITE (IO_UNIT%LEGO,'(A)') "vector per node:         1  1  Max-Velocity          fmaego.data.nvmx"
          WRITE (IO_UNIT%LEGO,'(A)') "vector per node:         1  1  Min-Velocity          fmaego.data.nvmn"
        ENDIF

        WRITE (IO_UNIT%LEGO,'(A)') "scalar per element:      1  1  Material              fmaego.data.emat"
        WRITE (IO_UNIT%LEGO,'(A)') "scalar per element:      1  1  Material-State        fmaego.data.esta"
        WRITE (IO_UNIT%LEGO,'(A)') "scalar per element:      1  1  Pressure              fmaego.data.eprs"
        WRITE (IO_UNIT%LEGO,'(A)') "scalar per element:      1  1  Bulk-Strain           fmaego.data.elnv"
        WRITE (IO_UNIT%LEGO,'(A)') "scalar per element:      1  1  Strain-Energy-Density fmaego.data.esed"
        WRITE (IO_UNIT%LEGO,'(A)') "scalar per element:      1  1  Effective-Stress      fmaego.data.edev"
!!!     WRITE (IO_UNIT%LEGO,'(A)') "scalar per element:      1  1  Max-Principal-Stress  fmaego.data.emxs"
!!!     WRITE (IO_UNIT%LEGO,'(A)') "scalar per element:      1  1  Min-Principal-Stress  fmaego.data.emns"
        WRITE (IO_UNIT%LEGO,'(A)') "vector per element:      1  1  Stress*Vel            fmaego.data.estv"
        WRITE (IO_UNIT%LEGO,'(A)') "tensor symm per element: 1  1  Stress                fmaego.data.estr"

        !                                                   ts
        WRITE (IO_UNIT%LEGO,'(A)') "TIME"
        WRITE (IO_UNIT%LEGO,'(A)') "time set:                1"
        WRITE (IO_UNIT%LEGO,'(A)') "number of steps:     "//NEGO
        WRITE (IO_UNIT%LEGO,'(A,12X,2(1PE24.12)/(3(1PE24.12)))') "time values:",(TIME_VALUES(n), n= 1,TStep)

        !                                                      fs
        WRITE (IO_UNIT%LEGO,'(A)') "FILE"
        WRITE (IO_UNIT%LEGO,'(A)') "file set:                   1"
        WRITE (IO_UNIT%LEGO,'(A)') "number of steps:        "//NEGO

        CLOSE (UNIT=IO_UNIT%LEGO, STATUS='KEEP')
      ENDIF

      RETURN

      CONTAINS
!!
!!======================================================================
!! Parts-Defined-by-Material Construction Function
!!======================================================================
!!
      FUNCTION ELEMENTS_AND_NODES_USED ( M )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Purpose: Mark elements and nodes used, count number of unique elements 
!! and nodes used by the elements, and create index translation arrays 
!! for material "M".
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: M  !  Material Index
!!
!! Function return value.
      INTEGER(4) :: ELEMENTS_AND_NODES_USED
!!
!! Local Variables.
      INTEGER(ITYPE), PARAMETER :: IUNITY(MXPHN) = 1  !  (/(1,1=1,MXPHN)/)
      INTEGER(ITYPE) :: K,N,I,Imax,MmtID,NPS,KPH
!!
!! Clear marker/index-sequence and translation arrays.
!!
      NELUSED = 0
      NPTUSED = 0
      NPTNOWI = 0
!!
!! Initialize the program's element-type counters
!! for this material.
!!
      MHX = 0;  MPX = 0;  MPY = 0;  MTX = 0;  
      MM3 = 0;  MP3 = 0;  MM4 = 0;  MP4 = 0;  
      MTR = 0;  MPH = 0;  MPG = 0;
!!
!! Scan all elements for the elements using MATERIAL(M).
!! Note: The element scan order here must match that in
!! the other extraction routines: NUMHX, NUMPX, NUMPY,
!! NUMTX, NUMPH, NUMM3, NUMP3, NUMM4, NUMP4, NUMTR,
!! NUMPH, NUMPG.
!!
      K = 0
      ELEMENTS_AND_NODES_USED = 0

      DO N = 1,NUMHX
        K = K + 1
        IF (HEXAH(N)%PAR%MatID .EQ. M) THEN
          MHX = MHX + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPTUSED(HEXAH(N)%PAR%IX(1:8)) = IUNITY(1:8)
        ELSE IF (HEXAH(N)%PAR%MatID .LT. 0) THEN
          MmtID = ABS(HEXAH(N)%PAR%MatID)
          Imax = MULTIMAT(MmtID)%NUMMX
          DO I = 1,Imax
            IF (MULTIMAT(MmtID)%MaterialID(I) .EQ. M) THEN
              MHX = MHX + 1
              NELUSED(K) = 1
              ELEMENTS_AND_NODES_USED = 1
              NPTUSED(HEXAH(N)%PAR%IX(1:8)) = IUNITY(1:8)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      DO N = 1,NUMPX
        K = K + 1
        IF (PENTA(N)%PAR%MatID .EQ. M) THEN
          MPX = MPX + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPTUSED(PENTA(N)%PAR%IX(1:6)) = IUNITY(1:6)
        ENDIF
      ENDDO

      DO N = 1,NUMPY
        K = K + 1
        IF (PYRAMID(N)%PAR%MatID .EQ. M) THEN
          MPY = MPY + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPTUSED(PYRAMID(N)%PAR%IX(1:5)) = IUNITY(1:5)
        ENDIF
      ENDDO

      DO N = 1,NUMTX
        K = K + 1
        IF (TETRA(N)%PAR%MatID .EQ. M) THEN
          MTX = MTX + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1

! ParaView does not know how to handle 8-node tetrahedrons.
!          Iform = SECTION_3D( TETRA(N)%PAR%SecID )%Iform
!          IF (Iform .EQ. Eight_Nodes) Then
!            NPTUSED(TETRA(N)%PAR%IX(1:8)) = IUNITY(1:8)
!          ELSE
!            NPTUSED(TETRA(N)%PAR%IX(1:4)) = IUNITY(1:4)
!==========ENDIF

          NPTUSED(TETRA(N)%PAR%IX(1:4)) = IUNITY(1:4)
        ENDIF
      ENDDO

      DO N = 1,NUMM3
        K = K + 1
        IF (MEMBT(N)%PAR%MatID .EQ. M) THEN
          MM3 = MM3 + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPTUSED(MEMBT(N)%PAR%IX(1:3)) = IUNITY(1:3)
        ENDIF
      ENDDO

      DO N = 1,NUMP3
        K = K + 1
        IF (PLATT(N)%PAR%MatID .EQ. M) THEN
          MP3 = MP3 + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPTUSED(PLATT(N)%PAR%IX(1:3)) = IUNITY(1:3)
        ENDIF
      ENDDO

      DO N = 1,NUMM4
        K = K + 1
        IF (MEMBQ(N)%PAR%MatID .EQ. M) THEN
          MM4 = MM4 + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPTUSED(MEMBQ(N)%PAR%IX(1:4)) = IUNITY(1:4)
        ENDIF
      ENDDO

      DO N = 1,NUMP4
        K = K + 1
        IF (PLATQ(N)%PAR%MatID .EQ. M) THEN
          MP4 = MP4 + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPTUSED(PLATQ(N)%PAR%IX(1:4)) = IUNITY(1:4)
        ENDIF
      ENDDO

      DO N = 1,NUMTR
        K = K + 1
        IF (TRUSS(N)%PAR%MatID .EQ. M) THEN
          MTR = MTR + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPTUSED(TRUSS(N)%PAR%IX(1:2)) = IUNITY(1:2)
        ENDIF
      ENDDO

      DO N = 1,NUMPH
        K = K + 1
        IF (POLYH(N)%PAR%MatID .EQ. M) THEN
          MPH = MPH + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPS = POLYH(N)%PAR%NPHNP
          NPTUSED(POLYH(N)%PAR%IX(1:NPS)) = IUNITY(1:NPS)
        ELSE IF (POLYH(N)%PAR%MatID .LT. 0) THEN
          MmtID = ABS(POLYH(N)%PAR%MatID)
          Imax = MULTIMAT(MmtID)%NUMMX
          DO I = 1,Imax
            IF (MULTIMAT(MmtID)%MaterialID(I) .EQ. M) THEN
              MPH = MPH + 1
              NELUSED(K) = 1
              ELEMENTS_AND_NODES_USED = 1
              NPS = POLYH(N)%PAR%NPHNP
              NPTUSED(POLYH(N)%PAR%IX(1:NPS)) = IUNITY(1:NPS)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      DO N = 1,NUMPG
        K = K + 1
        KPH = POLYG(N)%PAR%ParID  !! Use parent polyhedron material
        IF (POLYH(KPH)%PAR%MatID .EQ. M) THEN
          MPG = MPG + 1
          NELUSED(K) = 1
          ELEMENTS_AND_NODES_USED = 1
          NPS = POLYG(N)%PAR%NPGNP
          NPTUSED(POLYG(N)%PAR%IX(1:NPS)) = IUNITY(1:NPS)
        ELSE IF (POLYH(KPH)%PAR%MatID .LT. 0) THEN
          MmtID = ABS(POLYH(KPH)%PAR%MatID)
          Imax = MULTIMAT(MmtID)%NUMMX
          DO I = 1,Imax
            IF (MULTIMAT(MmtID)%MaterialID(I) .EQ. M) THEN
              MPG = MPG + 1
              NELUSED(K) = 1
              ELEMENTS_AND_NODES_USED = 1
              NPS = POLYG(N)%PAR%NPGNP
              NPTUSED(POLYG(N)%PAR%IX(1:NPS)) = IUNITY(1:NPS)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
!!
!! Get count of elements and nodes used by this material.
!!
      MATERIAL(M)%NElems = SUM ( NELUSED(1:NUMXL) )
      MATERIAL(M)%NNodes = SUM ( NPTUSED(1:NUMNP) )

!!    write (IO_UNIT%LELO,*) ' *** Material Number:',M
!!    write (IO_UNIT%LELO,*) '              NElems:',MATERIAL(M)%NElems
!!    write (IO_UNIT%LELO,*) '              NNodes:',MATERIAL(M)%NNodes
!!    write (IO_UNIT%LELO,*) '    NELUSED(1:NUMXL):',(NELUSED(n), n=1,NUMXL)
!!    write (IO_UNIT%LELO,*) '    NPTUSED(1:NUMNP):',(NPTUSED(n), n=1,NUMNP)

!!
!! Create nodal point index map NPTNOWI for use with 
!! element connectivity n-tuples. (Maps program index 
!! to index written to esg files for this material.)
!!
!! Also, convert NPTUSED into a sequential index map.
!!
      K = 0
      DO N = 1,NUMNP
        IF (NPTUSED(N) .EQ. 1) THEN
          K = K + 1
          NPTNOWI(N) = K
          NPTUSED(K) = N
       ENDIF
      ENDDO
!!
!! Convert NELUSED into a sequential index map.
!!
      K= 0 
      DO N = 1,NUMXL
        IF (NELUSED(N) .EQ. 1) THEN
          K = K + 1
          NELUSED(K) = N
        ENDIF
      ENDDO

      RETURN
      END FUNCTION ELEMENTS_AND_NODES_USED
!!
!!======================================================================
!! Datum "Extraction" Functions. These functions access the programs
!! object-oriented datum arrays and return datum-sets for the current
!! Part, A "Part" being defined by a material model MATERIAL(*)%Type
!!======================================================================
!!
      FUNCTION ELEMENT_ESG_TYPE( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      CHARACTER(8) :: ELEMENT_ESG_TYPE
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N

      CHARACTER(8), SAVE :: ESG_VERTEX      = "point"
      CHARACTER(8), SAVE :: ESG_LINE        = "bar2"
      CHARACTER(8), SAVE :: ESG_TRIANGLE    = "tria3"
      CHARACTER(8), SAVE :: ESG_QUAD        = "quad4"
      CHARACTER(8), SAVE :: ESG_TETRAHEDRON = "tetra4"
      CHARACTER(8), SAVE :: ESG_HEXAHEDRON  = "hexa8"
      CHARACTER(8), SAVE :: ESG_WEDGE       = "penta6"
      CHARACTER(8), SAVE :: ESG_PYRAMID     = "pyramid5"
      CHARACTER(8), SAVE :: ESG_POLYHEDRON  = "nfaced"
      CHARACTER(8), SAVE :: ESG_POLYGON     = "nsided"
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      ELEMENT_ESG_TYPE = " "
      IF (NEL .LE. NHX) THEN
        N = NEL
        ELEMENT_ESG_TYPE = ESG_HEXAHEDRON  ! HEXAH(N)%PAR%EleID
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        ELEMENT_ESG_TYPE = ESG_WEDGE       ! PENTA(N)%PAR%EleID
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        ELEMENT_ESG_TYPE = ESG_PYRAMID     ! PYRAMID(N)%PAR%EleID
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        ELEMENT_ESG_TYPE = ESG_TETRAHEDRON ! TETRA(N)%PAR%EleID
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        ELEMENT_ESG_TYPE = ESG_TRIANGLE    ! MEMBT(N)%PAR%EleID
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        ELEMENT_ESG_TYPE = ESG_TRIANGLE    ! PLATT(N)%PAR%EleID
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        ELEMENT_ESG_TYPE = ESG_QUAD        ! MEMBQ(N)%PAR%EleID
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        ELEMENT_ESG_TYPE = ESG_QUAD        ! PLATQ(N)%PAR%EleID
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        ELEMENT_ESG_TYPE = ESG_LINE        ! TRUSS(N)%PAR%EleID
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        ELEMENT_ESG_TYPE = ESG_POLYHEDRON  ! POLYH(N)%PAR%EleID
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        ELEMENT_ESG_TYPE = ESG_POLYGON     ! POLYG(N)%PAR%PlgID
      ENDIF

      RETURN
      END FUNCTION ELEMENT_ESG_TYPE
!!
      FUNCTION EleID_DATA ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 12-APR-1991 18:59:44 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      INTEGER(4) :: EleID_DATA
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      EleID_DATA = 0
      IF (NEL .LE. NHX) THEN
        N = NEL
        EleID_DATA = HEXAH(N)%PAR%EleID
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        EleID_DATA = PENTA(N)%PAR%EleID
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        EleID_DATA = PYRAMID(N)%PAR%EleID
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        EleID_DATA = TETRA(N)%PAR%EleID
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        EleID_DATA = MEMBT(N)%PAR%EleID
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        EleID_DATA = PLATT(N)%PAR%EleID
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        EleID_DATA = MEMBQ(N)%PAR%EleID
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        EleID_DATA = PLATQ(N)%PAR%EleID
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        EleID_DATA = TRUSS(N)%PAR%EleID
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        EleID_DATA = POLYH(N)%PAR%EleID
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        EleID_DATA = POLYG(N)%PAR%PlgID
      ENDIF

      RETURN
      END FUNCTION EleID_DATA
!!
      FUNCTION ELEMENT_N_TUPLE( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      INTEGER(ITYPE) :: ELEMENT_N_TUPLE
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE) :: M
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      IF (NEL .LE. NHX) THEN
        N = NEL
        ELEMENT_N_TUPLE = 8
        NTUPLE(1:8) = NPTNOWI(HEXAH(N)%PAR%IX(1:8))
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        ELEMENT_N_TUPLE = 6
        NTUPLE(1:6) = NPTNOWI(PENTA(N)%PAR%IX(1:6))
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        ELEMENT_N_TUPLE = 5
        NTUPLE(1:5) = NPTNOWI(PYRAMID(N)%PAR%IX(1:5))
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        ELEMENT_N_TUPLE = 4
        NTUPLE(1:4) = NPTNOWI(TETRA(N)%PAR%IX(1:4))
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        ELEMENT_N_TUPLE = 3
        NTUPLE(1:3) = NPTNOWI(MEMBT(N)%PAR%IX(1:3))
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        ELEMENT_N_TUPLE = 3
        NTUPLE(1:3) = NPTNOWI(PLATT(N)%PAR%IX(1:3))
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        ELEMENT_N_TUPLE = 4
        NTUPLE(1:4) = NPTNOWI(MEMBQ(N)%PAR%IX(1:4))
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        ELEMENT_N_TUPLE = 4
        NTUPLE(1:4) = NPTNOWI(PLATQ(N)%PAR%IX(1:4))
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        ELEMENT_N_TUPLE = 2
        NTUPLE(1:2) = NPTNOWI(TRUSS(N)%PAR%IX(1:2))
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        M = POLYH(N)%PAR%NPHNP
        ELEMENT_N_TUPLE = M
        NTUPLE(1:M) = NPTNOWI(POLYH(N)%PAR%IX(1:M))
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        M = POLYG(N)%PAR%NPGNP
        ELEMENT_N_TUPLE = M
        NTUPLE(1:M) = NPTNOWI(POLYG(N)%PAR%IX(1:M))
      ENDIF

      RETURN
      END FUNCTION ELEMENT_N_TUPLE
!!
      FUNCTION POLYHEDRON_FACET_COUNT ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-APR-2011
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      INTEGER(ITYPE) :: POLYHEDRON_FACET_COUNT
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: M
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      POLYHEDRON_FACET_COUNT = 0

      IF (NEL .LE. NHX) THEN
      ELSE IF (NEL .LE. NPX) THEN
      ELSE IF (NEL .LE. NPY) THEN
      ELSE IF (NEL .LE. NTX) THEN
      ELSE IF (NEL .LE. NM3) THEN
      ELSE IF (NEL .LE. NP3) THEN
      ELSE IF (NEL .LE. NM4) THEN
      ELSE IF (NEL .LE. NP4) THEN
      ELSE IF (NEL .LE. NTR) THEN
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        POLYHEDRON_FACET_COUNT = POLYH(N)%PAR%NPHPG
      ELSE IF (NEL .LE. NPG) THEN
      ENDIF

      END FUNCTION POLYHEDRON_FACET_COUNT
!!
      FUNCTION POLYGON_NP_COUNTERS ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-APR-2011
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth *Polyhedral* element
!!
!! Function return value.
      INTEGER(ITYPE) :: POLYGON_NP_COUNTERS
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: KPG,k,IGk
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      POLYGON_NP_COUNTERS = 0

      IF (NEL .LE. NHX) THEN
      ELSE IF (NEL .LE. NPX) THEN
      ELSE IF (NEL .LE. NPY) THEN
      ELSE IF (NEL .LE. NTX) THEN
      ELSE IF (NEL .LE. NM3) THEN
      ELSE IF (NEL .LE. NP3) THEN
      ELSE IF (NEL .LE. NM4) THEN
      ELSE IF (NEL .LE. NP4) THEN
      ELSE IF (NEL .LE. NTR) THEN
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        KPG = POLYH(N)%PAR%NPHPG
        POLYGON_NP_COUNTERS = KPG
        DO k = 1,KPG
          IGk = POLYH(N)%PAR%IG(k)
          NTUPLE(k)= POLYG(IGk)%PAR%NPGNP
        ENDDO
      ELSE IF (NEL .LE. NPG) THEN
      ENDIF

      END FUNCTION POLYGON_NP_COUNTERS
!!
      FUNCTION MatID_DATA ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 12-APR-1991 18:59:44 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      INTEGER(ITYPE) :: MatID_DATA
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      MatID_DATA = 0
      IF (NEL .LE. NHX) THEN
        N = NEL
        MatID_DATA = MatID_HEXAH(N)
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        MatID_DATA = PENTA(N)%PAR%MatID
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        MatID_DATA = PYRAMID(N)%PAR%MatID
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        MatID_DATA = TETRA(N)%PAR%MatID
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        MatID_DATA = MEMBT(N)%PAR%MatID
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        MatID_DATA = PLATT(N)%PAR%MatID
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        MatID_DATA = MEMBQ(N)%PAR%MatID
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        MatID_DATA = PLATQ(N)%PAR%MatID
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        MatID_DATA = TRUSS(N)%PAR%MatID
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        MatID_DATA = MatID_POLYG(N)
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        MatID_DATA = 0
      ENDIF

      RETURN
      END FUNCTION MatID_DATA
!!
      FUNCTION MATERIAL_STATE ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 16-NOV-2006 18:59:44 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      REAL :: MATERIAL_STATE
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: MatID,Ibgn,Ioff,M
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      MATERIAL_STATE = ZERO
      IF (NEL .LE. NHX) THEN
        N = NEL
        MatID = MatID_HEXAH(N)
        Ibgn = HEXAH(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        MatID = PENTA(N)%PAR%MatID
        Ibgn = PENTA(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        MatID = PYRAMID(N)%PAR%MatID
        Ibgn = PYRAMID(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        MatID = TETRA(N)%PAR%MatID
        Ibgn = TETRA(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        MatID = MEMBT(N)%PAR%MatID
        Ibgn = MEMBT(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        MatID = PLATT(N)%PAR%MatID
        Ibgn = PLATT(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        MatID = MEMBQ(N)%PAR%MatID
        Ibgn = MEMBQ(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        MatID = PLATQ(N)%PAR%MatID
        Ibgn = PLATQ(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        MatID = TRUSS(N)%PAR%MatID
        Ibgn = TRUSS(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        MatID = MatID_POLYH(N)
        Ibgn = POLYH(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        MatID = MatID_POLYH(M)
        Ibgn = POLYH(N)%PAR%Isv
        Ioff = MATERIAL(MatID)%Ipt - 1
        IF (Ioff .GE. 0) MATERIAL_STATE = STATE_VARIABLES(Ibgn+Ioff)
      ENDIF

      RETURN
      END FUNCTION MATERIAL_STATE
!!
      FUNCTION STRESS_DATA ( NEL, I )
!!
!! Copyright (c) by FMA Development, LLC, 8-APR-1991 20:33:30 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL  ! Nth element
      INTEGER(ITYPE), INTENT(IN) :: I    ! Stress component requested
!!
!! Function return value.
      REAL :: STRESS_DATA
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: NDEX(4),M
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      STRESS_DATA = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        STRESS_DATA = HEXAH(N)%RES%Stress(I)
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        STRESS_DATA = PENTA(N)%RES%Stress(I)
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        STRESS_DATA = PYRAMID(N)%RES%Stress(I)
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        Iform = SECTION_3D( TETRA(N)%PAR%SecID )%Iform
        IF (Iform .EQ. Nodal_Based) Then
          NDEX = TETRA(N)%RES%IN
          STRESS_DATA = P4TH * SUM( TETNP(NDEX)%Stress(I) )
        ELSE
          STRESS_DATA = TETRA(N)%RES%Stress(I)
        ENDIF
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        IF (I .LE. 2) THEN
          STRESS_DATA = MEMBT(N)%RES%Stress(I)
        ELSE IF (I .EQ. 4) THEN
          STRESS_DATA = MEMBT(N)%RES%Stress(3)
        ENDIF
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        STRESS_DATA = PLATT(N)%RES%Arst(I)
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        IF (I .LE. 2) THEN
          STRESS_DATA = MEMBQ(N)%RES%Stress(I)
        ELSE IF (I .EQ. 4) THEN
          STRESS_DATA = MEMBQ(N)%RES%Stress(3)
        ENDIF
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        STRESS_DATA = PLATQ(N)%RES%Arst(I)
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        IF (I .EQ. 1) THEN
          STRESS_DATA = TRUSS(N)%RES%Stress
        ELSE
          STRESS_DATA = 0.0
        ENDIF
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        STRESS_DATA = POLYH(N)%RES%Stress(I)
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        STRESS_DATA = POLYH(M)%RES%Stress(I)
      ENDIF

      RETURN
      END FUNCTION STRESS_DATA
!!
      FUNCTION STRESS_FLUX ( NEL, I )
!!
!! Copyright (c) by FMA Development, LLC, 18-FEB-2007 18:56:00 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL  ! Nth element
      INTEGER(ITYPE), INTENT(IN) :: I    ! Stress*Velocity component requested
!!
!! Function return value.
      REAL :: STRESS_FLUX
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: NDEX(MXPHN),K,M
      INTEGER(ITYPE), SAVE :: ISTR(3,3)

       DATA ISTR(1:3,1) /1,4,5/
       DATA ISTR(1:3,2) /4,2,6/
       DATA ISTR(1:3,3) /5,6,3/
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      STRESS_FLUX = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        NDEX(1:8) = HEXAH(N)%PAR%IX(1:8)
        Vx = P8TH * SUM( MOTION(NDEX(1:8))%Vx )
        Vy = P8TH * SUM( MOTION(NDEX(1:8))%Vy )
        Vz = P8TH * SUM( MOTION(NDEX(1:8))%Vz )
        STRESS_FLUX = DOT_PRODUCT( (/Vx,Vy,Vz/), HEXAH(N)%RES%Stress(ISTR(1:3,I)) )
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        NDEX(1:6) = PENTA(N)%PAR%IX(1:6)
        Vx = P6TH * SUM( MOTION(NDEX(1:6))%Vx )
        Vy = P6TH * SUM( MOTION(NDEX(1:6))%Vy )
        Vz = P6TH * SUM( MOTION(NDEX(1:6))%Vz )
        STRESS_FLUX = DOT_PRODUCT( (/Vx,Vy,Vz/), PENTA(N)%RES%Stress(ISTR(1:3,I)) )
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        NDEX(1:5) = PYRAMID(N)%PAR%IX(1:5)
        Vx = P5TH * SUM( MOTION(NDEX(1:5))%Vx )
        Vy = P5TH * SUM( MOTION(NDEX(1:5))%Vy )
        Vz = P5TH * SUM( MOTION(NDEX(1:5))%Vz )
        STRESS_FLUX = DOT_PRODUCT( (/Vx,Vy,Vz/), PYRAMID(N)%RES%Stress(ISTR(1:3,I)) )
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        Iform = SECTION_3D( TETRA(N)%PAR%SecID )%Iform
        IF (Iform .EQ. Nodal_Based) Then
          NDEX(1:4) = TETRA(N)%PAR%IX(1:4)
          Vx = P4TH * SUM( MOTION(NDEX(1:4))%Vx )
          Vy = P4TH * SUM( MOTION(NDEX(1:4))%Vy )
          Vz = P4TH * SUM( MOTION(NDEX(1:4))%Vz )
          NDEX(1:4) = TETRA(N)%RES%IN(1:4)
          S1 = P4TH * SUM( TETNP(NDEX(1:4))%Stress(ISTR(1,I)) )
          S2 = P4TH * SUM( TETNP(NDEX(1:4))%Stress(ISTR(2,I)) )
          S3 = P4TH * SUM( TETNP(NDEX(1:4))%Stress(ISTR(3,I)) )
          STRESS_FLUX = DOT_PRODUCT( (/Vx,Vy,Vz/),(/S1,S2,S3/) )
        ELSE
          NDEX(1:4) = TETRA(N)%PAR%IX(1:4)
          Vx = P4TH * SUM( MOTION(NDEX(1:4))%Vx )
          Vy = P4TH * SUM( MOTION(NDEX(1:4))%Vy )
          Vz = P4TH * SUM( MOTION(NDEX(1:4))%Vz )
          STRESS_FLUX = DOT_PRODUCT( (/Vx,Vy,Vz/), TETRA(N)%RES%Stress(ISTR(1:3,I)) )
        ENDIF
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        STRESS_FLUX = ZERO
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        STRESS_FLUX = ZERO
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        STRESS_FLUX = ZERO
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        STRESS_FLUX = ZERO
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        M = POLYH(N)%PAR%NPHNP
        NDEX(1:M) = POLYH(N)%PAR%IX(1:M)
        QMTH = PONE / REAL(M,KIND(0.0D+0))
        Vx = QMTH * SUM( MOTION(NDEX(1:M))%Vx )
        Vy = QMTH * SUM( MOTION(NDEX(1:M))%Vy )
        Vz = QMTH * SUM( MOTION(NDEX(1:M))%Vz )
        STRESS_FLUX = DOT_PRODUCT( (/Vx,Vy,Vz/), POLYH(N)%RES%Stress(ISTR(1:3,I)) )
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        K = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        M = POLYH(K)%PAR%NPHNP
        NDEX(1:M) = POLYH(K)%PAR%IX(1:M)
        QMTH = PONE / REAL(M,KIND(0.0D+0))
        Vx = QMTH * SUM( MOTION(NDEX(1:M))%Vx )
        Vy = QMTH * SUM( MOTION(NDEX(1:M))%Vy )
        Vz = QMTH * SUM( MOTION(NDEX(1:M))%Vz )
        STRESS_FLUX = DOT_PRODUCT( (/Vx,Vy,Vz/), POLYH(K)%RES%Stress(ISTR(1:3,I)) )
      ENDIF

      RETURN
      END FUNCTION STRESS_FLUX
!!
      FUNCTION MAX_PRINCIPAL_STRESS ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 8-APR-1991 20:33:30 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL  ! Nth element
!!
!! Function return value.
      REAL :: MAX_PRINCIPAL_STRESS
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: NDEX(4),i,M
      REAL(RTYPE)          :: SIGMA(6),EV(3),EVEC(3,3)
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      MAX_PRINCIPAL_STRESS = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        SIGMA(1:6) = HEXAH(N)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        SIGMA(1:6) = PENTA(N)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        SIGMA(1:6) = PYRAMID(N)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        Iform = SECTION_3D( TETRA(N)%PAR%SecID )%Iform
        IF (Iform .EQ. Nodal_Based) Then
          NDEX = TETRA(N)%RES%IN
          DO i = 1,6
            SIGMA(i) = P4TH * SUM( TETNP(NDEX)%Stress(i) )
          ENDDO
        ELSE
          SIGMA(1:6) = TETRA(N)%RES%Stress(1:6)
        ENDIF
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        SIGMA(1:3) = MEMBT(N)%RES%Stress(1:3)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:2),ZERO,SIGMA(3),ZERO,ZERO/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        SIGMA(1:6) = PLATT(N)%RES%Arst
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        SIGMA(1:3) = MEMBQ(N)%RES%Stress(1:3)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:2),ZERO,SIGMA(3),ZERO,ZERO/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        Ist = PLATQ(N)%PAR%Ist
        SIGMA(1:6) = STRESS(I,Ist)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        MAX_PRINCIPAL_STRESS = MAX( ZERO,TRUSS(N)%RES%Stress )
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        SIGMA(1:6) = POLYH(N)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        SIGMA(1:6) = POLYH(M)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MAX_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ENDIF
      RETURN
      END FUNCTION MAX_PRINCIPAL_STRESS
!!
      FUNCTION MIN_PRINCIPAL_STRESS ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 8-APR-1991 20:33:30 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL  ! Nth element
!!
!! Function return value.
      REAL :: MIN_PRINCIPAL_STRESS
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: NDEX(4),i,M
      REAL(RTYPE)          :: SIGMA(6),EV(3),EVEC(3,3)
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      MIN_PRINCIPAL_STRESS = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        SIGMA(1:6) = HEXAH(N)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MIN( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        SIGMA(1:6) = PENTA(N)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MIN( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        SIGMA(1:6) = PYRAMID(N)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MIN( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        Iform = SECTION_3D( TETRA(N)%PAR%SecID )%Iform
        IF (Iform .EQ. Nodal_Based) Then
          NDEX = TETRA(N)%RES%IN
          DO i = 1,6
            SIGMA(i) = P4TH * SUM( TETNP(NDEX)%Stress(i) )
          ENDDO
        ELSE
          SIGMA(1:6) = TETRA(N)%RES%Stress(1:6)
        ENDIF
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MIN( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        SIGMA(1:3) = MEMBT(N)%RES%Stress(1:3)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:2),ZERO,SIGMA(3),ZERO,ZERO/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MIN( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        SIGMA(1:6) = PLATT(N)%RES%Arst
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MIN( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        SIGMA(1:3) = MEMBQ(N)%RES%Stress(1:3)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:2),ZERO,SIGMA(3),ZERO,ZERO/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MIN( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        Ist = PLATQ(N)%PAR%Ist
        SIGMA(1:6) = STRESS(I,Ist)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MIN( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        MIN_PRINCIPAL_STRESS = MIN( ZERO,TRUSS(N)%RES%Stress )
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        SIGMA(1:6) = POLYH(N)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        SIGMA(1:6) = POLYH(M)%RES%Stress(1:6)
        CALL GET_EIGENVALVES_EIGENVECTORS ( (/SIGMA(1:4),SIGMA(6),SIGMA(5)/),EV,EVEC )
        MIN_PRINCIPAL_STRESS = MAX( EV(1),EV(2),EV(3) )
      ENDIF

      RETURN
      END FUNCTION MIN_PRINCIPAL_STRESS
!!
      FUNCTION BULK_STRAIN ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 12-APR-1991 18:59:44 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL  ! Nth element
!!
!! Function return value.
      REAL :: BULK_STRAIN
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: NDEX(4),M
      REAL(RTYPE)          :: AVERAGE_VOLUME
!!
!! Function reference.
      REAL(RTYPE), EXTERNAL :: LOG_POS
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      BULK_STRAIN = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        BULK_STRAIN = LOG_POS(HEXAH(N)%RES%Volume/HEXAH(N)%PAR%Volume)
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        BULK_STRAIN = LOG_POS(PENTA(N)%RES%Volume/PENTA(N)%PAR%Volume)
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        BULK_STRAIN = LOG_POS(PYRAMID(N)%RES%Volume/PYRAMID(N)%PAR%Volume)
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        Iform = SECTION_3D( TETRA(N)%PAR%SecID )%Iform
        IF (Iform .EQ. Nodal_Based) Then
          NDEX = TETRA(N)%RES%IN
          AVERAGE_VOLUME = P4TH * SUM( TETNP(NDEX)%Volume )
          BULK_STRAIN = LOG_POS(AVERAGE_VOLUME/TETRA(N)%PAR%Volume)
        ELSE
          BULK_STRAIN = LOG_POS(TETRA(N)%RES%Volume/TETRA(N)%PAR%Volume)
        ENDIF
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        BULK_STRAIN = LOG_POS(MEMBT(N)%RES%Area/MEMBT(N)%PAR%Area)
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        BULK_STRAIN = LOG_POS(PLATT(N)%RES%Area/PLATT(N)%PAR%Area)
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        BULK_STRAIN = LOG_POS(MEMBQ(N)%RES%Area/MEMBQ(N)%PAR%Area)
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        BULK_STRAIN = LOG_POS(PLATQ(N)%RES%Area/PLATQ(N)%PAR%Area)
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        BULK_STRAIN = LOG_POS(TRUSS(N)%RES%Length/TRUSS(N)%PAR%Length)
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        BULK_STRAIN = LOG_POS(POLYH(N)%RES%Volume/POLYH(N)%PAR%Volume)
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        BULK_STRAIN = LOG_POS(POLYH(M)%RES%Volume/POLYH(M)%PAR%Volume)
      ENDIF

      RETURN
      END FUNCTION BULK_STRAIN
!!
      FUNCTION STRAIN_ENERGY_DENSITY ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 12-APR-1991 18:59:44 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      REAL :: STRAIN_ENERGY_DENSITY
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: NDEX(4),M
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      STRAIN_ENERGY_DENSITY = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        STRAIN_ENERGY_DENSITY = HEXAH(N)%RES%Str_Eng
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        STRAIN_ENERGY_DENSITY = PENTA(N)%RES%Str_Eng
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        STRAIN_ENERGY_DENSITY = PYRAMID(N)%RES%Str_Eng
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        Iform = SECTION_3D( TETRA(N)%PAR%SecID )%Iform
        IF (Iform .EQ. Nodal_Based) Then
          NDEX = TETRA(N)%RES%IN
          STRAIN_ENERGY_DENSITY = P4TH * SUM( TETNP(NDEX)%Str_Eng )
        ELSE
          STRAIN_ENERGY_DENSITY = TETRA(N)%RES%Str_Eng
        ENDIF
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        STRAIN_ENERGY_DENSITY = MEMBT(N)%RES%Str_Eng
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        Iform = SECTION_2D( PLATT(N)%PAR%SecID )%Iform
        IF (Iform .EQ. Nodal_Based_P3EL) Then
          NDEX(1:3) = PLATT(N)%RES%IN
          STRAIN_ENERGY_DENSITY = P3RD * SUM( PLTNP(NDEX(1:3))%Str_Eng )
        ELSE
          STRAIN_ENERGY_DENSITY = PLATT(N)%RES%Str_Eng
        ENDIF
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        STRAIN_ENERGY_DENSITY = MEMBQ(N)%RES%Str_Eng
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        STRAIN_ENERGY_DENSITY = PLATQ(N)%RES%Str_Eng
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        STRAIN_ENERGY_DENSITY = TRUSS(N)%RES%Str_Eng
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        STRAIN_ENERGY_DENSITY = POLYH(N)%RES%Str_Eng
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        STRAIN_ENERGY_DENSITY = POLYH(M)%RES%Str_Eng
      ENDIF
!!
      RETURN
      END FUNCTION STRAIN_ENERGY_DENSITY
!!
      FUNCTION PRESSURE ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 8-APR-1991 20:33:30 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      REAL :: PRESSURE
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: NDEX(4),M
      REAL(RTYPE)          :: TEMP(4)
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      PRESSURE = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        PRESSURE = P3RD * SUM( HEXAH(N)%RES%Stress(1:3) )
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        PRESSURE = P3RD * SUM( PENTA(N)%RES%Stress(1:3) )
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        PRESSURE = P3RD * SUM( PYRAMID(N)%RES%Stress(1:3) )
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        Iform = SECTION_3D( TETRA(N)%PAR%SecID )%Iform
        IF (Iform .EQ. Nodal_Based) Then
          NDEX = TETRA(N)%RES%IN
          TEMP = (/ZERO,ZERO,ZERO,ZERO/)
          TEMP = TEMP + TETNP(NDEX)%Stress(1)
          TEMP = TEMP + TETNP(NDEX)%Stress(2)
          TEMP = TEMP + TETNP(NDEX)%Stress(3)
          PRESSURE = P12TH * SUM( TEMP ) 
        ELSE
          PRESSURE = P3RD * SUM( TETRA(N)%RES%Stress(1:3) )
        ENDIF
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        PRESSURE = P3RD * SUM( MEMBT(N)%RES%Stress(1:2) )
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        PRESSURE = P3RD * SUM( PLATT(N)%RES%Arst(1:2) )
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        PRESSURE = P3RD * SUM( MEMBQ(N)%RES%Stress(1:2) )
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        Ist = PLATQ(N)%PAR%Ist
        PRESSURE = P3RD * SUM( STRESS(1:2,Ist) )
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        PRESSURE = P3RD * TRUSS(N)%RES%Stress
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        PRESSURE = P3RD * SUM( POLYH(N)%RES%Stress(1:3) )
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        PRESSURE = P3RD * SUM( POLYH(M)%RES%Stress(1:3) )
      ENDIF

      RETURN
      END FUNCTION PRESSURE
!!
      FUNCTION EFFECTIVE_STRESS ( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 8-APR-1991 20:33:30 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      REAL :: EFFECTIVE_STRESS
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: NDEX(4),M
      REAL(RTYPE)          :: TEMP(4)
      REAL(RTYPE)          :: PRESSURE,S1,S2,S3,S4,S5,S6
!!
!! Function reference.
      REAL(RTYPE) :: EFF_STR
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      EFFECTIVE_STRESS = 0.0
      IF (NEL .LE. NHX) THEN
        N = NEL
        PRESSURE = P3RD * SUM( HEXAH(N)%RES%Stress(1:3) )
        S1 = HEXAH(N)%RES%Stress(1) - PRESSURE
        S2 = HEXAH(N)%RES%Stress(2) - PRESSURE
        S3 = HEXAH(N)%RES%Stress(3) - PRESSURE
        S4 = HEXAH(N)%RES%Stress(4) * HEXAH(N)%RES%Stress(4)
        S5 = HEXAH(N)%RES%Stress(5) * HEXAH(N)%RES%Stress(5) 
        S6 = HEXAH(N)%RES%Stress(6) * HEXAH(N)%RES%Stress(6)
        S4 = S4 + S5 + S6
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        PRESSURE = P3RD * SUM( PENTA(N)%RES%Stress(1:3) )
        S1 = PENTA(N)%RES%Stress(1) - PRESSURE
        S2 = PENTA(N)%RES%Stress(2) - PRESSURE
        S3 = PENTA(N)%RES%Stress(3) - PRESSURE
        S4 = PENTA(N)%RES%Stress(4) * PENTA(N)%RES%Stress(4)
        S5 = PENTA(N)%RES%Stress(5) * PENTA(N)%RES%Stress(5) 
        S6 = PENTA(N)%RES%Stress(6) * PENTA(N)%RES%Stress(6)
        S4 = S4 + S5 + S6
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        PRESSURE = P3RD * SUM( PYRAMID(N)%RES%Stress(1:3) )
        S1 = PYRAMID(N)%RES%Stress(1) - PRESSURE
        S2 = PYRAMID(N)%RES%Stress(2) - PRESSURE
        S3 = PYRAMID(N)%RES%Stress(3) - PRESSURE
        S4 = PYRAMID(N)%RES%Stress(4) * PYRAMID(N)%RES%Stress(4)
        S5 = PYRAMID(N)%RES%Stress(5) * PYRAMID(N)%RES%Stress(5) 
        S6 = PYRAMID(N)%RES%Stress(6) * PYRAMID(N)%RES%Stress(6)
        S4 = S4 + S5 + S6
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        Iform = SECTION_3D( TETRA(N)%PAR%SecID )%Iform
        IF (Iform .EQ. Nodal_Based) THEN
          NDEX = TETRA(N)%RES%IN
          TEMP = (/ZERO,ZERO,ZERO,ZERO/)
          TEMP = TEMP + TETNP(NDEX)%Stress(1)
          TEMP = TEMP + TETNP(NDEX)%Stress(2)
          TEMP = TEMP + TETNP(NDEX)%Stress(3)
          PRESSURE = P12TH * SUM( TEMP ) 
        ELSE
          PRESSURE = P3RD * SUM( TETRA(N)%RES%Stress(1:3) )
        ENDIF
        IF (Iform .EQ. Nodal_Based) THEN
          S1 =  P4TH * SUM( TETNP(NDEX)%Stress(1) ) - PRESSURE
          S2 =  P4TH * SUM( TETNP(NDEX)%Stress(2) ) - PRESSURE
          S3 =  P4TH * SUM( TETNP(NDEX)%Stress(3) ) - PRESSURE
          S4 = (P4TH * SUM( TETNP(NDEX)%Stress(4) ) )**2
          S5 = (P4TH * SUM( TETNP(NDEX)%Stress(5) ) )**2
          S6 = (P4TH * SUM( TETNP(NDEX)%Stress(6) ) )**2
          S4 = (S4 + S5 + S6)
        ELSE
          S1 = TETRA(N)%RES%Stress(1) - PRESSURE
          S2 = TETRA(N)%RES%Stress(2) - PRESSURE
          S3 = TETRA(N)%RES%Stress(3) - PRESSURE
          S4 = TETRA(N)%RES%Stress(4) * TETRA(N)%RES%Stress(4)
          S5 = TETRA(N)%RES%Stress(5) * TETRA(N)%RES%Stress(5) 
          S6 = TETRA(N)%RES%Stress(6) * TETRA(N)%RES%Stress(6)
          S4 = S4 + S5 + S6
        ENDIF
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        PRESSURE = P3RD * SUM ( MEMBT(N)%RES%Stress(1:2) )
        S1 = MEMBT(N)%RES%Stress(1) - PRESSURE
        S2 = MEMBT(N)%RES%Stress(2) - PRESSURE
        S3 =                        - PRESSURE
        S4 = MEMBT(N)%RES%Stress(3) * MEMBT(N)%RES%Stress(3)
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        PRESSURE = P3RD * SUM( PLATT(N)%RES%Arst(1:2) )
        S1 = PLATT(N)%RES%Arst(1) - PRESSURE
        S2 = PLATT(N)%RES%Arst(2) - PRESSURE
        S3 =                      - PRESSURE
        S4 = PLATT(N)%RES%Arst(4)**2
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        PRESSURE = P3RD * SUM( MEMBQ(N)%RES%Stress(1:2) )
        S1 = MEMBQ(N)%RES%Stress(1) - PRESSURE
        S2 = MEMBQ(N)%RES%Stress(2) - PRESSURE
        S3 =                        - PRESSURE
        S4 = MEMBQ(N)%RES%Stress(3) * MEMBQ(N)%RES%Stress(3)
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4 
        Ist = PLATQ(N)%PAR%Ist
        PRESSURE = P3RD * SUM( STRESS(1:2,Ist) )
        S1 = STRESS(1,Ist) - PRESSURE
        S2 = STRESS(2,Ist) - PRESSURE
        S3 =               - PRESSURE
        S4 = STRESS(4,Ist) * STRESS(4,Ist) 
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        PRESSURE = P3RD * TRUSS(N)%RES%Stress
        S1 = TRUSS(N)%RES%Stress - PRESSURE
        S2 =                     - PRESSURE
        S3 =                     - PRESSURE
        S4 = 0.0
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        PRESSURE = P3RD * SUM( POLYH(N)%RES%Stress(1:3) )
        S1 = POLYH(N)%RES%Stress(1) - PRESSURE
        S2 = POLYH(N)%RES%Stress(2) - PRESSURE
        S3 = POLYH(N)%RES%Stress(3) - PRESSURE
        S4 = POLYH(N)%RES%Stress(4) * POLYH(N)%RES%Stress(4)
        S5 = POLYH(N)%RES%Stress(5) * POLYH(N)%RES%Stress(5) 
        S6 = POLYH(N)%RES%Stress(6) * POLYH(N)%RES%Stress(6)
        S4 = S4 + S5 + S6
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        PRESSURE = P3RD * SUM( POLYH(M)%RES%Stress(1:3) )
        S1 = POLYH(M)%RES%Stress(1) - PRESSURE
        S2 = POLYH(M)%RES%Stress(2) - PRESSURE
        S3 = POLYH(M)%RES%Stress(3) - PRESSURE
        S4 = POLYH(M)%RES%Stress(4) * POLYH(M)%RES%Stress(4)
        S5 = POLYH(M)%RES%Stress(5) * POLYH(M)%RES%Stress(5) 
        S6 = POLYH(M)%RES%Stress(6) * POLYH(M)%RES%Stress(6)
        S4 = S4 + S5 + S6
        EFFECTIVE_STRESS = EFF_STR (S1,S2,S3,S4)
      ENDIF

      RETURN
      END FUNCTION EFFECTIVE_STRESS
!!
      FUNCTION ELEMENT_VOLUME( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 12-APR-1991 18:59:44 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      REAL(4) :: ELEMENT_VOLUME
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: M
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      ELEMENT_VOLUME= ZERO
      IF (NEL .LE. NHX) THEN
        N = NEL
        ELEMENT_VOLUME= HEXAH(N)%PAR%Volume
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        ELEMENT_VOLUME= PENTA(N)%PAR%Volume
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        ELEMENT_VOLUME= PYRAMID(N)%PAR%Volume
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        ELEMENT_VOLUME= TETRA(N)%PAR%Volume
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        ELEMENT_VOLUME= MEMBT(N)%PAR%Area
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        ELEMENT_VOLUME= PLATT(N)%PAR%Area
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        ELEMENT_VOLUME= MEMBQ(N)%PAR%Area
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        ELEMENT_VOLUME= PLATQ(N)%PAR%Area
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        ELEMENT_VOLUME= TRUSS(N)%PAR%Length
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        ELEMENT_VOLUME= POLYH(N)%PAR%Volume
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        ELEMENT_VOLUME= POLYH(M)%PAR%Volume
      ENDIF

      RETURN
      END FUNCTION ELEMENT_VOLUME
!!
      FUNCTION ELEMENT_CRITICAL_DT( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 12-APR-1991 18:59:44 
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Nth element
!!
!! Function return value.
      REAL(4) :: ELEMENT_CRITICAL_DT
!!
!! Local variables.
      LOGICAL,        SAVE :: FIRST = .TRUE.
      INTEGER(ITYPE), SAVE :: NHX,NPX,NPY,NTX,NM3,NP3,NM4,NP4,NTR,NPH,NPG,N
      INTEGER(ITYPE)       :: M
!!
      IF (FIRST) THEN
        NHX = NUMHX
        NPX = NHX + NUMPX
        NPY = NPX + NUMPY
        NTX = NPY + NUMTX
        NM3 = NTX + NUMM3
        NP3 = NM3 + NUMP3
        NM4 = NP3 + NUMM4
        NP4 = NM4 + NUMP4
        NTR = NP4 + NUMTR
        NPH = NTR + NUMPH
        NPG = NPH + NUMPG
        FIRST = .FALSE.
      ENDIF
!!
      ELEMENT_CRITICAL_DT= ZERO
      IF (NEL .LE. NHX) THEN
        N = NEL
        ELEMENT_CRITICAL_DT= HEXAH(N)%RES%DTelt
      ELSE IF (NEL .LE. NPX) THEN
        N = NEL - NHX
        ELEMENT_CRITICAL_DT= PENTA(N)%RES%DTelt
      ELSE IF (NEL .LE. NPY) THEN
        N = NEL - NPX
        ELEMENT_CRITICAL_DT= PYRAMID(N)%RES%DTelt
      ELSE IF (NEL .LE. NTX) THEN
        N = NEL - NPY
        ELEMENT_CRITICAL_DT= TETRA(N)%RES%DTelt
      ELSE IF (NEL .LE. NM3) THEN
        N = NEL - NTX
        ELEMENT_CRITICAL_DT= MEMBT(N)%RES%DTelt
      ELSE IF (NEL .LE. NP3) THEN
        N = NEL - NM3
        ELEMENT_CRITICAL_DT= PLATT(N)%RES%DTelt
      ELSE IF (NEL .LE. NM4) THEN
        N = NEL - NP3
        ELEMENT_CRITICAL_DT= MEMBQ(N)%RES%DTelt
      ELSE IF (NEL .LE. NP4) THEN
        N = NEL - NM4
        ELEMENT_CRITICAL_DT= PLATQ(N)%RES%DTelt
      ELSE IF (NEL .LE. NTR) THEN
        N = NEL - NP4
        ELEMENT_CRITICAL_DT= TRUSS(N)%RES%DTelt
      ELSE IF (NEL .LE. NPH) THEN
        N = NEL - NTR
        ELEMENT_CRITICAL_DT= POLYH(N)%RES%DTelt
      ELSE IF (NEL .LE. NPG) THEN
        N = NEL - NPH
        M = POLYG(N)%PAR%ParID  !! Use parent polyhedron
        ELEMENT_CRITICAL_DT= POLYH(M)%RES%DTelt
      ENDIF

      RETURN
      END FUNCTION ELEMENT_CRITICAL_DT
!!
      FUNCTION SEGMENTS_AND_NODES_USED ( M )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Purpose: Mark segments and nodes used, count number of unique segments
!! and nodes used by the segment set, and create index translation arrays 
!! for segment set "M".
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: M   !  Segment set ID
!!
!! Function return value.
      INTEGER(ITYPE) :: SEGMENTS_AND_NODES_USED
!!
!! Local Variables.
      INTEGER(ITYPE), SAVE :: IX(MXPGN),N,K,KNP
!!
!! Clear marker/index-sequence and translation arrays.
!!
      NELUSED = 0
      NPTUSED = 0
      NPTNOWI = 0
!!
!! Mark nodes used by this segment set.
!!
      SEGMENTS_AND_NODES_USED = 0

      N = 0
      DO WHILE (NEXT_SEG_ID(M,N))
        NELUSED(N) = 1
        KNP = SEGMENT(N)%PAR%KNP
        IX(1:KNP) = SEGMENT(N)%PAR%IX(1:KNP)
        IF (IX(4) .LE. 0) THEN
          NPTUSED(IX(1:3)) = (/1,1,1/)            
        ELSE
          NPTUSED(IX(1:KNP)) = (/(1,k=1,KNP)/)
        ENDIF
        SEGMENTS_AND_NODES_USED = 1
      ENDDO
!!
!! Create nodal point index map NPTNOWI for use with 
!! element connectivity n-tuples. (Maps program index 
!! to index written to esg files for this segment set.)
!!
!! Also, convert NPTUSED into a sequential index map.
!!
      K = 0
      DO N = 1,NUMNP
        IF (NPTUSED(N) .EQ. 1) THEN
          K = K + 1
          NPTNOWI(N) = K
          NPTUSED(K) = N
        ENDIF
      ENDDO
      NNodes = K
!!
!! Convert NELUSED into a sequential index map.
!!
      K= 0 
      DO N = 1,NUMSG
        IF (NELUSED(N) .EQ. 1) THEN
          K = K + 1
          NELUSED(K) = N
        ENDIF
      ENDDO
      NElems = K

      RETURN
      END FUNCTION SEGMENTS_AND_NODES_USED
!!
      FUNCTION SEGMENT_N_TUPLE( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Segment ID
!!
!! Function return value.
      INTEGER(ITYPE) :: SEGMENT_N_TUPLE
!!
!! Local variables.
      INTEGER(ITYPE), SAVE :: IX(MXPGN),KNP
!!
      KNP = SEGMENT(NEL)%PAR%KNP
      IX(1:KNP) = SEGMENT(NEL)%PAR%IX(1:KNP)
      IF (KNP.EQ.4 .AND. IX(4).LT.0) THEN
        !!
        !! SEGMENTS are defined as polygons, that is, 
        !! as "nsided" elements. Note that for this 
        !! case, KNP=4 has already been written to the
        !! EnSight-formatted file. Thus, we cannot 
        !! change it into a 3-sided polygon.
        !!
        IX(4) = IX(3)
        SEGMENT_N_TUPLE = 4
        NTUPLE(1:4) = NPTNOWI(IX(1:4))
      ELSE
        SEGMENT_N_TUPLE = KNP
        NTUPLE(1:KNP) = NPTNOWI(IX(1:KNP))
      ENDIF

      RETURN
      END FUNCTION SEGMENT_N_TUPLE
!!
      FUNCTION WEDGES_AND_NODES_USED ( M )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Purpose: Mark wedges and nodes used, count number of unique wedges
!! and nodes used by the wedge set, and create index translation arrays 
!! for wedge set "M".
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: M   ! Solid-to-solid interface ID
!!
!! Function return value.
      INTEGER(ITYPE) :: WEDGES_AND_NODES_USED
!!
!! Local Variables.
      INTEGER(ITYPE), SAVE :: IX(4),N,K
      LOGICAL              :: FOUND
!!
!! Clear marker/index-sequence and translation arrays.
!!
      NELUSED = 0
      NPTUSED = 0
      NPTNOWI = 0
!!
!! Mark wedges used by this solid-to-solid interface.
!!
      WEDGES_AND_NODES_USED = 0

      DO N = 1,NUMWX
        IF (WEDGE(N)%MPCID .EQ. M) THEN
          NELUSED(N) = 1
          WEDGES_AND_NODES_USED = 1
        ENDIF
      ENDDO
!!
!! Note: The wedges don't use the nodal points of the
!! mesh, but contain information for building their
!! own vertex coordinates, see WEDGE_NP_COORDINATE.
!!
!! Create nodal point index map NPTNOWI for use with 
!! element connectivity n-tuples. (Maps program index 
!! to index written to esg files for this wedge set.)
!!
!! Also, convert NPTUSED into a sequential index map.
!!
!!      K = 0
!!      DO N = 1,NUMNP
!!        IF (NPTUSED(N) .EQ. 1) THEN
!!          K = K + 1
!!          NPTNOWI(N) = K
!!          NPTUSED(K) = N
!!        ENDIF
!!      ENDDO
!!      NNodes = K
!!
!! Convert NELUSED into a sequential index map.
!!
      K= 0 
      DO N = 1,NUMWX
        IF (NELUSED(N) .EQ. 1) THEN
          K = K + 1
          NELUSED(K) = N
        ENDIF
      ENDDO

      NElems = K
      NNodes = 6 * NElems

      RETURN
      END FUNCTION WEDGES_AND_NODES_USED
!!
      FUNCTION WEDGE_N_TUPLE( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Wedge ID
!!
!! Function return value.
      INTEGER(ITYPE) :: WEDGE_N_TUPLE
!!
!! This routine supplies element connectivity n-tuples.
!!
      IOS = 6*NEL - 6

      NTUPLE(1:6) = (/1+IOS,2+IOS,3+IOS,4+IOS,5+IOS,6+IOS/)

      WEDGE_N_TUPLE = 6

      RETURN
      END FUNCTION WEDGE_N_TUPLE
!!
      FUNCTION WEDGE_NP_COORDINATE( I, J )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: I     ! Nodal counter
      INTEGER(ITYPE), INTENT(IN) :: J     ! 1/2/3=x/y/z
!!
!! Function return value.
      REAL(RTYPE) :: WEDGE_NP_COORDINATE
!!
!! Local variables.
      REAL(RTYPE)    :: PXYZ,QXYZ,W
      INTEGER(ITYPE) :: NEL,NPT,N
!!
!! Caution: This routine contains a hard-coded dependency.
!! This routine will generate nodal point coordinates for
!! the sequence of wedges contained in NELUSED(*). Each
!! wedge has its own six nodal points.
!!
      NEL = ((I-1)/6) + 1
      NPT = I - 6*((I-1)/6)

      MEL = NELUSED(NEL)

      PXYZ = ZERO

      DO k = 1,WEDGE(MEL)%K(NPT)
        N =    WEDGE(MEL)%N(k,NPT)
        W =    WEDGE(MEL)%W(k,NPT)
        SELECT CASE (J)
          CASE(1);  QXYZ = W * MOTION(N)%Px
          CASE(2);  QXYZ = W * MOTION(N)%Py
          CASE(3);  QXYZ = W * MOTION(N)%Pz
        END SELECT
        PXYZ = PXYZ + QXYZ
      ENDDO

      WEDGE_NP_COORDINATE = PXYZ

      RETURN
      END FUNCTION WEDGE_NP_COORDINATE
!!
      FUNCTION WEDGE_NP_VELOCITY( I, J )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: I     ! Nodal counter
      INTEGER(ITYPE), INTENT(IN) :: J     ! 1/2/3=x/y/z
!!
!! Function return value.
      REAL(RTYPE) :: WEDGE_NP_VELOCITY
!!
!! Local variables.
      REAL(RTYPE)    :: VXYZ,QXYZ,W
      INTEGER(ITYPE) :: NEL,NPT,N
!!
!! Caution: This routine contains a hard-coded dependency.
!! This routine will generate nodal point velocities for
!! the sequence of wedges contained in NELUSED(*). Each
!! wedge has its own six nodal points.
!!
      NEL = ((I-1)/6) + 1
      NPT = I - 6*((I-1)/6)

      MEL = NELUSED(NEL)

      VXYZ = ZERO

      DO k = 1,WEDGE(MEL)%K(NPT)
        N =    WEDGE(MEL)%N(k,NPT)
        W =    WEDGE(MEL)%W(k,NPT)
        SELECT CASE (J)
          CASE(1);  QXYZ = W * MOTION(N)%Vx
          CASE(2);  QXYZ = W * MOTION(N)%Vy
          CASE(3);  QXYZ = W * MOTION(N)%Vz
        END SELECT
        VXYZ = VXYZ + QXYZ
      ENDDO

      WEDGE_NP_VELOCITY = VXYZ

      RETURN
      END FUNCTION WEDGE_NP_VELOCITY
!!
      FUNCTION WEDGE_VOLUME( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Wedge ID
!!
!! Function return value.
      REAL(4) :: WEDGE_VOLUME
!!
!! This routine supplies the wedge element volume.
!!
      WEDGE_VOLUME = WEDGE(NEL)%Volume

      RETURN
      END FUNCTION WEDGE_VOLUME
!!
      FUNCTION WEDGE_CRITICAL_DT( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2004
!!
!! Supply the wedge element's slave parent's critical dt.
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Wedge ID
!!
!! Function return value.
      REAL(4) :: WEDGE_CRITICAL_DT
!!
!! Local variables.
      INTEGER(ITYPE) :: ParID, Ptype
!!
!! Find parennt element and type.
!!
      ParID = SEGMENT(WEDGE(NEL)%SlaID)%PAR%ParID
      Ptype = SEGMENT(WEDGE(NEL)%SlaID)%PAR%Ptype
!!
      WEDGE_CRITICAL_DT = ZERO

      SELECT CASE (Ptype)
        CASE ( Element_Type_HEXAH   );  WEDGE_CRITICAL_DT =   HEXAH(ParID)%RES%DTelt
        CASE ( Element_Type_PENTA   );  WEDGE_CRITICAL_DT =   PENTA(ParID)%RES%DTelt
        CASE ( Element_Type_PYRAMID );  WEDGE_CRITICAL_DT = PYRAMID(ParID)%RES%DTelt
        CASE ( Element_Type_TETRA   );  WEDGE_CRITICAL_DT =   TETRA(ParID)%RES%DTelt
        CASE ( Element_Type_MEMBT   );  WEDGE_CRITICAL_DT =   MEMBT(ParID)%RES%DTelt
        CASE ( Element_Type_PLATT   );  WEDGE_CRITICAL_DT =   PLATT(ParID)%RES%DTelt
        CASE ( Element_Type_MEMBQ   );  WEDGE_CRITICAL_DT =   MEMBQ(ParID)%RES%DTelt
        CASE ( Element_Type_PLATQ   );  WEDGE_CRITICAL_DT =   PLATQ(ParID)%RES%DTelt
      END SELECT

      RETURN
      END FUNCTION WEDGE_CRITICAL_DT
!!
      FUNCTION LINKED_PAIR_NODES_USED ( M )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2009
!!
!! Purpose: Mark linked-pair nodes used, count number of unique nodal
!! pairs and nodes used by the linked-pair set, and create index translation 
!! arrays for linked-pair set "M".
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: M   ! Node-to-node interface ID
!!
!! Function return value.
      INTEGER(ITYPE) :: LINKED_PAIR_NODES_USED
!!
!! Local Variables.
      INTEGER(ITYPE) :: N,K
!!
!! Clear marker/index-sequence and translation arrays.
!!
      NELUSED = 0
      NPTUSED = 0
      NPTNOWI = 0
!!
!! Mark linked-nodes used by this node-pair interface.
!!
      LINKED_PAIR_NODES_USED = 0

      IF (NUMC5 .GT. 0) THEN

        LINKED_PAIR_NODES_USED = 1

        DO N = 1,NUMC5
          NELUSED(N) = 1
          NPTUSED(LINKED_PAIR_INTERFACE(N)%Node_One) = 1
          NPTUSED(LINKED_PAIR_INTERFACE(N)%Node_Two) = 1
        ENDDO
!!
!! Create nodal point index map NPTNOWI for use with 
!! element connectivity n-tuples. (Maps program index 
!! to index written to esg files for the linked-pair set.)
!!
!! Also, convert NPTUSED into a sequential index map.
!!
        K = 0
        DO N = 1,NUMNP
          IF (NPTUSED(N) .EQ. 1) THEN
            K = K + 1
            NPTNOWI(N) = K
            NPTUSED(K) = N
          ENDIF
        ENDDO
        NNodes = K
!!
!! Convert NELUSED into a sequential index map.
!!
        K= 0 
        DO N = 1,NUMC5
          IF (NELUSED(N) .EQ. 1) THEN
            K = K + 1
            NELUSED(K) = N
          ENDIF
        ENDDO
        NElems = K

      ENDIF

      RETURN
      END FUNCTION LINKED_PAIR_NODES_USED
!!
      FUNCTION LINKED_PAIR_N_TUPLE( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2009
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! LINKED_PAIR_INTERFACE Index
!!
!! Function return value.
      INTEGER(ITYPE) :: LINKED_PAIR_N_TUPLE
!!
!! Local variables.
      INTEGER(ITYPE) :: I1,I2
!!
!! This routine supplies element connectivity n-tuples.
!!
      I1 = NPTNOWI(LINKED_PAIR_INTERFACE(NEL)%Node_One)
      I2 = NPTNOWI(LINKED_PAIR_INTERFACE(NEL)%Node_Two)

      NTUPLE(1:2) = (/I1,I2/)

      LINKED_PAIR_N_TUPLE = 2

      RETURN
      END FUNCTION LINKED_PAIR_N_TUPLE
!!
      FUNCTION LINKED_PAIR_VOLUME( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2009
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Linked-pair interface ID
!!
!! Function return value.
      REAL(4) :: LINKED_PAIR_VOLUME
!!
!! Local variables.
      INTEGER(ITYPE) :: I1,I2
      REAL(RTYPE)    :: dX,dY,dZ
!!
!! This routine supplies the linked-pair element separation.
!!
      N1 = LINKED_PAIR_INTERFACE(NEL)%Node_One
      N2 = LINKED_PAIR_INTERFACE(NEL)%Node_Two

      dX = MOTION(N2)%Px - MOTION(N1)%Px      
      dY = MOTION(N2)%Py - MOTION(N1)%Py      
      dZ = MOTION(N2)%Pz - MOTION(N1)%Pz      

      LINKED_PAIR_VOLUME = SQRT( dX*dX + dY*dY + dZ*dZ )

      RETURN
      END FUNCTION LINKED_PAIR_VOLUME
!!
      FUNCTION LINKED_PAIR_CRITICAL_DT( NEL )
!!
!! Copyright (c) by FMA Development, LLC, 26-MAY-2009
!!
!! Supply the linked-pair element's critical dt a zero.
!!
!! Arguments.
      INTEGER(ITYPE), INTENT(IN) :: NEL   ! Linked-pair interface ID
!!
!! Function return value.
      REAL(4) :: LINKED_PAIR_CRITICAL_DT
!!
      LINKED_PAIR_CRITICAL_DT = ZERO

      RETURN
      END FUNCTION LINKED_PAIR_CRITICAL_DT

      END SUBROUTINE WRITE_TO_ESG_RESULTS_FILE
