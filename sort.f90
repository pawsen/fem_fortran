module sort_array

  private
  public ::  heapsort,remove_dups, binarySearch, brute_force_search

  public :: heapsort_real, Bubble_Sort_real, Shell_Sort_real

contains

  subroutine heapsort(a)
    ! Fra http://rosettacode.org/wiki/Sorting_algorithms/Heapsort
    ! Sorterer integer arrays i stigende orden

    integer, intent(in out) :: a(0:)
    integer :: start, n, bottom
    integer :: temp

    n = size(a)
    do start = (n - 2) / 2, 0, -1
       call siftdown(a, start, n);
    end do

    do bottom = n-1, 1, -1
       temp = a(0)
       a(0) = a(bottom)
       a(bottom) = temp;
       call siftdown(a, 0, bottom)
    end do

  end subroutine heapsort

  subroutine siftdown(a, start, bottom)

    integer, intent(in out) :: a(0:)
    integer, intent(in) :: start, bottom
    integer :: child, root
    integer :: temp

    root = start
    do while(root*2 + 1 < bottom)
       child = root * 2 + 1

       if ((child + 1 < bottom) .and. (a(child) < a(child+1))) then
          child = child + 1
       end if

       if (a(root) < a(child)) then
          temp = a(child)
          a(child) = a (root)
          a(root) = temp
          root = child
       else
          return
       end if
    end do


  end subroutine siftdown

  subroutine heapsort_real(a)

    real(8), intent(in out) :: a(0:)

    integer :: start, n, bottom
    real(8) :: temp

    n = size(a)
    do start = (n - 2) / 2, 0, -1
       call siftdown_real(a, start, n)
    end do

    do bottom = n - 1, 1, -1
       temp = a(0)
       a(0) = a(bottom)
       a(bottom) = temp;
       call siftdown_real(a, 0, bottom)
    end do

  end subroutine heapsort_real

  subroutine siftdown_real(a, start, bottom)

    real(8), intent(in out) :: a(0:)
    integer, intent(in) :: start, bottom
    integer :: child, root
    real(8) :: temp

    root = start
    do while(root*2 + 1 < bottom)
       child = root * 2 + 1

       if ((child + 1 < bottom) .and. (a(child) < a(child+1))) then
          child = child + 1
       end if

       if (a(root) < a(child)) then
          temp = a(child)
          a(child) = a (root)
          a(root) = temp
          root = child
       else
          return
       end if
    end do

  end subroutine siftdown_real

  SUBROUTINE Bubble_Sort_real(a)
    !http://rosettacode.org/wiki/Sorting_algorithms/Bubble_sort#Fortran
    !bubble sort is generally considered to be the simplest sorting algorithm
    !Because of its abysmal O(n^2) performance, it is not used often for large (or even medium-sized) datasets

    REAL(8), INTENT(in out), DIMENSION(:) :: a
    REAL(8) :: temp
    INTEGER :: i, j
    LOGICAL :: swapped = .TRUE.

    DO j = SIZE(a)-1, 1, -1
       swapped = .FALSE.
       DO i = 1, j
          IF (a(i) > a(i+1)) THEN
             temp = a(i)
             a(i) = a(i+1)
             a(i+1) = temp
             swapped = .TRUE.
          END IF
       END DO
       IF (.NOT. swapped) EXIT
    END DO
  END SUBROUTINE Bubble_Sort_real

  SUBROUTINE Shell_Sort_real(a,index_list)
    !http://rosettacode.org/wiki/Sorting_algorithms/Shell_sort#Fortran
    ! bedre end buble sort

    IMPLICIT NONE
    INTEGER :: i, j, increment
    REAL(8) :: temp
    REAL(8), INTENT(in out) :: a(:)
    integer, intent(out) :: index_list(:)

    increment = SIZE(a) / 2
    DO WHILE (increment > 0)
       DO i = increment+1, SIZE(a)
          j = i
          temp = a(i)
          DO WHILE ((j >= increment+1) .AND. (a(MAX(j-increment,1)) > temp))
             a(j) = a(j-increment)
             index_list(j) = j-increment
             j = j - increment
          END DO
          a(j) = temp
          index_list(j) = i
       END DO
       IF (increment == 2) THEN
          increment = 1
       ELSE
          increment = increment * 5 / 11
       END IF
    END DO

  END SUBROUTINE Shell_Sort_real


  SUBROUTINE remove_dups(in_vec,out_vec,k)
    ! sort an integer array

    integer, intent(in) :: in_vec(:)
    integer, pointer :: out_vec(:)
    integer, allocatable :: res_vec(:)
    integer, intent(out) :: k                   ! The number of unique elements
    integer :: i, j

    allocate(res_vec(size(in_vec)))

    k = 1
    res_vec(1) = in_vec(1)
    outer: do i=2,size(in_vec)
       do j=1,k
          if (res_vec(j) == in_vec(i)) then
             ! Found a match so start looking again
             cycle outer
          end if
       end do
       ! No match found so add it to the output
       k = k + 1
       res_vec(k) = in_vec(i)
    end do outer

    ALLOCATE(out_vec(k))
    out_vec = rec_vec

  end subroutine remove_dups


  recursive function binarySearch (a, find_value) result (bsresult)
    ! Binary search, Recursive function, eg call itself. Retuns index of matching value

    !http://rosettacode.org/wiki/Binary_search#Fortran
    integer, intent(in) :: a(:), find_value
    integer          :: bsresult, mid

    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
    else if (a(mid) > find_value) then
       bsresult= binarySearch(a(:mid-1), find_value)
    else if (a(mid) < find_value) then
       bsresult = binarySearch(a(mid+1:), find_value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binarySearch


  subroutine brute_force_search(a, find_value,bsresult)
    ! search trough array by a for-loop. Returns index of matching value

    integer, intent(in) :: a(:), find_value
    integer, intent(out) :: bsresult

    do i = 1, num_elements
       if (a(i) .eq. target_value) then
          bsresult = i
          exit
       endif
    end do

  end subroutine brute_force_search

end module sort_array


!!$
!!$! http://www.andrewduey.com/cscs252d.htm
!!$
!!$! This was my final Project for my Fortran 90 course
!!$! Because I had such a hard time finding examples of various sorting algorithms on the internet
!!$! I have decided to make my program available to anyone who is interested.  
!!$! I make no promises or warranties other than it seems to work pretty well.  If for some reason 
!!$! you need to get a hold of me contact information can be found on my website at http://www.andrewduey.com
!!$!
!!$! This program will read from the file called datain.dat and save data to dataout.dat
!!$! This program will load and sort the information in datain.dat, then prompt the user
!!$! to select which type of sort they would like to perform.  After the program has sorted 
!!$! the list it will display the output on the screen as well as to the data file (without the 
!!$! headers) and then exit.
!!$! Programs is designed to work with lists up to 10,000 records long
!!$! Written By Andrew Duey for CS252d
!!$! Program last modified 5-7-04
!!$
!!$MODULE types
!!$  !Defines the person variable type
!!$  IMPLICIT NONE
!!$  TYPE :: person !Defines the person type
!!$     CHARACTER (LEN=12) :: first !a place to store the first name
!!$     CHARACTER (LEN=12) :: last !a place to store the last name
!!$     CHARACTER (LEN=11) :: ssn !a place to store the social security number
!!$     INTEGER :: account_Num !a place to store the account number 
!!$     REAL :: amount_owed  !a place to store the amount owed on the account
!!$     TYPE (person), POINTER :: next_value  !the pointer which we use for the linked list as part of the insertion sort
!!$  END TYPE person
!!$END MODULE types
!!$
!!$PROGRAM final_project
!!$  USE types !Use the custom variable type we defined above
!!$  IMPLICIT NONE
!!$  !Declare the variables we will be using here
!!$
!!$  INTEGER, PARAMETER :: MAX_SIZE = 10000 !delcare the max number of records for eash changing
!!$  TYPE(person), DIMENSION(MAX_SIZE) :: customer_list !Declare the customer list where all customer data will be stored
!!$  INTEGER :: option ! Delecare the variable to read the user response into
!!$  INTEGER :: nvals = 0 ! The number of values from records we're going to sort
!!$  INTEGER :: status !Used to check the status of read and write operations
!!$  INTEGER :: i = 0 !delcare our counter
!!$  ! Write out the main menu to the screen
!!$  WRITE(*,*) ' Welcome to the Duey List Sorter'
!!$  WRITE(*,*) ''
!!$  WRITE(*,*) 'Please select from the following options'
!!$  WRITE(*,*) ''
!!$  WRITE(*,*) '1) Sort using Bubble Sort'
!!$  WRITE(*,*) '2) Sort using Shell Sort'
!!$  WRITE(*,*) '3) Sort using Selection Sort'
!!$  WRITE(*,*) '4) Sort using Quick Sort'
!!$  WRITE(*,*) '5) Sort using Insertion Sort'
!!$  WRITE(*,*) ''
!!$  WRITE(*,*) 'Please type the number of the sort you wish to perform and press enter'
!!$  READ(*,*) option !read the desired action from the user
!!$  OPEN (UNIT=3, FILE='datain.dat', STATUS='OLD', ACTION='READ', IOSTAT=status) !opens the data file
!!$  openif: IF (status == 0 ) THEN
!!$     !Open worked
!!$     readloop: DO	!loop through the records and read them into the array
!!$        READ (3, 1001, IOSTAT=status) customer_list(nvals)%first, customer_list(nvals)%last,customer_list(nvals)%ssn,&
!!$             customer_list(nvals)%account_Num, customer_list(nvals)%account_Num, customer_list(nvals)%amount_owed
!!$1001    FORMAT (1X, A12, 1X, A12, 1X, A11, 1X, I10, 1X, F10.2)
!!$        IF ( status /=0 ) EXIT
!!$        nvals = nvals + 1
!!$     END DO readloop
!!$     readif: IF ( status > 0 ) THEN ! if there was a problem reading the file tell the user
!!$        WRITE (*,*) 'An error occured while reading line ', nvals + 1, ' ' , status
!!$     ELSE
!!$        !WRITE (*,*) nvals, ' was read sucessfully'
!!$     END IF readif
!!$
!!$
!!$     !now that we have the data loaded into the array we need to sort it how the user wants
!!$     SelectOption: SELECT CASE (option)
!!$     CASE (1) SelectOption
!!$        WRITE (*,*) 'Bubble Sort selected'
!!$        CALL sort_bubble (customer_list, nvals)
!!$     CASE (2) SelectOption
!!$        WRITE (*,*) 'Shell Sort selected'
!!$        CALL sort_shell (customer_list, nvals)
!!$     CASE (3) SelectOption
!!$        WRITE (*,*) 'Selection Sort selected'
!!$        CALL sort_selection (customer_list, nvals)
!!$     CASE (4) SelectOption
!!$        WRITE (*,*) 'Quick Sort selected'
!!$        CALL sort_quick (customer_list, nvals)
!!$     CASE (5) SelectOption
!!$        WRITE (*,*) 'Insertion Sort selected'
!!$        CALL sort_insertion (customer_list, nvals)
!!$     CASE DEFAULT
!!$        WRITE (*,*) 'Sort option not recogonized, original data written to output file'
!!$     END SELECT SelectOption
!!$     !when we are done sorting the array we will dump it to the screen and to a file
!!$     OPEN (UNIT=4, FILE='dataout.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=status) !Opens the data file for writing
!!$     WRITE (*,*) '' !puts a couple of blank lines at the top of the screen for estetics
!!$     WRITE (*,*) ''
!!$     WRITE (*,*) 'First Name   Last Name    SSN              Acct#    Balance'
!!$     WRITE (*,*) '-----------------------------------------------------------'
!!$     outputloop: DO i=0, nvals - 1, 1
!!$        WRITE (4, 1000) customer_list(i)%first, customer_list(i)%last, customer_list(i)%ssn, customer_list(i)%account_Num, & 
!!$             customer_list(i)%amount_owed !writes to the file
!!$        WRITE (*,1000) customer_list(i)%first, customer_list(i)%last, customer_list(i)%ssn, customer_list(i)%account_Num, & 
!!$             customer_list(i)%amount_owed !writes to the screen
!!$1000    FORMAT (1X, A12, 1X, A12, 1X, A11, 1X, I10, 1X, F10.2)
!!$     END DO outputloop
!!$     CLOSE (UNIT=4) !closes output files
!!$  ELSE openif
!!$     WRITE (*,*) 'error reading input file with error code=', status !reports file errors if any
!!$  END IF openif
!!$  !close the file here
!!$  CLOSE (UNIT=3)
!!$  STOP !ends program
!!$END PROGRAM final_project
!!$
!!$
!!$SUBROUTINE sort_selection (customer_list, nvals) 
!!$  !implements selection sort
!!$  USE types
!!$  LOGICAL, EXTERNAL :: gt_person !funtion to tell which person goes first
!!$  INTEGER :: i = 0
!!$  INTEGER :: j = 0
!!$  TYPE(person):: temp_person2 !temp variable
!!$  TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list !Define the customer list to handle whatever size is sent
!!$  INTEGER, INTENT(IN) :: nvals !grab the number of values from the calling code
!!$  ! Here is where we will do the selection sort
!!$  WRITE (*,*) 'Now doing selection sort'
!!$  sortloop1: DO i = 0, nvals - 2, 1
!!$     !this loop increments i which is our starting point for the comparison
!!$     sortloop2:DO j = i+1, nvals -1, 1
!!$        !this loop increments j which is the ending point for the comparison		
!!$        swapposition: IF ( gt_person(customer_list(i),customer_list(j)) )  THEN
!!$           !WRITE (*,*) 'WE SWAPED ', customer_list(i)%last , ' and ', customer_list(j)%last
!!$           !swap the name here
!!$           temp_person2 = customer_list(i)
!!$           customer_list(i) = customer_list(j)
!!$           customer_list(j) = temp_person2
!!$        END IF swapposition
!!$     END DO sortloop2
!!$  END DO sortloop1
!!$
!!$END SUBROUTINE sort_selection
!!$
!!$
!!$SUBROUTINE sort_quick (customer_list, nvals) 
!!$  !Sets up for the quick sort recursive method
!!$  USE types
!!$  LOGICAL, EXTERNAL :: gt_person !funtion to tell which person goes first
!!$  INTEGER :: i = 0
!!$  INTEGER :: j = 0
!!$  TYPE(person):: temp_person2
!!$  TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list !Define the customer list to handle whatever size is sent
!!$  INTEGER, INTENT(IN) :: nvals !grab the number of values from the calling code
!!$  ! Here is where we will do the selection sort
!!$  WRITE (*,*) 'Now doing Quick sort'
!!$  CALL qsRecursive(0, nvals-1, customer_list) !kicks off the recursive process
!!$
!!$END SUBROUTINE sort_quick
!!$
!!$
!!$RECURSIVE SUBROUTINE qsRecursive (lo, hi, customer_list)
!!$  !This is the actualy recursive portion of the quicksort
!!$  USE types
!!$  INTEGER :: pivotPoint
!!$  INTEGER, INTENT(IN) :: lo
!!$  INTEGER, INTENT(IN) :: hi
!!$  TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list
!!$  pivotPoint = qsPartition(lo, hi, customer_list); !basically all we do is find the pivot point, adjust elements, then call it again
!!$  IF (lo < pivotPoint) CALL qsRecursive(lo, pivotPoint -1, customer_list)
!!$  IF (pivotPoint < hi) CALL qsRecursive(pivotPoint + 1, hi, customer_list)
!!$END SUBROUTINE qsRecursive
!!$FUNCTION qsPartition (loin, hiin, customer_list)
!!$  !The partition portios of the Quick Sort is the must involved part
!!$  USE types
!!$  LOGICAL, EXTERNAL :: gt_person !funtion to tell which person goes first
!!$  TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list
!!$  INTEGER, INTENT(IN) :: loin
!!$  INTEGER:: lo !variable so we can manipulate the hi and lo values without changing things elsewhere in the program by reference
!!$  INTEGER, INTENT(IN) :: hiin
!!$  INTEGER:: hi !variable so we can manipulate the hi and lo values without changing things elsewhere in the program by reference
!!$  TYPE(person)::pivot !the temp location for the pivitoal element to which everything will be compaired
!!$  hi = hiin
!!$  lo = loin
!!$  pivot = customer_list(lo)
!!$  DO
!!$     IF (lo >= hi) EXIT !exit the loop when done
!!$     DO !move in from the right
!!$        IF ((gt_person(pivot, customer_list(hi))) .OR. (lo >= hi)) EXIT
!!$        hi = hi - 1
!!$     END DO
!!$     IF (hi /= lo) then !move the entry indexed by hi to left side of partition
!!$        customer_list(lo) = customer_list(hi) 
!!$        lo = lo + 1
!!$     END IF
!!$     DO !move in from the left
!!$        IF ((gt_person(customer_list(lo),pivot)) .OR. (lo >= hi)) EXIT
!!$        lo = lo + 1
!!$     END DO
!!$     IF (hi /= lo) then !move the entry indexed by hi to left side of partition
!!$        customer_list(hi) = customer_list(lo) 
!!$        hi = hi - 1
!!$     END IF
!!$  END DO
!!$  customer_list(hi) = pivot !put the pivot element back when we're done
!!$  qsPartition = hi !return the correct position of the pivot element
!!$END FUNCTION qsPartition
!!$
!!$
!!$SUBROUTINE sort_insertion (customer_list, nvals) 
!!$  !the sub that handles the insertion sort
!!$  USE types
!!$  LOGICAL, EXTERNAL :: gt_person !funtion to tell which person goes first
!!$  LOGICAL, EXTERNAL:: sp !Function call to tell if the it's the same person
!!$  TYPE(person):: temp_person2
!!$  TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list !Define the customer list to handle whatever size is sent
!!$  INTEGER, INTENT(IN) :: nvals !grab the number of values from the calling code
!!$  TYPE (person), POINTER :: ptr !declare the pointers we use to build and maintain the linked list
!!$  TYPE (person), POINTER :: ptr1
!!$  TYPE (person), POINTER :: ptr2
!!$  TYPE (person), POINTER :: tail
!!$  TYPE (person), POINTER :: head
!!$  INTEGER :: istat !a variable to hold status flags in
!!$  INTEGER :: counter !Used to count the number of records we've read
!!$  INTEGER :: i = 0 !ah, gold ole i
!!$  ! Here is where we will do the selection sort
!!$  WRITE (*,*) 'Now doing Inerstion sort using pointers'
!!$  input: DO !Input the values from the existing customer array
!!$     if (counter==nvals+1) EXIT
!!$     temp_person2 = customer_list(counter)
!!$     counter = counter + 1
!!$     ALLOCATE (ptr, STAT=istat)
!!$     ptr = temp_person2
!!$     !Now we find where to put it in the list
!!$
!!$     new: IF (.NOT. ASSOCIATED(head)) THEN !check to see if we need to start the list
!!$        !ADD to front of list
!!$        head => ptr	!place at front
!!$        tail => head	!tail points to new value
!!$        NULLIFY (ptr%next_value) !Nullify next ptr
!!$     ELSE !if the list already exists
!!$        !Values already in list. Check for location.
!!$        front: IF (gt_person(head,ptr )) THEN !if it belongs at the start of the list
!!$           !Add to front of list
!!$           ptr%next_value => head
!!$           head => ptr
!!$        ELSE IF ( gt_person(ptr, tail) .OR. (sp(ptr,tail)) ) THEN !if it belongs at the end of the list do that
!!$           !Add at end of list
!!$           tail%next_value => ptr
!!$           tail => ptr
!!$           NULLIFY ( tail%next_value)
!!$        ELSE !otherwise figure out where in the list it belongs
!!$           !Find Place to add value
!!$           ptr1 => head
!!$           ptr2 => ptr1%next_value
!!$           search: DO
!!$              IF ( (gt_person(ptr,ptr1) .OR. (sp(ptr,ptr1))) .AND. (gt_person(ptr2,ptr))) THEN
!!$                 !Insert Value Here
!!$                 ptr%next_value => ptr2
!!$                 ptr1%next_value => ptr
!!$                 EXIT search
!!$              END IF
!!$              ptr1 => ptr2
!!$              ptr2 => ptr2%next_value
!!$           END DO search
!!$        End IF front
!!$     END IF new
!!$  END DO input
!!$  !Now write out the data
!!$  ptr => head
!!$  i = -1
!!$  output: DO !this writes the data back to the same customer_list array for output by common functions
!!$     IF ( .NOT. ASSOCIATED(ptr)) EXIT
!!$     customer_list(i)=ptr	
!!$     ptr => ptr%next_value
!!$     i = i + 1
!!$  END DO output
!!$
!!$END SUBROUTINE sort_insertion
!!$
!!$
!!$SUBROUTINE sort_shell (customer_list, nvals) !-----needs numvals and custlist
!!$  !This is where we do the shell sort
!!$  USE types
!!$  LOGICAL, EXTERNAL :: gt_person !funtion to tell which person goes first
!!$  INTEGER :: i = 0
!!$  INTEGER :: j = 0
!!$  INTEGER :: increment = 3 !this is the increment which can be adjusted up or down depending on condition and size of dataset
!!$  TYPE(person):: temp_person2
!!$  TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list !Define the customer list to handle whatever size is sent
!!$  INTEGER, INTENT(IN) :: nvals !grab the number of values from the calling code
!!$  ! Here is where we will do the selection sort
!!$  WRITE (*,*) 'Now doing Shell sort'
!!$
!!$  outloop: DO
!!$     IF (increment==0) EXIT !check to make sure it's not time to end
!!$     sortloop1: DO i = 0, nvals - 1, 1
!!$        !this loop increments i which is our starting point for the comparison
!!$        j=i
!!$        temp_person2 = customer_list(i)
!!$        sortloop2:DO !here in the inner loop is where the comparisons happen
!!$           IF ((j<increment) .OR. (gt_person(temp_person2,customer_list(j-increment)))) EXIT
!!$           !this loop increments j which is the ending point for the comparison		
!!$           customer_list(j) = customer_list(j - increment)
!!$           j=j-increment
!!$        END DO sortloop2
!!$        customer_list(j)=temp_person2
!!$     END DO sortloop1
!!$     IF ((increment/2) /= 0) THEN !make adjustments up and down to the increment
!!$        increment = increment/2
!!$     ELSE IF	(increment==1) then
!!$        increment = 0
!!$     ELSE
!!$        increment=1;
!!$     END IF
!!$  END DO outloop
!!$END SUBROUTINE sort_shell
!!$
!!$
!!$SUBROUTINE sort_bubble (customer_list, nvals) !-----needs numvals and custlist
!!$  !this is where teh bubble sort is done
!!$  USE types
!!$  LOGICAL, EXTERNAL :: gt_person !funtion to tell which person goes first
!!$  INTEGER :: i = 0
!!$  INTEGER :: j = 0
!!$  TYPE(person):: temp_person2
!!$  TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list !Define the customer list to handle whatever size is sent
!!$  INTEGER, INTENT(IN) :: nvals !grab the number of values from the calling code
!!$  ! Here is where we will do the selection sort
!!$  WRITE (*,*) 'Now doing a bubble sort'
!!$  sortloop1: DO i = nvals -1, 0, -1 !basically we just loop through every element to compare it against every other element
!!$     !this loop increments i which is our starting point for the comparison
!!$     sortloop2:DO j = 1, i, 1
!!$        !this loop increments j which is the ending point for the comparison		
!!$        swapposition: IF ( gt_person(customer_list(j-1),customer_list(j)) )  THEN
!!$           !swap the name here
!!$           temp_person2 = customer_list(j-1)
!!$           customer_list(j-1) = customer_list(j)
!!$           customer_list(j) = temp_person2
!!$        END IF swapposition
!!$     END DO sortloop2
!!$  END DO sortloop1
!!$
!!$END SUBROUTINE sort_bubble
!!$
!!$
!!$LOGICAL FUNCTION gt_person (a, b) !Greater Than Person is what it stands for
!!$  ! This function takes a person as the argument and figurs out which sorts out first
!!$  !this is used as part of every sorting method we use
!!$  USE types
!!$  IMPLICIT NONE
!!$  TYPE(person), INTENT(IN)::a,b !grab the arguments and format them to make Fortran happy
!!$  gt_person = .FALSE.!if no other conditions are met then the second person comes first
!!$  !WRITE (*,*) a%last, ' ', b%last
!!$  last: IF (a%last==b%last) THEN !check the last name
!!$     first: IF (a%first==b%first) then !check first name if last is same
!!$        ssn: IF (a%ssn>b%ssn) then  !check SSN if both first and last are same
!!$           gt_person = .TRUE.
!!$           !WRITE (*,*) 'swapped because of ssn'
!!$        END IF ssn
!!$     else
!!$        first2: IF (LLT ( b%first , a%first)) then !checking first if last matches
!!$           gt_person = .TRUE.
!!$           !WRITE (*,*) 'swapped because of first'
!!$        END IF first2
!!$     END IF first
!!$  ELSE
!!$     last2: IF (LLT ( b%last , a%last)) THEN !if nothing else we just check last (LLT adjusts for case sensitivy)
!!$        gt_person = .TRUE.
!!$        !WRITE (*,*) 'swapped because of last'
!!$     END IF last2
!!$  END IF last
!!$END FUNCTION gt_person
!!$
!!$LOGICAL FUNCTION sp (a, b)
!!$  ! This function takes a person as the argument and figurs out if they are the same person
!!$  ! This function is only used by insertion sort where points make the == operator not work
!!$  USE types
!!$  IMPLICIT NONE
!!$  TYPE(person), INTENT(IN)::a,b !grab the arguments and format them
!!$  IF ((a%first==b%first) .AND. (a%last==b%last) .AND. (a%ssn==b%ssn)) THEN !check to see if the first, last, and ssn are the same
!!$     sp = .TRUE.
!!$  ELSE
!!$     sp = .FALSE.
!!$  END IF
!!$END FUNCTION sp
