! Provides input parser in Fortran 90/95/03/08
! Inputs from command line are 'key=value'-pairs, and inputs from file are read in Simantics format.
! Can be used to set default values for input parameters, automatically process command line or input file
! parameters, and get the values of the parameters.
!
! If you use this code in a publication, please make a reference to:
! A. Penttilä, Fortran 95 Input Parser implementation (computer code),
! http://wiki.helsinki.fi/display/~aipentti@helsinki.fi/Collection+of+codes (2012).
!
! Update 2015, removed usage of ErrorCapsule-module. Errors are not handeled or
! issued in any way, but are still detected in the code. Add your own error handling if needed.
!
! Antti Penttilä
! 2012, 2015
! Department of Physics, University of Helsinki
!
! Modifications 2015 by Timo Väisänen
! -replaced some code
! -added errorHandler
! -commenting is done with # and not //
!
! NOTIFICATION!
! Not all errors go through errorHandler! This can result that
! program will run without knowing that something failed while
! reading input.
!
! Copyright (C) 2016 Antti Penttilä, Timo Väisänen and University of Helsinki
! All rights reserved.
! The new BSD License is applied to this software, see LICENSE.txt


MODULE InputParser
  ! Hashtable for storing key/value pairs of input
  USE Hashtable
    use error_handler
    use constants
  IMPLICIT NONE
 
  PRIVATE

  ! Maximum string lenght for input keyword, inherited from Hastable
  INTEGER, PARAMETER :: inp_key_length = hash_key_length
  ! Inherited variable type codes, inherited from Hastable
  INTEGER, PARAMETER :: inp_type_int = type_int, inp_type_real = type_real, inp_type_char = type_char, &
    inp_type_logical = type_logical, inp_type_string = type_string, inp_type_file_et_line = type_file_et_line
  ! Maximum string length for string value, inherited from Hastable
  INTEGER, PARAMETER :: inp_max_string_length = max_string_length
  ! Error codes
  INTEGER, PARAMETER :: inp_error_defined = -131, inp_error_wrong_type = error_wrong_type, inp_error_file_open = -132, &
    inp_command_line = -133

  ! Public interface, setup input key and give default value. This must be done for every possible input key
  ! before the input file is read.
  ! Call setup_input(key, default_value) where
  ! key is input keyword (string), and default value is integer, real, logical or character
  ! or
  ! setup_input(key, default_string, str_len)
  ! with string value and string length
  ! or
  ! setup_input(key, default_integer, default_string, str_len)
  ! with both string value and integer value (typically for file name and line number).
  INTERFACE setup_value
    MODULE PROCEDURE setup_int, setup_real, setup_char, setup_logical, setup_string, setup_file_et_line
  END INTERFACE

  ! Public interface, setup command line input key and give default value. This must be done for
  ! every possible input key before the command line is processed.
  ! Call setup_cmd_arg(key, default_value) where
  ! key is command line keyword (string), and default value is integer, real, logical or character
  ! or
  ! setup_cmd_arg(key, default_string, str_len)
  ! with string value and string length.
  ! Command line arguments must be given in 'key=value'-format.
  INTERFACE setup_cmd_arg
    MODULE PROCEDURE setup_cmd_int, setup_cmd_real, setup_cmd_char, setup_cmd_logical, &
      setup_cmd_string
  END INTERFACE
  
  ! Public interface, get the current value of the input key.
  ! Call get_value(key, value_data) where
  ! key is input keyword (string), and value_data is integer, real, logical or character
  ! or
  ! get_value(key, value_string, str_len)
  ! with string value and string length
  ! or
  ! get_value(key, value_integer, value_string, str_len)
  ! with both string value and integer value (typically for file name and line number).
  INTERFACE get_value
    MODULE PROCEDURE get_int, get_real, get_char, get_logical, get_string, get_file_et_line
  END INTERFACE

  ! Public interface, get the current value of the command line key.
  ! Call get_cmd_arg(key, value_data) where
  ! key is input keyword (string), and value_data is integer, real, logical or character
  ! or
  ! get_cmd_arg(key, value_string, str_len)
  ! with string value and string length.
  INTERFACE get_cmd_arg
    MODULE PROCEDURE get_cmd_int, get_cmd_real, get_cmd_char, get_cmd_logical, get_cmd_string
  END INTERFACE
  
  ! Public interface, read and process input file. All the possible
  ! input keys in the file must been set a default value beforehand.
  ! Call read_input(fn, DEBUG0) where
  ! fn is the filename of the input file, and if optional logical
  ! DEGUG is set .TRUE. (default is .FALSE.), debug information is printed to std_out
  ! while processing the input file
  ! or
  ! read_input(ch, DEBUG0, FN) where
  ! ch is read-cabable open unit. FN is optional string value that is given to
  ! file_et_line-type integer input as the file name if read from input.
  INTERFACE read_input
    MODULE PROCEDURE read_input_fn, read_input_ch
  END INTERFACE

  INTEGER :: i, j, k, astat
  REAL(kind=rk) :: x, y, z
  LOGICAL :: b, is_init = .FALSE.
  CHARACTER :: c
  CHARACTER(LEN=1024) :: err_str
  CHARACTER(LEN=inp_max_string_length) :: inp_str
  CHARACTER(LEN=1024) :: str
  
  PUBLIC :: init_input, setup_value, get_value, read_input, print_input_values, &
    setup_cmd_arg, get_cmd_arg, read_cmd_arg, &
    ! Codes
    inp_key_length, &
    ! Type codes
    inp_type_int, inp_type_real, inp_type_char, inp_type_logical, inp_type_string, inp_type_file_et_line, &
    ! Other
    inp_max_string_length, &
    ! Error codes
    inp_error_defined, inp_error_wrong_type, inp_error_file_open


CONTAINS
  
!  PUBLIC!!!!!!!!!!!!!!!!!!!!!
! Initialize input parser by calling this before
! anything else from this module is called.
SUBROUTINE init_input(n)
  IMPLICIT NONE
  ! Input, hashtable that is used to store input/value-pairs is
  ! inited with this number. Optimal case would be that n is somewhat larger,
  ! e.g. twice as large, as the number of inputs. Will work also if n is smaller,
  ! but not so efficiently.
  INTEGER, INTENT(IN) :: n

  CALL init_hash(n)
  is_init = .TRUE.
    
END SUBROUTINE init_input


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checks if input is initiated.
SUBROUTINE check_init()
  IMPLICIT NONE
  
  IF (is_init) RETURN
  
  WRITE(*,*) "Error from InputParser: subroutine 'init_input'", &
    " has to be called before any other routine is called."
  WRITE(*,*) "***Emergency call stop***"
  !call terminateProgram("InputParser: "//trim(err_str))
  
END SUBROUTINE check_init


! INTERFACE!!!!!!!!!!!!!!
! For generic interface setup_value
SUBROUTINE setup_int(key, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(IN) :: value_data

  CALL check_init()
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up the input value '", TRIM(key), "'"      
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, .TRUE.)
  END IF
END SUBROUTINE setup_int

! INTERFACE!!!!!!!!!!!!!!
! For generic interface setup_value
SUBROUTINE setup_real(key, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  REAL(kind=rk), INTENT(IN) :: value_data

  CALL check_init()
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up the input value '", TRIM(key), "'"    
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, .TRUE.)
  END IF
END SUBROUTINE setup_real

! INTERFACE!!!!!!!!!!!!!!
! For generic interface setup_value
SUBROUTINE setup_char(key, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER, INTENT(IN) :: value_data

  CALL check_init()
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up the input value '", TRIM(key), "'"      
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, .TRUE.)
  END IF
END SUBROUTINE setup_char

! INTERFACE!!!!!!!!!!!!!!!!
! For generic interface setup_value
SUBROUTINE setup_logical(key, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  LOGICAL, INTENT(IN) :: value_data

  CALL check_init()
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up the input value '", TRIM(key), "'"    
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, .TRUE.)
  END IF
END SUBROUTINE setup_logical


! INTERFACE!!!!!!!!!!!!!!!!!!!!
! For generic interface setup_value
SUBROUTINE setup_string(key, value_data, str_len)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER(LEN=*), INTENT(IN) :: value_data
  INTEGER, INTENT(IN) :: str_len

  CALL check_init()
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up the input value '", TRIM(key), "'"  
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, str_len, .TRUE.)
  END IF
END SUBROUTINE setup_string

! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For generic interface setup_value
SUBROUTINE setup_file_et_line(key, value_data, value_data2, str_len)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(IN) :: value_data
  CHARACTER(LEN=*), INTENT(IN) :: value_data2
  INTEGER, INTENT(IN) :: str_len

  CALL check_init()
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up the input value '", TRIM(key), "'"    
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, value_data2, str_len, .TRUE.)
  END IF
END SUBROUTINE setup_file_et_line


! INTERFACE!!!!!!!!!!!!!!!!!
! For generic interface setup_cmd_arg
SUBROUTINE setup_cmd_int(cmdkey, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  INTEGER, INTENT(IN) :: value_data
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up command line argument '", TRIM(key), "'"    
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, .TRUE.)
  END IF
END SUBROUTINE setup_cmd_int

! INTERFACE!!!!!!!!!!!!!!!!!!
! For generic interface setup_cmd_arg
SUBROUTINE setup_cmd_real(cmdkey, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  REAL(kind=rk), INTENT(IN) :: value_data
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up command line argument '", TRIM(key), "'"  
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, .TRUE.)
  END IF
END SUBROUTINE setup_cmd_real

! INTERFACE!!!!!!!!!!!!!!!!!!
! For generic interface setup_cmd_arg
SUBROUTINE setup_cmd_char(cmdkey, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  CHARACTER, INTENT(IN) :: value_data
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up command line argument '", TRIM(key), "'"
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, .TRUE.)
  END IF
END SUBROUTINE setup_cmd_char

! INTERFACE!!!!!!!!!!!!!!!!!!!
! For generic interface setup_cmd_arg
SUBROUTINE setup_cmd_logical(cmdkey, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  LOGICAL, INTENT(IN) :: value_data
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up command line argument '", TRIM(key), "'"  
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, .TRUE.)
  END IF
END SUBROUTINE setup_cmd_logical

! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!
! For generic interface setup_cmd_arg
SUBROUTINE setup_cmd_string(cmdkey, value_data, str_len)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  CHARACTER(LEN=*), INTENT(IN) :: value_data
  INTEGER, INTENT(IN) :: str_len
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  ! Already exists ?
  IF(i > 0) THEN
    ! Error
    WRITE(err_str,*) "Error in setting up command line argument '", TRIM(key), "'"    
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = add_hash_value(key, value_data, str_len, .TRUE.)
  END IF
END SUBROUTINE setup_cmd_string


! INTERFACE!!!!!!!!!!!!!
! For generic interface get_value
SUBROUTINE get_int(key, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(INOUT) :: value_data

  CALL check_init()
  
  i = hash_value_type(key)
  IF(i /= inp_type_int) THEN
    call addError("inputparser: Wrong or nonexistent value when asked for int: "//key)
    !call terminateProgram("InputParser: "//trim(err_str))
    ! Error "Wrong or nonexistent value when asked for int"
  ELSE
    i = get_hash_value(key,value_data)
  END IF

END SUBROUTINE get_int

! INTERFACE!!!!!!!!!!!!!
! For generic interface get_value
SUBROUTINE get_real(key, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  REAL(KIND=rk), INTENT(INOUT) :: value_data

  CALL check_init()
  
  i = hash_value_type(key)
  IF(i /= inp_type_real) THEN
    call addError("inputparser: Wrong or nonexistent value when asked for real: "//key)
    !call terminateProgram("InputParser: "//trim(err_str))
    !  Error "Wrong or nonexistent value when asked for real"
  ELSE
    i = get_hash_value(key,value_data)
  END IF

END SUBROUTINE get_real

! INTERFACE!!!!!!!!!!!!!
! For generic interface get_value
SUBROUTINE get_char(key, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER, INTENT(INOUT) :: value_data

  CALL check_init()
  
  i = hash_value_type(key)
  IF(i /= inp_type_char) THEN
    call addError("inputparser: Wrong or nonexistent value when asked for char: "//key)
    !call terminateProgram("InputParser: "//trim(err_str))
    !  Error "Wrong or nonexistent value when asked for char"
  ELSE
    i = get_hash_value(key,value_data)
  END IF

END SUBROUTINE get_char

! INTERFACE!!!!!!!!!!!!!!!
! For generic interface get_value
SUBROUTINE get_logical(key, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  LOGICAL, INTENT(INOUT) :: value_data

  CALL check_init()
  
  i = hash_value_type(key)
  IF(i /= inp_type_logical) THEN
    ! Error
    WRITE(str, *) "Wrong or nonexistent value when asked for logical for key '", TRIM(key), &
      "'. Status code was ", i
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = get_hash_value(key,value_data)
  END IF

END SUBROUTINE get_logical

! INTERFACE!!!!!!!!!!!!!!!!!!!
! For generic interface get_value
SUBROUTINE get_string(key, value_data, str_len)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER(LEN=*), INTENT(INOUT) :: value_data
  INTEGER, INTENT(INOUT) :: str_len

  CALL check_init()

  i = hash_value_type(key)
  IF(i /= inp_type_string) THEN     
    call addError("inputparser: Wrong or nonexistent value when asked for string: "//key)
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = get_hash_value(key, value_data, str_len)
  END IF

END SUBROUTINE get_string

! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For generic interface get_value
SUBROUTINE get_file_et_line(key, value_data, value_data2, str_len)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(INOUT) :: value_data
  CHARACTER(LEN=*), INTENT(INOUT) :: value_data2
  INTEGER, INTENT(INOUT) :: str_len
    
  CALL check_init()
  
  i = hash_value_type(key)
  IF(i /= inp_type_file_et_line) THEN
    call addError("inputparser: Wrong or nonexistent value when asked for string: "//key)
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = get_hash_value(key, value_data, value_data2, str_len)
  END IF

END SUBROUTINE get_file_et_line


! INTERFACE!!!!!!!!!!!!!!!!!
! For generic interface get_cmd_arg
SUBROUTINE get_cmd_int(cmdkey, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  INTEGER, INTENT(INOUT) :: value_data
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  IF(i /= inp_type_int) THEN
    call addError("inputparser: Wrong or nonexistent value when asked for int: "//key)
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = get_hash_value(key,value_data)
  END IF

END SUBROUTINE get_cmd_int

! INTERFACE!!!!!!!!!!!!!!!!!
! For generic interface get_cmd_arg
SUBROUTINE get_cmd_real(cmdkey, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  REAL(KIND=rk), INTENT(INOUT) :: value_data
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  IF(i /= inp_type_real) THEN
    call addError("inputparser: Wrong or nonexistent value when asked for real: "//key)
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = get_hash_value(key,value_data)
  END IF

END SUBROUTINE get_cmd_real

! INTERFACE!!!!!!!!!!!!!!!!!
! For generic interface get_cmd_arg
SUBROUTINE get_cmd_char(cmdkey, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  CHARACTER, INTENT(INOUT) :: value_data
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  IF(i /= inp_type_char) THEN
    call addError("inputparser: Wrong or nonexistent value when asked for char"//key)
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = get_hash_value(key,value_data)
  END IF

END SUBROUTINE get_cmd_char

! INTERFACE!!!!!!!!!!!!!!!!!!
! For generic interface get_cmd_arg
SUBROUTINE get_cmd_logical(cmdkey, value_data)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  LOGICAL, INTENT(INOUT) :: value_data
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  IF(i /= inp_type_logical) THEN
    ! Error
    WRITE(str, *) "Wrong or nonexistent value when asked for logical for key '", TRIM(key), &
      "'. Status code was ", i
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = get_hash_value(key,value_data)
  END IF

END SUBROUTINE get_cmd_logical

! INTERFACE!!!!!!!!!!!!!!!!!!!!!!
! For generic interface get_cmd_arg
SUBROUTINE get_cmd_string(cmdkey, value_data, str_len)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: cmdkey
  CHARACTER(LEN=*), INTENT(INOUT) :: value_data
  INTEGER, INTENT(INOUT) :: str_len
  CHARACTER(LEN=inp_key_length) :: key

  CALL check_init()
  
  WRITE(key, '(A,A)') TRIM(cmdkey), "-CMDARG"
  
  i = hash_value_type(key)
  IF(i /= inp_type_string) THEN
    call addError("inputparser: Wrong or nonexistent value when asked for string: "//key)
    !call terminateProgram("InputParser: "//trim(err_str))
  ELSE
    i = get_hash_value(key, value_data, str_len)
  END IF

END SUBROUTINE get_cmd_string


! INTERFACE!!!!!!!!!!!!!
! For generic interface read_input
SUBROUTINE read_input_fn(fn, DEBUG0)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: fn
  LOGICAL, INTENT(IN), OPTIONAL :: DEBUG0
  INTEGER :: ch, ln
  CHARACTER, PARAMETER :: tab = CHAR(9)
  CHARACTER(LEN=hash_key_length) :: key
  LOGICAL :: db
  
  CALL check_init()
  
  IF(PRESENT(DEBUG0)) THEN
    db = DEBUG0
  ELSE
    db = .FALSE.
  END IF
  
  ch = 40004
  OPEN(ch, FILE=TRIM(fn), ACTION='read', STATUS='old', IOSTAT=astat)
  IF(astat /= 0) THEN
    ! Error
    WRITE(err_str,*) "Error in opening file '", TRIM(fn), "' for reading"
    !call terminateProgram("InputParser: "//trim(err_str))
  END IF
  
  CALL read_input_ch(ch, DEBUG0=db)
  
  CLOSE(ch)
  
END SUBROUTINE read_input_fn

! INTERFACE!!!!!!!!!!!!!!!
! For generic interface read_input
SUBROUTINE read_input_ch(ch, DEBUG0, FN)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ch
  LOGICAL, INTENT(IN), OPTIONAL :: DEBUG0
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: FN
  INTEGER :: ln
  CHARACTER, PARAMETER :: tab = CHAR(9)
  CHARACTER(LEN=inp_max_string_length) :: lfn
  CHARACTER(LEN=hash_key_length) :: key
  LOGICAL :: db
  
  CALL check_init()
  
  IF(PRESENT(DEBUG0)) THEN
    db = DEBUG0
  ELSE
    db = .FALSE.
  END IF
  IF(PRESENT(FN)) THEN
    lfn = FN
  ELSE
    lfn = "unknown"
  END IF
  
  READ(ch, '(A)', IOSTAT=astat) inp_str
  ln = 0
  DO WHILE(astat == 0)
    ln = ln+1
    IF(db) THEN
      WRITE(*,*) "line ", ln, ":", TRIM(inp_str)
    END IF
    WRITE(inp_str, '(A)') ADJUSTL(inp_str)
    ! Find comments
    i = SCAN(inp_str, "#")
    IF(i == 1) THEN
      READ(ch, '(A)', IOSTAT=astat) inp_str
      CYCLE
    ELSE IF(i > 1) THEN
      WRITE(inp_str, '(A)') inp_str(1:i-1)
    END IF
    ! Find other
    CALL trim_string(inp_str, (/ "{", "}", ",", tab /))
    IF(LEN_TRIM(inp_str) == 0) THEN
      READ(ch, '(A)', IOSTAT=astat) inp_str
      CYCLE
    END IF
    ! Read key
    i = SCAN(inp_str, "=")
    ! No key?
    IF(i == 0) THEN
      READ(ch, '(A)', IOSTAT=astat) inp_str
      CYCLE
    END IF
    WRITE(key, '(A)') inp_str(1:i-1)
    ! Value type
    j = hash_value_type(key)
    SELECT CASE(j)
    CASE(inp_type_int)
      READ(inp_str(i+1:), *) k
      astat = change_hash_value(key, k, .TRUE.)
    CASE(inp_type_real)
      READ(inp_str(i+1:), *) x
      astat = change_hash_value(key, x, .TRUE.)
    CASE(inp_type_char)
      READ(inp_str(i+1:), '(A1)') c
      astat = change_hash_value(key, c, .TRUE.)
    CASE(inp_type_logical)
      READ(inp_str(i+1:), *) b
      astat = change_hash_value(key, b, .TRUE.)
    CASE(inp_type_string)
      READ(inp_str(i+1:), '(A)') str
      WRITE(str, '(A)') ADJUSTL(str)
      astat = change_hash_value(key, str, LEN_TRIM(str), .TRUE.)
    CASE(inp_type_file_et_line)
      WRITE(str, '(A)') ADJUSTL(str)
      astat = change_hash_value(key, ln, lfn, LEN_TRIM(lfn))
    CASE DEFAULT
      ! Error
      WRITE(err_str, *) "Wrong variable type: ", TRIM(key)
      !call terminateProgram("InputParser: "//trim(err_str))
    END SELECT
    IF(astat /= 0) THEN
      ! Error
      WRITE(err_str, *) "Wrong variable type: ", TRIM(key)
      !call terminateProgram("InputParser: "//trim(err_str))
    END IF

    READ(ch, '(A)', IOSTAT=astat) inp_str
  END DO
  
END SUBROUTINE read_input_ch


!  PUBLIC!!!!!!!!!!!!!
! Process command line arguments. All the possible command line keys
! must have been set default values. Command line keys adn values must
! be given as 'key=value'. If optional argument DEBUG0 is .TRUE. (default is .FALSE.),
! debugging information is printed in std_out.
SUBROUTINE read_cmd_arg(DEBUG0)
  IMPLICIT NONE
  ! Optional, set .TRUE. for debugging.
  LOGICAL, INTENT(IN), OPTIONAL :: DEBUG0
  INTEGER :: ai, ac, cal, i, ti, rstat
  CHARACTER, PARAMETER :: tab = CHAR(9)
  CHARACTER(LEN=inp_max_string_length) :: ca
  CHARACTER(LEN=hash_key_length) :: key
  LOGICAL :: db
  
  CALL check_init()
  
  IF(PRESENT(DEBUG0)) THEN
    db = DEBUG0
  ELSE
    db = .FALSE.
  END IF
  
  ac = COMMAND_ARGUMENT_COUNT()
  
  DO ai=1,ac
    CALL GET_COMMAND_ARGUMENT(ai,ca)
    cal = LEN_TRIM(ca)
    
    IF(db) WRITE(*,*) "processing '", ca(1:cal), "'..."

    i = SCAN(ca, "=")
    IF(i == 0) THEN
      ! Error
      WRITE(err_str, *) "Use '=' in command-line arguments:", TRIM(ca)
      !call terminateProgram("InputParser: "//trim(err_str))
    END IF

    WRITE(key, '(A,A)') ca(1:i-1), "-CMDARG"
    ti = hash_value_type(key)
    SELECT CASE(ti)
    CASE(inp_type_int)
      READ(ca(i+1:cal), *, IOSTAT=rstat) k
      astat = change_hash_value(key, k, .TRUE.)
    CASE(inp_type_real)
      READ(ca(i+1:cal), *, IOSTAT=rstat) x
      astat = change_hash_value(key, x, .TRUE.)
    CASE(inp_type_char)
      READ(ca(i+1:cal), *, IOSTAT=rstat) c
      astat = change_hash_value(key, c, .TRUE.)
    CASE(inp_type_logical)
      READ(ca(i+1:cal), *, IOSTAT=rstat) b
      astat = change_hash_value(key, b, .TRUE.)
    CASE(inp_type_string)
      READ(ca(i+1:cal), *, IOSTAT=rstat) str
      astat = change_hash_value(key, str, cal-1+1, .TRUE.)
    CASE DEFAULT
      ! Error
      WRITE(err_str, *) "Wrong variable type: ", TRIM(ca)
      !call terminateProgram("InputParser: "//trim(err_str))
    END SELECT
    IF(astat /= 0 .OR. rstat /= 0) THEN
      ! Error
      WRITE(err_str, *) "Wrong variable type: ", TRIM(ca)
      !call terminateProgram("InputParser: "//trim(err_str))
    END IF
    
  END DO
  
END SUBROUTINE read_cmd_arg


!  PUBLIC!!!!!!!!!!!!!!!!!!!!!
! Print all input keys and their current values,
! use for debugging.
SUBROUTINE print_input_values()
  IMPLICIT NONE
  CHARACTER(LEN=inp_key_length), DIMENSION(:), POINTER :: keyl
  INTEGER :: kn, typ, local_i, local_j, local_k
  
  CALL check_init()
  
  keyl => list_hash_keys()
  
  WRITE(*,*) ''
  WRITE(*,*) "---input_values---"
  kn = SIZE(keyl, 1)
  DO local_i=1,kn
    typ = hash_value_type(keyl(local_i))
    SELECT CASE(typ)
    CASE(inp_type_int)
      CALL get_value(keyl(local_i), local_j)
      WRITE(*,*) TRIM(keyl(local_i)), " (int): ", local_j
    CASE(inp_type_real)
      CALL get_value(keyl(local_i), x)
      WRITE(*,*) TRIM(keyl(local_i)), " (real): ", x
    CASE(inp_type_char)
      CALL get_value(keyl(local_i), c)
      WRITE(*,*) TRIM(keyl(local_i)), " (char): ", c
    CASE(inp_type_logical)
      CALL get_value(keyl(local_i), b)
      WRITE(*,*) TRIM(keyl(local_i)), " (boolean): ", b
    CASE(inp_type_string)
      CALL get_value(keyl(local_i), inp_str, local_j)
      WRITE(*,*) TRIM(keyl(local_i)), " (string): ", inp_str(1:local_j)
    CASE(inp_type_file_et_line)
      CALL get_value(keyl(local_i), local_k, inp_str, local_j)
      WRITE(*,*) TRIM(keyl(local_i)), " (fileetline): ", local_k, "(", inp_str(1:local_j), ")"
    END SELECT
  END DO

  WRITE(*,*) "---end input_values---"
  WRITE(*,*) ''
  
END SUBROUTINE print_input_values


!  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE trim_string(str, chars)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(INOUT) :: str
  CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: chars
  INTEGER :: ncodes, loci, ii, pn

  ncodes = SIZE(chars, 1)
  DO loci=1,ncodes
    pn = LEN(chars(loci))
    ii = SCAN(str, chars(loci))
    DO WHILE (ii > 0)
      IF(ii == 1) THEN
        WRITE(str, '(A)') str(2:)
      ELSE
        WRITE(str, '(A,A)') str(1:ii-1), str(ii+pn:)
      END IF
      ii = SCAN(str, chars(loci))
    END DO
  END DO
  WRITE(str, '(A)') TRIM(ADJUSTL(str))

END SUBROUTINE trim_string


END MODULE InputParser
