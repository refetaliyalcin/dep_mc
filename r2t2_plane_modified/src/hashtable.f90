MODULE Hashtable
!! Provides hashtable in Fortran 90/95/03/08
!! Can store, modify and read values from
!! key/value pairs, where key is a string, and
!! value can have different datatypes.
!!
!! Antti Penttilä
!! 2012, 2015
!! Department of Physics, University of Helsinki
!!
!! Modifications 2015 by Timo Väisänen
!! -replaced some code(toLower(),kinds)
!!
!! Copyright (C) 2016 Antti Penttilä, Timo Väisänen and University of Helsinki
!! All rights reserved.
!! The new BSD License is applied to this software, see LICENSE.txt


    use constants
  IMPLICIT NONE

  PRIVATE
  ! max hash key length
  INTEGER, PARAMETER :: hash_key_length = 128
  ! max length for string value
  INTEGER, PARAMETER :: max_string_length = 1024
  ! type codes for hash node value
  INTEGER, PARAMETER :: type_int = 1, type_real = 2, type_char = 3, type_logical = 4, type_string = 5, type_file_et_line = 6
  ! error codes
  INTEGER, PARAMETER :: error_not_exists = -111, error_wrong_type = -112, error_allocation = -113, &
    error_use_bef_allocation = -114, error_string_too_long = -115, error_string_length_mismatch = -116, &
    stat_empty = 1
  
  ! Adds new key-value pairs to hashtable. Supports various data types.
  ! Always results status value of the operation, 0 if success
  ! Call interface with
  !   add_hash_value(key, value_data, do_empty) RESULT(stat)
  ! where
  !   key  --  key of the hash value, string, max. length given in hash_key_length
  !   value_data -- integer, real(rk), character or logical
  !   do_empty -- optional, but if given will create empty valued node
  ! or with
  !   add_hash_value(key, value_data, str_len, do_empty) RESULT(stat)
  ! where
  !   value_data is string and str_len is its length
  ! or with
  !   add_hash_value(key, value_data, value_data2, str_len, do_empty) RESULT(stat)
  ! where
  !  value_data is integer, value_data2 is string and str_len is its length. 
  INTERFACE add_hash_value
    MODULE PROCEDURE add_hash_int, add_hash_real, add_hash_char, add_hash_logical, add_hash_string, add_hash_file_et_line
  END INTERFACE

  ! Gets the value of key-value pair from hashtable. Supports various data types.
  ! Always results status value of the operation, 0 if success and stat_empty if empty node
  ! Call interface with
  !   get_hash_value(key, value_data) RESULT(stat)
  ! where
  !   key  --  key of the hash value, string, max. length given in hash_key_length
  !   value_data -- integer, real(kind=rk), character or logical
  ! or with
  !   get_hash_value(key, value_data, str_len) RESULT(stat)
  ! where
  !   value_data is string and str_len will be its length
  ! or with
  !   get_hash_value(key, value_data, value_data2, str_len) RESULT(stat)
  ! where
  !  value_data is integer, value_data2 is string and str_len is its length.
  INTERFACE get_hash_value
    MODULE PROCEDURE get_hash_int, get_hash_real, get_hash_char, get_hash_logical, get_hash_string, get_hash_file_et_line
  END INTERFACE

  ! Change the value of existing key-value pair in hashtable. Supports various data types.
  ! Always results status value of the operation, 0 if succes. If same_type is .TRUE. and the value data is not of the same
  ! type as the previous value, error error_wrong_type will be given to stat.
  ! Call interface with
  !   change_hash_value(key, value_data, same_type) RESULT(stat)
  ! where
  !   key  --  key of the hash value, string, max. length given in hash_key_length
  !   value_data -- integer, real(kind=rk), character or logical
  ! or with
  !   change_hash_value(key, value_data, str_len, same_type) RESULT(stat)
  ! where
  !   value_data is string and str_len will be its length
  ! or with
  !   change_hash_value(key, value_data, value_data2, str_len) RESULT(stat)
  ! where
  !  value_data is integer, value_data2 is string and str_len is its length.
  INTERFACE change_hash_value
    MODULE PROCEDURE change_hash_int, change_hash_real, change_hash_char, change_hash_logical, change_hash_string, &
      change_hash_file_et_line
  END INTERFACE

  INTEGER, PARAMETER :: stat_ok = 0
  INTEGER :: hash_table_length, astat, i, j, k
  LOGICAL :: inited = .FALSE.
  CHARACTER(LEN=512) :: err_msg
  
  TYPE hash_node
    ! Atomic types
    LOGICAL :: empty
    INTEGER :: int_value
    REAL(kind=rk) :: real_value
    CHARACTER :: char_value
    LOGICAL :: log_value
    CHARACTER(max_string_length) :: string_value
    ! Node internal variables
    CHARACTER(LEN=hash_key_length) :: hash_key
    INTEGER :: node_type, string_length
    TYPE(hash_node), POINTER :: next_node
  END TYPE hash_node
  
  TYPE hash_unit
    TYPE(hash_node), POINTER :: hash_object
  END TYPE hash_unit
  
  TYPE(hash_unit), DIMENSION(:), ALLOCATABLE :: hash_table
 

  PUBLIC :: init_hash, add_hash_value, get_hash_value, change_hash_value, &
    hash_value_type, delete_hash_node, list_hash_keys, &
    error_not_exists, error_wrong_type, error_allocation, error_use_bef_allocation, &
    error_string_too_long, error_string_length_mismatch, &
    type_int, type_real, type_char, type_logical, type_string, type_file_et_line, &
    hash_key_length, max_string_length, stat_empty
  

CONTAINS
  
! PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the hash table before using it.
! Call with maximum number of different hash keys n.
! Hashtable can store more than n key/value pairs, but
! it is more efficient when n is larger than the final
! number of pairs.
SUBROUTINE init_hash(n)
  ! Input, max. number of different hash keys
  INTEGER, INTENT(IN) :: n
  
  hash_table_length = n
  ALLOCATE(hash_table(hash_table_length), STAT=astat)
  IF (astat /= 0) THEN
    ! "Error in allocating Hashtable."
    RETURN
  END IF
  
  DO i=1,hash_table_length
    NULLIFY(hash_table(i)%hash_object)
  END DO
  
  inited = .TRUE.
  
END SUBROUTINE init_hash


! PRIVATE!!!!!!!!!!!!!!!!!!!!!!!
! Add a new node to first hash node
SUBROUTINE add_new_node_to_chain(prev_node, new_node)
  IMPLICIT NONE
  TYPE(hash_node), POINTER, INTENT(INOUT) :: prev_node, new_node

  ! Travel to the end of chain
  DO WHILE(ASSOCIATED(prev_node%next_node))
    prev_node => prev_node%next_node
  END DO
  prev_node%next_node => new_node
  
END SUBROUTINE add_new_node_to_chain


! PRIVATE!!!!!!!!!!!!!!!!
! Hash function for mapping the keys
FUNCTION hash_index(key_str) RESULT(key)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key_str
  INTEGER :: key
  INTEGER :: key_length
  CHARACTER :: c
  CHARACTER(LEN=hash_key_length) :: wkey_str
  INTEGER, DIMENSION(hash_key_length) :: key_iarray
  
  key_iarray(:) = 0
  wkey_str(:) = ''
  key_length = LEN_TRIM(key_str)
  IF(key_length > hash_key_length) THEN
    key_length = hash_key_length
  END IF
    
  WRITE(wkey_str,*) toLower(key_str(1:key_length))
  
  DO i=1,key_length
    c = wkey_str(i+1:i+1)
    key_iarray(i) = IACHAR(c)
  END DO
  
  key = MOD(SUM(key_iarray), hash_table_length)+1

END FUNCTION hash_index


! PRIVATE!!!!!!!!!!!!!!!!!!!!
! Find and travel to the node with the key
FUNCTION find_node(key, curr_node) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  TYPE(hash_node), POINTER, INTENT(INOUT) :: curr_node
  INTEGER :: stat, hash_i
  
  hash_i = hash_index(key)

  stat = error_not_exists
  curr_node => hash_table(hash_i)%hash_object
  ! Does node exists?
  IF( .NOT. ASSOCIATED(curr_node)) THEN
    RETURN
  END IF
  ! Check the key value
  DO WHILE(TRIM(curr_node%hash_key) /= TRIM(toLower(key)))
    IF( .NOT. ASSOCIATED(curr_node%next_node)) THEN
      RETURN
    END IF
    curr_node => curr_node%next_node
  END DO
  ! Now in right place
  stat = 0
  
END FUNCTION find_node


! PUBLIC!!!!!!!!!!!!!!!!!!!
! Deletes (and deallocates) the node with a
! give key
FUNCTION delete_hash_node(key) RESULT(stat)
  IMPLICIT NONE
  ! key-value of the pair, input
  CHARACTER(LEN=*), INTENT(IN) :: key
  ! result, status of operation
  INTEGER :: stat
  INTEGER :: hash_i
  TYPE(hash_node), POINTER :: curr_node, prev_node
  
  hash_i = hash_index(key)

  stat = error_not_exists
  curr_node => hash_table(hash_i)%hash_object
  prev_node => NULL()
  ! Does node exists?
  IF( .NOT. ASSOCIATED(curr_node)) THEN
    RETURN
  END IF
  ! Is only node?
  IF(ASSOCIATED(curr_node%next_node)) THEN
    DO WHILE(curr_node%hash_key /= toLower(TRIM(key)))
      IF( .NOT. ASSOCIATED(curr_node%next_node)) THEN
        RETURN
      END IF
      prev_node => curr_node
      curr_node => curr_node%next_node
    END DO
  END IF
  ! Now in right place
  
  ! Four possibilities - find the correct one
  ! First in hash chain
  IF(.NOT. ASSOCIATED(prev_node)) THEN
    ! The only node in hash chain
    IF(.NOT. ASSOCIATED(curr_node%next_node)) THEN
      NULLIFY(hash_table(hash_i)%hash_object)
    ! The first but not only in chain
    ELSE
      hash_table(hash_i)%hash_object => curr_node%next_node
    END IF
  ! Not first in chain
  ELSE
    ! Not the first nor the last
    IF(ASSOCIATED(curr_node%next_node)) THEN
      prev_node%next_node => curr_node%next_node
    ! The last in chain
    ELSE
      prev_node%next_node => NULL()
    END IF
  END IF
  
  DEALLOCATE(curr_node)
  stat = 0

END FUNCTION delete_hash_node


! PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gives the type of the hash value with a given key.
! Type codes are integers and are public variables in this
! module. Optional argument can be given to check
! if the hash value exists, but has empty value.
! Will result 0 hashvalue with key does not exist.
FUNCTION hash_value_type(key, is_empty) RESULT(stat)
  IMPLICIT NONE
  ! Key of the hashvalue
  CHARACTER(LEN=*), INTENT(IN) :: key
  ! If given, will be .TRUE. on exit if the key exists
  ! but the value is empty
  LOGICAL, INTENT(OUT), OPTIONAL :: is_empty
  ! Result, the type code of the value with the key
  ! or 0 if no such key/value
  INTEGER :: stat
  LOGICAL :: ch_empty = .FALSE.

  TYPE(hash_node), POINTER :: curr_node
  
  IF(PRESENT(is_empty)) THEN
    ch_empty = .TRUE.
  END IF

  stat = find_node(key, curr_node)

  ! No such node
  IF(stat /= 0) RETURN
  ! Node exists
  stat = curr_node%node_type
  IF(ch_empty) THEN
    is_empty = curr_node%empty
  END IF
  
END FUNCTION hash_value_type


! PUBLIC!!!!!!!!!!!!!!!!!!
! List all the hash keys that are in use.
! mainly for debugging
FUNCTION list_hash_keys() RESULT(keylist)
  IMPLICIT NONE
  ! Result, array of key-values (strings) that are in the hashtable
  CHARACTER(LEN=hash_key_length), DIMENSION(:), POINTER :: keylist
  TYPE(hash_node), POINTER :: curr_node
  
  i = 0
  ! Travel one time through to get the size
  DO j=1,hash_table_length
    IF(.NOT. ASSOCIATED(hash_table(j)%hash_object)) CYCLE
    i = i+1
    curr_node => hash_table(j)%hash_object
    ! Is only one in chain?
    DO WHILE(ASSOCIATED(curr_node%next_node))
      i = i+1
      curr_node => curr_node%next_node
    END DO
  END DO
  
  ALLOCATE(keylist(i))
  ! Travel second time through to get the keys
  i = 0
  DO j=1,hash_table_length
    IF(.NOT. ASSOCIATED(hash_table(j)%hash_object)) CYCLE
    i = i+1
    curr_node => hash_table(j)%hash_object
    keylist(i) = curr_node%hash_key
    ! Is only one in chain?
    DO WHILE(ASSOCIATED(curr_node%next_node))
      i = i+1
      curr_node => curr_node%next_node
      keylist(i) = curr_node%hash_key
    END DO
  END DO
  
END FUNCTION list_hash_keys

  
! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface add_hash_value
FUNCTION add_hash_int(key, value_data, do_empty) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(IN) :: value_data
  LOGICAL, INTENT(IN), OPTIONAL :: do_empty
  INTEGER :: stat
  INTEGER :: hash_i
  TYPE(hash_node), POINTER :: new_node, prev_node
  
  hash_i = hash_index(key)
  
  ALLOCATE(new_node)
  IF(LEN_TRIM(key) > hash_key_length) THEN
    new_node%hash_key = toLower(key(1:hash_key_length))
  ELSE
    new_node%hash_key = toLower(TRIM(key))
  END IF
  new_node%node_type = type_int
  new_node%int_value = value_data
  IF(PRESENT(do_empty)) THEN
    new_node%empty = .TRUE.
  ELSE
    new_node%empty = .FALSE.
  END IF
  new_node%next_node => NULL()
  ! Is node the first in chain?
  IF(.NOT. ASSOCIATED(hash_table(hash_i)%hash_object)) THEN
    hash_table(hash_i)%hash_object => new_node
  ELSE
    prev_node => hash_table(hash_i)%hash_object
    CALL add_new_node_to_chain(prev_node, new_node)
  END IF
  
  stat = 0
  
END FUNCTION add_hash_int


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface add_hash_value
FUNCTION add_hash_real(key, value_data, do_empty) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  REAL(kind=rk), INTENT(IN) :: value_data
  LOGICAL, INTENT(IN), OPTIONAL :: do_empty
  INTEGER :: stat
  INTEGER :: hash_i
  TYPE(hash_node), POINTER :: new_node, prev_node
  
   hash_i = hash_index(key)
  
  ALLOCATE(new_node)
  IF(LEN_TRIM(key) > hash_key_length) THEN
    new_node%hash_key = toLower(key(1:hash_key_length))
  ELSE
    new_node%hash_key = toLower(TRIM(key))
  END IF
  new_node%node_type = type_real
  new_node%real_value = value_data
  IF(PRESENT(do_empty)) THEN
    new_node%empty = .TRUE.
  ELSE
    new_node%empty = .FALSE.
  END IF
  new_node%next_node => NULL()
  ! Is node the first in chain?
  IF(.NOT. ASSOCIATED(hash_table(hash_i)%hash_object)) THEN
    hash_table(hash_i)%hash_object => new_node
  ELSE
    prev_node => hash_table(hash_i)%hash_object
    CALL add_new_node_to_chain(prev_node, new_node)
  END IF
  
  stat = 0
  
END FUNCTION add_hash_real


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface add_hash_value
FUNCTION add_hash_char(key, value_data, do_empty) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER, INTENT(IN) :: value_data
  LOGICAL, INTENT(IN), OPTIONAL :: do_empty
  INTEGER :: stat
  INTEGER :: hash_i
  TYPE(hash_node), POINTER :: new_node, prev_node
  
   hash_i = hash_index(key)
  
  ALLOCATE(new_node)
  IF(LEN_TRIM(key) > hash_key_length) THEN
    new_node%hash_key = toLower(key(1:hash_key_length))
  ELSE
    new_node%hash_key = toLower(TRIM(key))
  END IF
  new_node%node_type = type_char
  new_node%char_value = value_data
  IF(PRESENT(do_empty)) THEN
    new_node%empty = .TRUE.
  ELSE
    new_node%empty = .FALSE.
  END IF
  new_node%next_node => NULL()
  ! IS node the first in chain?
  IF(.NOT. ASSOCIATED(hash_table(hash_i)%hash_object)) THEN
    hash_table(hash_i)%hash_object => new_node
  ELSE
    prev_node => hash_table(hash_i)%hash_object
    CALL add_new_node_to_chain(prev_node, new_node)
  END IF
  
  stat = 0
  
END FUNCTION add_hash_char


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface add_hash_value
FUNCTION add_hash_logical(key, value_data, do_empty) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  LOGICAL, INTENT(IN) :: value_data
  LOGICAL, INTENT(IN), OPTIONAL :: do_empty
  INTEGER :: stat
  INTEGER :: hash_i
  TYPE(hash_node), POINTER :: new_node, prev_node
  
   hash_i = hash_index(key)
  
  ALLOCATE(new_node)
  IF(LEN_TRIM(key) > hash_key_length) THEN
    new_node%hash_key = toLower(key(1:hash_key_length))
  ELSE
    new_node%hash_key = toLower(TRIM(key))
  END IF
  new_node%node_type = type_logical
  new_node%log_value = value_data
  IF(PRESENT(do_empty)) THEN
    new_node%empty = .TRUE.
  ELSE
    new_node%empty = .FALSE.
  END IF
  new_node%next_node => NULL()
  ! Is node the first in chain?
  IF(.NOT. ASSOCIATED(hash_table(hash_i)%hash_object)) THEN
    hash_table(hash_i)%hash_object => new_node
  ELSE
    prev_node => hash_table(hash_i)%hash_object
    CALL add_new_node_to_chain(prev_node, new_node)
  END IF
  
  stat = 0
  
END FUNCTION add_hash_logical


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface add_hash_value
FUNCTION add_hash_string(key, value_data, str_len, do_empty) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER(LEN=*), INTENT(IN) :: value_data
  INTEGER, INTENT(IN) :: str_len
  LOGICAL, INTENT(IN), OPTIONAL :: do_empty
  INTEGER :: stat
  INTEGER :: hash_i
  TYPE(hash_node), POINTER :: new_node, prev_node
  
   hash_i = hash_index(key)
  
  IF(str_len > max_string_length) THEN
    ! Error
    WRITE(err_msg,*) "String length ", str_len, " is too large for string hash value. Must be ", &
      max_string_length, " at max"
    RETURN
  END IF
  IF(str_len > LEN(value_data)) THEN
    ! Error
    WRITE(err_msg,*) "String length ", str_len, " larger than the actual string"
    RETURN
  END IF
  
  ALLOCATE(new_node)
  IF(LEN_TRIM(key) > hash_key_length) THEN
    new_node%hash_key = toLower(key(1:hash_key_length))
  ELSE
    new_node%hash_key = toLower(TRIM(key))
  END IF
  new_node%node_type = type_string
  new_node%string_value = value_data(1:str_len)
  new_node%string_length = str_len
  IF(PRESENT(do_empty)) THEN
    new_node%empty = .TRUE.
  ELSE
    new_node%empty = .FALSE.
  END IF
  new_node%next_node => NULL()
  ! IS node the first in chain?
  IF(.NOT. ASSOCIATED(hash_table(hash_i)%hash_object)) THEN
    hash_table(hash_i)%hash_object => new_node
  ELSE
    prev_node => hash_table(hash_i)%hash_object
    CALL add_new_node_to_chain(prev_node, new_node)
  END IF
  
  stat = 0
  
END FUNCTION add_hash_string


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface add_hash_value
FUNCTION add_hash_file_et_line(key, value_data, value_data2, str_len, do_empty) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(IN) :: value_data
  CHARACTER(LEN=*), INTENT(IN) :: value_data2
  INTEGER, INTENT(IN) :: str_len
  LOGICAL, INTENT(IN), OPTIONAL :: do_empty
  INTEGER :: stat
  INTEGER :: hash_i
  TYPE(hash_node), POINTER :: new_node, prev_node
  
   hash_i = hash_index(key)
  
  IF(str_len > max_string_length) THEN
    ! Error
    WRITE(err_msg,*) "String length ", str_len, " is too large for string hash value. Must be ", &
      max_string_length, " at max"
    RETURN
  END IF
  IF(str_len > LEN(value_data2)) THEN
    ! Error
    WRITE(err_msg,*) "String length ", str_len, " larger than the actual string"
    RETURN
  END IF
  
  ALLOCATE(new_node)
  IF(LEN_TRIM(key) > hash_key_length) THEN
    new_node%hash_key = toLower(key(1:hash_key_length))
  ELSE
    new_node%hash_key = toLower(TRIM(key))
  END IF
  new_node%node_type = type_file_et_line
  new_node%int_value = value_data
  new_node%string_value = value_data2(1:str_len)
  new_node%string_length = str_len
  IF(PRESENT(do_empty)) THEN
    new_node%empty = .TRUE.
  ELSE
    new_node%empty = .FALSE.
  END IF
  new_node%next_node => NULL()
  ! IS node the first in chain?
  IF(.NOT. ASSOCIATED(hash_table(hash_i)%hash_object)) THEN
    hash_table(hash_i)%hash_object => new_node
  ELSE
    prev_node => hash_table(hash_i)%hash_object
    CALL add_new_node_to_chain(prev_node, new_node)
  END IF
  
  stat = 0
  
END FUNCTION add_hash_file_et_line


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface get_hash_value
FUNCTION get_hash_int(key, value_data) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(INOUT) :: value_data
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node

    ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(curr_node%node_type /= type_int) THEN
    stat = error_wrong_type
    RETURN
  END IF
  ! Everything all right
  value_data = curr_node%int_value
  IF(curr_node%empty) THEN
    stat = stat_empty
  ELSE
    stat = 0
  END IF
  
END FUNCTION get_hash_int


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface get_hash_value
FUNCTION get_hash_real(key, value_data) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  REAL(kind=rk), INTENT(INOUT) :: value_data
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(curr_node%node_type /= type_real) THEN
    stat = error_wrong_type
    RETURN
  END IF
  ! Everything all right
  value_data = curr_node%real_value
  IF(curr_node%empty) THEN
    stat = stat_empty
  ELSE
    stat = 0
  END IF
  
END FUNCTION get_hash_real


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface get_hash_value
FUNCTION get_hash_char(key, value_data) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER, INTENT(INOUT) :: value_data
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(curr_node%node_type /= type_char) THEN
    stat = error_wrong_type
    RETURN
  END IF
  ! Everything all right
  value_data = curr_node%char_value
  IF(curr_node%empty) THEN
    stat = stat_empty
  ELSE
    stat = 0
  END IF
  
END FUNCTION get_hash_char


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface get_hash_value
FUNCTION get_hash_logical(key, value_data) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  LOGICAL, INTENT(INOUT) :: value_data
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(curr_node%node_type /= type_logical) THEN
    stat = error_wrong_type
    RETURN
  END IF
  ! Everything all right
  value_data = curr_node%log_value
  IF(curr_node%empty) THEN
    stat = stat_empty
  ELSE
    stat = 0
  END IF
  
END FUNCTION get_hash_logical


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface get_hash_value
FUNCTION get_hash_string(key, value_data, str_len) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER(LEN=*), INTENT(INOUT) :: value_data
  INTEGER, INTENT(INOUT) :: str_len
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(curr_node%node_type /= type_string) THEN
    stat = error_wrong_type
    RETURN
  END IF
  
  str_len = curr_node%string_length
  
  IF(LEN(value_data) < str_len) THEN
    ! Error
    WRITE(err_msg,*) "Length of string variable is not large enough"
    RETURN
  END IF

  ! Everything all right
  value_data = curr_node%string_value(1:str_len)
  IF(curr_node%empty) THEN
    stat = stat_empty
  ELSE
    stat = 0
  END IF
  
END FUNCTION get_hash_string


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface get_hash_value
FUNCTION get_hash_file_et_line(key, value_data, value_data2, str_len) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(INOUT) :: value_data
  CHARACTER(LEN=*), INTENT(INOUT) :: value_data2
  INTEGER, INTENT(INOUT) :: str_len
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(curr_node%node_type /= type_file_et_line) THEN
    stat = error_wrong_type
    RETURN
  END IF
  
  str_len = curr_node%string_length
  
  IF(LEN(value_data2) < str_len) THEN
    ! Return
    WRITE(err_msg,*) "Length of string variable is not large enough"
    RETURN
  END IF

  ! Everything all right
  value_data = curr_node%int_value
  value_data2 = curr_node%string_value(1:str_len)
  IF(curr_node%empty) THEN
    stat = stat_empty
  ELSE
    stat = 0
  END IF
  
END FUNCTION get_hash_file_et_line


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface change_hash_value
FUNCTION change_hash_int(key, value_data, same_type) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(IN) :: value_data
  LOGICAL, INTENT(IN) :: same_type
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(same_type) THEN
    IF(curr_node%node_type /= type_int) THEN
      stat = error_wrong_type
      RETURN
    END IF
  END IF
  ! change the value
  curr_node%int_value = value_data
  curr_node%node_type = type_int
  curr_node%empty = .FALSE.
  
  stat = 0    
  
END FUNCTION change_hash_int


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface change_hash_value
FUNCTION change_hash_real(key, value_data, same_type) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  REAL(kind=rk), INTENT(IN) :: value_data
  LOGICAL, INTENT(IN) :: same_type
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(same_type) THEN
    IF(curr_node%node_type /= type_real) THEN
      stat = error_wrong_type
      RETURN
    END IF
  END IF
  ! change the value
  curr_node%real_value = value_data
  curr_node%node_type = type_real
  curr_node%empty = .FALSE.
  
  stat = 0  
  
END FUNCTION change_hash_real


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface change_hash_value
FUNCTION change_hash_char(key, value_data, same_type) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER, INTENT(IN) :: value_data
  LOGICAL, INTENT(IN) :: same_type
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(same_type) THEN
    IF(curr_node%node_type /= type_char) THEN
      stat = error_wrong_type
      RETURN
    END IF
  END IF
  ! change the value
  curr_node%char_value = value_data
  curr_node%node_type = type_char
  curr_node%empty = .FALSE.
  
  stat = 0    
  
END FUNCTION change_hash_char


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface change_hash_value
FUNCTION change_hash_logical(key, value_data, same_type) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  LOGICAL, INTENT(IN) :: value_data
  LOGICAL, INTENT(IN) :: same_type
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(same_type) THEN
    IF(curr_node%node_type /= type_logical) THEN
      stat = error_wrong_type
      RETURN
    END IF
  END IF
  ! change the value
  curr_node%log_value = value_data
  curr_node%node_type = type_logical
  curr_node%empty = .FALSE.
  
  stat = 0
  
END FUNCTION change_hash_logical


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface change_hash_value
FUNCTION change_hash_string(key, value_data, str_len, same_type) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  CHARACTER(LEN=*), INTENT(IN) :: value_data
  INTEGER, INTENT(IN) :: str_len
  LOGICAL, INTENT(IN) :: same_type
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(same_type) THEN
    IF(curr_node%node_type /= type_string) THEN
      stat = error_wrong_type
      RETURN
    END IF
  END IF
  
  ! Test the value
  IF(str_len > max_string_length) THEN
    ! Return
    WRITE(err_msg,*) "String length ", str_len, " is too large for string hash value. Must be ", &
      max_string_length, " at max"
    RETURN
  END IF
  IF(str_len > LEN(value_data)) THEN
    ! Error
    WRITE(err_msg,*) "String length ", str_len, " larger than the actual string"
    RETURN
  END IF

  ! change the value
  curr_node%string_value = value_data(1:str_len)
  curr_node%string_length = str_len
  curr_node%node_type = type_string
  curr_node%empty = .FALSE.
  
  stat = 0
  
END FUNCTION change_hash_string


! INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Different versions of the generic interface change_hash_value
FUNCTION change_hash_file_et_line(key, value_data, value_data2, str_len) RESULT(stat)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: key
  INTEGER, INTENT(IN) :: value_data
  CHARACTER(LEN=*), INTENT(IN) :: value_data2
  INTEGER, INTENT(IN) :: str_len
  INTEGER :: stat
  TYPE(hash_node), POINTER :: curr_node
  
  ! Find correct place in the hash
  stat = find_node(key, curr_node)
  IF(stat /= 0) RETURN
  ! Now in correct place
  IF(curr_node%node_type /= type_file_et_line) THEN
    stat = error_wrong_type
    RETURN
  END IF
  
  ! Test the value
  IF(str_len > max_string_length) THEN
    ! Error
    WRITE(err_msg,*) "String length ", str_len, " is too large for string hash value. Must be ", &
      max_string_length, " at max"
    RETURN
  END IF
  IF(str_len > LEN(value_data2)) THEN
    ! Error
    WRITE(err_msg,*) "String length ", str_len, " larger than the actual string"
    RETURN
  END IF

  ! change the value
  curr_node%int_value = value_data
  curr_node%string_value = value_data2(1:str_len)
  curr_node%string_length = str_len
  curr_node%node_type = type_file_et_line
  curr_node%empty = .FALSE.
  
  stat = 0  
  
END FUNCTION change_hash_file_et_line



    !toLower
    !Convert uppercase letters to lowercase
    function toLower(str) result(retStr)
        character(*), intent(in) :: str
        character(len(str)) :: retStr
        integer :: i
        retStr=str
        do i=1,len(str)
            if(iachar(str(i:i))>=65 .and. iachar(str(i:I))<=90) then
                retStr(i:i) = achar(iachar(str(i:i))+32)
            else
                retStr(i:i)=str(i:i)
            endif
        enddo
    end function


END MODULE Hashtable
