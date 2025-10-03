MODULE cable_mpi_stub_types_mod
  !* Stubs for MPI datatypes when compiling without MPI
  !
  ! Adapted from https://github.com/open-mpi/ompi/blob/v3.0.0/ompi/mpi/fortran/use-mpi-f08/mod/mpi-f08-types.F90
  IMPLICIT NONE

  TYPE, BIND(c) :: MPI_Status
    INTEGER :: MPI_SOURCE
    INTEGER :: MPI_TAG
    INTEGER :: MPI_ERROR
  END TYPE MPI_Status

  TYPE, bind(c) :: MPI_Comm
    INTEGER :: MPI_VAL
  END TYPE MPI_Comm

  TYPE, bind(c) :: MPI_Datatype
    INTEGER :: MPI_VAL
  END TYPE MPI_Datatype

  TYPE, bind(c) :: MPI_Errhandler
    INTEGER :: MPI_VAL
  END TYPE MPI_Errhandler

  TYPE, bind(c) :: MPI_File
    INTEGER :: MPI_VAL
  END TYPE MPI_File

  TYPE, bind(c) :: MPI_Group
    INTEGER :: MPI_VAL
  END TYPE MPI_Group

  TYPE, bind(c) :: MPI_Info
    INTEGER :: MPI_VAL
  END TYPE MPI_Info

  TYPE, bind(c) :: MPI_Message
    INTEGER :: MPI_VAL
  END TYPE MPI_Message

  TYPE, bind(c) :: MPI_Op
    INTEGER :: MPI_VAL
  END TYPE MPI_Op

  TYPE, bind(c) :: MPI_Request
    INTEGER :: MPI_VAL
  END TYPE MPI_Request

  TYPE, bind(c) :: MPI_Session
    INTEGER :: MPI_VAL
  END TYPE MPI_Session

  TYPE, bind(c) :: MPI_Win
    INTEGER :: MPI_VAL
  END TYPE MPI_Win

  TYPE(MPI_Comm), PARAMETER       :: MPI_COMM_NULL       = MPI_Comm(0)
  TYPE(MPI_Datatype), PARAMETER   :: MPI_DATATYPE_NULL   = MPI_Datatype(0)
  TYPE(MPI_Errhandler), PARAMETER :: MPI_ERRHANDLER_NULL = MPI_Errhandler(0)
  TYPE(MPI_File), PARAMETER       :: MPI_FILE_NULL       = MPI_File(0)
  TYPE(MPI_Group),  PARAMETER     :: MPI_GROUP_NULL      = MPI_Group(0)
  TYPE(MPI_Info), PARAMETER       :: MPI_INFO_NULL       = MPI_Info(0)
  TYPE(MPI_Message), PARAMETER    :: MPI_MESSAGE_NULL    = MPI_Message(0)
  TYPE(MPI_Op), PARAMETER         :: MPI_OP_NULL         = MPI_Op(0)
  TYPE(MPI_Request), PARAMETER    :: MPI_REQUEST_NULL    = MPI_Request(0)
  TYPE(MPI_Win), PARAMETER        :: MPI_WIN_NULL        = MPI_Win(0)

  TYPE(MPI_Datatype), PARAMETER :: MPI_INTEGER           = MPI_Datatype(0)
  TYPE(MPI_Datatype), PARAMETER :: MPI_INTEGER8          = MPI_Datatype(0)
  TYPE(MPI_Datatype), PARAMETER :: MPI_DOUBLE_PRECISION  = MPI_Datatype(0)
  TYPE(MPI_Datatype), PARAMETER :: MPI_DOUBLE_COMPLEX    = MPI_Datatype(0)
  TYPE(MPI_Datatype), PARAMETER :: MPI_2DOUBLE_PRECISION = MPI_Datatype(0)
  TYPE(MPI_Datatype), PARAMETER :: MPI_CHARACTER         = MPI_Datatype(0)
  TYPE(MPI_Datatype), PARAMETER :: MPI_LOGICAL           = MPI_Datatype(0)

  TYPE(MPI_Op), PARAMETER :: MPI_SUM      = MPI_Op(0)
  TYPE(MPI_Op), PARAMETER :: MPI_MINLOC   = MPI_Op(0)
  TYPE(MPI_Op), PARAMETER :: MPI_MAXLOC   = MPI_Op(0)
  TYPE(MPI_Op), PARAMETER :: MPI_LOR      = MPI_Op(0)
  TYPE(MPI_Op), PARAMETER :: MPI_LAND     = MPI_Op(0)
  TYPE(MPI_Op), PARAMETER :: MPI_MAX      = MPI_Op(0)
  TYPE(MPI_Op), PARAMETER :: MPI_MIN      = MPI_Op(0)
  TYPE(MPI_Op), PARAMETER :: MPI_IN_PLACE = MPI_Op(0)

  INTERFACE OPERATOR (.EQ.)
    MODULE PROCEDURE oct_mpi_comm_op_eq
    MODULE PROCEDURE oct_mpi_datatype_op_eq
    MODULE PROCEDURE oct_mpi_errhandler_op_eq
    MODULE PROCEDURE oct_mpi_file_op_eq
    MODULE PROCEDURE oct_mpi_group_op_eq
    MODULE PROCEDURE oct_mpi_info_op_eq
    MODULE PROCEDURE oct_mpi_message_op_eq
    MODULE PROCEDURE oct_mpi_op_op_eq
    MODULE PROCEDURE oct_mpi_request_op_eq
    MODULE PROCEDURE oct_mpi_win_op_eq
  END INTERFACE OPERATOR (.EQ.)

  INTERFACE OPERATOR (.NE.)
    MODULE PROCEDURE oct_mpi_comm_op_ne
    MODULE PROCEDURE oct_mpi_datatype_op_ne
    MODULE PROCEDURE oct_mpi_errhandler_op_ne
    MODULE PROCEDURE oct_mpi_file_op_ne
    MODULE PROCEDURE oct_mpi_group_op_ne
    MODULE PROCEDURE oct_mpi_info_op_ne
    MODULE PROCEDURE oct_mpi_message_op_ne
    MODULE PROCEDURE oct_mpi_op_op_ne
    MODULE PROCEDURE oct_mpi_request_op_ne
    MODULE PROCEDURE oct_mpi_win_op_ne
  END INTERFACE OPERATOR (.NE.)

CONTAINS

  LOGICAL FUNCTION oct_mpi_comm_op_eq(a, b)
    TYPE(MPI_Comm), intent(in) :: a, b
    oct_mpi_comm_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_comm_op_eq

  LOGICAL FUNCTION oct_mpi_datatype_op_eq(a, b)
    TYPE(MPI_Datatype), intent(in) :: a, b
    oct_mpi_datatype_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_datatype_op_eq

  LOGICAL FUNCTION oct_mpi_errhandler_op_eq(a, b)
    TYPE(MPI_Errhandler), intent(in) :: a, b
    oct_mpi_errhandler_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_errhandler_op_eq

  LOGICAL FUNCTION oct_mpi_file_op_eq(a, b)
    TYPE(MPI_File), intent(in) :: a, b
    oct_mpi_file_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_file_op_eq

  LOGICAL FUNCTION oct_mpi_group_op_eq(a, b)
    TYPE(MPI_Group), intent(in) :: a, b
    oct_mpi_group_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_group_op_eq

  LOGICAL FUNCTION oct_mpi_info_op_eq(a, b)
    TYPE(MPI_Info), intent(in) :: a, b
    oct_mpi_info_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_info_op_eq

  LOGICAL FUNCTION oct_mpi_message_op_eq(a, b)
    TYPE(MPI_Message), intent(in) :: a, b
    oct_mpi_message_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_message_op_eq

  LOGICAL FUNCTION oct_mpi_op_op_eq(a, b)
    TYPE(MPI_Op), intent(in) :: a, b
    oct_mpi_op_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_op_op_eq

  LOGICAL FUNCTION oct_mpi_request_op_eq(a, b)
    TYPE(MPI_Request), intent(in) :: a, b
    oct_mpi_request_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_request_op_eq

  LOGICAL FUNCTION oct_mpi_win_op_eq(a, b)
    TYPE(MPI_Win), intent(in) :: a, b
    oct_mpi_win_op_eq = (a%MPI_VAL .EQ. b%MPI_VAL)
  END FUNCTION oct_mpi_win_op_eq

  LOGICAL FUNCTION oct_mpi_comm_op_ne(a, b)
    TYPE(MPI_Comm), intent(in) :: a, b
    oct_mpi_comm_op_ne = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_comm_op_ne

  LOGICAL FUNCTION oct_mpi_datatype_op_ne(a, b)
    TYPE(MPI_Datatype), intent(in) :: a, b
    oct_mpi_datatype_op_ne = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_datatype_op_ne

  LOGICAL FUNCTION oct_mpi_errhandler_op_ne(a, b)
    TYPE(MPI_Errhandler), intent(in) :: a, b
    oct_mpi_errhandler_op_ne = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_errhandler_op_ne

  LOGICAL FUNCTION oct_mpi_file_op_ne(a, b)
    TYPE(MPI_File), intent(in) :: a, b
    oct_mpi_file_op_ne = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_file_op_ne

  LOGICAL FUNCTION oct_mpi_group_op_ne(a, b)
    TYPE(MPI_Group), intent(in) :: a, b
    oct_mpi_group_op_ne = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_group_op_ne

  LOGICAL FUNCTION oct_mpi_info_op_ne(a, b)
    TYPE(MPI_Info), intent(in) :: a, b
    oct_mpi_info_op_ne = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_info_op_ne

  LOGICAL FUNCTION oct_mpi_message_op_ne(a, b)
    TYPE(MPI_Message), intent(in) :: a, b
    oct_mpi_message_op_ne = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_message_op_ne

  LOGICAL FUNCTION oct_mpi_op_op_ne(a, b)
    TYPE(MPI_Op), intent(in) :: a, b
    oct_mpi_op_op_ne = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_op_op_ne

  LOGICAL FUNCTION oct_mpi_request_op_ne(a, b)
    TYPE(MPI_Request), intent(in) :: a, b
    oct_mpi_request_op_ne  = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_request_op_ne

  LOGICAL FUNCTION oct_mpi_win_op_ne(a, b)
    TYPE(MPI_Win), intent(in) :: a, b
    oct_mpi_win_op_ne = (a%MPI_VAL .NE. b%MPI_VAL)
  END FUNCTION oct_mpi_win_op_ne

END MODULE
