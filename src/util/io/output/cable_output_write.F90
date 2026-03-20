! TODO(Sean): The preprocessor define ENFORCE_SINGLE_PRECISION is enabled
! temporarily to restore bitwise reproducibility with the previous output module
! which enforces double precision variables to be sampled using single precision
! arrays, and enforces writing both double and single precision data as single
! precision.
#define ENFORCE_SINGLE_PRECISION

submodule (cable_output_mod:cable_output_internal) cable_output_write

  use aggregator_mod, only: aggregator_int32_0d_t
  use aggregator_mod, only: aggregator_int32_1d_t
  use aggregator_mod, only: aggregator_int32_2d_t
  use aggregator_mod, only: aggregator_int32_3d_t
  use aggregator_mod, only: aggregator_real32_0d_t
  use aggregator_mod, only: aggregator_real32_1d_t
  use aggregator_mod, only: aggregator_real32_2d_t
  use aggregator_mod, only: aggregator_real32_3d_t
  use aggregator_mod, only: aggregator_real64_0d_t
  use aggregator_mod, only: aggregator_real64_1d_t
  use aggregator_mod, only: aggregator_real64_2d_t
  use aggregator_mod, only: aggregator_real64_3d_t
  use aggregator_mod, only: SINGLE_PRECISION_ENFORCED_KIND

  use cable_netcdf_mod, only: cable_netcdf_decomp_t

  use cable_grid_reductions_mod, only: grid_cell_average
  use cable_grid_reductions_mod, only: first_patch_in_grid_cell

  implicit none

contains

  module subroutine cable_output_write_variable(output_stream, output_variable, patch, landpt, frame, restart)
    type(cable_output_stream_t), intent(inout) :: output_stream
    type(cable_output_variable_t), intent(inout), target :: output_variable
    type(patch_type), intent(in), optional :: patch(:)
    type(land_type), intent(in), optional :: landpt(:)
    integer, intent(in), optional :: frame
    logical, intent(in), optional :: restart

    class(cable_netcdf_decomp_t), pointer :: decomp
    integer :: i, ndims
    logical :: restart_local
    character(128) :: variable_name

    integer(kind=int32), pointer :: write_buffer_int32_0d
    integer(kind=int32), pointer :: write_buffer_int32_1d(:)
    integer(kind=int32), pointer :: write_buffer_int32_2d(:, :)
    integer(kind=int32), pointer :: write_buffer_int32_3d(:, :, :)
    real(kind=real32),   pointer :: write_buffer_real32_0d
    real(kind=real32),   pointer :: write_buffer_real32_1d(:)
    real(kind=real32),   pointer :: write_buffer_real32_2d(:, :)
    real(kind=real32),   pointer :: write_buffer_real32_3d(:, :, :)
    real(kind=SINGLE_PRECISION_ENFORCED_KIND), pointer :: write_buffer_real64_0d
    real(kind=SINGLE_PRECISION_ENFORCED_KIND), pointer :: write_buffer_real64_1d(:)
    real(kind=SINGLE_PRECISION_ENFORCED_KIND), pointer :: write_buffer_real64_2d(:, :)
    real(kind=SINGLE_PRECISION_ENFORCED_KIND), pointer :: write_buffer_real64_3d(:, :, :)

    decomp => null()

    write_buffer_int32_0d  => null()
    write_buffer_int32_1d  => null()
    write_buffer_int32_2d  => null()
    write_buffer_int32_3d  => null()
    write_buffer_real32_0d => null()
    write_buffer_real32_1d => null()
    write_buffer_real32_2d => null()
    write_buffer_real32_3d => null()
    write_buffer_real64_0d => null()
    write_buffer_real64_1d => null()
    write_buffer_real64_2d => null()
    write_buffer_real64_3d => null()

    restart_local = .false.
    if (present(restart)) restart_local = restart

    if (.not. restart_local .and. output_variable%reduction_method /= "none") then
      if (.not. present(patch) .or. .not. present(landpt)) then
        call cable_abort("Optional arguments patch and landpt must be present for grid reductions", __FILE__, __LINE__)
      end if
    end if

    variable_name = output_variable%get_netcdf_name()
    if (restart_local) variable_name = output_variable%field_name

    select type (aggregator => output_variable%aggregator)
    type is (aggregator_int32_0d_t)
      if (output_variable%reduction_method /= "none") then
        call cable_abort("Grid cell reductions are not supported for scalar variables", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_abort("Distributed writes are not supported for scalar variables", __FILE__, __LINE__)
      end if
      write_buffer_int32_0d => aggregator%aggregated_data
      if (restart_local) write_buffer_int32_0d => aggregator%source_data
      if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_int32_0d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_int32_0d)
      end if
    type is (aggregator_int32_1d_t)
      if (restart_local) then
        write_buffer_int32_1d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_int32_1d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_abort("Reduction method grid_cell_average is not supported for integer variables", __FILE__, __LINE__)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_int32_1d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_int32_1d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
        call output_stream%output_file%write_darray( &
              var_name=variable_name, &
              values=write_buffer_int32_1d, &
              decomp=decomp, &
              fill_value=CABLE_OUTPUT_FILL_VALUE_INT32, &
              frame=frame)
      else if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_int32_1d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_int32_1d)
      end if
    type is (aggregator_int32_2d_t)
      if (restart_local) then
        write_buffer_int32_2d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_int32_2d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_abort("Reduction method grid_cell_average is not supported for integer variables", __FILE__, __LINE__)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_int32_2d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_int32_2d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
        call output_stream%output_file%write_darray( &
              var_name=variable_name, &
              values=write_buffer_int32_2d, &
              decomp=decomp, &
              fill_value=CABLE_OUTPUT_FILL_VALUE_INT32, &
              frame=frame)
      else if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_int32_2d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_int32_2d)
      end if
    type is (aggregator_int32_3d_t)
      if (restart_local) then
        write_buffer_int32_3d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_int32_3d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_abort("Reduction method grid_cell_average is not supported for integer variables", __FILE__, __LINE__)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_int32_3d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_int32_3d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
        call output_stream%output_file%write_darray( &
              var_name=variable_name, &
              values=write_buffer_int32_3d, &
              decomp=decomp, &
              fill_value=CABLE_OUTPUT_FILL_VALUE_INT32, &
              frame=frame)
      else if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_int32_3d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_int32_3d)
      end if
    type is (aggregator_real32_0d_t)
      if (output_variable%reduction_method /= "none") then
        call cable_abort("Grid cell reductions are not supported for scalar variables", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_abort("Distributed writes are not supported for scalar variables", __FILE__, __LINE__)
      end if
      write_buffer_real32_0d => aggregator%aggregated_data
      if (restart_local) write_buffer_real32_0d => aggregator%source_data
      if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real32_0d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real32_0d)
      end if
    type is (aggregator_real32_1d_t)
      if (restart_local) then
        write_buffer_real32_1d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_real32_1d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real32_1d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_1d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real32_1d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_1d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
        call output_stream%output_file%write_darray( &
              var_name=variable_name, &
              values=write_buffer_real32_1d, &
              decomp=decomp, &
              fill_value=CABLE_OUTPUT_FILL_VALUE_REAL32, &
              frame=frame)
      else if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real32_1d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real32_1d)
      end if
    type is (aggregator_real32_2d_t)
      if (restart_local) then
        write_buffer_real32_2d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_real32_2d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real32_2d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_2d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real32_2d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_2d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
        call output_stream%output_file%write_darray( &
              var_name=variable_name, &
              values=write_buffer_real32_2d, &
              decomp=decomp, &
              fill_value=CABLE_OUTPUT_FILL_VALUE_REAL32, &
              frame=frame)
      else if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real32_2d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real32_2d)
      end if
    type is (aggregator_real32_3d_t)
      if (restart_local) then
        write_buffer_real32_3d => aggregator%source_data
      else if (output_variable%reduction_method == "none") then
        write_buffer_real32_3d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real32_3d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_3d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real32_3d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real32_3d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
        call output_stream%output_file%write_darray( &
              var_name=variable_name, &
              values=write_buffer_real32_3d, &
              decomp=decomp, &
              fill_value=CABLE_OUTPUT_FILL_VALUE_REAL32, &
              frame=frame)
      else if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real32_3d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real32_3d)
      end if
    type is (aggregator_real64_0d_t)
      if (output_variable%reduction_method /= "none") then
        call cable_abort("Grid cell reductions are not supported for scalar variables", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
        call cable_abort("Distributed writes are not supported for scalar variables", __FILE__, __LINE__)
      end if
      write_buffer_real64_0d => aggregator%aggregated_data
#ifdef ENFORCE_SINGLE_PRECISION
      if (restart_local) then
        call output_stream%output_file%put_var(variable_name, real(aggregator%source_data, kind=SINGLE_PRECISION_ENFORCED_KIND))
        return
      end if
#else
      if (restart_local) write_buffer_real64_0d => aggregator%source_data
#endif
      if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real64_0d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real64_0d)
      end if
    type is (aggregator_real64_1d_t)
      if (restart_local) then
#ifdef ENFORCE_SINGLE_PRECISION
        if (output_variable%distributed) then
          call cable_output_decomp_associate(output_stream, output_variable, decomp)
          call output_stream%output_file%write_darray( &
                var_name=variable_name, &
                values=real(aggregator%source_data, kind=SINGLE_PRECISION_ENFORCED_KIND), &
                decomp=decomp, &
                fill_value=CABLE_OUTPUT_FILL_VALUE_REAL32, &
                frame=frame)
        else
          call output_stream%output_file%put_var(variable_name, real(aggregator%source_data, kind=SINGLE_PRECISION_ENFORCED_KIND))
        end if
        return
#else
        write_buffer_real64_1d => aggregator%source_data
#endif
      else if (output_variable%reduction_method == "none") then
        write_buffer_real64_1d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real64_1d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_1d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real64_1d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_1d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
#ifdef ENFORCE_SINGLE_PRECISION
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
#else
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
#endif
        call output_stream%output_file%write_darray( &
              var_name=variable_name, &
              values=write_buffer_real64_1d, &
              decomp=decomp, &
#ifdef ENFORCE_SINGLE_PRECISION
              fill_value=CABLE_OUTPUT_FILL_VALUE_REAL32, &
#else
              fill_value=CABLE_OUTPUT_FILL_VALUE_REAL64, &
#endif
              frame=frame)
      else if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real64_1d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real64_1d)
      end if
    type is (aggregator_real64_2d_t)
      if (restart_local) then
#ifdef ENFORCE_SINGLE_PRECISION
        if (output_variable%distributed) then
          call cable_output_decomp_associate(output_stream, output_variable, decomp)
          call output_stream%output_file%write_darray( &
                var_name=variable_name, &
                values=real(aggregator%source_data, kind=SINGLE_PRECISION_ENFORCED_KIND), &
                decomp=decomp, &
                fill_value=CABLE_OUTPUT_FILL_VALUE_REAL32, &
                frame=frame)
        else
          call output_stream%output_file%put_var(variable_name, real(aggregator%source_data, kind=SINGLE_PRECISION_ENFORCED_KIND))
        end if
        return
#else
        write_buffer_real64_2d => aggregator%source_data
#endif
      else if (output_variable%reduction_method == "none") then
        write_buffer_real64_2d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real64_2d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_2d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real64_2d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_2d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
#ifdef ENFORCE_SINGLE_PRECISION
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
#else
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
#endif
        call output_stream%output_file%write_darray( &
              var_name=variable_name, &
              values=write_buffer_real64_2d, &
              decomp=decomp, &
#ifdef ENFORCE_SINGLE_PRECISION
              fill_value=CABLE_OUTPUT_FILL_VALUE_REAL32, &
#else
              fill_value=CABLE_OUTPUT_FILL_VALUE_REAL64, &
#endif
              frame=frame)
      else if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real64_2d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real64_2d)
      end if
    type is (aggregator_real64_3d_t)
      if (restart_local) then
#ifdef ENFORCE_SINGLE_PRECISION
        if (output_variable%distributed) then
          call cable_output_decomp_associate(output_stream, output_variable, decomp)
          call output_stream%output_file%write_darray( &
                var_name=variable_name, &
                values=real(aggregator%source_data, kind=SINGLE_PRECISION_ENFORCED_KIND), &
                decomp=decomp, &
                fill_value=CABLE_OUTPUT_FILL_VALUE_REAL32, &
                frame=frame)
        else
          call output_stream%output_file%put_var(variable_name, real(aggregator%source_data, kind=SINGLE_PRECISION_ENFORCED_KIND))
        end if
        return
#else
        write_buffer_real64_3d => aggregator%source_data
#endif
      else if (output_variable%reduction_method == "none") then
        write_buffer_real64_3d => aggregator%aggregated_data
      else if (output_variable%reduction_method == "grid_cell_average") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real64_3d)
        call grid_cell_average( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_3d, &
              landpt=landpt, &
              patch=patch)
      else if (output_variable%reduction_method == "first_patch_in_grid_cell") then
        call cable_output_reduction_buffers_associate(output_variable, write_buffer_real64_3d)
        call first_patch_in_grid_cell( &
              input_array=aggregator%aggregated_data, &
              output_array=write_buffer_real64_3d, &
              landpt=landpt)
      else
        call cable_abort("Invalid reduction method", __FILE__, __LINE__)
      end if
      if (output_variable%distributed) then
#ifdef ENFORCE_SINGLE_PRECISION
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
#else
        call cable_output_decomp_associate(output_stream, output_variable, decomp)
#endif
        call output_stream%output_file%write_darray( &
              var_name=variable_name, &
              values=write_buffer_real64_3d, &
              decomp=decomp, &
#ifdef ENFORCE_SINGLE_PRECISION
              fill_value=CABLE_OUTPUT_FILL_VALUE_REAL32, &
#else
              fill_value=CABLE_OUTPUT_FILL_VALUE_REAL64, &
#endif
              frame=frame)
      else if (present(frame)) then
        call output_stream%output_file%inq_var_ndims(variable_name, ndims)
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real64_3d, &
              start=[(1, i = 1, ndims - 1), frame])
      else
        call output_stream%output_file%put_var( &
              var_name=variable_name, &
              values=write_buffer_real64_3d)
      end if
    class default
      call cable_abort("Unexpected aggregator type", __FILE__, __LINE__)
    end select

  end subroutine cable_output_write_variable

end submodule
