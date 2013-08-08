!    Copyright (C) 2006 Imperial College London and others.
!
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module arb_fun_diagnostics

  use diagnostic_source_fields
  use field_options
  use fields_manipulation
  use initialise_fields_module
  use fields
  use fldebug
  use global_parameters, only : timestep, OPTION_PATH_LEN, current_time
  use spud
  use state_fields_module
  use state_module
  use arbitrary_function



  implicit none
  
  private

  public :: calculate_arbitrary_function

  contains

  
  subroutine calculate_arbitrary_function(states,state_index,s_field)

    type(state_type), dimension(:) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout):: s_field
    integer :: i, stat, diagnostic_count, n
    type (function_data) :: fun_data


    call initialize_function_data(states,trim(s_field%option_path)//'/diagnostic/algorithm',&
         fun_data)

    call calculate_function_at_nodes(fun_data,s_field)

    call finalize_function_data(fun_data)
  end subroutine calculate_arbitrary_function

end module arb_fun_diagnostics
