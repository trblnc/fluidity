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


  module granular_flow

  use state_module
  use fields
  use field_options
  use spud
  use global_parameters, only: option_path_len
  use futils, only: int2str

  implicit none

  private

 real,parameter ::  DEFAULT_E=0.9, T_MIN=1.0e-3

  public :: update_distribution_function, calculate_solid_phase_pressure, clean_solid_phase_pressure,&
       calculate_solid_shear_viscosity, calculate_solid_energy_source, calculate_solid_absorbtion,&
       calculate_solid_diffusion, set_div_u




  contains

    subroutine calculate_solid_shear_viscosity(state,PhaseVolumeFractions,PhaseVolumeFractionsOld,u_diffusion)
      type(state_type), dimension(:), intent(in) ::  state
      real, dimension(:), target, intent(in) :: PhaseVolumeFractions, PhaseVolumeFractionsOld
      real, dimension(:,:,:,:), intent(inout) :: u_diffusion
      
      real :: d0
      real, dimension(:), allocatable :: mu
      type(scalar_field), pointer :: distribution_function, temp, rho
      real, dimension(:), pointer :: PhaseVolumeFraction,PhaseVolumeFractionOld
   
      character(len=OPTION_PATH_LEN) :: dist_type, option_path
      type(tensor_field), pointer ::viscosity
      type(mesh_type), pointer :: mesh
      real :: CoefficientOfRestitution


      integer, dimension(:), pointer :: Tnodes, MatNodes
      integer :: phase, idim, nodes,ele

      mesh=>extract_mesh(state(1),"PressureMesh_Discontinuous")

      do phase=1, size(state)
         option_path='/material_phase['//int2str(phase-1)//']/scalar_field::RadialDistributionFunction'
         if (have_option(trim(option_path))) then
            option_path='/material_phase['//int2str(phase-1)//']/multiphase_properties/drag/diameter'
            call get_option(trim(option_path),d0,default=0.001)
            option_path='/material_phase['//int2str(phase-1)//']/multiphase_properties/coefficient_of_restitution'
            call get_option(trim(option_path),&
                 CoefficientOfRestitution,default=DEFAULT_E)
            distribution_function=>extract_scalar_field(state(phase),"RadialDistributionFunction")
            temp=>extract_scalar_field(state(phase),"Temperature")
            rho=> extract_scalar_field(state(phase),"Density")
            viscosity=> extract_tensor_field(state(phase),"Viscosity")
            nodes=node_count(distribution_function)
            
            PhaseVolumeFraction=>PhaseVolumeFractions((phase-1)*nodes+1:phase*nodes)
            PhaseVolumeFractionOld=>PhaseVolumeFractionsOld((phase-1)*nodes+1:phase*nodes)

            u_diffusion(:,:,:,phase)=0.0
            call zero(viscosity)

            do ele=1,ele_count(temp)
               allocate(mu(ele_loc(temp,ele)))

               Tnodes=>ele_nodes(temp,ele)
               Matnodes=>ele_nodes(mesh,ele)

            mu(:)=4.0/3.0*rho%val(Tnodes)*(0.5*PhaseVolumeFraction(Tnodes)+0.5*PhaseVolumeFractionOld(Tnodes))*distribution_function%val(Tnodes)*d0*(1.0+CoefficientOfRestitution)*sqrt(abs(temp%val(Tnodes))/3.1415927)
            do idim=1,size(u_diffusion,2)
               u_diffusion(MatNodes,idim,idim,phase)=mu
               viscosity%val(idim,idim,Tnodes)=min(mu,1e6)
            end do

            deallocate(mu)
         end do
         END if
      end do

    end subroutine calculate_solid_shear_viscosity

    subroutine update_distribution_function(state, PhaseVolumeFractions,PhaseVolumeFractionsOld)

      type(state_type), dimension(:), intent(in) ::  state
      real, dimension(:), target, intent(in) :: PhaseVolumeFractions, PhaseVolumeFractionsOld

      integer :: phase, nodes
      type(scalar_field), pointer :: distribution_function
      real, dimension(:), pointer :: PhaseVolumeFraction,PhaseVolumeFractionOld

      character(len=OPTION_PATH_LEN) :: dist_type, option_path
      real :: MaxVolFrac


      do phase=1, size(state)
         option_path='/material_phase['//int2str(phase-1)//']/scalar_field::RadialDistributionFunction'
         if (have_option(trim(option_path))) then
            distribution_function=>extract_scalar_field(state(phase),"RadialDistributionFunction")
            nodes=node_count(distribution_function)
            PhaseVolumeFraction=>PhaseVolumeFractions((phase-1)*nodes+1:phase*nodes)
            PhaseVolumeFractionOld=>PhaseVolumeFractionsOld((phase-1)*nodes+1:phase*nodes)

            call get_option(trim(option_path)//'/function/name',dist_type)
            select case(dist_type)

            case('Gidaspow')

               call get_option(trim(option_path)//'/maximum_phase_volume_fraction',MaxVolFrac)
               distribution_function%val=(3.0/5.0)*1.0/(max(1.0e-16,1.0-(min(abs(0.0*PhaseVolumeFraction+1.0*PhaseVolumeFractionOld)/MaxVolFrac,1.0))**(1.0/3.0)))

            end select
         end if
      end do

    end subroutine update_distribution_function


    subroutine clean_solid_phase_pressure(solid_phase_pressure)

      type(scalar_field), dimension(:) :: solid_phase_pressure
      integer :: phase
      character(len=OPTION_PATH_LEN) :: option_path

      do phase=1,size(solid_phase_pressure)
         option_path='/material_phase['//int2str(phase-1)//']'
         if (have_option(trim(option_path)//'/scalar_field::RadialDistributionFunction')) then
            call deallocate(solid_phase_pressure(phase))
         end if
      end do
    end subroutine clean_solid_phase_pressure
      

    subroutine calculate_solid_phase_pressure(state,solid_phase_pressure,&
         PhaseVolumeFractions,PhaseVolumeFractionsOld)
      type(state_type), dimension(:) :: state
      type(scalar_field), dimension(:) :: solid_phase_pressure
      real, dimension(:), target, intent(in) :: PhaseVolumeFractions, PhaseVolumeFractionsOld
      type(scalar_field), pointer :: distribution_function, density, temperature, PVF,OldPVF
      type(scalar_field), pointer :: spp
      integer :: phase, nodes
      real, dimension(:), pointer :: PhaseVolumeFraction,OldPhaseVolumeFraction
      real :: CoefficientOfRestitution
      character(len=OPTION_PATH_LEN) :: option_path


      do phase=1,size(state)
         option_path='/material_phase['//int2str(phase-1)//']'

         if (have_option(trim(option_path)//'/scalar_field::RadialDistributionFunction')) then

            call get_option(trim(option_path)//'/multiphase_properties/coefficient_of_restitution',&
                 CoefficientOfRestitution,default=DEFAULT_E)
            distribution_function=>extract_scalar_field(state(phase),"RadialDistributionFunction")
            Temperature=>extract_scalar_field(state(phase),"Temperature")
            Density=>extract_scalar_field(state(phase),"Density")
            nodes=node_count(distribution_function)
            PVF=>extract_scalar_field(state(phase),"PhaseVolumeFraction")
            OldPVF=>extract_scalar_field(state(phase),"OldPhaseVolumeFraction")
            PhaseVolumeFraction=>PhaseVolumeFractions((phase-1)*nodes+1:phase*nodes)
            OldPhaseVolumeFraction=>PhaseVolumeFractionsOld((phase-1)*nodes+1:phase*nodes)

         !   PVF%val=PhaseVolumeFraction
          !  OldPVF%val=OldPhaseVolumeFraction
            call allocate(solid_phase_pressure(phase),Temperature%mesh,"SolidPhasePressure")

            solid_phase_pressure(phase)%val=density%val*temperature%val*&
                 (1.0+2.0*(1.0+CoefficientOfRestitution)*distribution_function%val*PhaseVolumeFraction)

            if (has_scalar_field(state(phase),"SolidPhasePressure")) then
               spp=>extract_scalar_field(state(phase),"SolidPhasePressure")
               spp%val=solid_phase_pressure(phase)%val
            end if

         end if

      end do
      

    end subroutine calculate_solid_phase_pressure

    subroutine calculate_solid_energy_source(state,phase,t_absorb,t_source)

      type(state_type),dimension(:) :: state
      integer :: phase
      real, dimension(:) ::  t_absorb,t_source

      type(tensor_field), pointer :: viscosity
      type(vector_field), pointer :: velocity, positions
      type(scalar_field), pointer :: T,  source_out

      type(scalar_field) :: mass, source

      type(element_type), pointer :: vel_shape, T_shape

      real, allocatable, dimension(:,:,:) :: dvel_shape,mu

      real, allocatable, dimension(:) :: detwei

      real, dimension(:,:), allocatable :: U
      real, dimension(:,:,:), allocatable :: DU

      integer :: idim, jdim,kdim, gi, ele
      integer, dimension(:), pointer :: nodes
      real :: vol

      velocity=>extract_vector_field(state(phase),"Velocity")
      viscosity=>extract_tensor_field(state(phase),"Viscosity")
      T=>extract_scalar_field(state,"Temperature")
      positions=>extract_vector_field(state,"Coordinate")

      allocate(detwei(ele_ngi(T,1)))
      allocate(dvel_shape(ele_loc(velocity, 1), ele_ngi(velocity,1), mesh_dim(t)))
      allocate(U(mesh_dim(t),ele_loc(velocity,1)))
      allocate(DU(mesh_dim(t),ele_ngi(T,ele),mesh_dim(t)))
      allocate(mu(mesh_dim(t),mesh_dim(t),ele_loc(t,1)))
      

      call allocate(mass,T%mesh,"CVMass")
      call allocate(source,T%mesh,"Source")

      call zero(mass)
      call zero(source)

      vel_shape=>ele_shape(velocity,1)
      T_shape=>ele_shape(T,1)

      do ele=1, ele_count(T)
         u=ele_val(velocity,ele)
         mu=ele_val(viscosity,ele)
         nodes=>ele_nodes(T,ele)
         call transform_to_physical(positions, ele, vel_shape, &
              & dshape = dvel_shape, detwei = detwei)

         forall (gi=1:ele_ngi(velocity,ele),idim=1:mesh_dim(t))
            DU(:,gi,idim)=matmul(U,dvel_shape(:,gi,idim))
         end forall

         mass%val(nodes)=mass%val(nodes)+matmul(T_shape%n,detwei)
         do idim=1,mesh_dim(t)
            do jdim=1,mesh_dim(t)
               do kdim=1,mesh_dim(t)
                  source%val(nodes)=source%val(nodes)&
                       +matmul(T_shape%n,detwei)*(2.0/3.0)&
                       *(DU(Idim,1,jdim)+DU(jdim,1,idim))&
                       *mu(jdim,kdim,:)*du(idim,1,kdim)
               end do
            end do
         end do
      end do

      source%val=source%val/mass%val+t_absorb*T_MIN
      t_source=source%val

      if (has_scalar_field(state(phase),"TemperatureSource")) then
         source_out=>extract_scalar_field(state(phase),"TemperatureSource")
         source_out%val=source%val
      end if

      call deallocate(source)
      call deallocate(mass)

      deallocate(detwei)
      deallocate(dvel_shape)
      deallocate(U)
      deallocate(DU)
      deallocate(mu)


    end subroutine calculate_solid_energy_source

    subroutine calculate_solid_absorbtion(state,TNew,t_absorb,solid_source)

      type(state_type) :: state
      real, dimension(:) ::  tnew,t_absorb
      type(scalar_field_pointer) :: solid_source

      type(scalar_field), pointer :: distribution_function, density,&
           T, epsilon, solid_pressure,absorb_out
      real :: CoefficientOfRestitution, d0
      character(len=OPTION_PATH_LEN) :: option_path


      allocate(solid_source%ptr)

      distribution_function=>extract_scalar_field(state,"RadialDistributionFunction")
      density=>extract_scalar_field(state,"Density")
      T=>extract_scalar_field(state,"Temperature")
      epsilon=>extract_scalar_field(state,"PhaseVolumeFraction")
      solid_pressure=>extract_scalar_field(state,"SolidPhasePressure")

      call allocate(solid_source%ptr,T%mesh,"SourceWithDivU")
      call zero(solid_source%ptr)

      T%val=TNew

      call get_option(trim(state%option_path)//'/multiphase_properties/coefficient_of_restitution',&
                 CoefficientOfRestitution,default=DEFAULT_E)
      option_path=trim(state%option_path)//'/multiphase_properties/drag/diameter'
            call get_option(trim(option_path),d0,default=0.001)

      t_absorb=2.0*(1.0-CoefficientOfRestitution**2)*abs(epsilon%val*density%val*distribution_function%val)*&
           (4.0/d0*sqrt(abs(T%val)/3.1415927))
         ! misssing a div u term here!
!                    + 2.0*drag_coeff
      
      if (has_scalar_field(state,"TemperatureAbsorption")) then
         absorb_out=>extract_scalar_field(state,"TemperatureAbsorption")
         absorb_out%val(:)=t_absorb
      end if


      solid_source%ptr%val=-2.0*(1.0-CoefficientOfRestitution**2)*epsilon%val*density%val*distribution_function%val*T%val+solid_pressure%val

    end subroutine calculate_solid_absorbtion


    subroutine calculate_solid_diffusion(state,phase,tdiffusion)
      type(state_type), dimension(:) :: state
      integer ::phase
      real, dimension(:,:,:) :: tdiffusion
      
      type(scalar_field), pointer :: distribution_function, density, T, epsilon
      real :: CoefficientOfRestitution, d0
      character(len=OPTION_PATH_LEN) :: option_path

      integer :: ele, idim
      integer , dimension(:), pointer :: mat_nodes, t_nodes
      
      type(mesh_type), pointer :: mat_mesh


      distribution_function=>extract_scalar_field(state(phase),"RadialDistributionFunction")
      density=>extract_scalar_field(state(phase),"Density")
      T=>extract_scalar_field(state(phase),"Temperature")
      epsilon=>extract_scalar_field(state(phase),"PhaseVolumeFraction")

      mat_mesh=>extract_mesh(state(1),"PressureMesh_Discontinuous")

      call get_option(trim(state(phase)%option_path)//'/multiphase_properties/coefficient_of_restitution',&
           CoefficientOfRestitution,default=DEFAULT_E)
      option_path=trim(state(phase)%option_path)//'/multiphase_properties/drag/diameter'
            call get_option(trim(option_path),d0,default=0.001)

      do ele=1,ele_count(T)
         t_nodes=>ele_nodes(T,ele)
         mat_nodes=>ele_nodes(mat_mesh,ele)

         do idim=1,mesh_dim(t)
            tdiffusion(mat_nodes,idim,idim)=4.0/3.0*abs(density%val(t_nodes)&
                 *epsilon%val(T_nodes)*d0*(1.0+CoefficientOfRestitution)&
                 *distribution_function%val(T_nodes))*sqrt(abs(T%val(t_nodes))/3.1415927)
         end do

      end do



    end subroutine calculate_solid_diffusion

    subroutine set_div_u(state,U, V, W, CV_NONODS, U_NONODS, NDIM, NPHASE, &
         CT, NCOLCT, FINDCT, COLCT)
      type(state_type), dimension(:) :: state
      INTEGER, intent( in ) :: CV_NONODS, U_NONODS, NDIM, NPHASE, NCOLCT
      REAL, DIMENSION( U_NONODS * NPHASE ), intent( in ) :: U, V, W  
      INTEGER, DIMENSION( CV_NONODS + 1 ), intent( in ) :: FINDCT
      INTEGER, DIMENSION( NCOLCT ), intent( in ) :: COLCT
      REAL, DIMENSION( NCOLCT * NDIM * NPHASE ), intent( in ) :: CT

    end subroutine set_div_u
             
  end module granular_flow
               
               
