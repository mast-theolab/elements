module propinfo
    !! Module providing basic information on supported quantities
    use workdata, only: PropertyDB
    use exception, only: BaseException, ArgumentError

    implicit none

    interface load_property_info
        module procedure load_propinfo_from_id, load_propinfo_from_tag
    end interface load_property_info

    private :: load_propinfo_from_id, load_propinfo_from_tag

contains

! ======================================================================

subroutine load_propinfo_from_id(property, identifier)
    !! Load property information corresponding to an integer identifier.
    !!
    !! Loads basic property information related to the identifier
    !! provided as integer in input.
    !!
    !! NOTE: Information contained on property are overwritten.
    !!
    !! CAUTION: If the identifier is not found, empty information is
    !! provided.  This should be checked carefully by calling routine.
    class(PropertyDB), intent(inout) :: property
    !! Property database.
    integer, intent(in) :: identifier
    !! Identifier of the property of interest.

    select case(identifier)
    case(101)
        property%label = 'eldip'
        property%name = 'electric dipole'
        property%unit = 'e.a0'
        property%pdim = [3]
        property%known = .true.
    case(111)
        property%label = 'eldip_v'
        property%name = 'electric dipole (vel.)'
        property%unit = 'e.a0'
        property%pdim = [3]
        property%known = .true.
    case(102)
        property%label = 'magdip'
        property%name = 'magnetic dipole'
        property%unit = 'hbar.e/m_e'
        property%pdim = [3]
        property%known = .true.
    case(103)
        property%label = 'poltens'
        property%name = 'polarizability tensor'
        property%unit = 'a0^3'
        property%pdim = [3, 3]
        property%known = .true.
    case(107)
        property%label = 'elquad'
        property%name = 'electric quadrupole'
        property%unit = 'e.a0^2'
        property%pdim = [3, 3]
        property%known = .true.
    case default
        property%label = 'N/A'
        property%name = 'unknown'
        property%unit = 'N/A'
        property%pdim = [1]
    end select

end subroutine load_propinfo_from_id

! ======================================================================

subroutine load_propinfo_from_tag(property, name, tag)
    !! Load property information based on its name (and tag).
    !!
    !! Loads basic property information related to the identifier
    !! provided as a string or a group of string in input.
    !! `name` can be directly the name of the property of itnerest
    !! or a group.  In the latter case, a `tag` name should be
    !! provided to fully characterize the property of interest.
    !!
    !! NOTE: Information contained on property are overwritten.
    !!
    !! CAUTION: If the identifier is not found, empty information is
    !! provided.  This should be checked carefully by calling routine.
    class(PropertyDB), intent(inout) :: property
    !! Property database.
    character(len=*), intent(in) :: name
    !! Name/group name of the property of interest.
    character(len=*), intent(in), optional :: tag
    !! Tag of the property of interest within group.

    select case(name)
    case default
        property%label = 'N/A'
        property%name = 'unknown'
        property%unit = 'N/A'
        property%pdim = [1]
    end select
    
end subroutine load_propinfo_from_tag

! ======================================================================

    
end module propinfo
