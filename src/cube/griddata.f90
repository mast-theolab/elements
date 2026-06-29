module griddata
    use iso_fortran_env, only:real64

    character(len=42), dimension(*), parameter :: grid_keys = [ &
        'Grid type                         ', &  !  1.
        'Grid origin                       ', &  !  2.
        'N points                          ', &  !  3.
        'Step                              '  &  !  4.
        ]

end module griddata