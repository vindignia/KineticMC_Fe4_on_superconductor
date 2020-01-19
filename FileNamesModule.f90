module FileNames
    implicit none
    character(len=40),parameter   :: directory='output_directory/'
    character(len=40),parameter   :: label='ZZZ.dat'
    character(len=40),parameter   :: directory_field_values='field_data/'	! used to read nominal/true fields from file
    character(len=80)   :: filename_fields=trim(adjustl(directory_field_values))//'fields_prova.dat' ! short '...input_2.dat'
    character(len=80)   :: filename_1=trim(adjustl(directory))//'magnetization_curves/Fe4_'//trim(adjustl(label))
    character(len=80)   :: filename_eigenvalues=trim(adjustl(directory))//'eigenvalues/eig_'//trim(adjustl(label))
    character(len=80)   :: filename_transition_rates=trim(adjustl(directory))//'transition_rates/W_'//trim(adjustl(label))
end module FileNames