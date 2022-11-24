# Guidelines for the scientific documentation

The scientific documentation must be written inline in the source code. Each of the following components:

- subroutine
- function
- interface
- derived type declaration
- non-dummy variable or array declaration

must be accompanied by a comment describing the related science.

The inline documentation is processed by [FORD][FORD] which supports formatting with [Markdown][markdown].

## Minimum requirements

At a minimum a module, subroutine, function, or interface must contain the following standard header template 
just underneath the first defining line:

    !*## Purpose
    !
    ! Just underneath the first line defining a procedure, start with a 
    ! description of the procedure's purpose in English. 
    !
    !## Method
    !
    ! Include what quantity(ies) the procedure is calculating and a short 
    ! description of the relationship between the outputs and the inputs.
    ! It is fine to refer to a published paper as explanation of the 
    ! relationship (e.g. "following Wang and Leuning 1998"). 
    !
    ! Mention all scientific methods used in the procedure.
    !
    ! Avoid using abbreviations (e.g. write "maximum" not "max.").
    ! 
    ! Remember to keep all lines within 80 columns width.
    !
    !## References
    !
    ! Add references to the published papers containing the description of 
    ! the science and equations used in the procedure. 
    !
    ! The references should show the name of the authors and the publication 
    ! year in the rendered text (e.g. "Wang and Leuning 1998") and be a hyperlink
    ! to the paper using the DOI link as follow:
    !
    ! [Wang and Leuning 1998](https://doi.org/10.1016/S0168-1923(98)00061-6)
    !
    ! Do not use the citation text! It is too long and unnecessary with a hyperlink.
    !

### Description procedure arguments and public module variables

- Add a description for each argument at the end of the line or on the line after the argument declaration
- If the argument is the direct equivalent of a variable in an equation in a published paper, add the symbol used in the equation (using LaTeX Maths) and the reference, if it's ambiguous, after the description within parentheses
- Add the unit at the end of the description within parentheses. Use LaTex Maths to format superscripts and subscripts. The unit should be the last item in the description line

### Example of minimum requirements

```fortran
MODULE MIN_DOC_EXAMPLE
    !! This module is to show how to write a minimum documentation
    !! for CABLE.

    IMPLICIT NONE

CONTAINS
    SUBROUTINE SUB1_EXAMPLE(myarg)
    !!## Purpose
    !!
    !! This subroutine only serves to show the minimum documentation.  
    !!
    !!## Method
    !!
    !! It writes the argument to the standard output.
    !! 
    !!## References
    !!
    !! The minimum requirements for documentation are taken from:  
    !! [CABLE developer guide](https://cable-lsm.github.io/CABLE/developer_guide/doc_guide/science_doc/)

    ! Arguments
    INTEGER, INTENT(IN) :: myarg !! Example argument to show the documentation (-)

    ! Write the argument as-is to standard output
    PRINT*, myarg

    END SUBROUTINE SUB1_EXAMPLE
END MODULE MIN_DOC_EXAMPLE
```

## Optional information to include

There is obviously a lot of additional information that can be included with the documentation and it highly depends on the procedure or derived type that is being documented. Below are some ideas but the list is not exhaustive.

- **Headings and special formatting**. It is possible to use Markdown to format the documentation. For example, one can add headings or lists.
  
  For *procedures and derived types*, **do not use level-1 headings** as these are used for the procedure declaration line or the derived type declaration.

- **Equations**. You can insert equations with [LaTeX Mathematics][latex-maths] notations
- **Order of procedure**. You may want to explain the successive steps of the code. When doing so it might be useful to:
    - use unordered or ordered lists
    - place the list items for each step throughout the code *before* the code it refers too. Doing so, the code comments 
can then easily be used for the documentation as well. 

### Example of more complex documentation

```fortran
MODULE EXTENSIVE_DOC_EXAMPLE
    !*# Overview 
    ! This module is to show how to write more extensive documentation for CABLE.  
    ! It is using Markdown to format the documentation which is supported by FORD.

    IMPLICIT NONE

CONTAINS
    SUBROUTINE SUB1_EXAMPLE(myarg)
    !*## Purpose
    !
    ! This subroutine only serves to show some more extensive documentation.  
    ! 
    !## Method
    !
    ! The subroutine takes a float as input, prints its value then adds 2/3, prints it and returns it.  
    ! The equation used is:
    ! \[ myarg = myarg + !frac{2.}{3.} \]
    !
    !## References
    !
    ! The guidelines for documentation can be found in:  
    ! [CABLE developer guide](https://cable-lsm.github.io/CABLE/developer_guide/doc_guide/science_doc/)
    !

    ! Arguments
    REAL, INTENT(INOUT) :: myarg !! Example argument to show the documentation (-)

    !>### Order of procedure
    !> 1. Start an ordered list for my documentation.
    PRINT*, myarg

    !> 2. Continue the same list. The numbering is important here.
    myarg = myarg + 2./3.

    !> 3. Place each comment before the line of code it refers to.
    PRINT*, myarg

    END SUBROUTINE SUB1_EXAMPLE
END MODULE MIN_DOC_EXAMPLE
```

[FORD]: https://forddocs.readthedocs.io/en/latest/index.html
[markdown]: https://www.markdownguide.org/cheat-sheet/
[latex-maths]: https://en.wikibooks.org/wiki/LaTeX/Mathematics