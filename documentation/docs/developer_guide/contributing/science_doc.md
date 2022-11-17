# Guidelines for the scientific documentation

The scientific documentation is written in the source code. Each procedure and/or module and derived type declaration will contain some information about the science it solves. A procedure is:

 - a subroutine
 - a function
 - an interface

The inline documentation is processed by [FORD][FORD] which supports formatting with [Markdown][markdown].

## Minimum information per procedure

### Example of minimum information

```fortran
MODULE MIN_DOC_EXAMPLE
    !! This module is to show how to write a minimum documentation
    !! for CABLE.

    IMPLICIT NONE

CONTAINS
    SUBROUTINE SUB1_EXAMPLE(myarg)
    !! This subroutine only serves to show the minimum documentation.  
    !! It writes the argument to the standard output.
    !! 
    !! The minimum requirements for documentation are taken from:  
    !! [CABLE developer guide](https://cable-lsm.github.io/CABLE/developer_guide/science_doc_guidelines/)

    ! Arguments
    INTEGER, INTENT(IN) :: myarg !! Example argument to show the documentation (-)

    ! Write the argument as-is to standard output
    PRINT*, myarg

    END SUBROUTINE SUB1_EXAMPLE
END MODULE MIN_DOC_EXAMPLE
```

### Purpose

Just underneath the first line defining a procedure, start with a description of the **procedure's purpose** in English. 

Include what quantity(ies) the procedure is calculating and a short description of the relationship between the outputs and the inputs. It is fine to refer to a published paper as explanation of the relationship (e.g. "following Wang and Leuning 1998"). 

**Mention all scientific methods** used in the procedure.

**Avoid using abbreviations** (e.g. "maximum" not "max.").

### References

**Add references** to the published papers containing the description of the science and equations used in the procedure. 

The references should show the name of the authors and the publication year in the rendered text (e.g. "Wang and Leuning 1998") and be a hyperlink to the paper using the DOI link as follow:

```
    [Wang and Leuning 1998](https://doi.org/10.1016/S0168-1923(98)00061-6)
```

*Do not use the citation text!* It is too long and unnecessary with a hyperlink.

### Description of the arguments

**Add a description for each argument** at the end of the line or on the line *after* the argument declaration.

If the argument is the direct equivalent of a variable in an equation in a published paper, **add the symbol used in the equation** (using LaTeX Maths) and the *reference* after the description within parentheses.

**Add the unit** *at the end* of the description within parentheses. Use LaTex Maths to format superscripts and subscripts. The unit should be the last item in the description line.

## Optional information to include

There are obviously a lot of additional information that can be included with the documentation and it highly depends on the procedure or derived type that is being documented. Below are some ideas but the list is not exhaustive.

### Example of more complex documentation

```fortran
MODULE EXTENSIVE_DOC_EXAMPLE
    !*# Overview 
    ! This module is to show how to write more extensive documentation for CABLE.  
    ! It is using Markdown to format the documentation which is supported by FORD.

    IMPLICIT NONE

CONTAINS
    SUBROUTINE SUB1_EXAMPLE(myarg)
    !* This subroutine only serves to show some more extensive documentation.  
    ! 
    ! The guilines for documentation can be found in:  
    ! [CABLE developer guide](https://cable-lsm.github.io/CABLE/developer_guide/science_doc_guidelines/)
    !
    ! The subroutine takes an float as input, print its value then adds 2/3 and prints it and returns it.  
    ! The equation used is:
    ! \[ myarg = myarg + !frac{2.}{3.} \]
    !

    ! Arguments
    REAL, INTENT(INOUT) :: myarg !! Example argument to show the documentation (-)

    !>### Order of procedure
    !> 1. Start an ordered list for my documentation.
    PRINT*, myarg

    !> 2. Continue the same list. The numbering is important here.
    myarg = myarg + 2./3.

    !> 3. Final step
    PRINT*, myarg

    END SUBROUTINE SUB1_EXAMPLE
END MODULE MIN_DOC_EXAMPLE
```

### Headings and special formatting

For longer documentation, it is possible to use Markdown to format the documentation. For example, one can add headings or lists.

For *procedures and derived types*, **do not use level-1 headings** as these are used for the procedure declaration line or the derived type declaration.

### Equations

FORD supports the use of LaTeX Maths notations to insert equations in the documentation

### Order of procedure

It is possible to add an order of procedure to explain the successive steps of the code. When doing so it might be useful to:

 - use unordered or ordered lists
 - place the list items for each step throughout the code *before* the code it refers too. Doing so, the code comments can then easily be used for the documentation as well. 

[FORD]: https://forddocs.readthedocs.io/en/latest/index.html
[markdown]: https://www.markdownguide.org/cheat-sheet/