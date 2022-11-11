# Guidelines for the scientific documentation

The scientific documentation is written in the source code. Each procedure and/or module and derived type declaration will contain some information about the science is solves. A procedure is:

 - a subroutine
 - a function
 - an interface

## Minimum information per procedure

Just underneath the first line defining a procedure, please include:

### Purpose

Start with a description of the procedure's purpose in English. 

Include what quantity(ies) the procedure is calculating and a short description of the relationship between the outputs and the inputs. It is fine to refer to a published paper as explanation of the relationship (e.g. "following Wang and Leuning 1998"). 

Describe all scientific methods used in the procedure.

Avoid using abbreviations (e.g. "maximum" not "max.").

### References

Add references to the published papers containing the description of the science and equations used in the procedure. 

The should show the name of the authors and the publication year in the rendered text (e.g. "Wang and Leuning 1998") and be a hyperlink to the paper using the DOI link as follow:

```
    [Wang and Leuning 1998](https://doi.org/10.1016/S0168-1923(98)00061-6)
```

Do not use the citation text! It is too long and unnecessary with a hyperlink.

### Description of the arguments

Add a *description* for each argument at the end of the line or on the line **after** the argument declaration.

If the argument is the direct equivalent of a variable in an equation in a published paper, add the *symbol* used in the equation (using LaTeX Maths) and the *reference* after the description within parentheses.

Add the *unit* at the end of the description within parentheses. Use LaTex Maths to format superscripts and subscripts. The unit should be the last item in the description line.

